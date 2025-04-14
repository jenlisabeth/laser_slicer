# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

bl_info = {
    "name": "Laser Slicer",
    "author": "Ryan Southall",
    "version": (1, 0, 0),
    "blender": (2, 80, 0),
    "location": "3D View > Tools Panel",
    "description": "Creates cross-sections and exports as SVG for laser cutting",
    "warning": "",
    "wiki_url": "tba",
    "tracker_url": "https://github.com/rgsouthall/laser_slicer/issues",
    "category": "Object"
}

import bpy, os, bmesh, numpy as np, math, logging
from bpy.props import FloatProperty, BoolProperty, EnumProperty, IntProperty, StringProperty, FloatVectorProperty

# Set up logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
logger.addHandler(logging.StreamHandler())

# UI row helper
def newrow(layout, label, obj, prop_name):
    row = layout.row()
    row.label(text=label)
    row.prop(obj, prop_name)

def perform_axis_slicing(aob, axis, settings, scale_mm):
    scene = bpy.context.scene
    dp = bpy.context.evaluated_depsgraph_get()

    # Axis selection logic
    axis_index = {'X': 0, 'Y': 1, 'Z': 2}[axis]
    normal_vector = (int(axis == 'X'), int(axis == 'Y'), int(axis == 'Z'))

    # Determine projection axes based on slicing axis
    proj_axes = [i for i in range(3) if i != axis_index]

    temp_mesh = aob.evaluated_get(dp).to_mesh()
    bm = bmesh.new()
    bm.from_mesh(temp_mesh)
    bm.transform(aob.matrix_world)
    aob.evaluated_get(dp).to_mesh_clear()
    aob.select_set(False)

    # Settings
    mat_thickness = settings.laser_slicer_material_thick / scale_mm
    cut_spacing = settings.laser_slicer_cut_spacing / scale_mm
    use_polygons = settings.laser_slicer_accuracy

    # Get slice range
    coords = [v.co[axis_index] for v in bm.verts]
    minv = min(coords)
    maxv = max(coords)
    slice_pos = minv + mat_thickness * 0.5

    slice_vpos = []
    vertex_blocks = []
    edge_blocks = []
    v_total, e_total = 0, 0
    vcount_list, ecount_list = [0], [0]
    vertex_index = np.empty(0, dtype=np.int8)

    # Output mesh object
    cob = None
    cob_exists = False
    for obj in bpy.context.scene.objects:
        if obj.get('Slices_' + axis):
            cob = obj
            cob_exists = True
            me = cob.data
            bpy.context.view_layer.objects.active = cob
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.mesh.delete(type='VERT')
            bpy.ops.object.mode_set(mode='OBJECT')
            break
    if not cob_exists:
        me = bpy.data.meshes.new('Slices_' + axis)
        cob = bpy.data.objects.new('Slices_' + axis, me)
        cob['Slices_' + axis] = 1

    # Slice loop
    slice_positions = []
    while slice_pos < maxv:
        slice_positions.append(slice_pos)
        axis_vector = tuple([slice_pos*x for x in normal_vector])
        cbm = bm.copy()
        cut = bmesh.ops.bisect_plane(
            cbm, 
            geom=cbm.edges[:] + cbm.faces[:], 
            dist=0, 
            plane_co=axis_vector,
            plane_no=normal_vector,
            clear_outer=False, 
            clear_inner=False
        )['geom_cut']

        verts = [v for v in cut if isinstance(v, bmesh.types.BMVert)]
        edges = [e for e in cut if isinstance(e, bmesh.types.BMEdge)]

        voffset = min(v.index for v in verts)

        v_block = np.array([v.co for v in verts]).flatten()
        e_index = np.array([[v.index - voffset + v_total for v in e.verts] for e in edges]).flatten()

        slice_vpos = np.append(slice_vpos, v_block)
        vertex_blocks.append([
            [(v.co - cob.location)[proj_axes[0]], (v.co - cob.location)[proj_axes[1]]] for v in verts
        ])
        edge_blocks.append([
            [[(v.co - cob.location)[proj_axes[0]], (v.co - cob.location)[proj_axes[1]]] for v in e.verts]
            for e in edges
        ])
        vertex_index = np.append(vertex_index, e_index)

        v_total += len(verts)
        e_total += len(edges)
        vcount_list.append(v_total)
        ecount_list.append(e_total)
        slice_pos += mat_thickness + cut_spacing
        cbm.free()

    bm.free()

    # Finalize mesh
    me.vertices.add(v_total)
    me.vertices.foreach_set('co', slice_vpos)
    me.edges.add(e_total)
    me.edges.foreach_set('vertices', vertex_index)

    if use_polygons:
        # Regenerate edges and sort them into continuous paths
        polygon_paths = []
        for i in range(len(vcount_list) - 1):
            vstart, vend = vcount_list[i], vcount_list[i+1]
            estart, eend = ecount_list[i], ecount_list[i+1]
            slice_edges = me.edges[estart:eend]
            edge_verts = [v for e in slice_edges for v in e.vertices]
            singles = [v for v in edge_verts if edge_verts.count(v) == 1]
            used = []
            path = []
            edge_chain = []

            if singles:
                e = [e for e in slice_edges if e.vertices[0] in singles or e.vertices[1] in singles][0]
                path.extend(e.vertices)
                edge_chain.append(e)
            else:
                e = slice_edges[0]
                path.extend(e.vertices)
                edge_chain.append(e)

            while len(edge_chain) < len(slice_edges):
                for e in [e for e in slice_edges if e not in edge_chain]:
                    if e.vertices[0] == path[-1] and e.vertices[1] not in path:
                        path.append(e.vertices[1])
                        edge_chain.append(e)
                        break
                    elif e.vertices[1] == path[-1] and e.vertices[0] not in path:
                        path.append(e.vertices[0])
                        edge_chain.append(e)
                        break
                    elif e.vertices[0] in path and e.vertices[1] in path:
                        edge_chain.append(e)
                        break
                else:
                    path.append(-1)

            polygon_paths.append([
                ([me.vertices[v].co[proj_axes[0]], me.vertices[v].co[proj_axes[1]]] if v >= 0 else v)
                for v in path
            ])
        vertex_blocks = polygon_paths

    if not cob_exists:
        bpy.context.scene.collection.objects.link(cob)

    bpy.context.view_layer.objects.active = cob
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.object.mode_set(mode='OBJECT')
    aob.select_set(True)
    bpy.context.view_layer.objects.active = aob

    logger.debug(f"[DEBUG] Axis {axis} slicing complete. Slices in world space: {slice_positions}")
    return {
        "vpos": slice_vpos,
        "vertex_blocks": vertex_blocks,
        "edge_blocks": edge_blocks,
        "vertex_index": vertex_index,
        "v_total": v_total,
        "e_total": e_total,
        "vcount_list": vcount_list,
        "ecount_list": ecount_list,
        "slice_positions": slice_positions
    }

def generate_svg_files(data, aob, axis, settings, scale_mm, all_slice_positions):
    vertex_blocks = data['vertex_blocks']
    edge_blocks = data['edge_blocks']

    # Settings
    mat_width = settings.laser_slicer_material_width
    mat_height = settings.laser_slicer_material_height
    cut_thickness = settings.laser_slicer_cut_thickness / scale_mm
    sep_files = settings.laser_slicer_separate_files
    use_polygons = settings.laser_slicer_accuracy
    dpi = settings.laser_slicer_dpi
    svg_pos = settings.laser_slicer_svg_position
    output_path = settings.laser_slicer_ofile
    cut_color = settings.laser_slicer_cut_colour
    line_thick = settings.laser_slicer_cut_line

    dpi_scale = dpi / 25.4
    scale = scale_mm * dpi_scale
    logger.debug(f"[DEBUG] Overall scale: {scale}")

    # Attempt to retrieve the slice object (cob) from the scene.
    slice_obj = None
    for obj in bpy.context.scene.objects:
        if obj.get('Slices_' + axis):
            slice_obj = obj
            break
    if slice_obj is None:
        # Fallback to aob if no slice object is found.
        slice_obj = aob

    # Setup SVG filenames
    if sep_files:
        if os.path.isdir(bpy.path.abspath(output_path)):
            filepaths = [os.path.join(bpy.path.abspath(output_path), f"{aob.name}_{axis}_{i}.svg") for i in range(len(vertex_blocks))]
        else:
            out_dir = os.path.dirname(bpy.path.abspath(output_path or bpy.data.filepath))
            name_base = bpy.path.display_name_from_filepath(output_path or aob.name)
            filepaths = [os.path.join(out_dir, f"{name_base}_{axis}_{i}.svg") for i in range(len(vertex_blocks))]
    else:
        if os.path.isdir(bpy.path.abspath(output_path)):
            filepath = os.path.join(bpy.path.abspath(output_path), f"{aob.name}_{axis}.svg")
        else:
            filepath = os.path.join(os.path.dirname(bpy.data.filepath), f"{aob.name}_{axis}.svg") if not output_path else bpy.path.abspath(output_path)

    svg_doc = ""
    x_offset, y_offset = 0, 0
    row_height = 0
    last_xmax = 0

    for i, vertices in enumerate(vertex_blocks):
        if sep_files or i == 0:
            svg_group = ""
        valid_verts = [v for v in vertices if isinstance(v, (list, tuple))]
        xmin = min(v[0] for v in valid_verts)
        xmax = max(v[0] for v in valid_verts)
        ymin = min(v[1] for v in valid_verts)
        ymax = max(v[1] for v in valid_verts)
        logger.debug(f"[DEBUG] Slice {i}: Polygon bounds: xmin={xmin}, xmax={xmax}, ymin={ymin}, ymax={ymax}")

        cxsize, cysize = xmax - xmin + cut_thickness, ymax - ymin + cut_thickness

        if (sep_files and svg_pos == '0') or (sep_files and i == 0 and svg_pos == '1'):
            x_offset, y_offset = -xmin + cut_thickness, -ymin + cut_thickness
        elif (sep_files and svg_pos == '1') or not sep_files:
            if scale_mm * (last_xmax + cxsize) <= mat_width:
                x_offset = last_xmax - xmin + cut_thickness
                y_offset = y_offset - ymin + cut_thickness
                row_height = max(row_height, cysize)
                last_xmax += cxsize
            else:
                y_offset += row_height
                x_offset = -xmin + cut_thickness
                last_xmax = cxsize
                row_height = cysize
        elif sep_files and svg_pos == '2':
            x_offset = mat_width / (2 * scale_mm) - 0.5 * cxsize - xmin
            y_offset = mat_height / (2 * scale_mm) - 0.5 * cysize - ymin

        logger.debug(f"[DEBUG] Slice {i}: x_offset={x_offset}, y_offset={y_offset}, row_height={row_height}, last_xmax={last_xmax}")
        svg_group += '<g>\n'

        if settings.laser_slicer_labels:
            slice_label = f"{axis}-{i:02d}"
            if valid_verts:
                label_x = scale * (x_offset + (xmin + xmax) / 2)
                label_y = scale * (y_offset + (ymin + ymax) / 2)
                svg_group += f'<text x="{label_x:.2f}" y="{label_y:.2f}" font-size="12" fill="black" text-anchor="middle">{slice_label}</text>\n'

        svg_points = []
        if not use_polygons:
            for edge in edge_blocks[i]:
                p1 = (scale * (x_offset + edge[0][0]), scale * (y_offset + edge[0][1]))
                p2 = (scale * (x_offset + edge[1][0]), scale * (y_offset + edge[1][1]))
                svg_points.append((p1, p2))
                svg_group += f'<line x1="{p1[0]}" y1="{p1[1]}" x2="{p2[0]}" y2="{p2[1]}" style="stroke:rgb({int(255*cut_color[0])},{int(255*cut_color[1])},{int(255*cut_color[2])});stroke-width:{line_thick}" />\n'
        else:
            points = ""
            for v in vertices:
                if v in (-1, -2):
                    # Add the points as tuples to the svg_points list
                    for point in points.split():
                        svg_points.append(tuple(map(float, point.split(','))))
                    tag = 'polygon' if v == -1 else 'polyline'
                    svg_group += f'<{tag} points="{points.strip()}" style="fill:none;stroke:rgb({int(255*cut_color[0])},{int(255*cut_color[1])},{int(255*cut_color[2])});stroke-width:{line_thick}" />\n'
                    points = ""
                else:
                    x, y = scale * (x_offset + v[0]), scale * (y_offset + v[1])
                    svg_points.append((x, y))
                    points += f"{x:.4f},{y:.4f} "
            if points:
                svg_group += f'<polygon points="{points.strip()}" style="fill:none;stroke:rgb({int(255*cut_color[0])},{int(255*cut_color[1])},{int(255*cut_color[2])});stroke-width:{line_thick}" />\n'
        logger.debug(f"[DEBUG] SVG points for slice {i}: {svg_points}")

        # Check for identical points and give only the identicals
        if len(svg_points) > 1:
            dupe_points = [p for p in svg_points if svg_points.count(p) > 1]
            logger.debug(f"[DEBUG] Duplicate points for slice {i}: {dupe_points}")

        if settings.laser_slicer_intersections:
            # Retrieve the slice object (which holds the polygon data).
            slice_obj = None
            for obj in bpy.context.scene.objects:
                if obj.get('Slices_' + axis):
                    slice_obj = obj
                    break
            if slice_obj is None:
                slice_obj = aob  # Fallback.

            logger.debug(f"[DEBUG] Processing intersections for slice {i} with slice_obj location: {slice_obj.location}")

            for other_axis, positions in all_slice_positions.items():
                if other_axis == axis:
                    continue  # Skip the current slicing axis.

                other_axis_index = {'X': 0, 'Y': 1, 'Z': 2}[other_axis]
                current_axis_index = {'X': 0, 'Y': 1, 'Z': 2}[axis]
                proj_axes = [i for i in range(3) if i != current_axis_index]

                if other_axis_index in proj_axes:
                    # For vertical lines, we use the x coordinate; for horizontal, the y coordinate.
                    dir_2d = proj_axes.index(other_axis_index)
                    logger.debug(f"[DEBUG] Other Axis: {other_axis} | dir_2d: {dir_2d}")
                    for pos in positions:
                        logger.debug(f"[DEBUG] Other Axis: {other_axis} | pos: {pos}")
                        local_pos = pos - slice_obj.location[other_axis_index]
                        if dir_2d == 0:
                            local_x_pos = scale * (x_offset + local_pos)
                            logger.debug(f"[DEBUG] Calculated local coordinate for {other_axis}: {local_x_pos}")
                            min_y = None  # To store the minimum y intersection
                            max_y = None  # To store the maximum y intersection

                            for p1, p2 in svg_points:
                                # Check if the vertical line at local_x_pos might intersect the segment.
                                if min(p1[0], p2[0]) <= local_x_pos <= max(p1[0], p2[0]):
                                    # Avoid division by zero: if the segment is vertical, you might skip or handle separately.
                                    if p1[0] == p2[0]:
                                        logger.debug(f"[DEBUG] Skipping vertical segment between {p1} and {p2}")
                                        continue  # Or handle this case as needed.

                                    # Calculate the slope.
                                    slope = (p2[1] - p1[1]) / (p2[0] - p1[0])
                                    # Calculate the y-value of the intersection.
                                    y_int = slope * (local_x_pos - p1[0]) + p1[1]
                                    logger.debug(f"[DEBUG] Intersection found at ({local_x_pos}, {y_int}) on segment between {p1} and {p2}")

                                    # Update the minimum and maximum y intersection values.
                                    if min_y is None or y_int < min_y:
                                        min_y = y_int
                                    if max_y is None or y_int > max_y:
                                        max_y = y_int

                            if min_y is None or max_y is None:
                                logger.debug(f"[DEBUG] No intersections found for slice {i} on axis {other_axis}")
                                continue
                            else:
                                # After the loop, min_y holds the lowest intersection y-value
                                # and max_y holds the highest intersection y-value.
                                logger.debug(f"[DEBUG] Final min_y: {min_y}, max_y: {max_y}")

                                # Draw min and max points
                                r = 3
                                svg_group += f'<circle cx="{local_x_pos:.2f}" cy="{min_y:.2f}" r="{r}" fill="red" />\n'
                                svg_group += f'<circle cx="{local_x_pos:.2f}" cy="{max_y:.2f}" r="{r}" fill="red" />\n'
                                # Draw the intersection line
                                svg_group += f'<line x1="{local_x_pos:.2f}" y1="{min_y:.2f}" x2="{local_x_pos:.2f}" y2="{max_y:.2f}" style="stroke:gray;stroke-width:0.5;stroke-dasharray:5,5" />\n'
                                if settings.laser_slicer_labels:
                                    label = f"{other_axis}-{positions.index(pos):02d}"
                                    svg_group += f'<text x="{local_x_pos:.2f}" y="{min_y-5:.2f}" font-size="8" fill="gray" text-anchor="middle">{label} start</text>\n'
                                    svg_group += f'<text x="{local_x_pos:.2f}" y="{max_y+12:.2f}" font-size="8" fill="gray" text-anchor="middle">{label} end</text>\n'
                            
                        else:
                            local_y_pos = scale * (y_offset + local_pos)
                            logger.debug(f"[DEBUG] Calculated local coordinate for {other_axis}: {local_y_pos}")
                            min_x = None  # To store the minimum x intersection
                            max_x = None  # To store the maximum x intersection

                            for p1, p2 in svg_points:
                                # Check if the horizontal line at local_y_pos might intersect the segment
                                if min(p1[1], p2[1]) <= local_y_pos <= max(p1[1], p2[1]):
                                    # Avoid division by zero: if the segment is horizontal, handle separately or skip.
                                    if p1[1] == p2[1]:
                                        # If the entire segment lies on local_y_pos, you may want to consider a different logic,
                                        # such as taking one of the endpoints, but here we skip it.
                                        logger.debug(f"[DEBUG] Skipping horizontal segment between {p1} and {p2}")
                                        continue

                                    # Calculate the slope of the segment.
                                    slope = (p2[1] - p1[1]) / (p2[0] - p1[0])
                                    # Compute the x coordinate of the intersection.
                                    x_int = p1[0] + (local_y_pos - p1[1]) / slope
                                    logger.debug(f"[DEBUG] Intersection found at ({x_int}, {local_y_pos}) on segment between {p1} and {p2}")

                                    # Update the minimum and maximum x intersection values.
                                    if min_x is None or x_int < min_x:
                                        min_x = x_int
                                    if max_x is None or x_int > max_x:
                                        max_x = x_int

                            if min_x is None or max_x is None:
                                logger.debug(f"[DEBUG] No intersections found for slice {i} on axis {other_axis}")
                                continue
                            else:
                                # After the loop, min_x holds the leftmost intersection x-value 
                                # and max_x holds the rightmost intersection x-value.
                                logger.debug(f"[DEBUG] Final min_x: {min_x}, max_x: {max_x}")

                                # Draw min and max points
                                r = 3
                                svg_group += f'<circle cx="{min_x:.2f}" cy="{local_y_pos:.2f}" r="{r}" fill="red" />\n'
                                svg_group += f'<circle cx="{max_x:.2f}" cy="{local_y_pos:.2f}" r="{r}" fill="red" />\n'
                                # Draw the intersection line
                                svg_group += f'<line x1="{min_x:.2f}" y1="{local_y_pos:.2f}" x2="{local_y_pos:.2f}" y2="{max_x:.2f}" style="stroke:gray;stroke-width:0.5;stroke-dasharray:5,5" />\n'
                                if settings.laser_slicer_labels:
                                    label = f"{other_axis}-{positions.index(pos):02d}"
                                    svg_group += f'<text x="{min_x-5:.2f}" y="{local_y_pos:.2f}" font-size="8" fill="gray" text-anchor="middle">{label} start</text>\n'
                                    svg_group += f'<text x="{max_x+12:.2f}" y="{local_y_pos:.2f}" font-size="8" fill="gray" text-anchor="middle">{label} end</text>\n'                 


        svg_group += '</g>\n'

        if sep_files:
            logger.debug(f"[DEBUG] Saving separate SVG file: {filepaths[i]}")
            with open(filepaths[i], 'w') as f:
                f.write('<?xml version="1.0"?>\n')
                f.write('<svg xmlns="http://www.w3.org/2000/svg" version="1.1" ')
                f.write(f'width="{mat_width*dpi_scale}" height="{mat_height*dpi_scale}" viewBox="0 0 {mat_width*dpi_scale} {mat_height*dpi_scale}">\n')
                f.write(svg_group)
                f.write('</svg>\n')
            logger.debug(f"[DEBUG] File saved: {filepaths[i]}")
        else:
            svg_doc += svg_group

    if not sep_files:
        logger.debug(f"[DEBUG] Saving combined SVG file: {filepath}")
        with open(filepath, 'w') as f:
            f.write('<?xml version="1.0"?>\n')
            f.write('<svg xmlns="http://www.w3.org/2000/svg" version="1.1" ')
            f.write(f'width="{mat_width*dpi_scale}" height="{mat_height*dpi_scale}" viewBox="0 0 {mat_width*dpi_scale} {mat_height*dpi_scale}">\n')
            f.write(svg_doc)
            f.write('</svg>\n')
        logger.debug(f"[DEBUG] File saved: {filepath}")

# === Main Slicer ===

def slicer(settings):
    scene = bpy.context.scene
    aob = bpy.context.active_object
    scale_mm = 1000 * scene.unit_settings.scale_length

    axes = []
    if settings.laser_slicer_axis_x:
        axes.append('X')
    if settings.laser_slicer_axis_y:
        axes.append('Y')
    if settings.laser_slicer_axis_z:
        axes.append('Z')

    all_slices = {}
    results_per_axis = {}

    for axis in axes:
        result = perform_axis_slicing(aob, axis, settings, scale_mm)
        results_per_axis[axis] = result
        all_slices[axis] = result['slice_positions']

    for axis in axes:
        generate_svg_files(results_per_axis[axis], aob, axis, settings, scale_mm, all_slices)

# Operator
class OBJECT_OT_Laser_Slicer(bpy.types.Operator):
    bl_label = "Laser Slicer"
    bl_idname = "object.laser_slicer"

    def execute(self, context):
        log_path = bpy.path.abspath(context.scene.slicer_settings.laser_slicer_log_path)
        log_handler = logging.FileHandler(log_path)

        try:
            logger.addHandler(log_handler)
            slicer(context.scene.slicer_settings)
            return {'FINISHED'}

        except Exception as e:
            logger.error(e)
            return {'CANCELLED'}

        finally:
            logger.removeHandler(log_handler)

# UI Panel
class OBJECT_PT_Laser_Slicer_Panel(bpy.types.Panel):
    bl_label = "Laser Slicer Panel"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Laser"
    bl_context = "objectmode"

    def draw(self, context):
        layout = self.layout
        settings = context.scene.slicer_settings

        box = layout.box()
        box.label(text="Material Dimensions:")
        newrow(box, "Thickness (mm):", settings, 'laser_slicer_material_thick')
        newrow(box, "Width (mm):", settings, 'laser_slicer_material_width')
        newrow(box, "Height (mm):", settings, 'laser_slicer_material_height')

        box = layout.box()
        box.label(text="Axis to Slice:")
        row = box.row()
        row.prop(settings, "laser_slicer_axis_x", toggle=True)
        row.prop(settings, "laser_slicer_axis_y", toggle=True)
        row.prop(settings, "laser_slicer_axis_z", toggle=True)
        newrow(box, "Show intersections:", settings, 'laser_slicer_intersections')

        box = layout.box()
        box.label(text="Cut Settings:")
        newrow(box, "DPI:", settings, 'laser_slicer_dpi')
        newrow(box, "Line Colour:", settings, 'laser_slicer_cut_colour')
        newrow(box, "Thickness (px):", settings, 'laser_slicer_cut_line')
        newrow(box, "Separate Files:", settings, 'laser_slicer_separate_files')
        if settings.laser_slicer_separate_files:
            newrow(box, "Cut Position:", settings, 'laser_slicer_svg_position')
        newrow(box, "Cut spacing (mm):", settings, 'laser_slicer_cut_spacing')
        newrow(box, "Laser kerf (mm):", settings, 'laser_slicer_cut_thickness')
        newrow(box, "SVG Polygons:", settings, 'laser_slicer_accuracy')
        newrow(box, "Export Path:", settings, 'laser_slicer_ofile')
        newrow(box, "Enable Labels:", settings, 'laser_slicer_labels')

        
        box = layout.box()
        box.label(text="Log Settings:")
        newrow(box, "Log Path:", settings, 'laser_slicer_log_path')

        if context.active_object and context.active_object.select_get():
            aob = context.active_object
            dp = context.evaluated_depsgraph_get()
            temp_mesh = aob.evaluated_get(dp).to_mesh()
            bm = bmesh.new()
            bm.from_mesh(temp_mesh)
            bm.transform(aob.matrix_world)
            aob.evaluated_get(dp).to_mesh_clear()

            axes = []
            if settings.laser_slicer_axis_x:
                axes.append('X')
            if settings.laser_slicer_axis_y:
                axes.append('Y')
            if settings.laser_slicer_axis_z:
                axes.append('Z')

            slices = 0
            for axis in axes:
                axis_index = {'X': 0, 'Y': 1, 'Z': 2}[axis]
                coords = [v.co[axis_index] for v in bm.verts]
                minv = min(coords)
                maxv = max(coords)

                scale_mm = 1000 * context.scene.unit_settings.scale_length
                mat_thickness = settings.laser_slicer_material_thick / scale_mm
                spacing = settings.laser_slicer_cut_spacing / scale_mm

                total_length = maxv - minv
                slices += 1 + math.ceil(total_length / (mat_thickness + spacing))

            bm.free()

            row = layout.row()
            row.label(text=f"No. of slices: {slices}")
            if bpy.data.filepath or settings.laser_slicer_ofile:
                layout.operator("object.laser_slicer", text="Slice the object")

# Properties
class Slicer_Settings(bpy.types.PropertyGroup):
    laser_slicer_material_thick: FloatProperty(name="", description="Material thickness (mm)", min=0.1, max=50, default=2)
    laser_slicer_material_width: FloatProperty(name="", description="Material width (mm)", min=1, max=5000, default=450)
    laser_slicer_material_height: FloatProperty(name="", description="Material height (mm)", min=1, max=5000, default=450)
    laser_slicer_dpi: IntProperty(name="", description="DPI", min=50, max=500, default=96)
    laser_slicer_separate_files: BoolProperty(name="", description="Separate SVG files", default=False)
    laser_slicer_svg_position: EnumProperty(
        items=[('0', 'Top left', ''), ('1', 'Staggered', ''), ('2', 'Centre', '')],
        name="", description="SVG layout", default='0')
    laser_slicer_cut_spacing: FloatProperty(name="", description="Cut spacing (mm)", min=0, max=500, default=1)
    laser_slicer_cut_thickness: FloatProperty(name="", description="Laser kerf (mm)", min=0, max=5, default=1)
    laser_slicer_ofile: StringProperty(name="", description="Output file", subtype="FILE_PATH", default="")
    laser_slicer_accuracy: BoolProperty(name="", description="Use accurate polygons", default=False)
    laser_slicer_cut_colour: FloatVectorProperty(size=3, name="", subtype='COLOR', min=0, max=1, default=(1.0, 0.0, 0.0))
    laser_slicer_cut_line: FloatProperty(name="", description="SVG line width (px)", min=0, max=5, default=1)
    laser_slicer_axis_x: BoolProperty(name="X", description="Slice along X axis", default=False)
    laser_slicer_axis_y: BoolProperty(name="Y", description="Slice along Y axis", default=False)
    laser_slicer_axis_z: BoolProperty(name="Z", description="Slice along Z axis", default=True)
    laser_slicer_labels: BoolProperty(name="", description="Enable labels in svg", default=True)
    laser_slicer_intersections: BoolProperty(name="", description="Show intersections in svg", default=True)
    laser_slicer_log_path: StringProperty(name="", description="Log file path", subtype="FILE_PATH", default="")

# Registering
classes = (OBJECT_PT_Laser_Slicer_Panel, OBJECT_OT_Laser_Slicer, Slicer_Settings)

def register():
    for cls in classes:
        bpy.utils.register_class(cls)
    bpy.types.Scene.slicer_settings = bpy.props.PointerProperty(type=Slicer_Settings)

def unregister():
    for cls in classes:
        bpy.utils.unregister_class(cls)
    del bpy.types.Scene.slicer_settings
