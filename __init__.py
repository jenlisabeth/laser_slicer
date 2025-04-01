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
    "version": (0, 9, 2),
    "blender": (2, 80, 0),
    "location": "3D View > Tools Panel",
    "description": "Creates cross-sections and exports as SVG for laser cutting",
    "warning": "",
    "wiki_url": "tba",
    "tracker_url": "https://github.com/rgsouthall/laser_slicer/issues",
    "category": "Object"
}

import bpy, os, bmesh, numpy as np
from bpy.props import FloatProperty, BoolProperty, EnumProperty, IntProperty, StringProperty, FloatVectorProperty

# UI row helper
def newrow(layout, label, obj, prop_name):
    row = layout.row()
    row.label(text=label)
    row.prop(obj, prop_name)

# Main slicing logic
def slicer(settings):
    scene = bpy.context.scene
    dp = bpy.context.evaluated_depsgraph_get()
    aob = bpy.context.active_object
    scale_mm = 1000 * scene.unit_settings.scale_length
    temp_mesh = aob.evaluated_get(dp).to_mesh()
    bm = bmesh.new()
    bm.from_mesh(temp_mesh)
    bm.transform(aob.matrix_world)
    aob.evaluated_get(dp).to_mesh_clear()
    aob.select_set(False)

    # Settings
    mat_thickness = settings.laser_slicer_material_thick / scale_mm
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

    # Get slice range
    minz = min(v.co.z for v in bm.verts)
    maxz = max(v.co.z for v in bm.verts)
    slice_z = minz + mat_thickness * 0.5

    # Setup variables
    dpi_scale = dpi / 25.4
    scale = scale_mm * dpi_scale
    slice_vpos = []
    vertex_blocks = []
    edge_blocks = []
    svg_elements = []
    v_total, e_total = 0, 0
    vcount_list, ecount_list = [0], [0]
    vertex_index = np.empty(0, dtype=np.int8)

    # Output mesh object
    cob = None
    cob_exists = False
    for obj in bpy.context.scene.objects:
        if obj.get('Slices'):
            cob = obj
            cob_exists = True
            me = cob.data
            bpy.context.view_layer.objects.active = cob
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.mesh.delete(type='VERT')
            bpy.ops.object.mode_set(mode='OBJECT')
            break
    if not cob_exists:
        me = bpy.data.meshes.new('Slices')
        cob = bpy.data.objects.new('Slices', me)
        cob['Slices'] = 1

    # Slice loop
    while slice_z < maxz:
        cbm = bm.copy()
        cut = bmesh.ops.bisect_plane(
            cbm, 
            geom=cbm.edges[:] + cbm.faces[:], 
            dist=0, 
            plane_co=(0, 0, slice_z), 
            plane_no=(0, 0, 1),
            clear_outer=False, 
            clear_inner=False
        )['geom_cut']

        verts = [v for v in cut if isinstance(v, bmesh.types.BMVert)]
        edges = [e for e in cut if isinstance(e, bmesh.types.BMEdge)]
        voffset = min(v.index for v in verts)

        v_block = np.array([v.co for v in verts]).flatten()
        e_index = np.array([[v.index - voffset + v_total for v in e.verts] for e in edges]).flatten()

        slice_vpos = np.append(slice_vpos, v_block)
        vertex_blocks.append([(v.co - cob.location)[0:2] for v in verts])
        edge_blocks.append([[(v.co - cob.location)[0:2] for v in e.verts] for e in edges])
        vertex_index = np.append(vertex_index, e_index)

        v_total += len(verts)
        e_total += len(edges)
        vcount_list.append(v_total)
        ecount_list.append(e_total)
        slice_z += mat_thickness
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

            polygon_paths.append([(me.vertices[v].co[:2] if v >= 0 else v) for v in path])
        vertex_blocks = polygon_paths

    # Setup SVG filenames
    if sep_files:
        if os.path.isdir(bpy.path.abspath(output_path)):
            filepaths = [os.path.join(bpy.path.abspath(output_path), f"{aob.name}{i}.svg") for i in range(len(vertex_blocks))]
        else:
            out_dir = os.path.dirname(bpy.path.abspath(output_path or bpy.data.filepath))
            name_base = bpy.path.display_name_from_filepath(output_path or aob.name)
            filepaths = [os.path.join(out_dir, f"{name_base}{i}.svg") for i in range(len(vertex_blocks))]
    else:
        if os.path.isdir(bpy.path.abspath(output_path)):
            filepath = os.path.join(bpy.path.abspath(output_path), f"{aob.name}.svg")
        else:
            filepath = os.path.join(os.path.dirname(bpy.data.filepath), f"{aob.name}.svg") if not output_path else bpy.path.abspath(output_path)

    # Generate SVG
    svg_doc = ""
    x_offset, y_offset = 0, 0
    row_height = 0
    last_xmax = 0

    for i, vertices in enumerate(vertex_blocks):
        if sep_files or i == 0:
            svg_group = ""

        valid_verts = [v for v in vertices if isinstance(v, (list, tuple))]
        xmin, xmax = min(v[0] for v in valid_verts), max(v[0] for v in valid_verts)
        ymin, ymax = min(v[1] for v in valid_verts), max(v[1] for v in valid_verts)
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

        if not use_polygons:
            svg_group += '<g>\n'
            for edge in edge_blocks[i]:
                p1 = (scale * (x_offset + edge[0][0]), scale * (y_offset + edge[0][1]))
                p2 = (scale * (x_offset + edge[1][0]), scale * (y_offset + edge[1][1]))
                svg_group += f'<line x1="{p1[0]}" y1="{p1[1]}" x2="{p2[0]}" y2="{p2[1]}" style="stroke:rgb({int(255*cut_color[0])},{int(255*cut_color[1])},{int(255*cut_color[2])});stroke-width:{line_thick}" />\n'
            svg_group += '</g>\n'
        else:
            svg_group += '<g>\n'
            points = ""
            for v in vertices:
                if v in (-1, -2):
                    tag = 'polygon' if v == -1 else 'polyline'
                    svg_group += f'<{tag} points="{points.strip()}" style="fill:none;stroke:rgb({int(255*cut_color[0])},{int(255*cut_color[1])},{int(255*cut_color[2])});stroke-width:{line_thick}" />\n'
                    points = ""
                else:
                    x, y = scale * (x_offset + v[0]), scale * (y_offset + v[1])
                    points += f"{x:.4f},{y:.4f} "
            if points:
                svg_group += f'<polygon points="{points.strip()}" style="fill:none;stroke:rgb({int(255*cut_color[0])},{int(255*cut_color[1])},{int(255*cut_color[2])});stroke-width:{line_thick}" />\n'
            svg_group += '</g>\n'

        if sep_files:
            with open(filepaths[i], 'w') as f:
                f.write('<?xml version="1.0"?>\n')
                f.write('<svg xmlns="http://www.w3.org/2000/svg" version="1.1" ')
                f.write(f'width="{mat_width*dpi_scale}" height="{mat_height*dpi_scale}" viewBox="0 0 {mat_width*dpi_scale} {mat_height*dpi_scale}">\n')
                f.write(svg_group)
                f.write('</svg>\n')
        else:
            svg_doc += svg_group

    if not sep_files:
        with open(filepath, 'w') as f:
            f.write('<?xml version="1.0"?>\n')
            f.write('<svg xmlns="http://www.w3.org/2000/svg" version="1.1" ')
            f.write(f'width="{mat_width*dpi_scale}" height="{mat_height*dpi_scale}" viewBox="0 0 {mat_width*dpi_scale} {mat_height*dpi_scale}">\n')
            f.write(svg_doc)
            f.write('</svg>\n')

    if not cob_exists:
        bpy.context.scene.collection.objects.link(cob)

    bpy.context.view_layer.objects.active = cob
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.object.mode_set(mode='OBJECT')
    aob.select_set(True)
    bpy.context.view_layer.objects.active = aob

# Operator
class OBJECT_OT_Laser_Slicer(bpy.types.Operator):
    bl_label = "Laser Slicer"
    bl_idname = "object.laser_slicer"

    def execute(self, context):
        slicer(context.scene.slicer_settings)
        return {'FINISHED'}

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

        layout.label(text="Material Dimensions:")
        newrow(layout, "Thickness (mm):", settings, 'laser_slicer_material_thick')
        newrow(layout, "Width (mm):", settings, 'laser_slicer_material_width')
        newrow(layout, "Height (mm):", settings, 'laser_slicer_material_height')

        layout.label(text="Cut Settings:")
        newrow(layout, "DPI:", settings, 'laser_slicer_dpi')
        newrow(layout, "Line Colour:", settings, 'laser_slicer_cut_colour')
        newrow(layout, "Thickness (px):", settings, 'laser_slicer_cut_line')
        newrow(layout, "Separate Files:", settings, 'laser_slicer_separate_files')

        if settings.laser_slicer_separate_files:
            newrow(layout, "Cut Position:", settings, 'laser_slicer_svg_position')

        newrow(layout, "Cut Spacing (mm):", settings, 'laser_slicer_cut_thickness')
        newrow(layout, "SVG Polygons:", settings, 'laser_slicer_accuracy')
        newrow(layout, "Export Path:", settings, 'laser_slicer_ofile')

        if context.active_object and context.active_object.select_get():
            row = layout.row()
            slices = context.active_object.dimensions[2] * 1000 * context.scene.unit_settings.scale_length / settings.laser_slicer_material_thick
            row.label(text=f"No. of slices: {slices:.0f}")
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
    laser_slicer_cut_thickness: FloatProperty(name="", description="Laser kerf (mm)", min=0, max=5, default=1)
    laser_slicer_ofile: StringProperty(name="", description="Output file", subtype="FILE_PATH", default="")
    laser_slicer_accuracy: BoolProperty(name="", description="Use accurate polygons", default=False)
    laser_slicer_cut_colour: FloatVectorProperty(size=3, name="", subtype='COLOR', min=0, max=1, default=(1.0, 0.0, 0.0))
    laser_slicer_cut_line: FloatProperty(name="", description="SVG line width (px)", min=0, max=5, default=1)

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
