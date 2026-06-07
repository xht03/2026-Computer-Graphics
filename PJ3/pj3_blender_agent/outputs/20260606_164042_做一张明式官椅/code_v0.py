import bpy
import math

# Clear scene
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

# Create wood material
wood_mat = bpy.data.materials.new(name="Wood_Material")
wood_mat.use_nodes = True
nodes = wood_mat.node_tree.nodes
nodes.clear()
bsdf = nodes.new('ShaderNodeBsdfPrincipled')
bsdf.inputs['Base Color'].default_value = (0.4, 0.25, 0.15, 1.0)
bsdf.inputs['Roughness'].default_value = 0.6
output = nodes.new('ShaderNodeOutputMaterial')
wood_mat.node_tree.links.new(bsdf.outputs['BSDF'], output.inputs['Surface'])

def add_mesh(name, verts, faces, location, size):
    mesh = bpy.data.meshes.new(name)
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.collection.objects.link(obj)
# [sanitized]     mesh.from_pydata(verts, [], faces)
    obj.location = location
    obj.scale = size
    bpy.ops.object.select_all(action='DESELECT')
    obj.select_set(True)
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY', center='BOUNDS')
    obj.data.materials.append(wood_mat)
    return obj

# COMPONENT: leg
leg_positions = [
    (-0.225, -0.19, 0.24),
    (0.225, -0.19, 0.24),
    (-0.225, 0.19, 0.24),
    (0.225, 0.19, 0.24)
]
for i, pos in enumerate(leg_positions):
    bpy.ops.mesh.primitive_cube_add(size=1, location=pos)
    leg = bpy.context.active_object
    leg.name = f"leg_{i+1}"
    leg.scale = (0.1, 0.1, 0.48)
    leg.data.materials.append(wood_mat)

# COMPONENT: seat
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, 0.705))
seat = bpy.context.active_object
seat.name = "seat"
seat.scale = (0.55, 0.48, 0.45)
seat.data.materials.append(wood_mat)

# COMPONENT: back_panel
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, -0.29, 0.6))
back_panel = bpy.context.active_object
back_panel.name = "back_panel"
back_panel.scale = (0.55, 0.1, 1.1)
back_panel.data.materials.append(wood_mat)

# COMPONENT: armrest
armrest_positions = [
    (-0.225, -0.05, 0.9),
    (0.225, -0.05, 0.9)
]
for i, pos in enumerate(armrest_positions):
    bpy.ops.mesh.primitive_cube_add(size=1, location=pos)
    armrest = bpy.context.active_object
    armrest.name = f"armrest_{i+1}"
    armrest.scale = (0.1, 0.48, 0.1)
    armrest.data.materials.append(wood_mat)

# COMPONENT: top_rail
# Create curved top rail using bezier curve
curve_data = bpy.data.curves.new('top_rail_curve', type='CURVE')
curve_data.dimensions = '3D'
spline = curve_data.splines.new('BEZIER')
spline.resolution_u = 32
spline.bezier_points.add(2)

# Control points for gentle forward curve
points = spline.bezier_points
points[0].co = (-0.275, -0.29, 1.15)
points[0].handle_right = (-0.14, -0.29, 1.15)
points[0].handle_left = (-0.275, -0.29, 1.15)

points[1].co = (0, -0.18, 1.15)
points[1].handle_left = (-0.14, -0.2, 1.15)
points[1].handle_right = (0.14, -0.2, 1.15)

points[2].co = (0.275, -0.29, 1.15)
points[2].handle_left = (0.14, -0.29, 1.15)
points[2].handle_right = (0.275, -0.29, 1.15)

# Set bevel for thickness
curve_data.bevel_depth = 0.05
curve_data.bevel_resolution = 4
curve_data.fill_mode = 'FULL'

top_rail = bpy.data.objects.new("top_rail", curve_data)
bpy.context.collection.objects.link(top_rail)
top_rail.data.materials.append(wood_mat)