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

# COMPONENT: leg (continuous posts)
# Back legs extending to top rail
bpy.ops.mesh.primitive_cube_add(size=1, location=(-0.225, -0.19, 0.575))
leg_1 = bpy.context.active_object
leg_1.name = "leg_1"
leg_1.scale = (0.08, 0.08, 1.15)
leg_1.data.materials.append(wood_mat)

bpy.ops.mesh.primitive_cube_add(size=1, location=(0.225, -0.19, 0.575))
leg_2 = bpy.context.active_object
leg_2.name = "leg_2"
leg_2.scale = (0.08, 0.08, 1.15)
leg_2.data.materials.append(wood_mat)

# Front legs extending to armrest height
bpy.ops.mesh.primitive_cube_add(size=1, location=(-0.225, 0.19, 0.45))
leg_3 = bpy.context.active_object
leg_3.name = "leg_3"
leg_3.scale = (0.08, 0.08, 0.9)
leg_3.data.materials.append(wood_mat)

bpy.ops.mesh.primitive_cube_add(size=1, location=(0.225, 0.19, 0.45))
leg_4 = bpy.context.active_object
leg_4.name = "leg_4"
leg_4.scale = (0.08, 0.08, 0.9)
leg_4.data.materials.append(wood_mat)

# COMPONENT: stretchers (horizontal supports between legs)
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0.19, 0.25))
front_stretcher = bpy.context.active_object
front_stretcher.name = "front_stretcher"
front_stretcher.scale = (0.34, 0.06, 0.06)
front_stretcher.data.materials.append(wood_mat)

bpy.ops.mesh.primitive_cube_add(size=1, location=(0, -0.19, 0.25))
back_stretcher = bpy.context.active_object
back_stretcher.name = "back_stretcher"
back_stretcher.scale = (0.34, 0.06, 0.06)
back_stretcher.data.materials.append(wood_mat)

bpy.ops.mesh.primitive_cube_add(size=1, location=(-0.225, 0, 0.25))
left_stretcher = bpy.context.active_object
left_stretcher.name = "left_stretcher"
left_stretcher.scale = (0.06, 0.38, 0.06)
left_stretcher.data.materials.append(wood_mat)

bpy.ops.mesh.primitive_cube_add(size=1, location=(0.225, 0, 0.25))
right_stretcher = bpy.context.active_object
right_stretcher.name = "right_stretcher"
right_stretcher.scale = (0.06, 0.38, 0.06)
right_stretcher.data.materials.append(wood_mat)

# COMPONENT: seat (thinner panel)
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, 0.55))
seat = bpy.context.active_object
seat.name = "seat"
seat.scale = (0.5, 0.4, 0.08)
seat.data.materials.append(wood_mat)

# COMPONENT: back_panel (proper proportions between seat and top rail)
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, -0.19, 0.87))
back_panel = bpy.context.active_object
back_panel.name = "back_panel"
back_panel.scale = (0.4, 0.05, 0.28)
back_panel.data.materials.append(wood_mat)

# COMPONENT: armrest (curved bezier curves)
# Left armrest curving outward and downward
curve_data_l = bpy.data.curves.new('armrest_left_curve', type='CURVE')
curve_data_l.dimensions = '3D'
spline_l = curve_data_l.splines.new('BEZIER')
spline_l.resolution_u = 32
spline_l.bezier_points.add(2)

points_l = spline_l.bezier_points
points_l[0].co = (-0.225, -0.19, 0.9)
points_l[0].handle_right = (-0.26, -0.1, 0.9)
points_l[0].handle_left = (-0.19, -0.28, 0.9)

points_l[1].co = (-0.3, 0, 0.875)
points_l[1].handle_left = (-0.3, -0.1, 0.875)
points_l[1].handle_right = (-0.3, 0.1, 0.875)

points_l[2].co = (-0.225, 0.19, 0.85)
points_l[2].handle_left = (-0.26, 0.1, 0.85)
points_l[2].handle_right = (-0.19, 0.28, 0.85)

curve_data_l.bevel_depth = 0.04
curve_data_l.bevel_resolution = 4
curve_data_l.fill_mode = 'FULL'

armrest_left = bpy.data.objects.new("armrest_left", curve_data_l)
bpy.context.collection.objects.link(armrest_left)
armrest_left.data.materials.append(wood_mat)

# Right armrest (mirrored)
curve_data_r = bpy.data.curves.new('armrest_right_curve', type='CURVE')
curve_data_r.dimensions = '3D'
spline_r = curve_data_r.splines.new('BEZIER')
spline_r.resolution_u = 32
spline_r.bezier_points.add(2)

points_r = spline_r.bezier_points
points_r[0].co = (0.225, -0.19, 0.9)
points_r[0].handle_right = (0.19, -0.28, 0.9)
points_r[0].handle_left = (0.26, -0.1, 0.9)

points_r[1].co = (0.3, 0, 0.875)
points_r[1].handle_left = (0.3, -0.1, 0.875)
points_r[1].handle_right = (0.3, 0.1, 0.875)

points_r[2].co = (0.225, 0.19, 0.85)
points_r[2].handle_left = (0.19, 0.28, 0.85)
points_r[2].handle_right = (0.26, 0.1, 0.85)

curve_data_r.bevel_depth = 0.04
curve_data_r.bevel_resolution = 4
curve_data_r.fill_mode = 'FULL'

armrest_right = bpy.data.objects.new("armrest_right", curve_data_r)
bpy.context.collection.objects.link(armrest_right)
armrest_right.data.materials.append(wood_mat)

# COMPONENT: top_rail (extended horizontally beyond back posts)
curve_data = bpy.data.curves.new('top_rail_curve', type='CURVE')
curve_data.dimensions = '3D'
spline = curve_data.splines.new('BEZIER')
spline.resolution_u = 32
spline.bezier_points.add(2)

points = spline.bezier_points
points[0].co = (-0.4, -0.29, 1.15)
points[0].handle_right = (-0.2, -0.29, 1.15)
points[0].handle_left = (-0.4, -0.29, 1.15)

points[1].co = (0, -0.18, 1.15)
points[1].handle_left = (-0.2, -0.2, 1.15)
points[1].handle_right = (0.2, -0.2, 1.15)

points[2].co = (0.4, -0.29, 1.15)
points[2].handle_left = (0.2, -0.29, 1.15)
points[2].handle_right = (0.4, -0.29, 1.15)

curve_data.bevel_depth = 0.05
curve_data.bevel_resolution = 4
curve_data.fill_mode = 'FULL'

top_rail = bpy.data.objects.new("top_rail", curve_data)
bpy.context.collection.objects.link(top_rail)
top_rail.data.materials.append(wood_mat)