import bpy
import math

bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

wood = bpy.data.materials.new(name="Huanghuali")
wood.use_nodes = True
bsdf = wood.node_tree.nodes["Principled BSDF"]
bsdf.inputs["Base Color"].default_value = (0.45, 0.25, 0.1, 1.0)
bsdf.inputs["Roughness"].default_value = 0.6

def add_bevel(obj, width=0.004):
    mod = obj.modifiers.new("Bevel", type='BEVEL')
    mod.width = width
    mod.segments = 2

# COMPONENT: front_left_leg
bpy.ops.mesh.primitive_cube_add(size=1, location=(-0.255, 0.20, 0.26))
front_left_leg = bpy.context.active_object
front_left_leg.name = "front_left_leg"
front_left_leg.scale = (0.04, 0.04, 0.52)
add_bevel(front_left_leg)
front_left_leg.data.materials.append(wood)

# COMPONENT: front_right_leg
bpy.ops.mesh.primitive_cube_add(size=1, location=(0.255, 0.20, 0.26))
front_right_leg = bpy.context.active_object
front_right_leg.name = "front_right_leg"
front_right_leg.scale = (0.04, 0.04, 0.52)
add_bevel(front_right_leg)
front_right_leg.data.materials.append(wood)

# COMPONENT: back_left_leg
bpy.ops.mesh.primitive_cube_add(size=1, location=(-0.255, -0.20, 0.575))
back_left_leg = bpy.context.active_object
back_left_leg.name = "back_left_leg"
back_left_leg.scale = (0.04, 0.04, 1.15)
add_bevel(back_left_leg)
back_left_leg.data.materials.append(wood)

# COMPONENT: back_right_leg
bpy.ops.mesh.primitive_cube_add(size=1, location=(0.255, -0.20, 0.575))
back_right_leg = bpy.context.active_object
back_right_leg.name = "back_right_leg"
back_right_leg.scale = (0.04, 0.04, 1.15)
add_bevel(back_right_leg)
back_right_leg.data.materials.append(wood)

# COMPONENT: front_seat_rail
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0.20, 0.52))
front_seat_rail = bpy.context.active_object
front_seat_rail.name = "front_seat_rail"
front_seat_rail.scale = (0.55, 0.04, 0.04)
add_bevel(front_seat_rail)
front_seat_rail.data.materials.append(wood)

# COMPONENT: back_seat_rail
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, -0.20, 0.52))
back_seat_rail = bpy.context.active_object
back_seat_rail.name = "back_seat_rail"
back_seat_rail.scale = (0.55, 0.04, 0.04)
add_bevel(back_seat_rail)
back_seat_rail.data.materials.append(wood)

# COMPONENT: side_seat_rail_left
bpy.ops.mesh.primitive_cube_add(size=1, location=(-0.255, 0, 0.52))
side_seat_rail_left = bpy.context.active_object
side_seat_rail_left.name = "side_seat_rail_left"
side_seat_rail_left.scale = (0.04, 0.40, 0.04)
add_bevel(side_seat_rail_left)
side_seat_rail_left.data.materials.append(wood)

# COMPONENT: side_seat_rail_right
bpy.ops.mesh.primitive_cube_add(size=1, location=(0.255, 0, 0.52))
side_seat_rail_right = bpy.context.active_object
side_seat_rail_right.name = "side_seat_rail_right"
side_seat_rail_right.scale = (0.04, 0.40, 0.04)
add_bevel(side_seat_rail_right)
side_seat_rail_right.data.materials.append(wood)

# COMPONENT: seat_panel
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, 0.53))
seat_panel = bpy.context.active_object
seat_panel.name = "seat_panel"
seat_panel.scale = (0.47, 0.40, 0.02)
add_bevel(seat_panel, width=0.002)
seat_panel.data.materials.append(wood)

# COMPONENT: front_stretcher
bpy.ops.mesh.primitive_cylinder_add(radius=0.015, depth=0.47, location=(0, 0.20, 0.15))
front_stretcher = bpy.context.active_object
front_stretcher.name = "front_stretcher"
front_stretcher.rotation_euler = (0, math.pi/2, 0)
front_stretcher.data.materials.append(wood)

# COMPONENT: back_stretcher
bpy.ops.mesh.primitive_cylinder_add(radius=0.015, depth=0.47, location=(0, -0.20, 0.25))
back_stretcher = bpy.context.active_object
back_stretcher.name = "back_stretcher"
back_stretcher.rotation_euler = (0, math.pi/2, 0)
back_stretcher.data.materials.append(wood)

# COMPONENT: side_stretcher_left
bpy.ops.mesh.primitive_cylinder_add(radius=0.015, depth=0.48, location=(-0.255, 0, 0.20))
side_stretcher_left = bpy.context.active_object
side_stretcher_left.name = "side_stretcher_left"
side_stretcher_left.rotation_euler = (math.pi/2, 0, 0)
side_stretcher_left.data.materials.append(wood)

# COMPONENT: side_stretcher_right
bpy.ops.mesh.primitive_cylinder_add(radius=0.015, depth=0.48, location=(0.255, 0, 0.20))
side_stretcher_right = bpy.context.active_object
side_stretcher_right.name = "side_stretcher_right"
side_stretcher_right.rotation_euler = (math.pi/2, 0, 0)
side_stretcher_right.data.materials.append(wood)

# COMPONENT: back_splat
cu_splat = bpy.data.curves.new("back_splat_curve", type='CURVE')
cu_splat.dimensions = '3D'
cu_splat.bevel_depth = 0.02
sp_splat = cu_splat.splines.new('BEZIER')
sp_splat.bezier_points.add(2)
pts_splat = [(-0.06, -0.20, 0.52), (0, 0.0, 0.815), (0.06, -0.20, 1.07)]
for bp, co in zip(sp_splat.bezier_points, pts_splat):
    bp.co = co
    bp.handle_left_type = 'AUTO'
    bp.handle_right_type = 'AUTO'
back_splat = bpy.data.objects.new("back_splat", cu_splat)
bpy.context.collection.objects.link(back_splat)
back_splat.data.materials.append(wood)

# COMPONENT: top_rail
cu_top = bpy.data.curves.new("top_rail_curve", type='CURVE')
cu_top.dimensions = '3D'
cu_top.bevel_depth = 0.04
sp_top = cu_top.splines.new('BEZIER')
sp_top.bezier_points.add(2)
pts_top = [(-0.30, -0.20, 1.10), (0, -0.20, 1.15), (0.30, -0.20, 1.10)]
for bp, co in zip(sp_top.bezier_points, pts_top):
    bp.co = co
    bp.handle_left_type = 'AUTO'
    bp.handle_right_type = 'AUTO'
top_rail = bpy.data.objects.new("top_rail", cu_top)
bpy.context.collection.objects.link(top_rail)
top_rail.data.materials.append(wood)

# COMPONENT: armrest_left
cu_al = bpy.data.curves.new("armrest_left_curve", type='CURVE')
cu_al.dimensions = '3D'
cu_al.bevel_depth = 0.02
sp_al = cu_al.splines.new('BEZIER')
sp_al.bezier_points.add(1)
pts_al = [(-0.255, 0.25, 0.52), (-0.255, -0.25, 1.15)]
for bp, co in zip(sp_al.bezier_points, pts_al):
    bp.co = co
    bp.handle_left_type = 'AUTO'
    bp.handle_right_type = 'AUTO'
armrest_left = bpy.data.objects.new("armrest_left", cu_al)
bpy.context.collection.objects.link(armrest_left)
armrest_left.data.materials.append(wood)

# COMPONENT: armrest_right
cu_ar = bpy.data.curves.new("armrest_right_curve", type='CURVE')
cu_ar.dimensions = '3D'
cu_ar.bevel_depth = 0.02
sp_ar = cu_ar.splines.new('BEZIER')
sp_ar.bezier_points.add(1)
pts_ar = [(0.255, 0.25, 0.52), (0.255, -0.25, 1.15)]
for bp, co in zip(sp_ar.bezier_points, pts_ar):
    bp.co = co
    bp.handle_left_type = 'AUTO'
    bp.handle_right_type = 'AUTO'
armrest_right = bpy.data.objects.new("armrest_right", cu_ar)
bpy.context.collection.objects.link(armrest_right)
armrest_right.data.materials.append(wood)