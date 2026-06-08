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

# COMPONENT: back_splat
cu_splat = bpy.data.curves.new("back_splat_curve", type='CURVE')
cu_splat.dimensions = '3D'
cu_splat.extrude = 0.0075
sp_splat = cu_splat.splines.new('BEZIER')
sp_splat.bezier_points.add(3)
pts_splat = [(-0.235, -0.20, 0.545), (0.235, -0.20, 0.545), (0.235, -0.20, 1.15), (-0.235, -0.20, 1.15)]
for bp, co in zip(sp_splat.bezier_points, pts_splat):
    bp.co = co
    bp.handle_left_type = 'VECTOR'
    bp.handle_right_type = 'VECTOR'
sp_splat.use_cyclic_u = True
back_splat = bpy.data.objects.new("back_splat", cu_splat)
bpy.context.collection.objects.link(back_splat)
back_splat.data.materials.append(wood)

# COMPONENT: front_apron
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0.21, 0.47))
front_apron = bpy.context.active_object
front_apron.name = "front_apron"
front_apron.scale = (0.47, 0.02, 0.10)
add_bevel(front_apron)
front_apron.data.materials.append(wood)

# COMPONENT: top_rail
cu_top = bpy.data.curves.new("top_rail_curve", type='CURVE')
cu_top.dimensions = '3D'
cu_top.bevel_depth = 0.04
sp_top = cu_top.splines.new('BEZIER')
sp_top.bezier_points.add(2)
pts_top = [(-0.30, -0.25, 1.10), (0, -0.20, 1.15), (0.30, -0.25, 1.10)]
for bp, co in zip(sp_top.bezier_points, pts_top):
    bp.co = co
    bp.handle_left_type = 'AUTO'
    bp.handle_right_type = 'AUTO'
top_rail = bpy.data.objects.new("top_rail", cu_top)
bpy.context.collection.objects.link(top_rail)
top_rail.data.materials.append(wood)

# COMPONENT: left_seat_rail
bpy.ops.mesh.primitive_cube_add(size=1, location=(-0.255, 0, 0.52))
left_seat_rail = bpy.context.active_object
left_seat_rail.name = "left_seat_rail"
left_seat_rail.scale = (0.04, 0.40, 0.04)
add_bevel(left_seat_rail)
left_seat_rail.data.materials.append(wood)

# COMPONENT: right_seat_rail
bpy.ops.mesh.primitive_cube_add(size=1, location=(0.255, 0, 0.52))
right_seat_rail = bpy.context.active_object
right_seat_rail.name = "right_seat_rail"
right_seat_rail.scale = (0.04, 0.40, 0.04)
add_bevel(right_seat_rail)
right_seat_rail.data.materials.append(wood)

# COMPONENT: left_gooseneck
bpy.ops.mesh.primitive_cube_add(size=1, location=(-0.255, 0.20, 0.62))
left_gooseneck = bpy.context.active_object
left_gooseneck.name = "left_gooseneck"
left_gooseneck.scale = (0.04, 0.04, 0.20)
add_bevel(left_gooseneck)
left_gooseneck.data.materials.append(wood)

# COMPONENT: right_gooseneck
bpy.ops.mesh.primitive_cube_add(size=1, location=(0.255, 0.20, 0.62))
right_gooseneck = bpy.context.active_object
right_gooseneck.name = "right_gooseneck"
right_gooseneck.scale = (0.04, 0.04, 0.20)
add_bevel(right_gooseneck)
right_gooseneck.data.materials.append(wood)

# COMPONENT: left_armrest
cu_lar = bpy.data.curves.new("left_armrest_curve", type='CURVE')
cu_lar.dimensions = '3D'
cu_lar.bevel_depth = 0.02
sp_lar = cu_lar.splines.new('BEZIER')
sp_lar.bezier_points.add(1)
pts_lar = [(-0.255, -0.20, 0.72), (-0.275, 0.35, 0.72)]
for bp, co in zip(sp_lar.bezier_points, pts_lar):
    bp.co = co
    bp.handle_left_type = 'AUTO'
    bp.handle_right_type = 'AUTO'
left_armrest = bpy.data.objects.new("left_armrest", cu_lar)
bpy.context.collection.objects.link(left_armrest)
left_armrest.data.materials.append(wood)

# COMPONENT: right_armrest
cu_rar = bpy.data.curves.new("right_armrest_curve", type='CURVE')
cu_rar.dimensions = '3D'
cu_rar.bevel_depth = 0.02
sp_rar = cu_rar.splines.new('BEZIER')
sp_rar.bezier_points.add(1)
pts_rar = [(0.255, -0.20, 0.72), (0.275, 0.35, 0.72)]
for bp, co in zip(sp_rar.bezier_points, pts_rar):
    bp.co = co
    bp.handle_left_type = 'AUTO'
    bp.handle_right_type = 'AUTO'
right_armrest = bpy.data.objects.new("right_armrest", cu_rar)
bpy.context.collection.objects.link(right_armrest)
right_armrest.data.materials.append(wood)

# COMPONENT: left_apron
bpy.ops.mesh.primitive_cube_add(size=1, location=(-0.265, 0, 0.47))
left_apron = bpy.context.active_object
left_apron.name = "left_apron"
left_apron.scale = (0.02, 0.40, 0.08)
add_bevel(left_apron)
left_apron.data.materials.append(wood)

# COMPONENT: right_apron
bpy.ops.mesh.primitive_cube_add(size=1, location=(0.265, 0, 0.47))
right_apron = bpy.context.active_object
right_apron.name = "right_apron"
right_apron.scale = (0.02, 0.40, 0.08)
add_bevel(right_apron)
right_apron.data.materials.append(wood)

# COMPONENT: left_stretcher
bpy.ops.mesh.primitive_cylinder_add(radius=0.015, depth=0.48, location=(-0.255, 0, 0.20))
left_stretcher = bpy.context.active_object
left_stretcher.name = "left_stretcher"
left_stretcher.rotation_euler = (math.pi/2, 0, 0)
left_stretcher.data.materials.append(wood)

# COMPONENT: right_stretcher
bpy.ops.mesh.primitive_cylinder_add(radius=0.015, depth=0.48, location=(0.255, 0, 0.20))
right_stretcher = bpy.context.active_object
right_stretcher.name = "right_stretcher"
right_stretcher.rotation_euler = (math.pi/2, 0, 0)
right_stretcher.data.materials.append(wood)