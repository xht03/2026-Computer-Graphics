import bpy
import math

# Clear scene
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

# Material setup
wood = bpy.data.materials.new(name="Ming_Hardwood")
wood.use_nodes = True
bsdf = wood.node_tree.nodes["Principled BSDF"]
bsdf.inputs["Base Color"].default_value = (0.35, 0.18, 0.08, 1.0)
bsdf.inputs["Roughness"].default_value = 0.4

# Common dimensions
seat_w, seat_d, seat_h = 0.55, 0.48, 0.50
seat_thick = 0.04
leg_size = 0.03  # Reduced from 0.04 for slender proportions
back_leg_h = 1.05  # Reduced from 1.15 for better proportions
front_leg_h = 0.70
splay = 0.025  # radians for leg splay

x_left, x_right = -seat_w/2, seat_w/2
y_front, y_back = -seat_d/2, seat_d/2

def add_bevel(obj, width=0.003):
    mod = obj.modifiers.new("Bevel", type='BEVEL')
    mod.width = width
    mod.segments = 2

# Bevel profile for curved aprons (cross-section: 0.02 thick x 0.10 high)
apron_profile = bpy.data.meshes.new("ApronProfile")
apron_profile_obj = bpy.data.objects.new("apron_profile_obj", apron_profile)
bpy.context.collection.objects.link(apron_profile_obj)
# Profile as thin box to avoid degenerate geometry
apron_profile.from_pydata(
    [(0,0,0), (0.02,0,0), (0.02,0.10,0), (0,0.10,0), (0,0,0.001), (0.02,0,0.001), (0.02,0.10,0.001), (0,0.10,0.001)],
    [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)],
    [(0,1,2,3),(4,7,6,5),(0,4,5,1),(1,5,6,2),(2,6,7,3),(3,7,4,0)]
)
apron_profile_obj.hide_set(True)

# Profile for back splat (wide panel: 0.12 wide x 0.015 thick)
splat_profile = bpy.data.meshes.new("SplatProfile")
splat_profile_obj = bpy.data.objects.new("splat_profile_obj", splat_profile)
bpy.context.collection.objects.link(splat_profile_obj)
# Profile as thin box to avoid degenerate geometry
splat_profile.from_pydata(
    [(0,0,0), (0.12,0,0), (0.12,0.015,0), (0,0.015,0), (0,0,0.001), (0.12,0,0.001), (0.12,0.015,0.001), (0,0.015,0.001)],
    [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)],
    [(0,1,2,3),(4,7,6,5),(0,4,5,1),(1,5,6,2),(2,6,7,3),(3,7,4,0)]
)
splat_profile_obj.hide_set(True)

# COMPONENT: back_left_leg
bpy.ops.mesh.primitive_cube_add(size=1, location=(x_left, y_back, back_leg_h/2))
back_left_leg = bpy.context.active_object
back_left_leg.name = "back_left_leg"
back_left_leg.scale = (leg_size, leg_size, back_leg_h)
back_left_leg.rotation_euler = (splay, -splay, 0)  # splay outward
add_bevel(back_left_leg)
back_left_leg.data.materials.append(wood)

# COMPONENT: back_right_leg
bpy.ops.mesh.primitive_cube_add(size=1, location=(x_right, y_back, back_leg_h/2))
back_right_leg = bpy.context.active_object
back_right_leg.name = "back_right_leg"
back_right_leg.scale = (leg_size, leg_size, back_leg_h)
back_right_leg.rotation_euler = (splay, splay, 0)
add_bevel(back_right_leg)
back_right_leg.data.materials.append(wood)

# COMPONENT: front_left_leg
bpy.ops.mesh.primitive_cube_add(size=1, location=(x_left, y_front, front_leg_h/2))
front_left_leg = bpy.context.active_object
front_left_leg.name = "front_left_leg"
front_left_leg.scale = (leg_size, leg_size, front_leg_h)
front_left_leg.rotation_euler = (-splay, -splay, 0)
add_bevel(front_left_leg)
front_left_leg.data.materials.append(wood)

# COMPONENT: front_right_leg
bpy.ops.mesh.primitive_cube_add(size=1, location=(x_right, y_front, front_leg_h/2))
front_right_leg = bpy.context.active_object
front_right_leg.name = "front_right_leg"
front_right_leg.scale = (leg_size, leg_size, front_leg_h)
front_right_leg.rotation_euler = (-splay, splay, 0)
add_bevel(front_right_leg)
front_right_leg.data.materials.append(wood)

# COMPONENT: seat_frame
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, seat_h))
seat_frame = bpy.context.active_object
seat_frame.name = "seat_frame"
seat_frame.scale = (seat_w, seat_d, seat_thick)
add_bevel(seat_frame, 0.003)
seat_frame.data.materials.append(wood)

# COMPONENT: top_rail (Hat-shaped crest rail with extended ends)
cu = bpy.data.curves.new("top_rail_curve", type='CURVE')
cu.dimensions = '3D'
cu.bevel_depth = 0.018
cu.resolution_u = 12
sp = cu.splines.new('BEZIER')
sp.bezier_points.add(2)
# Extended beyond rear posts for official's hat style
points = [(-0.40, y_back, back_leg_h+0.02), (0, y_back, back_leg_h-0.02), (0.40, y_back, back_leg_h+0.02)]
for bp, co in zip(sp.bezier_points, points):
    bp.co = co
    bp.handle_left_type = bp.handle_right_type = 'AUTO'
top_rail = bpy.data.objects.new("top_rail", cu)
bpy.context.collection.objects.link(top_rail)
top_rail.data.materials.append(wood)

# COMPONENT: back_splat (S-curved panel connecting seat to top rail)
cu = bpy.data.curves.new("back_splat_curve", type='CURVE')
cu.dimensions = '3D'
cu.bevel_mode = 'OBJECT'
cu.bevel_object = splat_profile_obj
cu.resolution_u = 12
sp = cu.splines.new('BEZIER')
sp.bezier_points.add(2)
# S-curve in X-Z plane at back, connecting seat to top rail
pts = [(0, y_back-0.02, seat_h+0.02), (0.04, y_back-0.02, seat_h+0.25), (-0.04, y_back-0.02, seat_h+0.40), (0, y_back-0.02, back_leg_h-0.03)]
for bp, co in zip(sp.bezier_points, pts):
    bp.co = co
    bp.handle_left_type = bp.handle_right_type = 'AUTO'
back_splat = bpy.data.objects.new("back_splat", cu)
bpy.context.collection.objects.link(back_splat)
back_splat.data.materials.append(wood)

# COMPONENT: left_armrest (Goose neck curve connecting front and back legs)
cu = bpy.data.curves.new("left_armrest_curve", type='CURVE')
cu.dimensions = '3D'
cu.bevel_depth = 0.012
sp = cu.splines.new('BEZIER')
sp.bezier_points.add(2)
# Connects front leg top to back leg top with graceful outward curve
pts = [(x_left-0.05, y_front-0.03, front_leg_h), (x_left-0.03, y_front+0.08, front_leg_h+0.05), (x_left-0.01, (y_front+y_back)/2, (front_leg_h+back_leg_h)/2), (x_left, y_back, back_leg_h-0.05)]
for bp, co in zip(sp.bezier_points, pts):
    bp.co = co
    bp.handle_left_type = bp.handle_right_type = 'AUTO'
left_armrest = bpy.data.objects.new("left_armrest", cu)
bpy.context.collection.objects.link(left_armrest)
left_armrest.data.materials.append(wood)

# COMPONENT: right_armrest
cu = bpy.data.curves.new("right_armrest_curve", type='CURVE')
cu.dimensions = '3D'
cu.bevel_depth = 0.012
sp = cu.splines.new('BEZIER')
sp.bezier_points.add(2)
pts = [(x_right+0.05, y_front-0.03, front_leg_h), (x_right+0.03, y_front+0.08, front_leg_h+0.05), (x_right+0.01, (y_front+y_back)/2, (front_leg_h+back_leg_h)/2), (x_right, y_back, back_leg_h-0.05)]
for bp, co in zip(sp.bezier_points, pts):
    bp.co = co
    bp.handle_left_type = bp.handle_right_type = 'AUTO'
right_armrest = bpy.data.objects.new("right_armrest", cu)
bpy.context.collection.objects.link(right_armrest)
right_armrest.data.materials.append(wood)

# COMPONENT: front_apron (More pronounced curve)
cu = bpy.data.curves.new("front_apron_curve", type='CURVE')
cu.dimensions = '3D'
cu.bevel_mode = 'OBJECT'
cu.bevel_object = apron_profile_obj
sp = cu.splines.new('BEZIER')
sp.bezier_points.add(1)
# Curved bottom edge, deeper concave upward
pts = [(x_left, y_front+0.01, 0.38), (x_left/2, y_front+0.01, 0.28), (x_right/2, y_front+0.01, 0.28), (x_right, y_front+0.01, 0.38)]
for bp, co in zip(sp.bezier_points, pts):
    bp.co = co
    bp.handle_left_type = bp.handle_right_type = 'AUTO'
front_apron = bpy.data.objects.new("front_apron", cu)
bpy.context.collection.objects.link(front_apron)
front_apron.data.materials.append(wood)

# COMPONENT: left_apron
cu = bpy.data.curves.new("left_apron_curve", type='CURVE')
cu.dimensions = '3D'
cu.bevel_mode = 'OBJECT'
cu.bevel_object = apron_profile_obj
sp = cu.splines.new('BEZIER')
sp.bezier_points.add(1)
pts = [(x_left-0.01, y_front, 0.40), (x_left-0.01, (y_front+y_back)/2, 0.32), (x_left-0.01, y_back, 0.40)]
for bp, co in zip(sp.bezier_points, pts):
    bp.co = co
    bp.handle_left_type = bp.handle_right_type = 'AUTO'
left_apron = bpy.data.objects.new("left_apron", cu)
bpy.context.collection.objects.link(left_apron)
left_apron.data.materials.append(wood)

# COMPONENT: right_apron
cu = bpy.data.curves.new("right_apron_curve", type='CURVE')
cu.dimensions = '3D'
cu.bevel_mode = 'OBJECT'
cu.bevel_object = apron_profile_obj
sp = cu.splines.new('BEZIER')
sp.bezier_points.add(1)
pts = [(x_right+0.01, y_front, 0.40), (x_right+0.01, (y_front+y_back)/2, 0.32), (x_right+0.01, y_back, 0.40)]
for bp, co in zip(sp.bezier_points, pts):
    bp.co = co
    bp.handle_left_type = bp.handle_right_type = 'AUTO'
right_apron = bpy.data.objects.new("right_apron", cu)
bpy.context.collection.objects.link(right_apron)
right_apron.data.materials.append(wood)

# COMPONENT: back_apron
cu = bpy.data.curves.new("back_apron_curve", type='CURVE')
cu.dimensions = '3D'
cu.bevel_mode = 'OBJECT'
cu.bevel_object = apron_profile_obj
sp = cu.splines.new('BEZIER')
sp.bezier_points.add(1)
pts = [(x_left, y_back-0.01, 0.40), (x_left/2, y_back-0.01, 0.36), (x_right/2, y_back-0.01, 0.36), (x_right, y_back-0.01, 0.40)]
for bp, co in zip(sp.bezier_points, pts):
    bp.co = co
    bp.handle_left_type = bp.handle_right_type = 'AUTO'
back_apron = bpy.data.objects.new("back_apron", cu)
bpy.context.collection.objects.link(back_apron)
back_apron.data.materials.append(wood)

# COMPONENT: front_stretcher (scaled to connect legs with slender proportions)
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, y_front+0.015, 0.15))
front_stretcher = bpy.context.active_object
front_stretcher.name = "front_stretcher"
front_stretcher.scale = (0.54, 0.02, 0.02)  # Increased X to connect, reduced YZ for slender look
add_bevel(front_stretcher)
front_stretcher.data.materials.append(wood)

# COMPONENT: back_stretcher
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, y_back-0.015, 0.15))
back_stretcher = bpy.context.active_object
back_stretcher.name = "back_stretcher"
back_stretcher.scale = (0.54, 0.02, 0.02)
add_bevel(back_stretcher)
back_stretcher.data.materials.append(wood)

# COMPONENT: left_stretcher
bpy.ops.mesh.primitive_cube_add(size=1, location=(x_left-0.015, 0, 0.18))
left_stretcher = bpy.context.active_object
left_stretcher.name = "left_stretcher"
left_stretcher.scale = (0.02, 0.46, 0.02)  # Slightly longer to connect slender legs
add_bevel(left_stretcher)
left_stretcher.data.materials.append(wood)

# COMPONENT: right_stretcher
bpy.ops.mesh.primitive_cube_add(size=1, location=(x_right+0.015, 0, 0.18))
right_stretcher = bpy.context.active_object
right_stretcher.name = "right_stretcher"
right_stretcher.scale = (0.02, 0.46, 0.02)
add_bevel(right_stretcher)
right_stretcher.data.materials.append(wood)

# CONNECTIONS list as required
CONNECTIONS = [
    ("seat_frame", "back_left_leg"),
    ("seat_frame", "back_right_leg"),
    ("seat_frame", "front_left_leg"),
    ("seat_frame", "front_right_leg"),
    ("seat_frame", "front_apron"),
    ("seat_frame", "left_apron"),
    ("seat_frame", "right_apron"),
    ("seat_frame", "back_apron"),
    ("seat_frame", "back_splat"),
    ("back_left_leg", "top_rail"),
    ("back_right_leg", "top_rail"),
    ("back_left_leg", "back_splat"),
    ("back_right_leg", "back_splat"),
    ("back_left_leg", "left_armrest"),
    ("back_right_leg", "right_armrest"),
    ("front_left_leg", "left_armrest"),
    ("front_right_leg", "right_armrest"),
    ("back_left_leg", "back_apron"),
    ("back_right_leg", "back_apron"),
    ("front_left_leg", "front_apron"),
    ("front_right_leg", "front_apron"),
    ("back_left_leg", "left_apron"),
    ("front_left_leg", "left_apron"),
    ("back_right_leg", "right_apron"),
    ("front_right_leg", "right_apron"),
    ("front_left_leg", "front_stretcher"),
    ("front_right_leg", "front_stretcher"),
    ("back_left_leg", "back_stretcher"),
    ("back_right_leg", "back_stretcher"),
    ("front_left_leg", "left_stretcher"),
    ("back_left_leg", "left_stretcher"),
    ("front_right_leg", "right_stretcher"),
    ("back_right_leg", "right_stretcher"),
]