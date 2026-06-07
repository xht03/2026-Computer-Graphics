import bpy
import math

bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

wood = bpy.data.materials.new(name="MingWood")
bsdf = wood.node_tree.nodes["Principled BSDF"]
bsdf.inputs['Base Color'].default_value = (0.45, 0.28, 0.12, 1.0)
bsdf.inputs['Roughness'].default_value = 0.5

def make_box(name, sx, sy, sz, x, y, z):
    bpy.ops.mesh.primitive_cube_add(size=1, location=(x, y, z))
    obj = bpy.context.active_object
    obj.name = name
    obj.scale = (sx, sy, sz)
    obj.data.materials.append(wood)
    return obj

seat_w, seat_d = 0.55, 0.48
leg_h = 0.55
leg_x = (seat_w/2) - 0.05
leg_y_front = (seat_d/2) - 0.05
leg_y_back = seat_d/2

# COMPONENT: leg_front_left
leg_fl = make_box("leg_front_left", 0.1, 0.1, leg_h, -leg_x, leg_y_front, leg_h/2)

# COMPONENT: leg_front_right
leg_fr = make_box("leg_front_right", 0.1, 0.1, leg_h, leg_x, leg_y_front, leg_h/2)

# COMPONENT: leg_back_left
leg_bl = make_box("leg_back_left", 0.1, 0.1, leg_h, -leg_x, -leg_y_back, leg_h/2)

# COMPONENT: leg_back_right
leg_br = make_box("leg_back_right", 0.1, 0.1, leg_h, leg_x, -leg_y_back, leg_h/2)

# COMPONENT: seat
seat_thick = 0.04
seat = make_box("seat", seat_w, seat_d, seat_thick, 0, 0, leg_h + seat_thick/2)

# COMPONENT: back_panel
# Positioned between back legs, centered behind seat
back_h = 0.6
back_thick = 0.04
back_y = -(seat_d/2 + back_thick/2)
back_panel = make_box("back_panel", seat_w, back_thick, back_h, 0, back_y, back_h/2)

# COMPONENT: armrest_left
# Extends from back panel to front leg
arm_w = 0.1
arm_h = 0.1
arm_z = leg_h
arm_y_start = back_y
arm_y_end = leg_y_front
arm_len = arm_y_end - arm_y_start
arm_y_center = (arm_y_start + arm_y_end) / 2
arm_x = -(seat_w/2 + arm_w/2)
armrest_l = make_box("armrest_left", arm_w, arm_len, arm_h, arm_x, arm_y_center, arm_z)

# COMPONENT: armrest_right
armrest_r = make_box("armrest_right", arm_w, arm_len, arm_h, -arm_x, arm_y_center, arm_z)

# COMPONENT: top_rail
# Curved bezier curve, beveled to 0.1m diameter
bpy.ops.curve.primitive_bezier_curve_add(location=(0, back_y, back_h))
top_rail = bpy.context.active_object
top_rail.name = "top_rail"
curve = top_rail.data
curve.dimensions = '3D'
curve.resolution_u = 64
curve.fill_mode = 'FULL'
curve.bevel_depth = 0.05

# Configure 3-point curve for forward arch (hat shape)
spline = curve.splines[0]
spline.bezier_points.add(1)

# Left end (local coordinates)
spline.bezier_points[0].co = (-seat_w/2, 0, 0)
spline.bezier_points[0].handle_right = (-0.09, 0, 0)
spline.bezier_points[0].handle_left_type = 'FREE'
spline.bezier_points[0].handle_right_type = 'FREE'

# Center (forward by 0.1m)
spline.bezier_points[1].co = (0, 0.1, 0)
spline.bezier_points[1].handle_left = (-0.09, 0.1, 0)
spline.bezier_points[1].handle_right = (0.09, 0.1, 0)
spline.bezier_points[1].handle_left_type = 'FREE'
spline.bezier_points[1].handle_right_type = 'FREE'

# Right end
spline.bezier_points[2].co = (seat_w/2, 0, 0)
spline.bezier_points[2].handle_left = (0.09, 0, 0)
spline.bezier_points[2].handle_left_type = 'FREE'
spline.bezier_points[2].handle_right_type = 'FREE'

# Convert to mesh for material application
bpy.ops.object.convert(target='MESH')
top_rail = bpy.context.active_object
top_rail.scale = (1.2, 1.0, 0.35)
bpy.ops.object.transform_apply(scale=True)
top_rail.data.materials.append(wood)

CONNECTIONS = [
    ("leg_front_left", "seat"),
    ("leg_front_right", "seat"),
    ("leg_back_left", "seat"),
    ("leg_back_right", "seat"),
    ("leg_back_left", "back_panel"),
    ("leg_back_right", "back_panel"),
    ("seat", "back_panel"),
    ("back_panel", "top_rail"),
    ("back_panel", "armrest_left"),
    ("back_panel", "armrest_right"),
    ("armrest_left", "leg_front_left"),
    ("armrest_right", "leg_front_right")
]