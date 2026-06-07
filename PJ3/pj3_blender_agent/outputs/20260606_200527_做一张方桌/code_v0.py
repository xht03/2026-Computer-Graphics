import bpy
import math

# Clear scene
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

# Dimensions
total_width = 0.9
total_depth = 0.9
total_height = 0.78
top_thickness = 0.02

# Leg dimensions (diameter 0.08m, radius 0.04m)
leg_radius = 0.04
# Apron height from specification
apron_height = 0.3
# Derive leg height: must fill space from floor to apron bottom
# Total = leg_height + apron_height + top_thickness
leg_height = total_height - top_thickness - apron_height  # 0.46m

# Positions
top_z = total_height - top_thickness / 2  # 0.77
leg_z = leg_height / 2  # 0.23
apron_z = leg_height + apron_height / 2  # 0.61

# Leg inset from center (half width minus radius)
leg_inset = total_width / 2 - leg_radius  # 0.41

# Wood material
wood_mat = bpy.data.materials.new(name="Wood")
wood_mat.use_nodes = True
bsdf = wood_mat.node_tree.nodes["Principled BSDF"]
bsdf.inputs["Base Color"].default_value = (0.45, 0.25, 0.1, 1.0)
bsdf.inputs["Roughness"].default_value = 0.6

# COMPONENT: top_board
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, top_z))
top_board = bpy.context.active_object
top_board.name = "top_board"
top_board.scale = (total_width, total_depth, top_thickness)
top_board.data.materials.append(wood_mat)

# COMPONENT: apron (frame composed of 4 boards, joined into one object)
apron_thickness = 0.02

# Front apron
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, leg_inset + apron_thickness/2, apron_z))
apron_front = bpy.context.active_object
apron_front.scale = (total_width, apron_thickness, apron_height)

# Back apron
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, -leg_inset - apron_thickness/2, apron_z))
apron_back = bpy.context.active_object
apron_back.scale = (total_width, apron_thickness, apron_height)

# Left apron (shortened to avoid overlapping front/back at corners)
bpy.ops.mesh.primitive_cube_add(size=1, location=(-leg_inset - apron_thickness/2, 0, apron_z))
apron_left = bpy.context.active_object
apron_left.scale = (apron_thickness, total_depth - 2*apron_thickness, apron_height)

# Right apron
bpy.ops.mesh.primitive_cube_add(size=1, location=(leg_inset + apron_thickness/2, 0, apron_z))
apron_right = bpy.context.active_object
apron_right.scale = (apron_thickness, total_depth - 2*apron_thickness, apron_height)

# Join all apron parts into single object
bpy.ops.object.select_all(action='DESELECT')
apron_front.select_set(True)
apron_back.select_set(True)
apron_left.select_set(True)
apron_right.select_set(True)
bpy.context.view_layer.objects.active = apron_front
bpy.ops.object.join()
apron = bpy.context.active_object
apron.name = "apron"
apron.data.materials.append(wood_mat)

# COMPONENT: leg (4 cylinders at corners)
leg_positions = [
    (leg_inset, leg_inset),
    (leg_inset, -leg_inset),
    (-leg_inset, leg_inset),
    (-leg_inset, -leg_inset)
]

for i, (x, y) in enumerate(leg_positions):
    bpy.ops.mesh.primitive_cylinder_add(radius=leg_radius, depth=leg_height, location=(x, y, leg_z))
    leg = bpy.context.active_object
    leg.name = f"leg_{i+1}"
    leg.data.materials.append(wood_mat)

# Declare connections for assembly verification
CONNECTIONS = [
    ("apron", "top_board"),
    ("leg_1", "apron"),
    ("leg_2", "apron"),
    ("leg_3", "apron"),
    ("leg_4", "apron")
]