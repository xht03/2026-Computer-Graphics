import bpy
import math

# Clear scene
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

# Create wood material
wood_mat = bpy.data.materials.new(name="Wood")
wood_mat.use_nodes = True
principled = wood_mat.node_tree.nodes["Principled BSDF"]
principled.inputs['Base Color'].default_value = (0.5, 0.3, 0.1, 1.0)
principled.inputs['Roughness'].default_value = 0.5

# COMPONENT: top_board
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, 0.77))
top_board = bpy.context.active_object
top_board.name = "top_board"
top_board.scale = (0.9, 0.9, 0.02)
top_board.data.materials.append(wood_mat)

# COMPONENT: apron
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, 0.75))
apron = bpy.context.active_object
apron.name = "apron"
apron.scale = (0.9, 0.9, 0.02)
apron.data.materials.append(wood_mat)

# COMPONENT: legs
# Table half-width is 0.45m, leg radius is 0.04m, so leg centers are at 0.45 - 0.04 = 0.41m
leg_offsets = 0.41
leg_height = 0.45
leg_radius = 0.04

leg_positions = [
    (leg_offsets, leg_offsets),
    (leg_offsets, -leg_offsets),
    (-leg_offsets, leg_offsets),
    (-leg_offsets, -leg_offsets)
]

for i, (x, y) in enumerate(leg_positions):
    bpy.ops.mesh.primitive_cylinder_add(
        radius=leg_radius, 
        depth=leg_height, 
        location=(x, y, leg_height / 2)
    )
    leg = bpy.context.active_object
    leg.name = f"leg_{i+1}"
    leg.data.materials.append(wood_mat)

# CONNECTIONS
CONNECTIONS = [("apron", "top_board")]