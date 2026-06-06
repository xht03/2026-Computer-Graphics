import bpy
import math

# Clear the scene
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

# Create a wood-like material
mat = bpy.data.materials.new(name="WoodMat")
mat.use_nodes = True
bsdf = mat.node_tree.nodes["Principled BSDF"]
bsdf.inputs["Base Color"].default_value = (0.45, 0.25, 0.1, 1.0)
bsdf.inputs["Roughness"].default_value = 0.6

# COMPONENT: top_board
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, 0.76))
top_board = bpy.context.active_object
top_board.name = "top_board"
top_board.scale = (0.9 / 2, 0.9 / 2, 0.05)
top_board.data.materials.append(mat)

# COMPONENT: leg
leg_radius = 0.05
leg_height = 0.75
leg_positions = [(0.45, 0.45), (-0.45, 0.45), (-0.45, -0.45), (0.45, -0.45)]

leg_z = top_board.location.z - (top_board.dimensions.z / 2) - (leg_height / 2)

for i, (x, y) in enumerate(leg_positions):
    bpy.ops.mesh.primitive_cylinder_add(radius=leg_radius, depth=leg_height, location=(x, y, leg_z))
    leg = bpy.context.active_object
    leg.name = f"leg_{i+1}"
    leg.scale = (1, 1, 1)
    leg.data.materials.append(mat)