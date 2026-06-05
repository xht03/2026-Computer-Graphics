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
bpy.ops.mesh.primitive_cube_add(size=1, scale=(1.0, 1.0, 0.05), location=(0, 0, 0.75))
top_board = bpy.context.active_object
top_board.name = "top_board"
top_board.scale = (0.9, 0.9, 0.05)
top_board.data.materials.append(mat)

# COMPONENT: apron
bpy.ops.mesh.primitive_cube_add(size=1, scale=(0.9, 0.3, 0.05), location=(0, 0, 0.78 - 0.05 / 2))
apron = bpy.context.active_object
apron.name = "apron"
apron.scale = (0.9, 0.3, 0.05)
apron.data.materials.append(mat)

# COMPONENT: leg
leg_radius = 0.04
leg_height = 1.5
leg_curve = 0.02  # Slight curve at the base
for i in range(4):
    bpy.ops.mesh.primitive_cube_add(size=1, scale=(0.08, 0.08, leg_height), location=(-0.4 * math.cos(i * math.pi / 2), -0.4 * math.sin(i * math.pi / 2), 0.375))
    leg = bpy.context.active_object
    leg.name = f"leg_{i+1}"
# [sanitized]     leg.scale.z -= leg_curve
# [sanitized]     leg.location.z += leg_curve
    leg.data.materials.append(mat)