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
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, 0.78 / 2))
top_board = bpy.context.active_object
top_board.name = "top_board"
top_board.scale = (0.9 / 2, 0.9 / 2, 0.02 / 2)
top_board.data.materials.append(mat)

# COMPONENT: apron
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, 0.78 - 0.02 / 2))
apron = bpy.context.active_object
apron.name = "apron"
apron.scale = (0.9 / 2, 0.3 / 2, 0.02 / 2)
apron.data.materials.append(mat)

# COMPONENT: leg
leg_radius = 0.08 / 2
leg_height = 0.7 / 2
leg_curve = 0.02  # Slight curve at the base
leg_curve_angle = math.radians(10)  # 10 degrees

for i in range(4):
    bpy.ops.mesh.primitive_cylinder_add(radius=leg_radius, depth=leg_height, location=(0.45 * math.cos(math.radians(i * 90)), 0.45 * math.sin(math.radians(i * 90)), 0.78 - leg_height - leg_curve))
    leg = bpy.context.active_object
    leg.name = f"leg_{i+1}"
    leg.scale = (leg_radius, leg_radius, leg_height)
    leg.data.materials.append(mat)
    # Create a slight curve at the base
# [sanitized]     leg.location.rotate((1, 0, 0), leg_curve_angle)