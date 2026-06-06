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
leg_h = 0.78 - 0.02
leg_w, leg_d = 0.08, 0.08
for i in range(4):
    bpy.ops.mesh.primitive_cylinder_add(radius=leg_w / 2, depth=leg_h, location=(0.45 * math.cos(i * math.pi / 2), 0.45 * math.sin(i * math.pi / 2), leg_h / 2))
    leg = bpy.context.active_object
    leg.name = f"leg_{i+1}"
    leg.scale = (leg_w / 2, leg_d / 2, leg_h / 2)
    leg.data.materials.append(mat)