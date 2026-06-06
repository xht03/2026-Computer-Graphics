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
apron.scale = (0.9 / 2, 0.9 / 2, 0.02 / 2)
apron.data.materials.append(mat)

# COMPONENT: leg
leg_height = 0.78 * 0.65
bpy.ops.mesh.primitive_cylinder_add(radius=0.08 / 2, depth=leg_height / 2, location=(0.45, 0.45, leg_height / 2))
leg1 = bpy.context.active_object
leg1.name = "leg_1"
leg1.scale = (0.08 / 2, 0.08 / 2, leg_height / 2)
leg1.data.materials.append(mat)

bpy.ops.mesh.primitive_cylinder_add(radius=0.08 / 2, depth=leg_height / 2, location=(-0.45, -0.45, leg_height / 2))
leg2 = bpy.context.active_object
leg2.name = "leg_2"
leg2.scale = (0.08 / 2, 0.08 / 2, leg_height / 2)
leg2.data.materials.append(mat)

bpy.ops.mesh.primitive_cylinder_add(radius=0.08 / 2, depth=leg_height / 2, location=(-0.45, 0.45, leg_height / 2))
leg3 = bpy.context.active_object
leg3.name = "leg_3"
leg3.scale = (0.08 / 2, 0.08 / 2, leg_height / 2)
leg3.data.materials.append(mat)

bpy.ops.mesh.primitive_cylinder_add(radius=0.08 / 2, depth=leg_height / 2, location=(0.45, -0.45, leg_height / 2))
leg4 = bpy.context.active_object
leg4.name = "leg_4"
leg4.scale = (0.08 / 2, 0.08 / 2, leg_height / 2)
leg4.data.materials.append(mat)