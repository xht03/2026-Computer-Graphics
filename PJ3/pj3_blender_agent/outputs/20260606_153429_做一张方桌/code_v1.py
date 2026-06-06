import bpy
import math

# Clear the scene
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

# Material
mat = bpy.data.materials.new(name="WoodMat")
mat.use_nodes = True
bsdf = mat.node_tree.nodes["Principled BSDF"]
bsdf.inputs["Base Color"].default_value = (0.45, 0.25, 0.1, 1.0)
bsdf.inputs["Roughness"].default_value = 0.6

# Constants
leg_height = 0.75
top_thickness = 0.02

# COMPONENT: top_board
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, leg_height + top_thickness / 2))
top_board = bpy.context.active_object
top_board.name = "top_board"
top_board.scale = (0.9 / 2, 0.9 / 2, top_thickness / 2)
top_board.data.materials.append(mat)

# COMPONENT: apron (positioned as lower shelf)
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, leg_height * 0.25))
apron = bpy.context.active_object
apron.name = "apron"
apron.scale = (0.9 / 2, 0.9 / 2, 0.02 / 2)
apron.data.materials.append(mat)

# COMPONENT: leg
leg_radius = 0.04

for i in range(4):
    x = 0.45 * math.cos(i * math.pi / 2)
    y = 0.45 * math.sin(i * math.pi / 2)
    bpy.ops.mesh.primitive_cylinder_add(radius=leg_radius, depth=leg_height, location=(x, y, leg_height / 2))
    leg = bpy.context.active_object
    leg.name = f"leg_{i+1}"
    leg.data.materials.append(mat)