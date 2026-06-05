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
top_board.scale = (0.9 / 2, 0.9 / 2, 0.05 / 2)
top_board.data.materials.append(mat)

# COMPONENT: apron
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, 0.78 - 0.05 / 2))
apron = bpy.context.active_object
apron.name = "apron"
apron.scale = (0.9 / 2, 0.3 / 2, 0.05 / 2)
apron.data.materials.append(mat)

# COMPONENT: leg_base
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, 0))
leg_base = bpy.context.active_object
leg_base.name = "leg_base"
leg_base.scale = (0.9 / 2, 0.9 / 2, 0.05 / 2)
leg_base.data.materials.append(mat)

# Create four vertical legs
leg_radius = 0.05 / 2
leg_height = 0.7 / 2
leg_curve = 0.02  # Slight curve at the base
leg_curve_angle = math.radians(10)  # 10 degrees

# Leg positions
leg_positions = [
    (0.4, 0.4, 0.35),
    (-0.4, 0.4, 0.35),
    (0.4, -0.4, 0.35),
    (-0.4, -0.4, 0.35)
]

# Create legs
for pos in leg_positions:
    bpy.ops.mesh.primitive_cube_add(size=1, location=pos)
    leg = bpy.context.active_object
    leg.name = f"leg_{leg_positions.index(pos)}"
    leg.scale = (0.08, 0.08, 0.7)
    leg.data.materials.append(mat)

# Reposition the top_board to rest on the legs
top_board.location.z = leg_height + 0.05 / 2

# Create structural legs connecting the tabletop to floor level
for pos in leg_positions:
    bpy.ops.mesh.primitive_cylinder_add(radius=0.04, depth=0.7, location=pos)
    leg = bpy.context.active_object
    leg.name = f"leg_structural_{leg_positions.index(pos)}"
    leg.scale = (0.08, 0.08, 0.7)
    leg.data.materials.append(mat)