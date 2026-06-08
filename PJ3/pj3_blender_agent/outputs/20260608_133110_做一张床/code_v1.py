import bpy
import math

bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

# Material: Dark Hardwood (Hongmu/Zitan style)
wood_mat = bpy.data.materials.new(name="DarkHardwood")
wood_mat.use_nodes = True
bsdf = wood_mat.node_tree.nodes["Principled BSDF"]
bsdf.inputs["Base Color"].default_value = (0.12, 0.06, 0.03, 1.0)
bsdf.inputs["Roughness"].default_value = 0.5

# Material: Cotton Mattress
mattress_mat = bpy.data.materials.new(name="CottonMattress")
mattress_mat.use_nodes = True
bsdf2 = mattress_mat.node_tree.nodes["Principled BSDF"]
bsdf2.inputs["Base Color"].default_value = (0.9, 0.85, 0.75, 1.0)
bsdf2.inputs["Roughness"].default_value = 0.9

# Post positions (centers): FrontLeft, FrontRight, BackLeft, BackRight
# Bed dimensions: X=1.5 (width), Y=2.0 (length)
# Posts inset slightly so canopy overhangs
post_coords = [(-0.70, -0.95), (0.70, -0.95), (-0.70, 0.95), (0.70, 0.95)]
post_names = []

# COMPONENT: corner_post
for i, (x, y) in enumerate(post_coords):
    bpy.ops.mesh.primitive_cube_add(size=1, location=(x, y, 0.5))
    post = bpy.context.active_object
    post.name = f"corner_post_{i+1}"
    post.scale = (0.15, 0.15, 1.0)
    post.data.materials.append(wood_mat)
    # Bevel for rounded edges
    mod = post.modifiers.new("Bevel", type='BEVEL')
    mod.width = 0.005
    mod.segments = 2
    post_names.append(post.name)

# COMPONENT: bed_frame
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, 0.45))
bed_frame = bpy.context.active_object
bed_frame.name = "bed_frame"
bed_frame.scale = (1.5, 2.0, 0.08)
bed_frame.data.materials.append(wood_mat)
mod = bed_frame.modifiers.new("Bevel", type='BEVEL')
mod.width = 0.004
mod.segments = 2

# COMPONENT: bed_board
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, 0.50))
bed_board = bpy.context.active_object
bed_board.name = "bed_board"
bed_board.scale = (1.46, 1.96, 0.02)
bed_board.data.materials.append(wood_mat)

# COMPONENT: mattress
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, 0.57))
mattress = bpy.context.active_object
mattress.name = "mattress"
mattress.scale = (1.4, 1.9, 0.12)
mattress.data.materials.append(mattress_mat)

# COMPONENT: headboard_panel
# At head side (y=0.95), between back posts (3 and 4)
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0.91, 0.8))
headboard = bpy.context.active_object
headboard.name = "headboard_panel"
headboard.scale = (1.4, 0.04, 0.6)
headboard.data.materials.append(wood_mat)
mod = headboard.modifiers.new("Bevel", type='BEVEL')
mod.width = 0.003
mod.segments = 2

# COMPONENT: footboard_panel
# At foot side (y=-0.95), between front posts (1 and 2)
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, -0.91, 0.8))
footboard = bpy.context.active_object
footboard.name = "footboard_panel"
footboard.scale = (1.4, 0.04, 0.6)
footboard.data.materials.append(wood_mat)
mod = footboard.modifiers.new("Bevel", type='BEVEL')
mod.width = 0.003
mod.segments = 2

# COMPONENT: canopy_rails
# Connect posts at top for visual stability and detail
bpy.ops.mesh.primitive_cube_add(size=1, location=(-0.70, 0, 1.0))
left_rail = bpy.context.active_object
left_rail.name = "left_rail"
left_rail.scale = (0.08, 1.9, 0.06)
left_rail.data.materials.append(wood_mat)

bpy.ops.mesh.primitive_cube_add(size=1, location=(0.70, 0, 1.0))
right_rail = bpy.context.active_object
right_rail.name = "right_rail"
right_rail.scale = (0.08, 1.9, 0.06)
right_rail.data.materials.append(wood_mat)

bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0.95, 1.0))
head_rail = bpy.context.active_object
head_rail.name = "head_rail"
head_rail.scale = (1.46, 0.08, 0.06)
head_rail.data.materials.append(wood_mat)

bpy.ops.mesh.primitive_cube_add(size=1, location=(0, -0.95, 1.0))
foot_rail = bpy.context.active_object
foot_rail.name = "foot_rail"
foot_rail.scale = (1.46, 0.08, 0.06)
foot_rail.data.materials.append(wood_mat)

# COMPONENT: blanket
# Thick comforter, axis-aligned, distinct warm color (not grey)
blanket_mat = bpy.data.materials.new(name="Blanket")
blanket_mat.use_nodes = True
bsdf3 = blanket_mat.node_tree.nodes["Principled BSDF"]
bsdf3.inputs["Base Color"].default_value = (0.6, 0.25, 0.2, 1.0)
bsdf3.inputs["Roughness"].default_value = 0.8

bpy.ops.mesh.primitive_cube_add(size=1, location=(0, -0.1, 0.66))
blanket = bpy.context.active_object
blanket.name = "blanket"
blanket.scale = (1.38, 1.4, 0.06)
blanket.data.materials.append(blanket_mat)

# COMPONENT: pillows
# Two pillows at head of bed
for px in [-0.35, 0.35]:
    bpy.ops.mesh.primitive_cube_add(size=1, location=(px, 0.75, 0.68))
    pillow = bpy.context.active_object
    pillow.name = f"pillow_{'left' if px < 0 else 'right'}"
    pillow.scale = (0.3, 0.2, 0.08)
    pillow.data.materials.append(mattress_mat)
    mod = pillow.modifiers.new("Bevel", type='BEVEL')
    mod.width = 0.015
    mod.segments = 3