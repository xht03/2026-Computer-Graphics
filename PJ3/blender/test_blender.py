"""
Minimal bpy test script — run inside Blender to verify the environment:
    blender --background --python blender/test_blender.py

Generates a red cube and renders it to outputs/test_render.png
"""
import bpy
import os
import math

out_path = os.path.join(os.path.dirname(__file__), "..", "outputs", "test_render.png")
out_path = os.path.abspath(out_path)
os.makedirs(os.path.dirname(out_path), exist_ok=True)

# Clear scene
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

# Add red cube
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, 0.5))
cube = bpy.context.active_object

mat = bpy.data.materials.new(name="RedMat")
mat.use_nodes = True
bsdf = mat.node_tree.nodes.get("Principled BSDF")
bsdf.inputs["Base Color"].default_value = (0.8, 0.05, 0.05, 1.0)
cube.data.materials.append(mat)

# Light
bpy.ops.object.light_add(type='SUN', location=(4, -4, 8))
sun = bpy.context.active_object
sun.rotation_euler = (math.radians(45), 0, math.radians(45))
bpy.context.active_object.data.energy = 3.0

# Camera
bpy.ops.object.camera_add(location=(3, -4, 2.5))
cam = bpy.context.active_object
import mathutils
direction = mathutils.Vector((0, 0, 0.5)) - cam.location
cam.rotation_euler = direction.to_track_quat('-Z', 'Y').to_euler()
bpy.context.scene.camera = cam

# World background
world = bpy.context.scene.world
if not world:
    world = bpy.data.worlds.new("World")
    bpy.context.scene.world = world
world.use_nodes = True
bg = world.node_tree.nodes.get("Background")
if bg:
    bg.inputs[0].default_value = (0.1, 0.1, 0.1, 1)

# Render
scene = bpy.context.scene
scene.render.engine = 'CYCLES'
scene.cycles.samples = 32
scene.render.resolution_x = 512
scene.render.resolution_y = 512
scene.render.image_settings.file_format = 'PNG'
scene.render.filepath = out_path

bpy.ops.render.render(write_still=True)
print(f"TEST_RENDER_OK: {out_path}")
