"""
Provides build_render_script() which returns a bpy code snippet.
Appended to generated model code before running in Blender.
"""

_RENDER_SNIPPET = r"""
# ── Multi-view render (injected by render.py) ──────────────────────────────
import bpy
import os
import math

_OUTPUT_DIR = {output_dir!r}
os.makedirs(_OUTPUT_DIR, exist_ok=True)

scene = bpy.context.scene
# Use EEVEE for reliable headless rendering (CYCLES needs GPU driver in headless)
scene.render.engine = 'BLENDER_EEVEE'
scene.eevee.taa_render_samples = 128
scene.render.image_settings.file_format = 'PNG'
scene.render.resolution_x = 1024
scene.render.resolution_y = 1024
scene.render.film_transparent = False

# Ensure a world background exists
if scene.world is None:
    scene.world = bpy.data.worlds.new("World")
scene.world.use_nodes = True
bg = scene.world.node_tree.nodes.get("Background")
if bg:
    bg.inputs[0].default_value = (0.05, 0.05, 0.05, 1)

# Remove existing cameras and lights to start clean
for obj in list(bpy.data.objects):
    if obj.type in ('CAMERA', 'LIGHT'):
        bpy.data.objects.remove(obj, do_unlink=True)

# Add key + fill lights
sun_data = bpy.data.lights.new("Sun", type='SUN')
sun_data.energy = 3.0
sun_obj = bpy.data.objects.new("Sun", sun_data)
scene.collection.objects.link(sun_obj)
sun_obj.location = (4, -4, 8)
sun_obj.rotation_euler = (math.radians(45), 0, math.radians(45))

fill_data = bpy.data.lights.new("Fill", type='AREA')
fill_data.energy = 200
fill_obj = bpy.data.objects.new("Fill", fill_data)
scene.collection.objects.link(fill_obj)
fill_obj.location = (-4, 4, 4)

# Add camera
cam_data = bpy.data.cameras.new("Camera")
cam_data.lens = 50
cam_obj = bpy.data.objects.new("Camera", cam_data)
scene.collection.objects.link(cam_obj)
scene.camera = cam_obj

def _look_at(cam, target=(0, 0, 0.8)):
    import mathutils
    direction = mathutils.Vector(target) - cam.location
    rot = direction.to_track_quat('-Z', 'Y')
    cam.rotation_euler = rot.to_euler()

_VIEWS = [
    ("front", (0, -4.5, 1.2)),
    ("side",  (4.5, 0,  1.2)),
    ("top",   (0,   0,  6.0)),
    ("iso",   (3.5, -3.5, 3.0)),
]

for name, loc in _VIEWS:
    cam_obj.location = loc
    _look_at(cam_obj)
    scene.render.filepath = os.path.join(_OUTPUT_DIR, name)
    bpy.ops.render.render(write_still=True)
    print(f"RENDER_OK: {{name}}.png")

print("ALL_RENDERS_COMPLETE")
"""


def build_render_script(output_dir: str) -> str:
    """Return bpy render code snippet with output_dir interpolated."""
    return _RENDER_SNIPPET.format(output_dir=output_dir)
