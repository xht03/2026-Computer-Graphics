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

def _look_at(cam, target):
    import mathutils
    direction = mathutils.Vector(target) - cam.location
    rot = direction.to_track_quat('-Z', 'Y')
    cam.rotation_euler = rot.to_euler()

# ── Auto-frame: fit the camera distance to the actual scene bounding box so the
# object FILLS the frame regardless of its real size (small parts no longer render
# as tiny specks). Lights are left fixed — only the camera moves.
import mathutils as _mu

_mesh_objs = [o for o in bpy.data.objects if o.type == 'MESH']
if _mesh_objs:
    _corners = []
    for _o in _mesh_objs:
        for _v in _o.bound_box:
            _corners.append(_o.matrix_world @ _mu.Vector(_v))
    _minx = min(c.x for c in _corners); _maxx = max(c.x for c in _corners)
    _miny = min(c.y for c in _corners); _maxy = max(c.y for c in _corners)
    _minz = min(c.z for c in _corners); _maxz = max(c.z for c in _corners)
    _center = _mu.Vector(((_minx + _maxx) / 2, (_miny + _maxy) / 2, (_minz + _maxz) / 2))
    _diag = _mu.Vector((_maxx - _minx, _maxy - _miny, _maxz - _minz)).length
    _radius = max(_diag / 2.0, 0.05)
else:
    _center = _mu.Vector((0.0, 0.0, 0.5))
    _radius = 0.5

# Distance so the bounding sphere fits the field of view (with margin).
_dist = _radius / math.tan(cam_data.angle / 2.0) * 1.25

# View directions (unit-ish vectors from the object centre); scaled by _dist.
_VIEWS = [
    ("front", _mu.Vector((0.0, -1.0, 0.25))),
    ("side",  _mu.Vector((1.0,  0.0, 0.25))),
    ("top",   _mu.Vector((0.0,  0.0, 1.0))),
    ("iso",   _mu.Vector((1.0, -1.0, 0.7))),
]

for name, _dir in _VIEWS:
    cam_obj.location = _center + _dir.normalized() * _dist
    _look_at(cam_obj, _center)
    scene.render.filepath = os.path.join(_OUTPUT_DIR, name)
    bpy.ops.render.render(write_still=True)
    print(f"RENDER_OK: {{name}}.png")

print("ALL_RENDERS_COMPLETE")
"""


def build_render_script(output_dir: str) -> str:
    """Return bpy render code snippet with output_dir interpolated."""
    return _RENDER_SNIPPET.format(output_dir=output_dir)
