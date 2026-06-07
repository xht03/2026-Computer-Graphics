import bpy
import math

bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

def make_material(name, color, roughness=0.5, metallic=0.0):
    mat = bpy.data.materials.new(name=name)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes["Principled BSDF"]
    bsdf.inputs["Base Color"].default_value = (*color, 1.0)
    bsdf.inputs["Roughness"].default_value = roughness
    bsdf.inputs["Metallic"].default_value = metallic
    return mat

wood_dark = make_material("WoodDark", (0.15, 0.06, 0.03), 0.6)
wood_bamboo = make_material("Bamboo", (0.65, 0.45, 0.18), 0.7)
metal_brass = make_material("Brass", (0.4, 0.3, 0.1), 0.3, 0.8)
mat_bulb = make_material("BulbGlass", (0.95, 0.95, 0.9), 0.1)
mat_shade = make_material("RicePaper", (0.92, 0.88, 0.75), 0.9)

# COMPONENT: base
bpy.ops.mesh.primitive_cylinder_add(radius=0.10, depth=0.04, location=(0, 0, 0.02))
base = bpy.context.active_object
base.name = "base"
base.data.materials.append(wood_dark)
mod = base.modifiers.new("Bevel", type='BEVEL')
mod.width = 0.004
mod.segments = 2

# COMPONENT: stem
bpy.ops.mesh.primitive_cylinder_add(radius=0.015, depth=0.40, location=(0, 0, 0.24))
stem = bpy.context.active_object
stem.name = "stem"
stem.data.materials.append(wood_dark)

# COMPONENT: shade_frame_bottom
# Spindle torus geometry to ensure physical contact with centered stem
bpy.ops.mesh.primitive_torus_add(major_radius=0.01, minor_radius=0.115, location=(0, 0, 0.20), major_segments=64, minor_segments=16)
shade_frame_bottom = bpy.context.active_object
shade_frame_bottom.name = "shade_frame_bottom"
shade_frame_bottom.data.materials.append(wood_bamboo)

# COMPONENT: shade
bpy.ops.mesh.primitive_cylinder_add(radius=0.125, depth=0.22, location=(0, 0, 0.31))
shade = bpy.context.active_object
shade.name = "shade"
shade.data.materials.append(mat_shade)

# COMPONENT: shade_frame_top
bpy.ops.mesh.primitive_torus_add(major_radius=0.01, minor_radius=0.115, location=(0, 0, 0.42), major_segments=64, minor_segments=16)
shade_frame_top = bpy.context.active_object
shade_frame_top.name = "shade_frame_top"
shade_frame_top.data.materials.append(wood_bamboo)

# COMPONENT: socket
bpy.ops.mesh.primitive_cylinder_add(radius=0.02, depth=0.06, location=(0, 0, 0.41))
socket = bpy.context.active_object
socket.name = "socket"
socket.data.materials.append(metal_brass)

# COMPONENT: bulb
bpy.ops.mesh.primitive_uv_sphere_add(radius=0.03, segments=32, ring_count=16, location=(0, 0, 0.35))
bulb = bpy.context.active_object
bulb.name = "bulb"
bulb.data.materials.append(mat_bulb)

# COMPONENT: finial
bpy.ops.mesh.primitive_uv_sphere_add(radius=0.02, segments=16, ring_count=8, location=(0, 0, 0.44))
finial = bpy.context.active_object
finial.name = "finial"
finial.data.materials.append(wood_dark)

CONNECTIONS = [
    ("base", "stem"),
    ("stem", "socket"),
    ("socket", "bulb"),
    ("stem", "shade_frame_bottom"),
    ("shade_frame_bottom", "shade"),
    ("shade", "shade_frame_top"),
    ("shade_frame_top", "finial")
]