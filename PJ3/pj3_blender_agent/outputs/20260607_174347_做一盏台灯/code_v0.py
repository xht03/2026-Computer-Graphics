import bpy
import math

# Clear scene
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

# Material setup
def make_material(name, color, roughness, metallic=0.0):
    mat = bpy.data.materials.new(name=name)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes["Principled BSDF"]
    bsdf.inputs["Base Color"].default_value = (*color, 1.0)
    bsdf.inputs["Roughness"].default_value = roughness
    bsdf.inputs["Metallic"].default_value = metallic
    return mat

wood_mat = make_material("DarkWood", (0.15, 0.08, 0.05), 0.3)
brass_mat = make_material("Brass", (0.8, 0.6, 0.2), 0.3, metallic=1.0)
fabric_mat = make_material("Silk", (0.9, 0.85, 0.7), 0.9)
glass_mat = make_material("Glass", (1.0, 1.0, 0.9), 0.1)

# COMPONENT: base
bpy.ops.mesh.primitive_cylinder_add(radius=0.08, depth=0.05, location=(0, 0, 0.025))
base = bpy.context.active_object
base.name = "base"
base.data.materials.append(wood_mat)
mod = base.modifiers.new("Bevel", type='BEVEL')
mod.width = 0.003
mod.segments = 2

# COMPONENT: stem
bpy.ops.mesh.primitive_cylinder_add(radius=0.0125, depth=0.35, location=(0, 0, 0.225))
stem = bpy.context.active_object
stem.name = "stem"
stem.data.materials.append(wood_mat)
mod = stem.modifiers.new("Bevel", type='BEVEL')
mod.width = 0.001
mod.segments = 2

# COMPONENT: bulb_socket
bpy.ops.mesh.primitive_cylinder_add(radius=0.015, depth=0.04, location=(0, 0, 0.42))
bulb_socket = bpy.context.active_object
bulb_socket.name = "bulb_socket"
bulb_socket.data.materials.append(brass_mat)

# COMPONENT: light_bulb
bpy.ops.mesh.primitive_uv_sphere_add(radius=0.025, location=(0, 0, 0.44))
light_bulb = bpy.context.active_object
light_bulb.name = "light_bulb"
light_bulb.data.materials.append(glass_mat)

# COMPONENT: shade_bottom_ring
bpy.ops.mesh.primitive_torus_add(mode='MAJOR_MINOR', major_radius=0.11, minor_radius=0.01, location=(0, 0, 0.45), major_segments=48, minor_segments=12)
shade_bottom_ring = bpy.context.active_object
shade_bottom_ring.name = "shade_bottom_ring"
shade_bottom_ring.data.materials.append(wood_mat)
mod = shade_bottom_ring.modifiers.new("Bevel", type='BEVEL')
mod.width = 0.002
mod.segments = 2

# COMPONENT: shade_fabric
bpy.ops.mesh.primitive_cone_add(radius1=0.11, radius2=0.06, depth=0.18, location=(0, 0, 0.54), vertices=32)
shade_fabric = bpy.context.active_object
shade_fabric.name = "shade_fabric"
shade_fabric.data.materials.append(fabric_mat)

# COMPONENT: shade_top_ring
bpy.ops.mesh.primitive_torus_add(mode='MAJOR_MINOR', major_radius=0.06, minor_radius=0.01, location=(0, 0, 0.62), major_segments=48, minor_segments=12)
shade_top_ring = bpy.context.active_object
shade_top_ring.name = "shade_top_ring"
shade_top_ring.data.materials.append(wood_mat)
mod = shade_top_ring.modifiers.new("Bevel", type='BEVEL')
mod.width = 0.002
mod.segments = 2

# COMPONENT: finial
bpy.ops.mesh.primitive_uv_sphere_add(radius=0.0125, location=(0, 0, 0.6425))
finial = bpy.context.active_object
finial.name = "finial"
finial.data.materials.append(brass_mat)

CONNECTIONS = [
    ("base", "stem"),
    ("stem", "bulb_socket"),
    ("bulb_socket", "light_bulb"),
    ("bulb_socket", "shade_bottom_ring"),
    ("shade_bottom_ring", "shade_fabric"),
    ("shade_fabric", "shade_top_ring"),
    ("shade_top_ring", "finial")
]