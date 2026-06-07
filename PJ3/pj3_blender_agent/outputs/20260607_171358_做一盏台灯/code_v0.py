import bpy
import math

# Clear scene
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

# Material setup
def create_material(name, base_color, roughness=0.5, metallic=0.0, emission=0.0):
    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes["Principled BSDF"]
    bsdf.inputs["Base Color"].default_value = (*base_color, 1.0)
    bsdf.inputs["Roughness"].default_value = roughness
    bsdf.inputs["Metallic"].default_value = metallic
    if emission > 0:
        bsdf.inputs["Emission Color"].default_value = (*base_color, 1.0)
        bsdf.inputs["Emission Strength"].default_value = emission
    return mat

wood_mat = create_material("Mahogany", (0.12, 0.04, 0.02), roughness=0.3)
brass_mat = create_material("Brass", (0.75, 0.55, 0.15), roughness=0.3, metallic=1.0)
paper_mat = create_material("RicePaper", (0.95, 0.9, 0.75), roughness=0.9)
jade_mat = create_material("Jade", (0.25, 0.6, 0.4), roughness=0.1)
bulb_mat = create_material("BulbLight", (1.0, 0.95, 0.8), roughness=0.0, emission=2.0)

# Derived heights for assembly
total_h = 0.52
base_plate_h = 0.04
base_pedestal_h = 0.08
finial_r = 0.02
central_stem_h = total_h - base_plate_h - base_pedestal_h - finial_r * 2  # 0.36

# COMPONENT: base_plate
bpy.ops.mesh.primitive_cylinder_add(radius=0.09, depth=base_plate_h, location=(0, 0, base_plate_h/2))
base_plate = bpy.context.active_object
base_plate.name = "base_plate"
base_plate.data.materials.append(wood_mat)
mod = base_plate.modifiers.new("Bevel", type='BEVEL')
mod.width = 0.003
mod.segments = 2

# COMPONENT: base_pedestal
pedestal_z = base_plate_h + base_pedestal_h/2
bpy.ops.mesh.primitive_cylinder_add(radius=0.06, depth=base_pedestal_h, location=(0, 0, pedestal_z))
base_pedestal = bpy.context.active_object
base_pedestal.name = "base_pedestal"
base_pedestal.data.materials.append(wood_mat)
mod = base_pedestal.modifiers.new("Bevel", type='BEVEL')
mod.width = 0.003
mod.segments = 2

# COMPONENT: central_stem
stem_bottom = base_plate_h + base_pedestal_h
stem_center = stem_bottom + central_stem_h/2
bpy.ops.mesh.primitive_cylinder_add(radius=0.015, depth=central_stem_h, location=(0, 0, stem_center))
central_stem = bpy.context.active_object
central_stem.name = "central_stem"
central_stem.data.materials.append(brass_mat)

# COMPONENT: socket_cup
socket_z = 0.23  # ~0.20 from bottom as specified
bpy.ops.mesh.primitive_cylinder_add(radius=0.02, depth=0.06, location=(0, 0, socket_z))
socket_cup = bpy.context.active_object
socket_cup.name = "socket_cup"
socket_cup.data.materials.append(brass_mat)

# COMPONENT: light_bulb
bpy.ops.mesh.primitive_uv_sphere_add(radius=0.03, location=(0, 0, socket_z))
light_bulb = bpy.context.active_object
light_bulb.name = "light_bulb"
light_bulb.data.materials.append(bulb_mat)

# COMPONENT: shade_support
# Positioned so shade sits in upper portion of stem
shade_support_z = stem_center + central_stem_h/2 - 0.20 - 0.02 - 0.005  # 0.265
bpy.ops.mesh.primitive_cylinder_add(radius=0.03, depth=0.02, location=(0, 0, shade_support_z))
shade_support = bpy.context.active_object
shade_support.name = "shade_support"
shade_support.data.materials.append(brass_mat)

# COMPONENT: shade_bottom_ring
ring_minor_r = 0.005
bottom_ring_z = shade_support_z + 0.01 + ring_minor_r  # 0.28
bpy.ops.mesh.primitive_torus_add(location=(0, 0, bottom_ring_z), major_radius=0.125, minor_radius=ring_minor_r, major_segments=64, minor_segments=12)
shade_bottom_ring = bpy.context.active_object
shade_bottom_ring.name = "shade_bottom_ring"
shade_bottom_ring.data.materials.append(brass_mat)

# COMPONENT: shade_top_ring
shade_h = 0.20
top_ring_z = bottom_ring_z + shade_h  # 0.48
bpy.ops.mesh.primitive_torus_add(location=(0, 0, top_ring_z), major_radius=0.075, minor_radius=ring_minor_r, major_segments=64, minor_segments=12)
shade_top_ring = bpy.context.active_object
shade_top_ring.name = "shade_top_ring"
shade_top_ring.data.materials.append(brass_mat)

# COMPONENT: shade_fabric
fabric_center = (bottom_ring_z + top_ring_z) / 2  # 0.38
bpy.ops.mesh.primitive_cone_add(radius1=0.125, radius2=0.075, depth=shade_h, location=(0, 0, fabric_center))
shade_fabric = bpy.context.active_object
shade_fabric.name = "shade_fabric"
shade_fabric.data.materials.append(paper_mat)

# COMPONENT: finial
finial_z = stem_bottom + central_stem_h + finial_r  # 0.50
bpy.ops.mesh.primitive_uv_sphere_add(radius=finial_r, location=(0, 0, finial_z))
finial = bpy.context.active_object
finial.name = "finial"
finial.data.materials.append(jade_mat)

# CONNECTIONS
CONNECTIONS = [
    ("base_plate", "base_pedestal"),
    ("base_pedestal", "central_stem"),
    ("central_stem", "socket_cup"),
    ("socket_cup", "light_bulb"),
    ("central_stem", "shade_support"),
    ("shade_support", "shade_bottom_ring"),
    ("shade_bottom_ring", "shade_fabric"),
    ("shade_fabric", "shade_top_ring"),
    ("central_stem", "shade_top_ring"),
    ("shade_top_ring", "finial"),
    ("central_stem", "finial")
]