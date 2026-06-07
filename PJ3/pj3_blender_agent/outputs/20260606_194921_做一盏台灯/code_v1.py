import bpy
import math

# Clear scene
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

# Materials
# Dark wood for base and stem
wood_mat = bpy.data.materials.new(name="DarkWood")
wood_mat.use_nodes = True
bsdf = wood_mat.node_tree.nodes["Principled BSDF"]
bsdf.inputs["Base Color"].default_value = (0.35, 0.15, 0.08, 1.0)
bsdf.inputs["Roughness"].default_value = 0.5

# Rice paper for shade with subtle pattern suggestion via geometry
paper_mat = bpy.data.materials.new(name="RicePaper")
paper_mat.use_nodes = True
bsdf_p = paper_mat.node_tree.nodes["Principled BSDF"]
bsdf_p.inputs["Base Color"].default_value = (0.92, 0.88, 0.75, 1.0)
bsdf_p.inputs["Roughness"].default_value = 0.8

# Ceramic material for socket
ceramic_mat = bpy.data.materials.new(name="Ceramic")
ceramic_mat.use_nodes = True
bsdf_cer = ceramic_mat.node_tree.nodes["Principled BSDF"]
bsdf_cer.inputs["Base Color"].default_value = (0.9, 0.9, 0.85, 1.0)
bsdf_cer.inputs["Roughness"].default_value = 0.3

# Glass material for bulb
bulb_mat = bpy.data.materials.new(name="BulbGlass")
bulb_mat.use_nodes = True
bsdf_bulb = bulb_mat.node_tree.nodes["Principled BSDF"]
bsdf_bulb.inputs["Base Color"].default_value = (0.98, 0.98, 0.9, 1.0)
bsdf_bulb.inputs["Roughness"].default_value = 0.1
bsdf_bulb.inputs["Transmission"].default_value = 0.9

# Dimensions (total height 0.5m)
total_h = 0.5
base_h = 0.04
base_r = 0.14
shade_h = 0.22
shade_r = 0.125
stem_h = total_h - base_h - shade_h  # 0.24m
stem_r = 0.025

# COMPONENT: base
bpy.ops.mesh.primitive_cylinder_add(radius=base_r, depth=base_h, location=(0, 0, base_h/2))
base = bpy.context.active_object
base.name = "base"
base.data.materials.append(wood_mat)

# Add intricate carvings to base (decorative rings and motifs)
# Carved bands
for i in range(3):
    z_pos = 0.008 + i * 0.012
    bpy.ops.mesh.primitive_torus_add(major_radius=base_r + 0.001, minor_radius=0.004, location=(0, 0, z_pos))
    ring = bpy.context.active_object
    ring.name = f"base_carving_band_{i}"
    ring.data.materials.append(wood_mat)

# Chinese motif details (spherical decorations at cardinal directions)
for angle in [0, math.pi/2, math.pi, 3*math.pi/2]:
    x = math.cos(angle) * (base_r - 0.02)
    y = math.sin(angle) * (base_r - 0.02)
    bpy.ops.mesh.primitive_uv_sphere_add(radius=0.012, location=(x, y, base_h - 0.01))
    motif = bpy.context.active_object
    motif.name = f"base_motif_{int(angle*2)}"
    motif.data.materials.append(wood_mat)

# COMPONENT: stem
stem_z = base_h + stem_h/2
bpy.ops.mesh.primitive_cylinder_add(radius=stem_r, depth=stem_h, location=(0, 0, stem_z))
stem = bpy.context.active_object
stem.name = "stem"
stem.data.materials.append(wood_mat)

# COMPONENT: shade (open-ended tube)
shade_z = base_h + stem_h + shade_h/2
bpy.ops.mesh.primitive_cylinder_add(radius=shade_r, depth=shade_h, location=(0, 0, shade_z), end_fill_type='NOTHING')
shade = bpy.context.active_object
shade.name = "shade"
shade.data.materials.append(paper_mat)

# Add Chinese pattern frame to shade (wooden trim and vertical elements)
# Top and bottom wooden rings
for z_offset in [-shade_h/2 + 0.01, shade_h/2 - 0.01]:
    bpy.ops.mesh.primitive_torus_add(major_radius=shade_r + 0.003, minor_radius=0.006, location=(0, 0, shade_z + z_offset))
    trim = bpy.context.active_object
    trim.name = f"shade_frame_{'bottom' if z_offset < 0 else 'top'}"
    trim.data.materials.append(wood_mat)

# Vertical frame elements (suggesting traditional lattice pattern)
num_slats = 8
for i in range(num_slats):
    angle = (2 * math.pi * i) / num_slats
    x = math.cos(angle) * (shade_r + 0.003)
    y = math.sin(angle) * (shade_r + 0.003)
    bpy.ops.mesh.primitive_cylinder_add(radius=0.004, depth=shade_h - 0.02, location=(x, y, shade_z))
    slat = bpy.context.active_object
    slat.name = f"shade_slat_{i}"
    slat.data.materials.append(wood_mat)

# COMPONENT: internal hardware (socket and bulb)
# Socket at top of stem
socket_h = 0.03
socket_r = 0.02
socket_z = base_h + stem_h + socket_h/2
bpy.ops.mesh.primitive_cylinder_add(radius=socket_r, depth=socket_h, location=(0, 0, socket_z))
socket_obj = bpy.context.active_object
socket_obj.name = "socket"
socket_obj.data.materials.append(ceramic_mat)

# Bulb above socket
bulb_r = 0.025
bulb_z = socket_z + socket_h/2 + bulb_r
bpy.ops.mesh.primitive_uv_sphere_add(radius=bulb_r, location=(0, 0, bulb_z))
bulb_obj = bpy.context.active_object
bulb_obj.name = "bulb"
bulb_obj.data.materials.append(bulb_mat)

# CONNECTIONS: List of physically touching parts
CONNECTIONS = [
    ("base", "stem"),
    ("stem", "shade"),
    ("stem", "socket"),
    ("socket", "bulb")
]