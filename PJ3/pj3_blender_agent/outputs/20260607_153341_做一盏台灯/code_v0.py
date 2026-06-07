import bpy
import math

bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

# Materials
wood_mat = bpy.data.materials.new(name="DarkHardwood")
wood_mat.use_nodes = True
bsdf = wood_mat.node_tree.nodes["Principled BSDF"]
bsdf.inputs["Base Color"].default_value = (0.12, 0.04, 0.02, 1.0)
bsdf.inputs["Roughness"].default_value = 0.5

brass_mat = bpy.data.materials.new(name="Brass")
brass_mat.use_nodes = True
bsdf = brass_mat.node_tree.nodes["Principled BSDF"]
bsdf.inputs["Base Color"].default_value = (0.75, 0.55, 0.15, 1.0)
bsdf.inputs["Metallic"].default_value = 0.9
bsdf.inputs["Roughness"].default_value = 0.25

glass_mat = bpy.data.materials.new(name="Glass")
glass_mat.use_nodes = True
bsdf = glass_mat.node_tree.nodes["Principled BSDF"]
bsdf.inputs["Base Color"].default_value = (1.0, 1.0, 0.95, 1.0)
bsdf.inputs["Roughness"].default_value = 0.05
bsdf.inputs["IOR"].default_value = 1.45

fabric_mat = bpy.data.materials.new(name="RicePaper")
fabric_mat.use_nodes = True
bsdf = fabric_mat.node_tree.nodes["Principled BSDF"]
bsdf.inputs["Base Color"].default_value = (0.92, 0.88, 0.78, 1.0)
bsdf.inputs["Roughness"].default_value = 0.7

# COMPONENT: base
bpy.ops.mesh.primitive_cylinder_add(radius=0.10, depth=0.05, location=(0.0, 0.0, 0.025))
base = bpy.context.active_object
base.name = "base"
base.data.materials.append(wood_mat)

# COMPONENT: stem
bpy.ops.mesh.primitive_cylinder_add(radius=0.02, depth=0.28, location=(0.0, 0.0, 0.19))
stem = bpy.context.active_object
stem.name = "stem"
stem.data.materials.append(wood_mat)

# COMPONENT: socket_holder
bpy.ops.mesh.primitive_cylinder_add(radius=0.04, depth=0.02, location=(0.0, 0.0, 0.34))
socket_holder = bpy.context.active_object
socket_holder.name = "socket_holder"
socket_holder.data.materials.append(wood_mat)

# COMPONENT: socket
bpy.ops.mesh.primitive_cylinder_add(radius=0.015, depth=0.04, location=(0.0, 0.0, 0.37))
socket = bpy.context.active_object
socket.name = "socket"
socket.data.materials.append(brass_mat)

# COMPONENT: bulb
bpy.ops.mesh.primitive_uv_sphere_add(radius=0.025, location=(0.0, 0.0, 0.37))
bulb = bpy.context.active_object
bulb.name = "bulb"
bulb.data.materials.append(glass_mat)

# COMPONENT: shade_bottom_ring
bpy.ops.mesh.primitive_torus_add(major_radius=0.115, minor_radius=0.01, location=(0.0, 0.0, 0.36))
shade_bottom_ring = bpy.context.active_object
shade_bottom_ring.name = "shade_bottom_ring"
shade_bottom_ring.data.materials.append(wood_mat)

# COMPONENT: shade_top_ring
bpy.ops.mesh.primitive_torus_add(major_radius=0.08, minor_radius=0.01, location=(0.0, 0.0, 0.56))
shade_top_ring = bpy.context.active_object
shade_top_ring.name = "shade_top_ring"
shade_top_ring.data.materials.append(wood_mat)

# COMPONENT: shade_diffuser
bpy.ops.mesh.primitive_cone_add(radius1=0.125, radius2=0.09, depth=0.20, location=(0.0, 0.0, 0.46))
shade_diffuser = bpy.context.active_object
shade_diffuser.name = "shade_diffuser"
shade_diffuser.data.materials.append(fabric_mat)

# COMPONENT: finial
bpy.ops.mesh.primitive_uv_sphere_add(radius=0.02, location=(0.0, 0.0, 0.59))
finial = bpy.context.active_object
finial.name = "finial"
finial.data.materials.append(wood_mat)

CONNECTIONS = [
    ("base", "stem"),
    ("stem", "socket_holder"),
    ("socket_holder", "socket"),
    ("socket_holder", "shade_bottom_ring"),
    ("shade_bottom_ring", "shade_diffuser"),
    ("shade_diffuser", "shade_top_ring"),
    ("shade_top_ring", "finial"),
    ("socket", "bulb")
]