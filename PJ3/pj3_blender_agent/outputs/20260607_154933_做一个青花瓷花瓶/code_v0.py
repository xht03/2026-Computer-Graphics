import bpy
import math

# Clear scene
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

# Material helper
def make_material(name, color, roughness=0.05):
    mat = bpy.data.materials.new(name=name)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes["Principled BSDF"]
    bsdf.inputs["Base Color"].default_value = (*color, 1.0)
    bsdf.inputs["Roughness"].default_value = roughness
    return mat

white_mat = make_material("WhitePorcelain", (0.97, 0.97, 0.99), 0.05)
blue_mat = make_material("CobaltBlue", (0.08, 0.22, 0.55), 0.1)

# COMPONENT: foot_ring
bpy.ops.mesh.primitive_torus_add(major_radius=0.06, minor_radius=0.01, location=(0, 0, 0.01))
foot_ring = bpy.context.active_object
foot_ring.name = "foot_ring"
foot_ring.data.materials.append(white_mat)
bpy.ops.object.shade_smooth()

# COMPONENT: lower_body
bpy.ops.mesh.primitive_cone_add(radius1=0.05, radius2=0.10, depth=0.10, location=(0, 0, 0.07))
lower_body = bpy.context.active_object
lower_body.name = "lower_body"
lower_body.data.materials.append(white_mat)
bpy.ops.object.shade_smooth()

# COMPONENT: belly
bpy.ops.mesh.primitive_uv_sphere_add(radius=0.10, location=(0, 0, 0.12))
belly = bpy.context.active_object
belly.name = "belly"
belly.data.materials.append(white_mat)
bpy.ops.object.shade_smooth()

# COMPONENT: shoulder
bpy.ops.mesh.primitive_cone_add(radius1=0.10, radius2=0.03, depth=0.10, location=(0, 0, 0.27))
shoulder = bpy.context.active_object
shoulder.name = "shoulder"
shoulder.data.materials.append(white_mat)
bpy.ops.object.shade_smooth()

# COMPONENT: neck
bpy.ops.mesh.primitive_cylinder_add(radius=0.03, depth=0.10, location=(0, 0, 0.37))
neck = bpy.context.active_object
neck.name = "neck"
neck.data.materials.append(white_mat)
bpy.ops.object.shade_smooth()

# COMPONENT: rim
bpy.ops.mesh.primitive_torus_add(major_radius=0.05, minor_radius=0.01, location=(0, 0, 0.43))
rim = bpy.context.active_object
rim.name = "rim"
rim.data.materials.append(white_mat)
bpy.ops.object.shade_smooth()

# Decorative blue bands (underglaze decoration style)
# Base wave pattern
bpy.ops.mesh.primitive_torus_add(major_radius=0.051, minor_radius=0.003, location=(0, 0, 0.02))
base_decoration = bpy.context.active_object
base_decoration.name = "base_decoration"
base_decoration.data.materials.append(blue_mat)
bpy.ops.object.shade_smooth()

# Waist lotus band (junction of lower_body and belly)
bpy.ops.mesh.primitive_torus_add(major_radius=0.101, minor_radius=0.004, location=(0, 0, 0.12))
waist_band = bpy.context.active_object
waist_band.name = "waist_band"
waist_band.data.materials.append(blue_mat)
bpy.ops.object.shade_smooth()

# Cloud-collar at shoulder (junction of belly and shoulder)
bpy.ops.mesh.primitive_torus_add(major_radius=0.101, minor_radius=0.004, location=(0, 0, 0.22))
shoulder_band = bpy.context.active_object
shoulder_band.name = "shoulder_band"
shoulder_band.data.materials.append(blue_mat)
bpy.ops.object.shade_smooth()

# Neck bamboo/ruyi motifs
bpy.ops.mesh.primitive_torus_add(major_radius=0.031, minor_radius=0.002, location=(0, 0, 0.35))
neck_band1 = bpy.context.active_object
neck_band1.name = "neck_band1"
neck_band1.data.materials.append(blue_mat)
bpy.ops.object.shade_smooth()

bpy.ops.mesh.primitive_torus_add(major_radius=0.031, minor_radius=0.002, location=(0, 0, 0.39))
neck_band2 = bpy.context.active_object
neck_band2.name = "neck_band2"
neck_band2.data.materials.append(blue_mat)
bpy.ops.object.shade_smooth()

CONNECTIONS = [
    ("foot_ring", "lower_body"),
    ("lower_body", "belly"),
    ("belly", "shoulder"),
    ("shoulder", "neck"),
    ("neck", "rim")
]