import bpy

# Clear scene
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

# Material setup
def create_material(name, base_color, roughness=0.5, metallic=0.0):
    mat = bpy.data.materials.new(name=name)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes["Principled BSDF"]
    bsdf.inputs["Base Color"].default_value = (*base_color, 1.0)
    bsdf.inputs["Roughness"].default_value = roughness
    bsdf.inputs["Metallic"].default_value = metallic
    return mat

wood_mat = create_material("Wood", (0.45, 0.25, 0.1), roughness=0.7)
gold_mat = create_material("Gold", (0.85, 0.65, 0.15), roughness=0.3, metallic=0.9)
red_mat = create_material("Red", (0.7, 0.1, 0.1), roughness=0.6)

# COMPONENT: base
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, 0.02))
base = bpy.context.active_object
base.name = "base"
base.scale = (0.2, 0.2, 0.04)
base.data.materials.append(wood_mat)

# Base decorative studs (gold corners)
stud_positions = [(0.08, 0.08), (-0.08, 0.08), (-0.08, -0.08), (0.08, -0.08)]
for i, (x, y) in enumerate(stud_positions):
    bpy.ops.mesh.primitive_cube_add(size=1, location=(x, y, 0.05))
    stud = bpy.context.active_object
    stud.name = f"base_stud_{i+1}"
    stud.scale = (0.02, 0.02, 0.02)
    stud.data.materials.append(gold_mat)

# Base center inlay (red)
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, 0.0425))
base_inlay = bpy.context.active_object
base_inlay.name = "base_inlay"
base_inlay.scale = (0.08, 0.08, 0.005)
base_inlay.data.materials.append(red_mat)

# COMPONENT: stem
bpy.ops.mesh.primitive_cylinder_add(radius=0.015, depth=0.4, location=(0, 0, 0.24))
stem = bpy.context.active_object
stem.name = "stem"
stem.data.materials.append(wood_mat)

# Stem decorative rings (gold accents)
ring_heights = [0.14, 0.24, 0.34]
for i, z in enumerate(ring_heights):
    bpy.ops.mesh.primitive_cylinder_add(radius=0.018, depth=0.01, location=(0, 0, z))
    ring = bpy.context.active_object
    ring.name = f"stem_ring_{i+1}"
    ring.data.materials.append(gold_mat)

# COMPONENT: shade
bpy.ops.mesh.primitive_cylinder_add(radius=0.125, depth=0.22, location=(0, 0, 0.55))
shade = bpy.context.active_object
shade.name = "shade"
shade.data.materials.append(red_mat)

# Shade gold trim (top and bottom)
trim_heights = [0.44, 0.66]
for i, z in enumerate(trim_heights):
    bpy.ops.mesh.primitive_torus_add(location=(0, 0, z), major_radius=0.125, minor_radius=0.005)
    trim = bpy.context.active_object
    trim.name = f"shade_trim_{i+1}"
    trim.data.materials.append(gold_mat)

# Shade vertical gold patterns (4 cardinal directions)
pattern_offsets = [(0.13, 0), (0, 0.13), (-0.13, 0), (0, -0.13)]
for i, (x, y) in enumerate(pattern_offsets):
    bpy.ops.mesh.primitive_cube_add(size=1, location=(x, y, 0.55))
    pattern = bpy.context.active_object
    pattern.name = f"shade_pattern_{i+1}"
    pattern.scale = (0.01, 0.03, 0.11)
    pattern.data.materials.append(gold_mat)

# CONNECTIONS: structural assembly connections
CONNECTIONS = [
    ("base", "stem"),
    ("stem", "shade")
]