import bpy
import math

# Clear scene
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

# Material creation helper
def create_material(name, color, metallic=0.0, roughness=0.5):
    mat = bpy.data.materials.new(name=name)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes["Principled BSDF"]
    bsdf.inputs['Base Color'].default_value = (*color, 1.0)
    bsdf.inputs['Metallic'].default_value = metallic
    bsdf.inputs['Roughness'].default_value = roughness
    return mat

wood_mat = create_material("Wood", (0.4, 0.22, 0.08), metallic=0.0, roughness=0.65)
bamboo_mat = create_material("Bamboo", (0.88, 0.78, 0.45), metallic=0.0, roughness=0.75)
metal_mat = create_material("Metal", (0.12, 0.12, 0.12), metallic=0.85, roughness=0.35)

# COMPONENT: wooden_base
bpy.ops.mesh.primitive_cylinder_add(radius=0.15, depth=0.1, location=(0, 0, 0.05))
base = bpy.context.active_object
base.name = "wooden_base"
base.scale = (1.0, 0.6666667, 1.0)
bpy.ops.object.transform_apply(scale=True)
base.data.materials.append(wood_mat)

# COMPONENT: carving
carving_parts = []
num_details = 16
for i in range(num_details):
    angle = (2 * math.pi * i) / num_details
    x = 0.15 * math.cos(angle)
    y = 0.1 * math.sin(angle)
    bpy.ops.mesh.primitive_uv_sphere_add(radius=0.01, location=(x, y, 0.05))
    detail = bpy.context.active_object
    carving_parts.append(detail)

bpy.ops.object.select_all(action='DESELECT')
for obj in carving_parts:
    obj.select_set(True)
bpy.context.view_layer.objects.active = carving_parts[0]
bpy.ops.object.join()
carving = bpy.context.active_object
carving.name = "carving"
carving.data.materials.append(wood_mat)

# COMPONENT: bamboo_shade
# [sanitized] bpy.ops.mesh.primitive_cone_add(radius1=0.1, radius2=0.05, depth=0.3, location=(0, 0, 0.25))
shade = bpy.context.active_object
shade.name = "bamboo_shade"
shade.data.materials.append(bamboo_mat)

# COMPONENT: light_bulb_holder
bpy.ops.mesh.primitive_cylinder_add(radius=0.025, depth=0.1, location=(0, 0, 0.4))
holder = bpy.context.active_object
holder.name = "light_bulb_holder"
holder.data.materials.append(metal_mat)

CONNECTIONS = [
    ("wooden_base", "bamboo_shade"),
    ("wooden_base", "carving"),
    ("bamboo_shade", "light_bulb_holder")
]