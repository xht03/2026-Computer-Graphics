import bpy
import math

# Clear scene
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

# Create wood material
wood_mat = bpy.data.materials.new(name="Wood")
wood_mat.use_nodes = True
nodes = wood_mat.node_tree.nodes
nodes.clear()
bsdf = nodes.new('ShaderNodeBsdfPrincipled')
output = nodes.new('ShaderNodeOutputMaterial')
wood_mat.node_tree.links.new(bsdf.outputs[0], output.inputs[0])
bsdf.inputs['Base Color'].default_value = (0.45, 0.28, 0.12, 1.0)
bsdf.inputs['Roughness'].default_value = 0.7
bsdf.inputs['Specular'].default_value = 0.3

# COMPONENT: top_board
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0, 0.77))
top_board = bpy.context.active_object
top_board.name = "top_board"
top_board.scale = (0.9, 0.9, 0.02)
top_board.data.materials.append(wood_mat)

# COMPONENT: apron (constructed as four boards joined into one object)
# Front
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, 0.43, 0.61))
apron_front = bpy.context.active_object
apron_front.scale = (0.9, 0.02, 0.3)
# Back
bpy.ops.mesh.primitive_cube_add(size=1, location=(0, -0.43, 0.61))
apron_back = bpy.context.active_object
apron_back.scale = (0.9, 0.02, 0.3)
# Right
bpy.ops.mesh.primitive_cube_add(size=1, location=(0.43, 0, 0.61))
apron_right = bpy.context.active_object
apron_right.scale = (0.02, 0.9, 0.3)
# Left
bpy.ops.mesh.primitive_cube_add(size=1, location=(-0.43, 0, 0.61))
apron_left = bpy.context.active_object
apron_left.scale = (0.02, 0.9, 0.3)

# Join apron parts
bpy.ops.object.select_all(action='DESELECT')
apron_front.select_set(True)
apron_back.select_set(True)
apron_right.select_set(True)
apron_left.select_set(True)
bpy.context.view_layer.objects.active = apron_front
bpy.ops.object.join()
apron = bpy.context.active_object
apron.name = "apron"
apron.data.materials.append(wood_mat)

# COMPONENT: leg (x4)
leg_positions = [
    (0.41, 0.41),
    (0.41, -0.41),
    (-0.41, 0.41),
    (-0.41, -0.41)
]

for i, pos in enumerate(leg_positions):
    # Create cylinder with slight flare at base for traditional curve
    bpy.ops.mesh.primitive_cylinder_add(radius=0.04, depth=0.7, location=(pos[0], pos[1], 0.35), vertices=24)
    leg = bpy.context.active_object
    leg.name = f"leg_{i+1}"
    
    # Modify bottom vertices to create curve (horse hoof style)
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action='DESELECT')
    bpy.ops.object.mode_set(mode='OBJECT')
    
    # Select bottom vertices (local z < -0.3)
    for vert in leg.data.vertices:
        if vert.co.z < -0.3:
            vert.select = True
    
    bpy.ops.object.mode_set(mode='EDIT')
    # Scale bottom outward for curve
    bpy.ops.transform.resize(value=(1.25, 1.25, 1.0), constraint_axis=(True, True, False))
    bpy.ops.object.mode_set(mode='OBJECT')
    
    leg.data.materials.append(wood_mat)

# CONNECTIONS list for assembly verification
CONNECTIONS = [
    ("leg_1", "apron"),
    ("leg_2", "apron"),
    ("leg_3", "apron"),
    ("leg_4", "apron"),
    ("apron", "top_board")
]