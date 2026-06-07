import bpy
import math

# Clear scene
bpy.ops.object.select_all(action='SELECT')
bpy.ops.object.delete(use_global=False)

# Materials
def make_material(name, color, roughness):
    mat = bpy.data.materials.new(name=name)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes["Principled BSDF"]
    bsdf.inputs["Base Color"].default_value = (*color, 1.0)
    bsdf.inputs["Roughness"].default_value = roughness
    return mat

mat_white = make_material("porcelain_white", (0.98, 0.98, 1.0), 0.05)
mat_blue = make_material("cobalt_blue", (0.05, 0.15, 0.5), 0.2)
mat_clay = make_material("unglazed_clay", (0.6, 0.4, 0.25), 0.9)

# COMPONENT: vase_body
# Profile for lathe (radius, z)
profile = [
    (0.0, 0.02),      # Center at foot top
    (0.035, 0.02),    # Foot top radius (inside foot ring)
    (0.11, 0.10),     # Belly lower (max radius 0.11)
    (0.11, 0.26),     # Belly upper
    (0.04, 0.36),     # Neck
    (0.04, 0.38),     # Rim top
    (0.0, 0.38),      # Close top
]
verts = [(r, 0.0, z) for (r, z) in profile]
edges = [(i, i+1) for i in range(len(verts)-1)]
me_body = bpy.data.meshes.new("vase_body_mesh")
obj_body = bpy.data.objects.new("vase_body", me_body)
bpy.context.collection.objects.link(obj_body)
me_body.from_pydata(verts, edges, [])
me_body.update()

screw = obj_body.modifiers.new("Screw", type='SCREW')
screw.axis = 'Z'
screw.angle = math.radians(360)
screw.steps = 64
screw.use_merge_vertices = True

# Bake modifier
dg = bpy.context.evaluated_depsgraph_get()
obj_body.data = bpy.data.meshes.new_from_object(obj_body.evaluated_get(dg))
obj_body.modifiers.clear()
obj_body.data.materials.append(mat_white)

# COMPONENT: foot_ring
bpy.ops.mesh.primitive_torus_add(
    major_radius=0.05,
    minor_radius=0.01,
    major_segments=48,
    minor_segments=16,
    location=(0, 0, 0.01)
)
obj_foot = bpy.context.active_object
obj_foot.name = "foot_ring"
obj_foot.data.materials.append(mat_clay)

# COMPONENT: rim_ring
bpy.ops.mesh.primitive_torus_add(
    major_radius=0.04,
    minor_radius=0.0075,
    major_segments=48,
    minor_segments=16,
    location=(0, 0, 0.38)
)
obj_rim = bpy.context.active_object
obj_rim.name = "rim_ring"
obj_rim.data.materials.append(mat_blue)

# COMPONENT: left_handle
cu_left = bpy.data.curves.new("left_handle_curve", type='CURVE')
cu_left.dimensions = '3D'
cu_left.bevel_depth = 0.01
cu_left.resolution_u = 12
sp_left = cu_left.splines.new('BEZIER')
sp_left.bezier_points.add(2)  # 3 points total
pts_left = [(-0.10, 0, 0.20), (-0.13, 0, 0.245), (-0.04, 0, 0.29)]
for bp, co in zip(sp_left.bezier_points, pts_left):
    bp.co = co
    bp.handle_left_type = 'AUTO'
    bp.handle_right_type = 'AUTO'
obj_left = bpy.data.objects.new("left_handle", cu_left)
bpy.context.collection.objects.link(obj_left)
obj_left.data.materials.append(mat_blue)

# COMPONENT: right_handle
cu_right = bpy.data.curves.new("right_handle_curve", type='CURVE')
cu_right.dimensions = '3D'
cu_right.bevel_depth = 0.01
cu_right.resolution_u = 12
sp_right = cu_right.splines.new('BEZIER')
sp_right.bezier_points.add(2)
pts_right = [(0.10, 0, 0.20), (0.13, 0, 0.245), (0.04, 0, 0.29)]
for bp, co in zip(sp_right.bezier_points, pts_right):
    bp.co = co
    bp.handle_left_type = 'AUTO'
    bp.handle_right_type = 'AUTO'
obj_right = bpy.data.objects.new("right_handle", cu_right)
bpy.context.collection.objects.link(obj_right)
obj_right.data.materials.append(mat_blue)

# Required connections
CONNECTIONS = [
    ("vase_body", "foot_ring"),
    ("vase_body", "rim_ring"),
    ("vase_body", "left_handle"),
    ("vase_body", "right_handle"),
]