# TITLE: 回纹装饰（程序化贴图）
# DESCRIPTION: 用程序化节点生成中式回纹（雷纹/方格螺旋纹），贴在面板表面
# TAGS: 回纹,雷纹,装饰,纹样,程序化,贴图,纹理,材质,中式

import bpy


def make_thunder_pattern_material(name: str = "ThunderPattern",
                                  base_color=(0.35, 0.18, 0.06, 1.0),
                                  pattern_color=(0.85, 0.70, 0.40, 1.0),
                                  scale: float = 20.0) -> bpy.types.Material:
    """创建回纹（雷纹）程序化材质。

    使用方波噪波模拟方格螺旋纹：
    - 波形纹理 + Checker 纹理叠加
    - 两种颜色交替形成回纹效果
    """
    mat = bpy.data.materials.new(name=name)
    mat.use_nodes = True
    tree = mat.node_tree
    nodes = tree.nodes
    links = tree.links
    nodes.clear()

    # Output
    out = nodes.new("ShaderNodeOutputMaterial")
    out.location = (700, 0)

    bsdf = nodes.new("ShaderNodeBsdfPrincipled")
    bsdf.location = (500, 0)
    bsdf.inputs["Roughness"].default_value = 0.70
    links.new(bsdf.outputs["BSDF"], out.inputs["Surface"])

    # MixRGB：将两种颜色按 checker 混合
    mix = nodes.new("ShaderNodeMixRGB")
    mix.location = (300, 0)
    mix.inputs["Color1"].default_value = base_color
    mix.inputs["Color2"].default_value = pattern_color
    links.new(mix.outputs["Color"], bsdf.inputs["Base Color"])

    # Checker 纹理（主格）
    checker1 = nodes.new("ShaderNodeTexChecker")
    checker1.location = (0, 100)
    checker1.inputs["Scale"].default_value = scale
    links.new(checker1.outputs["Fac"], mix.inputs["Fac"])

    # Wave 纹理（叠加条纹，模拟回纹线条）
    wave = nodes.new("ShaderNodeTexWave")
    wave.location = (0, -100)
    wave.inputs["Scale"].default_value = scale * 4
    wave.inputs["Distortion"].default_value = 2.0
    wave.inputs["Detail"].default_value = 3.0
    wave.wave_type = "BANDS"

    # Math: multiply checker × wave → 产生回纹格栅效果
    math_node = nodes.new("ShaderNodeMath")
    math_node.operation = "MULTIPLY"
    math_node.location = (150, 0)
    links.new(checker1.outputs["Fac"], math_node.inputs[0])
    links.new(wave.outputs["Fac"],     math_node.inputs[1])
    links.new(math_node.outputs["Value"], mix.inputs["Fac"])

    # Texture Coordinate
    coord = nodes.new("ShaderNodeTexCoord")
    coord.location = (-200, 0)
    links.new(coord.outputs["UV"], checker1.inputs["Vector"])
    links.new(coord.outputs["UV"], wave.inputs["Vector"])

    return mat


# 用法示例：
thunder_mat = make_thunder_pattern_material("ThunderPattern_Example")
# 赋给任意面板或装饰条：
# obj.data.materials.append(thunder_mat)
