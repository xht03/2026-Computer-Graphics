# TITLE: 木质材质 Shader 节点
# DESCRIPTION: 红木/花梨木质感材质，含纹理噪波+粗糙度，适用于所有家具部件
# TAGS: 材质,木材,红木,花梨木,shader,纹理,节点

import bpy

def make_wood_material(name: str = "RedWood",
                       base_color=(0.38, 0.18, 0.06, 1.0),
                       roughness: float = 0.55) -> bpy.types.Material:
    """创建逼真木质材质，包含随机噪波纹理模拟木纹。"""
    mat = bpy.data.materials.new(name=name)
    mat.use_nodes = True
    tree = mat.node_tree
    nodes = tree.nodes
    links = tree.links

    # 清除默认节点
    nodes.clear()

    # Principled BSDF
    bsdf = nodes.new("ShaderNodeBsdfPrincipled")
    bsdf.location = (400, 0)
    bsdf.inputs["Roughness"].default_value = roughness
    bsdf.inputs["Specular IOR Level"].default_value = 0.2

    # 噪波纹理（模拟木纹）
    noise = nodes.new("ShaderNodeTexNoise")
    noise.location = (0, 100)
    noise.inputs["Scale"].default_value = 18.0
    noise.inputs["Detail"].default_value = 8.0
    noise.inputs["Roughness"].default_value = 0.6
    noise.inputs["Distortion"].default_value = 0.3

    # 颜色渐变（将噪波映射到木纹颜色）
    ramp = nodes.new("ShaderNodeValToRGB")
    ramp.location = (200, 100)
    ramp.color_ramp.elements[0].color = (base_color[0] * 0.7, base_color[1] * 0.7, base_color[2] * 0.5, 1.0)
    ramp.color_ramp.elements[1].color = base_color

    # 连接
    links.new(noise.outputs["Fac"], ramp.inputs["Fac"])
    links.new(ramp.outputs["Color"], bsdf.inputs["Base Color"])

    # Output
    output = nodes.new("ShaderNodeOutputMaterial")
    output.location = (600, 0)
    links.new(bsdf.outputs["BSDF"], output.inputs["Surface"])

    return mat


# 示例用法：创建红木材质并赋给当前活动对象
wood_mat = make_wood_material("RedWood_Example")
# obj.data.materials.append(wood_mat)   # 注释掉示例调用，实际使用时取消注释
