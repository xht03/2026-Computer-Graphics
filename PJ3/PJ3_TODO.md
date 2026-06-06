# 复旦图形学 PJ3 任务清单
## 基于 LL3M 范式的多创新扩展:面向中式家具的开源 3D 代码生成 Agent

> 本文档是一份完整的项目执行简报,可直接交给 Coding Agent(如 Claude Code、Codex、Cursor Agent 等)按计划执行。

---

## 一、项目背景

- **课程**:复旦大学 2026 春《计算机图形学》Project 3
- **主题方向**:生成式 AI Agent 在计算机图形学中的应用
- **截止日期**:**2026 年 6 月 10 日 23:59**(提交到 elearning)
- **汇报时间**:6 月 11 日(第 15 周)或 6 月 18 日(第 16 周),每组约 10 分钟
- **团队**:单人独立完成
- **资源**:有本地 GPU 可用

### 评分规则(总分 20)
- 基础分 10:按时提交报告 PDF + 汇报 PPT,按时汇报
- 质量分 10:报告创新性 + 汇报展示综合评分
- **满分条件**:需包含"**原创性、有效的创新**"

### 提交物
1. 英文 PDF 报告(推荐 ICLR 模板)
2. 汇报 PPT
3. 两者打包为 zip,命名:`2026图形学Project3 [本人姓名].zip`

---

## 二、选题与定位

### 项目名称(英文)
**Open-Source Multi-Agent System for Code-Based 3D Generation of Chinese Traditional Furniture**

### 核心思路
扩展 LL3M(Lu et al., 2025)的"Code-as-Generator"范式,做三方面有针对性的改进:**程序化几何验证**、**中文领域专精**、**全开源 LLM 量化评测**。让多个 LLM Agent 协作编写 Blender Python (`bpy`) 代码生成 3D 模型,所有产物可读、可编辑、可文化定制。

### 与 LL3M 的明确区分(写报告必备)

| 维度 | LL3M(已有) | 本工作(差异化) |
|---|---|---|
| Critic 机制 | VLM 视觉 + 代码自我批评 | **+ 程序化几何验证器(α)** |
| 输入语言/领域 | 英文,通用物体 | **中文,中式家具专精(β)** |
| LLM 选择 | 大概率 GPT-4 等闭源 | **全开源(GLM-4 + Qwen2-VL)(γ)** |
| 评测方式 | 定性 case 展示 | **+ Eval3DAIGC-198 量化基线(γ)** |

---

## 三、系统架构

```
[中文用户描述]
    ↓
[规划 Agent] —— 拆解任务:子部件列表 + 空间关系
    ↓
[编码 Agent + 领域 RAG] —— 检索中式家具示例 → 生成 bpy 代码
    ↓
[Blender 执行环境] —— 运行代码 + 渲染 4 视角(前/侧/顶/45°)
    ↓
[评论 Agent (VLM)]     [几何验证 Agent (程序化)]
        ↓                       ↓
        └────── 合并问题列表 ─────┘
                    ↓
              [修复 Agent] —— 定位代码、修改参数
                    ↓ (循环 2-3 轮)
[输出:3D 模型 (.blend/.obj) + 可读代码 + 渲染图]
```

### 各 Agent 的职责
| Agent | 输入 | 输出 | 实现方式 |
|---|---|---|---|
| 规划 | 中文用户描述 | JSON 子任务列表 | LLM + 结构化输出 prompt |
| 编码 | 子任务 + RAG 检索结果 | 完整 bpy 代码(带组件注释) | LLM + few-shot 模板 |
| **VLM Critic** | 4 张渲染图 + 用户描述 | 语义/视觉问题列表 | VLM |
| **几何验证 Agent** | Blender mesh 数据 | 几何问题列表(JSON) | 程序化 bpy 检查脚本 |
| 修复 | 合并问题 + 当前代码 | 修订代码 | LLM diff 改写 |

---

## 四、三个创新点(满分关键)

### 创新 α:程序化几何验证器 ⭐⭐⭐⭐(技术贡献)

**问题**:LL3M 的 critic 靠"看图 + 看代码",**漏掉几何细微缺陷**:
- 部件间未真正连接(看着连了,mesh 没合并)
- Non-manifold 边(无法 3D 打印/仿真)
- 部件穿模/碰撞
- 比例失衡(扶手低于坐板等结构不合理)

**方案**:用 `bpy` API 直接读 mesh 做程序化检查:
- `mesh.is_manifold` 检查
- bounding box 重叠检测(判断穿模 / 未连接)
- 距离阈值检查(部件间距是否合理)
- 比例约束(如:椅腿长度 / 椅高 ∈ [0.45, 0.55])

**输出**:JSON 问题列表,与 VLM critic 的输出合并送给修复 Agent。

**实现复杂度**:中。重点是 bpy 几何 API 调用 + 阈值调优。

### 创新 β:中文 + 中式家具领域专精 ⭐⭐⭐(应用贡献)

**问题**:LL3M 用英文 prompt 做通用物体,**缺乏文化与领域专精**。

**方案**:三个组件
1. **中文输入支持**:prompt 模板支持中文,planner 先做中→结构化 JSON 翻译
2. **领域 RAG 库**(12-15 条 bpy 片段):
   - 简化榫卯结构示意
   - 弧形椅圈构造函数
   - 中式回纹/卷云纹图案 procedural texture
   - 木质材质 shader 节点
   - 明式家具典型比例参数
3. **目标物体清单**(由易到难):
   - 简单:长凳、案几、方桌、屏风(几何规整)
   - 中等:八仙桌、太师椅、官帽椅
   - 困难:明式圈椅、罗汉床(留作"limitation discussion")

**实现复杂度**:中-低。RAG 用简单 embedding 检索(BGE-small)即可。

### 创新 γ:全开源 LLM + 量化评测 ⭐⭐⭐(方法学贡献)

**问题**:LL3M 大概率用 GPT-4 等闭源模型,且只做定性展示。

**方案**:
1. **全开源栈**:LLM 用 GLM-4-Flash(免费 API),VLM 用本地 Qwen2-VL-7B
2. **量化基线**:在 Eval3DAIGC-198(子集 30-50 个 prompt)上跑系统,报告:
   - 代码一次执行成功率
   - 平均迭代轮数(收敛速度)
   - 最终 CLIP score(图像与文本对齐度)
   - 几何验证通过率(α 创新的直接证明)
3. **消融实验**:对比"有 α / 无 α"、"有 β RAG / 无 β RAG"、"有迭代 / 无迭代"

**实现复杂度**:低-中。主要是跑实验和收集数据。

---

## 五、技术栈

### 核心组件
- **3D 引擎**:Blender 4.x(命令行无头模式:`blender --background --python script.py`)
- **Python 库**:`bpy`、`Pillow`、`numpy`、`requests`、`sentence-transformers`
- **LLM**:智谱 GLM-4-Flash(免费、中文友好);备选 Gemini 2.5 Flash AI Studio
- **VLM**:本地 Qwen2-VL-7B(用户有 GPU);备选 Gemini 2.5 Flash
- **RAG**:`bge-small-zh-v1.5` 做中文 embedding + FAISS / 简单 cosine 检索
- **编排**:自写 Python 状态机
- **演示界面**:Gradio

### 项目目录结构
```
pj3_blender_agent/
├── agents/
│   ├── planner.py        # 规划 Agent(中文输入)
│   ├── coder.py          # 编码 Agent(集成 RAG)
│   ├── vlm_critic.py     # VLM 视觉评论
│   ├── geom_verifier.py  # 程序化几何验证(α)
│   ├── fixer.py          # 修复 Agent
│   └── rag.py            # 领域 RAG 检索(β)
├── blender/
│   ├── runner.py         # Blender 子进程封装
│   └── render.py         # 多视角渲染
├── examples/             # 中式家具 RAG 示例库(12-15 条 bpy 片段)
├── prompts/              # 各 Agent 的 prompt 模板(中文)
├── eval/
│   ├── eval3daigc.py     # Eval3DAIGC-198 评测脚本
│   └── metrics.py        # CLIP score、几何通过率等
├── outputs/              # 生成的 .blend / 渲染图 / 代码
├── demo_app.py           # Gradio demo
├── main.py               # 主流程入口
├── requirements.txt
└── README.md
```

---

## 六、7 天执行计划(6 月 4 日 — 6 月 10 日)

### Day 1(6/4)环境搭建 + 最小可行流程
- [ ] 安装 Blender 4.x、配置无头模式调用
- [ ] 测试简单 `bpy` 脚本:生成立方体并渲染
- [ ] 申请 GLM-4 API key、测试调用
- [ ] 本地部署 Qwen2-VL-7B,测试图像理解
- [ ] 写 `main.py` 骨架:文本 → LLM → bpy 代码 → Blender 执行 → 渲染一张图
- **里程碑**:输入"a red cube",输出红色立方体渲染图

### Day 2(6/5)规划 + 编码 Agent(中文支持)
- [ ] 实现 `planner.py`:**支持中文输入**,LLM 输出结构化 JSON 子任务
- [ ] 实现 `coder.py`:基于子任务生成 bpy 代码,**在代码中加 `# COMPONENT: xxx` 注释**(便于 fixer 定位)
- [ ] 多视角渲染:`render.py`,从前/侧/顶/45° 渲染 4 张图
- [ ] 在 3 个中文测试用例上跑通:"长凳"、"方桌"、"屏风"
- **里程碑**:中文输入家具描述能产出可执行代码 + 4 视角渲染图

### Day 3(6/6)VLM Critic + Fixer 闭环
- [ ] 实现 `vlm_critic.py`:Qwen2-VL 输入 4 张图 + 描述,输出 JSON 问题列表
- [ ] 实现 `fixer.py`:接收问题 + 当前代码,LLM 输出修订代码
- [ ] 在 `main.py` 接入迭代循环(上限 3 轮)
- [ ] 测试 3 个 case,记录每轮改进
- **里程碑**:能看到第 2 轮渲染明显比第 1 轮好

### Day 4(6/7)**创新 α:程序化几何验证器**
- [ ] 实现 `geom_verifier.py`,集成 5 个核心检查:
  - Manifold 性检查
  - 部件间 bounding box 重叠 / 间隙检测
  - 整体比例约束(可配置规则)
  - Mesh 完整性(无悬空顶点)
  - 几何合理性(如重心稳定性,可选)
- [ ] 在 `main.py` 里把几何验证结果与 VLM 评论 **合并** 后送给 fixer
- [ ] 构造 3-5 个"VLM 看不出但几何有问题"的对比 case(消融实验素材)
- **里程碑**:能展示"VLM 漏掉、α 抓到"的具体案例

### Day 5(6/8)**创新 β:领域 RAG + 中式家具**
- [ ] 构建中式家具 bpy 示例库(`examples/`,约 12-15 条):
  - 至少包含:椅圈、桌腿、屏风框架、回纹纹理、木质 shader、典型比例参数
  - 每条带详细中文注释 + 元数据(适用物品类型)
- [ ] 实现 `rag.py`:用 `bge-small-zh-v1.5` 做 embedding,对编码 Agent 请求检索 top-3 相关片段
- [ ] 把检索结果作为 few-shot 加入编码 Agent 的 prompt
- [ ] 测试 5-6 个中式家具 case:长凳 / 案几 / 八仙桌 / 太师椅 / 屏风 / 简化版圈椅
- **里程碑**:中式家具的生成质量明显优于 baseline(无 RAG 的情况)

### Day 6(6/9)**创新 γ:量化评测 + 报告 + PPT**
- [ ] **上午**:Eval3DAIGC-198 评测
  - 从 HuggingFace 下载,挑 30 个适配你系统的 prompt
  - 跑全流程,收集指标:执行成功率、平均迭代轮数、CLIP score、几何通过率
- [ ] **上午**:消融实验(每组跑 5-8 个 case)
  - 完整系统 vs 无 α vs 无 β-RAG vs 无迭代
- [ ] **下午**:写英文报告初稿(ICLR 模板,详见第七节)
- [ ] **晚上**:做 PPT 初稿(详见第八节)
- **里程碑**:报告和 PPT 完整可用版本

### Day 7(6/10)缓冲 + 提交
- [ ] 上午:精修报告、润色英文表达、补图
- [ ] 下午:PPT 演练,卡时间在 10 分钟内
- [ ] **录制 demo 备份视频**(中式案几现场生成 + 修复 + RAG 编辑全过程)
- [ ] 整理 GitHub 仓库,完善 README,在报告中放链接
- [ ] **18:00 前提交 elearning**,zip 命名格式正确
- [ ] 提交后再确认一遍上传成功

---

## 七、报告结构(英文,6-8 页 ICLR 格式)

1. **Abstract**(150-200 词,明确列出三条贡献)
2. **Introduction**
   - 文生 3D 的可编辑性问题
   - LL3M 范式优势
   - 本工作三个具体改进
3. **Related Work**
   - 3.1 Text-to-3D Generation(引 Wen 2025 survey、Liu 2024 survey)
   - 3.2 Code-as-Generator(引 LL3M、CAD-Assistant)
   - 3.3 Multi-Agent LLM for 3D(引 Worldcraft、Idea23D、SAGE)
4. **Method**
   - 4.1 System Overview(架构图)
   - 4.2 Planner & Coder with Domain RAG(**创新 β**)
   - 4.3 Hybrid Critic: VLM + Geometric Verifier(**创新 α**)
   - 4.4 Fixer with Component Localization
5. **Experiments**
   - 5.1 Implementation Details(**创新 γ**:GLM-4 + Qwen2-VL,开源栈)
   - 5.2 Quantitative Evaluation on Eval3DAIGC-198(**创新 γ**)
   - 5.3 Chinese Furniture Case Studies(**创新 β**)
   - 5.4 Ablation:VLM-only vs +Geometric Verifier(**创新 α**)
6. **Discussion**:失败案例(明式圈椅)、局限、未来方向
7. **Conclusion**
8. **References**
9. **附录**:代码仓库链接、Prompt 模板、RAG 示例库样例

---

## 八、PPT 结构(10 张,10 分钟)

| 页 | 内容 | 时长 |
|---|---|---|
| 1 | 标题页(姓名+学号) | 30s |
| 2 | 背景:文生 3D 不可编辑 → Code-as-Generator(LL3M 引入) | 1 min |
| 3 | 与 LL3M 区分的三个创新点(总览图) | 1 min |
| 4 | 系统架构(突出 α 几何验证 + β RAG 模块) | 1 min |
| 5 | **创新 α**:几何 critic 抓到 VLM 漏掉的缺陷(对比图) | 1.5 min |
| 6 | **创新 β**:中式家具 case 展示 + RAG 库样例 | 1.5 min |
| 7 | **创新 γ**:Eval3DAIGC-198 量化结果 + 消融实验表 | 1.5 min |
| 8 | **现场 Live Demo**(中式案几) | 2 min |
| 9 | 失败案例 + 局限讨论(显学术诚实) | 30s |
| 10 | 总结 + 致谢 + 仓库链接 | 30s |

---

## 九、现场 Demo 设计(三创新一次性展示)

1. **打开 Gradio 界面**,中文输入框
2. **输入**:"做一张明式案几"(展示 β 中文输入)
3. **系统跑**:规划 → RAG 检索"案几"相关示例(屏幕展示检索结果)→ 生成代码 → 渲染
4. **几何 critic 报错**:屏幕展示"右前腿与桌面有 2mm 间隙,VLM 未检测到"(展示 α)
5. **自动修复 → 重新渲染** → 最终结果
6. **邀请评委说**:"加上回纹装饰"(展示 β RAG)→ 系统检索回纹片段 → 修改代码 → 再渲染
7. **口播带过**:整套系统跑在 GLM-4 + Qwen2-VL 上,完全免费(γ)

**应急**:每个步骤都预录视频,任何环节卡住立刻切到视频继续讲。

---

## 十、风险与应急方案

| 风险 | 应对 |
|---|---|
| `bpy` 代码总报错 | 加 try-except 重试;给编码 Agent 投喂 Blender API 文档片段 |
| VLM 评论质量差 | 简化任务到家具类;靠 α 几何 critic 兜底 |
| 几何验证误报多 | 阈值调宽松些;允许"可疑"等级而非直接报错 |
| 中文 RAG 检索不准 | 预先人工标注示例的关键词,加 keyword 检索作为补充 |
| LLM API 限流 | 准备 GLM + Gemini 双备选,自动切换 |
| 迭代不收敛 | 硬性上限 3 轮;给修复 Agent 加"最小改动"指令 |
| 现场 demo 翻车 | 预录视频备份;只演示已验证过的指令集 |
| 时间不够 | **优先级:α > β > γ**。γ 可降级为小样本评测;β 的难物体砍掉 |

---

## 十一、提交前最终 Checklist

- [ ] 英文报告 PDF(6-8 页,ICLR 模板,附代码链接)
- [ ] 汇报 PPT(10 页左右)
- [ ] 代码仓库 README 完整(安装、运行、demo 说明)
- [ ] 中式家具 RAG 示例库已整理好,放入仓库
- [ ] Eval3DAIGC-198 评测结果表格已加入报告
- [ ] Demo 备份视频已录制
- [ ] zip 命名格式:`2026图形学Project3 [本人姓名].zip`
- [ ] 在 6 月 10 日 23:59 之前提交到 elearning
- [ ] 准备好 6 月 11 日(或 18 日)现场汇报

---

## 十二、给执行 Agent 的特别说明

1. **本人可独立处理的部分**:申请 API key、提交 elearning、现场汇报、录 demo 视频
2. **需要 Agent 完成的核心部分**:全部代码实现、调试、跑实验、报告初稿、PPT 初稿
3. **关键设计决策(Agent 不得擅自简化)**:
   - 三个创新都必须有可见的、可演示的实现
   - 创新 α 必须有"VLM 漏掉、α 抓到"的具体对比 case
   - 创新 β 的 RAG 示例库不得少于 10 条
4. **代码风格**:模块化清晰,每个 Agent 一个文件,所有 prompt 放独立目录
5. **里程碑回报**:每天结束时,Agent 汇报当天完成项、未完成项、阻塞点

---

## 十三、关键参考文献(Related Work 必引)

1. **LL3M** — Lu et al., "LL3M: Large Language 3D Modelers", arXiv:2508.08228, 2025.(**直接对标**)·
2. **CAD-Assistant** — Mallis et al., ICCV 2025.(同范式工程方向)
3. **Worldcraft** — Liu et al., arXiv:2502.15601, 2025.(LLM agents 3D 定制)
4. **Idea23D** — Chen et al., COLING 2025.(多模态 Agent 3D 生成)
5. **3D scene generation: A survey** — Wen et al., arXiv:2505.05474, 2025.
6. **A comprehensive survey on 3D content generation** — Liu et al., arXiv:2402.01166, 2024.

**评测数据集**:
- **Eval3DAIGC-198** — `https://huggingface.co/yisuanwang/Idea23D`(创新 γ 评测用)

---

**文档版本**:v2(整合三创新)
**生成日期**:2026 年 6 月 4 日
**项目截止**:2026 年 6 月 10 日 23:59
