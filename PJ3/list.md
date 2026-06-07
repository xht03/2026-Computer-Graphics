# 待解决问题清单

## 原始想法

1. planner换模型，画台灯的时候明显组件太简单化了

2. vlm是不是和geom分工比较好，vlm不检查connection缝隙之类的，主要看构型，设计，风格，等几何检查看不出来的，感觉有些vlm的几何警告并不可靠

3. 这是我看最近一次生成台灯的案例里，第二次迭代中，vlm给出了组件间有间隙的警告，这个是真的有间隙吗？还是vlm的误报？因为没有看到geom的警告，所以我觉得可能是vlm的误报了

4. coder的prompt中约束其只用简单的几何体，是不是这样僵化了设计，有没有什么方案可以在可靠性和设计自由度之间找到平衡？

5. 能不能有那种3D可以转视角查看的方式，现在的四个视图还行，但是如果我想跟细致的检查，就没办法了

---

## 分析与方案（Claude）

### 1. Planner 换模型 + 组件太简单

**诊断**：Planner 用的是 `planner.py:13` 的 `glm-4-flash` —— 最便宜最弱的档。台灯只拆出 base/stem/shade 三件、全是 `cylinder`，过于贫瘠。原因两层：模型弱 + prompt 没鼓励它细分。

**方案**：
- 换 `glm-4-plus` 或用 Coder 同款 `kimi-k2.5`（已有 client，统一也好）。
- 改 `planner_system.txt`，要求对每个大组件再拆**子部件**（如灯罩 = 罩面 + 上下圈口 + 骨架），并允许 `shape` 字段超出 cylinder/cube（cone/torus/profile 等）。

**倾向**：模型升级 + prompt 鼓励层级化拆分一起做。

### 2 & 3. VLM 和 Geom 分工 + 那个"间隙"是不是误报

**误报问题：是误报。** Round 2 里 VLM 说"stem 由多段圆柱组成、段间有可见缝隙"，但实际代码里 stem 是**一根连续圆柱**（`radius=0.015, depth=0.4`），外面套了 3 个金色装饰环——视觉上像分段，其实是整根。同一轮 Geom 报 `0 issue`。这是 VLM 看图产生的几何幻觉。

**分工思路正确**：
- **Geom**（程序化、读真实坐标）→ 负责一切几何事实：连接缝隙、穿模、尺度、悬空。有 ground truth，零误报。
- **VLM**（看图）→ 只负责 Geom 看不到的：整体比例/重心观感、风格是否到位、设计是否单调、组件是否缺失、对称性、材质搭配。**明确禁止它报缝隙/穿模/接触类问题**。

**方案**：改写 `vlm_critic_system.txt`，加硬规则"不要评价部件间是否接触/缝隙/穿插——这些由几何模块负责；你只看造型、比例、风格、完整度"。

### 4. Coder 只许用简单几何体 → 设计僵化 vs 可靠性

**最值得改的地方。** 典型 bug：**灯罩永远是直筒圆柱，做不出"上窄下宽"的传统灯罩。**

根因：prompt 明说 "do not attempt a cone"、"never write radius1=/radius2="，而 `coder.py:26-27` 的 sanitizer 会把**任何**含 `radius1=`/`radius2=` 的行删掉。但 `primitive_cone_add(radius1=底, radius2=顶, depth=)` 是**完全合法**的 API——锥台正是灯罩的天然形状。这条正则本为 cylinder 不支持 radius1/2 而写，却误伤 cone，从结构上禁止了所有锥形。

**平衡方案**（可叠加）：
- **A. 精准放开 cone**：sanitizer 只在 `cylinder_add` 行里剥 radius1/2，放行 `cone_add`；prompt 改成"灯罩/锥形件用 `primitive_cone_add`"。低风险、立竿见影。
- **B. 安全地扩充工具箱**：允许 bevel/subdivision_surface modifier、`primitive_cone_add`、screw/lathe 做回转体（中式灯罩、花瓶、葫芦造型），不需 edit-mode、不易崩。
- **C. RAG 配方库**：把"锥形灯罩""带弧度回转体""镂空骨架"等做成已验证代码片段，检索注入。自由度交给配方，可靠性由"已验证"保证——最佳平衡点，但工作量最大。

**倾向**：先 A（零成本、直接解决灯罩），再 B 扩两三个安全构造，C 后续。

### 5. 可转视角的 3D 查看

**方案**：渲染后让 Blender 额外 `export_scene.gltf` 导出 `.glb`，再生成单文件 `view.html` 内嵌 Google 的 `<model-viewer>`（一个 `<script>` + 一个标签）。双击 html 即可在浏览器自由旋转/缩放，零依赖、零服务器。比转盘动画更适合细致检查。

---

## 落地顺序

3/2（VLM 分工，最快止血误报）→ 4A（放开 cone，解决灯罩）→ 1（Planner 升级）→ 5（glb 查看器）→ 4B/4C（扩工具箱）

> 注：实际按用户要求从 **问题 1（Planner）** 开始逐个修复。

---

## 修复进度

### ✅ 问题 1：Planner 换模型 + 组件太简单（已完成）

**效果**：组件数 3→8/9，形状从全 cylinder 扩展到 torus/cylinder/cone/sphere，连接数 2→11，装饰正确归入 style_notes，并补齐灯泡+灯座等功能件。

**改动**：
- `agents/planner.py`：ZhipuAI→OpenAI client；模型 `glm-4-flash`→`kimi-k2.5`；`max_tokens` 1024→16000（kimi 是推理模型，思考吃 token）；`temperature=1`（kimi 强制）；JSON 解析加“提取首个`{`到末个`}`”兜底；`plan_to_coder_description` 支持 `part_of`。
- `prompts/planner_system.txt`：schema 加 `cone/sphere/curve/custom` + 可选 `part_of`；新增 **DECOMPOSITION GRANULARITY**（6-12 个真实结构件、按形状选 cone/sphere/torus/curve、圆台=带两半径的 cone）；新增**功能件软规则**（灯具/电器需含灯泡 sphere + 灯座 cylinder，即使被遮挡）。

**衔接点**：planner 已把灯罩标为 `cone`，但 Coder 仍禁 cone（问题 4A 待解锁），故当前灯罩仍会退化成圆柱。

### ✅ 问题 2 & 3：VLM 与 Geom 分工，消除间隙误报（已完成）

**确认**：round 2 那条“stem 段间有缝隙”是 VLM 看图幻觉（实际是一根连续圆柱+装饰环），同轮 geom 报 0 issue。

**分工原则**：
- Geom（读真实坐标，零误报）→ 接触/缝隙/穿模/悬空/绝对尺度。
- VLM（看图）→ 只看几何检查器看不到的：缺件/多件、比例与视觉平衡、风格契合度、设计丰富度、对称性、材质配色、整体形态。**明令禁止报缝隙/接触/穿模/悬空/尺度。**

**改动**：重写 `prompts/vlm_critic_system.txt`（开头加 DIVISION OF LABOR，列❌禁报 / ✅应报清单）+ `prompts/vlm_critic_user.txt`（去掉“structural/geometric problems”措辞，改为聚焦设计）。代码 `vlm_critic.py` 无需改（仅读 prompt + 解析 JSON，GLM-4V 兜底共用同一 prompt）。

**实测**（对上次产生误报的渲染图）：旧 prompt 出 4 条含 3 条缝隙/接触/穿模误报；新 prompt 出 1 条比例/平衡观察，连跑 iter_1/iter_2 均**零几何事实泄漏**，且准确抓到 round 2 灯罩过宽的真实审美问题。

### ✅ 问题 4A：放开 cone，圆台灯罩生效（已完成）

**两道闸门**：① `coder.py` sanitizer 在任何行删 `radius1=/radius2=`；② `coder_system.txt` 明令禁 cone。

**改动**：
- `agents/coder.py`：sanitizer 的 radius1/2 规则改成只在 `cylinder_add` 行触发（cone 合法使用两半径做圆台，放行）。
- `prompts/coder_system.txt`：删掉 cone 禁令，新增 **Cone/frustum** 用法说明（`primitive_cone_add(radius1=底, radius2=顶, depth=)`，圆台用于灯罩/喇叭口底座/收腰腿/瓶身）。

**实测**：单元测试三种情况正确（cone 保留 / cylinder 保留 / cylinder 误用 radius1/2 拦截）。端到端跑 `做一盏台灯 --geom`：planner 出 9 件含 cone 灯罩 → coder 生成 `primitive_cone_add(radius1=0.125, radius2=0.09, ...)`，0 行被 sanitize，geom 0 问题，渲染出**标准圆台灯罩**（下宽上窄+顶钮），相比之前直筒鼓状灯罩质变。

**RAG 配方库（C）说明**：基础设施已存在（`examples/` 15 个手写片段 + `rag.py` BGE 检索，`--rag` 启用）。C = 往 examples/ 再写已验证的高级配方文件，工作量在“逐个写好并验证不崩”。本轮按用户要求跳过。

### ✅ 问题 4B：扩充工具箱 — 回转体 lathe + Bevel 倒角（已完成，跳过 subdivision）

**Bevel**：`coder_system.txt` 新增 Bevel modifier 用法（`obj.modifiers.new("Bevel", type='BEVEL')`），给硬边倒角，零风险（不碰任何 sanitizer 规则）。

**回转体 lathe（真 screw modifier）**：用户选了真 lathe 而非堆叠圆台，接受重开 `from_pydata`。
- `agents/coder.py`：移除 `from_pydata` 全面禁令（改为靠 prompt 约束 faces=[]）。
- `coder_system.txt`：新增 **Revolution bodies** recipe —— from_pydata 建轮廓线（verts+edges，faces 必须为 []）→ Screw modifier 旋转 → **depsgraph 烘焙成真实网格**（关键：不烘焙则 geom 读 bound_box 会因 y≈0 误判退化）。并强调“连续曲线体必须整体 lathe，禁止用圆台堆（留台阶）”。
- `planner_system.txt`：平滑曲线体（花瓶/葫芦/瓷身/收腰）是**一个** shape="curve" 组件，禁止切成 foot/belly/neck 多段。

**实测**：
- 机制单测：lathe 花瓶烘焙后尺寸 0.16×0.16×0.28（非退化）、258 真顶点、geom 零误报。
- 端到端 `做一个青花瓷花瓶`：修复前 planner 切 6 段→coder 堆圆台→**渲染有台阶**；加引导后 planner 出 1 个 `vase_body`→coder 用 Screw lathe（0 圆台）→**渲染平滑连续**，geom 零问题。

### 🔄 RAG 重构：bpy API 文档（LL3M BlenderRAG 路线）

**背景**：examples/ 是 AI 生成的成品样例（非优质库，可能含 bug）。参考 LL3M 的 BlenderRAG（arXiv 2508.08228）——知识库用**官方 API 文档**而非成品样例，且在编码前 + 报错时两处检索。用户决定：① 知识库换成 bpy API 文档；② **报错驱动检索（Fixer）优先做**。

**✅ 已完成：知识库构建 + 检索器 + Fixer 接入**
- `build_bpy_kb.py`：跑在 Blender 里**内省 RNA 元数据**（非 ML、非抓取，确定性提取），输出 `knowledge/bpy_api.jsonl`（561 条：493 operator + 55 modifier + 13 shader node），含精确参数签名/enum/属性/socket 名。**版本匹配本机 Blender 5.1**——抓到的 Principled socket 真名是 `Specular IOR Level`/`Transmission Weight`，与 prompt 里硬编码的“4.x”不同。
- `agents/api_rag.py` `ApiDocRetriever`：bge-small-en-v1.5 向量（磁盘缓存）+ 英文关键词兜底。`retrieve_for_error()` 专为 traceback 设计：提炼出错行+引号串+bpy路径+kwargs 做聚焦查询，关键词+语义混合排序，加停用词过滤。
- Fixer 接入：`fixer.py` 崩溃时按 error 文本检索 API 文档，注入 `fixer_user.txt` 的 `{api_docs}`；`fixer_system.txt` 加“信任 API REFERENCE”规则，并对齐 4A/4B（放行 cone、lathe 的 from_pydata faces=[]）。
- 三类崩溃检索验证：socket名错→Principled、cylinder误用radius1→cylinder+cone、modifier属性typo→ScrewModifier，均命中。端到端：fixer 把 `Specular`→`Specular IOR Level`、`Transmission`→`Transmission Weight`，用对 5.1 真名。
- `.gitignore`：忽略 embedding 二进制缓存（`*.npy`/`*.meta`，可重建）。

**✅ 已完成：Coder 接入（plan 驱动针对性检索）+ 删除旧 RAG**
- 关键洞察：plan 的 `shape` 字段已声明每个组件用什么构造，故检索**确定性驱动**——扫描计划用到的高级构造（cone/curve→lathe/torus/sphere、style_notes 提 bevel）才检索对应文档，cube/cylinder 不加噪；并**总是注入 Principled socket**（每个脚本都设材质，5.1 socket 名易错，主动预防）。
- `api_rag.py` 加 `retrieve_for_plan(plan)`；`coder.py` `generate(description, api_docs=)` 注入；`coder_user.txt` 占位符 `{rag_context}`→`{api_docs}`。
- `main.py`：删 `--rag`/`use_rag`，改为 plan 驱动 API-doc 检索（默认开）。
- **彻底删除** `examples/`（15 个 AI 生成样例）+ `agents/rag.py`（RAGRetriever）。
- 端到端 `做一盏台灯 --geom`：检索 3 条（sphere/torus/Principled）→ 生成代码全用安全 socket 名（Base Color/Metallic/Roughness，无 Specular/Transmission 崩溃）→ geom 0 问题、渲染成功。

**RAG 重构小结**：知识库从「AI 生成成品样例」彻底转为「内省本机 Blender 的 API 文档」（LL3M BlenderRAG 路线），两处接入——Coder 编码前（plan 驱动）+ Fixer 报错时（traceback 驱动）。

### ✅ VLM 拿回「宏观摆放错误」职责（已完成，修正问题 2/3 的矫枉过正）

**背景**：问题 2/3 把所有位置/结构都禁了，导致「灯泡露在灯罩外」这类**粗大摆放错误**没人管（geom 只查声明连接的缝隙，灯泡↔socket 相接即通过；VLM 被全面禁位置）。

**分界改为按"尺度"而非"是否碰位置"**：
- VLM 该报：粗大、一眼可见的**摆放/空间关系错误**（错位、错朝向、该藏的露在外/该包含的在外面）。
- VLM 不该报：细粒度接触（精确相接、毫米缝隙、微穿插、绝对尺度）——仍归 geom。

**改动**：`vlm_critic_system.txt`（DIVISION OF LABOR 重写为按尺度分工，✅清单加"GROSS placement/spatial-relationship errors"项）+ `vlm_critic_user.txt`。

**实测**：灯泡露出渲染——VLM 能准确报为 error（"bulb exposed below shade"），好模型仍 0 误报。但**检测率仅 1/3**，瓶颈是**渲染取景太小**（灯泡在画面里太小看不清）→ 指向下一步：渲染自动取景。

### ✅ 渲染自动取景（已完成）

**动机**：固定相机为 ~1m 家具标定，台灯/部件/花瓶都在画面里很小 → 害 VLM 检测、害人工检查（曾把正确八仙桌误看成断腿）。

**改动**：`blender/render.py` —— 渲染前算所有 MESH 的世界包围盒，取中心+半径，按相机 FOV（`cam_data.angle`）算距离 `radius/tan(fov/2)*1.25`，四视图改为「中心 + 单位方向×距离、look_at 中心」。**只动相机，灯光不变**（灯光在远处、范围大，缩放相机不影响照明，避开照明强度连锁问题）。空场景兜底默认距离。

**实测**：台灯填满画面，灯泡露出一目了然；VLM 灯泡检测 **1/3 → 2/3**，报错更精准。一改多利（VLM 检测↑、人工检查↑、问题5 查看器复用）。
