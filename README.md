# ReaxNet Lite

一个面向 LAMMPS ReaxFF 轨迹的反应分析核心库。它读取标准
`fix reax/c/bonds` 输出，识别每一帧中的分子、相邻帧之间的反应事件，生成
SMILES/InChI，并将重复的正逆反应聚合为反应网络。

2026-10-28 的解析与分子构建性能优化及逐事件等价验证见
[反应分析性能验证](docs/reaction_analysis_performance_2026-10-28.md)。

键连接数量、官能团拓扑候选数量及其时间折线图见
[逐帧数量分析](docs/bond_and_functional_group_analysis_2026-10-28.md)。

本项目参考了 ChemTraYzer 的处理思路，但代码是独立的 Python 3 实现，不依赖
ChemTraYzer 源码或数据库。当前版本定位为可运行、可测试、可继续接入 Web
服务的第一阶段 MVP。

## PRD 第一阶段算法入口

2026-09-17 新增独立入口 `python -m reaxnet.connectivity`，用于输入核查、
逐帧文件拓扑读取和分子实例输出。2026-09-22 已取消单／双阈值及元素对阈值，
所有入口直接采用 bonds 邻接关系，不按价态删边；保留原始 BO 与原子追踪依据。
输出 cutoff 仅为元数据，可未知；化学键类型/电子状态由后续结构模型解释。
变更和已确认的目标筛选边界见 [文件拓扑与目标匹配设计](docs/file_topology_and_target_matching_2026-09-22.md)。

实现规则、运行命令和真实数据小段验证见
[第一阶段算法实现与验证](docs/algorithm_stage1_2026-09-17.md)。
下文的 `python -m reaxnet` 是原有分析流程；新入口尚未接入其物种与事件算法。

## PRD 结构推断与物种注册入口

`python -m reaxnet.structures` 在第一阶段连接图后运行 CHO 结构候选搜索，
显式分配化学键类型、形式电荷和自由基电子，再检查并生成 SMILES/InChI。
必须指定片段总电荷范围与自由基电子数上限，例如 `--charges 0 --max-radicals 1`。

它输出结构候选、supported/ambiguous/unresolved 状态、物种清单、实例关联、
逐帧物种数量和 SQLite 查询索引。搜索达到上限时保留未决结果。
实现规则、示例及真实数据覆盖限制见
[结构推断实现与验证](docs/algorithm_stage2_2026-09-17.md)。
新结构结果尚未接入旧版事件/往返/网络算法。

## 已实现功能

- 逐帧解析标准 `# Timestep` 分隔的 ReaxFF bonds 文件。
- 直接读取文件邻接关系，不再使用 BO 阈值或价态/配位数限制删边。
- 校验双向邻接及 BO 一致性，保留原始 BO 证据。
- 通过真实键图的连通分量识别分子，包括孤立原子和自由基。
- 生成分子式、与原子编号无关的图指纹和局部键签名。
- 使用 RDKit 尽可能生成 canonical SMILES、InChI 和 InChIKey。
- 比较相邻帧的断键/成键，并将共享分子的变化合并为一个反应事件。
- 区分 association、dissociation、isomerization 和 exchange 事件。
- 使用原子 ID 与真实键状态识别两步和多步回穿闭环。
- 聚合正向次数、逆向次数、总次数和净次数，输出反应超图 JSON。

## 算法流程

```text
fix reax/c/bonds
        |
        v
文件邻接关系 + 输入一致性校验
        |
        v
真实键图 -> 连通分量 -> 分子物种
        |
        +----> RDKit -> SMILES / InChI / InChIKey
        |
        v
相邻帧断键/成键差分
        |
        v
共享分子合并 -> 回穿过滤 -> 反应网络聚合
```

## 快速运行

当前核心分析不要求第三方依赖。SMILES/InChI 是可选能力，需要 RDKit：

```powershell
python -m pip install -e ".[identifiers]"
```

直接分析仓库中的最小示例：

```powershell
python -m reaxnet `
  --bonds examples\ch4_dissociation.bonds `
  --elements "1:C,2:H" `
  --output-dir outputs\ch4_demo `
  --timestep-scale 0.001
```

也可以安装后使用命令行入口：

```powershell
reaxnet --bonds bonds.reaxc `
  --elements "1:C,2:H,3:O" `
  --output-dir outputs\run01 `
  --timestep-scale 0.0001 `
  --recrossing-window 0.01
```

`--timestep-scale` 表示一个 LAMMPS timestep 对应多少 ps。例如积分步长为
`0.1 fs` 时应填写 `0.0001`。

## Web 工作台与大文件

2026-10-27 已加入以候选反应事件为入口的浏览器：分页筛选、前后拓扑对照、
同号原子高亮、守恒证据和单事件 JSON 导出。可直接分析整个轨迹，无需先指定分子式。
实现范围、运行方式及真实小段验证见
[候选反应事件工作台](docs/reaction_event_workbench_2026-10-27.md)。

Web 服务使用同一套有界内存流式分析器。对于 4 GB 等本机轨迹，推荐直接填写
绝对路径，服务端会原地读取文件，不经过 HTTP 上传或临时复制：

```powershell
cd D:\smiletosvg
python -m pip install -e ".[identifiers]"
cd web
python -m pip install -e .
smiletosvg-web
```

浏览器打开 `http://127.0.0.1:8000`。需要跨机器提交时可切换为分块上传；每个
请求只包含 8 MiB 数据，支持断点续传，当前服务端上限为 16 GiB。分析过程中
只保留当前帧、回穿时间窗内的状态和去重物种/反应聚合信息，不会缓存全部轨迹帧。

元素类型必须提供映射；连接划分不再依赖元素价态表。CHO 结构推断的适用范围另见第二阶段说明。

## 输出文件

| 文件 | 内容 |
| --- | --- |
| `species_timeseries.csv` | 每帧各物种数量，以及 formula、SMILES、InChI、InChIKey |
| `reaction_events.csv` | 每次断键/成键事件、参与原子、反应方向与过滤状态 |
| `molecule_registry.csv` | 去重后的物种身份、示例原子、键签名和转换状态 |
| `reaction_network.json` | 物种节点、反应超边、正逆次数、净次数和时间范围 |
| `event_index.sqlite` | 流式入口生成的完整事件证据、原子映射和分页索引 |

标识符生成失败不会终止分析。失败原因会写入 `molecule_registry.csv` 的
`status` 和 `message` 字段，其他 CSV 和网络仍然会生成。

## 验证

运行单元测试：

```powershell
python -m unittest discover -s tests -v
```

当前测试覆盖：标准 bonds 解析、原始邻接保留、甲烷断键、SMILES/InChI、
回穿过滤、流式内存边界、乱序分块上传、本机路径 API 和完整输出契约。旧版曾使用 ChemTraYzer 2.1 附带的 101 帧
`CH4 + OH` 示例做过兼容性验证。以下为旧阈值版本的历史结果，尚未在新文件拓扑策略下复验：

```text
CH4 + [OH] <-> [CH3] + H2O
forward = 2, reverse = 1, net = 1
```

设置 `--recrossing-window 0.01` 后，一对窗口内的正逆事件会被标记为回穿并
过滤，保留 1 个净事件。

## 当前边界

- 多步回穿通过原子相交事件链和首尾真实键状态实现；复杂催化循环仍需要谨慎设置
  时间窗，避免把真实的快速可逆过程误判为回穿。
- SMILES/InChI 来自 ReaxFF 键连通和离散化键级，是 best-effort 结果。对于
  芳香体系、过渡金属和强电荷转移体系，建议结合坐标、形式电荷或结构优化验证。
- 图指纹用于稳定地区分连接关系，但不是严格的完整图同构证明；对外交换数据时
  应优先使用成功生成的 InChIKey。
- 当前未实现速率常数、Arrhenius 拟合、量化化学任务提交和 GML 绘图。

## 代码结构

```text
reaxnet/parser.py     bonds 文件解析与元素映射
reaxnet/chemistry.py  真实键筛选、分子识别与图指纹
reaxnet/identity.py   RDKit SMILES/InChI 转换与降级处理
reaxnet/analysis.py   反应事件、回穿过滤与网络聚合
reaxnet/outputs.py    CSV/JSON 输出契约
reaxnet/cli.py        命令行入口
```
