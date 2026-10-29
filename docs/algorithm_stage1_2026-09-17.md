# 第一阶段算法实现与验证

日期：2026-09-17；连接规则于 2026-09-22 更新。对应 PRD F01–F04。
当前规则以用户 2026-09-22 决定为准：直接使用文件邻接关系，删除单／双阈值判断。

这次交付将连续 BO 的连接判定与化学键类型、价态及物种身份推断分开。
新入口是 `python -m reaxnet.connectivity`；原有 `python -m reaxnet` 仍使用旧版完整分析流程。
本阶段不生成 SMILES、species_id、反应或路径，不能将其输出称为已支持的化学物种。

## 输入与错误行为

逐帧解析，保留 timestep、输出截断、原始原子字段和起止行号。
拒绝重复原子、负邻接数、自连接、重复邻居、未知邻居、非有限数值和原子行缺失。
连接入口要求明确的原子总数与元素映射；后续帧的 ID 到类型映射必须与第一帧相同。

默认遇到损坏帧、重复或倒退时间、原子增删/类型改变、间隔变化时立即报错。
本阶段没有跨错误帧生成事件，也没有自动分段恢复。建议显式提供 expected_interval；
未提供时以首两个有效帧的差值作为后续检查基准，因此不能识别一开始就缺帧，
也不能发现整个文件始终按同一倍数抽帧的情况。

缺少 timestep_scale 时 time_ps 为 JSON null，只按 timestep 工作。
本阶段暂不支持时间重置、变步长或坐标/PBC 对齐。

## BO 归并和连接状态

用排序后的两个原子 ID 标识原子对。两个方向的 BO 都保存。
差值不大于显式 tolerance 时取均值并保留说明；超过 tolerance 或仅有单向记录时停止。
默认 tolerance 为 0，即只接受一致值。这是保守的首版实现策略，尚未经全轨迹异常分布校准。

连接仅由当前帧的邻接记录决定：列出的边全部保留，不再比较 BO 数值，
没有 formation、breaking、元素对阈值，也没有历史滞后状态。
上一帧存在而当前帧不再列出的边，记录 merged_order=null、空双向观测和
not_listed_in_frame 原因；这不等于物理 BO=0。列出的边记录 listed_in_frame。

output_cutoff 仅记录写文件时的截断元数据，不参与过滤。头部未提供时可为 null；
显式配置与头部冲突仍报错。打印 BO 等于输出 cutoff 时也保留该邻接记录，
避免再次比较打印数值导致舍入附近的连接被删除。

这里不以价态删掉已接受的连接，也不将 BO 四舍五入成化学键类型。
例如接受五条 C-H 边会产生组成 CH5 的连接实例；其化学解释留给结构推断层。

## 分子实例与可追溯记录

对所有原子和接受的边求连通分量，保留孤立原子。实例 ID 形如
F00000000:M00000000，属于特定保存帧，不是跨帧持续分子身份。
每个实例包含原子集合、分子式和连接集合。跨帧对应依靠原子 ID。

- frames.jsonl：每帧一行，含原始原子字段、BO 观测、连接决策、源文件行号，
  以及对应实例在实例文件中的首行和条数。
- molecule_instances.jsonl：每个逐帧实例一行；通过 frame_id 回查完整证据。
- run_manifest.json：运行 ID、规则、显式元素映射、时间、输入定位及大小/mtime、
  Python/平台版本、运行状态、耗时和数据输出体积。

输出目录必须不存在，避免覆盖前一次运行。完整到 EOF 为 complete，主动抽取有限帧
为 partial，异常或中断为 incomplete。设置 max_frames 时即使恰好达到 EOF 也保守标记 partial。
运行中的 manifest 首先标记 incomplete；正常结束或捕获异常时更新统计。
输入身份目前使用路径、大小及 mtime，不声称已做内容哈希校验。

解析、连接及输出仅保留当前帧和上一帧的原子类型及连接状态，无整段轨迹缓存；
单帧记录会在序列化时有临时副本。连接分量遍历是 O(N+E)，确定性排序带来额外排序开销。
JSONL 是首版可追溯交换记录，目前没有 SQLite 等随机访问索引，不承诺固定内存上限。

## 运行示例

```powershell
cd D:\smiletosvg
python -m reaxnet.connectivity `
  --bonds examples\ch4_dissociation.bonds `
  --elements '1:C,2:H' `
  --output-dir outputs\connectivity_demo
```

真实样本小段：

```powershell
python -m reaxnet.connectivity `
  --bonds D:\cellulose\0.9d3000\bonds.reax `
  --elements '1:H,2:C,3:O' `
  --output-dir outputs\connectivity_real_sample `
  --expected-interval 1000 --timestep-scale 0.0002 `
  --max-frames 20
```

已删除的阈值参数会被命令行拒绝。当前 schema 为 connectivity-2，算法为 file-topology-2。
新策略真实样本与测试见 [2026-09-22 更新](file_topology_and_target_matching_2026-09-22.md)。

## 2026-09-17 历史验证（旧阈值版本）

`python -m unittest discover -s tests -v`：27 项通过，其中 14 项新增。
覆盖阈值等号、双阈值初始化/滞后、输出截断、缺失与零值、双向冲突、元素对覆盖、
孤立原子、帧内 ID、原子增删及类型变化、时间缺口、原始字段保留、异常状态与防覆盖。
现有事件、标识符和 Web 测试仍通过，这不代表旧分析器已满足新 PRD。

真实数据输入大小经本次核查为 4,424,628,012 字节。
配套 in.celluloseTD 确认元素顺序为 H C O，units real，积分步长 0.2 fs。
此次处理前 20 帧，timestep 0–19000，即 0–3.8 ps，每帧 10,330 个原子。
在此次显式 0.50/0.35 规则下，共保存 40,768 条逐帧实例，首帧 2015 个、末帧 2078 个，
没有发现双向 BO 差异。输出约 60.6 MB。
这些统计没有证明连接判据科学正确，逐帧实例条数也不等于独立分子总数。

本次没有重新扫描全部 5001 帧，没有测定全文件峰值内存，未完成 A01 全量检查和 A12
全文件运行验收，也未完成阈值敏感性验收。不要将这次小段结果外推为后期高反应区的验证。

## 后续衔接

下一阶段使用实例连接图和实际元素组成，显式声明电子状态假设，搜索化学键类型、
电荷及自由基候选，再核对实际氢数，生成 supported/ambiguous/unresolved 身份。
CHO 小案例优先覆盖 CH3、CH2O、乙醇/二甲醚；候选身份不能通过补氢或强制单键回退伪装成确定结果。
随后再接事件区间、往返标记、关键产物与观察路径，并评估实例索引及压缩存储。
