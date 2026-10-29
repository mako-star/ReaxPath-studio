"""Read count time series and backfill old jobs without rerunning reactions."""
from __future__ import annotations

import asyncio
import csv
import json
import math
from pathlib import Path
import time

from reaxnet.parser import parse_element_mapping

from . import job_manager as state_module
from .job_manager import JobStatus


COUNT_FILES = ("bond_counts", "bond_changes", "functional_groups",
               "functional_group_molecules", "count_metadata")
GROUP_LABELS = {
    "hydroxyl": "羟基候选", "ether": "醚候选", "aldehyde": "醛基候选",
    "ketone": "酮型羰基候选", "carboxyl": "羧基候选", "ester": "酯基候选",
}


class CountOutputError(ValueError):
    """A persisted count output is incomplete, corrupt or incompatible."""


def _metadata(snapshot):
    filename = snapshot.get("result_files", {}).get("count_metadata")
    if not filename:
        return None
    try:
        result = json.loads(Path(filename).read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise CountOutputError("计数元数据无法读取，请检查输出文件。") from exc
    if not isinstance(result, dict) or result.get("status") != "completed":
        raise CountOutputError("计数输出尚未完整提交，不能作为完整曲线显示。")
    count = result.get("frame_count")
    if type(count) is not int or count < 1:
        raise CountOutputError("计数元数据缺少有效帧数。")
    try:
        scale = float(result["timestep_scale_ps"])
        start_time = float(result["time_start_ps"])
        end_time = float(result["time_end_ps"])
    except (KeyError, TypeError, ValueError) as exc:
        raise CountOutputError("计数元数据缺少有效的时间换算或时间范围。") from exc
    if not all(math.isfinite(value) for value in (scale, start_time, end_time)) or scale <= 0 or start_time > end_time:
        raise CountOutputError("计数元数据的时间换算或时间范围无效。")
    configured_scale = snapshot.get("config", {}).get("timestep_scale")
    if configured_scale is not None:
        try:
            matches_scale = math.isclose(scale, float(configured_scale), rel_tol=1e-12, abs_tol=1e-15)
        except (TypeError, ValueError) as exc:
            raise CountOutputError("原任务的时间换算参数无效。") from exc
        if not matches_scale:
            raise CountOutputError("计数时间换算与原任务参数不一致。")
    expected = snapshot.get("summary", {}).get("frames")
    if type(expected) is int and expected != count:
        raise CountOutputError("计数帧数与原分析范围不一致。")
    return result


def _spec(kind, metric):
    choices = {
        ("bonds", "count"): ("bond_counts", ""),
        ("bonds", "formed"): ("bond_changes", "_formed"),
        ("bonds", "broken"): ("bond_changes", "_broken"),
        ("functional", "count"): ("functional_groups", ""),
        ("functional", "molecules"): ("functional_group_molecules", ""),
    }
    try:
        return choices[(kind, metric)]
    except KeyError as exc:
        raise ValueError("该分析类型不支持所选计数指标。") from exc


def _sample_csv(path, metadata, categories, suffix, max_points):
    """Keep the full range and each bucket's extrema, using bounded storage."""
    expected = metadata["frame_count"]
    keys = [item["key"] for item in categories]
    if max_points < 2 * len(keys) + 2:
        raise ValueError("max_points 太小，无法保留各条曲线的极值。")
    bucket_count = max(1, (max_points - 2) // (2 * len(keys)))
    buckets = [{} for _ in range(bucket_count)] if expected > max_points else None
    selected = {}
    previous_step = previous_time = None
    first_step = last_step = first_time = last_time = None
    count = 0
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            required = {"timestep", "time_ps", *(key + suffix for key in keys)}
            if not required.issubset(reader.fieldnames or ()):
                raise CountOutputError("计数 CSV 缺少必要列。")
            for index, original in enumerate(reader):
                if index >= expected:
                    raise CountOutputError("计数 CSV 帧数超过已记录的分析范围。")
                step = int(original["timestep"])
                physical_time = float(original["time_ps"])
                if not math.isfinite(physical_time):
                    raise CountOutputError("计数 CSV 包含非有限时间。")
                if not math.isclose(physical_time, step * float(metadata["timestep_scale_ps"]),
                                    rel_tol=1e-12, abs_tol=1e-9):
                    raise CountOutputError("计数 CSV 的时间与 step 换算不一致。")
                if previous_step is not None and (step <= previous_step or physical_time <= previous_time):
                    raise CountOutputError("计数 CSV 的时间或 step 顺序无效。")
                row = {"timestep": step, "time_ps": physical_time}
                for key in keys:
                    raw = original[key + suffix]
                    value = None if raw == "" else int(raw)
                    if value is None and (not suffix or index != 0):
                        raise CountOutputError("计数 CSV 存在缺失数量。")
                    if value is not None and value < 0:
                        raise CountOutputError("计数 CSV 包含负数量。")
                    row[key] = value
                first_step = step if first_step is None else first_step
                last_step = step
                first_time = physical_time if first_time is None else first_time
                last_time = physical_time
                previous_step, previous_time = step, physical_time
                count += 1
                if buckets is None or index in (0, expected - 1):
                    selected[index] = row
                if buckets is not None:
                    bucket = buckets[min(bucket_count - 1, index * bucket_count // expected)]
                    for key in keys:
                        value = row[key]
                        if value is None:
                            continue
                        for mode, improves in (("min", lambda x, y: x < y), ("max", lambda x, y: x > y)):
                            tag = (key, mode)
                            old = bucket.get(tag)
                            if old is None or improves(value, old[1][key]):
                                bucket[tag] = (index, row)
    except CountOutputError:
        raise
    except (OSError, TypeError, ValueError, KeyError, csv.Error) as exc:
        raise CountOutputError("计数 CSV 无法读取或数量格式无效。") from exc
    if count != expected or first_step != metadata.get("first_timestep") or last_step != metadata.get("last_timestep"):
        raise CountOutputError("计数 CSV 与元数据的帧数或 step 范围不一致。")
    if not math.isclose(first_time, float(metadata["time_start_ps"]), rel_tol=1e-12, abs_tol=1e-9) or not math.isclose(
            last_time, float(metadata["time_end_ps"]), rel_tol=1e-12, abs_tol=1e-9):
        raise CountOutputError("计数 CSV 与元数据的时间范围不一致。")
    if buckets is not None:
        for bucket in buckets:
            for index, row in bucket.values():
                selected[index] = row
    return [selected[index] for index in sorted(selected)]


def read_series(snapshot, *, kind="bonds", metric="count", max_points=5000):
    file_key, suffix = _spec(kind, metric)
    frame_count = snapshot.get("summary", {}).get("frames", 0)
    base = {"kind": kind, "metric": metric, "frame_count": frame_count,
            "returned_points": 0, "downsampled": False, "categories": [],
            "rows": [], "metadata": None, "csv_url": None}
    auxiliary = snapshot.get("count_analysis", {})
    if auxiliary.get("status") == "running":
        return {**base, **auxiliary, "status": "running"}
    metadata = _metadata(snapshot)
    if metadata is None:
        if auxiliary.get("status") == "failed":
            return {**base, **auxiliary, "status": "failed"}
        return {**base, "status": "missing", "reason": "此任务尚未生成键和官能团计数，可单独补算。"}
    category_key = "bond_categories" if kind == "bonds" else "functional_group_categories"
    raw_categories = metadata.get(category_key)
    if not isinstance(raw_categories, list) or not raw_categories:
        raise CountOutputError("计数元数据缺少曲线分类。")
    categories = []
    seen = set()
    for item in raw_categories:
        if not isinstance(item, dict) or not isinstance(item.get("key"), str) or not item["key"] or item["key"] in seen:
            raise CountOutputError("计数元数据的曲线分类无效。")
        seen.add(item["key"])
        categories.append({"key": item["key"], "label": str(item.get("label", item["key"]))})
    filename = snapshot.get("result_files", {}).get(file_key)
    if not filename:
        raise CountOutputError("此任务的计数 CSV 未登记。")
    rows = _sample_csv(Path(filename), metadata, categories, suffix, max_points)
    return {**base, "status": "completed", "frame_count": metadata["frame_count"],
            "returned_points": len(rows), "downsampled": len(rows) < metadata["frame_count"],
            "sampling": "bucket_extrema" if len(rows) < metadata["frame_count"] else "all_frames",
            "metadata": metadata, "categories": categories, "rows": rows,
            "csv_url": f"/api/jobs/{snapshot['job_id']}/outputs/{file_key}"}


def _backfill_parameters(snapshot):
    config = snapshot.get("config", {})
    filename = config.get("source_path")
    if not isinstance(filename, str) or not filename:
        raise ValueError("旧任务没有保存原始文件路径，无法补算；请重新提交原始文件。")
    source = Path(filename)
    if not source.is_absolute() or not source.is_file():
        raise ValueError("任务记录中的原始 bonds 文件不存在。")
    stat = source.stat()
    expected_size = config.get("source_size_bytes", snapshot.get("summary", {}).get("source_size_bytes"))
    expected_mtime = config.get("source_mtime_ns")
    if expected_size is not None and stat.st_size != expected_size:
        raise ValueError("原始文件大小与原任务记录不一致，停止补算。")
    if expected_mtime is not None and stat.st_mtime_ns != expected_mtime:
        raise ValueError("原始文件修改时间与原任务记录不一致，停止补算。")
    frame_limit = snapshot.get("summary", {}).get("frames")
    if type(frame_limit) is not int or frame_limit < 1:
        raise ValueError("旧任务缺少明确的已分析帧数，无法安全确定补算范围。")
    raw_mapping = snapshot.get("element_mapping", {})
    if not isinstance(raw_mapping, dict):
        raise ValueError("任务的元素映射无效。")
    mapping_text = raw_mapping.get("_raw")
    if mapping_text is None:
        mapping_text = ",".join(f"{key}:{value}" for key, value in raw_mapping.items())
    mapping = parse_element_mapping(mapping_text)
    scale = float(config.get("timestep_scale", 1.0))
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError("任务的时间换算无效。")
    interval = config.get("expected_interval")
    return source.resolve(), mapping, scale, interval, frame_limit, stat


async def prepare_backfill(manager, job_id):
    """Claim a job before scheduling a worker, including concurrent requests."""
    async with manager._lock:
        job = manager._jobs.get(job_id)
        if job is None:
            raise KeyError(job_id)
        if job.status != JobStatus.COMPLETED:
            raise ValueError("原反应分析尚未完成，暂不能补算。")
        snapshot = job.to_dict()
        if job.count_analysis.get("status") == "running":
            return dict(job.count_analysis), None
        if all(name in job.result_files and Path(job.result_files[name]).is_file() for name in COUNT_FILES):
            metadata = _metadata(snapshot)
            if metadata is not None:
                return {"status": "completed", "frame_count": metadata["frame_count"], "progress": 100.0}, None
        source, mapping, scale, interval, frame_limit, stat = _backfill_parameters(snapshot)
        output = state_module.JOB_ROOT / job_id / "count_analysis" / f"run-{time.time_ns()}"
        job.count_analysis = {"status": "running", "progress": 0.0, "frames_processed": 0,
                              "frame_limit": frame_limit, "started_at": time.time(),
                              "message": "正在逐帧补算键和官能团数量"}
        manager._persist(job)
        return dict(job.count_analysis), (source, mapping, scale, interval, frame_limit, output,
                                         (stat.st_size, stat.st_mtime_ns))


async def run_backfill(manager, job_id, parameters):
    source, mapping, scale, interval, frame_limit, output, input_stat = parameters
    loop = asyncio.get_running_loop()
    pending = None
    last_report = 0.0

    def progress(payload):
        nonlocal pending, last_report
        now = time.monotonic()
        frames = int(payload.get("frames_processed", 0))
        if now - last_report < 0.5 and frames < frame_limit:
            return
        last_report = now
        pending = asyncio.run_coroutine_threadsafe(manager.update(job_id, count_analysis={
            "status": "running", "progress": min(99.0, 100 * frames / frame_limit),
            "frames_processed": frames, "frame_limit": frame_limit,
            "message": f"已处理 {frames:,}/{frame_limit:,} 帧"}), loop)

    try:
        from reaxnet.counts import analyze_counts_file

        current_stat = source.stat()
        if (current_stat.st_size, current_stat.st_mtime_ns) != input_stat:
            raise ValueError("原始文件在补算启动前发生变化。")
        result = await asyncio.to_thread(analyze_counts_file, source, mapping, output,
                                        timestep_scale=scale, expected_interval=interval,
                                        max_frames=frame_limit, progress_callback=progress)
        if pending is not None:
            await asyncio.wrap_future(pending)
        if result.summary.get("frame_count") != frame_limit:
            raise ValueError("补算帧数与原任务范围不一致。")
        if result.metadata.get("status") != "completed":
            raise ValueError("补算输出未完整提交。")
        if not all(name in result.paths and Path(result.paths[name]).is_file() for name in COUNT_FILES):
            raise ValueError("补算结果缺少已完成的计数文件。")
        snapshot = manager.snapshot(job_id)
        if snapshot is None:
            return
        files = {**snapshot.get("result_files", {}),
                 **{key: str(Path(result.paths[key]).resolve()) for key in COUNT_FILES}}
        await manager.update(job_id, result_files=files, count_analysis={
            "status": "completed", "progress": 100.0, "frame_count": frame_limit,
            "frames_processed": frame_limit, "finished_at": time.time(),
            "message": "键和官能团计数完成"})
    except Exception as exc:
        if pending is not None:
            try:
                await asyncio.wrap_future(pending)
            except Exception:
                pass
        await manager.update(job_id, count_analysis={"status": "failed", "error": f"{type(exc).__name__}: {exc}",
                                                    "finished_at": time.time(), "message": "计数补算失败，原反应结果保留"})
