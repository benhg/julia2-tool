"""
Track submitted jobs and report their status.
"""

import csv
import logging
import os
import shutil
import subprocess
from collections import Counter
from datetime import datetime, timezone
from statistics import median

logger = logging.getLogger("julia2.job_tracking")

JOB_HEADERS = [
    "job_id",
    "job_name",
    "job_type",
    "index_id",
    "reads_sample_id",
    "sbatch_script",
    "stdout_path",
    "submitted_at",
    "completed_at",
    "state",
    "elapsed",
    "exit_code",
    "reason",
]

ACTIVE_STATES = {"PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "SUSPENDED"}
TERMINAL_STATES = {
    "COMPLETED",
    "FAILED",
    "CANCELLED",
    "TIMEOUT",
    "OUT_OF_MEMORY",
    "NODE_FAIL",
    "PREEMPTED",
}
FAILED_STATES = TERMINAL_STATES - {"COMPLETED"}


def job_status_path(project_config):
    return f"{project_config.project_dir}/output/job_status.csv"


def ensure_job_status_file(project_config):
    path = job_status_path(project_config)
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return path
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=JOB_HEADERS)
        writer.writeheader()
    return path


def append_job_record(project_config, row):
    ensure_job_status_file(project_config)
    with open(job_status_path(project_config), "a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=JOB_HEADERS)
        writer.writerow(row)


def load_job_records(project_config):
    path = ensure_job_status_file(project_config)
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for row in rows:
        for header in JOB_HEADERS:
            row.setdefault(header, "")
    return rows


def write_job_records(project_config, rows):
    path = ensure_job_status_file(project_config)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=JOB_HEADERS)
        writer.writeheader()
        writer.writerows(rows)


def remove_job_records(project_config, predicate):
    rows = load_job_records(project_config)
    kept_rows = [row for row in rows if not predicate(row)]
    write_job_records(project_config, kept_rows)
    return len(rows) - len(kept_rows)


def make_job_record(job_id, sbatch_name, job_type, sbatch_script, stdout_path,
                    index_id="", reads_sample_id="", state="SUBMITTED",
                    elapsed="", exit_code="", reason=""):
    return {
        "job_id": job_id,
        "job_name": sbatch_name,
        "job_type": job_type,
        "index_id": index_id,
        "reads_sample_id": reads_sample_id,
        "sbatch_script": sbatch_script,
        "stdout_path": stdout_path,
        "submitted_at": datetime.now(timezone.utc).isoformat(),
        "completed_at": "",
        "state": state,
        "elapsed": elapsed,
        "exit_code": exit_code,
        "reason": reason,
    }


def parse_sbatch_job_id(output_text):
    for token in output_text.strip().split():
        if token.isdigit():
            return token
    raise ValueError(f"Could not parse job id from sbatch output: {output_text!r}")


def _run_status_cmd(cmd):
    try:
        return subprocess.check_output(cmd, shell=True, text=True).splitlines()
    except subprocess.CalledProcessError as exc:
        logger.warning("Status command failed: %s", exc)
        return []


def _query_squeue(job_ids):
    if not job_ids:
        return {}
    ids = ",".join(job_ids)
    lines = _run_status_cmd(
        f'squeue -h -j {ids} -o "%i|%T|%M|%R"'
    )
    results = {}
    for line in lines:
        parts = line.split("|", 3)
        if len(parts) != 4:
            continue
        job_id, state, elapsed, reason = parts
        results[job_id] = {
            "state": state,
            "elapsed": elapsed,
            "reason": reason,
            "exit_code": "",
        }
    return results


def _query_sacct(job_ids):
    if not job_ids:
        return {}
    ids = ",".join(job_ids)
    lines = _run_status_cmd(
        f'sacct -n -X -P -j {ids} --format=JobIDRaw,State,Elapsed,ExitCode,Reason'
    )
    results = {}
    for line in lines:
        parts = line.split("|", 4)
        if len(parts) != 5:
            continue
        job_id, state, elapsed, exit_code, reason = parts
        if not job_id.isdigit() or "." in job_id:
            continue
        results[job_id] = {
            "state": state.split()[0],
            "elapsed": elapsed,
            "exit_code": exit_code,
            "reason": reason,
        }
    return results


def _query_sacct_node_data(job_ids):
    if not job_ids:
        return {}
    ids = ",".join(job_ids)
    lines = _run_status_cmd(
        f"sacct -n -X -P -j {ids} --format=JobIDRaw,NodeList,AllocNodes,AllocCPUS,ElapsedRaw"
    )
    results = {}
    for line in lines:
        parts = line.split("|", 4)
        if len(parts) != 5:
            continue
        job_id, node_list, alloc_nodes, alloc_cpus, elapsed_raw = parts
        if not job_id.isdigit() or "." in job_id:
            continue
        try:
            alloc_nodes_value = int(alloc_nodes) if alloc_nodes else 0
        except ValueError:
            alloc_nodes_value = 0
        try:
            alloc_cpus_value = int(alloc_cpus) if alloc_cpus else 0
        except ValueError:
            alloc_cpus_value = 0
        try:
            elapsed_raw_value = int(elapsed_raw) if elapsed_raw else None
        except ValueError:
            elapsed_raw_value = None
        results[job_id] = {
            "node_list": node_list.strip(),
            "alloc_nodes": alloc_nodes_value,
            "alloc_cpus": alloc_cpus_value,
            "elapsed_raw": elapsed_raw_value,
        }
    return results


def _expand_slurm_nodelist(node_list_text):
    if not node_list_text:
        return []
    node_list_text = node_list_text.strip()
    if not node_list_text or node_list_text in {"Unknown", "None", "N/A"}:
        return []

    scontrol = shutil.which("scontrol")
    if scontrol:
        try:
            expanded = subprocess.check_output(
                [scontrol, "show", "hostnames", node_list_text],
                text=True
            ).splitlines()
            return [node.strip() for node in expanded if node.strip()]
        except (subprocess.CalledProcessError, OSError):
            logger.warning("Failed to expand SLURM nodelist: %s", node_list_text)

    if "," in node_list_text and "[" not in node_list_text and "]" not in node_list_text:
        return [part.strip() for part in node_list_text.split(",") if part.strip()]
    return [node_list_text]


def _format_cpu_hours(cpu_hours):
    if cpu_hours is None:
        return "unknown"
    return f"{cpu_hours:.2f} CPU-hours"


def refresh_job_statuses(project_config, use_slurm):
    rows = load_job_records(project_config)
    slurm_ids = [row["job_id"] for row in rows if row["job_id"].isdigit()]
    now = datetime.now(timezone.utc).isoformat()

    if use_slurm and slurm_ids:
        live_status = _query_squeue(slurm_ids)
        acct_status = _query_sacct(slurm_ids)
    else:
        live_status = {}
        acct_status = {}

    updated_rows = []
    for row in rows:
        job_id = row["job_id"]
        previous_state = row["state"]
        status = live_status.get(job_id) or acct_status.get(job_id)
        if status:
            row["state"] = status.get("state", row["state"]) or row["state"]
            row["elapsed"] = status.get("elapsed", row["elapsed"]) or row["elapsed"]
            row["exit_code"] = status.get("exit_code", row["exit_code"]) or row["exit_code"]
            row["reason"] = status.get("reason", row["reason"]) or row["reason"]
        elif not use_slurm and row["state"] == "SUBMITTED":
            row["state"] = "RUNNING"
        if row["state"] in TERMINAL_STATES and not row.get("completed_at"):
            row["completed_at"] = _infer_completed_at(row, now)
        elif previous_state in TERMINAL_STATES and row["state"] not in TERMINAL_STATES:
            row["completed_at"] = ""
        updated_rows.append(row)

    write_job_records(project_config, updated_rows)
    return updated_rows


def summarize_jobs(rows):
    counts = Counter(row["state"] or "UNKNOWN" for row in rows)
    return counts


def _parse_elapsed_to_seconds(elapsed_text):
    if not elapsed_text:
        return None
    elapsed_text = elapsed_text.strip()
    if not elapsed_text or elapsed_text in {"Unknown", "INVALID"}:
        return None

    day_count = 0
    time_text = elapsed_text
    if "-" in elapsed_text:
        day_text, time_text = elapsed_text.split("-", 1)
        if not day_text.isdigit():
            return None
        day_count = int(day_text)

    parts = time_text.split(":")
    if len(parts) == 3:
        hours, minutes, seconds = parts
    elif len(parts) == 2:
        hours = "0"
        minutes, seconds = parts
    else:
        return None

    try:
        total_seconds = (day_count * 24 * 3600 + int(hours) * 3600 +
                         int(minutes) * 60 + int(seconds))
    except ValueError:
        return None
    return total_seconds


def _format_seconds(total_seconds):
    if total_seconds is None:
        return "unknown"
    total_seconds = int(round(total_seconds))
    days, remainder = divmod(total_seconds, 24 * 3600)
    hours, remainder = divmod(remainder, 3600)
    minutes, seconds = divmod(remainder, 60)
    if days:
        return f"{days}d {hours}h {minutes}m"
    if hours:
        return f"{hours}h {minutes}m"
    if minutes:
        return f"{minutes}m {seconds}s"
    return f"{seconds}s"


def _parse_timestamp(timestamp_text):
    if not timestamp_text:
        return None
    timestamp_text = timestamp_text.strip()
    if not timestamp_text:
        return None
    try:
        return datetime.fromisoformat(timestamp_text)
    except ValueError:
        return None


def _infer_completed_at(row, default_time):
    submitted_at = _parse_timestamp(row.get("submitted_at"))
    elapsed_seconds = _parse_elapsed_to_seconds(row.get("elapsed"))
    if submitted_at is not None and elapsed_seconds is not None:
        return datetime.fromtimestamp(
            submitted_at.timestamp() + elapsed_seconds, timezone.utc
        ).isoformat()
    return default_time


def _effective_completed_at(row):
    recorded_completed_at = _parse_timestamp(row.get("completed_at"))
    submitted_at = _parse_timestamp(row.get("submitted_at"))
    elapsed_seconds = _parse_elapsed_to_seconds(row.get("elapsed"))
    inferred_completed_at = None
    if submitted_at is not None and elapsed_seconds is not None:
        inferred_completed_at = datetime.fromtimestamp(
            submitted_at.timestamp() + elapsed_seconds, timezone.utc
        )

    if recorded_completed_at is None:
        return inferred_completed_at
    if inferred_completed_at is None:
        return recorded_completed_at

    # Older rows may have been backfilled with "now" when completed_at was introduced.
    # If the recorded completion is much later than the historical lower bound, prefer
    # the inferred timestamp for ETA calculations to avoid a fake burst of completions.
    if (recorded_completed_at - inferred_completed_at).total_seconds() > 3600:
        return inferred_completed_at
    return recorded_completed_at


def _estimate_remaining_time_from_completion_rate(rows):
    active_rows = [
        row for row in rows
        if row["state"] in ACTIVE_STATES or row["state"] == "SUBMITTED"
    ]
    if not active_rows:
        return None

    completed_rows = [
        row for row in rows
        if row["state"] == "COMPLETED" and _effective_completed_at(row)
    ]
    completed_rows.sort(key=_effective_completed_at)
    if len(completed_rows) < 3:
        return None

    earliest_completion = _effective_completed_at(completed_rows[0])
    if earliest_completion is None:
        return None

    now = datetime.now(timezone.utc)
    observation_seconds = (now - earliest_completion).total_seconds()
    if observation_seconds <= 0:
        return None

    completion_rate = len(completed_rows) / observation_seconds
    if completion_rate <= 0:
        return None

    return {
        "model": "completion_rate",
        "completed_samples": len(completed_rows),
        "window_seconds": observation_seconds,
        "completion_rate_per_hour": completion_rate * 3600,
        "active_jobs": len(active_rows),
        "estimated_wall_seconds": len(active_rows) / completion_rate,
    }


def _estimate_remaining_time_from_runtime_heuristic(rows):
    completed_seconds = [
        _parse_elapsed_to_seconds(row["elapsed"])
        for row in rows
        if row["state"] == "COMPLETED"
    ]
    completed_seconds = [value for value in completed_seconds if value is not None]

    active_rows = [
        row for row in rows
        if row["state"] in ACTIVE_STATES or row["state"] == "SUBMITTED"
    ]
    running_rows = [row for row in active_rows if row["state"] == "RUNNING"]
    pending_rows = [row for row in active_rows if row["state"] != "RUNNING"]

    if not active_rows or not completed_seconds:
        return None

    median_runtime = median(completed_seconds)
    running_elapsed = [
        _parse_elapsed_to_seconds(row["elapsed"])
        for row in running_rows
    ]
    running_remaining = [
        max(median_runtime - elapsed, 0)
        for elapsed in running_elapsed
        if elapsed is not None
    ]

    if running_rows:
        current_wave_remaining = max(running_remaining) if running_remaining else median_runtime
        current_wave_remaining = max(current_wave_remaining, 0)
    else:
        current_wave_remaining = 0

    concurrency = max(len(running_rows), 1)
    queued_jobs = len(pending_rows)
    queued_waves = (queued_jobs + concurrency - 1) // concurrency
    queued_wall_seconds = queued_waves * median_runtime
    estimated_wall_seconds = current_wave_remaining + queued_wall_seconds

    return {
        "model": "median_runtime",
        "completed_samples": len(completed_seconds),
        "median_runtime_seconds": median_runtime,
        "active_jobs": len(active_rows),
        "running_jobs": len(running_rows),
        "pending_jobs": len(pending_rows),
        "estimated_wall_seconds": estimated_wall_seconds,
        "estimated_serial_seconds": len(active_rows) * median_runtime,
    }


def _estimate_remaining_time(rows):
    return (
        _estimate_remaining_time_from_completion_rate(rows)
        or _estimate_remaining_time_from_runtime_heuristic(rows)
    )


def format_job_status_report(rows):
    counts = summarize_jobs(rows)
    ordered_states = [
        "SUBMITTED",
        "PENDING",
        "RUNNING",
        "COMPLETED",
        "FAILED",
        "CANCELLED",
        "TIMEOUT",
        "OUT_OF_MEMORY",
        "NODE_FAIL",
        "PREEMPTED",
    ]
    lines = [f"Total jobs: {len(rows)}"]
    for state in ordered_states:
        if counts.get(state):
            lines.append(f"{state}: {counts[state]}")
    other_states = sorted(set(counts) - set(ordered_states))
    for state in other_states:
        lines.append(f"{state}: {counts[state]}")

    eta = _estimate_remaining_time(rows)
    if eta:
        lines.append("")
        lines.append("Estimated remaining time:")
        if eta["model"] == "completion_rate":
            lines.append(
                f"Observed completion rate: {eta['completion_rate_per_hour']:.2f} jobs/hour "
                f"from {eta['completed_samples']} completed jobs "
                f"over {_format_seconds(eta['window_seconds'])}"
            )
            lines.append(
                f"Estimated wall time to drain active queue at the recent completion rate: "
                f"{_format_seconds(eta['estimated_wall_seconds'])}"
            )
        else:
            lines.append(
                f"Median completed runtime: {_format_seconds(eta['median_runtime_seconds'])} "
                f"from {eta['completed_samples']} completed jobs"
            )
            lines.append(
                f"Estimated wall time to drain active queue at current concurrency: "
                f"{_format_seconds(eta['estimated_wall_seconds'])}"
            )
            lines.append(
                f"Estimated aggregate compute time remaining: "
                f"{_format_seconds(eta['estimated_serial_seconds'])}"
            )
    elif any(row["state"] in ACTIVE_STATES or row["state"] == "SUBMITTED" for row in rows):
        lines.append("")
        lines.append("Estimated remaining time:")
        lines.append("Not enough completed jobs with timing history to estimate yet")

    failed_rows = [
        row for row in rows
        if row["state"] in TERMINAL_STATES - {"COMPLETED"}
    ]
    if failed_rows:
        lines.append("")
        lines.append("Failed jobs:")
        for row in failed_rows[:20]:
            lines.append(
                f'{row["job_id"]} {row["job_type"]} {row["index_id"]} {row["reads_sample_id"]} '
                f'{row["state"]} {row["reason"]} {row["stdout_path"]}'
            )
        if len(failed_rows) > 20:
            lines.append(f"... {len(failed_rows) - 20} more")

    active_rows = [
        row for row in rows
        if row["state"] in ACTIVE_STATES or row["state"] == "SUBMITTED"
    ]
    if active_rows:
        lines.append("")
        lines.append("Active jobs:")
        for row in active_rows[:20]:
            lines.append(
                f'{row["job_id"]} {row["job_type"]} {row["index_id"]} {row["reads_sample_id"]} '
                f'{row["state"]} {row["elapsed"]} {row["reason"]}'
            )
        if len(active_rows) > 20:
            lines.append(f"... {len(active_rows) - 20} more")

    return "\n".join(lines)


def format_job_time_by_node_report(rows, use_slurm):
    if not rows:
        return "No tracked jobs found"
    if not use_slurm:
        return "Node time breakdown is only available when SLURM mode is enabled"

    slurm_ids = [row["job_id"] for row in rows if row["job_id"].isdigit()]
    if not slurm_ids:
        return "No tracked SLURM jobs found"

    node_data = _query_sacct_node_data(slurm_ids)
    totals_by_node = Counter()
    cpu_hours_by_node = Counter()
    jobs_by_node = Counter()
    jobs_without_node = 0
    jobs_without_elapsed = 0
    jobs_without_cpu_count = 0

    for row in rows:
        job_id = row["job_id"]
        if job_id not in node_data:
            jobs_without_node += 1
            continue

        job_node_data = node_data[job_id]
        node_names = _expand_slurm_nodelist(job_node_data["node_list"])
        if not node_names:
            jobs_without_node += 1
            continue

        elapsed_seconds = job_node_data["elapsed_raw"]
        if elapsed_seconds is None:
            elapsed_seconds = _parse_elapsed_to_seconds(row.get("elapsed"))
        if elapsed_seconds is None:
            jobs_without_elapsed += 1
            continue

        seconds_per_node = elapsed_seconds / max(len(node_names), 1)
        alloc_cpus = job_node_data.get("alloc_cpus", 0)
        cpu_hours_per_node = ((elapsed_seconds * alloc_cpus) / 3600.0) / max(len(node_names), 1)
        if alloc_cpus <= 0:
            jobs_without_cpu_count += 1
        for node_name in node_names:
            totals_by_node[node_name] += seconds_per_node
            cpu_hours_by_node[node_name] += cpu_hours_per_node
            jobs_by_node[node_name] += 1

    if not totals_by_node:
        return "No node timing data available yet from sacct"

    lines = ["Job time by node:"]
    for node_name, total_seconds in sorted(
        totals_by_node.items(),
        key=lambda item: (-item[1], item[0])
    ):
        average_seconds = total_seconds / jobs_by_node[node_name]
        lines.append(
            f"{node_name}: {_format_seconds(total_seconds)} total, "
            f"{_format_seconds(average_seconds)} average/job across {jobs_by_node[node_name]} jobs, "
            f"{_format_cpu_hours(cpu_hours_by_node[node_name])} total"
        )

    lines.append("")
    lines.append(f"Nodes with timing data: {len(totals_by_node)}")
    lines.append(f"Jobs with node data: {sum(jobs_by_node.values())}")
    lines.append(f"Total CPU-hours: {_format_cpu_hours(sum(cpu_hours_by_node.values()))}")
    lines.append(
        f"Average time per node-attributed job: "
        f"{_format_seconds(sum(totals_by_node.values()) / sum(jobs_by_node.values()))}"
    )
    if jobs_without_node:
        lines.append(f"Jobs missing node data: {jobs_without_node}")
    if jobs_without_elapsed:
        lines.append(f"Jobs missing elapsed time: {jobs_without_elapsed}")
    if jobs_without_cpu_count:
        lines.append(f"Jobs missing CPU allocation data: {jobs_without_cpu_count}")

    return "\n".join(lines)
