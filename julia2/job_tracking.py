"""
Track submitted jobs and report their status.
"""

import csv
import logging
import os
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
