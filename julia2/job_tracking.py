"""
Track submitted jobs and report their status.
"""

import csv
import logging
import os
import subprocess
from collections import Counter
from datetime import datetime, timezone

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
        return list(csv.DictReader(fh))


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

    if use_slurm and slurm_ids:
        live_status = _query_squeue(slurm_ids)
        acct_status = _query_sacct(slurm_ids)
    else:
        live_status = {}
        acct_status = {}

    updated_rows = []
    for row in rows:
        job_id = row["job_id"]
        status = live_status.get(job_id) or acct_status.get(job_id)
        if status:
            row["state"] = status.get("state", row["state"]) or row["state"]
            row["elapsed"] = status.get("elapsed", row["elapsed"]) or row["elapsed"]
            row["exit_code"] = status.get("exit_code", row["exit_code"]) or row["exit_code"]
            row["reason"] = status.get("reason", row["reason"]) or row["reason"]
        elif not use_slurm and row["state"] == "SUBMITTED":
            row["state"] = "RUNNING"
        updated_rows.append(row)

    write_job_records(project_config, updated_rows)
    return updated_rows


def summarize_jobs(rows):
    counts = Counter(row["state"] or "UNKNOWN" for row in rows)
    return counts


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
