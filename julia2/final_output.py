"""
Postprocess the summary to extract and apply the thresholding we talked about
"""
import csv
import logging
import os

logger = logging.getLogger("julia2.final_output")

import julia2.job_tracking as job_tracking

# Database headers list is only stored in database.py file
from julia2.database import headers


def hopper_threshold(value):
    """
    Take a value, and return the "hopperness" associated with this value
    """
    if value <= 1:
        return "LIKELY_HOPPER"
    elif value < 20:
        return "MAYBE_HOPPER"
    else:
        return "NOT_HOPPER"


def _sequence_fasta_paths(project_config, sequence_name):
    canonical_name = _canonical_sequence_name(sequence_name)
    return [
        f"{project_config.project_dir}/indexes/{canonical_name}.fasta",
        f"{project_config.project_dir}/indexes/{canonical_name}_index/{canonical_name}.fasta",
    ]


def _canonical_sequence_name(sequence_name):
    if sequence_name.startswith("index_"):
        return sequence_name[len("index_"):]
    return sequence_name


def _sequence_length(project_config, sequence_name):
    for fasta_path in _sequence_fasta_paths(project_config, sequence_name):
        if not os.path.exists(fasta_path):
            continue
        with open(fasta_path, "r") as fh:
            return sum(len(line.strip()) for line in fh if line.strip() and not line.startswith(">"))
    return ""


def _latest_alignment_rows(project_config):
    latest = {}
    try:
        rows = job_tracking.load_job_records(project_config)
    except FileNotFoundError:
        return latest

    for row in rows:
        if row["job_type"] != "alignment":
            continue
        latest[(row["index_id"], row["reads_sample_id"])] = row
    return latest


def _tracking_summary(latest_alignment_rows, sequence_name):
    canonical_name = _canonical_sequence_name(sequence_name)
    matching_rows = [
        row for (index_id, _reads_sample_id), row in latest_alignment_rows.items()
        if index_id == canonical_name
    ]
    states = [row["state"] or "UNKNOWN" for row in matching_rows]
    completed = sum(1 for state in states if state == "COMPLETED")
    active = sum(
        1 for state in states
        if state in job_tracking.ACTIVE_STATES or state == "SUBMITTED")
    failed = sum(1 for state in states if state in job_tracking.FAILED_STATES)
    return {
        "tracked_alignment_jobs": len(matching_rows),
        "completed_alignment_jobs": completed,
        "active_alignment_jobs": active,
        "failed_alignment_jobs": failed,
    }


def create_final_output(project_config):

    input_file = f"{project_config.project_dir}/output/alignment_database.csv"
    output_file = f"{project_config.project_dir}/output/hopping_results.csv"

    sample_to_taxon = project_config.sample_to_taxon
    sample_to_taxon = project_config.sample_to_taxon_short
    latest_alignment_rows = _latest_alignment_rows(project_config)

    index_to_rows = {}

    f = open(input_file, 'r')
    reader = csv.DictReader(f, fieldnames=headers)

    for row in reader:
        existing_entry = index_to_rows.get(row["index_sample"], [])
        existing_entry.append(row)

        # Ignore the header row.
        if row["reads_sample"] != "reads_sample":
            index_to_rows[row["index_sample"]] = existing_entry

    out_headers = [
        "sequence_name", "source_index_sample_label", "source_index_taxon_labelled",
        "sequence_length", "tracked_alignment_jobs", "completed_alignment_jobs",
        "active_alignment_jobs", "failed_alignment_jobs", "hopper_status",
        "total_reads_mapped", "read_max_sample", "read_max_taxon", "pairtype_max_sample",
        "num_reads_max_sample",
        "percent_reads_from_max", "num_reads_from_labelled_true_auto",
        "percent_reads_from_labelled_true_auto", "threshold_metric"
    ]

    hdr_fh = open(output_file, 'w')
    header_writer = csv.writer(hdr_fh)
    header_writer.writerow(out_headers)
    hdr_fh.close()

    out_file = open(output_file, 'a')
    writer = csv.DictWriter(out_file, fieldnames=out_headers)

    for sequence, group in index_to_rows.items():
        sequence_sample_label = sequence.split("_")[1]
        #print(sequence_sample_label)
        total_single_aligned = 0
        max_aligned_num = 0
        max_aligned_read = ""
        max_aligned_taxon = ""
        max_aligned_pairtype = ""
        aligned_num_auto = 0
        aligned_taxon_auto = ""
        for row in group:
            total_single_aligned += int(row["num_aligned_once"])
            if int(row["num_aligned_once"]) >= max_aligned_num:
                max_aligned_num = int(row["num_aligned_once"])
                max_aligned_read = row["reads_sample"]
                max_aligned_taxon = row["reads_taxon"]
                max_aligned_pairtype = row["pairtype"]
            if sequence_sample_label in row["reads_sample"]:
                #print(row["reads_sample"])
                aligned_num_auto = int(row["num_aligned_once"])
                aligned_taxon_auto = row["reads_taxon"]
        if total_single_aligned == 0:
            percent_reads_from_max = 0
            percent_reads_from_auto = 0
            threshold_metric = 0
        else:
            percent_reads_from_max = (max_aligned_num / total_single_aligned) * 100
            percent_reads_from_auto = (aligned_num_auto /
                                       total_single_aligned) * 100
            if percent_reads_from_max == 0:
                threshold_metric = 0
            else:
                threshold_metric = (percent_reads_from_auto /
                                    percent_reads_from_max) * 100
        tracking_summary = _tracking_summary(latest_alignment_rows, sequence)
        summary_row = {
            "sequence_name": sequence,
            "source_index_sample_label": sequence_sample_label,
            "source_index_taxon_labelled":
            sample_to_taxon[sequence_sample_label],
            "sequence_length": _sequence_length(project_config, sequence),
            "tracked_alignment_jobs": tracking_summary["tracked_alignment_jobs"],
            "completed_alignment_jobs": tracking_summary["completed_alignment_jobs"],
            "active_alignment_jobs": tracking_summary["active_alignment_jobs"],
            "failed_alignment_jobs": tracking_summary["failed_alignment_jobs"],
            "total_reads_mapped": total_single_aligned,
            "read_max_sample": max_aligned_read,
            "read_max_taxon": max_aligned_taxon,
            "pairtype_max_sample": max_aligned_pairtype,
            "num_reads_max_sample": max_aligned_num,
            "percent_reads_from_max": percent_reads_from_max,
            "num_reads_from_labelled_true_auto": aligned_num_auto,
            "percent_reads_from_labelled_true_auto": percent_reads_from_auto,
            "threshold_metric": threshold_metric,
            "hopper_status": hopper_threshold(threshold_metric)
        }
        writer.writerow(summary_row)
    print(f"Generated output to {output_file}")
