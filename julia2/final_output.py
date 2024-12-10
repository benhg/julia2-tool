"""
Postprocess the summary to extract and apply the thresholding we talked about
"""
from glob import glob
import csv
import subprocess
import sys
import logging

logger = logging.getLogger("julia2.final_output")

import julia2.utils as utils

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


def create_final_output(project_config):

    input_file = f"{project_config.project_dir}/output/alignment_database.csv"
    output_file = f"{project_config.project_dir}/output/hopping_results.csv"

    sample_to_taxon = project_config.sample_to_taxon
    sample_to_taxon = project_config.sample_to_taxon_short

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
        "sequence_name", "hopper_status", "source_index_taxon_labelled",
        "total_reads_mapped", "read_max_sample", "read_max_taxon",
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
        aligned_num_auto = 0
        aligned_taxon_auto = ""
        for row in group:
            total_single_aligned += int(row["num_aligned_once"])
            if int(row["num_aligned_once"]) >= max_aligned_num:
                max_aligned_num = int(row["num_aligned_once"])
                max_aligned_read = row["reads_sample"]
                max_aligned_taxon = row["reads_taxon"]
            if sequence_sample_label in row["reads_sample"]:
                #print(row["reads_sample"])
                aligned_num_auto = int(row["num_aligned_once"])
                aligned_taxon_auto = row["reads_taxon"]
        percent_reads_from_max = (max_aligned_num / total_single_aligned) * 100
        percent_reads_from_auto = (aligned_num_auto /
                                   total_single_aligned) * 100
        threshold_metric = (percent_reads_from_auto /
                            percent_reads_from_max) * 100
        summary_row = {
            "sequence_name": sequence,
            "source_index_taxon_labelled":
            sample_to_taxon[sequence_sample_label],
            "total_reads_mapped": total_single_aligned,
            "read_max_sample": max_aligned_read,
            "read_max_taxon": max_aligned_taxon,
            "percent_reads_from_max": percent_reads_from_max,
            "num_reads_from_labelled_true_auto": aligned_num_auto,
            "percent_reads_from_labelled_true_auto": percent_reads_from_auto,
            "threshold_metric": threshold_metric,
            "hopper_status": hopper_threshold(threshold_metric)
        }
        writer.writerow(summary_row)
    print(f"Generated output to {output_file}")
