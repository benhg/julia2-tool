"""
Generate a Summary CSV based on the Bowtie alignment output(s)

Summary contains the information that Bowtie produced as well as where to get the full SAM file
"""
from glob import glob
import os
import csv
import subprocess
import sys
import logging
import multiprocessing
import fcntl

logger = logging.getLogger("julia2.database")

import julia2.utils as utils

headers = [
    "reads_sample", "reads_taxon", "index_sample", "index_taxon", "num_reads",
    "pairtype", "num_aligned_none", "num_aligned_once", "num_aligned_multiple",
    "none_alignment_rate", "single_alignment_rate", "multiple_alignment_rate",
    "num_aligned_any", "alignment_rate", "exec_time"
]

def combine_out_err_files(project_config):
    """
    In case SLURM out and error files get separated, we neeed to fold them all into the corresponding .out file
    """
    out_path = f"{project_config.project_dir}/output/alignment_database_data/slurm-*.*"
    # This should only ever return 2 files.
    # Combine the .out and the .err into .out, with the .out first

def update_database_single(vals_tuple):
    """
    Update the database for a single file as part of the multiprocessing Pool
    """
    file, output_file, data,sample_to_taxon,sample_to_taxon_short,row = vals_tuple
    with open(file) as fh2:
        try:
            # Name and metadata
            slurm_job_name = file.split("slurm-")[1].split(".out")[0]
            slurm_time_str = utils.run_cmd(
                f'sacct --format="Elapsed" -j {slurm_job_name}')
            slurm_time = slurm_time_str.split("\n")[-2].strip()

            # Text from stderr (which is summary info)
            data = fh2.readlines()
            data = [l for l in data if l.strip() != ""]

            #Sample IDs
            index_sample = data[0].split(" ")[0].split("_")[1].strip()
            index_sample_long = data[0].split(" ")[0]
            reads_sample = data[0].split(" ")[1].split("_")[1].strip()

            #print(
            #    f"Job {slurm_job_name} for index {index_sample} and reads {reads_sample} took {slurm_time}"
            #)
            # This is gonna be gross
            row = {
                "index_sample": index_sample_long,
                "index_taxon": sample_to_taxon[index_sample],
                "reads_sample": reads_sample,
                "reads_taxon": sample_to_taxon[reads_sample],
                "num_reads": int(data[1].split(" ")[0]),
                "num_aligned_none": int(data[3].split("(")[0].strip()),
                "num_aligned_once": int(data[4].split("(")[0].strip()),
                "num_aligned_multiple":
                int(data[5].split("(")[0].strip()),
                "exec_time": slurm_time
            }

            # Pair type
            try:
                if int(reads_sample.split("s")[1]) <= project_config.num_samples_per_lane - 1:
                    reads_lane = 1
                else:
                    reads_lane = 2
            except:
                reads_lane = 2

            try:
                if int(index_sample.split("s")[1]) <= project_config.num_samples_per_lane - 1:
                    index_lane = 1
                else:
                    index_lane = 2
            except:
                index_lane = 2

            if reads_sample == index_sample:
                row["pairtype"] = "True_Auto"
            elif sample_to_taxon_short[
                    reads_sample] == sample_to_taxon_short[
                        index_sample]:
                row["pairtype"] = "Taxon_Auto"
            elif reads_lane != index_lane:
                row["pairtype"] = "Other_Lane"
            else:
                row["pairtype"] = "Allo"

            # Alignment rates
            row["single_alignment_rate"] = row[
                "num_aligned_once"] / int(row["num_reads"])
            row["none_alignment_rate"] = row["num_aligned_none"] / row[
                "num_reads"]
            row["multiple_alignment_rate"] = row[
                "num_aligned_multiple"] / row["num_reads"]
            row["num_aligned_any"] = int(
                row["num_aligned_once"]) + int(
                    row["num_aligned_multiple"])
            row["alignment_rate"] = row["num_aligned_any"] / row[
                "num_reads"]
            with open(output_file, "a") as fh:
                fcntl.flock(g, fcntl.LOCK_EX)
                writer = csv.DictWriter(fh, fieldnames=headers)
                writer.writerow(row)
                fcntl.flock(g, fcntl.LOCK_UN)
        except Exception as e:
            if "list index out of range" in str(e):
                print(f"File {file} is still running")
                return

            print(f"failed for file {file}")

def update_database(project_config):
    """
    Update the intermediate alignment database checking for newly completed runs
    """

    path = f"{project_config.project_dir}/output/alignment_database_data/slurm-*.out"
    output_file = f"{project_config.project_dir}/output/alignment_database.csv"

    sample_to_taxon = project_config.sample_to_taxon
    sample_to_taxon_short = project_config.sample_to_taxon_short

    with open(output_file, "w") as fh:
        writer = csv.writer(fh)
        writer.writerow(headers)

    
    all_files = glob(path) 
    inputs = []
    for file in all_files:
        inputs.append((file, output_file, data,sample_to_taxon,sample_to_taxon_short,row))
   
    pool = multiprocessing.Pool()
    pool.map(update_database_single, inputs)

    while True:
        try:
            ready = [result.ready() for result in results]
            successful = [result.successful() for result in results]
            break
        except Exception:
            continue
       

    print(f"Generated database to {output_file}")
