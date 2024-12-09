"""
Generate a Summary CSV based on the Bowtie alignment output(s)

Summary contains the information that Bowtie produced as well as where to get the full SAM file
"""
from glob import glob
import csv
import subprocess
import sys

headers = [
    "reads_sample", "reads_taxon", "index_sample", "index_taxon", "num_reads",
    "pairtype", "num_aligned_none", "num_aligned_once", "num_aligned_multiple",
    "none_alignment_rate", "single_alignment_rate", "multiple_alignment_rate",
    "num_aligned_any", "alignment_rate", "exec_time"
]


def update_database(project_config):
    """
    Update the intermediate alignment database checking for newly completed runs
    """

    path = f"{project_config.project_path}/output/alignment_database_data/slurm-*.out"
    output_file = f"{project_config.project_path}/output/alignment_database_data/alignment_database.csv"

    sample_to_taxon = project_config.sample_to_taxon
    sample_to_taxon = project_config.sample_to_taxon_short

    with open(output_file, "w") as fh:
        writer = csv.writer(fh)
        writer.writerow(headers)

    with open(output_file, "a") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers)
        all_files = glob(path)
        for file in all_files:
            with open(file) as fh2:
                try:
                    # Name and metadata
                    slurm_job_name = file.split("slurm-")[1].split(".out")[0]
                    slurm_time_str = run_cmd(
                        f'sacct --format="Elapsed" -j {slurm_job_name}')
                    slurm_time = slurm_time_str.split("\n")[-2].strip()

                    # Text from stderr (which is summary info)
                    data = fh2.readlines()

                    #Sample IDs
                    index_sample = data[0].split(" ")[0].split("_")[1].strip()
                    index_sample_long = data[0].split(" ")[0]
                    reads_sample = data[0].split(" ")[1].split("_")[1].strip()

                    print(
                        f"Job {slurm_job_name} for index {index_sample} and reads {reads_sample} took {slurm_time}"
                    )
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
                        if int(reads_sample.split("s")[1]) <= 11:
                            reads_lane = 1
                        else:
                            reads_lane = 2
                    except:
                        reads_lane = 2

                    try:
                        if int(index_sample.split("s")[1]) <= 11:
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

                    writer.writerow(row)
                except Exception as e:

                    if "list index out of range":
                        print(f"File {file} is still running")
                        continue

                    print(f"failed for file {file}")
