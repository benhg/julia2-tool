"""
Run alignments on the small indexes

They'll go into /home/labs/binford/taxon_confirmation_indexes/s_{XYZ}/
"""

import subprocess
from glob import glob
import time

combined_files_dir = "/home/labs/binford/raw_reads_fasta_tagged_batched/combined_files/"
base_dir = "/home/labs/binford/taxon_confirmation_indexes/"

sbatch_template = """#!/bin/bash
#SBATCH --cpus-per-task=48

echo "index_{} read_s{}"

bowtie2 -f --threads 48 -x /home/labs/binford/taxon_confirmation_indexes/{}_index/{}_index -U {} > /home/labs/binford/taxon_confirmation_indexes/{}_index/index_{}_read_s{}.sam
"""


def run_alignment(reads_sample_id, index_id):
    if int(reads_sample_id) <= 11:
        lane = 1
    else:
        lane = 2

    print(lane, f"{combined_files_dir}/lane{lane}-s{reads_sample_id}*R1*")

    #print(f"{combined_files_dir}/lane{lane}-s{index_id}*R1*")
    dir_1_filename = glob(
        f"{combined_files_dir}/lane{lane}-s{reads_sample_id}*R1*")[0]

    job_filename = f"bowtie_cmds/unpaired_align_{index_id}_s{reads_sample_id}.sh"

    sbatch_text = sbatch_template.format(index_id, reads_sample_id, index_id,
                                         index_id, dir_1_filename, index_id,
                                         index_id, reads_sample_id)
    print(dir_1_filename, index_id, reads_sample_id)
    with open(job_filename, "w") as fh:
        fh.write(sbatch_text)
    time.sleep(0.1)
    print(subprocess.check_output(f"sbatch {job_filename}", shell=True))


def run_all_intra_lane_samples():
    # For each sequence
    for index in glob(f"{base_dir}/*.fasta"):
        index_id = index.split(".fasta")[0].split("/home/labs/binford/taxon_confirmation_indexes/")[1]
        # For each sample
        for i in range(1,23):
            reads_sample_id = str(i).zfill(3)
            print(reads_sample_id, index_id)
            run_alignment(reads_sample_id, index_id)

run_all_intra_lane_samples()
