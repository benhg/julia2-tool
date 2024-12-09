"""
Run alignments on the small indexes

They'll go into /home/labs/binford/taxon_confirmation_indexes/s_{XYZ}/
"""

import subprocess
from glob import glob
import time

import utils

combined_files_dir = "/home/labs/binford/raw_reads_fasta_tagged_batched/combined_files/"
base_dir = "/home/labs/binford/taxon_confirmation_indexes/"

def run_alignment(reads_sample_id, index_id, slurm_settings, project_config):
    # TODO: Handle more than 2 lanes
    if int(reads_sample_id) <= num_samples_per_lane - 1:
        lane = 1
    else:
        lane = 2

    sbatch_template, cpus = utils.create_sbatch_template(slurm_settings, project_config, cpus=True, align_index="ALIGN")

    print(lane, f"{combined_files_dir}/lane{lane}-s{reads_sample_id}*R1*")

    #print(f"{combined_files_dir}/lane{lane}-s{index_id}*R1*")
    dir_1_filename = glob(
        f"{combined_files_dir}/lane{lane}-s{reads_sample_id}*R1*")[0]

    sbatch_cmds = f"""
echo "index_{index_id} read_s{reads_sample_id}"

bowtie2 -f --threads {cpus} -x /home/labs/binford/taxon_confirmation_indexes/{index_id}_index/{index_id}_index -U {dir_1_filename} > {project_config.project_dir}/output/alignment_database_data/raw/index_{index_id}_read_s{reads_sample_id}.sam
"""
    
    sbatch_text = f"""
{sbatch_template}

{sbatch_cmds}
"""

    print(dir_1_filename, index_id, reads_sample_id)
    utils.run_slurm_job(sbatch_text, f"align_index_{index_id}_reads_{reads_sample_id}", project_config)



def run_all_intra_lane_samples(slurm_settings, project_config):
    # For each sequence
    for index in glob(f"{base_dir}/*.fasta"):
        index_id = index.split(".fasta")[0].split("/home/labs/binford/taxon_confirmation_indexes/")[1]
        # For each sample
        for i in range(1,project_config.num_samples + 1):
            reads_sample_id = str(i).zfill(3)
            print(reads_sample_id, index_id)
            run_alignment(reads_sample_id, index_id, slurm_settings, project_config)

run_all_intra_lane_samples()
