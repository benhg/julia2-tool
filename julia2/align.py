"""
Run alignments on the small indexes

They'll go into /home/labs/binford/taxon_confirmation_indexes/s_{XYZ}/
"""

import subprocess
from glob import glob
import time

import utils

def run_alignment(reads_sample_id, index_id, slurm_settings, project_config):
    # TODO: Handle more than 2 lanes
    if int(reads_sample_id) <= num_samples_per_lane - 1:
        lane = 1
    else:
        lane = 2

    sbatch_template, cpus = utils.create_sbatch_template(slurm_settings,
                                                         project_config,
                                                         cpus=True,
                                                         align_index="ALIGN")

    print(lane, f"{project_config.project_dir}/raw_reads/lane{lane}-s{reads_sample_id}*R1*")

    #print(f"{combined_files_dir}/lane{lane}-s{index_id}*R1*")
    dir_1_filename = glob(
        f"{project_config.project_dir}/raw_reads/lane{lane}-s{reads_sample_id}*R1*")[0]

    sbatch_cmds = f"""
echo "index_{index_id} read_s{reads_sample_id}"

bowtie2 -f --threads {cpus} -x /home/labs/binford/taxon_confirmation_indexes/{index_id}_index/{index_id}_index -U {dir_1_filename} > {project_config.project_dir}/output/alignment_database_data/raw/index_{index_id}_read_s{reads_sample_id}.sam
"""

    sbatch_text = f"""
{sbatch_template}

{sbatch_cmds}
"""

    print(dir_1_filename, index_id, reads_sample_id)
    utils.run_slurm_job(sbatch_text,
                        f"align_index_{index_id}_reads_{reads_sample_id}",
                        project_config)


def run_all_samples(slurm_settings, project_config, sequence_name_list):
    # For each sequence
    sequences = open(sequence_name_list).readlines()
    for index in sequences:
        index_id = index.split(".fasta")[0].split(
            "/home/labs/binford/taxon_confirmation_indexes/")[1]
        # For each sample
        for i in range(1, project_config.num_samples + 1):
            reads_sample_id = str(i).zfill(3)
            print(reads_sample_id, index_id)
            run_alignment(reads_sample_id, index_id, slurm_settings,
                          project_config)

def run_all_true_auto_samples(slurm_settings, project_config, sequence_name_list):
    sequences = open(sequence_name_list).readlines()
    for index in sequences:
        index_id = index.split(".fasta")[0].split(
            "/home/labs/binford/taxon_confirmation_indexes/")[1]
        # For each sample
        for i in range(1, project_config.num_samples + 1):
            reads_sample_id = str(i).zfill(3)
            if reads_sample_id in index_id:
                print(reads_sample_id, index_id)
                run_alignment(reads_sample_id, index_id, slurm_settings,
                              project_config)

def run_all_allo_samples(slurm_settings, project_config, sequence_name_list):
    sequences = open(sequence_name_list).readlines()
    for index in sequences:
        index_id = index.split(".fasta")[0].split(
            "/home/labs/binford/taxon_confirmation_indexes/")[1]
        # For each sample
        for i in range(1, project_config.num_samples + 1):
            reads_sample_id = str(i).zfill(3)
            if int(reads_sample_id) not in int(index_id):
                print(reads_sample_id, index_id)
                run_alignment(reads_sample_id, index_id, slurm_settings,
                              project_config)

def run_all_taxon_auto_samples(slurm_settings, project_config, sequence_name_list):
    sequences = open(sequence_name_list).readlines()
    for index in sequences:
        index_id = index.split(".fasta")[0].split(
            "/home/labs/binford/taxon_confirmation_indexes/")[1]
        # For each sample
        for i in range(1, project_config.num_samples + 1):
            reads_sample_id = str(i).zfill(3)
            index_sample_id = index_id.split("_")[0]
            # TODO: Remove HACK
            if "c74742" in index_sample_id:
                index_sample_id = "s020"
            if "c49446" in index_sample_id:
                index_sample_id = "s018"
            if project_config.index_to_taxon_short[int(index_sample_id)] == project_config.index_to_taxon_short[int(reads_sample_id)]:
                print(reads_sample_id, index_id)
                run_alignment(reads_sample_id, index_id, slurm_settings,
                              project_config)

def run_all_intra_lane_samples(slurm_settings, project_config, sequence_name_list):
    sequences = open(sequence_name_list).readlines()
    for index in sequences:
        index_id = index.split(".fasta")[0].split(
            "/home/labs/binford/taxon_confirmation_indexes/")[1]
        # For each sample
        for i in range(1, project_config.num_samples + 1):
            reads_sample_id = str(i).zfill(3)
            index_sample_id = index_id.split("_")[0]
            # TODO: Remove HACK
            if "c74742" in index_sample_id:
                index_sample_id = "s020"
            if "c49446" in index_sample_id:
                index_id = "s018"

            # Pair type
            if int(reads_sample.split("s")[1]) <= project_config.num_samples_per_lane - 1:
                reads_lane = 1
            else:
                reads_lane = 2
            if int(index_sample.split("s")[1]) <= project_config.num_samples_per_lane - 1:
                index_lane = 1
            else:
                index_lane = 2


            if index_lane == reads_lane:
                print(reads_sample_id, index_id)
                run_alignment(reads_sample_id, index_id, slurm_settings,
                              project_config)

def run_all_cross_lane_samples(slurm_settings, project_config, sequence_name_list):
    sequences = open(sequence_name_list).readlines()
    for index in sequences:
        index_id = index.split(".fasta")[0].split(
            "/home/labs/binford/taxon_confirmation_indexes/")[1]
        # For each sample
        for i in range(1, project_config.num_samples + 1):
            reads_sample_id = str(i).zfill(3)
            index_sample_id = index_id.split("_")[0]
            # TODO: Remove HACK
            if "c74742" in index_sample_id:
                index_sample_id = "s020"
            if "c49446" in index_sample_id:
                index_id = "s018"

            # Pair type
            if int(reads_sample.split("s")[1]) <= project_config.num_samples_per_lane - 1:
                reads_lane = 1
            else:
                reads_lane = 2
            if int(index_sample.split("s")[1]) <= project_config.num_samples_per_lane - 1:
                index_lane = 1
            else:
                index_lane = 2


            if index_lane != reads_lane:
                print(reads_sample_id, index_id)
                run_alignment(reads_sample_id, index_id, slurm_settings,
                              project_config)