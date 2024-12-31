"""
Run alignments on the small indexes

They'll go into the indexes directory for the project
"""

import subprocess
from glob import glob
import time
import logging
import os

logger = logging.getLogger("julia2.align")

import julia2.utils as utils

def run_alignment(reads_sample_id, index_id, system_config, project_config):
    # TODO: Handle more than 2 lanes
    if int(reads_sample_id) <= project_config.num_samples_per_lane - 1:
        lane = 1
    else:
        lane = 2
    logging.debug(f"Running alignment for reads {reads_sample_id} and index {index_id}")
    sbatch_template, cpus = utils.create_sbatch_template(system_config.slurm_settings,
                                                         project_config,
                                                         cpus=True,
                                                         align_index="ALIGN")
    dir_1_filename = glob(
        f"{project_config.project_dir}/raw_reads/lane{lane}-s{reads_sample_id}*R1*")[0]

    sbatch_cmds = f"""
echo "index_{index_id} read_s{reads_sample_id}"

bowtie2 -f --threads {cpus} -x {project_config.project_dir}/indexes/{index_id}_index/{index_id}_index -U {dir_1_filename} > {project_config.project_dir}/output/alignment_database_data/raw/index_{index_id}_read_s{reads_sample_id}.sam
"""

    sbatch_text = f"""{sbatch_template}

{sbatch_cmds}
"""
    utils.run_slurm_job(sbatch_text,
                        f"align_index_{index_id}_reads_{reads_sample_id}",
                        project_config, system_config)

def run_all_samples(system_config, project_config, sequence_name_list):
    # For each sequence
    sequences = open(sequence_name_list).readlines()
    for index in sequences:
        index_id = os.path.basename(index.split(".fasta")[0]).strip()
        logging.info(f"Running all alignments on sample {index_id}")
        # For each sample
        for i in range(1, project_config.num_samples + 1):
            reads_sample_id = str(i).zfill(3)
            run_alignment(reads_sample_id, index_id, system_config,
                          project_config)

def run_all_true_auto_samples(system_config, project_config, sequence_name_list):
    sequences = open(sequence_name_list).readlines()
    for index in sequences:
        index_id = os.path.basename(index.split(".fasta")[0]).strip()
        logging.info(f"Running true-auto alignment on sample {index_id}")
        # For each sample
        for i in range(1, project_config.num_samples + 1):
            reads_sample_id = str(i).zfill(3)
            if f"s{reads_sample_id}_" in index_id:
                run_alignment(reads_sample_id, index_id, system_config,
                              project_config)
                logging.debug(f"Running Job for {reads_sample_id}, {index_id}")


def run_all_allo_samples(system_config, project_config, sequence_name_list):
    sequences = open(sequence_name_list).readlines()
    for index in sequences:
        index_id = os.path.basename(index.split(".fasta")[0]).strip()
        logging.debug(f"Running allo alignment on sample {index_id}")
        # For each sample
        for i in range(1, project_config.num_samples + 1):
            reads_sample_id = str(i).zfill(3)
            index_sample_id = index_id.split("_")[0]
            if project_config.sample_to_taxon_short[index_sample_id] != project_config.sample_to_taxon_short[f"s{reads_sample_id}"]:
                run_alignment(reads_sample_id, index_id, system_config,
                              project_config)

def run_all_taxon_auto_samples(system_config, project_config, sequence_name_list):
    sequences = open(sequence_name_list).readlines()
    for index in sequences:
        index_id = os.path.basename(index.split(".fasta")[0]).strip()
        logging.info(f"Running taxon-auto alignment on sample {index_id}")
        # For each sample
        for i in range(1, project_config.num_samples + 1):
            reads_sample_id = str(i).zfill(3)
            index_sample_id = index_id.split("_")[0]
            # TODO: Remove HACK
            if "c74742" in index_sample_id:
                index_sample_id = "s020"
            if "c49446" in index_sample_id:
                index_sample_id = "s018"
            if project_config.sample_to_taxon_short[index_sample_id] == project_config.sample_to_taxon_short[f"s{reads_sample_id}"]:
                run_alignment(reads_sample_id, index_id, system_config,
                              project_config)
                logging.debug(f"Running Job for {reads_sample_id}, {index_sample_id}, {project_config.sample_to_taxon_short[index_sample_id]}, {project_config.sample_to_taxon_short['s'+str(reads_sample_id)] }")

def run_all_intra_lane_samples(system_config, project_config, sequence_name_list):

    sequences = open(sequence_name_list).readlines()
    for index in sequences:
        index_id = os.path.basename(index.split(".fasta")[0]).strip()
        logging.info(f"Running same-lane alignment on sample {index_id}")
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
            if int(reads_sample_id) <= project_config.num_samples_per_lane - 1:
                reads_lane = 1
            else:
                reads_lane = 2
            if int(index_id.split("s")[1].split("_")[0]) <= project_config.num_samples_per_lane - 1:
                index_lane = 1
            else:
                index_lane = 2


            if index_lane == reads_lane:
                run_alignment(reads_sample_id, index_id, system_config,
                              project_config)

def run_all_cross_lane_samples(system_config, project_config, sequence_name_list):
    sequences = open(sequence_name_list).readlines()
    for index in sequences:
        index_id = os.path.basename(index.split(".fasta")[0]).strip()
        logging.info(f"Running other-lane alignment on sample {index_id}")
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
            if int(reads_sample_id) <= project_config.num_samples_per_lane - 1:
                reads_lane = 1
            else:
                reads_lane = 2
            if int(index_id.split("s")[1].split("_")[0]) <= project_config.num_samples_per_lane - 1:
                index_lane = 1
            else:
                index_lane = 2


            if index_lane != reads_lane:
                run_alignment(reads_sample_id, index_id, system_config,
                              project_config)
