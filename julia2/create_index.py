"""
Make separate bowtie indexes for each sequence

They'll go into /home/labs/binford//home/labs/binford/taxon_confirmation_indexes/{NAME}/
"""

import subprocess
from Bio import SeqIO
import glob
import shutil

import utils


def split_fasta_file_into_indexes(new_fasta_path, project_config):
    """Turn each sequence in the fasta file into its own fasta file
    """
    with open(new_fasta_path, "r") as old_handle:
        sequences = SeqIO.parse(old_handle, "fasta")
        for record in sequences:
            with open(f"{project_config.project_path}/indexes/{record.id}.fasta", "w") as new_file:
                new_file.write(f">{record.id}\n")
                new_file.write(f"{record.seq}\n")


def submit_all_index_requests(project_config, slurm_settings):
    """
    Submit a whole bunch of SBatch scripts to create a lot of indexes

    Assumes that all .fasta files in the indexes directory should be split into indexes
    """
    base_dir = f"{project_config.project_path}/indexes/"
    files = glob.glob(f"{base_dir}/*.fasta")
    for file in files:
        name = file.split(".fasta")[0].split("/home/labs/binford/taxon_confirmation_indexes/")[1]
        create_index_slurm(name, project_config, slurm_settings)


def create_index_slurm(index_name, slurm_settings, project_config):
    """
    Create an individual index by submitting a SLURM command
    """
    sbatch_template, cpus = utils.create_sbatch_template(slurm_settings, project_config, cpus=True, align_index="INDEX")
    sbatch_commands_template = f"""
mkdir -p /home/labs/binford/taxon_confirmation_indexes/{index_name}_index
chmod 777 /home/labs/binford/taxon_confirmation_indexes/{index_name}_index

bowtie2-build --threads {cpus} /home/labs/binford/taxon_confirmation_indexes/{index_name}.fasta /home/labs/binford/taxon_confirmation_indexes/{index_name}_index/{index_name}_index
"""
    sbatch_text  = f"""
{sbatch_template}

{sbatch_commands_template}
"""
    utils.run_slurm_job(sbatch_text, f"index_{index_name}", project_config)


def cleanup_index_fastas(project_config):
    """
    Find all the FASTA files that already have indexes associated with them
    Move them into their index directory
    """
    base_dir = f"{project_config.project_path}/indexes/"
    files = glob.glob(f"{base_dir}/*.fasta")
    for file in files:
        dir_name = f"{file.split(".fasta")[0]}"
        if os.path.isdir(dir_name):
            shutil.move(file, dir_name)


def create_all_indexes_for_new_fasta(new_fasta_path):
    """
    Create all the indexes and add them to the index directory
    Requires a new FASTA labelled with index names
    """
    cleanup_index_fastas(project_config)
    split_fasta_file_into_indexes(new_fasta_path, project_config)
    submit_all_index_requests(project_config, slurm_settings)
