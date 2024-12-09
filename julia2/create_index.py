"""
Make separate bowtie indexes for each sequence

They'll go into /home/labs/binford//home/labs/binford/taxon_confirmation_indexes/{NAME}/
"""

import subprocess
from Bio import SeqIO
import glob

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


def submit_all_index_requests():
    """
    Submit a whole bunch of SBatch scripts to create a lot of indexes
    """
    ## Then, submit a bunch of indexing jobs to make indexes
        sbatch_commands_template = """
mkdir -p /home/labs/binford/taxon_confirmation_indexes/{}_index
chmod 777 /home/labs/binford/taxon_confirmation_indexes/{}_index

bowtie2-build --threads {} /home/labs/binford/taxon_confirmation_indexes/{}.fasta /home/labs/binford/taxon_confirmation_indexes/{}_index/{}_index
"""
    base_dir = f"{project_config.project_path}/indexes/"
    files = glob.glob(f"{base_dir}/*.fasta")
    for file in files:
        name = file.split(".fasta")[0].split("/home/labs/binford/taxon_confirmation_indexes/")[1]
        sbatch_commands = sbatch_commands_template.format(name, name, name,
                                             name, name)
        create_index_slurm(sbatch_commands)


def create_index_slurm(sbatch_commands):
    """
    Create an individual index by submitting a SLURM command
    """
    sbatch_text = utils.create_sbatch_template(sbatch_commands, slurm_settings)
    utils.run_slurm_job(sbatch_text, sbatch_name, project_config)


def cleanup_index_fastas():
    """
    Find all the FASTA files that already have indexes associated with them
    Move them into their index directory
    """
    base_dir = f"{project_config.project_path}/indexes/"
    files = glob.glob(f"{base_dir}/*.fasta")
    for file in files:

def create_all_indexes_for_new_fasta(new_fasta_path):
    """
    Create all the indexes and add them to the index directory
    Requires a new FASTA labelled with index names
    """
    cleanup_index_fastas()
    split_fasta_file_into_indexes(new_fasta_path)
    submit_all_index_requests()
