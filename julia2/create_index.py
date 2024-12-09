"""
Make separate bowtie indexes for each sequence

They'll go into /home/labs/binford//home/labs/binford/taxon_confirmation_indexes/{NAME}/
"""

import subprocess
from Bio import SeqIO
import glob

import utils

base_dir = "/home/labs/binford/taxon_confirmation_indexes/"


def split_fasta_file_into_indexes():
    ## First, turn each sequence in the fasta file into its own fasta file
    with open("data/all_sequences.fasta", "r") as old_handle:
        sequences = SeqIO.parse(old_handle, "fasta")
        for record in sequences:
            with open(f"{base_dir}/{record.id}.fasta", "w") as new_file:
                new_file.write(f">{record.id}\n")
                new_file.write(f"{record.seq}\n")



def submit_all_index_requests():
    """
    Submit a whole bunch of SBatch scripts to create a lot of indexes
    """
    ## Then, submit a bunch of indexing jobs to make indexes
    files = glob.glob(f"{base_dir}/*.fasta")
    for file in files:
        name = file.split(".fasta")[0].split("/home/labs/binford/taxon_confirmation_indexes/")[1]
        sbatch_text = sbatch_template.format(name, name, name,
                                             name, name)
        with open(f"bowtie_cmds/gen_index_s{name}.sh", "w") as fh:
            fh.write(sbatch_text)
        print(
            subprocess.check_output(
                f"sbatch /home/glick/JULIA-Take-2/src/single_sequence_index/bowtie_cmds/gen_index_s{name}.sh",
               shell=True))


def create_index_slurm():
    """
    Create an individual index by submitting a SLURM command
    """

    sbatch_commands = """
mkdir -p /home/labs/binford/taxon_confirmation_indexes/{}_index
chmod 777 /home/labs/binford/taxon_confirmation_indexes/{}_index

bowtie2-build --threads {} /home/labs/binford/taxon_confirmation_indexes/{}.fasta /home/labs/binford/taxon_confirmation_indexes/{}_index/{}_index
"""

    sbatch_text = create_sbatch_template(sbatch_commands, slurm_settings)
    run_slurm_job(sbatch_text, sbatch_name, project_config)


def cleanup_index_fastas():
    """
    Find all the FASTA files that already have indexes associated with them
    Move them into their index directory
    """



    """
    Create a new directory structure for a new project

    Takes project config (for project path and project name) as input
    
    File layout:
    .
    ├── project_config.json
    ├── indexes
    │   ├── seq1_index
    │   │   ├── seq1_index.btl2
    │   │   └── seq1_index.fasta
    │   └── seq2_index
    │       ├── seq2_index.btl2
    │       └── seq2_index.fasta
    ├── logs
    │   └── julia2.log
    ├── output
    │   ├── alignment_database.csv
    │   ├── alignment_database_data
    │   │   ├── index_01_sample_01.out
    │   │   ├── index_01_sample_02.out
    │   │   ├── index_02_sample_01.out
    │   │   └── index_02_sample_02.out
    │   └── hopping_results.csv
    ├── raw_reads
    │   ├── sample1.fasta
    │   └── sample2.fasta
    └── slurm_jobs
        ├── align1.sh
        ├── align2.sh
        ├── index.sh
        └── index1.sh 
    """