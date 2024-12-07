"""
Make separate bowtie indexes for each sequence

They'll go into /home/labs/binford//home/labs/binford/taxon_confirmation_indexes/{NAME}/
"""

import subprocess
from Bio import SeqIO
import glob

base_dir = "/home/labs/binford/taxon_confirmation_indexes/"

sbatch_template = """#!/bin/bash
#SBATCH --cpus-per-task=48

mkdir -p /home/labs/binford/taxon_confirmation_indexes/{}_index
chmod 777 /home/labs/binford/taxon_confirmation_indexes/{}_index

bowtie2-build --threads 48 /home/labs/binford/taxon_confirmation_indexes/{}.fasta /home/labs/binford/taxon_confirmation_indexes/{}_index/{}_index
"""

## First, turn each sequence in the fasta file into its own fasta file
with open("data/all_sequences.fasta", "r") as old_handle:
    sequences = SeqIO.parse(old_handle, "fasta")
    for record in sequences:
        with open(f"{base_dir}/{record.id}.fasta", "w") as new_file:
            new_file.write(f">{record.id}\n")
            new_file.write(f"{record.seq}\n")

    


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

