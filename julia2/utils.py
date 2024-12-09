"""
Utility functions for the project
"""
import subprocess

def run_cmd(cmd):
    return subprocess.check_output(cmd, shell=True).decode(sys.stdout.encoding)

def convert_to_fasta(input_file, output_file):
    """
    Convert gdoc text to FASTA

    Converts the google docs I was given into more usable fasta
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        label = None
        sequence = []
        in_block = False  # Flag to indicate we are in a block

        for line in infile:
            line = line.strip()

            if line.startswith(
                    ">"
            ) and not in_block:  # Only capture the first label in each block
                if label and sequence:  # If a sequence was already collected, write it to file
                    outfile.write(f">{label}\n")
                    outfile.write(''.join(sequence) + "\n")
                label = line[1:]  # Store the first label (without ">")
                sequence = []  # Reset the sequence
                in_block = True  # Indicate we are processing a block

            elif line.startswith(">") and in_block:
                continue  # Skip additional labels in the same block

            elif line.startswith("#"):  # Skip comment lines
                continue

            elif line:  # If it's a sequence line, add to the sequence
                sequence.append(line)

            else:  # Blank line or end of block, reset for the next block
                in_block = False

        # Write the last sequence to the file
        if label and sequence:
            outfile.write(f">{label}\n")
            outfile.write(''.join(sequence) + "\n")


def create_sbatch_template(slurm_settings,
                           project_config,
                           cpus=True,
                           align_index="ALIGN"):
    """
    Create an SBatch file - return the text and the number of CPUs
    """

    node = slurm_settings.current_node

    # Use all CPUs by default
    if cpus == True:
        cpus = slurm_settings.nodes[node]

    if align_index == "ALIGN":
        out_dir = f"{project_config.project_dir}/output/alignment_database_data"
    else:
        out_dir = f"{project_config.project_dir}/output/index_creation"

    sbatch_template = f"""#!/bin/bash
#SBATCH --cpus-per-task={cpus}
#SBATCH --partition {slurm_settings.partition}
#SBATCH --mail-user {slurm_settings.email}
#SBATCH --mail-type BEGIN
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
#SBATCH -e {out_dir}/slurm-%j.err
#SBATCH -o {out_dir}/slurm-%j.out
"""

    keyList = sorted(slurm_settings.nodes.keys())
    for i, v in enumerate(keyList):
        if v == node:
            if i + 1 >= len(keyList):
                i = 0
            slurm_settings.current_node = slurm_settings.nodes[keyList[i + 1]]
            break

    return sbatch_template, cpus


def run_slurm_job(sbatch_text, sbatch_name, project_config):
    """
    Submit a Slurm job with SBatch text passed in
    """
    with open(f"{project_config.project_dir}/slurm_jobs/{sbatch_name}.sh",
              "w") as fh:
        fh.write(sbatch_text)
    print(
        subprocess.check_output(
            f"sbatch {project_config.project_dir}/slurm_jobs/{sbatch_name}.sh",
            shell=True))
