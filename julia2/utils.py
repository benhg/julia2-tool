"""
Utility functions for the project
"""
import subprocess
import sys
import logging
import shutil
import os

logger = logging.getLogger("julia2.utils")


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
#SBATCH --partition {slurm_settings.partition_name}
#SBATCH --mail-user {slurm_settings.email}
#SBATCH --mail-type BEGIN
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
#SBATCH -e {out_dir}/slurm-%j.out
#SBATCH -o {out_dir}/slurm-%j.out
"""
    logging.debug(f"Running slurm job on node {node} with CPUs {cpus}")
    keyList = sorted(slurm_settings.nodes.keys())
    for i, v in enumerate(keyList):
        if v == node:
            if i + 1 >= len(keyList):
                i = 0
            slurm_settings.current_node = keyList[i + 1]
            break

    return sbatch_template, cpus


def run_slurm_job(sbatch_text, sbatch_name, project_config, system_config):
    """
    Submit a Slurm job with SBatch text passed in

    sbatch_name must not contain any /
    """
    logger.debug(f"Creating sbatch file at {project_config.project_dir}/slurm_jobs/{sbatch_name}.sh")
    with open(f"{project_config.project_dir}/slurm_jobs/{sbatch_name}.sh",
              "w") as fh:
        fh.write(sbatch_text)

    logging.debug(f"Command: [sbatch] {project_config.project_dir}/slurm_jobs/{sbatch_name}.sh")
    if system_config.use_slurm:
        logging.info(
        subprocess.check_output(
            f"sbatch {project_config.project_dir}/slurm_jobs/{sbatch_name}.sh",
            shell=True))
    else:
        logging.info(subprocess.check_output(
            f"/bin/bash {project_config.project_dir}/slurm_jobs/{sbatch_name}.sh &",
            shell=True))


def setup_logging(project_config, log_level: int = logging.ERROR):
    """
    Set up basic logging configuration.

    Args:
        log_file (str): The path to the log file.
        log_level (int): The logging level (e.g., logging.DEBUG, logging.INFO).
    # TODO: lower the log level to INFO when ready
    """
    # Define log format
    log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    date_format = "%Y-%m-%d %H:%M:%S"

    log_file = f"{project_config.project_dir}/logs/julia2.log"

    # Set up the root logger
    logging.basicConfig(
        level=log_level,
        format=log_format,
        datefmt=date_format,
        handlers=[
            logging.FileHandler(log_file),   # Log to a file
            logging.StreamHandler(sys.stdout)         # Log to console
        ]
    )

def move_with_exist_ok(src: str, dst: str, exist_ok: bool = False):
    """
    Move a file or directory with optional `exist_ok` behavior.

    Args:
        src (str): The source path.
        dst (str): The destination path.
        exist_ok (bool): If True, overwrite the destination if it exists.
                         If False, raise an error if the destination exists.
    """
    if os.path.exists(dst):
        if exist_ok:
            # Remove the existing destination if it's a file
            if os.path.isfile(dst) or os.path.islink(dst):
                os.remove(dst)
            # If it's a directory, remove it and its contents
            elif os.path.isdir(dst):
                raise FileExistsError(f"Destination '{dst}' already exists as dir. Not deleting directory")
                shutil.rmtree(dst)
        else:
            raise FileExistsError(f"Destination '{dst}' already exists.")
    
    # Perform the move
    shutil.move(src, dst)
