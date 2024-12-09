"""
Manage the data for a project
"""

import os

import logging

logger = logging.getLogger("julia2.manage_data")



def create_blank_project(project_config, system_config):
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

    # Project path should be handled by earlier steps when config is created
    project_dir = project_config.project_dir
    project_name = project_config.project_name

    if not os.path.exists(project_dir):
        os.makedirs(project_dir, exist_ok=True)
    else:
        res = input("WARNING: Project already exists. Continue anyways? y/N")
        if res.lower() != "y":
            return

    # Make project config
    with open(f"{project_dir}/project_config.json", "a"):
        os.utime(f"{project_dir}/project_config.json", None)

    os.makedirs(f"{project_dir}/indexes", exist_ok=True)
    os.makedirs(f"{project_dir}/raw_reads", exist_ok=True)
    os.makedirs(f"{project_dir}/slurm_jobs", exist_ok=True)
    os.makedirs(f"{project_dir}/logs", exist_ok=True)
    os.makedirs(f"{project_dir}/output", exist_ok=True)
    os.makedirs(f"{project_dir}/output/alignment_database_data",
                exist_ok=True)
    os.makedirs(f"{project_dir}/output/index_creation", exist_ok=True)

    # Make output files
    with open(f"{project_dir}/output/alignment_database.csv", "a"):
        os.utime(f"{project_dir}/output/alignment_database.csv", None)

    with open(f"{project_dir}/output/hoppingg_results.csv", "a"):
        os.utime(f"{project_dir}/output/alignment_database.csv", None)


def delete_project(project_config, system_config):
    """
    Delete files associated with a project
    """
    res = input(f"Deleting {project_config.project_name}. Are you sure? y/N")
    if res.lower() == "y":
        os.rmdirs(project_config.project_dir)
        system_config.projects[project_config.project_name] = None
