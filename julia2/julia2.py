#!/usr/bin/env python3
"""
julia2.py

The `main` executable for this project.
"""
import argparse
import os

import config
import manage_data
import create_index
import final_output
import database

"""
Public API
"""


def check_system_config():
    """
    Check if the system is configured. If so, load the config.
    """
    if not os.path.exists(config.system_config_file):
        return False
    else:
        system_config = config.load_system_config()

    return system_config


def get_project_config(system_config, args):
    """
    Load the project config object with details about the project
    """
    project_config = config.load_project_config(system_config, args.project)
    if project_config == False:
        print("WARNING: Loading project config failed")
    return project_config

def create_project(args, system_config, project_config):
    """
    Create a new blank project.
    """
    manage_data.create_blank_project(project_config, system_config)


def delete_project(args, system_config, project_config):
    """
    Delete a project
    """
    manage_data.delete_project(project_config)
    # Set the project config to None so update knows it's gone
    project_config = None


def create_indexes(args, system_config, project_config):
    """
    Given a FASTA file of sequences, create one index for each sequence in that FASTA
    """
    create_index.create_all_indexes_for_new_fasta(args.file, project_config,
                                                  system_config.slurm_settings)


def configure_system(args, system_config, project_config):
    """
    Configure the system interactively
    """
    system_config = config.congfigure_system()


def generate_output(args, system_config, project_config):
    """
    Generate the final output file with the annotated guesses about what is/isn't a hopper
    """
    final_output.create_final_output(project_config)


def update_database(args, system_config, project_config):
    """
    Update the alignment database with all jobs that have
    """
    database.update_database(project_config)

def create_index_fasta_from_raw_reads(args, system_config, project_config):
    """
    Create a new index .fasta file from a list of sequences, in a newline-delimited plain text list, passed into -f option
    """
    create_index.find_reads(args.file, args.raw_reads)

def run_alignments(args, system_config, project_config):
    func_to_align_set = {
        "all": align.run_all_samples,
        "taxon_auto": run_all_taxon_auto_samples,
        "true_auto": run_all_true_auto_samples,
        "allo": run_all_allo_samples,
        "same_lane": run_all_intra_lane_samples,
        "other_lane": run_all_cross_lane_samples
    }

    func_to_align_set[args.taxon_set](slurm_settings, project_config, args.file)

def _parse_args(parser: argparse.ArgumentParser):
    """
    Parse args for our program
    """
    parser.add_argument("-a",
                        "--action",
                        help="Action to perform.",
                        type=str,
                        choices=act_to_func.keys(),
                        required=True)
    parser.add_argument(
        "-p",
        "--project",
        help=
        "Project name to act on. Either absolute path or path relative to JULIA2-Projects dir in system config",
        type=str,
        required=True)
    parser.add_argument("-r",
                        "--raw-reads",
                        help="Raw Reads FASTA file, used for creating index fastas",
                        type=str)
    parser.add_argument("-f",
                        "--file",
                        help="File to operate on. Action specific behavior",
                        type=str)
    parser.add_argument(
        "-t",
        "--taxon-set",
        help=
        "For alignment, specify which type of raw read taxons to run alignments against",
        choices=[
            "all", "true_auto", "taxon_auto", "allo", "same_lane", "other_lane"
        ],
        default="all",
        type=str)
    args = parser.parse_args()
    return args


"""
Switchboard
"""
if __name__ == '__main__':
    act_to_func = {
        "create_project": create_project,
        "delete_project": delete_project,
        "configure_project": config.create_project_config,
        "configure_system": configure_system,
        "create_index_fasta_from_raw_reads": create_index_fasta_from_raw_reads,
        "create_indexes": create_indexes,
        "run_alignmnents": run_alignments,
        "update_database": update_database,
        "generate_output": generate_output
    }

    configured = check_system_config()
    if configured == False:
        text = input("WARNING: System not configured. Configure now? Y/n")
        if text != "n":
            configure_system()
        else:
            print("WARNING: System not configured. ")

    system_config = configured

    parser = argparse.ArgumentParser(
        description="View and manage paper portfolios")
    args = _parse_args(parser)

    if args.action != "create_project":
        project_config = get_project_config(system_config, args)
    else:
        project_config = config.ProjectConfig(project_dir=args.project, sample_to_taxon=None, sample_to_taxon_short=None, project_name=os.path.basename(args.project), num_samples=0, num_samples_per_lane=0)
    
    act_to_func[args.action](args, system_config, project_config)

    # Write back any changed configs to the files
    config.update_configs(system_config, project_config)