#!/usr/bin/env python3
"""
tool.py

The `main` executable for this project.
"""

import argparse
from data_types import *
import os

import config

global act_to_func

global project_config
global system_config

"""
Public API
"""

def check_system_config():
    if not os.path.exists(config.system_config_file):
        return False
    else:
        system_config = config.load_system_config()


    



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
    parser.add_argument("-p",
                        "--project",
                        help="Project name to act on",
                        type=str,
                        required=True)
    parser.add_argument("-r",
                        "--raw-reads",
                        help="Raw Reads FASTA directory",
                        type=str,
                        required=True)

    parser.add_argument(
        "-t",
        "--taxon-set",
        help=
        "For alignment, specify which type of raw read taxons to run alignments against",
        choices=["all", "true_auto", "taxon_auto", "allo", "same_lane", "other_lane"],
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
        "configure_project": configure_project,
        "configure_system": configure_system,
        "display_config": display_config,
        "create_index_fasta_from_raw_reads": create_index_fasta_from_raw_reads,
        "create_indexes": create_indexes,
        "run_alignmnents": run_alignments,
        "update_database": update_database,
        "generate_output": generate_output
    }

    configured = check_system_config()
    if not configured:
        text = input("WARNING: System not configured. Configure now? Y/n")
        if text != "n":
            configure_system()
        else:
            print("WARNING: System not configured. ")


    parser = argparse.ArgumentParser(
        description="View and manage paper portfolios")
    args = _parse_args(parser)
    act_to_func[args.action](args)
