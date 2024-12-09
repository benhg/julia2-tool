"""
Functions and stuff to handle configuration-related things
"""

import os

import json
from dataclasses import dataclass, field
from typing import Dict, List


@dataclass
class SlurmSettings:
    nodes: Dict[str, int]
    current_node: str
    partition_name: str
    account: str


@dataclass
class SystemConfig:
    slurm_settings: SlurmSettings
    use_slurm: bool
    project_dir: str
    projects: List[str]

    @staticmethod
    def from_json(json_data: str) -> 'SystemConfig':
        """Parses JSON text into a SystemConfig object."""
        data = json.loads(json_data)
        print(json.dumps(data, indent=2))
        slurm_settings = SlurmSettings(**data['slurm_settings'])
        return SystemConfig(slurm_settings=slurm_settings,
                            use_slurm=data['use_slurm'],
                            project_dir=data['project_dir'],
                            projects=data['projects'])

    def to_json(self) -> str:
        """Serializes the SystemConfig object into JSON text."""
        return json.dumps(
            {
                "slurm_settings": {
                    "nodes": self.slurm_settings.nodes,
                    "current_node": self.slurm_settings.current_node,
                    "partition_name": self.slurm_settings.partition_name,
                    "account": self.slurm_settings.account
                },
                "use_slurm": self.use_slurm,
                "project_dir": self.project_dir,
                "projects": self.projects
            },
            indent=4)


@dataclass
class ProjectConfig:
    sample_to_taxon: Dict[str, str]
    sample_to_taxon_short: Dict[str, str]
    project_dir: str
    project_name: str
    num_samples: int
    num_samples_per_lane: int

    @staticmethod
    def from_json(json_data: str) -> 'ProjectConfig':
        """Parses JSON text into a ProjectConfig object."""
        data = json.loads(json_data)
        return ProjectConfig(
            sample_to_taxon=data['sample_to_taxon'],
            sample_to_taxon_short=data['sample_to_taxon_short'],
            project_dir=data['project_dir'],
            project_name=data['project_name'],
            num_samples=int(data['num_samples']),  # Ensure it's an integer
            num_samples_per_lane=int(
                data['num_samples_per_lane'])  # Ensure it's an integer
        )

    def to_json(self) -> str:
        """Serializes the ProjectConfig object into JSON text."""
        return json.dumps(
            {
                "sample_to_taxon": self.sample_to_taxon,
                "sample_to_taxon_short": self.sample_to_taxon_short,
                "project_dir": self.project_dir,
                "project_name": self.project_name,
                "num_samples": self.num_samples,
                "num_samples_per_lane": self.num_samples_per_lane
            },
            indent=4)


# The system configuration file.
# Defaults to ~/.julia2/config.json. Change to an absolute path for shared installs
system_config_file = os.path.expanduser("~/.julia2/system_config.json")


def load_system_config():
    fh = open(system_config_file, "r")
    json_text = fh.read()
    config = SystemConfig.from_json(json_text)
    fh.close()
    return config


def create_system_config():
    print("Not implemented yet. For now, please write the JSON file yourself")
    raise NotImplemented


def load_project_config(system_config, project):

    if project[0] == "/":
        project_dir = project
    else:
        project_dir = f"{system_config.project_dir}/{project}"
    fh = open(f"{project_dir}/project_config.json", "r")
    json_text = fh.read()
    config = SystemConfig.from_json(json_text)
    fh.close()
    return config


def update_configs(system_config, project_config):
    project_config_text = project_config.to_json()
    system_config_text = system_config.to_json()
    with open(system_config_file, "w") as fh:
        fh.write(system_config_text)
    with open(f"{project_config.project_dir}/project_config.json", "w") as fh:
        fh.write(project_config_text)


def create_project_config():
    print("Not implemented yet. For now, please write the JSON file yourself")
    raise NotImplemented


def congfigure_system():
    """
	Interactively set the system configuration
	1. Slurm settings
	2. Projects base dir
	"""
    print("Not implemented yet. For now, please write the JSON file yourself")
    raise NotImplemented
