"""
Functions and stuff to handle configuration-related things
"""

import os
import json
import subprocess
from dataclasses import dataclass
from typing import Dict, List
import logging

logger = logging.getLogger("julia2.config")


@dataclass
class SlurmSettings:
    nodes: Dict[str, int]
    current_node: str
    partition_name: str
    account: str
    email: str = ""


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
        slurm_data = data['slurm_settings']
        if 'email' not in slurm_data:
            slurm_data['email'] = ""
        slurm_settings = SlurmSettings(**slurm_data)
        return SystemConfig(slurm_settings=slurm_settings,
                            use_slurm=data['use_slurm'],
                            project_dir=data.get('project_dir',
                                                 data.get('project_path', '')),
                            projects=data['projects'])

    def to_json(self) -> str:
        """Serializes the SystemConfig object into JSON text."""
        return json.dumps(
            {
                "slurm_settings": {
                    "nodes": self.slurm_settings.nodes,
                    "current_node": self.slurm_settings.current_node,
                    "partition_name": self.slurm_settings.partition_name,
                    "account": self.slurm_settings.account,
                    "email": self.slurm_settings.email
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


def _prompt_with_default(prompt_text, default_value=""):
    if default_value in (None, ""):
        text = input(f"{prompt_text}: ").strip()
    else:
        text = input(f"{prompt_text} [{default_value}]: ").strip()
    if text == "":
        return default_value
    return text


def _prompt_yes_no(prompt_text, default_value):
    default_hint = "Y/n" if default_value else "y/N"
    text = input(f"{prompt_text} [{default_hint}]: ").strip().lower()
    if text == "":
        return default_value
    return text in {"y", "yes"}


def _parse_node_cpu_map(raw_text):
    nodes = {}
    if not raw_text.strip():
        return nodes
    for item in raw_text.split(","):
        item = item.strip()
        if not item:
            continue
        name, cpus = item.split(":", 1)
        nodes[name.strip()] = int(cpus.strip())
    return nodes


def _format_node_cpu_map(nodes):
    return ",".join(f"{name}:{cpus}" for name, cpus in sorted(nodes.items()))


def _detect_slurm_nodes_from_sinfo(partition_name=""):
    cmd = 'sinfo -N -h -o "%P|%n|%c"'
    try:
        lines = subprocess.check_output(cmd, shell=True, text=True).splitlines()
    except Exception as exc:
        logger.warning("Failed to query sinfo: %s", exc)
        return {}

    nodes = {}
    for line in lines:
        parts = line.split("|")
        if len(parts) != 3:
            continue
        partition, node_name, cpu_count = [part.strip() for part in parts]
        partition = partition.replace("*", "")
        if partition_name and partition != partition_name:
            continue
        try:
            nodes[node_name] = int(cpu_count)
        except ValueError:
            continue
    return nodes


def create_system_config():
    return congfigure_system()


def load_project_config(system_config, project):

    if project[0] == "/":
        project_dir = project
    else:
        project_dir = f"{system_config.project_dir}/{project}"
    fh = open(f"{project_dir}/project_config.json", "r")
    json_text = fh.read()
    config = ProjectConfig.from_json(json_text)
    fh.close()
    return config


def update_configs(system_config, project_config):
    if system_config is None or project_config is None:
        return
    project_config_text = project_config.to_json()
    system_config_text = system_config.to_json()
    with open(system_config_file, "w") as fh:
        fh.write(system_config_text)
    with open(f"{project_config.project_dir}/project_config.json", "w") as fh:
        fh.write(project_config_text)


def create_project_config():
    logger.error("Not implemented yet. For now, please write the JSON file yourself")
    raise NotImplemented


def congfigure_system():
    """
	Interactively set the system configuration
	1. Slurm settings
	2. Projects base dir
	"""
    existing = load_system_config() if os.path.exists(system_config_file) else None

    os.makedirs(os.path.dirname(system_config_file), exist_ok=True)

    use_slurm = _prompt_yes_no(
        "Use SLURM for job submission",
        existing.use_slurm if existing else True)

    existing_slurm = existing.slurm_settings if existing else SlurmSettings(
        nodes={},
        current_node="",
        partition_name="",
        account="",
        email="")

    partition_name = _prompt_with_default("SLURM partition name",
                                          existing_slurm.partition_name)

    detected_nodes = {}
    if use_slurm:
        detect_nodes = _prompt_yes_no("Auto-detect node CPU counts from sinfo",
                                      bool(existing_slurm.nodes))
        if detect_nodes:
            detected_nodes = _detect_slurm_nodes_from_sinfo(partition_name)
            if detected_nodes:
                print("Detected nodes:")
                for node_name, cpu_count in sorted(detected_nodes.items()):
                    print(f"  {node_name}: {cpu_count}")
            else:
                print("No nodes detected from sinfo; falling back to manual entry.")

    default_nodes = detected_nodes or existing_slurm.nodes
    node_map_text = _prompt_with_default(
        "Node CPU map as node:cpus pairs separated by commas",
        _format_node_cpu_map(default_nodes))
    node_map = _parse_node_cpu_map(node_map_text)

    current_node_default = existing_slurm.current_node
    if not current_node_default and node_map:
        current_node_default = sorted(node_map.keys())[0]
    current_node = _prompt_with_default("Current node to start round-robin from",
                                        current_node_default)
    account = _prompt_with_default("SLURM account", existing_slurm.account)
    email = _prompt_with_default("SLURM email for notifications",
                                 existing_slurm.email)

    project_dir = _prompt_with_default(
        "Base directory for julia2 projects",
        existing.project_dir if existing else "")
    projects_default = ",".join(existing.projects) if existing else ""
    projects_text = _prompt_with_default(
        "Known project names separated by commas",
        projects_default)
    projects = [project.strip() for project in projects_text.split(",") if project.strip()]

    system_config = SystemConfig(
        slurm_settings=SlurmSettings(nodes=node_map,
                                     current_node=current_node,
                                     partition_name=partition_name,
                                     account=account,
                                     email=email),
        use_slurm=use_slurm,
        project_dir=project_dir,
        projects=projects)

    with open(system_config_file, "w") as fh:
        fh.write(system_config.to_json())

    print(f"Wrote system config to {system_config_file}")
    return system_config
