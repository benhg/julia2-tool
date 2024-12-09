"""
Functions and stuff to handle configuration-related things
"""

import os

class SlurmSettings():
	def __init__(self):
		# Map of node to number of CPUs
		self.nodes = {}
		self.account = ""
		self.partition_name = ""

class SystemConfig:
	def __init__(self):
		# Map of name to path
		self.projects = {}
		self.project_path = ""
		self.slurm_settings = SlurmSettings()
		self.email = ""

class ProjectConfig:
	def __init__(self, project_name, project_path):
		self.project_name = project_name
		self.project_path = user_path if user_path[0] == "/" else f"{system_config.project_path}/{self.project_name}"

# The system configuration file.
# Defaults to ~/.julia2/config.json. Change to an absolute path for shared installs
system_config_file = os.path.expanduser("~/.julia2/config.json")


def load_system_config():
	pass

def create_system_config():
	pass

def load_project_config():
	pass

def create_project_config():
	pass


def congfigure_system():
	"""
	Interactively set the system configuration
	1. Slurm settings
	2. Projects base dir
	"""
	pass
