"""
Functions and stuff to handle configuration-related things
"""

import os

class SystemConfig:
	def __init__(self):
		pass

class ProjectConfig:
	def __init__(self):
		pass

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
