import os
import sys
import glob
import json

project_dir = "/home/labs/binford/index_hopping_project/output/alignment_database_data/"
failed = "/home/glick/julia2_tool_inputs/failed_aligns.txt"

failed_details = {}

with open(failed) as fh:
	failed_files = fh.readlines()
	for file in failed_files:
		file = file.split("failed for file ")[1].strip()
		data = open(file).readlines()
		index = data[0].split(" ")[0].split("index_")[1]
		if failed_details.get(index) is None:
			failed_details[index] = 1
		else:
			failed_details[index] += 1

# print(json.dumps(failed_details, indent=2))
print("\n".join(failed_details.keys()))

