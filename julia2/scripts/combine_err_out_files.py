import os
import glob

project_dir = "/home/labs/binford/index_hopping_project/output/alignment_database_data/"

files = glob.glob(f"{project_dir}/slurm-*.out")
for file in files:
	err_file = file.replace(".out", ".err")
	# If .out and .err both exist, combine them to be <out file> </n> <errfile>
	# Put this into the .out path
	# Put the old files into a backup dir somewhere else hopefully to never be seen again
	if os.path.isfile(err_file):
		out_file = open(file)
		err_file = open(err_file)
		out_text = out_file.read()
		err_text = err_file.read()
		new_text = f"{out_text}\n{err_text}"

		out_file.close()
		err_file.close()

		with open(file, "w") as fh:
			fh.write(new_text) 

