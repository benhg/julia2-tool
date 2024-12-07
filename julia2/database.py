"""
Generate a Summary CSV based on the Bowtie alignment output(s)

Summary contains the information that Bowtie produced as well as where to get the full SAM file
"""
from glob import glob
import csv
import subprocess
import sys

sample_to_taxon = {
    "s001": "Drymusa_serrana",
    "s002": "Loxo_arizonica",
    "s003": "Loxo_arizonica",
    "s004": "Loxo_arizonica",
    "s005": "Hexophthalma",
    "s006": "Hexophthalma",
    "s007": "Hexophthalma",
    "s008": "Periegops_MP_VG",
    "s009": "Periegops_MP_VG",
    "s010": "Periegops_MP_VG",
    "s011": "Periegops_MP_WB",
    "s012": "Periegops_VG_H",
    "s013": "Physocyclus",
    "s014": "Plectreurys",
    "s015": "Loxo_reclusa",
    "s016": "Zephryarchea",
    "s017": "Zephryarchea",
    "s018": "Scytodes",
    "s019": "Loxo_rufescens",
    "s020": "Loxo_spinulosa",
    "s021": "Periegops_MP_WB",
    "s022": "Periegops_MP_WB",
    "c74742": "Loxo_spinulosa",
    "c49446": "Scytodes"
}

# Handles variations like WB vs VG
sample_to_taxon_short = {
    "s001": "Drymusa_serrana",
    "s002": "Loxo_arizonica",
    "s003": "Loxo_arizonica",
    "s004": "Loxo_arizonica",
    "s005": "Hexophthalma",
    "s006": "Hexophthalma",
    "s007": "Hexophthalma",
    "s008": "Periegops",
    "s009": "Periegops",
    "s010": "Periegops",
    "s011": "Periegops",
    "s012": "Periegops",
    "s013": "Physocyclus",
    "s014": "Plectreurys",
    "s015": "Loxo_reclusa",
    "s016": "Zephryarchea",
    "s017": "Zephryarchea",
    "s018": "Scytodes",
    "s019": "Loxo_rufescens",
    "s020": "Loxo_spinulosa",
    "s021": "Periegops",
    "s022": "Periegops",
    "c74742": "Loxo_spinulosa",
    "c49446": "Scytodes"
}

# Collected from
# f'bowtie2-inspect --large-index /home/labs/binford/single_sample_indexes/{index_sample}_index/{index_sample}_index
# see get_sizes.py
sample_to_transcript_count = {
    "s001": 146829,
    "s002": 165083,
    "s003": 176396,
    "s004": 177026,
    "s005": 178906,
    "s006": 136577,
    "s007": 157499,
    "s008": 143469,
    "s009": 167627,
    "s010": 183810,
    "s011": 191292,
    "s012": 167942,
    "s013": 223989,
    "s014": 229674,
    "s015": 196020,
    "s016": 161846,
    "s017": 292600,
    "s018": 148310,
    "s019": 191686,
    "s020": 167585,
    "s021": 245515,
    "s022": 170940,
    "c74742": 167585,
    "c49446": 148310
}


def run_cmd(cmd):
    return subprocess.check_output(cmd, shell=True).decode(sys.stdout.encoding)


headers = [
    "reads_sample", "reads_taxon", "index_sample", "index_taxon", "num_reads",
    "num_transcripts", "pairtype", "num_aligned_none", "num_aligned_once",
    "num_aligned_multiple", "none_alignment_rate", "single_alignment_rate",
    "multiple_alignment_rate", "num_aligned_any", "alignment_rate",
    "reads_per_transcript_none", "reads_per_transcript_one",
    "reads_per_transcript_multiple", "reads_per_transcript_any", "exec_time"
]

path = "/home/glick/JULIA-Take-2/src/single_sequence_index/slurm-*.out"
output_file = "/home/labs/binford/taxon_confirmation_indexes/summary.csv"

with open(output_file, "w") as fh:
    writer = csv.writer(fh)
    writer.writerow(headers)

with open(output_file, "a") as fh:
    writer = csv.DictWriter(fh, fieldnames=headers)
    all_files = glob(path)
    for file in all_files:
        with open(file) as fh2:
            try:
                # Name and metadata
                slurm_job_name = file.split("slurm-")[1].split(".out")[0]
                slurm_time_str = run_cmd(
                    f'sacct --format="Elapsed" -j {slurm_job_name}')
                slurm_time = slurm_time_str.split("\n")[-2].strip()

                # Text from stderr (which is summary info)
                data = fh2.readlines()

                #Sample IDs
                index_sample = data[0].split(" ")[0].split("_")[1].strip()
                index_sample_long = data[0].split(" ")[0]
                reads_sample = data[0].split(" ")[1].split("_")[1].strip()

                print(
                    f"Job {slurm_job_name} for index {index_sample} and reads {reads_sample} took {slurm_time}"
                )
                # This is gonna be gross
                row = {
                    "index_sample": index_sample_long,
                    "index_taxon": sample_to_taxon[index_sample],
                    "reads_sample": reads_sample,
                    "reads_taxon": sample_to_taxon[reads_sample],
                    "num_reads": int(data[1].split(" ")[0]),
                    "num_transcripts":
                    sample_to_transcript_count[index_sample],
                    "num_aligned_none": int(data[3].split("(")[0].strip()),
                    "num_aligned_once": int(data[4].split("(")[0].strip()),
                    "num_aligned_multiple": int(data[5].split("(")[0].strip()),
                    "exec_time": slurm_time
                }

                # Pair type

                try:
                    if int(reads_sample.split("s")[1]) <= 11:
                        reads_lane = 1
                    else:
                        reads_lane = 2
                except:
                    reads_lane = 2

                try:
                    if int(index_sample.split("s")[1]) <= 11:
                        index_lane = 1
                    else:
                        index_lane = 2
                except:
                    index_lane = 2

                if reads_sample == index_sample:
                    row["pairtype"] = "True_Auto"
                elif sample_to_taxon_short[reads_sample] == sample_to_taxon_short[
                        index_sample]:
                    row["pairtype"] = "Taxon_Auto"
                elif reads_lane != index_lane:
                    row["pairtype"] = "Other_Lane"
                else:
                    row["pairtype"] = "Allo"

                # Alignment rates
                row["single_alignment_rate"] = row["num_aligned_once"] / int(
                    row["num_reads"])
                row["none_alignment_rate"] = row["num_aligned_none"] / row[
                    "num_reads"]
                row["multiple_alignment_rate"] = row[
                    "num_aligned_multiple"] / row["num_reads"]
                row["num_aligned_any"] = int(row["num_aligned_once"]) + int(
                    row["num_aligned_multiple"])
                row["alignment_rate"] = row["num_aligned_any"] / row[
                    "num_reads"]

                # Reads/transcript scores
                row["reads_per_transcript_none"] = row[
                    "num_aligned_none"] / sample_to_transcript_count[
                        index_sample]
                row["reads_per_transcript_one"] = row[
                    "num_aligned_once"] / sample_to_transcript_count[
                        index_sample]
                row["reads_per_transcript_multiple"] = row[
                    "num_aligned_multiple"] / sample_to_transcript_count[
                        index_sample]
                row["reads_per_transcript_any"] = row[
                    "num_aligned_any"] / sample_to_transcript_count[
                        index_sample]

                writer.writerow(row)
            except Exception as e:

                if "list index out of range":
                    print(f"File {file} is still running")
                    continue


                print(f"failed for file {file}")
