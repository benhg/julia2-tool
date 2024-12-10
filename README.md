# JULIA2

This repository contains the tool `julia2` which may be used to detect index hopping in a dataset.

## Introduction

This repository contains everything needed to go from having the following items:

1. A bunch of raw reads, annotated with whatever we think they came from
2. Specific sequences of interest from within the set of samples
3. A mapping of samples to 

## General Process

The general process that needs to be followed to start answering questions like "Did this sequence come from the sample it's labelled as or a different one?" looks like this:

1. Identify sequences of interest, and the whole dataset
   a. usually a subset ~20%
   b. Use the proteome alignment to suggest which are suspicious
2. Run against raw reads from "true-auto" taxon
3. Inspect 2. for suspicious levels of true-auto alignment (or lack thereof)
   a. apply automatic threshold here
4. Run (1 \ 3) against all data (or at least all of one lane)
5. Generate summary metrics and likely decision from before against results of 4.
   a. Denominator of this metric requires the total number of reads - could use from that lane, but may need to re-set the threshold for that

## Data Objects Created

Each step above produces and requires different data. This is a list:

- First, all the raw reads should be in FASTA files somewhere
- All the sequences of interest should be split out into a FASTA file as well
- We break up the FASTA files of interest into their own single-sequence FASTA
- After that, we create a Bowtie2 index from each one.
- Organize the directories like this:
```
.
├── config.json
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
│   │   ├── index_02_sample_02.out
│   │   └── raw
│   ├── hopping_results.csv
│   └── index_creation
│       └── index_01.out
├── raw_reads
│   ├── sample1.fasta
│   └── sample2.fasta
└── slurm_jobs
    ├── align1.sh
    ├── align2.sh
    ├── index.sh
    └── index1.sh
```

- Intermediate output (from the running of the alignments) will go into `alignment_database_data`
- Summary rows are added to `alignment_database.csv`
- Final answers get put in `hopping_results.csv`
- SLURM job batch files are stored in `slurm_jobs/`

## Steps to Detect Index Hopping

1. Configuration details - store to JSON or use configurator tool or both:
    a. Location of Bowtie2
    b. Location of raw reads (and format - only allow FASTA for now?)
    c. Number of nodes to use for runs
    d. Number of cores per node
    e. Amount of memory per node
    f. Slurm account to use
    g. Working directory. Subdirectories:
        i. output/
        ii. output/alignment_database.csv (currently called summary.csv)
        iii. output/alignment_database_data (the SBatch output files for the alignment runs)
        iv. output/hopping_results.csv (currently called post_summary.csv)
        v. sbatch/ (the SBatch files are created here)
        vi. logs/
        vii. indexes/
        viii. indexes/<sequence_name>/<bowtie index files>
        ix. indexes/<sequence_name>/sequence_name.fasta
        x. raw_reads/
        xi. raw_reads/<sample_name>/<sample_name>.fasta
2. Data related config:
    a. Taxon of each sample
    b. Size of each sample
    c. Other mappings 
2. Create indexes
    a. Take list of sample names
    b. Create directory in indexes/ for that sequence
    c. Create individual FASTA file for that sequence.
    d. Create BOWTIE index in that directory
3. Orchestrate runs
    a. Run all reads from reads of interest against indexes. Options:
        i. All reads
        ii. All reads from taxon-auto
        iii. All reads from true-auto
        iv. All reads from same lane/group
        vi. All reads from opposite lane/group
    b. Put output in output/alignment_database_data
4. Generate database
    a. Run tool currently called "generate_summary.py"
5. Check for hoppers
    a. Run tool currently called "postprocess_summary.py"

## Installation

### Prerequisites:

1. Bowtie2
2. SLURM cluster (local mode to follow)
3. Biopython

### Dev installation:

```
git clone https://github.com/benhg/julia2-tool
cd julia2-tool
python3 setup.py install
```

### `pip` installation:

```
pip install julia2-tool
```

## Usage

The help menu for the project is as follows:

```
usage: tool.py [-h] -a
               {create_project,delete_project,configure_project,configure_system,create_index_fasta_from_raw_reads,create_indexes,run_alignmnents,update_database,generate_output}
               -p PROJECT [-r RAW_READS] [-f FILE] [-t {all,true_auto,taxon_auto,allo,same_lane,other_lane}]

View and manage paper portfolios

optional arguments:
  -h, --help            show this help message and exit
  -a {create_project,delete_project,configure_project,configure_system,create_index_fasta_from_raw_reads,create_indexes,run_alignmnents,update_database,generate_output}, --action {create_project,delete_project,configure_project,configure_system,create_index_fasta_from_raw_reads,create_indexes,run_alignmnents,update_database,generate_output}
                        Action to perform.
  -p PROJECT, --project PROJECT
                        Project name to act on. Either absolute path or path relative to JULIA2-Projects dir in system
                        config
  -r RAW_READS, --raw-reads RAW_READS
                        Raw Reads FASTA file, used for creating index fastas
  -f FILE, --file FILE  File to operate on. Action specific behavior
  -t {all,true_auto,taxon_auto,allo,same_lane,other_lane}, --taxon-set {all,true_auto,taxon_auto,allo,same_lane,other_lane}
                        For alignment, specify which type of raw read taxons to run alignments against
```

Generally, follow the steps in the section "steps to detect index hopping". Example project and system config files can be found in the JSON files in this directory. The system config matches the HPC system at Lewis & Clark College and the project config matches a specific dataset we have been working with there.

0. Install BioPython and Bowtie2
1. Get access to a SLURM cluster (or wait for us to implement local run, and prepare to run your code for a LONG time)
2. Configure the system config to match your SLURM project in `~/.julia2/system_config.json`
3. Create a blank project with the `create_project` option
4. Create the config file for the project you just created
5. Put your raw reads in the `raw_reads` folder in the project you created
6. Prepare your dataset to match what the tool looks for - Raw read files must start with `sXXX_{name}.fasta`.
7. Either use the `create_index_fasta_from_raw_reads` option of this tool, or otherwise create a single FASTA of all the sequences you want to make an index out of
8. Create indexes using the `create_indexes` options
9. Run alignments of raw reads against your indexes. Make sure you are providing a newline-separated list of sequence names in the `--file` option
10. Create a database from the alignments you ran using `update_database`
11. Create an index hopping report using `generate_output`


## Example usage(s)

Run alignment on the raw read that a sequence (that is an index) is labelled as

```
cat ../alignment_test.txt 
s001_c39995_g1_i1_m_10141_SIC_TR_SIC_LRE_LAZ_LRUmerged
 ./julia2.py -f ../alignment_test.txt -a run_alignments -t true_auto -p /home/labs/binford/index_hopping_project/
```
To run other kinds of alignments, just change the `-t` flag.

Create indexes from a FASTA file that contains sequences you wnat to turn into indexes:
```
./julia2.py -p /home/labs/binford/index_hopping_project -a create_indexes -f ~/JULIA-Take-2/src/single_sequence_index/data/*.fasta
```

Clean up the index directory to match FASTAs to the indexes based on them:
```
python3 julia2.py -p /home/labs/binford/index_hopping_project -a cleanup_index_fastas
```

Create a FASTA containing the sequences with names listed in the input file (must be in raw reads)
```
cat ../test_sequences.txt

J00107:16:H2LFJBBXX:1:2106:12845:13236
./julia2.py -f ../test_sequences.txt -a create_index_fasta_from_raw_reads -o ../output_test.fasta

cat ../output_test.fasta 

>s001_J00107:16:H2LFJBBXX:1:2106:12845:13236
CAAATCTCGTGTTATTCTGGCAAAAACTTATCAGTATCAAGCTGAAATTAATGTGATTGAAAAAGACTTCCTCAAATTTCCCCGATAGCATCAGCTTCTGCGATGCATTCTGGATATTTAGGTTAATTTCGCTGTTTAGTATAACGTATTA
```

Generate the alignment database after alignments have been run:
```
./julia2.py -p /home/labs/binford/index_hopping_project/ -a update_database
````

Generate the final output after the database has been generated:
```
./julia2.py -p /home/labs/binford/index_hopping_project/ -a generate_output
```

