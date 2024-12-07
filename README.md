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
│   │   └── index_02_sample_02.out
│   └── hopping_results.csv
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

## Steps with future tool - tool interface design

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

