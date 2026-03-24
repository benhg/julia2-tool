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
├── project_config.json // The project configuration file. See example_project_config.json for explanation
├── indexes // The indexes directory
│   ├── seq1_index // Each sequence gets an index subdirectory
│   │   ├── seq1_index.btl2 // Each index has its own Bowtie2 index files
│   │   └── seq1_index.fasta // The sequence's FASTA file is also placed here
│   └── seq2_index
│       ├── seq2_index.btl2
│       └── seq2_index.fasta
├── logs
│   └── julia2.log
├── output
│   ├── alignment_database.csv // This intermediate database contains information about how many sequences from each sample mapped to indexes
│   ├── alignment_database_data // This is where the data for the alignment database goes
│   │   ├── index_01_sample_01.out
│   │   ├── index_01_sample_02.out
│   │   ├── index_02_sample_01.out
│   │   ├── index_02_sample_02.out
│   │   └── raw // This is where the raw SAM files that come out of Bowtie2 go in case they are needed for detailed analysis
│   ├── job_status.csv // One row per submitted job with tracked SLURM state
│   ├── hopping_results.csv // The results from the analysis come here
│   └── index_creation // Index creation output files go here. The indexes go to the indexes/ directory though.
│       └── index_01.out
├── raw_reads // Place your raw reads here
│   ├── sample1.fasta
│   └── sample2.fasta
|── assembled_untranslated_transcripts // Place assembled, untranslated transcripts from Trinity here
|   ├──transcript1.fasta
|   ├──transcript2.fasta
└── slurm_jobs // Slurm submit scripts go here
    ├── align1.sh
    ├── align2.sh
    ├── index.sh
    └── index1.sh
```

- Intermediate output (from the running of the alignments) will go into `alignment_database_data`
- Summary rows are added to `alignment_database.csv`
- Final answers get put in `hopping_results.csv`
- Submitted jobs and their current tracked state go into `job_status.csv`
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
               {create_project,delete_project,configure_project,configure_system,create_index_fasta_from_raw_reads,create_index_fasta_many_reads,create_indexes,run_alignments,update_database,generate_output,cleanup_index_fastas,check_index_creation_err,job_status,job_time_by_node}
               -p PROJECT [-r RAW_READS] [-i READ_ID] [-f FILE] [-o OUTPUT_FILE]
               [-t {all,true_auto,taxon_auto,allo,same_lane,other_lane}]
               [--resume | --reset]

View and manage paper portfolios

optional arguments:
  -h, --help            show this help message and exit
  -a {create_project,delete_project,configure_project,configure_system,create_index_fasta_from_raw_reads,create_index_fasta_many_reads,create_indexes,run_alignments,update_database,generate_output,cleanup_index_fastas,check_index_creation_err,job_status,job_time_by_node}, --action {create_project,delete_project,configure_project,configure_system,create_index_fasta_from_raw_reads,create_index_fasta_many_reads,create_indexes,run_alignments,update_database,generate_output,cleanup_index_fastas,check_index_creation_err,job_status,job_time_by_node}
                        Action to perform.
  -p PROJECT, --project PROJECT
                        Project name to act on. Either absolute path or path relative to JULIA2-Projects dir in system
                        config
  -r RAW_READS, --raw-reads RAW_READS
                        Raw Reads FASTA file, used for creating index fastas
  -i READ_ID, --read-id READ_ID
                        Read ID that goes with the -f option for create_index_fasta_from_raw_reads
  -f FILE, --file FILE  File to operate on. Action specific behavior
  -o OUTPUT_FILE, --output-file OUTPUT_FILE
                        File to output to. Action specific behavior
  -t {all,true_auto,taxon_auto,allo,same_lane,other_lane}, --taxon-set {all,true_auto,taxon_auto,allo,same_lane,other_lane}
                        For alignment, specify which type of raw read taxons to run alignments against
  --resume              Only submit target alignments that are not already completed or active
  --reset               Forget tracked state for the target alignments, then resubmit the full target set
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
   a. Use `-t same_lane` for an intra-lane run
   b. Use `--resume` to continue a partial run without resubmitting completed or active work
   c. Use `--reset` to discard tracked state for the requested target set and resubmit it cleanly
10. Check `job_status` while jobs are in flight instead of reconstructing progress from `squeue` and individual `slurm-*.out` files
11. Create a database from the alignments you ran using `update_database`
12. Create an index hopping report using `generate_output`

## Job Tracking and Resumable Runs

Each submitted index build or alignment job is recorded in `output/job_status.csv`. The tracker stores the SLURM job id, job type, target index, target read sample, submit time, stdout path, and the latest known SLURM state.

Use `job_status` to refresh and print the current tracked state:

```
python3 -m julia2.julia2 -p /home/labs/binford/index_hopping_project -a job_status
```

Use `job_time_by_node` to aggregate tracked job elapsed time by SLURM node:

```
python3 -m julia2.julia2 -p /home/labs/binford/index_hopping_project -a job_time_by_node
```

This report uses `sacct` node allocation data for tracked jobs. When a job spans multiple nodes, its elapsed wall time and allocated CPU count are split evenly across the allocated nodes, and the report includes total time, average time per job, and total CPU-hours for each node.

When enough jobs have already completed, `job_status` also prints an estimated remaining wall time for the active queue based on the recent observed completion rate of completed jobs. If there is not enough completion history yet, it falls back to a median-runtime heuristic. This is still a heuristic and depends on current cluster concurrency and scheduling.

## Interactive System Configuration

`configure_system` now supports interactive setup of `~/.julia2/system_config.json`.

It can:

- reuse values from an existing system config
- prompt for `use_slurm`, partition, account, email, project base directory, and known projects
- prompt for the node CPU map used for round-robin submission
- optionally query `sinfo` and prefill the node CPU map from the active SLURM cluster

Example:

```
python3 -m julia2.julia2 -p /home/labs/binford/index_hopping_project -a configure_system
```

If `sinfo` is available, choose the autodetect prompt to seed the node map. You can still edit the detected values before the file is written.

For large alignment batches, prefer `--resume`:

```
python3 -m julia2.julia2 -p /home/labs/binford/index_hopping_project -a run_alignments -f /home/glick/julia2_tool_inputs/proteome_full_alignment_names.txt -t same_lane --resume
```

`--resume` behavior:

- skip alignments already tracked as `COMPLETED`
- skip alignments already tracked as `RUNNING`, `PENDING`, or `SUBMITTED`
- resubmit alignments tracked as failed or cancelled
- submit alignments that do not yet appear in the ledger

Use `--reset` only when you want to rerun a specific requested target set from scratch:

```
python3 -m julia2.julia2 -p /home/labs/binford/index_hopping_project -a run_alignments -f /home/glick/julia2_tool_inputs/proteome_full_alignment_names.txt -t same_lane --reset
```

`--reset` removes tracked alignment records for the requested workset and then resubmits that full target set. It does not delete SLURM output files or database rows by itself.

## Database Field Meanings

This section explains the fields in the database file

| Field                   | Description                                                                                                                                                                                                                                             |
|-------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| reads_sample            | The sample ID that the reads came from (e.g. "s022")                                                                                                                                                                                                    |
| reads_taxon             | The name of the taxon that the reads came from (e.g. "Drymusa Serrana")                                                                                                                                                                                 |
| index_sample            | The sample ID that the index was generated from - note that determining the accuracy of this label is the point of this project                                                                                                                         |
| index_taxon             | The taxon that the sequence is labelled as coming from - note that determining the accuracy of this label is the point of this project                                                                                                                  |
| num_reads               | The number of raw reads in the raw_reads file                                                                                                                                                                                                           |
| pairtype                | The "type" of pair between the taxon of the reads and the presumptive taxon of the index. The options are "True auto" (they come from the same sample), "Taxon auto" (the same taxon, but a different sample), and "Allo" (Different taxon and sample)  |
| num_aligned_none        | The number of reads that do not align to the index                                                                                                                                                                                                      |
| num_aligned_once        | The number of reads that align exactly once to the index                                                                                                                                                                                                |
| num_aligned_multiple    | The number of reads that align multiple times to the index                                                                                                                                                                                              |
| none_alignment_rate     | The alignment rate of reads that do not align at all (the number that didn't align divided by the total number of reads)                                                                                                                                |
| single_alignment_rate   | The alignment rate of reads that align exactly once (the number that aligned once divided by the total number of reads)                                                                                                                                 |
| multiple_alignment_rate | The alignment rate of reads that align more than once (the number that aligned multiply divided by the total number of reads)                                                                                                                           |
| num_aligned_any         | The sum of reads that aligned once and multiply                                                                                                                                                                                                         |
| alignment_rate          | The num_aligned_any column divided by the num_reads column                                                                                                                                                                                              |
| exec_time               | The time the alignment run took to run to completion                                                                                                                                                                                                    |
|                         |                                                                                                                                                                                                                                                         |

## Final Output Field Meanings

This section explains the fields in the final output file which indicates if hopping likely occurred on this sequence

| Field                                 | Description                                                                                                                                                                                          |
|---------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| sequence_name                         | The name of the sequence of interest (the index)                                                                                                                                                     |
| hopper_status                         | Whether this sequence is likely to be an index hopper or not. Options: LIKELY_HOPPER, MAYBE_HOPPER (there is some evidence but not enough to be confident), and NOT_HOPPER                           |
| source_index_taxon_labelled           | The taxon that this sequence is labelled as                                                                                                                                                          |
| total_reads_mapped                    | The number of reads that mapped to this sequence from all reads                                                                                                                                      |
| num_reads_max_sample                  | The number of mapped reads from the sample that mapped the most times                                                                                                                                |
| read_max_sample                       | The sample ID of the set of raw reads that had the most reads map to this sequence. If the sequence is a hopper, it likely hopped from this sample.                                                  |
| read_max_taxon                        | The taxon name of the set of raw reads that had the most reads map to this sequence. If the sequence is a hopper, it likely hopped from this taxon.                                                  |
| percent_reads_from_max                | The percentage of reads that came from the sample that had the most reads (the read_max_taxon coumn / the total_reads_mapped column)                                                                 |
| num_reads_from_labelled_true_auto     | The number of reads that came from the sample labelled as taxon auto                                                                                                                                 |
| percent_reads_from_labelled_true_auto | The percentage of reads that came from the sample labelled as taxon auto (the num_reads_from_labelled_true_auto column / the total_reads_mapped column)                                              |
| threshold_metric                      | The ratio of the percentage of reads from the sample that mapped the most to the percentage of reads from the labelled true auto. This is used to determine whether the sequence is likely a hopper. |


## Example usage(s)

Run alignment on the raw read that a sequence (that is an index) is labelled as

```
cat ../alignment_test.txt 
s001_c39995_g1_i1_m_10141_SIC_TR_SIC_LRE_LAZ_LRUmerged
 ./julia2.py -f ../alignment_test.txt -a run_alignments -t true_auto -p /home/labs/binford/index_hopping_project/
```
To run other kinds of alignments, just change the `-t` flag.

Resume a partial intra-lane run without duplicating completed or active work:

```
python3 -m julia2.julia2 -p /home/labs/binford/index_hopping_project -a run_alignments -f /home/glick/julia2_tool_inputs/proteome_full_alignment_names.txt -t same_lane --resume
```

Check tracked job status:

```
python3 -m julia2.julia2 -p /home/labs/binford/index_hopping_project -a job_status
```

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
