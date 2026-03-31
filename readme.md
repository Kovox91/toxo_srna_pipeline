# Pipeline README

This README provides instructions for running the sRNA-seq pipeline. The pipeline processes raw FASTQ files, performs quality filtering, adapter trimming, read mapping, and quantification using various bioinformatics tools. Below are the steps required to execute the pipeline successfully.

## Step 1: Prepare Input Data

1. **Raw FASTQ Files:**
   - Place the raw FASTQ files into the `in` folder.
   - Ensure that the files are in gzip-compressed format (`*.fastq.gz`).
   - Name the files according to the following pattern: `<sample_name>_R1_001.fastq.gz` or `<sample_name>_R2_001.fastq.gz` respectively.

2. **Define Samples:**
   - Edit `samples.csv` with one sample and corresponding condition per line.
   - This file specifies the sample names used in the pipeline.

## Step 2: Execute the Pipeline

**Run Pipeline Script:**
   - Execute the `pipeline.sh` script in the terminal
   alternatively:
   - Right-click on the 'pipeline.sh' file and select "Run as Program" from the context menu.

## About `pipeline.sh`:

The `pipeline.sh` script automates the setup and execution of the sRNA-seq pipeline. It performs the following tasks:

1. **Check Dependencies:**
   - Verifies if Conda is installed. If not, installs it.
   - Checks for the required Conda environment. If not present, creates it.

2. **Create Conda Environment:**
   - Sets up a Conda environment.

3. **Execute Snakemake Pipeline:**
   - Runs the Snakemake pipeline within the created Conda environment.
   - Utilizes all available CPU cores for efficient processing.