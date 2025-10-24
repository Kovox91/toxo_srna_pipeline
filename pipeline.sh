#!/bin/bash

# Function to install conda
install_conda() {
    echo "Conda not found. Installing Miniconda..."
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        CONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        CONDA_URL="https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"
    else
        echo "Unsupported OS. Exiting."
        exit 1
    fi
    
    curl -L $CONDA_URL -o Miniconda3-latest.sh
    bash Miniconda3-latest.sh -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
    source $HOME/miniconda/etc/profile.d/conda.sh
    conda init
    echo "Conda installed successfully."
}

# Function to update conda
update_conda() {
    echo "Updating Conda..."
    conda update -n base -c defaults conda -y
    echo "Conda updated successfully."
}

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    install_conda
else
    echo "Conda is already installed."
    update_conda
fi

# Ensure that the environment is created
ENV_NAME="toxo_sRNA_seq"
ENV_FILE="environment.yaml"

if ! conda info --envs | grep -q $ENV_NAME; then
    echo "Conda environment $ENV_NAME does not exist. Creating it..."
    conda env create -f $ENV_FILE
else
    echo "Conda environment $ENV_NAME already exists"
fi

# Activate the environment using conda
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

# (Re-) Create the necessary reference files
cd references/pseudo_genome
./build_index.sh
cd ../..

# Run the Snakemake pipeline
echo "Running the Snakemake pipeline..."
snakemake --cores all

conda deactivate
