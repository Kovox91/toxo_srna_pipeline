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

# Ask for a descriptive run name
read -p "Enter a short name for this run (e.g. toxo_srna_seq_run1): " RUN_NAME

# Sanitize name (replace spaces or slashes)
RUN_NAME=$(echo "$RUN_NAME" | tr ' /' '__')

# Generate unique topic using name + random suffix
RAND_HEX=$(openssl rand -hex 6)
export NTFY_TOPIC="${RUN_NAME}"

# Optional: use custom server instead of public ntfy.sh
# export NTFY_URL="https://your.self.hosted.ntfy"

# Show user info
echo
echo "ntfy topic generated:  $NTFY_TOPIC"
echo "Subscribe to:          https://ntfy.sh/$NTFY_TOPIC"
echo
echo "You can paste this link into the ntfy app or open it in a browser."
echo

# quick test (you should get a ping on your phone)
curl -H "Title: test" -d "hello from server" "https://ntfy.sh/$NTFY_TOPIC"

# run the pipeline
snakemake -c all

conda deactivate
