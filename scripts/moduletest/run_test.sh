#!/bin/bash
#SBATCH --partition=kopp_1
#SBATCH --mem=128000
#SBATCH --time=01:00:00

# Activate the conda environment. Replace env by your environment (e.g. base)
env="facts" ; source ~/.bashrc ;conda activate $env

TESTSCRIPT_DIR=`pwd`

# FACTS main Directory
cd ../

BASEDIR=`pwd`

# Set log file paths
LOG_STDOUT=$BASEDIR/test/${MODULE_SET}.${MODULE}/test.out
LOG_STDERR=$BASEDIR/test/${MODULE_SET}.${MODULE}/test.err

echo "Launching test script.
Logging output to :
$LOG_STDOUT
$LOG_STDERR."

source test.sh > $LOG_STDOUT 2> $LOG_STDERR
