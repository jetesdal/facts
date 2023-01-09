#!/bin/bash

# Take inputs
module_set=$1
module=$2

# Navigate to parent directory of test and run python command runFACTS.py
cd ../..
echo "Create test.sh shell script for ${module_set}.${module}"
python3 runFACTS.py --shellscript test/${module_set}.${module} > test.sh

# Navigate back to moduletest.sh directory
cd -

# Define environment variables
export MODULE_SET="${module_set}"
export MODULE="${module}"

# Run sbatch command
sbatch run_test.sh
