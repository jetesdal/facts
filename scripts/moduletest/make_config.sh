#!/bin/bash

# Initialize variables
nsamps=""
scenario=""
pyear_start=""
pyear_end=""
pyear_step=""
baseyear=""

# Parse command line arguments
while getopts ":h-:n--:s--:p--:e--:t--:b--" opt; do 
  case $opt in
    h) # Display help message and exit
       cat <<EOF
make_config.sh - Create directories and config.yml files for a given list of module_set and module

usage: make_config.sh [-h] [-n NSAMPS] [-s SCENARIO] [-p PYEAR_START]
                      [-e PYEAR_END] [-t PYEAR_STEP] [-b BASEYEAR]
                      modules.txt

optional arguments:
  -h, --help            show this help message and exit
  -n NSAMPS, --nsamps NSAMPS
                        number of samples (default: 500)
  -s SCENARIO, --scenario SCENARIO
                        scenario (default: ssp585)
  -p PYEAR_START, --pyear_start PYEAR_START
                        projection year start (default: 2020)
  -e PYEAR_END, --pyear_end PYEAR_END
                        projection year end (default: 2100)
  -t PYEAR_STEP, --pyear_step PYEAR_STEP
                        projection year step (default: 10)
  -b BASEYEAR, --baseyear BASEYEAR
                        base year (default: 2005)

positional arguments:
  modules.txt           file with module_set and module values
EOF
       exit 0
       ;;
    n) nsamps="$OPTARG"
       ;;
    s) scenario="$OPTARG"
       ;;
    p) pyear_start="$OPTARG"
       ;;
    e) pyear_end="$OPTARG"
       ;;
    t) pyear_step="$OPTARG"
       ;;
    b) baseyear="$OPTARG"
       ;;
    :) # Display error message and exit for missing required argument
       echo "Error: Option -$OPTARG requires an argument." >&2
       exit 1
       ;;
    \?) # Display error message and exit for invalid option
       echo "Error: Invalid option -$OPTARG" >&2
       exit 1
       ;;
  esac
  shift # Shift positional arguments to the left
done

# Shift positional arguments
shift $((OPTIND-1))

# Get the remaining positional arguments
modules_file="$1"

# Check for required modules.txt file
if [ -z "$modules_file" ]; then
  echo "Error: Missing required modules.txt file" >&2
  exit 1
fi

# Set default values if not specified by user
nsamps=${nsamps:-500}
scenario=${scenario:-ssp585}
pyear_start=${pyear_start:-2020}
pyear_end=${pyear_end:-2100}
pyear_step=${pyear_step:-10}
baseyear=${baseyear:-2005}

# Read module_set and module values from modules.txt
while IFS=, read -r module_set module; do
  # Ignore header line
  if [ "$module_set" = "module_set" ]; then
    continue
  fi

  # Create directory and config.yml file
  mkdir -p "../../test/${module_set}.${module}"
  cat > "../../test/${module_set}.${module}/config.yml" <<EOF
global-options:
  nsamps: $nsamps
  scenario: $scenario
  pyear_start: $pyear_start
  pyear_end: $pyear_end
  pyear_step: $pyear_step
  baseyear: $baseyear

sealevel_step:
  onemodule:
    module_set: "$module_set"
    module: "$module"
EOF

  # Add scenario option for older modules (e.g., kopp14)
  case "$module_set" in
    "kopp14" | "ipccar5")
      cat >> "../../test/${module_set}.${module}/config.yml" <<EOF
    options:
      scenario: "rcp85"
EOF
      ;;
  esac

  # Copy first 10 rows of location.lst to directory
  head -n 10 "../../input_files/location.lst" > "../../test/${module_set}.${module}/location.lst"

done < "$modules_file"
