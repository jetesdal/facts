#!/bin/bash

# check if the required command line arguments were provided
if [ "$#" -ne 2 ]; then
  echo "Error: missing argument(s). Usage: $0 <modules_dir> <output_dir>"
  exit 1
fi

# get the path to the modules directory and the csv file from the command line arguments
modules_dir=$1
output_dir=$2

# create the CSV file
touch "$output_dir"
echo "module_set,module" >> "$output_dir"

# go to the modules directory
cd "$modules_dir"

# for each module set folder
for d in */ ; do
  # skip the facts folder
  if [ "$d" = "facts/" ]; then
    continue
  fi

  # get the module set name
  module_set=${d%/}
  # for each module folder
  for subd in $d/*/ ; do
    # get the module name
    module=${subd%/}
    module=${module##*/}
    # append the module set and module names to the CSV file
    echo "$module_set,$module" >> "$output_dir"
  done
done

# go back to the original directory
cd -
