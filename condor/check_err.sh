#!/bin/bash

# Check if the directory is provided as an argument
if [ -z "$1" ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

# Directory containing the .err files
log_dir="$1"

# Output file
output_file="output_err.txt"
echo "Err File, Number of files to process, LooperOutput" > $output_file

# Loop through all .err files in the specified directory
for err_file in "$log_dir"/*.err; do
    # Check if the .err file contains "segmentation violation"
    if grep -q "segmentation violation" "$err_file"; then
        # Extract the base name (without the extension) to find the corresponding .out file
        base_name="${err_file%.err}"
        out_file="${base_name}.out"

        # Check if the corresponding .out file exists
        if [ -f "$out_file" ]; then
            # Extract the required information from the .out file
            number_of_files=$(grep "Number of files to process" "$out_file" | awk '{print $NF}')
            looper_output=$(grep "LooperOutput " "$out_file" | awk '{$1=""; print $0}')

            # Append the information to the output file
            echo "$err_file, $number_of_files, $looper_output" >> $output_file
        fi
    fi
done

