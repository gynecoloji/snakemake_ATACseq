#!/usr/bin/env python3

import sys

def process_sam(input_sam, output_sam):
    # Open files
    with open(input_sam, 'r') as infile, open(output_sam, 'w') as outfile:
        prev_name = None
        line1 = None

        # Process line by line
        for line in infile:
            line = line.strip()

            # Handle header lines
            if line.startswith('@'):
                outfile.write(line + '\n')
                continue

            # Get read name (first field)
            fields = line.split('\t')
            read_name = fields[0]

            # For non-header lines, pair reads
            if prev_name is None:
                # Store first read of potential pair
                prev_name = read_name
                line1 = line
            else:
                if prev_name == read_name:
                    # Found a pair - write both reads
                    outfile.write(line1 + '\n')
                    outfile.write(line + '\n')
                    prev_name = None
                    line1 = None
                else:
                    # Not a pair - start over with this read
                    prev_name = read_name
                    line1 = line

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input_sam output_sam")
        sys.exit(1)

    input_sam = sys.argv[1]
    output_sam = sys.argv[2]

    process_sam(input_sam, output_sam)
