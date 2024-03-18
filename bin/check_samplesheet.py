#!/usr/bin/env python3

import csv
import sys


def process_file(input_file, output_file):
    with open(input_file, "r", newline="") as infile:
        with open(output_file, "w", newline="") as outfile:
            reader = csv.reader(infile, delimiter="\t")
            writer = csv.writer(outfile, delimiter="\t")

            header = next(reader)
            writer.writerow(header)  # Write header to output file

            for row in reader:
                if row and not row[0].startswith("#"):  # Skip lines starting with #
                    # Check if necessary columns are not empty
                    if all(
                        row[header.index(col)]
                        for col in [
                            "flowcell",
                            "library",
                            "sample_type",
                            "genome",
                            "read1",
                            "read2",
                        ]
                    ):
                        writer.writerow(row)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python check_samplesheet.py input_file output_file")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    process_file(input_file, output_file)
