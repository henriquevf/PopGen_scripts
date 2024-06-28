import argparse
import logging
import os
import subprocess
import sys

#The output file is now an organized tsv

def run_angsd(bam_file, output_prefix, anc_file, ref_file, scaffold, num_threads):
    # Build and run the ANGSD command
    cmd = f"angsd -P {num_threads} -i {bam_file} -anc {anc_file} -dosaf 1 -gl 1 -C 50 -minQ 20 -minmapq 30 -fold 1 -out {output_prefix}.{scaffold} -ref {ref_file} -r {scaffold}"
    subprocess.run(cmd, shell=True, check=True)


def run_realsfs(output_prefix, scaffold):
    # Build and run the realSFS command
    cmd = f"realSFS -nsites 200000 {output_prefix}.{scaffold}.saf.idx > {output_prefix}.{scaffold}.est.ml"
    subprocess.run(cmd, shell=True, check=True)


def calculate_heterozygosity(output_prefix, scaffold):
    # Calculate heterozygosity for each window in the est.ml file
    with open(f"{output_prefix}.{scaffold}.est.ml", "r") as ml_file:
        window_heterozygosities = []
        for line in ml_file:
            ml_data = line.split()
            heterozygous_sites = float(ml_data[1])
            total_sites = sum(float(x) for x in ml_data)
            heterozygosity = heterozygous_sites / total_sites
            window_heterozygosities.append(heterozygosity)

    return window_heterozygosities


def main():
    # Define and parse command-line arguments
    parser = argparse.ArgumentParser(description="Estimate heterozygosity using ANGSD and realSFS.")
    parser.add_argument("bam_file", help="Input BAM file.")
    parser.add_argument("output_prefix", help="Prefix for output files.")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use (default: 1).")
    parser.add_argument("-s", "--scaffold_list", required=True, help="File with list of scaffolds or chromosomes to analyze.")
    args = parser.parse_args()

    # Set up logging
    logging.basicConfig(level=logging.INFO)

    # Load the list of scaffolds from the provided file
    try:
        with open(args.scaffold_list, "r") as f:
            scaffolds = [line.strip() for line in f]
    except FileNotFoundError:
        logging.error("Could not find the scaffold list file. Please provide a file with the list of scaffolds or chromosomes to analyze using the -s option.")
        sys.exit(1)

    # Define reference and ancestral files
    anc_file = "/path/to/anc.fasta"
    ref_file = "/path/to/ref.fasta"

    # Extract the sample name from the BAM file path
    sample_name = os.path.splitext(os.path.basename(args.bam_file))[0]

    # Define the output file for heterozygosity
    heterozygosity_output_file = args.output_prefix + "_heterozygosity.txt"
    with open(heterozygosity_output_file, "w") as out_file:
        # Iterate over scaffolds and estimate heterozygosity for each window
        for scaffold in scaffolds:
            run_angsd(args.bam_file, args.output_prefix, anc_file, ref_file, scaffold, args.threads)
            run_realsfs(args.output_prefix, scaffold)
            window_heterozygosities = calculate_heterozygosity(args.output_prefix, scaffold)

            # Write the results for each scaffold and window to the output file in the desired format
            for i, heterozygosity in enumerate(window_heterozygosities):
                # Remove the specified suffix from the sample name
                sample_name_cleaned = sample_name.replace("_aligned_sorted_dedup_RG", "")
                # Write the data to the file in the desired format
                out_file.write(f"{scaffold}\t{sample_name_cleaned}\t{i + 1}\t{heterozygosity}\n")


if __name__ == "__main__":
    main()
