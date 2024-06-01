import multiprocessing
from skbio import DNA
from tqdm import tqdm
import argparse as ap
import os
import sys
from functools import partial
from varDetect_utils import *

# Main
if __name__ == '__main__':
    # Set up to parse command line arguments
    parser = ap.ArgumentParser(usage="python varDetect.py --mpileup-file [.mpileup file] --output-vcf [.vcf file] [OPTIONS]", 
                               description="VarDetect: A variant detection tool")

    parser.add_argument('--ref-genome', '-r', required=True, type = str, metavar="REF.fa", help = 'Reference genome')
    parser.add_argument('--mpileup-file', '-m', required=True,type = str, metavar="NAME.mpileup", nargs='+', help = 'The mpileup file(s) to run variant detection on.')
    parser.add_argument('--output-vcfss', '-o', required=True, type = str, metavar="NAME.vcfss", help = 'Output file in VCFss (VCF subset) format; provide the prefix: [prefix]_[filename].vcfss')
    parser.add_argument('--min-coverage', type = int, default = 8, metavar="INT", help = 'Minimum coverage at a position to make a variant call; default: 8')
    parser.add_argument('--min-reads', type = int, default = 2, help = 'Min supporting reads at a position to call variant; default: 2')
    parser.add_argument('--min-avg-qual', type = int, default = 15, help = 'Min average base quality at a position; default: 15')
    parser.add_argument('--min-var-freq', type = float, default = 0.01, help = 'Min variant allele freq threshold; default: 0.01')
    parser.add_argument('--min-freq-for-hom', type = float, default = 0.75, help = 'Minimum frequency to call homozygote; default: 0.75')
    parser.add_argument('--min-diff-threshold', type = float, default = 0.05, help = 'Minimum difference between variant allele and reference allele frequencies for calling variants; default: 0.05')
    args = parser.parse_args()
    args_dict = vars(args)

    # Retrieve command args
    pile_ups = args_dict["mpileup_file"]
    ref_genome = args_dict["ref_genome"]
    min_cov = args_dict["min_coverage"]
    min_reads = args_dict["min_reads"]
    min_avg_qual = args_dict["min_avg_qual"]
    min_var_freq = args_dict["min_var_freq"]
    min_freq_for_hom = args_dict["min_freq_for_hom"]
    min_diff_thres = args_dict["min_diff_threshold"]
    output_vcfss = args_dict["output_vcfss"]
    print("Mpileup file(s): ", pile_ups)
    print("Reference genome: %s\nMinimum coverage: %s\nMinimum supporting reads: %s\nMinimum average quality score: %s\nMinimum variant frequency: %s\nMinimum frequency for non-reference homozygous: %s\nMinimum difference threshold: %s" 
          % (ref_genome,min_cov,min_reads, min_avg_qual, min_var_freq, min_freq_for_hom, min_diff_thres))
    # Sanity checks for passed arguments, if does not pass, set to default
    if min_cov < 0:
        print("Minimum coverage cannot be negative, defaulting to 8")
        min_cov = 8
    if min_reads < 0: 
        print("Minimum reads cannot be negative, defaulting to 2")
        min_reads = 2
    if min_avg_qual < 0:
        print("Minimum average quality cannot be negative, defaulting to 15")
        min_avg_qual = 15
    if min_var_freq < 0:
        print("Minimum variant frequency cannot be negative, defaulting to 0.01")
        min_var_freq = 0.01
    if min_freq_for_hom < 0:
        print("Minimum frequency to call homozygote cannot be negative, defaulting to 0.75")
        min_freq_for_hom = 0.75
    if min_diff_thres < 0:
        print("Minimum threshold for calling variants cannot be negative, defaulting to 0")
        min_diff_thres = 0.05

    # Check if all specified files exist, then loop through and run variant calling on each individually 
    pcheck = [file for file in pile_ups if not os.path.exists(file)]
    if pcheck:
        print(f"The following file could not be found: {pcheck}")
        sys.exit(1)

    # Get the base frequencies from the genome of interest 
    print("Getting nucleotide frequencies from reference genome...")
    if not os.path.exists(ref_genome):
        print("Reference genome file '%s' could not be found" % (ref_genome))
        sys.exit(1)

    refg = DNA.read(ref_genome, lowercase=True)
    refgbf = refg.frequencies()
    totalb = sum(refgbf.values())
    basef = {base: (count / totalb) for base, count in refgbf.items()}
    
    outfs = []
    # Loop through all of the input files 
    for inp in pile_ups: 
        # Read mpileup file
        print("Reading in %s..." %(os.path.basename(inp)))
        with open(inp) as f:
            pileup_data = f.readlines()
    
        # Split pileup file into chunks
        print("Chunking pileup...")
        num_chunks = multiprocessing.cpu_count()
        chunk_size = len(pileup_data) // num_chunks
        if chunk_size == 0:
            chunks = [pileup_data[i:i + 1] for i in range(0, len(pileup_data), 1)]
        else:
            chunks = [pileup_data[i:i + chunk_size] for i in range(0, len(pileup_data), chunk_size)]
        if len(chunks) == 1:
            print("Data split into 1 chunk")
        else:
            print("Data split into %d chunks..." %(len(chunks)))

        # Perform variant calling utilizing multiprocessing
        print("Finding variants in each chunk...")
        
        with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
            partialpc = partial(process_chunk, refbf=basef, mc = min_cov, mr = min_reads, mvf = min_var_freq,
                                maq = min_avg_qual, mffh=min_freq_for_hom, mdt=min_diff_thres)
            results = list(tqdm(pool.imap(partialpc, chunks), total=len(chunks)))

        print("All variants found!")

        # Flatten the list of results
        flat_out_data = [item for sublist in results for item in sublist]

        # Print to output file 
        outfname = output_vcfss + "_" + os.path.splitext(os.path.basename(inp))[0] + ".VCFss"
        print("Writing variants to %s..." %(outfname))
        with open(outfname, "w") as f:
            f.write("CHROM\tPOS\tREF\tALT\tQUAL\tGT\n")
            num_vars = 0
            for entry in flat_out_data:
                if entry["REF"] != entry.get("ALT", "") and entry.get("ALT", "") != "":
                    f.write(f"{entry['CHROM']}\t{entry['POS']}\t{entry['REF']}\t{entry.get('ALT', '')}\t{entry['QUAL']}\t{entry['GT']}\n")
                    num_vars += 1
        print("%d variants found in %s.mpileup." % (num_vars, os.path.splitext(os.path.basename(inp))[0]))
        outfs.append(outfname)
    
    # Find common variants across output files if applicable 
    if len(outfs) > 1:
        print("Looking for shared variants...")
        shared_vars(outfs, output_vcfss)
        print("Writing shared variants to %s..." %(output_vcfss+"_shared.VCFss"))
