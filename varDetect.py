import multiprocessing
from skbio import DNA
from tqdm import tqdm
import argparse as ap
import time
import os
import sys
from functools import partial

# TODO: tqdm prog bars(one for overall, one for curr file), call for multiple input files(process one by one, don't need to parallel), validation?

# Main analysis function for calling variants
def var_call(crchrom, crpos, refb, readc, bpile, crqual, refbf, mffh, mt):
    refb = refb.upper()
    mostsig = 1
    samplef = {"A": 0, "C": 0, "G": 0, "T": 0}
    tba = {"CHROM": crchrom, "POS": crpos, "REF": refb, "GT": "0:1"} # Note that the default here is heterozygous, as we check for homozygous and if the sample allele matches the base, we don't write the variant to output 

    ubpile = bpile.upper()  
    i = 0 

    # Process the aligned reads(https://en.wikipedia.org/wiki/Pileup_format)
    while i < len(ubpile):
        cb = ubpile[i]
        if cb == "." or cb == ",":
            samplef[refb] += 1
        elif cb in samplef:
            samplef[cb] += 1
        elif cb == "+" or cb == "-":
            i += int(ubpile[i+1])
        elif cb == "^":
            i += 1  
        i += 1

    # Determine and report the most likely alternate allele 
    for key, value in samplef.items():
        samplef[key] = value / readc
        if (value - mt) > refbf[key] and value - refbf[key] > mostsig:
            mostsig = value - refbf[key]
            tba["ALT"] = key 

    # Mean phred is given as qual score 
    tba["QUAL"] = 0
    for q in crqual:
        tba["QUAL"] += (ord(q)-33)
    tba["QUAL"] /= readc

    # Check if homozygous 
    if "ALT" in tba and samplef.get(tba["ALT"], 0) >= mffh:
        tba["GT"] = "1:1"

    return tba

# Function to process each chunk 
def process_chunk(chunk, refbf, mffh, mt):
    out_data = []
    for line in chunk:
        rdata = line.split("\t")
        cr_chrom = rdata[0]
        cr_pos = int(rdata[1])
        cr_ref = rdata[2]
        cr_reads = int(rdata[3])
        cr_bases = rdata[4]
        cr_qual = rdata[5]
            
        # Skip rows with no reads 
        if cr_reads == 0:
            continue
            
        # Perform variant calling for the row and append it to 
        tba = var_call(cr_chrom, cr_pos, cr_ref, cr_reads, cr_bases, cr_qual, refbf, mffh, mt)
        out_data.append(tba)
    return out_data

# Main
if __name__ == '__main__':
    # Set up to parse command line arguments
    parser = ap.ArgumentParser(usage="python varDetect.py --mpileup-file [.mpileup file] --output-vcf [.vcf file] [OPTIONS]", 
                               description="VarDetect: A variant detection tool")

    parser.add_argument('--ref-genome', '-r', required=True, type = str, metavar="REF.fa", help = 'Reference genome')
    parser.add_argument('--mpileup-file', '-m', required=True,type = str, metavar="NAME.mpileup", nargs='+', help = 'The mpileup file(s) to run variant detection on.')
    parser.add_argument('--output-vcfss', '-o', required=True, type = str, metavar="NAME.vcfss", help = 'Output file in VCFss (VCF subset) format; provide the prefix: [prefix]_[filename].vcfss')
    parser.add_argument('--min-coverage', type = int, default = 8, metavar="INT", help = 'Minimum coverage at a position to make a variant call; default: 8')
    parser.add_argument('--min-reads', type = int, default = 2, help = 'Min supporting reads at a position; default: 2')
    parser.add_argument('--min-avg-qual', type = int, default = 15, help = 'Min average base quality at a position; default: 15')
    parser.add_argument('--min-var-freq', type = float, default = 0.01, help = 'Min variant allele freq threshold; default: 0.02')
    parser.add_argument('--min-freq-for-hom', type = float, default = 0.75, help = 'Minimum frequency to call homozygote; default: 0.75')
    parser.add_argument('--min-threshold', type = float, default = 0, help = 'Minimum percent threshold for calling variants; default: 0')
    args = parser.parse_args()
    args_dict = vars(args)

    # Retrieve command args
    pile_up = args_dict["mpileup_file"]
    ref_genome = args_dict["ref_genome"]
    min_cov = args_dict["min_coverage"]
    min_reads = args_dict["min_reads"]
    min_avg_qual = args_dict["min_avg_qual"]
    min_var_freq = args_dict["min_var_freq"]
    min_freq_for_hom = args_dict["min_freq_for_hom"]
    min_thres = args_dict["min_threshold"]
    output_vcfss = args_dict["output_vcfss"]
    
    # Sanity checks for passed arguments
    if min_cov < 0:
        print("Minimum coverage cannot be negative, defaulting to 8")
    if min_reads < 0: 
        print("Minimum reads cannot be negative, defaulting to 2")
    if min_avg_qual < 0:
        print("Minimum average quality cannot be negative, defaulting to 15")
    if min_var_freq < 0:
        print("Minimum variant frequency cannot be negative, defaulting to 0.01")
    if min_freq_for_hom < 0:
        print("Minimum frequency to call homozygote cannot be negative, defaulting to 0.75")
    if min_thres < 0:
        print("Minimum threshold for calling variants cannot be negative, defaulting to 0")

    # Check if all specified files exist, then loop through and run variant calling on each individually 
    pcheck = [file for file in pile_up if not os.path.exists(file)]
    if pcheck:
        print(f"The following file could not be found: {pcheck}")
        sys.exit(1)

    # Get the base frequencies from the genome of interest 
    print("Getting nucleotide frequencies from reference genome...")
    if not os.path.exists(ref_genome):
        print("Reference genome file could not be found")
        sys.exit(1)

    refg = DNA.read(ref_genome, lowercase=True)
    refgbf = refg.frequencies()
    totalb = sum(refgbf.values())
    basef = {base: (count / totalb) for base, count in refgbf.items()}
    
    # Loop through all of the input files 
    for inp in pile_up: 
        # Read mpileup file
        print("Reading in %s..." %(os.path.basename(inp)))
        with open(inp) as f:
            pileup_data = f.readlines()
    
        # Split pileup file into chunks
        print("Chunking pileup...")
        chunk_n = multiprocessing.cpu_count()
        chunk_s = len(pileup_data) // chunk_n
        chunks = [pileup_data[i:i + chunk_s] for i in range(0, len(pileup_data), chunk_s)]
        if len(chunks) == 1:
            print("Data splt into 1 chunk")
        else:
            print("Data split into %d chunks..." %(len(chunks)))

        # Perform variant calling utilizing multiprocessing
        print("Finding variants in each chunk...")
        
        with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
            partialpc = partial(process_chunk, refbf=basef, mffh=min_freq_for_hom, mt=min_thres)
            results = list(tqdm(pool.imap(partialpc, chunks), total=len(chunks)))

        print("All variants found!")

        # Flatten the list of results
        flat_out_data = [item for sublist in results for item in sublist]

        # Print to output file 
        outfname = output_vcfss + "_" + os.path.splitext(os.path.basename(inp))[0] + ".VCFss"
        print("Writing variants to %s..." %(outfname))
        with open(outfname, "w") as f:
            f.write("CHROM\tPOS\tREF\tALT\tQUAL\tGT\n")
            for entry in flat_out_data:
                if entry["REF"] != entry.get("ALT", "") and entry.get("ALT", "") != "":
                    f.write(f"{entry['CHROM']}\t{entry['POS']}\t{entry['REF']}\t{entry.get('ALT', '')}\t{entry['QUAL']}\t{entry['GT']}\n")
