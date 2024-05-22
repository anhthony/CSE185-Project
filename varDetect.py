import pandas as pd 
import sys
import os 
from skbio import DNA
from tqdm import tqdm

# refa: reference allele 
# bpile: bases from pileup file 
# note: ignoring all insertions and deletions and directionality
def var_call(crchrom, crpos, refb, readc, bpile, crqual, refbf):
    refb = refb.upper()
    mostsig = 1
    samplef = {"A": 0, "C": 0, "G": 0, "T": 0}
    tba = {}
    tba["CHROM"] = crchrom
    tba["POS"] = crpos
    tba["REF"] = refb

    for i in range(len(bpile)):
        cb = bpile[i]
        cb = cb.upper()
        if cb == "." or cb == ",":
            samplef[refb] += 1
        elif cb in samplef:
            samplef[cb] += 1
        elif cb == "+" or cb == "-":
            i += int(bpile[i+1])
        elif cb == "^":
            i += 1

    for key in samplef:
        samplef[key] /= readc
        if samplef[key] > refbf[key] and (samplef[key] - refbf[key] > mostsig):
            mostsig = samplef[key] - refbf[key]
            tba["ALT"] = key 

    # placeholder for qual score 
    tba["QUAL"] = 1

    return tba


# TODO: remove redundant imports from bfCalc, update requirements.txt
# TODO: if multithreading flag is passed, find base frequencies and call vars in parallel
# TODO: maybe seperate bf calc functions into a seperate file for readability
# Main, handle passed arguments and make calls to analysis functions
if __name__ =='__main__':
    
    # TODO: argument handling, paths are hardcoded right now  
    print("Reading in pileup...")
    pile_up = "/Users/ethanxu/scratch/185_proj/CSE185-Project/data/NA18555.mpileup"
    ref_genome = "/Users/ethanxu/Downloads/hg38.fa"
    columns = ['Chromosome', 'Position', 'Reference_Base', 'Coverage', 'Bases', 'Qualities']
    df = pd.read_csv(pile_up, sep='\t', header=None, names=columns)

    refg = DNA.read(ref_genome, lowercase=True)
    refgbf = refg.frequencies()
    totalb = sum(refgbf.values())
    basef = {base: (count / totalb) for base, count in refgbf.items()}
    outdf = pd.DataFrame(columns=["CHROM", "POS", "REF", "ALT", "QUAL"])

    # Parse the pileup row by row 
    print("Finding variants...")
    for index, row in tqdm(df.iterrows(), total=len(df), desc="Parsing pileup"):
        cr_chrom = row["Chromosome"]
        cr_pos = row["Position"]
        cr_ref = row["Reference_Base"]
        cr_reads = row["Coverage"]
        cr_bases = row["Bases"]
        cr_qual = row["Qualities"]
        tba = var_call(cr_chrom, cr_pos, cr_ref, cr_reads, cr_bases, cr_qual, basef)
        outdf = pd.concat([outdf, pd.Series(tba)], ignore_index=True)

    # TODO: flag for name of output, edit print and out file name
    print("Writing variants to output file...")
    outdf.to_csv("analysis.vcfss", sep='\t', index=False)
