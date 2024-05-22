import sys
import os 
import multiprocessing
from skbio import DNA
from tqdm import tqdm

def var_call(crchrom, crpos, refb, readc, bpile, crqual, refbf):
    refb = refb.upper()
    mostsig = 1
    samplef = {"A": 0, "C": 0, "G": 0, "T": 0}
    tba = {"CHROM": crchrom, "POS": crpos, "REF": refb}

    for cb in bpile.upper():
        if cb == "." or cb == ",":
            samplef[refb] += 1
        elif cb in samplef:
            samplef[cb] += 1
    for key, value in samplef.items():
        samplef[key] = value / readc
        if value > refbf[key] and value - refbf[key] > mostsig:
            mostsig = value - refbf[key]
            tba["ALT"] = key 
    tba["QUAL"] = 1
    return tba

def process_chunk(chunk, refbf):
    out_data = []
    for line in chunk:
        parts = line.split("\t")
        cr_chrom = parts[0]
        cr_pos = int(parts[1])
        cr_ref = parts[2]
        cr_reads = int(parts[3])
        cr_bases = parts[4]
        cr_qual = parts[5]
        
        # Skip processing if cr_reads is zero
        if cr_reads == 0:
            continue
        
        tba = var_call(cr_chrom, cr_pos, cr_ref, cr_reads, cr_bases, cr_qual, refbf)
        if tba.get("ALT") is None or tba["ALT"] == cr_ref:
            continue
        out_data.append(tba)
    return out_data


if __name__ == '__main__':
    print("Reading in pileup...")
    pile_up = "/Users/ethanxu/scratch/185_proj/CSE185-Project/data/NA18555.mpileup"
    ref_genome = "/Users/ethanxu/Downloads/hg38.fa"
    with open(pile_up) as f:
        pileup_data = f.readlines()

    refg = DNA.read(ref_genome, lowercase=True)
    refgbf = refg.frequencies()
    totalb = sum(refgbf.values())
    basef = {base: (count / totalb) for base, count in refgbf.items()}
    
    # Split pileup_data into chunks
    print("Chunking pileup...")
    num_chunks = multiprocessing.cpu_count()
    chunk_size = len(pileup_data) // num_chunks
    chunks = [pileup_data[i:i+chunk_size] for i in range(0, len(pileup_data), chunk_size)]

    print("Finding variants...")
    pool = multiprocessing.Pool(processes=num_chunks)
    results = [pool.apply_async(process_chunk, args=(chunk, basef)) for chunk in chunks]
    out_data = [result.get() for result in results]

    flattened_out_data = [item for sublist in out_data for item in sublist]

    print("Writing variants to output file...")
    with open("analysis.vcfss", "w") as f:
        f.write("CHROM\tPOS\tREF\tALT\tQUAL\n")
        for entry in flattened_out_data:
            f.write(f"{entry['CHROM']}\t{entry['POS']}\t{entry['REF']}\t{entry.get('ALT', '')}\t{entry['QUAL']}\n")
