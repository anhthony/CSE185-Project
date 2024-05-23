import multiprocessing
from skbio import DNA
from tqdm import tqdm

# TODO: tqdm prog bars(one for overall, one for curr file), call for multiple input files(process one by one, don't need to parallel), validation?

# Main analysis function for calling variants
def var_call(crchrom, crpos, refb, readc, bpile, crqual, refbf):
    refb = refb.upper()
    mostsig = 1
    samplef = {"A": 0, "C": 0, "G": 0, "T": 0}
    tba = {"CHROM": crchrom, "POS": crpos, "REF": refb}

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
        if value > refbf[key] and value - refbf[key] > mostsig:
            mostsig = value - refbf[key]
            tba["ALT"] = key 

    # Mean phred is given as qual score 
    tba["QUAL"] = 0
    for q in crqual:
        tba["QUAL"] += (ord(q)-33)
    tba["QUAL"] /= readc

    return tba

# Function to process each chunk 
def process_chunk(chunk, refbf):
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
        tba = var_call(cr_chrom, cr_pos, cr_ref, cr_reads, cr_bases, cr_qual, refbf)
        out_data.append(tba)
    return out_data

# Main
if __name__ == '__main__':
    # Paths are hardcoded!
    print("Reading in pileup...")
    pile_up = "/Users/ethanxu/scratch/185_proj/CSE185-Project/data/NA18555.mpileup"
    ref_genome = "/Users/ethanxu/Downloads/hg38.fa"
    with open(pile_up) as f:
        pileup_data = f.readlines()

    # Get the base frequencies from the genome of interest 
    refg = DNA.read(ref_genome, lowercase=True)
    refgbf = refg.frequencies()
    totalb = sum(refgbf.values())
    basef = {base: (count / totalb) for base, count in refgbf.items()}
    
    # Split pileup file into chunks
    print("Chunking pileup...")
    chunk_n = multiprocessing.cpu_count()
    chunk_s = len(pileup_data) // chunk_n
    chunks = [pileup_data[i:i + chunk_s] for i in range(0, len(pileup_data), chunk_s)]

    # Perform variant calling utilizing multiprocessing
    print("Finding variants...")
    pool = multiprocessing.Pool(processes=chunk_n)
    results = []
    for chunk in chunks:
        result = pool.apply_async(process_chunk, args=(chunk, basef))
        chunk_d = result.get()
        results.append(chunk_d)

    # Flatten the list of results
    flat_out_data = [item for sublist in results for item in sublist]

    # Print to output file 
    print("Writing variants to output file...")
    with open("analysis.vcfss", "w") as f:
        f.write("CHROM\tPOS\tREF\tALT\tQUAL\n")
        for entry in flat_out_data:
            if entry["REF"] != entry.get("ALT", "") and entry.get("ALT", "") != "":
                f.write(f"{entry['CHROM']}\t{entry['POS']}\t{entry['REF']}\t{entry.get('ALT', '')}\t{entry['QUAL']}\n")
