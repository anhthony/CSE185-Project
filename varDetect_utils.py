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

# Function to compile variants that are shared across all outputs 
def shared_vars(fpaths, prefix):
    def read_file(fpaths):
        with open(fpaths, 'r') as f:
            return {line.split('\t', 1)[0]: line.strip() for line in f}

    def write_file(fpaths, data):
        with open(fpaths, 'w') as f:
            for row in data:
                f.write(row + '\n')

    fin = [read_file(fpaths) for fpaths in fpaths]

    sk = set(fin[0].keys())
    for data in fin[1:]:
        sk.intersection_update(data.keys())

    sr = [fin[0][key] for key in sk]

    write_file(prefix+"_shared.VCFss", sr)