# Main analysis function for calling variants
""" Function to check for variant at a position
    
    Arguments:
    crchrom - chromosome number
    crpos - position in chromosome
    refb - reference base
    readc - number of reads at this position
    bpile - informatiomn about reads at this position
    crqual - quality score that corresponds to each reach
    refbf - frequencies of bases in ref. genome
    mr - minimum supporting reads at position to call variant
    mvf - minimum variant allele frequency to call variant
    mffh - minimum variant allele frequency to call homozygous non-ref
    mt - minimum threshold to call variant
"""
def var_call(crchrom, crpos, refb, readc, bpile, crqual, refbf, mr, mvf, mffh, mt):
    refb = refb.upper()
    mostsig = 1
    samplef = {"A": 0, "C": 0, "G": 0, "T": 0} #Count of bases in reads
    tba = {"CHROM": crchrom, "POS": crpos, "REF": refb, "GT": "0:1"} # Note that the default here is heterozygous, as we check for homozygous and if the sample allele matches the base, we don't write the variant to output 

    ubpile = bpile.upper()  #Make info about reads all uppercase
    i = 0 
    # Process the aligned reads(https://en.wikipedia.org/wiki/Pileup_format)
    while i < len(ubpile):
        cb = ubpile[i]
        if cb == "." or cb == ",": 
            samplef[refb] += 1 #Increment count for reference base
        elif cb in samplef:
            samplef[cb] += 1 #Increment count of non-reference base
        elif cb == "+" or cb == "-": #Found indel, skip over
            i += 1 + int(ubpile[i+1]) + 1
            continue
        elif cb == "^": #Skip start-of-read character and the base quality that follows
            i += 2
            continue
        elif cb == "$": #Skip end-of-read character
            i +=1 
            continue
        i += 1 #Next read
    # Determine and report the most likely alternate allele 
    for key, value in samplef.items():
        reads_base_freq = value / readc
        samplef[key] = reads_base_freq #Frequency of variant in all reads at the position
        #Checking difference between count of each base and mt(?) and making sure there're enough supporting reads and high enough variant allele freq.
        if (reads_base_freq - mt) > refbf[key] and reads_base_freq - refbf[key] > mostsig and value >= mr and reads_base_freq >= mvf:
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
    """ Function to process "chunks" of mpileup file
    
        Arguments:
        chunk - mpileup file
        refbf - frequencies of bases of the reference genome
        mc - minimum coverage at a position to call variant
        mr - minimum SUPPORTING reads at a position to call a variant (i.e must have x reads of variant allele)
        mvf - minimum variant frequency to call variant at a position
        maq - minimum average quality score of reads at a position to call variant
        mffh - minimum frequency to consider a variant non-reference homozygous
        mt - minimum threshold to call a variant
    """
def process_chunk(chunk, refbf, mc, mr, mvf, maq, mffh, mt):
    out_data = []
    for line in chunk: 
        rdata = line.split("\t")
        cr_chrom = rdata[0]
        cr_pos = int(rdata[1])
        cr_ref = rdata[2]
        cr_reads = int(rdata[3])
        cr_bases = rdata[4]
        cr_qual = rdata[5]
            
        # Skip rows with no reads or not enough reads
        if cr_reads == 0 or cr_reads < mc:
            continue
        
        #Skip row if average quality score is lesser than minimum avg. score
        sum_qual = 0
        for score in cr_qual:
            sum_qual += (ord(score)-33)
        if (sum_qual/len(cr_qual) < maq):
            continue
        
        # Perform variant calling for the row and append it to 
        tba = var_call(cr_chrom, cr_pos, cr_ref, cr_reads, cr_bases, cr_qual, refbf, mr, mvf, mffh, mt)
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