# imports for util functions
import csv
from collections import defaultdict

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
    mc - minimum read coverage to call variant
    maq - minimum base quality to count a read 
    mr - minimum supporting reads at position to call variant
    mvf - minimum variant allele frequency to call variant
    mffh - minimum variant allele frequency to call homozygous non-ref
    mt - minimum threshold to call variant
"""
def var_call(crchrom, crpos, refb, readc, bpile, crqual, refbf, mc, maq, mr, mvf, mffh, mdt):
    refb = refb.upper()
    mostsig = 0 #Difference in frequencies between observed bases and reference base
    samplef = {"A": 0, "C": 0, "G": 0, "T": 0} #Count of bases in reads
    tba = {"CHROM": crchrom, "POS": crpos, "REF": refb, "GT": "0/1"} # Note that the default here is heterozygous, as we check for homozygous and if the sample allele matches the base, we don't write the variant to output 

    ubpile = bpile.upper()  #Make info about reads all uppercase
    i = 0 #Index of reads string (bpile)
    j = 0 #Index of quality scores string (crqual)
    filt_readc = readc  #Read count after filtering out bad quality reads
    sum_qs = 0 #Sum of quality scores at a position
    
    # Process the aligned reads(https://en.wikipedia.org/wiki/Pileup_format)
    while i < len(ubpile):
        cb = ubpile[i]
        q = crqual[j]
        qs = (ord(q) - 33)
        if filt_readc < mc: #Stop variant calling if filtered read counts is lesser than minimum coverage
            return None
        if cb == "." or cb == ",": 
            if qs >= maq: #Check that quality score for the read pass the minimum quality score allowed
                samplef[refb] += 1 #Increment count for reference base
                sum_qs += qs 
                j += 1
            else: #Decrease read count if read does not pass quality score check
                filt_readc -= 1
        elif cb in samplef: 
            if qs >= maq:
                samplef[cb] += 1 #Increment count of non-reference base
                sum_qs += qs
                j += 1
            else:
                filt_readc -= 1
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
        if (reads_base_freq - mdt) > refbf[key] and reads_base_freq - refbf[key] > mostsig and value >= mr and reads_base_freq >= mvf:
            mostsig = reads_base_freq - refbf[key]
            tba["ALT"] = key 


    # Mean phred is given as qual score 
    tba["QUAL"] = sum_qs/filt_readc

    # Check if homozygous 
    if "ALT" in tba and samplef.get(tba["ALT"], 0) >= mffh:
        tba["GT"] = "1/1"

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
def process_chunk(chunk, refbf, mc, mr, mvf, maq, mffh, mdt):
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
        
        # Perform variant calling for the row and append it to 
        tba = var_call(cr_chrom, cr_pos, cr_ref, cr_reads, cr_bases, cr_qual, refbf, mc, maq, mr, mvf, mffh, mdt)
        if tba != None:
            out_data.append(tba)
    return out_data

# Function to compile variants that are shared across all outputs 
def shared_vars(fpaths, prefix):
    data = []
    header = None
    fcount = len(fpaths)
    
    data_out = defaultdict(lambda: {'count': 0, 'col5': [], 'col6': []})


    for file in fpaths:
        with open(file, 'r', newline='') as f:
            reader = csv.reader(f, delimiter='\t')
            if header is None:
                header = next(reader)  
            else:
                next(reader)  
            for row in reader:
                key = tuple(row[:4])
                data_out[key]['count'] += 1
                data_out[key]['col5'].append(row[4])
                data_out[key]['col6'].append(row[5])
    
    with open(prefix+"_shared.VCFss", 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(header)
        for key, value in data_out.items():
            if value['count'] == fcount:
                c5comb = ';'.join(value['col5'])
                c6comb = ';'.join(value['col6'])
                writer.writerow(list(key) + [c5comb, c6comb])