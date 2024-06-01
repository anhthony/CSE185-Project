import sys
import os
if __name__ == "__main__":
    file1 = sys.argv[1] #vcf file from VarScan
    file2 = sys.argv[2] #vcf file from VarDetect
    num_f1_vars = 0
    var_pos = []
    with open(file1, "r") as f1:
        vcf_data = f1.readlines()
        for line in vcf_data:
            if not line.startswith("#"):
                data = line.split("\t")
                var_pos.append(data[1])
                num_f1_vars += 1
    f1.close()
    num_f2_vars = 0
    num_match = 0
    mismatch = num_f2_vars - num_match
    with open(file2, "r") as f2:
        vcfss_data = f2.readlines()
        for line in vcfss_data:
            if "CHROM" not in line: #Skip header line
                data = line.split("\t")
                if data[1] in var_pos:
                    num_match += 1
                num_f2_vars += 1
    f2.close()
    match_perc = (num_match/num_f1_vars) *100
    fileName1 = os.path.basename(file1)
    fileName2 = os.path.basename(file2)

    print("%d variants reported in %s" % (num_f1_vars, fileName1))
    print("%d variants reported in %s" % (num_f2_vars, fileName2))
    print("Of %d variants in %s, %d of those variants are also found in %s" % (num_f2_vars, fileName2, num_match, fileName1))
    print("%s contains %f%% of variants found in %s" %(fileName2, match_perc, fileName1))