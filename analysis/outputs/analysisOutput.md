# Analysis Output Of NA12878_child (chr6:128405804-128605805; 200,002 bases)
```
(py312) (base) Anthony@Anthonys-Macbook-Pro-2 analysis % ./testTimeAndAccuracy.sh /mpileup/NA12878_child.mpileup ./data/hg19.fa 
Running VarScan on /mpileup/NA12878_child.mpileup...
Only SNPs will be reported
Min coverage:   8
Min reads2:     2
Min var freq:   0.01
Min avg qual:   15
P-value thresh: 0.99
Reading input from ./mpileup/NA12878_child.mpileup
200002 bases in pileup file
120 variant positions (95 SNP, 25 indel)
13 were failed by the strand-filter
82 variant positions reported (82 SNP, 0 indel)

real    0m16.936s
user    0m27.988s
sys     0m0.458s
Running VarDetect on /mpileup/NA12878_child.mpileup...
Mpileup file(s):  ['./analysis/mpileup/NA12878_child.mpileup']
Reference genome: ./data/hg19.fa
Minimum coverage: 8
Minimum supporting reads: 2
Minimum average quality score: 15
Minimum variant frequency: 0.01
Minimum frequency for non-reference homozygous: 0.75
Minimum difference threshold: 0.05
Getting nucleotide frequencies from reference genome...
Reading in NA12878_child.mpileup...
Chunking pileup...
Data split into 9 chunks...
Finding variants in each chunk...
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 9/9 [00:03<00:00,  2.94it/s]
All variants found!
Writing variants to variants_NA12878_child.vcfss...
11 variants found in NA12878_child.mpileup.

real    0m47.727s
user    0m53.288s
sys     0m11.601s
82 variants reported in NA12878_child.vcf
11 variants reported in variants_NA12878_child.vcfss
Of 11 variants in variants_NA12878_child.vcfss, 11 of those variants are also found in NA12878_child.vcf
variants_NA12878_child.vcfss contains 13.414634% of variants found in NA12878_child.vcf
```
# Analysis Output Of NA12891_father (chr6:128405804-128605805; 200,002 bases)
```
(py312) (base) Anthony@Anthonys-Macbook-Pro-2 analysis % ./testTimeAndAccuracy.sh /mpileup/NA12891_father.mpileup ./data/hg19.fa
Running VarScan on /mpileup/NA12891_father.mpileup...
Only SNPs will be reported
Min coverage:   8
Min reads2:     2
Min var freq:   0.01
Min avg qual:   15
P-value thresh: 0.99
Reading input from ./mpileup/NA12891_father.mpileup
200002 bases in pileup file
123 variant positions (92 SNP, 31 indel)
6 were failed by the strand-filter
86 variant positions reported (86 SNP, 0 indel)

real    0m15.909s
user    0m29.650s
sys     0m0.386s
Running VarDetect on /mpileup/NA12891_father.mpileup...
Mpileup file(s):  ['./analysis/mpileup/NA12891_father.mpileup']
Reference genome: ./data/hg19.fa
Minimum coverage: 8
Minimum supporting reads: 2
Minimum average quality score: 15
Minimum variant frequency: 0.01
Minimum frequency for non-reference homozygous: 0.75
Minimum difference threshold: 0.05
Getting nucleotide frequencies from reference genome...
Reading in NA12891_father.mpileup...
Chunking pileup...
Data split into 9 chunks...
Finding variants in each chunk...
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 9/9 [00:03<00:00,  2.75it/s]
All variants found!
Writing variants to variants_NA12891_father.vcfss...
16 variants found in NA12891_father.mpileup.

real    0m52.946s
user    0m54.770s
sys     0m13.332s
86 variants reported in NA12891_father.vcf
16 variants reported in variants_NA12891_father.vcfss
Of 16 variants in variants_NA12891_father.vcfss, 16 of those variants are also found in NA12891_father.vcf
variants_NA12891_father.vcfss contains 18.604651% of variants found in NA12891_father.vcf
```

# Analysis Output Of NA12892_mother (chr6:128405804-128605805; 200,002 bases)
```
(py312) (base) Anthony@Anthonys-Macbook-Pro-2 analysis % ./testTimeAndAccuracy.sh /mpileup/NA12892_mother.mpileup ./data/hg19.fa
Running VarScan on /mpileup/NA12892_mother.mpileup...
Only SNPs will be reported
Min coverage:   8
Min reads2:     2
Min var freq:   0.01
Min avg qual:   15
P-value thresh: 0.99
Reading input from ./mpileup/NA12892_mother.mpileup
200002 bases in pileup file
133 variant positions (105 SNP, 28 indel)
17 were failed by the strand-filter
88 variant positions reported (88 SNP, 0 indel)

real    0m18.944s
user    0m30.907s
sys     0m0.705s
Running VarDetect on /mpileup/NA12892_mother.mpileup...
Mpileup file(s):  ['./analysis/mpileup/NA12892_mother.mpileup']
Reference genome: ./data/hg19.fa
Minimum coverage: 8
Minimum supporting reads: 2
Minimum average quality score: 15
Minimum variant frequency: 0.01
Minimum frequency for non-reference homozygous: 0.75
Minimum difference threshold: 0.05
Getting nucleotide frequencies from reference genome...
Reading in NA12892_mother.mpileup...
Chunking pileup...
Data split into 9 chunks...
Finding variants in each chunk...
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 9/9 [00:03<00:00,  2.69it/s]
All variants found!
Writing variants to variants_NA12892_mother.vcfss...
13 variants found in NA12892_mother.mpileup.

real    1m3.771s
user    0m54.637s
sys     0m18.973s
88 variants reported in NA12892_mother.vcf
13 variants reported in variants_NA12892_mother.vcfss
Of 13 variants in variants_NA12892_mother.vcfss, 13 of those variants are also found in NA12892_mother.vcf
variants_NA12892_mother.vcfss contains 14.772727% of variants found in NA12892_mother.vcf
```
# Analysis Output Of NA18555 (chr10:10000000-10200000; 200001 bases)
```
(py312) (base) Anthony@Anthonys-Macbook-Pro-2 analysis % ./testTimeAndAccuracy.sh /mpileup/NA18555.mpileup ./data/hg38.fa       
Running VarScan on /mpileup/NA18555.mpileup...
Only SNPs will be reported
Min coverage:   8
Min reads2:     2
Min var freq:   0.01
Min avg qual:   15
P-value thresh: 0.99
Reading input from ./mpileup/NA18555.mpileup
200001 bases in pileup file
314 variant positions (261 SNP, 53 indel)
0 were failed by the strand-filter
261 variant positions reported (261 SNP, 0 indel)

real    0m17.487s
user    0m22.932s
sys     0m0.400s
Running VarDetect on /mpileup/NA18555.mpileup...
Mpileup file(s):  ['./analysis/mpileup/NA18555.mpileup']
Reference genome: ./data/hg38.fa
Minimum coverage: 8
Minimum supporting reads: 2
Minimum average quality score: 15
Minimum variant frequency: 0.01
Minimum frequency for non-reference homozygous: 0.75
Minimum difference threshold: 0.05
Getting nucleotide frequencies from reference genome...
Reading in NA18555.mpileup...
Chunking pileup...
Data split into 9 chunks...
Finding variants in each chunk...
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 9/9 [00:02<00:00,  3.19it/s]
All variants found!
Writing variants to variants_NA18555.vcfss...
136 variants found in NA18555.mpileup.

real    0m31.307s
user    0m36.718s
sys     0m8.773s
261 variants reported in NA18555.vcf
136 variants reported in variants_NA18555.vcfss
Of 136 variants in variants_NA18555.vcfss, 136 of those variants are also found in NA18555.vcf
variants_NA18555.vcfss contains 52.107280% of variants found in NA18555.vcf
```
# Analysis Output Of HG02603 (chr10:10000000-10200000; 200001 bases)
```
(py312) (base) Anthony@Anthonys-Macbook-Pro-2 analysis % ./testTimeAndAccuracy.sh /mpileup/HG02603.mpileup ./data/hg38.fa
Running VarScan on /mpileup/HG02603.mpileup...
Only SNPs will be reported
Min coverage:   8
Min reads2:     2
Min var freq:   0.01
Min avg qual:   15
P-value thresh: 0.99
Reading input from ./mpileup/HG02603.mpileup
200001 bases in pileup file
342 variant positions (288 SNP, 54 indel)
2 were failed by the strand-filter
286 variant positions reported (286 SNP, 0 indel)

real    0m16.606s
user    0m21.543s
sys     0m0.347s
Running VarDetect on /mpileup/HG02603.mpileup...
Mpileup file(s):  ['./analysis/mpileup/HG02603.mpileup']
Reference genome: ./data/hg38.fa
Minimum coverage: 8
Minimum supporting reads: 2
Minimum average quality score: 15
Minimum variant frequency: 0.01
Minimum frequency for non-reference homozygous: 0.75
Minimum difference threshold: 0.05
Getting nucleotide frequencies from reference genome...
Reading in HG02603.mpileup...
Chunking pileup...
Data split into 9 chunks...
Finding variants in each chunk...
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 9/9 [00:02<00:00,  3.25it/s]
All variants found!
Writing variants to variants_HG02603.vcfss...
119 variants found in HG02603.mpileup.

real    0m30.984s
user    0m36.375s
sys     0m9.270s
286 variants reported in HG02603.vcf
119 variants reported in variants_HG02603.vcfss
Of 119 variants in variants_HG02603.vcfss, 119 of those variants are also found in HG02603.vcf
variants_HG02603.vcfss contains 41.608392% of variants found in HG02603.vcf
```
# Analysis Output Of HG00145 (chr10:10000000-10200000; 200001 bases)
```
(py312) (base) Anthony@Anthonys-Macbook-Pro-2 analysis % ./testTimeAndAccuracy.sh /mpileup/HG00145.mpileup ./data/hg38.fa
Running VarScan on /mpileup/HG00145.mpileup...
Only SNPs will be reported
Min coverage:   8
Min reads2:     2
Min var freq:   0.01
Min avg qual:   15
P-value thresh: 0.99
Reading input from ./mpileup/HG00145.mpileup
200001 bases in pileup file
542 variant positions (459 SNP, 83 indel)
1 were failed by the strand-filter
458 variant positions reported (458 SNP, 0 indel)

real    0m8.438s
user    0m14.599s
sys     0m0.503s
Running VarDetect on /mpileup/HG00145.mpileup...
Mpileup file(s):  ['./analysis/mpileup/HG00145.mpileup']
Reference genome: ./data/hg38.fa
Minimum coverage: 8
Minimum supporting reads: 2
Minimum average quality score: 15
Minimum variant frequency: 0.01
Minimum frequency for non-reference homozygous: 0.75
Minimum difference threshold: 0.05
Getting nucleotide frequencies from reference genome...
Reading in HG00145.mpileup...
Chunking pileup...
Data split into 9 chunks...
Finding variants in each chunk...
100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 9/9 [00:02<00:00,  3.24it/s]
All variants found!
Writing variants to variants_HG00145.vcfss...
235 variants found in HG00145.mpileup.

real    0m34.293s
user    0m37.564s
sys     0m10.535s
458 variants reported in HG00145.vcf
235 variants reported in variants_HG00145.vcfss
Of 235 variants in variants_HG00145.vcfss, 235 of those variants are also found in HG00145.vcf
variants_HG00145.vcfss contains 51.310044% of variants found in HG00145.vcf
```

# Analysis Summary (Runtime (user + sys) + Accuracy):

## NA12878_child
VarScan: 27.988s + 0.458s = 28.446s  
VarDetect: 53.288s + 11.601s = 64.889  
VarDetect is about ~2.28x slower  
VarDetect found 13.41% of the 'true' variants
## NA12891_father
VarScan: 29.650s + 0.386s = 30.036s  
VarDetect: 54.770s + 13.332s = 68.102s  
VarDetect is about ~2.27x slower  
VarDetect found 18.60% of the 'true' variants
## NA12892_mother
VarScan: 30.907s + 0.705s = 31.612s  
VarDetect: 54.637s + 18.973s = 73.610s  
VarDetect is about ~2.33x slower  
VarDetect found 14.77% of the 'true' variants
## NA18555
VarScan: 22.932s + 0.400s = 23.332s  
VarDetect: 36.718s + 8.773s = 45.491s  
VarDetect is about ~1.95x slower  
VarDetect found 52.11% of the 'true' variants
## HG02603
VarScan: 21.543s + 0.347s = 21.890s  
VarDetect: 36.375s + 9.270s = 45.645s  
VarDetect is about ~2.09x slower  
VarDetect found 41.61% of the 'true' variants
## HG00145
VarScan: 14.599s + 0.503s = 15.102s  
VarDetect: 37.564s + 10.535s = 48.009  
VarDetect is about ~3.18x slower  
VarDetect found 51.31% of the 'true' variants