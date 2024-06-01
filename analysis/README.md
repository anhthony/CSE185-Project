# testTimeAndAccuracy.sh

`testTimeAndAccuracy.sh` is a bash script that takes in an input `.mpileup` file and the reference genome `.fa` file (for the `.mpileup` file) and compare the runtime when ran with `VarScan` vs. when ran with `VarDetect` using the same optional flags. The script also uses the python script `compareVCF.py` to compare the output variants using both tools. Here, we consider the `.vcf` output file from `VarScan` as the "ground truth" and see how accurate the `.vcfss` file from `VarDetect` is compared to the `.vcf` file from `VarScan`. The `.vcf` and `.vcfss` output files from each tool will be stored in the `vcf` and `vcfss` folder, respectively.

## Usage
`./testTimeAndAccuracy.sh [path-to-mpileup] [path-to-ref-genome]`

`[path-to-mpileuop]` - path to the `.mpileup` file. Make sure that the path does not include '.' and is accessible from the `analysis` folder.

`[path-to-ref-genome]` - path to the reference genome. Include `.` in the path and make sure the path is accessible from the `CSE185-Project` folder.

Example usage with provided files in `mpileup`:

`./testTimeAndAccuracy.sh /mpileup/HG02603.mpileup ./data/hg38.fa`