# testTimeAndAccuracy.sh

`testTimeAndAccuracy.sh` is a bash script that takes in an input `.mpileup` file and the reference genome `.fa` file (for the `.mpileup` file) and compare the runtime when ran with `VarScan` vs. when ran with `VarDetect` using the same optional flags. The script also uses the python script `compareVCF.py` to compare the output variants using both tools. Here, we consider the `.vcf` output file from `VarScan` as the "ground truth" and see how accurate the `.vcfss` file from `VarDetect` is compared to the `.vcf` file from `VarScan`. The `.vcf` and `.vcfss` output files from each tool will be stored in the `vcf` and `vcfss` folder, respectively.

## Usage
`./testTimeAndAccuracy.sh [path-to-mpileup] [path-to-ref-genome]`

`[path-to-mpileuop]` - path to the `.mpileup` file. Make sure that the path does not include '.' and is accessible from the `analysis` folder.

`[path-to-ref-genome]` - path to the reference genome. Include `.` in the path and make sure the path is accessible from the `CSE185-Project` folder.

Example usage with provided files in `mpileup`:

`./testTimeAndAccuracy.sh /mpileup/HG02603.mpileup ./data/hg38.fa`

# makePileUpSizes.sh

`makePileUpSizes.sh` is a bash script that generates 20 different `.mpileup` files with varying sizes from the same dataset. The sizes vary by 100,000 bases, from 100,000 to 2,000,000 bases. The script takes in a `ftp` link `.cram` file and a reference genome that is located in the parent directory (due to `cd ..` being used). It will create 20 `.mpileup` files starting at position chr10:10000000, each containing an increment of 100,000 bases. (i.e chr10:10000000-10100000, chr10:10000000-10200000,...). All `.mpileup` files generated are stored in `/mpileupSizes/`.

Usage: `./makePileUpSizes.sh [ftp link for .cram file] [ref-genome]`

Example usage: `./makePileUpSizes.sh ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239361/NA18555.final.cram ./data/hg38.fa`

Note: `/mpileupSizes/` is not found since the folder contains big files.

# runtimeBySizeVarD.sh / runtimeBySizeVarS.sh

`runtimeBySizeVarD.sh` is a bash script that runs the 20 `.mpileup` files generated from `.makePileUpSizes.sh` using `VarDetect` with default parameters to get the runtime for each.

`runtimeBySizeVarS.sh` is a bash script that runs the 20 `.mpileup` files generated from `.makePileUpSizes.sh` using `VarScan` with default parameters to get the runtime for each.

To run, simply execute the script(s), `./runtimeBySizeVarD.sh` or `./runtimeBySizeVarS.sh`.

# analysisOutput.md

`analysisOutput.md` contains output to terminal from running `testTimeAndAccuracy.sh` and further analysis.

# runtimeBySizeVarDetect.txt / runtimeBySizeVarScan.txt

Each file contains output to terminal from `runtimeBySizeVarD.sh` and `runtimeBySizeVarS.sh`.

# graph.ipynb

`graph.ipynb` contains code to graph the results (runtime based on size) from `runtimeBySizeVarD.sh` and `runtimeBySizeVarS.sh`.