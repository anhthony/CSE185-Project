# VarDetect
`VarDetect` is a Python package that detects single nucleotide polymorphism (SNPs) between 1 or more datasets against a reference genome. It identifies if there are variants at a positions of a region of a chromosome as specified by the user.

# Prerequisites
Python version 3.9 or newer is required for VarDetect. Installation instructions can be found [here](https://www.python.org/downloads/).

# Installation
Install VarDetect with the following command. 
```
git clone https://github.com/anhthony/CSE185-Project.git
cd CSE185-Project
```

If you don't have `git`, refer [here](https://github.com/git-guides/install-git).

After installing VarDetect, run the following command to install additional required Python packages:

```
pip install -r requirements.txt
```

# Data Generation
`VarDetect` requires two inputs, the `.mpileup` file and the reference genome of the data in the `.mpileup` file. 
## Mpileup File Generation

****THIS STEP IS OPTIONAL****

Provided with VarDetect is a bash script named `makePileUp.sh` to help convert `.sam`/`.cram` files into a `.mpileup` file. When running it for the first time, run `chmod +x makePileUp.sh` to make it an executable. To run the script, `samtools` is required, installation instructions can be found [here](https://www.htslib.org/download/).

To run the script, the usage is `./makePileUp.sh [.sam OR .cram file] [reference genome] [region]`.

Arguments:\
`[.sam OR .cram file]`:  a `.sam`/`.cram` file. You can also put in a ftp link, i.e `ftp://...` as the file.

`[reference genome]`: a `.fa` reference genome file.

`[region]`: region of the chromosome to look at. The format is `chrNUM:start-end` and the format is 1-indexed. For example, `chr4:10500-11000` means the 10500th nucleotide to the 11000th nucleotide on chromosome 4.

## Reference Genome
The reference genome file can be obtained online by referring to metadata about the `.sam`/`.cram`/`.mpileup` file. Ensure that you select the correct reference genome for the data that you are studying. 
# Usage
```
python varDetect.py --ref-genome [.fa ref-genome file] --mpileup-file [.mpileup file] --output-vcfss [.vcfss file] [OPTIONS]
```

# Usage Options
## Required Arguments 
    
`-r`/`--ref-genome` `[.fa file]`: specify a FASTA-formatted reference genome file. More info here [here](https://zhanggroup.org/FASTA/#:~:text=FASTA%20format%20is%20a%20text,by%20lines%20of%20sequence%20data.).   
`-m`/ `--mpileup-file` `[.mpileup file]`: a `.mpileup` file that contains pileups of reads at a single genomic position. More info [here](https://www.htslib.org/doc/samtools-mpileup.html).  
`-o`/`--output-vcfss` `[name_here.vcfss]`: a `.vcfss` file to write all variants to.

## Optional Arguments
There are also many optional arguments to provide additional specifications towards variant calling:  
   
```--min-coverage```: minimum coverage at a position to make a variant call (default: 8)   
```--min-reads```: minimum supporting reads at a position to call variants (default: 2)  
```--min-avg-qual```: minimum base quality ([Phred score](https://www.drive5.com/usearch/manual/quality_score.html)) of a read to count read (default: 15)  
```--min-freq-for-hom```: minimum frequency to call variant a homozygote (default: 0.75)   
```--min-var-freq```: specifiy the minimum variant allele frequency threshold (default: 0.01)   
```--min-threshold```:  minimum percent threshold for calling variants (default: 0)

# File Formats
```ref_genome.fa```   
   
The reference genome is in FASTA format and it contains the following information:
```
>chr[name]
[chromosome sequence]
```

      
```input.mpileup```   
   
The mpileup file is what ```VarDetect``` takes in to run variant detection on. It is a tab delimited text file that represents the alignment of sequence reads to a reference genome. It contains the following data
```
Chromosome    Position    Ref Base    Coverage    Read Bases    Quality Score
```
    
```output.vcfss```   
    
The vcfss output file is a tab delimited text file that stores the gene sequence variants identified from  ```VarDetect```. It contains the following data:
```
Chromosome    Position    Ref Base    Alternative Base    Quality Score
```

# Example Usage 

We have provided 2 sample datasets in`testData`. Dataset `NA18555.mpileup` is pulled from [here](https://www.internationalgenome.org/data-portal/sample/NA18555) (alignment data type under 1000 Genomes 30x on GRCh38) and contains the region chr10:10000000-10200000 (200,000 sequences). Dataset `HG00145.mpileup` is pulled from [here](https://www.internationalgenome.org/data-portal/sample/HG00145) (alignment data type under 1000 Genomes 30x on GRCh38) and contains the region chr10:9000000-10000000 (1,000,000 sequences). Note that hg38.fa, the reference genome fasta file for both datasets, is not provided. You should change the path below to match where you have saved the reference genome locally. Here we use hg38, which can be found on datahub at `/home/your_user/public/genomes/hg38.fa`. We can then perform variant calling in the following manner:
```
python varDetect.py -m /testData/NA18555.mpileup -r /data/hg38.fa -o NA18555.vcfss
```
We expect one output file named `NA18555.vcfss` from running the above command. 

If you want to try generating the `.mpileup` file using the bash script, use the `ftp` link found when selecting "Alignment" and "1000 Genomes 30x on GRCh38" [here](https://www.internationalgenome.org/data-portal/sample/NA18555). Copy the ftp link and you can run it as follow:
`./makePileUp.sh ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240150/HG00145.final.cram /data/hg38.fa chr10:10000000-10200000`

If "permission denied" is thrown when trying to run the bash script, run `chmod +x makePileUp.sh` and try again.

