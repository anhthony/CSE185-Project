# VarDetect
`VarDetect` is a Python package that detects single nucleotide polymorphism (SNPs) between 1 or more datasets against the human hg19 reference genome. It identifies if there are variants at a positions of a region of a chromosome as specified by the user. (?)

# Prerequisites
Required packages can be installed with the included requirements file:
```
pip install -r requirements.txt
```
# Installation
Once all the prerequisitive packages are installed, you can now install VarDetect with the following command: (?)
```
git clone https://github.com/anhthony/CSE185-Project.git
cd CSE185-Project
```
# Usage
```
python varDetect.py --ref-genome [.fa ref genome file] --mpileup-file [.mpileup file] --output-vcf [.vcf file] [OPTIONS]
```

# Usage Options
```VarDetect``` has the following required arguments:   
    
```-r``` ```--ref-genome```: specify the reference genome  
```-m``` ```--mpileup-file```: a text file that represents the alignment of sequence reads to a reference genome    
```-o``` ```--output-vcf```: a text file that stores the gene sequence variants indentified from  ```VarDetect```  
   
There are also many optional arguments to provide additional specifications toward variant calling:  
  
```--min-coverage```: specify the minimum coverage at a position (default is 8)  
```--min-reads```: specify the number of reads at a position (default 3)  
```--min-avg-qual``` specify the minimum average base quality (default 15)  
```--min-var-freq```: specifiy the minumum variant allele frequency threshold (default 0.02)  
```--p-value```: specify the p-value threshold for calling variants  
