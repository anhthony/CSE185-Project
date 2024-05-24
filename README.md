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
python vardetect (?)
```

# Usage Options
```--output-vcf```: 
```--min-coverage```: specify the minimum coverage at a position (default is 8)
```--min-reads```: specify the number of reads at a position (default 3)
```--min-avg-qual 20``` specify the minimum average base quality (default 15)
```--min-var-freq 0.05```: specifiy the minumum variant allele frequency threshold (default 0.02)
```--p-value```: specify the p-value threshold for calling variants
