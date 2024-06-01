#Run as ./testTime.sh [path to mpileup without '.'] [ref-genome]
#NOTE: path to ref-genome should start in the same directory as varDetect.py
#'.mpileup' file should be in the 'analysis' folder

FILE=$1
NAME=$(basename "$FILE" | cut -d'.' -f1)
REFGENOME=$2
echo "Running VarScan on $FILE..."
time java -jar VarScan.jar mpileup2snp .$FILE --min-var-frequency 0.01 --min-freq-for-hom 0.75 \
--p-value 0.99 --variants --output-vcf 1 > $NAME.vcf
mkdir -p vcf
mv $NAME.vcf vcf
mkdir -p vcfss

cd ..
echo "Running VarDetect on $FILE..."
time python varDetect.py -m ./analysis$FILE -r $REFGENOME -o variants
mv variants* ./analysis/vcfss

cd analysis
python compareVCF.py ./vcf/$NAME.vcf ./vcfss/variants_$NAME.vcfss