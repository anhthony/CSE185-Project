for f in mpileupSizes/*.mpileup;
do
    NAME=$(basename "$f" | cut -d'.' -f1)
    SIZE="${NAME#*-}"
    echo "Runtime for mpileup file with $SIZE bases"
    time java -jar VarScan.jar mpileup2snp $f --min-var-freq 0.01 --min-freq-for-hom 0.75 \
    --p-value 0.99 --variants --output-vcf 1 > /dev/null
done