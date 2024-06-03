cd ..

for f in analysis/mpileupSizes/*.mpileup;
do
    NAME=$(basename "$f" | cut -d'.' -f1)
    SIZE="${NAME#*-}"
    echo "Runtime for mpileup file with $SIZE bases"
    time python varDetect.py -m $f -r ./data/hg19.fa -o variants > /dev/null
done
rm variants*