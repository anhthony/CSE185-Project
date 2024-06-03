FILE=$1
REFGENOME=$2
NAME=$(basename "$FILE" | cut -d'.' -f1)
REGION=$3
CHR="${REGION%:*}"

mkdir -p mpileupSizes
cd ..

for num in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do 
    numbases="$((num*100000))"
    endpos="$((10000000+(numbases)))"

    samtools view -b -T $REFGENOME $FILE chr10:10000001-$endpos > $NAME-$numbases.bam
    cp $NAME.final.cram.crai $NAME-$numbases.final.cram.crai
    rm $NAME.final.cram.crai
    samtools sort $NAME-$numbases.bam > $NAME-$numbases.sorted.bam
    samtools index $NAME-$numbases.sorted.bam
    samtools mpileup -r chr10:10000001-$endpos -f $REFGENOME $NAME-$numbases.sorted.bam > ./analysis/mpileupSizes/$NAME-$numbases.mpileup
    rm $NAME-$numbases.final.cram.crai
    rm $NAME-$numbases.bam
    rm $NAME-$numbases.sorted.bam
    rm $NAME-$numbases.sorted.bam.bai
done