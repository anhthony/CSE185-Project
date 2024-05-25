#Run 'chmod +x ./makePileUp.sh' to make it executable (only need to run once)

FILE=$1
REFGENOME=$2
NAME=$(basename "$FILE" | cut -d'.' -f1)
REGION=$3
CHR="${REGION%:*}"

#Check if it's a .sam file
if [[ "$FILE" == *.sam ]]; then
    samtools view -b $FILE > $NAME.bam
    samtools sort $NAME.bam > $NAME.sorted.bam
    samtools index $NAME.sorted.bam
    mkdir -p mpileup
    samtools mpileup -r $REGION -f $REFGENOME $NAME.sorted.bam > ./mpileup/$NAME.mpileup

#Check that it's a .cram file
elif [[ "$FILE" == *.cram ]]; then
    samtools view -b -T $REFGENOME $FILE $REGION > $NAME.bam
    samtools sort $NAME.bam > $NAME.sorted.bam
    samtools index $NAME.sorted.bam
    mkdir -p mpileup
    samtools mpileup -r $REGION -f $REFGENOME $NAME.sorted.bam > ./mpileup/$NAME.mpileup
    rm $NAME.final.cram.crai    #Comment this out if you want to keep the cram index file
fi
##Comment these out if you want to keep the intermediate .bam and index files

rm $NAME.bam
rm $NAME.sorted.bam
rm $NAME.sorted.bam.bai