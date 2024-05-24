#Run 'chmod +x ./makePileUp.sh' to make it executable (only need to run once)

LINK=$1
REFGENOME=$2
NAME=$(basename "$LINK" | cut -d'.' -f1)
REGION=$3
CHR="${REGION%:*}"

samtools view -b -T $REFGENOME $LINK $CHR > $NAME.bam
samtools sort $NAME.bam > $NAME.sorted.bam
samtools index $NAME.sorted.bam
mkdir -p mpileup
samtools mpileup -r $REGION -f $REFGENOME $NAME.sorted.bam > ./mpileup/$NAME.mpileup

##Comment these out if you want to keep the intermediate .bam files
rm $NAME.final.cram.crai
rm $NAME.bam
rm $NAME.sorted.bam
rm $NAME.sorted.bam.bai
