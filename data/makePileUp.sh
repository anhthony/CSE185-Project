LINK=$1
REFGENOME=$2
NAME=$(basename "$LINK" | cut -d'.' -f1)
REGION=$3
CHR="${REGION%:*}"

samtools view -b -T $REFGENOME $LINK $CHR > $NAME.bam
samtools sort $NAME.bam > $NAME.sorted.bam
samtools index $NAME.sorted.bam
samtools mpileup -r $REGION -f $REFGENOME $NAME.sorted.bam > $NAME.mpileup
