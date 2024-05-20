FQ1=$1
FQ2=$2

bwa mem ./hg19.fa ./$FQ1 ./$FQ2 > data.sam
samtools view -S -b data.sam > data.bam
samtools sort data.bam > data.sorted.bam
samtools index data.sorted.bam
samtools mpileup -r chr6:128405804-128605805 -f ./hg19.fa data.sorted.bam > data.mpileup

echo "bwa mem ./hg19.fa ./$FQ1 ./$FQ2"
