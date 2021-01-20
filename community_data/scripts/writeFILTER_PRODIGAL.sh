#!/bin/bash


BASEDIR=/bio/abchase/MG/elevation_transplant/finaltime

cd $BASEDIR

for f in *_R1.fastq.gz
do

	sampleID=$(echo $f | cut -f1 -d'_')
	fread=$f
	rread=${sampleID}_R2.fastq.gz

	# submit the HPC job as array job to BLASTp in smaller chunks
	echo "#!/bin/bash
#$ -N f_${sampleID}
#$ -m a
#$ -q bio,mic
#$ -pe openmp 4

module load BBMap/37.50
module load bwa/0.7.8
module load picard-tools/1.96
module load samtools/1.3

REFGEN=/bio/abchase/refDB/filtered_euk
REFGENP=\$REFGEN/lolium_perenne_grass.fasta
REFGENF=\$REFGEN/pyrenophora_teres_fungi.fasta
REF=${sampleID}
REFBASE=${BASEDIR}
FFILE=\$REFBASE/${fread}
RFILE=\$REFBASE/${rread}

cd \$REFBASE

bbduk.sh qtrim=rl trimq=10 -Xmx20g threads=4 \\
in1=\$FFILE in2=\$RFILE \\
out1=\$REF.clean1.fq out2=\$REF.clean2.fq

# bwa index -a is \$REFGENP

bwa aln \$REFGENP \$REF.clean1.fq -t 4 > \$REF.clean1.sai
bwa aln \$REFGENP \$REF.clean2.fq -t 4 > \$REF.clean2.sai

bwa samse \$REFGENP \\
\$REF.clean1.sai \$REF.clean1.fq > \$REF.bwaP.R1.sam
bwa samse \$REFGENP \\
\$REF.clean2.sai \$REF.clean2.fq > \$REF.bwaP.R2.sam
rm \$REF.*.sai
rm \$REF.clean1.fq
rm \$REF.clean2.fq

samtools view \$REF.bwaP.R1.sam -b -o \$REF.bwaP.R1.bam
samtools view \$REF.bwaP.R2.sam -b -o \$REF.bwaP.R2.bam
rm \$REF.*.sam

reformat.sh in=\$REF.bwaP.R1.bam out=\$REF.bwaP.R1.fq unmappedonly
reformat.sh in=\$REF.bwaP.R2.bam out=\$REF.bwaP.R2.fq unmappedonly
rm \$REF.*.bam

# bwa index -a is \$REFGENF

bwa aln \$REFGENF \$REF.bwaP.R1.fq -t 4 > \$REF.bwaP.R1.sai
bwa aln \$REFGENF \$REF.bwaP.R2.fq -t 4 > \$REF.bwaP.R2.sai

bwa samse \$REFGENF \\
\$REF.bwaP.R1.sai \$REF.bwaP.R1.fq > \$REF.filter.R1.sam
bwa samse \$REFGENF \\
\$REF.bwaP.R2.sai \$REF.bwaP.R2.fq > \$REF.filter.R2.sam
rm \$REF.*.sai
samtools view \$REF.filter.R1.sam -b -o \$REF.filter.R1.bam
samtools view \$REF.filter.R2.sam -b -o \$REF.filter.R2.bam
rm \$REF.*.sam

reformat.sh in=\$REF.filter.R1.bam out=\$REF.filter.R1.fq unmappedonly
reformat.sh in=\$REF.filter.R2.bam out=\$REF.filter.R2.fq unmappedonly
rm \$REF.*.bam

repair.sh in=\$REF.filter.R1.fq in2=\$REF.filter.R2.fq \\
out=\$REF.filter.clean.R1.fq.gz out2=\$REF.filter.clean.R2.fq.gz

bbmerge.sh \\
in1=\$REF.filter.clean.R1.fq.gz in2=\$REF.filter.clean.R2.fq.gz \\
out=\$REF.filter.clean.merged.fq.gz outu=\$REF.filter.clean.unmerged.fq.gz

reformat.sh in=\$REF.filter.clean.merged.fq.gz out=\$REF.filter.clean.merged.fa
reformat.sh in=\$REF.filter.clean.unmerged.fq.gz out=\$REF.filter.clean.unmerged.fa

cat \$REF.filter.clean.merged.fa \$REF.filter.clean.unmerged.fa > \$REF.filter.total.fa
rm \$REF.filter.clean.merged.fq.gz
rm \$REF.filter.clean.unmerged.fq.gz
rm \$REF.filter.clean.merged.fa
rm \$REF.filter.clean.unmerged.fa

prodigal -i \$REF.filter.total.fa \\
-a \$REF.filter.total.faa -q \\
-f gff -p meta > \$REF.gff

rm -f \$REF.bwa*
rm -f \$REF.gff
rm -f \$REF.filter.R1.fq
rm -f \$REF.filter.R2.fq
rm -f \$REF.filter.total.fa



	" > $sampleID.sh

	qsub $sampleID.sh
	sleep 15

done

