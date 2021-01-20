#!/bin/bash
#$ -N MMLR14002
#$ -m a
#$ -q mic
#$ -pe openmp 16

echo "Job started on `hostname` at `date`"

module load samtools/1.3
module load BBMap/37.50
module load bwa/0.7.8 
module load pilon

REF=MMLR14002
READBASE=/bio/abchase/genomes/LRBACE/rawdata
FFILE=$READBASE/UCI_2_1.fastq
RFILE=$READBASE/UCI_2_2.fastq

BBMAPDIR=/data/users/abchase

OUTBASE=/bio/abchase/genomes/pacbio/
OUTDIR=$OUTBASE/spades/$REF
THREAD=16

bbduk.sh in1=$FFILE in2=$RFILE \
out1=$OUTBASE/$REF.R1.clean.fq.gz out2=$OUTBASE/$REF.R2.clean.fq.gz \
minlen=25 qtrim=rl trimq=10 ktrim=r k=25 ref=$BBMAPDIR/nextera.fa.gz hdist=1

# interleave for the PACBIO correction
reformat.sh \
in1=$OUTBASE/$REF.R1.clean.fq.gz \
in2=$OUTBASE/$REF.R2.clean.fq.gz \
out=$OUTBASE/$REF.R1R2.clean.fq.gz

repair.sh \
in=$OUTBASE/$REF.R1R2.clean.fq.gz \
out=$OUTBASE/$REF.R1R2.clean.fix.fq \
outsingle=$OUTBASE/$REF.R1R2.clean.single.fq.gz

# BWA spits out error from reads not being names the same so fix the suffix of read header
sed -e '/^[@]/ s/\// /' $OUTBASE/$REF.R1R2.clean.fix.fq > $OUTBASE/$REF.R1R2.clean.fix2.fq


#############################################################################
####### now run the pacbio correction
#############################################################################

PACBIO_ASSEMBLY=$OUTBASE/MMLR14002_pacbio.fasta

cd $OUTBASE

bwa index $PACBIO_ASSEMBLY

bwa mem -t $THREAD -p $PACBIO_ASSEMBLY $OUTBASE/$REF.R1R2.clean.fix2.fq | \
samtools sort -T tmp -o $OUTBASE/$REF.sorted.bam -

samtools index $OUTBASE/$REF.sorted.bam

# use pilon which will generate a new consensus sequence for the assembly
java -jar /data/apps/pilon/1.22/bin/pilon-1.22.jar \
--genome $PACBIO_ASSEMBLY \
--frags $OUTBASE/$REF.sorted.bam \
--output $OUTBASE/corrected


echo "Job Ended at `date`"









