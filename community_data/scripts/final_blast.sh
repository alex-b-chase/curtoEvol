#!/bin/bash
### $ -N submit
### $ -m a
### $ -q bio

module load BBMap/37.50

# take a metagenomic file, split it into smaller chunks, submit array job to run small BLASTp jobs in parallel

BASEDIR=/bio/abchase/MG/elevation_transplant/
MGDIR=$BASEDIR/finaltime

cd $MGDIR

for f in *.filter.total.faa
do
	# split the file into 100 smaller files 
	output=${f%.filter.total.faa}

	partition.sh in=$f out=${output}_%.faa ways=100

	# submit the HPC job as array job to BLASTp in smaller chunks
	echo "#!/bin/bash
#$ -N T${output}.bP
#$ -m a
#$ -q bio,mic,pub8i,pub64
#$ -pe openmp 4
#$ -t 1-100

module load enthought_python/7.3.2
module load blast/2.2.31

MGDIR=$MGDIR
BLASTDB=/bio/abchase/refDB 

cd \${MGDIR}

i=\$(expr \$SGE_TASK_ID - 1)
if [ ! -e \${MGDIR}/${output}_\${i}.blast.txt ]
then
	blastp -query ${output}_\${i}.faa -db \$BLASTDB/total_markergene \\
	-outfmt 6 -max_target_seqs 2 -evalue .00001 -num_threads 4 > ${output}_\${i}.blast.txt

	rm ${output}_\${i}.faa
fi

	" > $output.sh

	qsub $output.sh
done

