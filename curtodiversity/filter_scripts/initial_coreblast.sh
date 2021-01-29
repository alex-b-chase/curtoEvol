#!/bin/bash

module load BBMap/37.50

# take a metagenomic file, split it into smaller chunks, submit array job to run small BLASTp jobs in parallel

BASEDIR=/bio/abchase/MG/elevation_transplant
MGDIR=$BASEDIR/initialtime

cd $MGDIR

for f in *.filter.total.faa
do
	# split the file into 100 smaller files 
	output=${f%.filter.total.faa}

	partition.sh in=$f out=${output}_%.faa ways=100

	# submit the HPC job as array job to BLASTp in smaller chunks
	echo "#!/bin/bash
#$ -N ${output}.bP
#$ -q pub8i,bio,mic,free*
#$ -pe openmp 4
#$ -ckpt restart
#$ -t 1-100
#$ -j y
# RESTART_EMAIL_SUMMARY

module load enthought_python/7.3.2
module load blast/2.2.31

MGDIR=$MGDIR
BLASTDB=/bio/abchase/refDB/coregenes
OUTDIR=$BASEDIR/curto_coregenes

cd \${MGDIR}

i=\$(expr \$SGE_TASK_ID - 1)

if [ -e \$MGDIR/${output}_\${i}.faa ]
then

	rm -f \${OUTDIR}/${output}_\${i}.blast.txt

	blastp -query ${output}_\${i}.faa -db \$BLASTDB/total_coregenes \\
	-outfmt 6 -max_target_seqs 2 -evalue 1e-20 -num_threads 4 > \${OUTDIR}/${output}_\${i}.blast.txt

	sleep 10

	rm -f ${output}_\${i}.faa

	echo \"done with \${i}\"

else 
	echo \"already done!\"
fi

	" > $output.sh

	qsub $output.sh
done

