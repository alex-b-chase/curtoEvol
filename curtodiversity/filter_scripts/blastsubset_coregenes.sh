#!/bin/bash

# take the output from BLAST results and subset into marker gene specific files
# for some reason, (think the memory alottment), the HPC runs this SUPER slowly
# use locally and then do the HMMER and PPLACER analyses on the HPC

# first, subset from BIG MG files to make rest easier each time reading in

BASEDIR=/Volumes/Bennett_BACKUP/Research/curtobacterium/curto-evolution/transplant_communitybags/metagenomic_analyses
MGDIR=$BASEDIR/prodigal_reads

CURTODIR=$BASEDIR/curto_coregenes
OUTDIR=$CURTODIR/blast_output

cd $OUTDIR

# subset the BLAST results from the MG libraries

# for f in *.blast.total.txt
# do

# 	metagenome=${f%.blast.total.txt}

# 	if [[ $f == Z* ]]
# 		then
# 			READDIR=$MGDIR/time0
# 		else 
# 			READDIR=$MGDIR/time1
# 	fi

# 	cut -f1 $f | sort -u | sed 's/[ ]*$//' > ${metagenome}.temp.txt

# 	echo "Processing ${metagenome}..."


# 	# filter each MG by marker gene using BBMap software (WAY faster than my own scripts)
# 	filterbyname.sh \
# 	in=$READDIR/$metagenome.filter.total.faa \
# 	out=${metagenome}.curtomarkers.faa \
# 	names=${metagenome}.temp.txt ow=t include=t 2>/dev/null
# 	rm ${metagenome}.temp.txt

	# # rename each read to reflect the MGID on the fasta header and index each sequence
	# cat $f | sed -e "s/^>/>${metagenome}_/" | sed '/^>/ s/ .*//' > ${metagenome}.temp.faa

# done

# cat *.temp.faa > total.markers.faa 
# rm -f *.temp.faa


# for f in *.blast.total.txt
# do

# 	metagenome=${f%.blast.total.txt}

# 	cat $f | sed -e "s/^/${metagenome}_/" | cut -f1-2 > $metagenome.blast.temp.txt

# done

# cat *.blast.temp.txt > total.blastcore.txt 
# rm -f *.blast.temp.txt

# now subset by each marker gene

cd $OUTDIR

while read protein
do

	echo "Processing ${protein}..."

	
	if [ ! -f "${protein}.blast.faa" ]
	then

		grep "${protein}\b" total.blastcore.txt | cut -f1 | sort -u | sed 's/[ ]*$//' > ${protein}.temp.txt

		# filter each MG by marker gene using BBMap software (WAY faster than my own scripts)
		filterbyname.sh \
		in=total.markers.faa  \
		out=${protein}.temp.faa \
		names=${protein}.temp.txt ow=t include=t 2>/dev/null
		rm ${protein}.temp.txt

		# rename each read to index each sequence
		cat ${protein}.temp.faa | awk '/^>/ {$0=NR"_"$0} 1' | sed 's/ //g; s/>//g' | awk '/^[1-9]/ {$0=">"$0} 1' > ${protein}.blast.faa
		rm ${protein}.temp.faa

	else 
		continue 
	fi

done < $CURTODIR/coregenes.txt



