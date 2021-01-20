#!/bin/bash

# so this is the script used to process the metagenome files after the BLAST results are outputted
# BLAST output files will be *.blast.total.txt

# for this, we used Prodigal -meta to translate and used the refDB from curated genomic dataset of >5000 genomes

# additionally, you will have had to create a local database using the genomes that are included in the PPlacer reference tree
# reference database can be found in $REFDB
# this is the output of blast_hmmbuild.sh


# the reference DB we aare going to use is the one I constructed for the elevation gradient
# specifically, we will use the PATRIC representative genomes + Curto isolates + Sydney genomes from LR and elevation


initial_start=$(date +%s)

# source the appropriate directories for the files
BASEDIR=/Volumes/Bennett_BACKUP/Research/curtobacterium/curto-evolution/transplant_communitybags
MGDIR=$BASEDIR/metagenomic_analyses
AAREADS=$MGDIR/prodigal_reads
BLASTOUT=$MGDIR/blast_output

# reference the marker genes used and the refDB
REFDB=/Volumes/Bennett_BACKUP/Research/reference_db
BLASTDB=$REFDB/blastDB
ESSEN=$REFDB/essen_singlecopy
HMMPROF=$REFDB/hmmprofiles

# subset the marker genes used 
refmgs=$ESSEN/Essen_SingleCopy.txt
cat $refmgs | awk '{if ($6 == "1") print $0;}' | \
cut -f7 | cut -f2-3 -d' ' | sed 's/ /_/g' | cut -f1 -d'-' > $MGDIR/marker_genes.txt

markergene=$MGDIR/marker_genes.txt
mappingfile=$BASEDIR/mappingMGFile.txt

# rm -f $MGDIR/summary_mgs.txt
# rm -rf $MGDIR/markergene_subset
mkdir -p $MGDIR/markergene_subset

################################################################################################
#### subset by each protein by looping through each MG library and subsetting protein BLAST hits
################################################################################################

cd $BLASTOUT

ethresh=$1


# generate a summary file for the number of marker genes per library
echo "metagenome:num.AAreads:num.markers.BLAST" | tr ':' '\t' > $MGDIR/summary_mgs.txt

while read protein
do

	echo "############################################################################"
	echo "Searching for protein ${protein} in the Metagenomes."

	initial_start2=$(date +%s)

	for f in *.blast.total.txt
	do

		metagenome=${f%.blast.total.txt}

		if [[ $f == Z* ]]
		then
			READDIR=$AAREADS/time0
		else 
			READDIR=$AAREADS/time1
		fi

		# get basic stats for each library
		blastcount=$(cat $f | cut -f1 | sort | uniq | sed 's/[ ]*$//' | wc -l)
		readcount=$(cat $READDIR/$metagenome.filter.total.faa | grep '>' | wc -l)
		echo "${metagenome}:${readcount}:${blastcount}" | tr ':' '\t' >> $MGDIR/summary_mgs.txt

		cat $f | grep ${protein}p | cut -f1 | sort | uniq | sed 's/[ ]*$//' > ${metagenome}.${protein}p.temp.txt

		# filter each MG by marker gene using BBMap software (WAY faster than my own scripts)
		filterbyname.sh \
		in=$READDIR/$metagenome.filter.total.faa \
		out=${metagenome}.${protein}p.temp.faa \
		names=${metagenome}.${protein}p.temp.txt ow=t include=t 2>/dev/null
		rm ${metagenome}.${protein}p.temp.txt

		# rename each read to reflect the MGID on the fasta header
		cat ${metagenome}.${protein}p.temp.faa | sed -e "s/^>/>${metagenome}_/" > ${metagenome}.${protein}p.subset.faa
		rm ${metagenome}.${protein}p.temp.faa

	done

	# now, that we have gone through each MG for one marker gene, combine results
	cat *${protein}p.subset.faa > ${protein}p.blast.faa
	rm *${protein}p.subset.faa

	# secondary seach using HMMer for each marker gene using protein structure (remember HMM profiles are aligned)
	echo "Done with BLAST subset for ${protein}...starting secondary filter with HMMer"

	hmmsearch --tblout ${protein}p.hmm.txt -E 1e-${ethresh} --cpu 4 $HMMPROF/${protein}p.hmm ${protein}p.blast.faa > log.txt

	cat ${protein}p.hmm.txt | cut -f1 -d' ' | sort | uniq > ${protein}p.temp.txt

	# subset new HMMer filtered reads
	filterbyname.sh \
	in=${protein}p.blast.faa \
	out=total_${protein}p.${ethresh}.hmm.faa \
	names=${protein}p.temp.txt ow=t include=t 2>/dev/null
	rm ${protein}p.temp.txt
	rm ${protein}p.hmm.txt

	# move finished product to different folder
	# THESE ARE THE FILES READY FOR PPLACER!!!
	mv total_${protein}p.${ethresh}.hmm.faa $MGDIR/markergene_subset


	total_end2=$(date +%s)
	total_runtime2=$(echo "$total_end2 - $initial_start2" | bc -l)

	echo "Done with ${protein} - total time ${total_runtime2}"

done < $markergene

# remove duplicate looping of summaryfile
awk '!seen[$0]++' $MGDIR/summary_mgs.txt > $MGDIR/summaryMGs.txt
rm $MGDIR/summary_mgs.txt
rm log.txt


echo "################################################################################################"
echo ""
echo "DONE!"

echo ""
echo "The finished files are labeled 'total_<marker gene>.faa' and are in $MGDIR/markergene_subset"
echo "Each marker gene will have its own file with each sequence labeled by the metagenomic ID"
echo ""
echo "If you wish, you may want to screen these sequences at a higher e-value"
echo "If not, these are ready to be inputted into PPLACER once you align to reference packages"

total_end=$(date +%s)
total_runtime=$(echo "$total_end - $initial_start" | bc -l)
echo "################################################################################################"
echo ""
echo "Total time: $total_runtime seconds"









