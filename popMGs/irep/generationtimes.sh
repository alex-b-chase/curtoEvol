#!/bin/bash

BASE=/Volumes/Bennett_BACKUP/Research/curtobacterium/curto-evolution
BASEDIR=$BASE/MGprocess/irep 

GENOMEDIR=$BASE/reference_genome
REFGENOME=$GENOMEDIR/MMLR14002_pacbio-corrected.fa

### estimate the generation times using iRep for the metagenomic samples
### we have 2 sets of data
### 1. t3 for all transplanted isogenic bags (N=24) that were sequenced on a HiSeq (in MGprocess)
### 2. we deeply sequenced a single bag from all timepoints (1 bag x 5 sites x 3 timepoints)
### samples from #1 are site-MGN, samples from #2 are siteN-TN 

# iRep is a method for determining replication rates for bacteria from single time point metagenomics sequencing and draft-quality genomes
# bPTR -> measure replication rates using complete genome sequences (modified from Korem et al. Science 2015)

bPTR -f $REFGENOME \
-s /Volumes/JensenLabMGs/alex_alyssa/temp/*.sortedSH.sam \
-o elevation.irep.tsv \
-plot elevation.irep.pdf \
-m coverage


