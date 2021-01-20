 # write PBS shell scripts for SPAdes assembler to run on the HPC

GenomeMetadata = read.table("/Volumes/Bennett_BACKUP/Research/curtobacterium/curto-evolution/isolate_genomics_final/spades/reads.txt", header = T, sep = "\t", quote = "", comment.char = "", stringsAsFactors = F)

readFiles = GenomeMetadata[, c("FFile", "RFile")]
readIds = GenomeMetadata[, "fastq.id"]

# bbmap parameters
threads = 8

refBase = "/bio/abchase/genomes/evolution_complete"
outBase = "$READBASE/spades"
bin = "/data/users/abchase/bin"
blast = "/data/apps/commondata/blastdb"
tax = "/dfs1/bio/abchase/taxdump_ncbi"

setwd("/Volumes/Bennett_BACKUP/Research/curtobacterium/curto-evolution/isolate_genomics_final/spades/hpcjobs")

for(k in 1:nrow(readFiles)){
			jobName = paste(readIds[k], "_SPAdes", sep = "")
			if(nchar(jobName) > 15){
				jobFile = paste(jobName, ".sh", sep = "")
				jobName = substr(jobName, 1, 15)
			}else{
				jobFile = paste(jobName, ".sh", sep = "")
			}
			system(paste("rm -f ", jobFile, sep = ""))

			cat(k, jobFile, "\n")
			
			cat("#!/bin/bash", "\n", file = jobFile, sep = "", append = T)
			cat("#$ -N ", paste(jobName, sep = ""), "\n", file = jobFile, sep = "", append = T)
			cat("#$ -q mic,pub8i", "\n", file = jobFile, sep = "", append = T)
			cat("#$ -pe openmp ", threads, "\n", file = jobFile, sep = "", append = T)
	
			cat("\n", file = jobFile, sep = "", append = T)

			cat("echo \"Job started on `hostname` at `date`\"", "\n", file = jobFile, sep = "", append = T)
		
			cat("module load SPAdes/3.8.2", "\n", file = jobFile, sep = "", append = T)
			cat("module load blast/2.2.30", "\n", file = jobFile, sep = "", append = T)
			cat("module load bowtie2/2.2.7", "\n", file = jobFile, sep = "", append = T)
			cat("module load samtools/1.3", "\n", file = jobFile, sep = "", append = T)
			cat("module load BBMap/37.50", "\n", file = jobFile, sep = "", append = T)
		
			cat("\n", file = jobFile, sep = "", append = T)

			cat(paste("REF=", readIds[k], sep = ""), "\n", file = jobFile, sep = "", append = T)
			cat(paste("READBASE=", refBase, sep = ""), "\n", file = jobFile, sep = "", append = T)
			cat(paste("BBMAPDIR=/data/users/abchase", sep = ""), "\n", file = jobFile, sep = "", append = T)
			cat(paste("FFILE=$READBASE/genome_raw/", readFiles[k, "FFile"], sep = ""), "\n", file = jobFile, sep = "", append = T)
			cat(paste("RFILE=$READBASE/genome_raw/", readFiles[k, "RFile"], sep = ""), "\n", file = jobFile, sep = "", append = T)
			
			cat("\n", file = jobFile, sep = "", append = T)
			
			cat(paste("OUTDIR=", outBase, "/$REF", sep = ""), "\n", file = jobFile, sep = "", append = T)
			cat(paste("THREAD=", threads, sep = ""), "\n", file = jobFile, sep = "", append = T)
			
			cat("\n", file = jobFile, sep = "", append = T)
			
			line2write = "bbduk.sh in1=$FFILE in2=$RFILE \\"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			line2write = "out1=$REF.R1.clean.fq.gz out2=$REF.R2.clean.fq.gz \\"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			line2write = "minlen=25 qtrim=rl trimq=10 ktrim=r k=25 ref=$BBMAPDIR/nextera.fa.gz hdist=1"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			
			cat("\n", file = jobFile, sep = "", append = T)
			
			line2write = "rm -rf $OUTDIR"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			cat("\n", file = jobFile, sep = "", append = T)
			
			line2write = "spades.py -k 31,41,51,61,71,81,91 --careful -t $THREAD \\"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			
			line2write = "-1 $REF.R1.clean.fq.gz -2 $REF.R2.clean.fq.gz \\"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			
			line2write = "-o $OUTDIR"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			
			cat("\n", file = jobFile, sep = "", append = T)
			line2write = "rm -f $REF.R1.clean.fq.gz"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			line2write = "rm -f $REF.R2.clean.fq.gz"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			
			cat("\n", file = jobFile, sep = "", append = T)
			cat("# Now we will check the assembled contigs for any other contamination", "\n", file = jobFile, sep = "", append = T)
			cat("\n", file = jobFile, sep = "", append = T)
			
			cat(paste("BIN=", bin, sep = ""), "\n", file = jobFile, sep = "", append = T)
			cat(paste("BLASTDB=", blast, sep = ""), "\n", file = jobFile, sep = "", append = T)
			cat(paste("TAXDUMP=", tax, sep = ""), "\n", file = jobFile, sep = "", append = T)
			
			cat("\n", file = jobFile, sep = "", append = T)
			
			line2write = "cd $OUTDIR"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			
			cat("\n", file = jobFile, sep = "", append = T)
			
			line2write = "ASSEMBLY=scaffolds.fasta"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)

			cat("\n", file = jobFile, sep = "", append = T)
			
			line2write = "blastn -task megablast -query $ASSEMBLY -db $BLASTDB/nt -evalue 1e-5 -max_target_seqs 1 \\"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			line2write = "-num_threads $THREAD -outfmt '6 qseqid staxids' -out $REF.nt.1e-5.megablast"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			
			cat("\n", file = jobFile, sep = "", append = T)
			
			line2write = "bowtie2-build $ASSEMBLY $ASSEMBLY"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			
			cat("\n", file = jobFile, sep = "", append = T)
			
			line2write = "bowtie2 -x $ASSEMBLY --very-fast-local -k 1 -t -p $THREAD --reorder --mm \\"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			line2write = "-U <($BIN/shuffleSequences_fastx.pl 4 <(zcat $FFILE) <(zcat $RFILE)) \\"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			line2write = "| samtools view -S -b -T $ASSEMBLY - > $REF.bowtie2.bam"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			
			cat("\n", file = jobFile, sep = "", append = T)
			
			line2write = "$BIN/gc_cov_annotate.pl \\"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			line2write = "--blasttaxid $REF.nt.1e-5.megablast \\"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			line2write = "--assembly $ASSEMBLY --bam *.bam --out $REF.bowtie.blobplot.txt \\"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			line2write = "--taxdump $TAXDUMP --taxlist genus family phylum"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			
			cat("\n", file = jobFile, sep = "", append = T)
			
			line2write = "cp $ASSEMBLY $READBASE/spades/$REF.scaffolds.fasta"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			line2write = "cp $REF.bowtie.blobplot.txt $READBASE/spades/$REF.bowtie.blobplot.txt"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			
			cat("\n", file = jobFile, sep = "", append = T)
			
			line2write = "cd $READBASE/spades"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			line2write = "rm -rf $OUTDIR/"
			cat(line2write, "\n", file = jobFile, sep = "", append = T)
			
			cat("\n", file = jobFile, sep = "", append = T)
			cat("echo \"Job Ended at `date`\"", "\n", file = jobFile, sep = "", append = T)
			

}

