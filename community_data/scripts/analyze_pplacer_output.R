setwd('/Volumes/Bennett_BACKUP/Research/curtobacterium/curto-evolution/transplant_communitybags/metagenomic_analyses/community_markers/pplacer_output')

# read in pplacer output files .csv
# total them in terminal with and remove duplicate header row
# subset the dataframe to make it a little easier
# keep marker gene, read ID, node information, and taxonomy

# cat individual_genes/*.csv | cut -f1,2,4,11 -d',' | grep -v 'origin' > total_markers.fix.csv

# total_markers <- read.csv('total_markers.fix.removebad.csv')
total_markers <- read.csv('total_markers.fix.csv', sep = ',', header = F)

colnames(total_markers) <- c("origin","name","edge_num","tax_id")

# some sequences might have identical prefixed names from the .fastq files, BUT
# indexed each sequence so shouldn't be too big of an issue...

# add column for site so we can table it later
total_markers$bag <- lapply(strsplit(as.character(total_markers$name), "_"), "[", 2)


# $total_markers$bag is a list and will not write file
total_markers$bag <- vapply(total_markers$bag, paste, collapse = ", ", character(1L))

totalmarks <- t(table(total_markers$bag))
write.table(totalmarks, file = "totalmarkersXsite.txt", quote = F, sep = '\t', row.names = F)

# now read in the taxonomy file to match to classification ID
tax <- read.csv('/Volumes/Bennett_BACKUP/Research/curtobacterium/chapter3_biogeography/taxonomy/subtax.csv')
subtax <- tax[ , c("tax_id", "parent_id", "tax_name", "phylum", "family", "genus", "class", "order")]

# merge data
taxnew <- merge(total_markers, subtax, by = "tax_id")


# write out to file with all data
write.table(taxnew, file = "total_markers_wTAX.txt", quote = F, row.names = F, sep = '\t')


# check to see if some genes are biased towards phylogenies: subclade X origin
phylaXgene <- t(table(unlist(taxnew$phylum), unlist(taxnew$origin), exclude = NULL))

# create heatmap and don't reorder columns
library(gplots)

heatmap.2(phylaXgene, scale='none')



# table of phylum by site
phyla <- table(taxnew$bag, taxnew$phylum, exclude = NULL)
phylaA <- t(sweep(phyla, 1, rowSums(phyla), '/'))
write.table(phylaA, file = "siteXphylum.txt", quote = F, sep = '\t', col.names=NA)
write.table(phyla, file = "siteXphylum_raw.txt", quote = F, sep = '\t', col.names=NA)
phylaM <- table(taxnew$phylum, taxnew$origin, exclude = NULL)
phylatM <- sweep(phylaM, 1, rowSums(phylaM), '/')
write.table(phylatM, file = "markersXphylum.txt",  quote = F, sep = '\t', col.names=NA)

# table of genera by site
genus <- table(taxnew$bag, taxnew$genus, exclude = NULL)
genusA <- t(sweep(genus, 1, rowSums(genus), '/'))
write.table(genusA, file = "siteXgenus.txt", quote = F, sep = '\t', col.names=NA)
write.table(genus, file = "siteXgenus_raw.txt", quote = F, sep = '\t', col.names=NA)

# table of edge_num by site (for MDS plot?)
edge_num <- table(taxnew$bag, taxnew$edge_num, exclude = NULL)
edge_numA <- t(sweep(edge_num, 1, rowSums(edge_num), '/'))
write.table(edge_numA, file = "siteXedge_num.txt", quote = F, sep = '\t', col.names=NA)
write.table(edge_num, file = "siteXedge_num_raw.txt", quote = F, sep = '\t', col.names=NA)





