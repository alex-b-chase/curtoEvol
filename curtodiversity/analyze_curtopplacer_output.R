setwd('/Volumes/Bennett_BACKUP/Research/curtobacterium/curto-evolution/transplant_communitybags/metagenomic_analyses/curto_coregenes/pplacer_output')

# read in pplacer output files .csv
# total them in terminal with and remove duplicate header row
# subset the dataframe to make it a little easier
# keep marker gene, read ID, node information, and taxonomy

# system("cat individual_genes/*.csv | cut -f1,2,4,11 -d',' | awk -F , '$4 == \"2034\" { print }' > total_markers.fix.csv", intern = TRUE)
# cat individual_genes/*.csv | cut -f1,2,4,11 -d',' | awk -F , '$4 == "2034" { print }' > total_markers.fix.removebad.csv

# total_markers <- read.csv('total_markers.fix.removebad.csv')
total_markers <- read.csv('total_markers.fix.csv')

colnames(total_markers) <- c("origin","name","edge_num","classification")


# subset only the Curtobacterium taxID = '2034' 

newtotal <- total_markers[total_markers$classification == 2034 ,]

# some sequences might have identical prefixed names from the .fastq files, BUT
# indexed each sequence so shouldn't be too big of an issue...

# add column for site so we can table it later
newtotal$bag <- lapply(strsplit(as.character(newtotal$name), "_"), "[", 2)

# rename the column to match
names(newtotal)[names(newtotal) == 'classification'] <- 'tax_id'

# $newtotal$site is a list and will not write file
newtotal$bag <- vapply(newtotal$bag, paste, collapse = ", ", character(1L))

totalmarks <- t(table(newtotal$bag))
write.table(totalmarks, file = "totalmarkersXsite.txt", quote = F, sep = '\t', row.names = F)

# subset by clades of curtobacterium
# need to look these up in the .xml files
# ONLY NEED TO DO THIS ONCE SINCE THE TREE STRUCTURE SHOULD BE IDENTICAL

clade <- read.table('strain2clade.txt', header = T, sep = '\t')
names(clade)[names(clade) == 'nodeID'] <- 'edge_num'

clade2 <- merge(newtotal, clade, by = "edge_num", all = F)
cladest <- table(clade2$site, clade2$clade, exclude = NULL)
write.table(cladest, file = "siteXclade.txt", quote = F, sep = '\t', col.names=NA)

# table of clade by marker gene
clade2$subclade <- droplevels(clade2$subclade)
markerC <- table(clade2$clade, clade2$origin, exclude = NULL)
table(clade2$origin)


subcladeXsite <- table(clade2$subclade, unlist(clade2$site), exclude = T)
cladeXsite <- table(clade2$clade, unlist(clade2$site), exclude = T)

write.table(subcladeXsite, file = "siteXsubclade.txt", quote = F, sep = '\t', col.names=NA)
write.table(cladeXsite, file = "siteXclade.txt", quote = F, sep = '\t', col.names=NA)

# check to see if some genes are biased towards phylogenies: subclade X origin
clade2$subclade <- droplevels(clade2$subclade)
cladeXgene <- t(table(unlist(clade2$subclade), unlist(clade2$origin), exclude = T))
geneXclade <- t(table(unlist(clade2$origin), unlist(clade2$subclade), exclude = T))
genescaled <- as.matrix(scale(cladeXgene))
write.table(cladeXgene, file = "geneXsubclade.txt", quote = )
stdize = function(x, ...) {(x - min(x, ...)) / (max(x, ...) - min(x, ...))}
genescaled <- as.matrix(stdize(cladeXgene))  # lots of NA


# create heatmap and don't reorder columns

library(gplots)
pdf('geneXsubclade.pdf', height = 50, width = 20)

subclades <- as.data.frame.matrix(cladeXgene)
subclade.scale <- scale(subclades)
heatmap.2(subclade.scale)

dev.off()

pdf('subcladeXgene.pdf', height = 50, width = 20)

subclades <- as.data.frame.matrix(geneXclade)
subclade.scale <- scale(subclades)
heatmap(subclade.scale)

dev.off()

# example of a removed gene - almost all clade VA
# gapA.20 <- clade2[clade2$origin == "gapA.20" , ]

# remove bad genes and regenerate
