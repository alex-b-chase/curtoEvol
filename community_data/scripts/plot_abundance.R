setwd('/Volumes/Bennett_BACKUP/Research/curtobacterium/chapter5_curto-evolution/transplant_communitybags/metagenomic_analyses/community_markers/pplacer_output/pplacer_analysis')

library(ggplot2)
library(reshape2) # for melt

abund <- read.table('phylaXsite_abundance_matrix.txt', header = T, sep = '\t')

abundm <- melt(abund, "bag")

mappingfile <- read.table('../../../../mappingFile.txt', header = T, sep = '\t')

abundmtotal <- merge(abundm, mappingfile, by = "bag")

abundmtotal$siteF = factor(abundmtotal$site, levels=c('D', 'W', 'G', 'P', 'S'))
abundmtotal$inoculumF = factor(abundmtotal$inoculum, levels=c('desert', 'scrubland', 'grassland', 'oak-pine', 'subalpine'))

abundm$variable = factor(abundm$variable, levels = c("Proteobacteria", "Actinobacteria", "Bacteroidetes", "Acidobacteria", "other", "NA."))

# by each genera
gabund <- read.table('genusXsite_abundance_matrix.txt', header = T, sep = '\t')

gabundm <- melt(gabund, "site")

gabundm$gradient <- lapply(strsplit(as.character(gabundm$site), "-"), "[", 1)
gabundm$timepoint <- lapply(strsplit(as.character(gabundm$site), "-"), "[", 2)
gabundm$plot <- lapply(strsplit(as.character(gabundm$site), "-"), "[", 3)

gabundm$elev = factor(gabundm$gradient, levels=c('Subalpine', 'Pine', 'Woodland','Grassland', 'Desert', 'Salton'))
gabundm$variable = factor(gabundm$variable, levels = c("Sphingomonas..a.Proteo.", "Methylobacterium..a.Proteo.", "Massilia..b.Proteo.", 
                                                       "Pseudomonas..g.Proteo.", "Halomonas..g.Proteo.",
                                                     "Collinsella..Actino.", "Microbacterium..Actino.", "Modestobacter..Actino.", 
                                                    "Curtobacterium..Actino.", "Blastococcus..Actino..", "Nocardioides..Actino.",
                                                    "Pedobacter..Bact.", "Granulicella..Acido.", "Terriglobus..Acido.", "Other", "NA."))


# read in genus colors
gcolors <- read.table('genera_colors.txt', header =T, sep ='\t', stringsAsFactors = T, comment.char = '')
# create new color pallette for ggplot
gtaxcolors <- setNames(as.character(gcolors$color), gcolors$variable)

# read in phyla colors
colors <- read.table('phyla_colors.txt', header =T, sep ='\t', stringsAsFactors = T, comment.char = '')
# create new color pallette for ggplot
taxcolors <- setNames(as.character(colors$color), colors$variable)

meanabundm <- aggregate(value ~ variable + treatment + inoculumF + siteF, abundmtotal, mean)

# barplot of transplant and inoculum

pdf('phylaXtransplant.pdf', width = 14, height = 16)

ggplot(meanabundm, aes(y = value, x = treatment, fill = variable)) +
  geom_bar(stat = 'identity', position = 'stack') + facet_grid(inoculumF ~ siteF) +
  scale_fill_manual(values = taxcolors, name = "") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(x = "", y = "% Abundance") 

dev.off()

initialsamps <- meanabundm[meanabundm$treatment == "inoculum", ]
ggplot(initialsamps, aes(y = value, x = siteF, fill = variable)) +
  geom_bar(stat = 'identity', position = 'stack') +
  scale_fill_manual(values = taxcolors, name = "") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(x = "", y = "% Abundance") 

# PCO plot of abundances by inoculum and transplant
# read in the edgenum X site

abund <- read.table('edgenumXsite_abundance_matrix.txt', header = T, sep = '\t')

library(vegan)

# remove survey samples
abundsub <- abund[1:99,]

row.names(abund) <- abund$bag
abund$bag <- NULL

braycurt <- vegdist(abund, method = "bray")
hClustering <- hclust(braycurt, method = 'complete')
plot(hClustering, hang = -1)

row.names(abundsub) <- abundsub$bag
abundsub$bag <- NULL
dist.mat <- vegdist(abundsub, method = "bray")

mappingfile2 <- read.table('../../../../mappingFile.txt', header = T, sep = '\t', row.names = 1, comment.char = "")
mappingfile2sub <- mappingfile2[mappingfile2$treatment == "transplant", ]

perm.log <- adonis2(dist.mat ~ site * inoculum, data=mappingfile2sub, permutations = 999, method = "bray", strata = "PLOT")
perm.log

pops <- mappingfile2sub[,c(2,3)]
library(ecodist)
mgroup(dist.mat, pops, nperm=9999)
##       nclust   mantelr       pval
## inocula    5 0.0567048 0.00080008
## transplant 5 0.5383610 0.00010001

sol <- metaMDS(abund, distance = "bray", k = 2, trymax = 50)

row.names(mappingfile) <- mappingfile$bag

NMDS1 = data.frame(x = sol$point[, 1], y = sol$point[, 2])
NMDS2 = mappingfile[, c(2,3,4)]

NMDS <- merge(NMDS1, NMDS2, by = "row.names")

row.names(NMDS) <- NMDS$Row.names
NMDS$Row.names <- NULL

# subset only transplant bags
NMDSsub <- NMDS[NMDS$treatment == "transplant", ]


shape_values <- c(17, 4, 15, 3, 18)
color_values <- c("#ef8080", "#8fed8f", "#0000ff", "#9f20ef", "#edc800")

p1 <- ggplot(data=NMDSsub, aes(x, y, colour = site)) +
  geom_point(aes(shape = inoculum), size = 8, stroke = 2) +
  scale_shape_manual(values = shape_values) +
  scale_color_manual(values = color_values) +
  theme_bw()

p2 <- ggplot(data=NMDSsub, aes(x, y, colour = inoculum)) +
  geom_point(aes(shape = site), size = 6, stroke = 2) +
  scale_shape_manual(values = shape_values) +
  scale_color_manual(values = color_values) +
  theme_bw()

pdf("MDS-symbol.pdf", height = 8, width = 10)

p1

dev.off()



# add in the inoculum samples
# shape_values2 <- c(0, 2)
# color_values2 <- c("red", "green", "blue", "purple", "orange")
# 
# p2 <- ggplot(data=NMDS, aes(x, y, colour = site)) +
#   geom_point(aes(shape = treatment), size = 6.1, stroke = 2) +
#   scale_shape_manual(values = shape_values2) +
#   scale_color_manual(values = color_values2) +
#   theme_bw() +
#   ggtitle("Site") +
#   theme(legend.position="none") 


require(gridExtra)

pdf("MDSsubclade.pdf", height = 8, width = 14)

grid.arrange(p1, p2, ncol=2)

dev.off()

