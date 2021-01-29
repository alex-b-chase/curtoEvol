rm(list=ls())
setwd('/Volumes/Bennett_BACKUP/Research/curtobacterium/chapter5_curto-evolution/transplant_communitybags/metagenomic_analyses/curto_coregenes/pplacer_output/pplacer_analysis/')

library(ggplot2)
library(reshape2) # for melt

abund <- read.table('subcladeXsite_abundance_matrix_ABSOLUTE.txt', header = T, sep = '\t')
abundm <- melt(abund, "bag")

mappingfile <- read.table('../../../../mappingFile.txt', header = T, sep = '\t')
abundmtotal <- merge(abundm, mappingfile, by = "bag")
abundmtotal$siteF = factor(abundmtotal$site, levels=c('D', 'W', 'G', 'P', 'S'))
abundmtotal$inoculumF = factor(abundmtotal$inoculum, levels=c('desert', 'scrubland', 'grassland', 'oak-pine', 'subalpine'))

# read in colors
colors <- read.table('../clade_colors.txt', header =T, sep ='\t', stringsAsFactors = T, comment.char = '')
# create new color pallette for ggplot
taxcolors <- setNames(as.character(colors$color), colors$variable)

meanabundm <- aggregate(value ~ variable + treatment + inoculumF + siteF, abundmtotal, mean)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# get transplant by itself
countmetaTRANS <- abundmtotal[abundmtotal$treatment != "inoculum",]
countmetaTRANS$logcurto <- log10(countmetaTRANS$value)
summetaTRANS <- summarySE(countmetaTRANS, measurevar = "value", groupvars = c("siteF", "variable"))

# subset out the ecotypes
ecotypes <- c("IA", "IBC", "IVA", "IVB", "IVC", "VA")
summetaTRANSeco <- summetaTRANS[summetaTRANS$variable %in% ecotypes,]

pdf('ecotype_siteavg.pdf', height=5, width=8)
op <- par(mar=c(5, 6, 4, 2) + 0.1)
dodge <- position_dodge(width=0.5)

ggplot(summetaTRANSeco, aes(x = siteF, y = value, group = variable, color = variable)) + 
  geom_line(size = 2, position = dodge) + 
  geom_pointrange(aes(ymin = value - sd, ymax = value + sd), position = dodge, size = 2, shape = 15) +
  scale_y_log10() +
  theme_bw() +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank()) +
  labs(x = "Transplanted Site", y = "Log Cell Counts") +
  scale_fill_manual(values = taxcolors, name = "Curtobacterium Ecotype") +
  scale_color_manual(values = taxcolors, name = "Curtobacterium Ecotype") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(text = element_text(size=20)) +
  theme(legend.title = element_blank())

dev.off()

write.table(summetaTRANSeco, file = "curtoavgabund_postTRANS.txt", quote = F, row.names = F, sep = '\t')

# barplot of transplant and inoculum
pdf('subcladeXtransplant_ABSOLUTE.pdf', width = 14, height = 16)
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

##########################################
# get inocula by itself
##########################################
countmetaINO <- abundmtotal[abundmtotal$treatment == "inoculum",]
countmetaINO$logcurto <- log10(countmetaINO$value)
summetaINO <- summarySE(countmetaINO, measurevar = "value", groupvars = c("siteF", "variable"))
# subset out the ecotypes
summetaINOeco <- summetaINO[summetaINO$variable %in% ecotypes,]

pdf('ecotype_inoavg.pdf', height=5, width=8)
op <- par(mar=c(5, 6, 4, 2) + 0.1)
dodge <- position_dodge(width=0.5)

ggplot(summetaINOeco, aes(x = siteF, y = value, group = variable, color = variable)) + 
  geom_line(size = 2, position = dodge) + 
  geom_pointrange(aes(ymin = value - sd, ymax = value + sd), position = dodge, size = 2, shape = 15) +
  scale_y_log10() +
  theme_bw() +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank()) +
  labs(x = "Transplanted Site", y = "Log Cell Counts") +
  scale_fill_manual(values = taxcolors, name = "Curtobacterium Ecotype") +
  scale_color_manual(values = taxcolors, name = "Curtobacterium Ecotype") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(text = element_text(size=20)) +
  theme(legend.title = element_blank())

dev.off()

write.table(summetaINOeco, file = "curtoavgabund_postINO.txt", quote = F, row.names = F, sep = '\t')


##########################################
##########################################
# try and get the abundance by site
##########################################
##########################################
library(ggridges)
databas<-read.csv(text="site,sitenum
D,1
W,2
G,3
P,4
S,5")

totalabund <- merge(abundmtotal, databas, by = "site", all = T)
ecotypes <- c("IA", "IBC", "IVA", "IVB", "IVC", "VA")

library(dplyr)
totalabund2 <- totalabund %>%
  filter(variable %in% ecotypes)

totalabund2$sitebag <- as.numeric(paste(totalabund2$sitenum, totalabund2$bagnum, sep = '.'))
inoculumsub <- totalabund2[totalabund2$treatment == "inoculum", ]
transplantsub <- totalabund2[totalabund2$treatment == "transplant", ]

theme_set(theme_ridges())
require(scales)

pdf('density_subcladeXinoculum_ABSOLUTE.pdf', width = 14, height = 16)
ggplot(inoculumsub, aes(x = sqrt(value), y = siteF, color = variable, fill = variable)) + 
  geom_density_ridges(scale = .9, size = .1) +
  scale_fill_manual(values = c("#7ecfe8", "#004fff", "#eaa9ff", "#973eff", "#470099", "#00ff36")) +
  scale_color_manual(values = c("#7ecfe8", "#004fff", "#eaa9ff", "#973eff", "#470099", "#00ff36"))  +
  scale_x_log10()
dev.off()

pdf('density_subcladeXtransplant_ABSOLUTE.pdf', width = 14, height = 16)
ggplot(transplantsub, aes(x = sqrt(value), y = siteF, color = variable, fill = variable)) + 
  geom_density_ridges(scale = .9, size = .1) +
  scale_fill_manual(values = c("#7ecfe8", "#004fff", "#eaa9ff", "#973eff", "#470099", "#00ff36")) +
  scale_color_manual(values = c("#7ecfe8", "#004fff", "#eaa9ff", "#973eff", "#470099", "#00ff36"))
dev.off()

pdf('density_subcladeXtransplant-inocolum_ABSOLUTE.pdf', width = 14, height = 16)
ggplot(transplantsub, aes(x = sqrt(value), y = inoculumF, color = variable, fill = variable)) + 
  geom_density_ridges(scale = .9, size = .1) +
  scale_fill_manual(values = c("#7ecfe8", "#004fff", "#eaa9ff", "#973eff", "#470099", "#00ff36")) +
  scale_color_manual(values = c("#7ecfe8", "#004fff", "#eaa9ff", "#973eff", "#470099", "#00ff36")) +
  facet_grid(~ siteF) +
  scale_x_continuous(breaks = c(300, 1000, 3000, 10000, 20000))
dev.off()

rm(list=ls())
setwd('/Volumes/Bennett_BACKUP/Research/curtobacterium/chapter5_curto-evolution/transplant_communitybags/metagenomic_analyses/curto_coregenes/pplacer_output/pplacer_analysis/')

library(vegan)
library(ggplot2)
library(ecodist)

# read in abundnace table
abund <- read.table('subcladeXsite_abundance_matrix_ABSOLUTE.txt', header = T, sep = '\t')
row.names(abund) <- abund$bag
abund$bag <- NULL
# remove survey samples
abundsub <- abund[1:99,]

### data transformation and normalization 
tflex <- abundsub
logabundm <- log(tflex, 10)

### create the distance matrix
vegabundm <- vegdist(logabundm, method = "euclidean", diag = T)
vegabundmdist <- hclust(vegabundm)
plot(vegabundmdist)

hist(t(logabundm))

### test for additional data normalization
scaletflex <- scale(tflex)
hist(t(scaletflex))
vegabundm2 <- vegdist(scaletflex, method = "euclidean", diag = T)

# PERMANOVA stats analysis
# need to average features by strain before we can run
# read in metadata to plot the PCoA graph
metadata <- read.table("../../../../mappingFile.txt", header = T, sep = '\t', row.names = 1, comment.char = "")
metatrans <- metadata[metadata$treatment == "transplant",]
perm.log <- adonis2(vegabundm ~ site * inoculum, data=metatrans, permutations = 999, method = "euclidean", strata = "PLOT")
perm.log

perm.tflexscale <- adonis2(vegabundm2 ~ site * inoculum, data=metatrans, permutations = 999, method = "euclidean", strata = "PLOT")
perm.tflexscale

perm.notrans <- adonis2(tflex ~ site * inoculum, data=metatrans, permutations = 999, method = "bray", strata = "PLOT")
perm.notrans


################################################
##########        PCA plot       ###############
################################################
sitepca <- prcomp(t(logabundm), center = T, scale = F)
sitepcaresults <- summary(sitepca)
sitepcaresults$importance[3,1:3] 

sitepcadata <- as.data.frame(sitepcaresults$rotation)
sitepcadataPC12 <- sitepcadata[, c(1:3)]

sitepcameta <- merge(metadata, sitepcadataPC12, by = 0, all = F)
color_values <- c("red", "green", "blue", "purple", "orange")
shape_values <- c(17, 4, 15, 3, 18)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

summPC1 <- summarySE(sitepcameta, measurevar = "PC1", groupvars = c("site", "inoculum"))
summPC1$PC1sd <- summPC1$sd
summPC2 <- summarySE(sitepcameta, measurevar = "PC2", groupvars = c("site", "inoculum"))
summPC2$PC2sd <- summPC2$sd

totalPC12temp <- merge(summPC1, summPC2, by = c("site", "inoculum"))
totalPC12 <- totalPC12temp[, c("site", "PC1", "PC2", "PC1sd", "PC2sd", "inoculum")]

pdf("PCoA-absolute-averagesiteXinoculum.pdf", height = 10, width = 12)
ggplot(totalPC12, aes(PC1, PC2)) +
  geom_errorbarh(aes(xmin = PC1 - PC1sd, xmax = PC1 + PC1sd), height = 0) +
  geom_errorbar(aes(ymin = PC2 - PC2sd, ymax = PC2 + PC2sd), width = 0) +
  geom_point(aes(color = site, shape = inoculum), size = 10, stroke = 3) +
  theme_bw() +
  scale_colour_manual(values = color_values) +
  scale_shape_manual(values = shape_values)
dev.off()

################################################
##########        MDS plot        ##############
################################################
library(vegan)

braycurt <- vegdist(abundsub, method = "bray")
sol <- metaMDS(abundsub, distance = "bray", k = 2, trymax = 500, autotransform = T)

brayabundmdist <- hclust(braycurt)
plot(brayabundmdist)

hist(t(braycurt))

NMDS1 = data.frame(x = sol$point[, 1], y = sol$point[, 2])
NMDS2 = metadata[, c(1,2,3)]
NMDS <- merge(NMDS1, NMDS2, by = "row.names")
row.names(NMDS) <- NMDS$Row.names
NMDS$Row.names <- NULL

# only want confience intervals around the sample that stayed in the same site
desdes <- NMDS[NMDS$inoculum == "desert" & NMDS$site == "D",]
scrubscrub <- NMDS[NMDS$inoculum == "scrubland" & NMDS$site == "W",]
grassgrass <- NMDS[NMDS$inoculum == "grassland" & NMDS$site == "G",]
pinepine <- NMDS[NMDS$inoculum == "oak-pine" & NMDS$site == "P",]
subsub <- NMDS[NMDS$inoculum == "subalpine" & NMDS$site == "S",]

# Fit vectors to ordination
fit <- envfit(sol, abundsub, permutations = 999)
ef.df <- as.data.frame(fit$vectors$arrows * sqrt(fit$vectors$r))
ef.df$variable <- rownames(ef.df)

A <- as.list(fit$vectors)
#creating the dataframe
pvals <- as.data.frame(A$pvals)
arrows <- as.data.frame(A$arrows * sqrt(A$r))
C <- cbind(arrows, pvals)
#subset
Cred <- subset(C, pvals < 0.05)
Cred <- cbind(Cred, Species = rownames(Cred))

# relevant ecotypes
ecotypes <- c("IA", "IBC", "IVA", "IVB", "IVC", "VA")
arrow.p2 <- Cred %>%
  filter(Species %in% ecotypes)

shape_values <- c(17, 4, 15, 3, 18)
color_values <- c("red", "green", "blue", "purple", "orange")

pdf("MDS-absolute-siteXinoculum.pdf", height = 10, width = 12)
ggplot(NMDS, aes(x, y)) +
  geom_point(aes(shape = inoculum, color = site), size = 6, stroke = 2) + ##separates overlapping points
  stat_ellipse(data = desdes, aes(x, y), type = 't', size = 1, level = 0.8, color = "red") + ##draws 80% confidence interval ellipses
  stat_ellipse(data = scrubscrub, aes(x, y), type = 't', size = 1, level = 0.8, color = "orange") + 
  stat_ellipse(data = grassgrass, aes(x, y), type = 't', size = 1, level = 0.8, color = "green") + 
  stat_ellipse(data = pinepine, aes(x, y), type = 't', size = 1, level = 0.8, color = "blue") + 
  stat_ellipse(data = subsub, aes(x, y), type = 't', size = 1, level = 0.8, color = "purple") + 
  geom_segment(data = arrow.p2, aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = arrow(length = unit(0, "cm")), group = arrow.p2$Species,
               color = c("#7ECFE8", "#0000E8", "#E6C4FF", "#9300FF", "#4B0082", "#00FF36")) +
  theme_bw() +
  scale_shape_manual(values = shape_values) +
  scale_color_manual(values = color_values) +
  coord_fixed()
dev.off()

################################################
##########        PCoA plot        #############
################################################
library("FactoMineR")
library("factoextra")

## datawithoutVF = logabundm[ !(row.names(logabundm) %in% c("3W03")), ]
datawithoutVF = logabundm
scaletflex <- scale(datawithoutVF)
hist(t(scaletflex))

metaPCA <- metatrans
metaPCA$siteinoc <- paste(metaPCA$site, metaPCA$inoculum, sep = "_")
metaPCA$sampleID <- rownames(metaPCA)
metaPCA <- metaPCA[,c(6,5)]

logabundm2 <- merge(datawithoutVF, metaPCA, by = 0)
rownames(logabundm2) <- logabundm2$Row.names
logabundm2$Row.names <- NULL

site.pca <- PCA(logabundm2[,c(1:14)], graph = FALSE)
fviz_contrib(site.pca, choice = "var", axes = 1:2, top = 10)


pdf("PCA-absolute-siteXinoculum.pdf", height = 10, width = 18 )
fviz_pca_biplot(site.pca, geom.ind = "point",
                fill.ind = logabundm2$siteinoc, col.ind = "black",
                pointshape = 24, pointsize = 3,
                palette = c("#FF0000", "#cc4ec7", "#d9433a", "#c84b82", "#75da53", 
                            "#6c41c8", "#00FF00", "#d1db4b", "#649841", "#542c6e", 
                            "#67d5a7", "#c84b82", "#0000FF", "#c89e3d", "#7d7cd1",
                            "#d27a4e", "#80becf", "#833734", "#DEB0FF", "#cad09e", 
                            "#39252e", "#ce9cbc", "#3d5f3e", "#FFA500", "#87704a"), 
                addEllipses = TRUE, ellipse.level = 0.75, 
                select.var = list(name = c("IA", "IBC", "IVA", "IVB", "IVC", "VA")),
                label = "var", 
                repel = TRUE, mean.point = FALSE,
                alpha.var = "contrib", col.var = "contrib",
                legend.title = "Species") 
dev.off()

library(ape)

vegabundm <- vegdist(datawithoutVF, method = "euclidean", diag = T)
PCOA <- pcoa(vegabundm)
# plot the eigenvalues and interpret
barplot(PCOA$values$Relative_eig[1:10])
PCOA$values$Relative_eig[1:3]
# Can you also calculate the cumulative explained variance of the first 3 axes?

# Plot your results
biplot(PCOA)

# You see what`s missing? 
# Indeed, there are no species plotted on this biplot. 
# That's because we used a dissimilarity matrix (sites x sites) 
# as input for the PCOA function. 
# Hence, no species scores could be calculated. 
# However, we could work around this problem like this:
bipcoa <- biplot.pcoa(PCOA, datawithoutVF)


PCOAaxes <- PCOA$vectors[,c(1,2)]
PCOAmeta = metadata[, c(1,2,3)]

PCOAplot <- merge(PCOAaxes, PCOAmeta, by = "row.names")

# should look the same as the above one, just without the nice graphing stuff
ggplot(data = PCOAplot, aes(Axis.1, Axis.2, colour = site)) +
  geom_point(aes(shape = inoculum), size = 6, stroke = 2) +
  scale_shape_manual(values = shape_values) +
  scale_color_manual(values = color_values) +
  theme_bw()




