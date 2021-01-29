rm(list=ls())
setwd('/Volumes/Bennett_BACKUP/Research/curtobacterium/chapter5_curto-evolution/transplant_communitybags/metagenomic_analyses/curto_coregenes/ecotype_physio/')

traits <- read.table("pheno_assaysXecotype4.txt", header = T, sep = '\t', row.names = 1)

# took the average assay value regardless of temperature for each strain
# took the temperature each strain did the best at each assay

library(gplots)

subclade.scale <- scale(traits)

my_palette <- colorRampPalette(c("#4FAAA7", "#E8DF9C", "#C6460F"))(n = 299)

pdf("ecotype_physio.pdf", height = 6, width = 6)
heatmap.2(subclade.scale,
          density.info="none",  
          trace="none",         
          margins =c(12,9),    
          col=my_palette,       
          dendrogram='none',     
          Rowv = FALSE, Colv = FALSE)
dev.off()


na.rm <- TRUE
scaleCol <- function(x) {
  colMeans <- rm <- colMeans(x, na.rm = na.rm)
  x <- sweep(x, 2, rm)
  colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
  x <- sweep(x, 2, sx, "/")
  return(round(x, 6))
}

library(reshape2)
zscores <- melt(scaleCol(subclade.scale))
colnames(zscores) <- c("variable", "trait", "trait_value")


### community weight value means
### read in ecotype abundances
ecotypeabund <- read.table("../pplacer_output/pplacer_analysis/subcladeXsite_abundance_matrix_ABSOLUTE.txt", header = T, sep = '\t')
# subset only 6 ecotypes with trait data
ecotypelist <- c("bag", "IA", "IBC", "IVA", "IVB", "IVC", "VA")
ecotypeabund2 <- melt(ecotypeabund[, ecotypelist])

traits$variable <- rownames(traits)

transplant.cwm.trait <- merge(ecotypeabund2, traits, by = "variable")
transplant.cwm.trait$logabund <- sqrt(transplant.cwm.trait$value)

library(tidyr)
library(dplyr)

# Calculating CWM using dplyr and tidyr functions
summarize.trans.cwm <-   # New dataframe where we can inspect the result
  transplant.cwm.trait %>%   # First step in the next string of statements
  group_by(bag) %>%   # Groups the summary file by Plot number
  summarize(           # Coding for how we want our CWMs summarized
    Biofilm_cwm = weighted.mean(biofilm, value),   # Actual calculation of CWMs
    umax_cwm = weighted.mean(umax, value),
    Amax_cwm = weighted.mean(Amax, value),
    Cellulose_cwm = weighted.mean(cell_max, value),
    Xylan_cwm = weighted.mean(xy_max, value),
    Temp_cwm = weighted.mean(temp_avg, value)
  )

summarize.trans.cwm <- as.data.frame(summarize.trans.cwm)

## read in metadata
mappingFile <- read.table("../../../mappingFile.txt", sep = '\t', header = T)
mappingFile <- mappingFile[, c("bag", "treatment", "site")]
summarize.trans.cwmM <- merge(summarize.trans.cwm, mappingFile, by = "bag")

library(ggplot2)

summarize.trans.cwm2 <- melt(summarize.trans.cwmM)
summarize.trans.cwm2$siteF = factor(summarize.trans.cwm2$site, levels=c('D', 'W', 'G', 'P', 'S'))

traitcolors <- c("#00725C",
                  "#282828",
                  "#b3943f",
                  "#6295cd",
                  "#ca5d46",
                  "#c85990")
dodge <- position_dodge(width=0.5)

pdf("total.CWM.pdf", width = 12, height = 8)

ggplot(summarize.trans.cwm2, aes(x = siteF, y = value, group = treatment, color = treatment)) +
  geom_point(alpha = 0.4) +
  geom_smooth(size = 2, aes(linetype = treatment)) +
  facet_wrap(~ variable, scales = "free") +
  theme_bw() +
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank()) +
  labs(x = "Transplanted Site", y = "Community Weighted Means") +
  scale_fill_manual(values = traitcolors, name = "Traits") +
  scale_color_manual(values = traitcolors, name = "Traits") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(text = element_text(size=20)) +
  theme(legend.title = element_blank())

dev.off()

## add in MAP and MAT
siteF <- c("D", "W", "G", "P", "S")
mat <- c(28.9, 19.9, 20.9, 13.2, 13.3)
map <- c(213.5, 428.4, 569.4, 1415.8, 1376.5)
mapmat <- data.frame(siteF, map, mat)

totaltraits <- merge(summarize.trans.cwm2, mapmat, by = "siteF")
### r2 values for the traits across the gradient

biofilm <- totaltraits[totaltraits$variable == "Biofilm_cwm" & totaltraits$treatment == "inoculum",]
umax <- totaltraits[totaltraits$variable == "umax_cwm" & totaltraits$treatment == "inoculum",]
Amax <- totaltraits[totaltraits$variable == "Amax_cwm" & totaltraits$treatment == "inoculum",]
temperature <- totaltraits[totaltraits$variable == "Temp_cwm" & totaltraits$treatment == "inoculum",]
xylan <- totaltraits[totaltraits$variable == "Xylan_cwm" & totaltraits$treatment == "inoculum",]
cellulose <- totaltraits[totaltraits$variable == "Cellulose_cwm" & totaltraits$treatment == "inoculum",]

biofilmT <- totaltraits[totaltraits$variable == "Biofilm_cwm" & totaltraits$treatment == "transplant",]
umaxT <- totaltraits[totaltraits$variable == "umax_cwm" & totaltraits$treatment == "transplant",]
AmaxT <- totaltraits[totaltraits$variable == "Amax_cwm" & totaltraits$treatment == "transplant",]
temperatureT <- totaltraits[totaltraits$variable == "Temp_cwm" & totaltraits$treatment == "transplant",]
xylanT <- totaltraits[totaltraits$variable == "Xylan_cwm" & totaltraits$treatment == "transplant",]
celluloseT <- totaltraits[totaltraits$variable == "Cellulose_cwm" & totaltraits$treatment == "transplant",]

ggplot(temperature, aes(x = map, y = value)) + geom_point() + geom_smooth(method = "lm")
ggplot(umax, aes(x = mat, y = value)) + geom_point() + geom_smooth(method = "lm")

## get pearson correlation for each trait by mat and map

## survey samples first
cor.test(biofilm$value, biofilm$map, method = "pearson")
cor.test(biofilm$value, biofilm$mat, method = "pearson")

cor.test(umax$value, umax$map, method = "pearson")
cor.test(umax$value, umax$mat, method = "pearson")

cor.test(Amax$value, Amax$map, method = "pearson")
cor.test(Amax$value, Amax$mat, method = "pearson")

cor.test(temperature$value, temperature$map, method = "pearson")
cor.test(temperature$value, temperature$mat, method = "pearson")

cor.test(xylan$value, xylan$map, method = "pearson")
cor.test(xylan$value, xylan$mat, method = "pearson")

cor.test(cellulose$value, cellulose$map, method = "pearson")
cor.test(cellulose$value, cellulose$mat, method = "pearson")

## transplant samples second
cor.test(biofilmT$value, biofilmT$map, method = "pearson")
cor.test(biofilmT$value, biofilmT$mat, method = "pearson")

cor.test(umaxT$value, umaxT$map, method = "pearson")
cor.test(umaxT$value, umaxT$mat, method = "pearson")

cor.test(AmaxT$value, AmaxT$map, method = "pearson")
cor.test(AmaxT$value, AmaxT$mat, method = "pearson")

cor.test(temperatureT$value, temperatureT$map, method = "pearson")
cor.test(temperatureT$value, temperatureT$mat, method = "pearson")

cor.test(xylanT$value, xylanT$map, method = "pearson")
cor.test(xylanT$value, xylanT$mat, method = "pearson")

cor.test(celluloseT$value, celluloseT$map, method = "pearson")
cor.test(celluloseT$value, celluloseT$mat, method = "pearson")


summary(lm(value ~ map, data=biofilm))
summary(lm(value ~ mat, data=biofilm))
summary(lm(value ~ map, data=biofilmT))
summary(lm(value ~ mat, data=biofilmT))

summary(aov(value ~ map, data=xylan))





library(ape)
library(vegan)

vegabundm <- vegdist(scaletflex, method = "euclidean", diag = T)
PCOA <- pcoa(vegabundm)
# plot the eigenvalues and interpret
barplot(PCOA$values$Relative_eig[1:10])
PCOA$values$Relative_eig[1:3]
# Can you also calculate the cumulative explained variance of the first 3 axes?

# Plot your results
biplot(PCOA)
pdf("pcoa_cwm.pdf", width = 10, height = 10)
bipcoa <- biplot.pcoa(PCOA, scaletflex)
dev.off()

library(ggplot2)

curtoabund <- read.table("../pplacer_output/pplacer_analysis/curtoavgabund_postTRANS.txt", header = T, sep = '\t')
curtoabund2 <- curtoabund[, c(1,2,4)]

totaltraits <- merge(zscores, curtoabund2, by = "variable")
totaltraits$logabund <- log10(totaltraits$value)

ggplot(totaltraits, aes(x = siteF, y = logabund)) +
  geom_point(aes(color = variable, shape = siteF)) +
  facet_grid(~ trait_value) +
  theme_classic()

library(plotly)
plot_ly(totaltraits, x = ~siteF, y = ~logabund, z = ~trait_value, 
        type = "scatter3d", mode = "markers", color = ~variable, symbol = ~trait)


p <- plot_ly(iris, x = ~Sepal.Width, y = ~Sepal.Length) 
add_markers(p, color = ~Petal.Length, size = ~Petal.Length)
add_markers(p, color = ~Species)




### old stuff for 3D



scaleCWM <- summarize.trans.cwm
rownames(scaleCWM) <- paste(scaleCWM$siteF, scaleCWM$treatment, sep = "_")
scaleCWM$siteF <- NULL
scaleCWM$treatment <- NULL

scaleCWM <- as.data.frame(scale(scaleCWM))
scaleCWM$siteF <- rownames(scaleCWM)
scaleCWMm <- melt(scaleCWM)
colnames(scaleCWMm) <- c("siteF", "trait", "zscore")

scaleCWMm$functcat <- "growth"
scaleCWMm$functcat[scaleCWMm$trait == "Biofilm_cwm" | scaleCWMm$trait == "Temp_cwm"] <- "response"
scaleCWMm$functcat[scaleCWMm$trait == "Cellulose_cwm" | scaleCWMm$trait == "Xylan_cwm"] <- "carbon"

scaleCWMm$trait <- NULL
scaleCWM2 <- dcast(scaleCWMm, siteF ~ functcat, value.var = "zscore", fun.aggregate = mean)



scaleCWM2$site <- "zero"
scaleCWM2$site <- ifelse(grepl("D_", scaleCWM2$siteF, ignore.case = T), "D", 
                         ifelse(grepl("W_", scaleCWM2$siteF, ignore.case = T), "W", 
                                ifelse(grepl("S_", scaleCWM2$siteF, ignore.case = T), "S",
                                       ifelse(grepl("G_", scaleCWM2$siteF, ignore.case = T), "G", 
                                              ifelse(grepl("P_", scaleCWM2$siteF, ignore.case = T), "P", 
                                                     "Other")))))
scaleCWM2$treatment <- "zero"
scaleCWM2$treatment <- ifelse(grepl("_survey", scaleCWM2$siteF, ignore.case = T), "survey", 
                              ifelse(grepl("_trans", scaleCWM2$siteF, ignore.case = T), "transplant", 
                                     "Other"))

library(plotly)
plot_ly(scaleCWM2, x = ~growth, y = ~response, z = ~carbon, 
        type = "scatter3d", mode = "markers", color = ~site, symbol = ~treatment)


