setwd("/Users/student/Documents/GitHub/mxMcavSnp")

if (!require("pacman")) install.packages("pacman")
pacman::p_load("adegenet", "dendextend", "flextable", "gdata", "ggdendro", "hierfstat", "Imap", "kableExtra", "paletteer", "patchwork", "officer", "poppr", "RColorBrewer", "reshape2", "StAMPP", "tidyverse", "vcfR", "vegan", "WGCNA", "boa", "plyr", "rgdal", "broom", "rgeos", "ggmap", "moments", "car", "multcompView", "lsmeans", "pairwiseAdonis", "marmap", "ggsn")

options(stringsAsFactors = FALSE) 


devtools::install_github("dkahle/ggmap", force=TRUE)
mxSites = read.csv("mxSampleSites.csv", header=TRUE) 
box <- make_bbox(Lon, Lat, data = mxSites)
calc_zoom(box)


#ggmap::register_google(key = "AIzaSyAgi2MUDvA69dA5nwt1wtaIURGwr6l_xzU", write=TRUE)
baseMap=get_map(location = box, zoom = 9, scale = 2, maptype = c("satellite"), source = c("google"), crop = TRUE, color = c("color"))

#baseMap=get_googlemap(location = box, zoom = 9, size = c(640, 640), scale = 2, maptype = c("satellite"), color = c("color"))
x.min = -89.85185
y.min = 22.35036
x.max = -88.65330
y.max = 23.35262
mxMap=ggmap(baseMap)+ 
  geom_point(data = mxSites, size= 3, shape=21, aes(x = Lon, y = Lat, fill=Location)) + 
  scale_fill_manual(values = c("seagreen3", "blue"))+ 
  #scale_shape_manual(values = c(24,25), name = "Depth Zone", breaks=c("Shallow", "Mesophotic")) + # define shape/color scales 
  scale_y_continuous(label=function(x){return(paste(x,"°N"))}, expand = c(0, 0))+ 
  scale_x_continuous(label=function(x){return(paste(-x,"°W"))}, expand = c(0, 0))+ 
  coord_equal()+
  #ggsn::north(box, scale = 0.5, location= "topleft") + 
  #ggsn::scalebar(box, dist = 25, transform = TRUE, dist_unit="km",model = "WGS84",st.dist = 0.05, location="bottomright")+
  ggsn::north(x.min = -89.85185, y.min = 22.35036, x.max = -88.65330, y.max = 23.35262, scale = 0.25, location= "bottomright", anchor = c(x =  x.max+0.125, y = y.min-0.1)) + 
  ggsn::scalebar(x.min = -89.85185, y.min = 22.35036, x.max = -88.65330, y.max = 23.35262, dist = 25, transform = TRUE, 
                 dist_unit="km", model = "WGS84", st.dist = 0.05, location="bottomright", anchor = c(x =  x.max+0.2, y = y.min-0.175))+
  theme_bw()+ 
  theme(legend.position = "right", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())


mxMap

ggsave("mxMap.tiff", plot = mxMap, height = 5, width = 8, units = "in", dpi = 300)
#x.min = -89.85185, 
#y.min = 22.35036, 
#x.max = -88.65330, 
#y.max = 23.35262,
###########################################################################
#Test Merge Dendrograms
bams=read.table("bams")[,1] # list of bam files
goods=c(1:length(bams))

ma = as.matrix(read.table("mxMerged.ibsMat"))
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.5)  # this shows how similar clones are

ma=ma[goods,goods]
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.5) #
##########

bams=read.table("mxBamsClones")[,1] # list of bam files
goods=c(1:length(bams))

ma = as.matrix(read.table("mxMcavClones.ibsMat"))
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.5)  # this shows how similar clones are

ma=ma[goods,goods]
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.5) 

#####################################################
snpMa = as.matrix(read.table("mxMcavNoClones.ibsMat"))
mxMds = cmdscale(snpMa, eig = TRUE, x.ret = TRUE)

# Determine percent variation captured on each axis
# Calculate the eigenvalues so later we can figure out % variation shown on each Principal Coordinate
mxSnpPcoaVar = round(mxMds$eig/sum(mxMds$eig)*100, 1)
mxSnpPcoaVar

# Format data to plot
mxSnpPcoaValues = mxMds$points
mxSnpPcoaValues

snpI2P = read.csv("mxInds2PopsNoClones.csv") # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
row.names(snpI2P) = snpI2P[,1]
mxSnpPcoaValues=cbind(snpI2P, mxSnpPcoaValues)
mxSnpPcoaValues =as.data.frame(mxSnpPcoaValues, sample = rownames(mxSnpPcoaValues))
colnames(mxSnpPcoaValues)[5] <- "PCo1"
colnames(mxSnpPcoaValues)[6] <- "PCo2"
mxSnpPcoaValues

snpPCoA = merge(mxSnpPcoaValues, aggregate(cbind(mean.x=PCo1,mean.y=PCo2)~popDepth, mxSnpPcoaValues, mean), by="popDepth")

snpPCoA$depthZone=as.factor(snpPCoA$depthZone)
#snpPCoA$depthZone = factor(snpPCoA$depthZone, levels(snpPCoA$depthZone)[c(2,1)])

# SNP PCoA biplot
mxSnpPcoaPlotA = ggplot(snpPCoA, aes(x = PCo1, y = PCo2, color = pop, fill = pop, shape = depthZone, linetype = depthZone)) +
  geom_hline(yintercept = 0, color = "gray90", size = 0.5) +
  geom_vline(xintercept = 0, color = "gray90", size = 0.5) +
  #stat_ellipse(data = subset(snpPCoA, type = "n", geom = "polygon")) + #ellipse
  #scale_linetype_manual(values=c(1,2,4,3), name = "Depth Zone")+
  geom_point(aes(x = PCo1, y = PCo2, shape = depthZone), size = 3, alpha = 0.3, show.legend = FALSE, guide=FALSE) + #individual's points indicated by circles
  scale_shape_manual(values = c(24,22,23,25), breaks=c("10","15", "25", "35"), name = "Depth Zone") +
  geom_point(aes(x = mean.x, y = mean.y, shape = depthZone), size = 5, color = "black") + #population centroids indicated by triangles
  #scale_fill_paletteer_d("LaCroixColoR::PassionFruit", breaks=c("Alacranes", "Bajos del Norte"), name= "Population")+
  scale_fill_manual(values = c("seagreen3", "blue"), name= "Population")+
  #scale_color_paletteer_d("LaCroixColoR::PassionFruit", breaks=c("Alacranes", "Bajos del Norte"), name= "Population", guide=FALSE) +
  scale_color_manual(values = c("seagreen3", "blue"), name= "Population")+
  xlab(paste ("PCo 1 (", mxSnpPcoaVar[1],"%)", sep = "")) + #Prints percent variation explained by first axis
  ylab(paste ("PCo 2 (", mxSnpPcoaVar[2],"%)", sep = "")) + #Prints percent variation explained by second axis
  guides(shape = guide_legend(order = 2), fill = guide_legend(override.aes = list(shape = 22, size = 4, color = NA), order = 1))+ #linetype = guide_legend(order = 3),
  theme_bw()

mxSnpPcoaPlot = mxSnpPcoaPlotA +
  theme(axis.title.x = element_text(color = "black", size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "left",
        panel.border = element_rect(color = "black", size = 1.2),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



mxSnpPcoaPlot

ggsave("mxSnpPcoaPlot.tiff", plot = mxSnpPcoaPlot, height = 5, width = 7, units = "in", dpi = 300)

#################################################################################################
colPal2=c("#FFFFCC", "#41B6C4")
mxK2 <- read.csv("ngsAdmixK2.csv")
mxK2$sample = factor(mxK2$sample, levels = mxK2$sample[order(-mxK2$cluster2)])
mdat2 = melt(mxK2, id.vars=c("sample", "popDepth"), variable.name="Ancestry", value.name="Fraction")
mdat2$pop=as.factor(mdat2$pop)

p2 = ggplot(mdat2, aes(x=sample, y=Fraction, fill=Ancestry, order=sample)) +
  geom_bar(stat="identity", position="fill", width=1, colour="grey25") +
  facet_grid(.~pop, scales = "free", switch = "x", space = "free") +
  labs(x = "Population", y = "Ancestry") +
  ggtitle("K=2") +
  theme(plot.title = element_text(size=22),
        panel.grid=element_blank(),
        panel.background=element_rect(fill=NA, colour="grey25"),
        panel.spacing.x=grid:::unit(0, "lines"),
        panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        axis.text.x=element_blank(),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x=element_blank(),
        strip.background=element_blank(),
        strip.text=element_blank(),
        #strip.text=element_text(size=20, angle=90),
        legend.key=element_blank(),
        legend.position = "none",
        legend.title = element_blank()) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  #scale_fill_manual(values = brewer.pal(name = "YlGnBu", n = 4), name = "Cluster") +
  scale_fill_manual(values = colPal2, name = "Cluster") +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))
p2

#colPal3=c("#FFFFCC", "#41B6C4", "#225EA8")
colPal3=c("#41B6C4", "#225EA8", "#FFFFCC")
mxK3 <- read.csv("ngsAdmixK3.csv")
mxK3$sample = factor(mxK3$sample, levels = mxK3$sample[order(-mxK3$cluster2)])
mdat3 = melt(mxK3, id.vars=c("sample", "popDepth"), variable.name="Ancestry", value.name="Fraction")
mdat3$pop=as.factor(mdat3$pop)

p3 = ggplot(mdat3, aes(x=sample, y=Fraction, fill=Ancestry, order=sample)) +
  geom_bar(stat="identity", position="fill", width=1, colour="grey25") +
  facet_grid(.~pop, scales = "free", switch = "x", space = "free") +
  labs(x = "Population", y = "Ancestry") +
  ggtitle("K=3") +
  theme(plot.title = element_text(size=22),
        panel.grid=element_blank(),
        panel.background=element_rect(fill=NA, colour="grey25"),
        panel.spacing.x=grid:::unit(0, "lines"),
        panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        axis.text.x=element_blank(),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x=element_blank(),
        strip.background=element_blank(),
        strip.text=element_text(size=20, angle=90),
        legend.key=element_blank(),
        legend.position = "none",
        legend.title = element_blank()) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  #scale_fill_manual(values = brewer.pal(name = "YlGnBu", n = 4), name = "Cluster") +
  scale_fill_manual(values = colPal3, name = "Cluster") +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))
p3

combinedAdmix = (p2 / p3)

ggsave("combinedAdmix.tiff", plot = combinedAdmix, width = 50, height = 30, units = "cm", dpi = 300)

###################################################################################
mxVcf = read.vcfR("mxMcavNoClones.vcf")
mxGenlightPopDepth = vcfR2genlight(mxVcf, n.cores = 2) # Converts the vcf file into a file format that poppr uses the "genlight" format
locNames(mxGenlightPopDepth) = paste(mxVcf@fix[,1],mxVcf@fix[,2],sep="_")
popData = read.csv("mxInds2PopsNoClones.csv") # Reads in population data for each sample

strata(mxGenlightPopDepth) = data.frame(popData)
setPop(mxGenlightPopDepth) = ~popDepth
amova <- poppr.amova(mxGenlightPopDepth, ~popDepth) #Runs AMOVA
amova
set.seed(1999)
amovasignif <- randtest(amova, nrepet = 99) #Calculates significance levels of the AMOVA with 99 permutations
amovasignif
plot(amovasignif)

#mxGenlightPopDepth$pop=as.factor(mxGenlightPopDepth$pop)
#mxGenlightPopDepth$pop = factor(mxGenlightPopDepth$pop, levels(mxGenlightPopDepth$pop)[c(4,3,2,1,5,6,8,7)])

set.seed(694)
mx.fst <- stamppFst(mxGenlightPopDepth, nboots = 99, percent = 95, nclusters = 4) #99 permutations
mx.fst$Fsts
mx.fst$Pvalues

pca.1 <- glPca(mxGenlightPopDepth, nf=300, n.cores=1) 

# proportion of explained variance by first three axes
pca.1$eig[1]/sum(pca.1$eig) # proportion of variation explained by 1st axis
pca.1$eig[2]/sum(pca.1$eig) # proportion of variation explained by 2nd axis 
pca.1$eig[3]/sum(pca.1$eig) # proportion of variation explained by 3rd axis

#####K-Means Clustering and DAPC
grp <- find.clusters(mxGenlightPopDepth, max.n.clust=11, glPca = pca.1, perc.pca = 100, n.iter=1e6, n.start=1000) 


pop.order <- c("Alacranes-10", "Bajos del Norte-10", "Alacranes-15", "Bajos del Norte-15", "Alacranes-25", "Bajos del Norte-25","Alacranes-35", "Bajos del Norte-35")

# reads in fst matrix
snpFstMa <- as.matrix(mx.fst$Fsts)
upperTriangle(snpFstMa, byrow=TRUE) <- lowerTriangle(snpFstMa)
snpFstMa <- snpFstMa[,pop.order] %>%
  .[pop.order,]
snpFstMa[upper.tri(snpFstMa)] <- NA
snpFstMa <- as.data.frame(snpFstMa)

snpFstMa$Pop = factor(row.names(snpFstMa), levels = unique(pop.order))

snpQMa <- as.matrix(mx.fst$Pvalues)
upperTriangle(snpQMa, byrow=TRUE) <- lowerTriangle(snpQMa)
snpQMa <- snpQMa[,pop.order] %>%
  .[pop.order,]
snpQMa[upper.tri(snpQMa)] <- NA
snpQMa <- as.data.frame(snpQMa)
snpQMa$Pop = factor(row.names(snpQMa), levels = unique(pop.order))

snpFstMa$Pop = factor(row.names(snpFstMa), levels = unique(pop.order))
snpFst = melt(snpFstMa, id.vars = "Pop", value.name = "Fst", variable.name = "Pop2", na.rm = TRUE)
snpFst = snpFst[snpFst$Pop != snpFst$Pop2,]
snpFst$Fst = round(snpFst$Fst, 3)
snpFst = snpFst %>% mutate(Fst = replace(Fst, Fst < 0, 0))
head(snpFst)

snpQ = melt(snpQMa, id.vars = "Pop", value.name = "Pval", variable.name = "Pop2", na.rm = TRUE)
snpQ = snpQ[snpQ$Pop != snpQ$Pop2,]
snpQ$Qval = p.adjust(snpQ$Pval, method = "BH")
head(snpQ)

snpHeatmapA = ggplot(data = snpFst, aes(Pop, Pop2, fill = Fst))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "red", midpoint = 0, limit = c(0, 0.08),
                       space = "Lab", name = expression(paste(italic("F")[ST])))+
  geom_text(data = snpFst, aes(Pop, Pop2, label = Fst), color = ifelse(snpQ$Qval <= 0.05,"black", "darkgrey"), size = ifelse(snpQ$Qval < 0.05, 6, 5)) +
  guides(fill=guide_colorbar(barwidth = 1, barheight = 12, title.position = "top", title.hjust = 0.5))+                     
  scale_y_discrete(position = "right")+
  scale_x_discrete(labels = str_wrap(c("Bajos del Norte-10", "Alacranes-15", "Bajos del Norte-15", "Alacranes-25", "Bajos del Norte-25","Alacranes-35", "Bajos del Norte-35"), width = 6)) +
  #ggtitle("   SNP") +
  theme_minimal()

snpHeatmap = snpHeatmapA + theme(
  axis.text.x = element_text(vjust = 1, size = 16, hjust = 0.5, color = "black"),
  axis.text.y = element_text(size = 16, color = "black"),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.ticks = element_blank(),
  legend.position = "right",
  legend.direction = "vertical",
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 14),
  plot.title = element_text(size = 16)
)

snpHeatmap

ggsave("snpHeatMap.tiff", plot = snpHeatmap, width = 30, height = 15, units = "cm", dpi = 300)


#####Bayescan
source('plot_R.r')

#NewVersion
dat1 = read.table("mxMcav.baye_fst.txt",header=T)
head(dat1)
table(dat1[,"qval"]<0.1)
dat1$locus = c(1:nrow(dat1))
outs=which(dat1[,"qval"]<0.1)
plot_bayescan("mxMcav.baye_fst.txt",FDR=0.1,add_text=F,size=0.5,highlight=outs)

outliers= dat1 %>% filter(qval < 0.1)
outsQval=select(outliers, qval, locus)

#OldVersion
#dat2 = read.table("mxMcavInds.baye_fst.txt",header=T)
#head(dat2)
#table(dat2[,"qval"]<0.1)
#dat2$locus = c(1:nrow(dat2))
#outs=which(dat2[,"qval"]<0.1)
#plot_bayescan("mxMcavInds.baye_fst.txt",FDR=0.1,add_text=F,size=0.5,highlight=outs)


#Minor Allele Frequencies of Outliers Across Populations
# outAlleleFreq=read.csv("outlierAlleleFreqs.csv")
# outAlleleFreqPlot <- ggplot(outAlleleFreq, aes(x=PopDepth, y=MIF)) + 
#   geom_boxplot()+
#   theme_bw()
# outAlleleFreqPlot
# 
# outAlleleFreq.lm = lm(MIF ~ AvgDepth, data=outAlleleFreq)
# summary(outAlleleFreq.lm)$r.squared 
# 
# ggplot(outAlleleFreq, aes(x=AvgDepth, y=MIF)) + geom_point()+ geom_smooth(method=lm)

## SNPs potentially under selection
genes = read.table("mcav_gene_regions.tab")
names(genes) = c("chromo","start","end","gene")

# expand gene regions ± 2000 bp
genes$start = genes$start -2000
genes$end = genes$end +2000

gnames = read.table("mcav_cog.txt", sep = "\t", quote="", fill=FALSE)
names(gnames) = c("gene", "cog", "protein")

genes = full_join(genes, gnames, by = "gene")
genes$protein=as.character(genes$protein)
genes$protein[is.na(genes$protein)]="unknown"

#how many annotated genes do we have?
nrow(genes[genes$protein!="unknown",])
head(genes)

snpLoci = read.table("mxMcavNoClones.mafs.gz", header = TRUE)
snpLoci$locus = c(1:nrow(snpLoci))

outsByGene = snpLoci %>% dplyr::select(locus, chromo, position) %>% 
  filter(locus %in% outs) %>%
  #filter(locus %in% outsQval) %>%
  full_join(., outsQval, by = "locus") %>%
  full_join(., genes, by = "chromo") %>% 
  filter(position>start, position<end)

outsByGene = outsByGene[c(1:3,5:9,4)]

for(i in 1:nrow(outsByGene)) {
  if (outsByGene$qval[i] < .01) {
    outsByGene$sig[i] = "***"
  } else {
    if (outsByGene$qval[i] < .05) {
      outsByGene$sig[i] = "**"
    } else {
      outsByGene$sig[i] = "*"
    }
  }
}

outsByGene

#only show annotated genes
outsAnno = outsByGene %>% filter(protein != "unknown")
outsAnno

# df_uniq <- unique(outsAnno$locus)
# length(df_uniq)
# df_uniq1 <- unique(outsAnno$gene)
# length(df_uniq1)

write.csv(x=outsAnno, file="annotatedOutliers.csv")
###################################################################################
#Minor allele frequencies and depth?
mxMIF=read.csv("mxMIF.csv")
mxMIFPlot <- ggplot(mxMIF, aes(group=pop, x=depth, y=mif, color=site)) + 
geom_boxplot()+
theme_bw()
mxMIFPlot

mxMIF.lm = lm(mif ~ depth*site, data=mxMIF)
summary(mxMIF.lm)$r.squared 

as.factor(mxMIF$locus)
ggplot(mxMIF, aes(group=pop, x=depth, y=mif, color=site, shape=locus)) + geom_point()+ geom_smooth(method=lm)+scale_shape_identity()
###################################################################################
#####Bayescenve
source('plot_R.r')

#NewVersion
datDepth = read.table("mxMcavDepth.baye_fst.txt",header=T)
head(datDepth)
table(datDepth[,"qval_alpha"]<0.1)
datDepth$locus = c(1:nrow(datDepth))
outsDepth=which(datDepth[,"qval"]<0.1)
plot_bayescan("mxMcavDepth.baye_fst.txt",FDR=0.1,add_text=F,size=0.5,highlight=outs)

outliersDepth= datDepth %>% filter(qval < 0.1)
outsQvalDepth=select(outliersDepth, qval, locus)

genes = read.table("mcav_gene_regions.tab")
names(genes) = c("chromo","start","end","gene")

# expand gene regions ± 2000 bp
genes$start = genes$start -2000
genes$end = genes$end +2000

gnames = read.table("mcav_cog.txt", sep = "\t", quote="", fill=FALSE)
names(gnames) = c("gene", "cog", "protein")

genes = full_join(genes, gnames, by = "gene")
genes$protein=as.character(genes$protein)
genes$protein[is.na(genes$protein)]="unknown"

#how many annotated genes do we have?
nrow(genes[genes$protein!="unknown",])
head(genes)

snpLoci = read.table("mxMcavNoClones.mafs.gz", header = TRUE)
snpLoci$locus = c(1:nrow(snpLoci))

outsByGeneDepth = snpLoci %>% dplyr::select(locus, chromo, position) %>% 
  filter(locus %in% outsDepth) %>%
  #filter(locus %in% outsQval) %>%
  full_join(., outsQvalDepth, by = "locus") %>%
  full_join(., genes, by = "chromo") %>% 
  filter(position>start, position<end)

outsByGeneDepth = outsByGeneDepth[c(1:3,5:9,4)]

for(i in 1:nrow(outsByGeneDepth)) {
  if (outsByGeneDepth$qval[i] < .01) {
    outsByGeneDepth$sig[i] = "***"
  } else {
    if (outsByGeneDepth$qval[i] < .05) {
      outsByGeneDepth$sig[i] = "**"
    } else {
      outsByGeneDepth$sig[i] = "*"
    }
  }
}

outsByGeneDepth

#only show annotated genes
outsAnnoDepth = outsByGeneDepth %>% filter(protein != "unknown")
outsAnnoDepth

# df_uniq <- unique(outsAnno$locus)
# length(df_uniq)
# df_uniq1 <- unique(outsAnno$gene)
# length(df_uniq1)

write.csv(x=outsAnnoDepth, file="annotatedOutliersDepth.csv")
###################################################################################
zoox=read.csv("zooxAlignmentProp.csv", header=T)
set.seed(694)
zooxAdonis = adonis(zoox ~ depthZone*pop, data = zoox, permutations = 9999, method = "bray")
zooxAdonis

set.seed(694)
zooxPWAdonis = pairwise.adonis(zoox,
                                   factors = zoox$depthZone,
                                   sim.method = "bray", p.adjust.m = "BH", perm = 9999)

zooxPWAdonis

#make community matrix - extract columns with abundance information
com = zoox[,4:ncol(zoox)]
m_com = as.matrix(com)
set.seed(123)
nmds = metaMDS(m_com, distance = "bray")
nmds
plot(nmds)
data.scores = as.data.frame(scores(nmds))
data.scores$sample = zoox$sample
data.scores$pop = zoox$pop
data.scores$depthZone = zoox$depthZone
head(data.scores)

xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes( shape = as.factor(depthZone), colour = as.factor(pop)))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Population", y = "NMDS2", shape = "Depth Zone")  + 
  scale_colour_manual(values = c("#009E73", "#E69F00"))+
  scale_shape_manual(values = c(24,22,23,25))

xx

################################################################################################
#####################################################
snpZooxMa = as.matrix(read.table("mxZooxNoClones87Removed.ibsMat"))

mxZooxMds = cmdscale(snpZooxMa, eig = TRUE, x.ret = TRUE, k=2)

# Determine percent variation captured on each axis
# Calculate the eigenvalues so later we can figure out % variation shown on each Principal Coordinate
mxZooxSnpPcoaVar = round(mxZooxMds$eig/sum(mxZooxMds$eig)*100, 1)
mxZooxSnpPcoaVar

# Format data to plot
mxZooxSnpPcoaValues = mxZooxMds$points
mxZooxSnpPcoaValues

snpI2P = read.csv("mxZooxInds2PopsNoClones.csv") # 2-column tab-delimited table of individual assignments to populations; must be in the same order as samples in the bam list or vcf file.
row.names(snpI2P) = snpI2P[,1]
mxZooxSnpPcoaValues=cbind(snpI2P, mxZooxSnpPcoaValues)
mxZooxSnpPcoaValues =as.data.frame(mxZooxSnpPcoaValues, sample = rownames(mxZooxSnpPcoaValues))
colnames(mxZooxSnpPcoaValues)[5] <- "PCo1"
colnames(mxZooxSnpPcoaValues)[6] <- "PCo2"
mxZooxSnpPcoaValues

snpPCoA = merge(mxZooxSnpPcoaValues, aggregate(cbind(mean.x=PCo1,mean.y=PCo2)~popDepth, mxZooxSnpPcoaValues, mean), by="popDepth")

snpPCoA$depthZone=as.factor(snpPCoA$depthZone)
#snpPCoA$depthZone = factor(snpPCoA$depthZone, levels(snpPCoA$depthZone)[c(2,1)])

# SNP PCoA biplot
mxZooxSnpPcoaPlotA = ggplot(snpPCoA, aes(x = PCo1, y = PCo2, color = pop, fill = pop, shape = depthZone, linetype = depthZone)) +
  geom_hline(yintercept = 0, color = "gray90", size = 0.5) +
  geom_vline(xintercept = 0, color = "gray90", size = 0.5) +
  #stat_ellipse(data = subset(snpPCoA, type = "n", geom = "polygon")) + #ellipse
  #scale_linetype_manual(values=c(1,2,4,3), name = "Depth Zone")+
  geom_point(aes(x = PCo1, y = PCo2, shape = depthZone), size = 3, alpha = 0.3, show.legend = FALSE, guide=FALSE) + #individual's points indicated by circles
  scale_shape_manual(values = c(24,22,23,25), breaks=c("10","15", "25", "35"), name = "Depth Zone") +
  geom_point(aes(x = mean.x, y = mean.y, shape = depthZone), size = 5, color = "black") + #population centroids indicated by triangles
  #scale_fill_paletteer_d("LaCroixColoR::PassionFruit", breaks=c("Alacranes", "Bajos del Norte"), name= "Population")+
  scale_fill_manual(values = c("seagreen3", "blue"), name= "Population")+
  #scale_color_paletteer_d("LaCroixColoR::PassionFruit", breaks=c("Alacranes", "Bajos del Norte"), name= "Population", guide=FALSE) +
  scale_color_manual(values = c("seagreen3", "blue"), name= "Population")+
  xlab(paste ("PCo 1 (", mxZooxSnpPcoaVar[1],"%)", sep = "")) + #Prints percent variation explained by first axis
  ylab(paste ("PCo 2 (", mxZooxSnpPcoaVar[2],"%)", sep = "")) + #Prints percent variation explained by second axis
  guides(shape = guide_legend(order = 2), fill = guide_legend(override.aes = list(shape = 22, size = 4, color = NA), order = 1))+ #linetype = guide_legend(order = 3),
  theme_bw()

mxZooxSnpPcoaPlot = mxZooxSnpPcoaPlotA +
  theme(axis.title.x = element_text(color = "black", size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "left",
        panel.border = element_rect(color = "black", size = 1.2),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



mxZooxSnpPcoaPlot

ggsave("mxZooxSnpPcoaPlot.tiff", plot = mxZooxSnpPcoaPlot, height = 5, width = 7, units = "in", dpi = 300)
