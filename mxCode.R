setwd("/Users/student/Documents/GitHub/mxMcavSnp")

if (!require("pacman")) install.packages("pacman")
pacman::p_load("adegenet", "dendextend", "flextable", "gdata", "ggdendro", "hierfstat", "Imap", "kableExtra", 
               "paletteer", "patchwork", "officer", "poppr", "RColorBrewer", "reshape2", "StAMPP", 
               "tidyverse", "vcfR", "vegan", "WGCNA", "ggnewscale", "rnaturalearth", "sf", "ggspatial", 
               "cowplot", "marmap", "dplyr", "edgeR", "MCMC.OTU", "pairwiseAdonis", "distGeo")


pacman::p_load_gh("eliocamp/ggnewscale", "ropensci/rnaturalearthhires")

mxSites = read.csv("mxSampleSites.csv", header = TRUE)
mxSites$depthZone=as.factor(mxSites$depthZone)
# coerce to sf object

mxSites1 <- st_as_sf(x = mxSites, coords = c('lon', 'lat'), crs = 4326)
mxSitesDist=st_distance(mxSites1)
mxSitesDistKm=mxSitesDist/1000
rownames(mxSitesDistKm) <-c(mxSites$reef)
colnames(mxSitesDistKm) <-c(mxSites$reef)
#mxSites$depth = factor(mxSites$depth, levels = levels(mxSites$depth)[c(2,1)])
#mxSites$reef = factor(mxSites$reef, levels = levels(mxSites$reef)[c(1, 7, 8, 3, 6, 2, 5, 4, 11, 10, 9, 12)])
mxReefs = mxSites %>% group_by(reef) %>% summarize(lat = mean(lat), lon = mean(lon))
#mxReefs$depthShape=c("Mesophotic", "Shallow", "Shallow", "Shallow", "Shallow", "Shallow", "Shallow", "Shallow", "Paired", "Paired", "Paired", "Paired")
#mxReefs$depthShape=as.factor(mxReefs$depthShape)
#mxReefs$depthShape = factor(mxReefs$depthShape, levels = levels(mxReefs$depthShape)[c(2,3,1)])
gomCountries = st_as_sf(ne_countries(country = c("United States of America", "Cuba", "Mexico", "The Bahamas", "Belize", "Guatemala", "Cayman Islands"), scale = "Large"))

# gomBathy=st_read("/Users/student/Documents/GitHub/mxMcavSnp/GOM_isobath_100-1000m/w98e78n31s18_isobath_100-1000m.shp") %>% st_transform(54030)
gom5.100Bathy=st_read("/Users/student/Documents/GitHub/mxMcavSnp/gom_5-100m/w98e78n31s18_isobath_5-100m.shp") %>% st_transform(54030)
gom200Bathy=st_read("/Users/student/Documents/GitHub/mxMcavSnp/GOM_isobath_100-1000m/w98e78n31s18_isobath_100-1000m.shp") %>%
filter(CONTOUR == "-200") %>% st_transform(54030)
# gom150Bathy=st_read("/Users/student/Documents/GitHub/mxMcavSnp/gom_5-4000m/w98e78n31s18_isobath_selected_5-4000m.shp") %>%
#   filter(CONTOUR == "-50") %>% st_transform(54030)
#gom150Bathy=st_read("/Users/student/Documents/GitHub/mxMcavSnp/10mDepthContours/depth_contours_10m_intervals_GEBCO_2020_Gulf_of_Mexico_bathymetry.shp") %>% filter(contour == "-150") %>% st_transform(54030)
#gomBathy=st_read("/Users/student/Documents/GitHub/mxMcavSnp/gebco_2019_contours/gebco_2019_contours.shp")
#gomBathy= gomBathy %>% filter(DEPTH == "-150")

mxPal=c("seagreen3", "blue")

mxMapA = ggplot() +
  geom_sf(data = gomCountries, fill = "white", size = 0.25) +
  geom_sf(data = gom200Bathy, color = "red") +
  geom_point(data = mxSites, aes(x = lon, y = lat, shape = depthZone), size = 0) +
  geom_point(data = mxReefs, aes(x = lon, y = lat, fill = reef), size = 3, shape=21) +
  scale_fill_manual(values = mxPal, name = "Site") +
  scale_shape_manual(values = c(24,22,23,25), breaks=c("10","15", "25", "35"), name = "Depth Zone", labels=c("10 m","15 m", "25 m", "35 m")) +
  #scale_color_manual(values = c("blue"), labels = c("MPA Boundaries")) + #name = "Boundaries",
  coord_sf(xlim = c(-103, -80), ylim = c(18.5, 30)) +
  scale_x_continuous(breaks = c(seq(-103, -80, by = 2))) +
  scale_y_continuous(breaks = c(seq(18.5, 30, by = 2))) +
  #geom_rect(aes(xmin = -89.8, xmax = -89.6, ymin = 22.35, ymax = 22.6), color = "black", fill = NA, size = 0.4) +
  #geom_rect(aes(xmin = -88.74, xmax = -88.68, ymin = 23.24, ymax = 23.32), color = "black", fill = NA, size = 0.4) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "tr", which_north = "true", style = north_arrow_minimal()) +
  guides(fill = guide_legend(override.aes = list(shape = 22, color = NA, size = 4), ncol = 2, order = 1), shape = guide_legend(override.aes = list(size = 3), order = 2)) +
  #guides(fill = guide_legend(override.aes = list(shape = 21, size = 4), order = 1))+ #,
         #shape = guide_legend(order=2),
         #color = guide_legend(title=NULL, order = 3))+ #override.aes = list(fill=NA),
  theme_bw() +
  theme(plot.title = element_text(size = 9) ,
        panel.background = element_rect(fill = "aliceblue"),
        plot.background = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom")

mxMapA

#ggsave("mxMapNoInset.tiff", plot = mxMapA, width = 30, height = 20, units = "cm", dpi = 600)
set.seed(224)
inset = ggplot() +
  geom_segment(aes(x = -103, xend = -80, y = 18.5, yend = 30), color = "gray92", size = .55) +
  geom_sf(data = gomCountries, fill = "white", size = 0.25) +
  geom_point(data = subset(mxSites, subset = mxSites$reef == "Alacranes"), aes(x = lon, y = lat, fill = reef, shape=depthZone), size = 3, position=position_jitter(width=0.01, height=0.01)) +
  geom_point(data = subset(mxSites, subset = mxSites$reef == "Bajos del Norte"), aes(x = lon, y = lat, fill = reef, shape=depthZone), size = 3, position=position_jitter(width=0.002, height=0.002)) +
  scale_fill_manual(values = mxPal, name = "Reef") +
  scale_shape_manual(values = c(24,22,23,25), breaks=c("10","15", "25", "35"), name = "Depth Zone", labels=c("10 m","15 m", "25 m", "35 m")) +
  #geom_jitter()+
  #geom_sf(data = gom5.100Bathy, color = "gray") +
  guides(fill = guide_legend(override.aes = list(shape = 21, color = NA, size = 4))) +
  theme_bw() +
  theme(legend.title = element_text(size = 9, face = "bold"),
        axis.title = element_blank(),
        panel.background = element_rect(fill = "aliceblue"),
        legend.text = element_text(size = 9),
        legend.position = "none",
        legend.direction = "vertical",
        legend.box = "vertical",
        panel.grid=element_blank(),
        rect = element_blank())

# inset

alacranes = inset+
  ggtitle("Alacranes") +
  # geom_sf(data = gomCountries, fill = "white", size = 0.25) +
  annotation_scale(location = "bl") +
  coord_sf(xlim = c(-89.8, -89.6), ylim = c(22.35, 22.6)) +
  scale_x_continuous(breaks = c(seq(-89.8, -89.6, by = .1))) +
  scale_y_continuous(breaks = c(seq(22.35, 22.6, by = .1)))+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background=element_blank(),
        panel.border = element_rect(size = 1))

# alacranes

bajosDelNorte = inset +
  ggtitle("Bajos del Norte") +
  annotation_scale(location = "bl") +
  coord_sf(xlim = c(-88.74, -88.68), ylim = c(23.24, 23.32)) +
  scale_x_continuous(breaks = c(seq(-88.74, -88.68, by = .1))) +
  scale_y_continuous(breaks = c(seq(23.24, 23.32, by = .1)))+
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background=element_blank(),
        plot.title = element_text(vjust=-6.5, hjust=0.5),
        panel.border = element_rect(size = 1))

#bajosDelNorte

mxMap = ggdraw() +
  draw_plot(mxMapA) +
  #draw_plot(inset) +
  draw_plot(alacranes, x = 0.01, y = 0.25, width = 0.3, height = 0.3) +
  draw_plot(bajosDelNorte, x = 0.67, y = 0.25, width = 0.3, height = 0.3)


#mxMap

ggsave("mxMap.tiff", plot = mxMap, width = 25, height = 20, units = "cm", dpi = 600)

#############################################################################################

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

snpPCoA$depthZoneM=as.factor(snpPCoA$depthZoneM)
#snpPCoA$depthZone = factor(snpPCoA$depthZone, levels(snpPCoA$depthZone)[c(2,1)])

# SNP PCoA biplot
mxSnpPcoaPlotA = ggplot(snpPCoA, aes(x = PCo1, y = PCo2, color = site, fill = site, shape = depthZoneM, linetype = depthZoneM)) +
  geom_hline(yintercept = 0, color = "gray90", size = 0.5) +
  geom_vline(xintercept = 0, color = "gray90", size = 0.5) +
  #stat_ellipse(data = subset(snpPCoA, type = "n", geom = "polygon")) + #ellipse
  #scale_linetype_manual(values=c(1,2,4,3), name = "Depth Zone")+
  geom_point(aes(x = PCo1, y = PCo2, shape = depthZoneM), size = 3, alpha = 0.3, show.legend = FALSE, guide=FALSE) + #individual's points indicated by circles
  scale_shape_manual(values = c(24,22,23,25), breaks=c("10","15", "25", "35"), name = "Depth Zone", labels=c("10 m","15 m", "25 m", "35 m")) +
  geom_point(aes(x = mean.x, y = mean.y, shape = depthZoneM), size = 5, color = "black") + #population centroids indicated by triangles
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
paletteer_d("NineteenEightyR::miami2")
#colPal2=c("#FFFFCC", "#41B6C4")
colPal2=c("#3FB8AFFF", "#DAD8A7FF")
mxNgsAdmix <- read.csv("mxMcavNoClones_k2k3.csv") %>% arrange(site, depthZoneM, -cluster2.2)
# mxK2$sample = factor(mxK2$sample, levels = mxK2$sample[order(-mxK2$cluster2)])
sampleCounts = plyr::count(mxNgsAdmix, c('site','depthZoneM'))
meltedList = reshape2::melt(lapply(sampleCounts$freq,function(x){c(1:x)}))

mxNgsAdmix$barPlotOrder = meltedList$value
mxNgsAdmix$depthZoneM=as.factor(mxNgsAdmix$depthZoneM)
mdat = melt(mxNgsAdmix, id.vars=c("sample", "site", "depthZoneM", "barPlotOrder"), variable.name="Ancestry", value.name="Fraction")
#mdat2$site=as.factor(mdat2$site)
mdat$depthZoneM=as.factor(mdat$depthZoneM)
levels(mdat$depthZoneM)
levels(mdat$depthZoneM) = c("10 m", "15 m", "25 m", "35 m")

p2 = ggplot(data = subset(mdat, subset = mdat$Ancestry %in% c("cluster2.1", "cluster2.2")), aes(x=barPlotOrder, y=Fraction, fill=Ancestry, order=barPlotOrder)) +
  geom_bar(stat="identity", position="fill", width=1, colour="grey25") +
  #facet_grid(.~site, scales = "free", switch = "x", space = "free") +
  facet_grid(depthZoneM ~ site, scales = "free") + #faceting plots by Depth and Site
  labs(x = "Population", y = "Ancestry") +
  ggtitle(expression(paste(italic("K"), "=2"))) +
  theme(plot.title = element_text(size=30),
        panel.grid=element_blank(),
        panel.background=element_rect(fill=NA, colour="grey25"),
        panel.spacing.x=grid:::unit(0, "lines"),
        panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        axis.text.x=element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        strip.background = element_rect(fill = "white", size = 0.9, color="black"),
        strip.text.x = element_text(size = 30),
        strip.text.y = element_blank(),
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

###########
#colPal3=c("#FFFFCC", "#41B6C4", "#225EA8")
colPal3=c("#3FB8AFFF", "#DAD8A7FF", "#FF3D7FFF")
# mxK3 <- read.csv("mxMcavNoClones_K3.csv") %>% arrange(pop, depthZoneM, -cluster2)
# # mxK3$sample = factor(mxK3$sample, levels = mxK3$sample[order(-mxK3$cluster2)])
# sampleCountsK3 = plyr::count(mxK3, c('pop','depthZoneM'))
# meltedListK3 = reshape2::melt(lapply(sampleCountsK3$freq,function(x){c(1:x)}))
# 
# mxK3$barPlotOrderK3 = meltedListK3$value
# mdat3 = melt(mxK3, id.vars=c("sample", "pop", "depthZoneM", "barPlotOrderK3"), variable.name="Ancestry", value.name="Fraction")
# #mdat2$pop=as.factor(mdat3$pop)
# mdat3$depthZoneM=as.factor(mdat3$depthZoneM)
# levels(mdat3$depthZoneM)
# levels(mdat3$depthZoneM) = c("10 m", "15 m", "25 m", "35 m")


p3 = ggplot(data = subset(mdat, subset = mdat$Ancestry %in% c("cluster3.1", "cluster3.2", "cluster3.3")), aes(x=barPlotOrder, y=Fraction, fill=Ancestry, order=barPlotOrder)) +
  geom_bar(stat="identity", position="fill", width=1, colour="grey25") +
  #facet_grid(.~pop, scales = "free", switch = "x", space = "free") +
  facet_grid(depthZoneM ~ site, scales = "free") + #faceting plots by Depth and Site
  labs(x = "Population", y = "Ancestry") +
  ggtitle(expression(paste(italic("K"), "=3"))) +
  theme(plot.title = element_text(size=30),
        panel.grid=element_blank(),
        panel.background=element_rect(fill=NA, colour="grey25"),
        panel.spacing.x=grid:::unit(0, "lines"),
        panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        axis.text.x=element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        strip.background = element_rect(fill = "white", size = 0.9, color="black"),
        strip.text.x = element_text(size = 30),
        strip.text.y = element_text(size = 30),
        #strip.text=element_text(size=20, angle=90),
        legend.key=element_blank(),
        legend.position = "none",
        legend.title = element_blank()) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  #scale_fill_manual(values = brewer.pal(name = "YlGnBu", n = 4), name = "Cluster") +
  scale_fill_manual(values = colPal3, name = "Cluster") +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))
p3

combinedAdmix = (p2 | p3)
combinedAdmix

ggsave("combinedAdmix.tiff", plot = combinedAdmix, width = 70, height = 40, units = "cm", dpi = 300)

###################################################################################
mxVcf = read.vcfR("mxMcavNoClones.bcf")
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

# pca.1 <- glPca(mxGenlightPopDepth, nf=300, n.cores=1) 
# 
# # proportion of explained variance by first three axes
# pca.1$eig[1]/sum(pca.1$eig) # proportion of variation explained by 1st axis
# pca.1$eig[2]/sum(pca.1$eig) # proportion of variation explained by 2nd axis 
# pca.1$eig[3]/sum(pca.1$eig) # proportion of variation explained by 3rd axis
# 
# #####K-Means Clustering and DAPC
# grp <- find.clusters(mxGenlightPopDepth, max.n.clust=11, glPca = pca.1, perc.pca = 100, n.iter=1e6, n.start=1000) 
# 

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
  geom_text(data = snpFst, aes(Pop, Pop2, label = Fst), color = ifelse(snpQ$Qval <= 0.05,"black", "darkgrey"), size = ifelse(snpQ$Qval < 0.05, 6, 5), fontface = ifelse(snpQ$Qval <= 0.05, "bold", "plain")) +
  guides(fill=guide_colorbar(barwidth = 1, barheight = 12, title.position = "top", title.hjust = 0.5))+                     
  scale_y_discrete(position = "right", labels = c("Alacranes-10 m", "Bajos del Norte-10 m","Alacranes-15 m", "Bajos del Norte-15 m", "Alacranes-25 m", "Bajos del Norte-25 m","Alacranes-35 m"))+
  scale_x_discrete(labels = str_wrap(c("Bajos del Norte-10 m", "Alacranes-15 m", "Bajos del Norte-15 m", "Alacranes-25 m", "Bajos del Norte-25 m","Alacranes-35 m", "Bajos del Norte-35 m"), width = 4)) +
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

#####################################################################################
#popStats
mxMcavBreed = read.csv("./ngsRelate.csv")

mxMcavBreed2 = mxMcavBreed %>% group_by(a) %>% select("inbreed" = Fa)
mxMcavBreed3 = mxMcavBreed %>% group_by(b) %>% select("inbreed" = Fb)

mxMcavBreed = bind_rows(mxMcavBreed2, mxMcavBreed3) %>% group_by(a) %>% summarize("inbreed" = mean(inbreed)) 
mxMcavRelate = read.csv("./ngsRelate.csv")
mxMcavRelate2 = mxMcavRelate %>% group_by(a, b) %>% select("relate" = rab)

popData = read.csv("mxInds2PopsNoClones.csv", header=TRUE)
popData$a=c(0:96)
mxMcavRelate2 = mxMcavRelate2 %>% left_join(popData, by = "a") %>% left_join(popData, by = c("b" = "a"), suffix = c(".a", ".b")) 

mxMcavRelate = mxMcavRelate2 %>% filter(popDepth.a == popDepth.b) %>% rename(depthZoneM = depthZoneM.a, site = site.a)

mxMcavHet=read.csv("mxMcavHet.csv")
het = left_join(popData, mxMcavHet, by = "sample") %>% mutate("inbreed" = mxMcavBreed$inbreed)

hetAllSites= ggplot(data=het, aes(x=popDepth, y=observedHetAllSites, color=site)) + 
  geom_boxplot() +
  geom_point(position = "jitter")+
  scale_color_manual(values = c("seagreen3", "blue"))+
  #labs(subtitle="Boxplot with points using geom_point() with jitter")+
  theme_bw(base_size=16)+
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, vjust=0.5))

hetAllSites

hetVariantSites= ggplot(data=het, aes(x=popDepth, y=observedHetVariantSites, color=site)) + 
  geom_boxplot() +
  geom_point(position = "jitter")+
  scale_color_manual(values = c("seagreen3", "blue"))+
  #labs(subtitle="Boxplot with points using geom_point() with jitter")+
  theme_bw(base_size=16)+
  theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, vjust=0.5))

hetVariantSites

hetVariantSitesAnova <- aov(observedHetVariantSites ~ popDepth, data = het)
summary(hetVariantSitesAnova)

#####################################################################################

zoox = read.delim("zooxGenomeAlignmentRate", header = FALSE, check.names = FALSE)
head(zoox)

zoox$V2[is.na(zoox$V2)] <- as.character(zoox$V1[is.na(zoox$V2)])

zoox$V1 = gsub(pattern = "mx.*", "chr", zoox$V1)
zoox$V2 = gsub(".trim.*", "", zoox$V2)

zoox = zoox %>% filter(zoox$V1 != "*")

zooxLst = split(zoox$V2, as.integer(gl(length(zoox$V2), 20, length(zoox$V2))))

zooxMaps = NULL
for(i in zooxLst){
  zooxMaps = rbind(zooxMaps, data.frame(t(i)))
}
colnames(zooxMaps) = c("sample",zoox$V1[c(2:20)])

for(i in c(2:20)){
  zooxMaps[,i] = as.numeric(levels(zooxMaps[,i]))[zooxMaps[,i]]
}

str(zooxMaps)

zooxMaps$Symbiodinium = rowSums(zooxMaps[2:6])
zooxMaps$Breviolum = rowSums(zooxMaps[7:10])
zooxMaps$Cladocopium = rowSums(zooxMaps[11:16])
zooxMaps$Durusdinium = rowSums(zooxMaps[17:20])

zooxMaps = zooxMaps[,c(1, 21:24)]
zooxProp = zooxMaps

zooxProp$sum = apply(zooxProp[, c(2:length(zooxProp[1,]))], 1, function(x) {
  sum(x, na.rm = T)
})

zooxProp = cbind(zooxProp[, c(1)], (zooxProp[, c(2:(ncol(zooxProp)-1))]
                                    / zooxProp$sum))
head(zooxProp)

apply(zooxProp[, c(2:(ncol(zooxProp)))], 1, function(x) {
  sum(x, na.rm = T)
})

zooxProp=as.data.frame(zooxProp)
mean(zooxProp$Cladcopium)
###############################ITS2###################################################
its2Profs = read.csv("mxIts2ProfCounts.csv", header = TRUE, check.names = FALSE)
head(its2Profs)

its2Profs$depthZoneM = factor(its2Profs$depthZoneM, levels = c("10", "15", "25", "35"))
levels(its2Profs$depthZoneM)

its2Profs$site = factor(its2Profs$site, levels(its2Profs$site)[c(1,2)])
its2Profs = its2Profs[order(its2Profs$site, its2Profs$depthZoneM), ]
head(its2Profs)

its2Profs = its2Profs %>% arrange(site, depthZoneM,
                                  `G3b`,
                                  `G3am`,
                                  `C3/C3b`,
                                  `C3-C3fc-C40g-C3b-C3an-C21-C3fd-C3s-C3bb`,
                                  `C3an/C3-C21-C3bb`, 
                                  `C3-C21-C3gu-C3gy-C3gw-C3gx-C3gv-C3b`, 
                                  `C3-C3an-C3de-C21-C3bb-C3fc-C3b`,
                                  `C3-C3de-C3bb-C21ae-C21-C3an-C3s`, 
                                  `C3-C3fc-C21-C3b-C3an-C3fd-C3bb`, 
                                  `C3-C3de-C21-C3bb-C3an-C21ae-C65b-C3s`, 
                                  `C3-C21-C3an-C3fc-C3gz-C3b-C3bb`, 
                                  `C3-C21-C3fc-C3an-C65b-C3b-C3bb`)  
                      
                      

sampleCounts = plyr::count(its2Profs, c('site','depthZoneM'))
meltedList = reshape2::melt(lapply(sampleCounts$freq,function(x){c(1:x)}))
its2Profs$barPlotOrder = meltedList$value
its2Profs=its2Profs[c(1,ncol(its2Profs),2:(ncol(its2Profs)-1))]

##########TMM Normalized
its2ProfsTransposed = t(its2Profs[, 10:length(its2Profs[1, ])])
its2ProfsList = DGEList(counts = its2ProfsTransposed)
head(its2ProfsList$samples)

its2ProfsNorm =  calcNormFactors(its2ProfsList, method = "TMM")
head(its2ProfsNorm$samples)
its2TMM = t(cpm(its2ProfsNorm, normalized.lib.sizes = TRUE))
its2ProfsNorm = cbind(its2Profs[,c(1:9)], its2TMM)

#Prepping for bar plot
colOrder2Norm = order(colSums(its2ProfsNorm[10:length(its2ProfsNorm[1,])]), decreasing = TRUE)+9 #Add the number of columns of metadata

its2NormProfsPerc = cbind(its2ProfsNorm[,c(1:9)],its2ProfsNorm[,c(colOrder2Norm)])
its2NormProfsPerc$sum = apply(its2NormProfsPerc[, c(10:length(its2NormProfsPerc[1,]))], 1, function(x) {
  sum(x, na.rm = T)
})

its2NormProfsPerc = cbind(its2NormProfsPerc[, c(1:9)], (its2NormProfsPerc[, c(10:(ncol(its2NormProfsPerc)-1))] 
                                                / its2NormProfsPerc$sum))
head(its2NormProfsPerc)

apply(its2NormProfsPerc[, c(10:(ncol(its2NormProfsPerc)))], 1, function(x) {
  sum(x, na.rm = T)
})

head(its2NormProfsPerc)


gssNormProf = otuStack(its2NormProfsPerc, count.columns = c(10:length(its2NormProfsPerc[1, ])),
                   condition.columns = c(1:9)) # remove summ rows
gssNormProf = gssNormProf %>% filter(otu != "summ") %>% droplevels()

levels(gssNormProf$otu)

levels(gssNormProf$depthZoneM)
levels(gssNormProf$depthZoneM) = c("10 m", "15 m", "25 m", "35 m")
levels(gssNormProf$site)
levels(gssNormProf$site) = c("Alacranes", "Bajos del Norte")
levels(gssNormProf$depthZoneM)
levels(gssNormProf$site)

#gssNormProf %>% arrange(site, depthZoneM, count, otu)

#Construct bar plot

its2NormProfsPlotA = ggplot(gssNormProf, aes(x = barPlotOrder, y = count, fill = factor(otu))) +
  geom_bar(position = "stack", stat = "identity", color = "black", size = 0.25) + 
  ylab("Proportion") +
  scale_fill_paletteer_d("dichromat::Categorical.12")+
  labs(fill = expression(paste(italic("ITS2"), " type profile"))) +
  guides(fill = guide_legend(ncol = 2, reverse = FALSE)) +
  facet_grid(depthZoneM ~ site, scales = "free_x") + #faceting plots by Depth and Site
  theme_bw()

its2NormProfsPlot = its2NormProfsPlotA +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "bottom",
        legend.title = element_text(color = "black", size = 12, hjust = 0.5, angle = 90),
        legend.text = element_text(color = "black", size = 10),
        legend.key = element_blank(),
        legend.key.size = unit(0.75,"line"),
        legend.background = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid=element_blank(),
        plot.background = element_blank(),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        strip.background = element_rect(fill = "white", size = 0.9)
  )

its2NormProfsPlot

ggsave("its2NormProfsPlot.png", plot = its2NormProfsPlot, width = 8.5, height = 8.5, unit = "in", dpi = 600)
####################################################################################################
#Profile bar plot, no normalization
#Prepping for bar plot
colOrder2 = order(colSums(its2Profs[10:length(its2Profs[1,])]), decreasing = TRUE)+9 #Add the number of columns of metadata

its2ProfsPerc = cbind(its2Profs[,c(1:9)],its2Profs[,c(colOrder2)])
its2ProfsPerc$sum = apply(its2ProfsPerc[, c(10:length(its2ProfsPerc[1,]))], 1, function(x) {
  sum(x, na.rm = T)
})

its2ProfsPerc = cbind(its2ProfsPerc[, c(1:9)], (its2ProfsPerc[, c(10:(ncol(its2ProfsPerc)-1))] 
                                                        / its2ProfsPerc$sum))
head(its2ProfsPerc)

apply(its2ProfsPerc[, c(10:(ncol(its2ProfsPerc)))], 1, function(x) {
  sum(x, na.rm = T)
})

head(its2ProfsPerc)


gssProf = otuStack(its2ProfsPerc, count.columns = c(10:length(its2ProfsPerc[1, ])),
                       condition.columns = c(1:9)) # remove summ rows
gssProf = gssProf %>% filter(otu != "summ") %>% droplevels()

levels(gssProf$otu)

levels(gssProf$depthZoneM)
levels(gssProf$depthZoneM) = c("10 m", "15 m", "25 m", "35 m")
levels(gssProf$site)
levels(gssProf$site) = c("Alacranes", "Bajos del Norte")
levels(gssProf$depthZoneM)
levels(gssProf$site)

#gssProf %>% arrange(site, depthZoneM, count, otu)

#Construct bar plot

its2ProfsPlotA = ggplot(gssProf, aes(x = barPlotOrder, y = count, fill = factor(otu))) +
  geom_bar(position = "stack", stat = "identity", color = "black", size = 0.25) + 
  ylab("Proportion") +
  scale_fill_paletteer_d("dichromat::Categorical.12")+
  labs(fill = expression(paste(italic("ITS2"), " type profile"))) +
  guides(fill = guide_legend(ncol = 2, reverse = FALSE)) +
  facet_grid(depthZoneM ~ site, scales = "free_x") + #faceting plots by Depth and Site
  theme_bw()

its2ProfsPlot = its2ProfsPlotA +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "bottom",
        legend.title = element_text(color = "black", size = 12, hjust = 0.5, angle = 90),
        legend.text = element_text(color = "black", size = 10),
        legend.key = element_blank(),
        legend.key.size = unit(0.75,"line"),
        legend.background = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid=element_blank(),
        plot.background = element_blank(),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        strip.background = element_rect(fill = "white", size = 0.9)
  )

its2ProfsPlot

ggsave("its2ProfsPlot.png", plot = its2ProfsPlot, width = 8.5, height = 8.5, unit = "in", dpi = 600)

####################################################################################################
#Stats Analysis NOT NORMALIZED

set.seed(694) #setting seed allows repetition of randomized processes

its2dispS = betadisper(vegdist(its2Profs[, c(10:ncol(its2Profs))]), its2Profs$site)

anova(its2dispS)
#No significant effect of site on beta diversity

set.seed(694)

its2dispD = betadisper(vegdist(its2Profs[, c(10:ncol(its2Profs))]), its2Profs$depthZoneM)

anova(its2dispD)
#No significant effect of depth on beta diversity

#PERMANOVA
set.seed(694)
its2Adonis = adonis(its2Profs[, c(10:ncol(its2Profs))] ~ depthZoneM * site, 
                    data = its2Profs, permutations = 9999, method = "bray")

set.seed(694)

its2PWAdonis = pairwise.adonis(its2Profs[, c(10:ncol(its2Profs))],
                               factors = its2Profs$depthZoneM,
                               sim.method = "bray", p.adjust.m = "BH", perm = 9999)

its2PWAdonis 

###################Stats Analysis Normalized

set.seed(694) #setting seed allows repetition of randomized processes

its2dispNormS = betadisper(vegdist(its2ProfsNorm[, c(10:ncol(its2ProfsNorm))]), its2ProfsNorm$site)

anova(its2dispNormS)
#No significant effect of site on beta diversity

set.seed(694)

its2dispNormD = betadisper(vegdist(its2ProfsNorm[, c(10:ncol(its2ProfsNorm))]), its2ProfsNorm$depthZoneM)

anova(its2dispNormD)
#No significant effect of depth on beta diversity

#PERMANOVA
set.seed(694)
its2NormAdonis = adonis(its2ProfsNorm[, c(10:ncol(its2ProfsNorm))] ~ depthZoneM * site, 
                    data = its2ProfsNorm, permutations = 9999, method = "bray")

its2NormAdonis
#####################################################
#NON-NORMALIZED NMDS
Com = its2Profs[,10:ncol(its2Profs)]
Matrix = as.matrix(Com)
set.seed(123)
Nmds = metaMDS(Matrix, distance = "bray")
Nmds
DataScores = as.data.frame(scores(Nmds))
DataScores$sample = its2Profs$sample
DataScores$site = its2Profs$site
DataScores$depthZoneM = its2Profs$depthZoneM
head(DataScores)
DataScores$popDepth <- paste(DataScores$site, "-", DataScores$depthZoneM)
zooxNmds = merge(DataScores, aggregate(cbind(mean.x=NMDS1,mean.y=NMDS2)~popDepth, DataScores, mean), by="popDepth")

zooxNmdsPlot = ggplot(zooxNmds, aes(x = NMDS1, y = NMDS2, fill = site, shape = depthZoneM)) +
  geom_point(size = 3, alpha = 0.3, (aes(x = NMDS1, y = NMDS2, fill = site, shape = depthZoneM)), show.legend = FALSE, guide=FALSE)+
  geom_point(aes(x = mean.x, y = mean.y, shape = depthZoneM), size = 5, color = "black") +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 12),
        legend.text = element_text(size = 12, face ="bold", colour ="black"),
        legend.position = "left", axis.title.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) +
  labs(x = "NMDS1", fill = "Population", y = "NMDS2", shape = "Depth Zone")  +
  scale_fill_manual(values = c("seagreen3", "blue"))+
  scale_shape_manual(values = c(24,22,23,25), labels=c("10 m","15 m", "25 m", "35 m"))+
  guides(shape = guide_legend(order = 2), fill = guide_legend(override.aes = list(shape = 22, size = 4, color = NA), order = 1))
#guides(fill = guide_legend(override.aes=list(shape=22)))


zooxNmdsPlot

ggsave("zooxNmdsPlot.tiff", plot = zooxNmdsPlot, height = 5, width = 7, units = "in", dpi = 300)

#make NORMALIZED community matrix - extract columns with abundance information
normCom = its2ProfsNorm[,10:ncol(its2ProfsNorm)]
normMatrix = as.matrix(normCom)
set.seed(123)
normNmds = metaMDS(normMatrix, distance = "bray")
normNmds
normDataScores = as.data.frame(scores(normNmds))
normDataScores$sample = its2ProfsNorm$sample
normDataScores$site = its2ProfsNorm$site
normDataScores$depthZoneM = its2ProfsNorm$depthZoneM
head(normDataScores)
normDataScores$popDepth <- paste(normDataScores$site, "-", normDataScores$depthZoneM)
zooxNormNmds = merge(normDataScores, aggregate(cbind(mean.x=NMDS1,mean.y=NMDS2)~popDepth, normDataScores, mean), by="popDepth")

zooxNormNmdsPlot = ggplot(zooxNormNmds, aes(x = NMDS1, y = NMDS2, fill = site, shape = depthZoneM)) +
  geom_point(size = 3, alpha = 0.3, (aes(x = NMDS1, y = NMDS2, fill = site, shape = depthZoneM)), show.legend = FALSE, guide=FALSE)+
  geom_point(aes(x = mean.x, y = mean.y, shape = depthZoneM), size = 5, color = "black") +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 12),
        legend.text = element_text(size = 12, face ="bold", colour ="black"),
        legend.position = "left", axis.title.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) +
  labs(x = "NMDS1", fill = "Population", y = "NMDS2", shape = "Depth Zone")  +
  scale_fill_manual(values = c("seagreen3", "blue"))+
  scale_shape_manual(values = c(24,22,23,25), labels=c("10 m","15 m", "25 m", "35 m"))+
  guides(shape = guide_legend(order = 2), fill = guide_legend(override.aes = list(shape = 22, size = 4, color = NA), order = 1))
  #guides(fill = guide_legend(override.aes=list(shape=22)))


zooxNormNmdsPlot

ggsave("zooxNormNmdsPlot.tiff", plot = zooxNormNmdsPlot, height = 5, width = 7, units = "in", dpi = 300)

###############################Zoox ITS2 PCoA
zooxProfs = its2Profs[,10:ncol(its2Profs)]
zooxDist=vegdist(zooxProfs, method="bray")

zooxMds = cmdscale(zooxDist, eig = TRUE, x.ret = TRUE)

# Determine percent variation captured on each axis
# Calculate the eigenvalues so later we can figure out % variation shown on each Principal Coordinate
mxZooxPcoaVar = round(zooxMds$eig/sum(zooxMds$eig)*100, 1)
mxZooxPcoaVar

# Format data to plot
mxZooxPcoaValues = zooxMds$points
mxZooxPcoaValues

zooxI2P = its2Profs[,1:9]
row.names(zooxI2P) = zooxI2P[,1]
mxZooxPcoaValues=cbind(zooxI2P, mxZooxPcoaValues)
mxZooxPcoaValues =as.data.frame(mxZooxPcoaValues, sample = rownames(mxZooxPcoaValues))
colnames(mxZooxPcoaValues)[10] <- "PCo1"
colnames(mxZooxPcoaValues)[11] <- "PCo2"
mxZooxPcoaValues$popDepth<- paste(mxZooxPcoaValues$site, "-", mxZooxPcoaValues$depthZoneM)
head(mxZooxPcoaValues)

zooxPCoA = merge(mxZooxPcoaValues, aggregate(cbind(mean.x=PCo1,mean.y=PCo2)~popDepth, mxZooxPcoaValues, mean), by="popDepth")

#snpPCoA$depthZoneM=as.factor(snpPCoA$depthZoneM)
#snpPCoA$depthZoneM = factor(snpPCoA$depthZoneM, levels(snpPCoA$depthZoneM)[c(2,1)])

# SNP PCoA biplot
mxZooxPcoaPlotA = ggplot(zooxPCoA, aes(x = PCo1, y = PCo2, color = site, fill = site, shape = depthZoneM, linetype = depthZoneM)) +
  geom_hline(yintercept = 0, color = "gray90", size = 0.5) +
  geom_vline(xintercept = 0, color = "gray90", size = 0.5) +
  #stat_ellipse(data = subset(zooxPCoA, type = "n", geom = "polygon")) + #ellipse
  scale_linetype_manual(values=c(1,2,4,3), name = "Depth Zone")+
  geom_point(aes(x = PCo1, y = PCo2, shape = depthZoneM), size = 3, alpha = 0.3, show.legend = FALSE, guide=FALSE, position=position_jitter()) + #individual's points indicated by circles
  #geom_jitter(alpha = 0.3, show.legend=FALSE)+
  scale_shape_manual(values = c(24,22,23,25), breaks=c("10","15", "25", "35"), name = "Depth Zone", labels=c("10 m","15 m", "25 m", "35 m")) +
  geom_point(aes(x = mean.x, y = mean.y, shape = depthZoneM), size = 5, color = "black") + #population centroids indicated by triangles
  #scale_fill_paletteer_d("LaCroixColoR::PassionFruit", breaks=c("Alacranes", "Bajos del Norte"), name= "Population")+
  scale_fill_manual(values = c("seagreen3", "blue"), name= "Population")+
  #scale_color_paletteer_d("LaCroixColoR::PassionFruit", breaks=c("Alacranes", "Bajos del Norte"), name= "Population", guide=FALSE) +
  scale_color_manual(values = c("seagreen3", "blue"), name= "Population")+
  xlab(paste ("PCo 1 (", mxZooxPcoaVar[1],"%)", sep = "")) + #Prints percent variation explained by first axis
  ylab(paste ("PCo 2 (", mxZooxPcoaVar[2],"%)", sep = "")) + #Prints percent variation explained by second axis
  guides(shape = guide_legend(order = 2), fill = guide_legend(override.aes = list(shape = 22, size = 4, color = NA), order = 1))+ #linetype = guide_legend(order = 3),
  theme_bw()

mxZooxPcoaPlot = mxZooxPcoaPlotA +
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



mxZooxPcoaPlot

ggsave("mxZooxITS2PcoaPlot.tiff", plot = mxZooxPcoaPlot, height = 5, width = 7, units = "in", dpi = 300)

#####################Normalized PCoA
###############################Zoox ITS2 PCoA
zooxNormProfs = its2ProfsNorm[,10:ncol(its2ProfsNorm)]
zooxNormDist=vegdist(zooxNormProfs, method="bray")

zooxNormMds = cmdscale(zooxNormDist, eig = TRUE, x.ret = TRUE)

# Determine percent variation captured on each axis
# Calculate the eigenvalues so later we can figure out % variation shown on each Principal Coordinate
mxZooxNormPcoaVar = round(zooxNormMds$eig/sum(zooxNormMds$eig)*100, 1)
mxZooxNormPcoaVar

# Format data to plot
mxZooxNormPcoaValues = zooxNormMds$points
mxZooxNormPcoaValues

zooxNormI2P = its2ProfsNorm[,1:9]
row.names(zooxNormI2P) = zooxNormI2P[,1]
mxZooxNormPcoaValues=cbind(zooxNormI2P, mxZooxNormPcoaValues)
mxZooxNormPcoaValues =as.data.frame(mxZooxNormPcoaValues, sample = rownames(mxZooxNormPcoaValues))
colnames(mxZooxNormPcoaValues)[10] <- "PCo1"
colnames(mxZooxNormPcoaValues)[11] <- "PCo2"
mxZooxNormPcoaValues$popDepth<- paste(mxZooxNormPcoaValues$site, "-", mxZooxNormPcoaValues$depthZoneM)
head(mxZooxNormPcoaValues)

zooxNormPCoA = merge(mxZooxNormPcoaValues, aggregate(cbind(mean.x=PCo1,mean.y=PCo2)~popDepth, mxZooxNormPcoaValues, mean), by="popDepth")

#snpPCoA$depthZoneM=as.factor(snpPCoA$depthZoneM)
#snpPCoA$depthZoneM = factor(snpPCoA$depthZoneM, levels(snpPCoA$depthZoneM)[c(2,1)])

# SNP PCoA biplot
mxZooxNormPcoaPlotA = ggplot(zooxNormPCoA, aes(x = PCo1, y = PCo2, color = site, fill = site, shape = depthZoneM, linetype = depthZoneM)) +
  geom_hline(yintercept = 0, color = "gray90", size = 0.5) +
  geom_vline(xintercept = 0, color = "gray90", size = 0.5) +
  #stat_ellipse(data = subset(zooxPCoA, type = "n", geom = "polygon")) + #ellipse
  scale_linetype_manual(values=c(1,2,4,3), name = "Depth Zone")+
  geom_point(aes(x = PCo1, y = PCo2, shape = depthZoneM), size = 3, alpha = 0.3, show.legend = FALSE, guide=FALSE) + #individual's points indicated by circles
  scale_shape_manual(values = c(24,22,23,25), breaks=c("10","15", "25", "35"), name = "Depth Zone", labels=c("10 m","15 m", "25 m", "35 m")) +
  geom_point(aes(x = mean.x, y = mean.y, shape = depthZoneM), size = 5, color = "black") + #population centroids indicated by triangles
  #scale_fill_paletteer_d("LaCroixColoR::PassionFruit", breaks=c("Alacranes", "Bajos del Norte"), name= "Population")+
  scale_fill_manual(values = c("seagreen3", "blue"), name= "Population")+
  #scale_color_paletteer_d("LaCroixColoR::PassionFruit", breaks=c("Alacranes", "Bajos del Norte"), name= "Population", guide=FALSE) +
  scale_color_manual(values = c("seagreen3", "blue"), name= "Population")+
  xlab(paste ("PCo 1 (", mxZooxPcoaVar[1],"%)", sep = "")) + #Prints percent variation explained by first axis
  ylab(paste ("PCo 2 (", mxZooxPcoaVar[2],"%)", sep = "")) + #Prints percent variation explained by second axis
  guides(shape = guide_legend(order = 2), fill = guide_legend(override.aes = list(shape = 22, size = 4, color = NA), order = 1))+ #linetype = guide_legend(order = 3),
  theme_bw()

mxZooxNormPcoaPlot = mxZooxNormPcoaPlotA +
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



mxZooxNormPcoaPlot

################################################################
#Now read in the ITS2 Sequences (not profiles)

its2Seq = read.csv("mxIts2SeqCounts.csv", header = TRUE, check.names = FALSE)
head(its2Seq)

its2Seq$depthZoneM = factor(its2Seq$depthZoneM, levels = c("10", "15", "25", "35"))
levels(its2Seq$depthZoneM)

its2Seq$site = factor(its2Seq$site, levels(its2Seq$site)[c(1,2)])
its2Seq = its2Seq[order(its2Seq$site, its2Seq$depthZoneM), ]

head(its2Seq)

its2Seq$Symbiodinium = rowSums(zooxMaps[2:6])
its2Seq$Breviolum = rowSums(zooxMaps[7:10])
its2Seq$Cladocopium = rowSums(zooxMaps[11:16])
its2Seq$Durusdinium = rowSums(zooxMaps[17:20])

goods = purgeOutliers(its2Seq, count.columns = 9:length(its2Seq), otu.cut = 0.0001, sampleZcut = -5)

its2SeqTransposed = t(goods[, 9:length(goods[1, ])])
its2SeqList = DGEList(counts = its2SeqTransposed)
head(its2SeqList$samples)

its2SeqNorm =  calcNormFactors(its2SeqList, method = "TMM")
head(its2SeqNorm$samples)
its2SeqTMM = t(cpm(its2SeqNorm, normalized.lib.sizes = TRUE))

its2SeqGoods= its2Seq %>%
  filter(sample != 'MX112')

its2SeqNorm = cbind(its2SeqGoods[,c(1:8)], its2SeqTMM)

set.seed(694) #setting seed allows repetition of randomized processes

its2SeqdispS = betadisper(vegdist(its2SeqNorm[, c(9:ncol(its2SeqNorm))]), its2SeqNorm$site)

anova(its2SeqdispS)
#No significant effect of site on beta diversity

set.seed(694)

its2SeqdispD = betadisper(vegdist(its2SeqNorm[, c(9:ncol(its2SeqNorm))]), its2SeqNorm$depthZoneM)

anova(its2SeqdispD)
#No significant effect of depth on beta diversity

#PERMANOVA
set.seed(694)
its2SeqAdonis = adonis(its2SeqNorm[, c(9:ncol(its2SeqNorm))] ~ depthZoneM * site, 
                    data = its2SeqNorm, permutations = 9999, method = "bray")

its2SeqAdonis

colOrder = order(colSums(its2SeqNorm[9:length(its2SeqNorm[1,])]), decreasing = TRUE)+8 #Add the number of columns of metadata

its2SeqPerc = cbind(its2SeqGoods[,c(1:8)],its2SeqNorm[,c(colOrder)])
its2SeqPerc$sum = apply(its2SeqPerc[, c(9:length(its2SeqPerc[1,]))], 1, function(x) {
  sum(x, na.rm = T)
})

its2SeqPerc = cbind(its2SeqPerc[, c(1:8)], (its2SeqPerc[, c(9:(ncol(its2SeqPerc)-1))] 
                                                / its2SeqPerc$sum))

apply(its2SeqPerc[, c(9:(ncol(its2SeqPerc)))], 1, function(x) {
  sum(x, na.rm = T)
})

head(its2SeqPerc)

sampleCounts = plyr::count(its2SeqPerc, c('site','depthZoneM'))
meltedList = melt(lapply(sampleCounts$freq,function(x){c(1:x)}))
its2SeqPerc$barPlotOrder = meltedList$value
its2SeqPerc=its2SeqPerc[c(1,ncol(its2SeqPerc),2:(ncol(its2SeqPerc)-1))]


gssSeq = otuStack(its2SeqPerc, count.columns = c(10:length(its2SeqPerc[1, ])),
                   condition.columns = c(1:9)) # remove summ rows
gssSeq = gssSeq %>% filter(otu != "summ") %>% droplevels()

levels(gssSeq$otu)

levels(gssSeq$depthZoneM)
levels(gssSeq$depthZoneM) = c("10 m", "15 m", "25 m", "35 m")
levels(gssSeq$site)
levels(gssSeq$site) = c("Alacranes", "Bajos del Norte")
levels(gssSeq$depthZoneM)
levels(gssSeq$site)

#gssSeq %>% arrange(site, depthZoneM, count, otu)

#Construct bar plot
seqPal=paletteer_c("viridis::viridis", n = 88)

its2SeqPlotA = ggplot(gssSeq, aes(x = barPlotOrder, y = count, fill = factor(otu))) +
  geom_bar(position = "stack", stat = "identity", color = "black", size = 0.25) + 
  ylab("Proportion") +
  scale_fill_manual(values=seqPal)+
  labs(fill = expression(paste(italic("ITS2"), " type sequence"))) +
  guides(fill = guide_legend(ncol = 2, reverse = FALSE)) +
  facet_grid(depthZoneM ~ site, scales = "free_x") + #faceting plots by Depth and Site
  theme_bw()

its2SeqPlot = its2SeqPlotA +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "none",
        #legend.position = "bottom",
        legend.title = element_text(color = "black", size = 12, hjust = 0.5, angle = 90),
        legend.text = element_text(color = "black", size = 10),
        legend.key = element_blank(),
        legend.key.size = unit(0.75,"line"),
        legend.background = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        strip.background = element_rect(fill = "white", size = 0.9),
        panel.grid = element_blank()
  )

its2SeqPlot
##############
comSeq = its2SeqNorm[,9:ncol(its2SeqNorm)]
m_comSeq = as.matrix(comSeq)
set.seed(123)
nmdsSeq = metaMDS(m_comSeq, distance = "bray")
nmdsSeq
scoresSeq = as.data.frame(scores(nmdsSeq))
scoresSeq$sample = its2SeqNorm$sample
scoresSeq$site = its2SeqNorm$site
scoresSeq$depthZoneM = its2SeqNorm$depthZoneM
head(scoresSeq)
scoresSeq$popDepth <- paste(scoresSeq$site, "-", scoresSeq$depthZoneM)
zooxSeqNmds = merge(scoresSeq, aggregate(cbind(mean.x=NMDS1,mean.y=NMDS2)~popDepth, scoresSeq, mean), by="popDepth")

zooxSeqNmdsPlot = ggplot(zooxSeqNmds, aes(x = NMDS1, y = NMDS2, fill = site, shape = depthZoneM)) +
  geom_point(size = 3, alpha = 0.3, (aes(x = NMDS1, y = NMDS2, fill = site, shape = depthZoneM)), show.legend = FALSE, guide=FALSE)+
  geom_point(aes(x = mean.x, y = mean.y, shape = depthZoneM), size = 5, color = "black") +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 12),
        legend.text = element_text(size = 12, face ="bold", colour ="black"),
        legend.position = "left", axis.title.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) +
  labs(x = "NMDS1", fill = "Population", y = "NMDS2", shape = "Depth Zone")  +
  scale_fill_manual(values = c("seagreen3", "blue"))+
  scale_shape_manual(values = c(24,22,23,25), labels=c("10 m","15 m", "25 m", "35 m"))+
  guides(shape = guide_legend(order = 2), fill = guide_legend(override.aes = list(shape = 22, size = 4, color = NA), order = 1))
#guides(fill = guide_legend(override.aes=list(shape=22)))


zooxSeqNmdsPlot

ggsave("zooxSeqNmdsPlot.tiff", plot = zooxSeqNmdsPlot, height = 5, width = 7, units = "in", dpi = 300)

######################################################################################################
zooxSeqs = its2SeqNorm[,9:ncol(its2SeqNorm)]
zooxSeqDist=vegdist(zooxSeqs, method="bray")

zooxSeqMds = cmdscale(zooxSeqDist, eig = TRUE, x.ret = TRUE)

# Determine percent variation captured on each axis
# Calculate the eigenvalues so later we can figure out % variation shown on each Principal Coordinate
mxZooxSeqPcoaVar = round(zooxSeqMds$eig/sum(zooxSeqMds$eig)*100, 1)
mxZooxSeqPcoaVar

# Format data to plot
mxZooxSeqPcoaValues = zooxSeqMds$points
mxZooxSeqPcoaValues

zooxSeqI2P = its2SeqNorm[,1:8]
row.names(zooxSeqI2P) = zooxSeqI2P[,1]
mxZooxSeqPcoaValues=cbind(zooxSeqI2P, mxZooxSeqPcoaValues)
mxZooxSeqPcoaValues =as.data.frame(mxZooxSeqPcoaValues, sample = rownames(mmxZooxSeqPcoaValues))
colnames(mxZooxSeqPcoaValues)[9] <- "PCo1"
colnames(mxZooxSeqPcoaValues)[10] <- "PCo2"
mxZooxSeqPcoaValues$popDepth<- paste(mxZooxSeqPcoaValues$site, "-", mxZooxSeqPcoaValues$depthZoneM)
head(mxZooxSeqPcoaValues)

zooxSeqPCoA = merge(mxZooxSeqPcoaValues, aggregate(cbind(mean.x=PCo1,mean.y=PCo2)~popDepth, mxZooxSeqPcoaValues, mean), by="popDepth")

#snpPCoA$depthZoneM=as.factor(snpPCoA$depthZoneM)
#snpPCoA$depthZoneM = factor(snpPCoA$depthZoneM, levels(snpPCoA$depthZoneM)[c(2,1)])

# SNP PCoA biplot
mxZooxSeqPcoaPlotA = ggplot(zooxSeqPCoA, aes(x = PCo1, y = PCo2, color = site, fill = site, shape = depthZoneM, linetype = depthZoneM)) +
  geom_hline(yintercept = 0, color = "gray90", size = 0.5) +
  geom_vline(xintercept = 0, color = "gray90", size = 0.5) +
  #stat_ellipse(data = subset(zooxPCoA, type = "n", geom = "polygon")) + #ellipse
  scale_linetype_manual(values=c(1,2,4,3), name = "Depth Zone")+
  geom_point(aes(x = PCo1, y = PCo2, shape = depthZoneM), size = 3, alpha = 0.3, show.legend = FALSE, guide=FALSE) + #individual's points indicated by circles
  scale_shape_manual(values = c(24,22,23,25), breaks=c("10","15", "25", "35"), name = "Depth Zone", labels=c("10 m","15 m", "25 m", "35 m")) +
  geom_point(aes(x = mean.x, y = mean.y, shape = depthZoneM), size = 5, color = "black") + #population centroids indicated by triangles
  #scale_fill_paletteer_d("LaCroixColoR::PassionFruit", breaks=c("Alacranes", "Bajos del Norte"), name= "Population")+
  scale_fill_manual(values = c("seagreen3", "blue"), name= "Population")+
  #scale_color_paletteer_d("LaCroixColoR::PassionFruit", breaks=c("Alacranes", "Bajos del Norte"), name= "Population", guide=FALSE) +
  scale_color_manual(values = c("seagreen3", "blue"), name= "Population")+
  xlab(paste ("PCo 1 (", mxZooxPcoaVar[1],"%)", sep = "")) + #Prints percent variation explained by first axis
  ylab(paste ("PCo 2 (", mxZooxPcoaVar[2],"%)", sep = "")) + #Prints percent variation explained by second axis
  guides(shape = guide_legend(order = 2), fill = guide_legend(override.aes = list(shape = 22, size = 4, color = NA), order = 1))+ #linetype = guide_legend(order = 3),
  theme_bw()

mxZooxSeqPcoaPlot = mxZooxSeqPcoaPlotA +
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



mxZooxSeqPcoaPlot
ggsave("mxZooxSeqPcoaPlot.tiff", plot = mxZooxSeqPcoaPlot, height = 5, width = 7, units = "in", dpi = 300)

#####Combined NGSAdmix/ITS2 Plot
head(mxNgsAdmix) #Dataframe for NGSAdmix plot
head(its2NormProfsPerc) #Dataframe for ITS2 Profile plots

zooxAdmix=full_join(mxNgsAdmix, its2NormProfsPerc, by = c("sample", "site", "depthZoneM"), copy = FALSE, suffix = c(".ngs", ".prof"))

zooxAdmix=select(zooxAdmix, -c(lat, lon, collectionDate, collectionDepth, depthZone, barPlotOrder.ngs, barPlotOrder.prof))
zooxAdmix[is.na(zooxAdmix)] <- 0

zooxAdmix = zooxAdmix %>% arrange(site, depthZoneM, -cluster3.2)

sampleCounts = plyr::count(zooxAdmix, c('site','depthZoneM'))
meltedList = reshape2::melt(lapply(sampleCounts$freq,function(x){c(1:x)}))
zooxAdmix$barPlotOrder = meltedList$value

zooxAdmixMelt = melt(zooxAdmix, id.vars=c("sample", "site", "depthZoneM", "barPlotOrder"), variable.name="Cluster", value.name="Fraction")
#zooxAdmixMelt$Fraction[is.na(zooxAdmixMelt$Fraction)] <- 0

p3Stack = ggplot(data = subset(zooxAdmixMelt, subset = zooxAdmixMelt$Cluster %in% c("cluster3.1", "cluster3.2", "cluster3.3")), aes(x=barPlotOrder, y=Fraction, fill=Cluster, order=barPlotOrder)) +
  geom_bar(stat="identity", position="fill", width=1, colour="grey25") +
  #facet_grid(.~pop, scales = "free", switch = "x", space = "free") +
  facet_grid(depthZoneM ~ site, scales = "free") + #faceting plots by Depth and Site
  labs(x = "Population", y = "Ancestry") +
  ggtitle(expression(paste(italic("K"), "=3"))) +
  theme(plot.title = element_text(size=22),
        panel.grid=element_blank(),
        panel.background=element_rect(fill=NA, colour="grey25"),
        panel.spacing.x=grid:::unit(0, "lines"),
        panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        axis.text.x=element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        strip.background=element_blank(),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        #strip.text=element_text(size=20, angle=90),
        legend.key=element_blank(),
        legend.position = "none",
        legend.title = element_blank()) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  #scale_fill_manual(values = brewer.pal(name = "YlGnBu", n = 4), name = "Cluster") +
  scale_fill_manual(values = colPal3, name = "Cluster") +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))
p3Stack


its2NormProfsPlotStackA = ggplot(data = subset(zooxAdmixMelt, subset = zooxAdmixMelt$Cluster %in% c(
  "G3b",
                                                                                                "G3am",
                                                                                                "C3/C3b",
                                                                                                "C3-C3fc-C40g-C3b-C3an-C21-C3fd-C3s-C3bb",
                                                                                                "C3an/C3-C21-C3bb", 
                                                                                                "C3-C21-C3gu-C3gy-C3gw-C3gx-C3gv-C3b", 
                                                                                                "C3-C3an-C3de-C21-C3bb-C3fc-C3b",
                                                                                                "C3-C3de-C3bb-C21ae-C21-C3an-C3s", 
                                                                                                "C3-C3fc-C21-C3b-C3an-C3fd-C3bb", 
                                                                                                "C3-C3de-C21-C3bb-C3an-C21ae-C65b-C3s", 
                                                                                                "C3-C21-C3an-C3fc-C3gz-C3b-C3bb", 
                                                                                                "C3-C21-C3fc-C3an-C65b-C3b-C3bb")), aes(x=barPlotOrder, y=Fraction, fill=Cluster, order=barPlotOrder)) +
  geom_bar(stat="identity", position="fill", width=1, colour="grey25") +
  ylab("Proportion") +
  scale_fill_paletteer_d("ggsci::default_gsea")+
  labs(fill = expression(paste(italic("ITS2"), " type profile"))) +
  guides(fill = guide_legend(ncol = 2, reverse = FALSE)) +
  facet_grid(depthZoneM ~ site, scales = "free_x") + #faceting plots by Depth and Site
  theme_bw()

its2NormProfsPlotStack = its2NormProfsPlotStackA +
  # theme(axis.title.x = element_blank(),
  #       axis.text.x = element_blank(),
  #       axis.ticks.x = element_blank(),
  #       axis.title.y = element_text(color = "black", size = 12),
  #       axis.text.y = element_text(color = "black", size = 12),
  #       legend.position = "bottom",
  #       legend.title = element_text(color = "black", size = 12, hjust = 0.5, angle = 90),
  #       legend.text = element_text(color = "black", size = 10),
  #       legend.key = element_blank(),
  #       legend.key.size = unit(0.75,"line"),
  #       legend.background = element_blank(),
  #       panel.border = element_blank(),
  #       panel.background = element_rect(fill = "white"),
  #       plot.background = element_blank(),
  #       strip.text.x = element_text(size = 12),
  #       strip.text.y = element_text(size = 12),
  #       strip.background = element_rect(fill = "white", size = 0.9)
  # )

theme(plot.title = element_text(size=22),
      panel.grid=element_blank(),
      panel.background=element_rect(fill=NA, colour="grey25"),
      panel.spacing.x=grid:::unit(0, "lines"),
      panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
      axis.text.x=element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.x=element_blank(),
      axis.ticks.y=element_blank(),
      strip.background=element_blank(),
      strip.text.x = element_text(size = 12),
      strip.text.y = element_text(size = 12),
      #strip.text=element_text(size=20, angle=90),
      legend.key=element_blank(),
      legend.position = "none",
      legend.title = element_blank()) +
  scale_x_discrete(expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0)) +
  #scale_fill_manual(values = brewer.pal(name = "YlGnBu", n = 4), name = "Cluster") +
  #scale_fill_manual(values = colPal3, name = "Cluster") +
  guides(fill=guide_legend(override.aes=list(colour=NULL)))

its2NormProfsPlotStack

stackedPlot= (p3Stack/its2NormProfsPlotStack)
stackedPlot

