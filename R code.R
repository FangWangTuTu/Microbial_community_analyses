
#### Figure 1a Soil moisture content ####
library(ggplot2) 
library(patchwork)
data1 <- read.csv("Data/soil moisture content.csv", header=TRUE, row.names=1)
data1$Treatment <- factor(data1$Treatment, levels = c("CK", "N", "NPK", "NPKM", "NPK1.5M"))
a <- ggplot(data1, aes(x = Time, y = SMC, col = Time)) +  
  scale_color_gradient(low='blue', high='red') +
  geom_point(size = 0.8) +
  ylab("Soil moisture content (%)") +
  xlab("Air-drying and archiving time (h)") +
  theme(axis.title.x = element_blank()) + 
  facet_wrap(~ Treatment, ncol = 5) +
  labs(title = "Soil moisture content", tag = "a") +
  theme(plot.title = element_text(size = 15, colour = "black", hjust = 0.5, face = "bold")) 
a




#### Prokaryotic OTU Rarefaction ####
library(vegan)
data <- read.csv("Data/16s_OTU.csv", header=TRUE, row.names=1)
tdata <- t(data)
write.csv(tdata, file='Results/16s OTU_t.csv')

set.seed(88888) 
Bfn <- read.csv('Results/16s OTU_t.csv', header=TRUE, row.names=1)
S1 <- specnumber(Bfn)
(raremax1 <- min(rowSums(Bfn)))    
newOTU <- rrarefy(Bfn, raremax1)
write.csv(newOTU, file='Results/16s OTU_t_re.csv')

tnewOTU <- t(newOTU)
write.csv(tnewOTU, file='Results/16s OTU_t_re_tAA.csv')

# Transformation
data <- read.csv("Results/16s OTU_t_re_tAA.csv", row.names = 1, header = T) 
data <- decostand(data, MARGIN = 2, method = "total") 
write.csv(data, "Results/16s OTU_RA.csv")  

data <- read.csv("Results/16s OTU_RA.csv", row.names = 1, header = T) 
data <- data^0.5
write.csv(data, "Results/16s OTU_rootsign.csv") 


#### Figure S1 Prokaryotic community composition ####
# Prokaryotic Abundance at phylum level
OTU<-read.csv("Data/16s OTU_RA & taxa.csv", header=TRUE)
Phylum = aggregate(x = OTU[, 8:227], by=list(OTU$Phylum), FUN='sum')
write.csv(Phylum, file = "Results/16s Phyla.csv")

# plot
library(ggplot2)
library(ggalluvial)

taxon <- read.csv("Results/16s Phyla top14 for plot.csv", header=TRUE)  

taxon <- data.frame(subset(taxon, Treatment=="CK"))    #Choose different treatments
#taxon <- data.frame(subset(taxon, Treatment=="N"))    #Choose different treatments
#taxon <- data.frame(subset(taxon, Treatment=="NPK"))    #Choose different treatments
#taxon <- data.frame(subset(taxon, Treatment=="NPKM"))    #Choose different treatments
#taxon <- data.frame(subset(taxon, Treatment=="NPK1.5M"))    #Choose different treatments

Palette <- c("#F7FBFF", "#C6DBEF", "#6BAED6", "#4292C6", "#2171B5", "#08519C", 
             "#F7FCF5", "#C7E9C0", "#74C476", "#238B45", "#006D2C", 
             "#FFF5F0", "#FCBBA1", "#FB6A4A", "#CB181D", "#A50F15" )
taxon$Treatment <- factor(taxon$Treatment, 
                          levels = c("CK", "N", "NPK", "NPKM", "NPK1.5M"), 
                          ordered=TRUE)
taxon$Time <- factor(taxon$Time, 
                     levels = c("T0", "T16", "T32", "T64", "T128", "T256", "T512", "T1024", "T2048", "T4096", "T8192"), 
                     ordered=TRUE)
taxon$Taxon <- factor(taxon$Taxon, 
                      levels = c("Acidobacteria", "Proteobacteria", "Chloroflexi", "Thaumarchaeota", "Planctomycetes",
                                 "Actinobacteria", "Bacteroidetes", "Gemmatimonadetes", "Firmicutes", "Nitrospirae",
                                 "Armatimonadetes", "Cyanobacteria", "Latescibacteria", "Verrucomicrobia", "Other"), 
                      ordered=TRUE)
p1 <- ggplot(data = taxon, aes(x = Time, y = RA, alluvium = Taxon, stratum = Taxon)) +
  geom_alluvium(aes(fill = Taxon), alpha = 1, width = 0.2, show.legend = F) +  
  geom_stratum(aes(fill = Taxon), width = 0.2, show.legend = F) +
  scale_fill_manual(values = Palette) +
  scale_x_discrete(limits = taxon$Time) +
  ylab(label = "Relative Abundance") + 
  xlab(label = "") +
  labs(title = "Prokaryotic community") +
  theme(plot.title = element_text(size = 20, colour = "black", hjust = 0, vjust = 0.5, face = "bold"))
p1
 


#### Figure 1b Prokaryotic community similarity_time series ####
# Bray-Curtis dissimilarity
library(vegan)
library(ape)
library(ggplot2)
library(ggrepel)

data <- read.csv("Results/16s OTU_rootsign.csv", header=TRUE, row.names=1)
data <- t(data)
otu.bc <- vegdist(data, 'bray')                 
write.csv(as.matrix(otu.bc), 'Results/16s_Dissimilarity.csv')

# Bray-Curtis similarity
data <- read.csv("Results/16s_Dissimilarity.csv", row.names = 1, header = T) 
data <- 1-data 
write.csv(data,"Results/16s_Similarity.csv")

# Plot
library(ggplot2) 
library(patchwork)
data2 <- read.csv("Results/16s_Similarity Time series for plot.csv", header=TRUE, row.names=1)
data2$Treatment <- factor(data2$Treatment, levels = c("CK", "N", "NPK", "NPKM", "NPK1.5M"))
b <- ggplot(data2, aes(x = Time, y = Similarity, col = Time)) +   
  scale_color_gradient(low='blue', high='red') +
  geom_point(size = 0.8) +
  ylab("Bray-Curtis similarity") +
  xlab("Air-drying and archiving time (h)") +
  theme(axis.title.x = element_blank()) + 
  theme(plot.title = element_text(size = 15, colour = "black", hjust = 0.5, face = "bold")) +
  scale_y_continuous(limits=c(0.4, 0.8), breaks = c( 0.4, 0.6, 0.8), label = c("0.4", "0.6", "0.8")) + 
  facet_wrap(~ Treatment, ncol = 5) +
  labs(title = "Prokaryotic community", tag = "b")
  theme(plot.title = element_text(size = 20, colour = "black", hjust = 0, vjust = 0.5, face = "bold"))
b

#### Figure S3abc Prokaryotic community similarity_time series & replication & treatment effects ####
library(ggplot2) 
library(patchwork) 

# Prokaryotic community similarity_Air-drying and archiving time series
data <- read.csv("Results/16s_Similarity Time series for plot.csv", header=TRUE, row.names=1)
data$Treatment <- factor(data$Treatment, c("CK", "N", "NPK", "NPKM", "NPK1.5M"), ordered=TRUE)
colour5 <- c("#F7FBFF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#2171B5")
a <- ggplot(data, aes(x = Treatment, y = Similarity, fill = Treatment)) +  
  geom_boxplot() +
  geom_jitter(size = 1) +   
  scale_fill_manual(values = colour5) +  
  scale_y_continuous(limits=c(0.3, 0.9),    
                     breaks = c(0.4, 0.6, 0.8), 
                     label = c("0.4", "0.6", "0.8")) + 
  ylab("Bray-Curtis similarity") +
  xlab("Treatment") +
  theme(axis.title.x = element_blank()) +   
  theme(legend.position = "none") +  
  labs(title = "Prokaryotic community") + 
  labs(subtitle = "Air-drying and archiving time series", tag = "a") +
  theme(plot.subtitle = element_text(size = 12, colour = "black", hjust = 0, vjust = 0.5, face = "bold")) 
a

# Prokaryotic community similarity_Cryopreservation (T = 0 h): replication effect
data <- read.csv("Results/16s_Similarity replicate for plot.csv", header=TRUE, row.names=1)
data$Treatment <- factor(data$Treatment, levels = c("CK", "N", "NPK", "NPKM", "NPK1.5M"))
colour5 <- c("#F7FBFF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#2171B5")
data1 <- subset(data, Time=="0")   
b <- ggplot(data1, aes(x = Treatment, y = Similarity, fill = Treatment)) + 
  geom_boxplot() +
  geom_jitter(size = 1) +    
  scale_fill_manual(values = colour5) +  
  scale_y_continuous(limits=c(0.3, 0.9),    
                     breaks = c(0.4, 0.6, 0.8), 
                     label = c("0.4", "0.6", "0.8")) + 
  ylab("Bray-Curtis similarity") +
  theme(axis.title.y = element_blank()) +   
  xlab("Treatment") +
  theme(axis.title.x = element_blank()) +    
  theme(legend.position = "none") +          
  labs(title = "Prokaryotic community") + 
  labs(subtitle = "Cryopreservation (T = 0 h): replication effect", tag = "b") +
  theme(plot.subtitle = element_text(size = 12, colour = "black", hjust = 0, vjust = 0.5, face = "bold"))  
b

# Prokaryotic community similarity_Cryopreservation (T = 0 h): treatment effect
data <- read.csv("Results/16s_Similarity Treatment effects for plot.csv", header=TRUE, row.names=1)
data$Treatment <- factor(data$Treatment, c("N", "NPK", "NPKM", "NPK1.5M"), ordered=TRUE)
colour4 <- c("#C6DBEF", "#9ECAE1", "#6BAED6", "#2171B5")
c <- ggplot(data, aes(x = Treatment, y = Similarity, fill = Treatment)) +   
  geom_boxplot() +
  geom_jitter(size = 1) +    
  scale_fill_manual(values = colour4) + 
  scale_y_continuous(limits=c(0.3, 0.9),  
                     breaks = c(0.4, 0.6, 0.8), 
                     label = c("0.4", "0.6", "0.8")) + 
  ylab("Bray-Curtis similarity") +
  theme(axis.title.y = element_blank()) +  
  xlab("Treatment") +
  theme(axis.title.x = element_blank()) +   
  theme(legend.position = "none") +     
  labs(title = "Prokaryotic community") + 
  labs(subtitle = "Cryopreservation (T = 0 h): treatment effect", tag = "c") +
  theme(plot.subtitle = element_text(size = 12, colour = "black", hjust = 0, vjust = 0.5, face = "bold")) 
c

# wilcox.test: to compare similarity of time series & treatment effects
data <- read.csv("Results/16s_Similarity for test.csv", header=TRUE, row.names=1)
wilcox.test(data$TSB, data$TEB) 
wilcox.test(data$TSC, data$TEC) 
wilcox.test(data$TSD, data$TED) 
wilcox.test(data$TSE, data$TEE) 


#### Figure 2a Prokaryotic community structure_Air-drying and archiving time series ####
# NMDS
library(vegan)
OTU <- read.csv("Results/16s OTU_rootsign.csv", header=TRUE, row.names=1)
tOTU <- t(OTU)

set.seed(8888)  
tOTU_NMDS1 <- metaMDS(tOTU, distance = 'bray', k = 2, autotransform = FALSE, trymax = 2000)  
tOTU_NMDS <- metaMDS(tOTU, distance = "bray", k = 2, autotransform = FALSE, previous.best = tOTU_NMDS1)

Stress <- tOTU_NMDS$stress 
Stress

NMDS_scores <- scores(tOTU_NMDS, choices = c(1,2))  
write.csv(NMDS_scores, file="Results/16s NMDS Scores_Air-drying and archiving time series.csv")

# PERMANOVA/ADONIS
library(vegan)
otu <- read.csv("Results/16s OTU_rootsign.csv", header=TRUE, row.names=1)
otu <- data.frame(t(otu))  

group <- read.csv('Data/groups.csv', header=TRUE, row.names=1)

set.seed(8888)
adonis_result <- adonis(otu~Treatment, group, distance = 'bray', permutations = 2000)
summary(adonis_result)

R2 <- adonis_result$aov.tab$R2[1]
R2

P <- adonis_result$aov.tab$"Pr(>F)"[1]
P 

otuput <- data.frame(adonis_result$aov.tab, check.names = FALSE, stringsAsFactors = FALSE)
otuput <- cbind(rownames(otuput), otuput)
names(otuput) <- c('', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
write.table(otuput, file = 'Results/16s PERMANOVA_Air-drying and archiving time series.csv', row.names = FALSE, sep = ",")

# Plot
data <- read.csv("Results/16s NMDS Scores_Air-drying and archiving time series for plot.csv", header=TRUE)

data$Treatment <- factor(data$Treatment, c("CK", "N", "NPK", "NPKM", "NPK1.5M"), ordered=TRUE)

library(gridExtra) 
label1 <- c(paste("Air-drying and archiving time series"))
label2 <- c(paste("Stress ==", round(Stress, 3),  sep = ""))
label3 <- c(paste("PERMANOVA"),
            paste("R2 ==", round(R2, 3),  sep = ""),  
            paste("P < 0.001 ", sep = ""))

thm1 <- ttheme_minimal(core = list(fg_params = list(parse = TRUE, col = "black", hjust = 0, fontsize =10, x = 0.4)))
thm2 <- ttheme_minimal(core = list(fg_params = list(parse = TRUE, col = "black", hjust = 0,  x = 0.1)))

tbgrob1 <- tableGrob(as.matrix(label1), theme = thm1)
tbgrob2 <- tableGrob(as.matrix(label2), theme = thm2)
tbgrob3 <- tableGrob(as.matrix(label3), theme = thm2)

library(gridExtra) 
library(ggplot2) 
a <- ggplot(data = data) +
  geom_point(aes(x = NMDS1, y = NMDS2, col = Time, shape = Treatment), size=3) +
  scale_colour_gradient(low = "blue", high = "red") +
  scale_shape_manual(values = c(1, 3, 11, 4, 0)) +  
  scale_x_continuous(limits = c(-0.50, 0.50)) +  
  scale_y_continuous(limits = c(-0.30, 0.30)) + 
  xlab("NMDS1") +
  ylab("NMDS2") +
  labs(title = "Prokaryotic community", tag = "a") +
  theme(plot.title = element_text(size = 20, colour = "black", hjust = 0, vjust = 0.5, face = "bold")) +
  annotation_custom(tbgrob1, xmin = -1.38, ymax = 0.90) +   
  annotation_custom(tbgrob2, xmin = -1.1, ymax = -0.2) +    
  annotation_custom(tbgrob3, xmin = 0.12, ymax = -0.15)
a


#### Figure 2c Prokaryotic community structure_Cryopreservation (T=0h) ####
library(vegan)
library(ggplot2)  
library(dplyr) 
library(ggplot2)
library(ggrepel)  #geom_text_repel
library(dplyr)  #filter
library(cowplot) #拼图
library(gridExtra) 
# NMDS
data <- read.csv("Results/16s OTU_rootsign_sort by time.csv", header=TRUE, row.names=1)
#T=0h
data0 <- data[,c(1:20)] 
tOTU <- data.frame(t(data0))

set.seed(8888)  
tOTU_NMDS1 <- metaMDS(tOTU, distance = 'bray', k = 2, autotransform = FALSE, trymax = 2000)  
tOTU_NMDS <- metaMDS(tOTU, distance = "bray", k = 2, autotransform = FALSE, previous.best = tOTU_NMDS1)

Stress <- tOTU_NMDS$stress 
Stress

NMDS_scores=scores(tOTU_NMDS, choices = c(1,2))  
write.csv(NMDS_scores, file="Results/16s NMDS Scores_Cryopreservation (T=0h).csv")

# PERMANOVA/ADONIS
data <- read.csv("Results/16s OTU_rootsign_sort by time.csv", header=TRUE, row.names=1)
#T=0h
data0 <- data[,c(1:20)] 
tOTU <- data.frame(t(data0))

group_all <- read.csv('Data/groups.csv', header=TRUE, row.names=1)
group <- dplyr::filter(group_all, Time == "0")

set.seed(8888)
adonis_result <- adonis(tOTU~Treatment, group, distance = 'bray', permutations = 2000)
summary(adonis_result)

R2 <- adonis_result$aov.tab$R2[1]
R2

P <- adonis_result$aov.tab$"Pr(>F)"[1]
P 

otuput <- data.frame(adonis_result$aov.tab, check.names = FALSE, stringsAsFactors = FALSE)
otuput <- cbind(rownames(otuput), otuput)
names(otuput) <- c('', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
write.table(otuput, file = 'Results/16s PERMANOVA_Cryopreservation (T=0h).csv', row.names = FALSE, sep = ",")

# Plot
data <- read.csv("Results/16s NMDS Scores_Cryopreservation (T=0h) for plot.csv", header=TRUE)

data$Treatment <- factor(data$Treatment, c("CK", "N", "NPK", "NPKM", "NPK1.5M"), ordered=TRUE)

library(gridExtra) 
label1 <- c(paste("Time series"))
label2 <- c(paste("Stress ==", round(Stress, 3),  sep = ""))
label3 <- c(paste("PERMANOVA"),
            paste("R2 ==", round(R2, 3),  sep = ""),  
            paste("P < 0.001 ", sep = ""))

thm1 <- ttheme_minimal(core = list(fg_params = list(parse = TRUE, col = "black", hjust = 0, fontsize =10, x = 0.2)))
thm2 <- ttheme_minimal(core = list(fg_params = list(parse = TRUE, col = "black", hjust = 0,  x = 0.1)))

tbgrob1 <- tableGrob(as.matrix(label1), theme = thm1)
tbgrob2 <- tableGrob(as.matrix(label2), theme = thm2)
tbgrob3 <- tableGrob(as.matrix(label3), theme = thm2)

c <- ggplot(data = data) +
  geom_point(aes(x = NMDS1, y = NMDS2, shape = Treatment), col = "blue", size = 3) +
  scale_shape_manual(values = c(1, 3, 11, 4, 0)) +  
  scale_x_continuous(limits = c(-0.50, 0.50)) +  
  scale_y_continuous(limits = c(-0.30, 0.30)) + 
  xlab("NMDS1") +
  ylab("NMDS2") +
  labs(title = "Prokaryotic community", tag = "c") +
  theme(plot.title = element_text(size = 20, colour = "black", hjust = 0, vjust = 0.5, face = "bold")) +
  annotation_custom(tbgrob1, xmin = -1.38, ymax = 0.90) +   
  annotation_custom(tbgrob2, xmin = -1.1, ymax = -0.2) +    
  annotation_custom(tbgrob3, xmin = 0.12, ymax = -0.15)
c
 





#### Fungal OTU Rarefaction ####
library(vegan)
data <- read.csv("Data/ITS_OTU.csv", header=TRUE, row.names=1)

tdata <- t(data)
write.csv(tdata, file='Results/ITS OTU_t.csv')

set.seed(88888) 
Bfn <- read.csv('Results/ITS OTU_t.csv', header=TRUE, row.names=1)
S1 <- specnumber(Bfn)
(raremax1 <- min(rowSums(Bfn)))    
newOTU <- rrarefy(Bfn, raremax1)
write.csv(newOTU, file='Results/ITS OTU_t_re.csv')

tnewOTU <- t(newOTU)
write.csv(tnewOTU, file='Results/ITS OTU_t_re_tAA.csv')

# Transformation
data <- read.csv("Results/ITS OTU_t_re_tAA.csv", row.names = 1, header = T) 
data <- decostand(data, MARGIN = 2, method = "total") 
write.csv(data, "Results/ITS OTU_RA.csv")  

data <- read.csv("Results/ITS OTU_RA.csv", row.names = 1, header = T) 
data <- data^0.5
write.csv(data,"Results/ITS OTU_rootsign.csv") 

#### Figure S2 Fungal community composition  ####
# Fungal Abundance at phylum level
OTU<-read.csv("Data/ITS OTU_RA & taxa.csv", header=TRUE)
Phylum = aggregate(x = OTU[, 9:228], by=list(OTU$Phylum), FUN='sum')
write.csv(Phylum, file = "Results/ITS Phyla.csv")

# plot
library(ggplot2)
library(ggalluvial)

taxon <- read.csv("Results/ITS Phyla top9 for plot.csv", header=TRUE)

taxon <- data.frame(subset(taxon, Treatment=="CK"))    #Choose different treatments
#taxon <- data.frame(subset(taxon, Treatment=="N"))    #Choose different treatments
#taxon <- data.frame(subset(taxon, Treatment=="NPK"))    #Choose different treatments
#taxon <- data.frame(subset(taxon, Treatment=="NPKM"))    #Choose different treatments
#taxon <- data.frame(subset(taxon, Treatment=="NPK1.5M"))    #Choose different treatments

Palette <- c("#F7FBFF", "#C6DBEF", "#6BAED6", "#4292C6", "#2171B5", "#08519C", 
             "#74C476", "#238B45", "#006D2C")
taxon$Treatment <- factor(taxon$Treatment, 
                          levels = c("CK", "N", "NPK", "NPKM", "NPK1.5M"), 
                          ordered=TRUE)
taxon$Time <- factor(taxon$Time, 
                     levels = c("T0", "T16", "T32", "T64", "T128", "T256", "T512", "T1024", "T2048", "T4096", "T8192"), 
                     ordered=TRUE)
taxon$Taxon <- factor(taxon$Taxon, 
                      levels = c("Ascomycota", "Other", "Basidiomycota", "Mortierellomycota", "Mucoromycota",
                                 "Glomeromycota", "Chytridiomycota", "Kickxellomycota", "Aphelidiomycota"), 
                      ordered=TRUE)
p2 <- ggplot(data = taxon, aes(x = Time, y = RA, alluvium = Taxon, stratum = Taxon)) +
  geom_alluvium(aes(fill = Taxon), alpha = 1, width = 0.2, show.legend = F) +  
  geom_stratum(aes(fill = Taxon), width = 0.2, show.legend = F) +
  scale_fill_manual(values = Palette) +
  scale_x_discrete(limits = taxon$Time) +
  ylab(label = "Relative Abundance") + 
  xlab(label = "") +
  labs(title = "Fungal community") +
  theme(plot.title = element_text(size = 20, colour = "black", hjust = 0, vjust = 0.5, face = "bold"))
p2




#### Figure 1c Fungal community similarity_time series ####
# Bray-Curtis dissimilarity
library(vegan)
library(ape)
library(ggplot2)
library(ggrepel)

data <- read.csv("Results/ITS OTU_rootsign.csv", header=TRUE, row.names=1)
data <- t(data)
otu.bc <- vegdist(data, 'bray')                 
write.csv(as.matrix(otu.bc), 'Results/ITS_Dissimilarity.csv')

# Bray-Curtis similarity
data <- read.csv("Results/ITS_Dissimilarity.csv", row.names = 1, header = T) 
data <- 1-data 
write.csv(data,"Results/ITS_Similarity.csv")


# Plot
data3 <- read.csv("Results/ITS_Similarity Time series for plot.csv", header=TRUE, row.names=1)
data3$Treatment <- factor(data3$Treatment, levels = c("CK", "N", "NPK", "NPKM", "NPK1.5M"))
c <- ggplot(data3, aes(x = Time, y = Similarity, col = Time)) +  
  scale_color_gradient(low='blue', high='red') +
  geom_point(size = 0.8) +
  ylab("Bray-Curtis similarity") +
  xlab("Air-drying and archiving time (h)") +
  theme(plot.title = element_text(size = 15, colour = "black", hjust = 0.5, face = "bold")) +
  scale_y_continuous(limits=c(0.4, 0.8), breaks = c( 0.4, 0.6, 0.8), label = c("0.4", "0.6", "0.8")) + 
  facet_wrap(~ Treatment, ncol = 5) +
  labs(title = "Fungal community", tag = "c") +
  theme(plot.title = element_text(size = 20, colour = "black", hjust = 0, vjust = 0.5, face = "bold"))
c




#### Figure S3def Fungal community similarity_time series & replication & treatment effects ####
library(ggplot2) 
library(patchwork) 

# Fungal community similarity_Air-drying and archiving time series
data <- read.csv("Results/ITS_Similarity Time series for plot.csv", header=TRUE, row.names=1)
data$Treatment <- factor(data$Treatment, c("CK", "N", "NPK", "NPKM", "NPK1.5M"), ordered=TRUE)
colour5 <- c("#F7FBFF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#2171B5")
d <- ggplot(data, aes(x = Treatment, y = Similarity, fill = Treatment)) + 
  geom_boxplot() +
  geom_jitter(size = 1) +    
  scale_fill_manual(values = colour5) +  
  scale_y_continuous(limits=c(0, 0.8),    
                     breaks = c(0.2, 0.4, 0.6), 
                     label = c("0.2", "0.4", "0.6")) + 
  ylab("Bray-Curtis similarity") +
  xlab("Treatment") +
  theme(axis.title.x = element_blank()) +   
  theme(legend.position = "none") +    
  labs(title = "Fungal community") + 
  labs(subtitle = "Air-drying and archiving time series", tag = "d") +
  theme(plot.subtitle = element_text(size = 12, colour = "black", hjust = 0, vjust = 0.5, face = "bold")) 
d

# Fungal community similarity_Cryopreservation (T = 0 h): replication effect
data <- read.csv("Results/ITS_Similarity replicate for plot.csv", header=TRUE, row.names=1)
data$Treatment <- factor(data$Treatment, levels = c("CK", "N", "NPK", "NPKM", "NPK1.5M"))
colour5 <- c("#F7FBFF", "#C6DBEF", "#9ECAE1", "#6BAED6", "#2171B5")
data1 <- subset(data, Time=="0")   
e <- ggplot(data1, aes(x = Treatment, y = Similarity, fill = Treatment)) +   
  geom_boxplot() +
  geom_jitter(size = 1) +  
  scale_fill_manual(values = colour5) +  
  scale_y_continuous(limits=c(0, 0.8),   
                     breaks = c(0.2, 0.4, 0.6), 
                     label = c("0.2", "0.4", "0.6")) + 
  ylab("Bray-Curtis similarity") +
  theme(axis.title.y = element_blank()) +   
  xlab("Treatment") +
  theme(legend.position = "none") +        
  labs(title = "Fungal community") +  
  labs(subtitle = "Cryopreservation (T = 0 h): replication effect", tag = "e") +
  theme(plot.subtitle = element_text(size = 12, colour = "black", hjust = 0, vjust = 0.5, face = "bold"))  
e

# Fungal community similarity_Cryopreservation (T = 0 h): treatment effect
data <- read.csv("Results/ITS_Similarity Treatment effects for plot.csv", header=TRUE, row.names=1)
data$Treatment <- factor(data$Treatment, c("N", "NPK", "NPKM", "NPK1.5M"), ordered=TRUE)
colour4 <- c("#C6DBEF", "#9ECAE1", "#6BAED6", "#2171B5")
f <- ggplot(data, aes(x = Treatment, y = Similarity, fill = Treatment)) +  
  geom_boxplot() +
  geom_jitter(size = 1) +   
  scale_fill_manual(values = colour4) +  
  scale_y_continuous(limits=c(0, 0.8),   
                     breaks = c(0.2, 0.4, 0.6), 
                     label = c("0.2", "0.4", "0.6")) + 
  ylab("Bray-Curtis similarity") +
  theme(axis.title.y = element_blank()) +   
  xlab("Treatment") +
  theme(axis.title.x = element_blank()) +    
  theme(legend.position = "none") +    
  labs(title = "Fungal community") + 
  labs(subtitle = "Cryopreservation (T = 0 h): treatment effect", tag = "f") +
  theme(plot.subtitle = element_text(size = 12, colour = "black", hjust = 0, vjust = 0.5, face = "bold")) 
f

# wilcox.test: to compare similarity of time series & treatment effects
data <- read.csv("Results/ITS_Similarity for test.csv", header=TRUE, row.names=1)
wilcox.test(data$TSB, data$TEB) 
wilcox.test(data$TSC, data$TEC) 
wilcox.test(data$TSD, data$TED) 
wilcox.test(data$TSE, data$TEE) 








#### Figure 2b Fungal community structure_Air-drying and archiving time series ####
# NMDS
library(vegan)
OTU <- read.csv("Results/ITS OTU_rootsign.csv", header=TRUE, row.names=1)
tOTU <- t(OTU)

set.seed(8888)  
tOTU_NMDS1 <- metaMDS(tOTU, distance = 'bray', k = 2, autotransform = FALSE, trymax = 2000)  
tOTU_NMDS <- metaMDS(tOTU, distance = "bray", k = 2, autotransform = FALSE, previous.best = tOTU_NMDS1)

Stress <- tOTU_NMDS$stress 
Stress

NMDS_scores <- scores(tOTU_NMDS, choices = c(1,2))  
write.csv(NMDS_scores, file="Results/ITS NMDS Scores_Air-drying and archiving time series.csv")

# PERMANOVA/ADONIS
library(vegan)
otu <- read.csv("Results/ITS OTU_rootsign.csv", header=TRUE, row.names=1)
otu <- data.frame(t(otu)) 

group <- read.csv('Data/groups.csv', header=TRUE, row.names=1)

set.seed(8888)
adonis_result <- adonis(otu~Treatment, group, distance = 'bray', permutations = 2000)
summary(adonis_result)

R2 <- adonis_result$aov.tab$R2[1]
R2

P <- adonis_result$aov.tab$"Pr(>F)"[1]
P 

otuput <- data.frame(adonis_result$aov.tab, check.names = FALSE, stringsAsFactors = FALSE)
otuput <- cbind(rownames(otuput), otuput)
names(otuput) <- c('', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
write.table(otuput, file = 'Results/ITS PERMANOVA_Air-drying and archiving time series.csv', row.names = FALSE, sep = ",")

# Plot
data <- read.csv("Results/ITS NMDS Scores_Air-drying and archiving time series for plot.csv", header=TRUE)

data$Treatment <- factor(data$Treatment, c("CK", "N", "NPK", "NPKM", "NPK1.5M"), ordered=TRUE)

library(gridExtra) 
label1 <- c(paste("Air-drying and archiving time series"))
label2 <- c(paste("Stress ==", round(Stress, 3),  sep = ""))
label3 <- c(paste("PERMANOVA"),
            paste("R2 ==", round(R2, 3),  sep = ""),  
            paste("P < 0.001 ", sep = ""))

thm1 <- ttheme_minimal(core = list(fg_params = list(parse = TRUE, col = "black", hjust = 0, fontsize =10, x = 0.4)))
thm2 <- ttheme_minimal(core = list(fg_params = list(parse = TRUE, col = "black", hjust = 0,  x = 0.1)))

tbgrob1 <- tableGrob(as.matrix(label1), theme = thm1)
tbgrob2 <- tableGrob(as.matrix(label2), theme = thm2)
tbgrob3 <- tableGrob(as.matrix(label3), theme = thm2)

library(gridExtra) 
library(ggplot2) 
b <- ggplot(data = data) +
  geom_point(aes(x = NMDS1, y = NMDS2, col = Time, shape = Treatment), size=3) +
  scale_colour_gradient(low = "blue", high = "red") +
  scale_shape_manual(values = c(1, 3, 11, 4, 0)) +  
  scale_x_continuous(limits = c(-0.80, 0.80)) +  
  scale_y_continuous(limits = c(-0.60, 0.60)) + 
  xlab("NMDS1") +
  ylab("NMDS2") +
  labs(title = "Fungal community", tag = "b") +
  theme(plot.title = element_text(size = 20, colour = "black", hjust = 0, vjust = 0.5, face = "bold")) +
  annotation_custom(tbgrob1, xmin = -2.2, ymax = 1.8) +    
  annotation_custom(tbgrob2, xmin = -1.2, ymax = -0.4) +     
  annotation_custom(tbgrob3, xmin = 0.2, ymax = -0.3)     
b





#### Figure 2d Fungal community structure_Cryopreservation (T=0h) ####
library(vegan)
library(ggplot2)  
library(dplyr) 
library(ggplot2)
library(ggrepel)  #geom_text_repel
library(dplyr)  #filter
library(cowplot) #拼图
library(gridExtra) 
# NMDS
data <- read.csv("Results/ITS OTU_rootsign_sort by time.csv", header=TRUE, row.names=1)
#T=0h
data0 <- data[,c(1:20)] 
tOTU <- data.frame(t(data0))

set.seed(8888)  
tOTU_NMDS1 <- metaMDS(tOTU, distance = 'bray', k = 2, autotransform = FALSE, trymax = 2000)  
tOTU_NMDS <- metaMDS(tOTU, distance = "bray", k = 2, autotransform = FALSE, previous.best = tOTU_NMDS1)

Stress <- tOTU_NMDS$stress 
Stress

NMDS_scores=scores(tOTU_NMDS, choices = c(1,2))  
write.csv(NMDS_scores, file="Results/ITS NMDS Scores_Cryopreservation (T=0h).csv")

# PERMANOVA/ADONIS
data <- read.csv("Results/ITS OTU_rootsign_sort by time.csv", header=TRUE, row.names=1)
#T=0h
data0 <- data[,c(1:20)] 
tOTU <- data.frame(t(data0))

group_all <- read.csv('Data/groups.csv', header=TRUE, row.names=1)
group <- dplyr::filter(group_all, Time == "0")

set.seed(8888)
adonis_result <- adonis(tOTU~Treatment, group, distance = 'bray', permutations = 2000)
summary(adonis_result)

R2 <- adonis_result$aov.tab$R2[1]
R2

P <- adonis_result$aov.tab$"Pr(>F)"[1]
P 

otuput <- data.frame(adonis_result$aov.tab, check.names = FALSE, stringsAsFactors = FALSE)
otuput <- cbind(rownames(otuput), otuput)
names(otuput) <- c('', 'Df', 'Sums of squares', 'Mean squares', 'F.Model', 'Variation (R2)', 'Pr (>F)')
write.table(otuput, file = 'Results/ITS PERMANOVA_Cryopreservation (T=0h).csv', row.names = FALSE, sep = ",")

# Plot
data <- read.csv("Results/ITS NMDS Scores_Cryopreservation (T=0h) for plot.csv", header=TRUE)

data$Treatment <- factor(data$Treatment, c("CK", "N", "NPK", "NPKM", "NPK1.5M"), ordered=TRUE)

library(gridExtra) 
label1 <- c(paste("Time series"))
label2 <- c(paste("Stress ==", round(Stress, 3),  sep = ""))
label3 <- c(paste("PERMANOVA"),
            paste("R2 ==", round(R2, 3),  sep = ""),  
            paste("P < 0.001 ", sep = ""))

thm1 <- ttheme_minimal(core = list(fg_params = list(parse = TRUE, col = "black", hjust = 0, fontsize =10, x = 0.2)))
thm2 <- ttheme_minimal(core = list(fg_params = list(parse = TRUE, col = "black", hjust = 0,  x = 0.1)))

tbgrob1 <- tableGrob(as.matrix(label1), theme = thm1)
tbgrob2 <- tableGrob(as.matrix(label2), theme = thm2)
tbgrob3 <- tableGrob(as.matrix(label3), theme = thm2)

c <- ggplot(data = data) +
  geom_point(aes(x = NMDS1, y = NMDS2, shape = Treatment), col = "blue", size = 3) +
  scale_shape_manual(values = c(1, 3, 11, 4, 0)) +  
  scale_x_continuous(limits = c(-0.80, 0.80)) +  
  scale_y_continuous(limits = c(-0.60, 0.60)) + 
  xlab("NMDS1") +
  ylab("NMDS2") +
  labs(title = "Fungal community", tag = "d") +
  theme(plot.title = element_text(size = 20, colour = "black", hjust = 0, vjust = 0.5, face = "bold")) +
  annotation_custom(tbgrob1, xmin = -2.2, ymax = 1.8) +    
  annotation_custom(tbgrob2, xmin = -1.2, ymax = -0.4) +     
  annotation_custom(tbgrob3, xmin = 0.2, ymax = -0.3)   
c






#### Figure 2e Prokaryotic and fungal community_Mantel test ####
#Prokaryotic community_Mantel test
library(vegan)
comm1 <- t(read.csv("Results/16s OTU_rootsign_sort by time.csv", header=TRUE, row.names=1))

group.mantel <- function(a1, a2, n){
  coresult.out <- data.frame()
  g1 <- mantel(1-vegdist(comm1[a1:(a1+n-1),], method = 'bray'),
               1-vegdist(comm1[a2:(a2+n-1),], method = 'bray'), 
               method="spearman", permutations = 2000)
  coresult.out[1,1] <- a1
  coresult.out[1,2] <- (a1+n-1)
  coresult.out[1,3] <- a2
  coresult.out[1,4] <- (a2+n-1)
  coresult.out[1,5] <- g1$statistic
  coresult.out[1,6] <- g1$signif
  write.table(data.frame(coresult.out), 'Results/16s Mantel Test spearman_Similarity.csv', append= T, sep=',', row.names=F, col.names=F)
}

group.mantel(a1=1, a2=1, n=20)               
group.mantel(a1=1, a2=21, n=20)
group.mantel(a1=1, a2=41, n=20)
group.mantel(a1=1, a2=61, n=20)
group.mantel(a1=1, a2=81, n=20)
group.mantel(a1=1, a2=101, n=20)
group.mantel(a1=1, a2=121, n=20)
group.mantel(a1=1, a2=141, n=20)
group.mantel(a1=1, a2=161, n=20)
group.mantel(a1=1, a2=181, n=20)
group.mantel(a1=1, a2=201, n=20)

group.mantel(a1=21, a2=1, n=20)
group.mantel(a1=21, a2=21, n=20)
group.mantel(a1=21, a2=41, n=20)
group.mantel(a1=21, a2=61, n=20)
group.mantel(a1=21, a2=81, n=20)
group.mantel(a1=21, a2=101, n=20)
group.mantel(a1=21, a2=121, n=20)
group.mantel(a1=21, a2=141, n=20)
group.mantel(a1=21, a2=161, n=20)
group.mantel(a1=21, a2=181, n=20)
group.mantel(a1=21, a2=201, n=20)

group.mantel(a1=41, a2=1, n=20)
group.mantel(a1=41, a2=21, n=20)
group.mantel(a1=41, a2=41, n=20)
group.mantel(a1=41, a2=61, n=20)
group.mantel(a1=41, a2=81, n=20)
group.mantel(a1=41, a2=101, n=20)
group.mantel(a1=41, a2=121, n=20)
group.mantel(a1=41, a2=141, n=20)
group.mantel(a1=41, a2=161, n=20)
group.mantel(a1=41, a2=181, n=20)
group.mantel(a1=41, a2=201, n=20)

group.mantel(a1=61, a2=1, n=20)
group.mantel(a1=61, a2=21, n=20)
group.mantel(a1=61, a2=41, n=20)
group.mantel(a1=61, a2=61, n=20)
group.mantel(a1=61, a2=81, n=20)
group.mantel(a1=61, a2=101, n=20)
group.mantel(a1=61, a2=121, n=20)
group.mantel(a1=61, a2=141, n=20)
group.mantel(a1=61, a2=161, n=20)
group.mantel(a1=61, a2=181, n=20)
group.mantel(a1=61, a2=201, n=20)

group.mantel(a1=81, a2=1, n=20)
group.mantel(a1=81, a2=21, n=20)
group.mantel(a1=81, a2=41, n=20)
group.mantel(a1=81, a2=61, n=20)
group.mantel(a1=81, a2=81, n=20)
group.mantel(a1=81, a2=101, n=20)
group.mantel(a1=81, a2=121, n=20)
group.mantel(a1=81, a2=141, n=20)
group.mantel(a1=81, a2=161, n=20)
group.mantel(a1=81, a2=181, n=20)
group.mantel(a1=81, a2=201, n=20)

group.mantel(a1=101, a2=1, n=20)
group.mantel(a1=101, a2=21, n=20)
group.mantel(a1=101, a2=41, n=20)
group.mantel(a1=101, a2=61, n=20)
group.mantel(a1=101, a2=81, n=20)
group.mantel(a1=101, a2=101, n=20)
group.mantel(a1=101, a2=121, n=20)
group.mantel(a1=101, a2=141, n=20)
group.mantel(a1=101, a2=161, n=20)
group.mantel(a1=101, a2=181, n=20)
group.mantel(a1=101, a2=201, n=20)

group.mantel(a1=121, a2=1, n=20)
group.mantel(a1=121, a2=21, n=20)
group.mantel(a1=121, a2=41, n=20)
group.mantel(a1=121, a2=61, n=20)
group.mantel(a1=121, a2=81, n=20)
group.mantel(a1=121, a2=101, n=20)
group.mantel(a1=121, a2=121, n=20)
group.mantel(a1=121, a2=141, n=20)
group.mantel(a1=121, a2=161, n=20)
group.mantel(a1=121, a2=181, n=20)
group.mantel(a1=121, a2=201, n=20)

group.mantel(a1=141, a2=1, n=20)
group.mantel(a1=141, a2=21, n=20)
group.mantel(a1=141, a2=41, n=20)
group.mantel(a1=141, a2=61, n=20)
group.mantel(a1=141, a2=81, n=20)
group.mantel(a1=141, a2=101, n=20)
group.mantel(a1=141, a2=121, n=20)
group.mantel(a1=141, a2=141, n=20)
group.mantel(a1=141, a2=161, n=20)
group.mantel(a1=141, a2=181, n=20)
group.mantel(a1=141, a2=201, n=20)

group.mantel(a1=161, a2=1, n=20)
group.mantel(a1=161, a2=21, n=20)
group.mantel(a1=161, a2=41, n=20)
group.mantel(a1=161, a2=61, n=20)
group.mantel(a1=161, a2=81, n=20)
group.mantel(a1=161, a2=101, n=20)
group.mantel(a1=161, a2=121, n=20)
group.mantel(a1=161, a2=141, n=20)
group.mantel(a1=161, a2=161, n=20)
group.mantel(a1=161, a2=181, n=20)
group.mantel(a1=161, a2=201, n=20)

group.mantel(a1=181, a2=1, n=20)
group.mantel(a1=181, a2=21, n=20)
group.mantel(a1=181, a2=41, n=20)
group.mantel(a1=181, a2=61, n=20)
group.mantel(a1=181, a2=81, n=20)
group.mantel(a1=181, a2=101, n=20)
group.mantel(a1=181, a2=121, n=20)
group.mantel(a1=181, a2=141, n=20)
group.mantel(a1=181, a2=161, n=20)
group.mantel(a1=181, a2=181, n=20)
group.mantel(a1=181, a2=201, n=20)

group.mantel(a1=201, a2=1, n=20)
group.mantel(a1=201, a2=21, n=20)
group.mantel(a1=201, a2=41, n=20)
group.mantel(a1=201, a2=61, n=20)
group.mantel(a1=201, a2=81, n=20)
group.mantel(a1=201, a2=101, n=20)
group.mantel(a1=201, a2=121, n=20)
group.mantel(a1=201, a2=141, n=20)
group.mantel(a1=201, a2=161, n=20)
group.mantel(a1=201, a2=181, n=20)
group.mantel(a1=201, a2=201, n=20)

#Fungal community_Mantel test
library(vegan)
comm1 <- t(read.csv("Results/ITS OTU_rootsign_sort by time.csv", header=TRUE, row.names=1))

group.mantel <- function(a1, a2, n){
  coresult.out <- data.frame()
  g1 <- mantel(1-vegdist(comm1[a1:(a1+n-1),], method = 'bray'),
               1-vegdist(comm1[a2:(a2+n-1),], method = 'bray'), 
               method="spearman", permutations = 2000)
  coresult.out[1,1] <- a1
  coresult.out[1,2] <- (a1+n-1)
  coresult.out[1,3] <- a2
  coresult.out[1,4] <- (a2+n-1)
  coresult.out[1,5] <- g1$statistic
  coresult.out[1,6] <- g1$signif
  write.table(data.frame(coresult.out), 'Results/ITS Mantel Test spearman_Similarity.csv', append= T, sep=',', row.names=F, col.names=F)
}

group.mantel(a1=1, a2=1, n=20)                
group.mantel(a1=1, a2=21, n=20)
group.mantel(a1=1, a2=41, n=20)
group.mantel(a1=1, a2=61, n=20)
group.mantel(a1=1, a2=81, n=20)
group.mantel(a1=1, a2=101, n=20)
group.mantel(a1=1, a2=121, n=20)
group.mantel(a1=1, a2=141, n=20)
group.mantel(a1=1, a2=161, n=20)
group.mantel(a1=1, a2=181, n=20)
group.mantel(a1=1, a2=201, n=20)

group.mantel(a1=21, a2=1, n=20)
group.mantel(a1=21, a2=21, n=20)
group.mantel(a1=21, a2=41, n=20)
group.mantel(a1=21, a2=61, n=20)
group.mantel(a1=21, a2=81, n=20)
group.mantel(a1=21, a2=101, n=20)
group.mantel(a1=21, a2=121, n=20)
group.mantel(a1=21, a2=141, n=20)
group.mantel(a1=21, a2=161, n=20)
group.mantel(a1=21, a2=181, n=20)
group.mantel(a1=21, a2=201, n=20)

group.mantel(a1=41, a2=1, n=20)
group.mantel(a1=41, a2=21, n=20)
group.mantel(a1=41, a2=41, n=20)
group.mantel(a1=41, a2=61, n=20)
group.mantel(a1=41, a2=81, n=20)
group.mantel(a1=41, a2=101, n=20)
group.mantel(a1=41, a2=121, n=20)
group.mantel(a1=41, a2=141, n=20)
group.mantel(a1=41, a2=161, n=20)
group.mantel(a1=41, a2=181, n=20)
group.mantel(a1=41, a2=201, n=20)

group.mantel(a1=61, a2=1, n=20)
group.mantel(a1=61, a2=21, n=20)
group.mantel(a1=61, a2=41, n=20)
group.mantel(a1=61, a2=61, n=20)
group.mantel(a1=61, a2=81, n=20)
group.mantel(a1=61, a2=101, n=20)
group.mantel(a1=61, a2=121, n=20)
group.mantel(a1=61, a2=141, n=20)
group.mantel(a1=61, a2=161, n=20)
group.mantel(a1=61, a2=181, n=20)
group.mantel(a1=61, a2=201, n=20)

group.mantel(a1=81, a2=1, n=20)
group.mantel(a1=81, a2=21, n=20)
group.mantel(a1=81, a2=41, n=20)
group.mantel(a1=81, a2=61, n=20)
group.mantel(a1=81, a2=81, n=20)
group.mantel(a1=81, a2=101, n=20)
group.mantel(a1=81, a2=121, n=20)
group.mantel(a1=81, a2=141, n=20)
group.mantel(a1=81, a2=161, n=20)
group.mantel(a1=81, a2=181, n=20)
group.mantel(a1=81, a2=201, n=20)

group.mantel(a1=101, a2=1, n=20)
group.mantel(a1=101, a2=21, n=20)
group.mantel(a1=101, a2=41, n=20)
group.mantel(a1=101, a2=61, n=20)
group.mantel(a1=101, a2=81, n=20)
group.mantel(a1=101, a2=101, n=20)
group.mantel(a1=101, a2=121, n=20)
group.mantel(a1=101, a2=141, n=20)
group.mantel(a1=101, a2=161, n=20)
group.mantel(a1=101, a2=181, n=20)
group.mantel(a1=101, a2=201, n=20)

group.mantel(a1=121, a2=1, n=20)
group.mantel(a1=121, a2=21, n=20)
group.mantel(a1=121, a2=41, n=20)
group.mantel(a1=121, a2=61, n=20)
group.mantel(a1=121, a2=81, n=20)
group.mantel(a1=121, a2=101, n=20)
group.mantel(a1=121, a2=121, n=20)
group.mantel(a1=121, a2=141, n=20)
group.mantel(a1=121, a2=161, n=20)
group.mantel(a1=121, a2=181, n=20)
group.mantel(a1=121, a2=201, n=20)

group.mantel(a1=141, a2=1, n=20)
group.mantel(a1=141, a2=21, n=20)
group.mantel(a1=141, a2=41, n=20)
group.mantel(a1=141, a2=61, n=20)
group.mantel(a1=141, a2=81, n=20)
group.mantel(a1=141, a2=101, n=20)
group.mantel(a1=141, a2=121, n=20)
group.mantel(a1=141, a2=141, n=20)
group.mantel(a1=141, a2=161, n=20)
group.mantel(a1=141, a2=181, n=20)
group.mantel(a1=141, a2=201, n=20)

group.mantel(a1=161, a2=1, n=20)
group.mantel(a1=161, a2=21, n=20)
group.mantel(a1=161, a2=41, n=20)
group.mantel(a1=161, a2=61, n=20)
group.mantel(a1=161, a2=81, n=20)
group.mantel(a1=161, a2=101, n=20)
group.mantel(a1=161, a2=121, n=20)
group.mantel(a1=161, a2=141, n=20)
group.mantel(a1=161, a2=161, n=20)
group.mantel(a1=161, a2=181, n=20)
group.mantel(a1=161, a2=201, n=20)

group.mantel(a1=181, a2=1, n=20)
group.mantel(a1=181, a2=21, n=20)
group.mantel(a1=181, a2=41, n=20)
group.mantel(a1=181, a2=61, n=20)
group.mantel(a1=181, a2=81, n=20)
group.mantel(a1=181, a2=101, n=20)
group.mantel(a1=181, a2=121, n=20)
group.mantel(a1=181, a2=141, n=20)
group.mantel(a1=181, a2=161, n=20)
group.mantel(a1=181, a2=181, n=20)
group.mantel(a1=181, a2=201, n=20)

group.mantel(a1=201, a2=1, n=20)
group.mantel(a1=201, a2=21, n=20)
group.mantel(a1=201, a2=41, n=20)
group.mantel(a1=201, a2=61, n=20)
group.mantel(a1=201, a2=81, n=20)
group.mantel(a1=201, a2=101, n=20)
group.mantel(a1=201, a2=121, n=20)
group.mantel(a1=201, a2=141, n=20)
group.mantel(a1=201, a2=161, n=20)
group.mantel(a1=201, a2=181, n=20)
group.mantel(a1=201, a2=201, n=20)



# plot
library(ggcor)
library(ggplot2)

upper <- read.csv("Results/ITS Mantel Test spearman matarix.csv", header = TRUE, row.names = 1)      
lower <- read.csv("Results/16s Mantel Test spearman matarix.csv", header = TRUE, row.names = 1)

cor_df <- cor_tbl(upper, extra.mat = list(lower = lower))

quickcor(cor_df) +
  geom_circle2(data = get_data(.col.id + .row.id > 12)) +
  geom_circle2(aes(fill = lower, r0 = lower), 
               data = get_data(.col.id + .row.id < 12)) +
  geom_number(aes(num = r), nsmall = 3, colour = "white", size = 3,
              data = get_data(.col.id + .row.id > 12)) +
  geom_number(aes(num = lower), nsmall = 3, colour = "white", size = 3,
              data = get_data(.col.id + .row.id < 12)) +
  geom_text(aes(x = .col.id, y = .row.id), label = names(upper), 
            angle = 45, data = get_data(.col.id + .row.id == 12),
            inherit.aes = FALSE) +
  scale_fill_gradient2n(limits = c(-1, 1), colours = rev(red_blue())) +
  guides(fill = guide_colorbar("")) +
  remove_all_axis()


