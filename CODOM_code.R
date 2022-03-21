library(phyloseq)
library(agricolae)
library(tidyverse)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(emmeans)
library(car)
library(MuMIn)
library(agridat)
library(DHARMa)
library(viridis)
library(hillR)
library(GGally)
library(lavaan)
library(compositions)
library(vegan)
library(corrplot)
library(ggplot2)
library(ape)

#### Set working directory
setwd("~/Desktop/Codom/BF")

#### Import into R
## ASV table
otu_taxa <- read_csv("OTU_FinalFinal.csv")
head(otu_taxa)

## Mapping file
mapping<-read_csv("META_Final.csv")
head(mapping)

### reshape data
otu_long <- otu_taxa %>% pivot_longer(cols=c(2:41), names_to = "ID", values_to = 'abund') %>% full_join(mapping,.) %>% mutate(pres=if_else(abund>0,1,0))

## ALPHA DIVERSITY ##
##### Calculate diversity metrics for each trophic group, then combine into one dataset ###############################################################
## BF
otu_BF <- otu_long %>% filter(Trophic=='BF') %>% select(ID, Site, Plot, Plant_removal, N_addition, Trophic, OTU_ID, abund) %>%                      
  group_by(ID, Site, Plot, Plant_removal, N_addition) %>% pivot_wider(names_from = OTU_ID, values_from = abund)                      

otu_BF_div <- otu_BF[1:6]

otu_BF_div$rich <- hill_taxa(otu_BF[7:74], q=0)                                                                                               
otu_BF_div$shan <- hill_taxa(otu_BF[7:74], q=1)                                                                                               
otu_BF_div$simp <- hill_taxa(otu_BF[7:74], q=2)
## FF
otu_FF <- otu_long %>% filter(Trophic=='FF') %>% select(ID, Site, Plot, Plant_removal, N_addition, Trophic, OTU_ID, abund) %>%                      
  group_by(ID, Site, Plot, Plant_removal, N_addition) %>% pivot_wider(names_from = OTU_ID, values_from = abund)                      

otu_FF_div <- otu_FF[1:6]

otu_FF_div$rich <- hill_taxa(otu_FF[7:42], q=0)                                                                                               
otu_FF_div$shan <- hill_taxa(otu_FF[7:42], q=1)                                                                                               
otu_FF_div$simp <- hill_taxa(otu_FF[7:42], q=2)
## OM
otu_OM <- otu_long %>% filter(Trophic=='OM') %>% select(ID, Site, Plot, Plant_removal, N_addition, Trophic, OTU_ID, abund) %>%                      
  group_by(ID, Site, Plot, Plant_removal, N_addition) %>% pivot_wider(names_from = OTU_ID, values_from = abund)                      

otu_OM_div <- otu_OM[1:6]

otu_OM_div$rich <- hill_taxa(otu_OM[7:28], q=0)                                                                                               
otu_OM_div$shan <- hill_taxa(otu_OM[7:28], q=1)                                                                                               
otu_OM_div$simp <- hill_taxa(otu_OM[7:28], q=2)
## PP
otu_PP <- otu_long %>% filter(Trophic=='PP') %>% select(ID, Site, Plot, Plant_removal, N_addition, Trophic, OTU_ID, abund) %>%                      
  group_by(ID, Site, Plot, Plant_removal, N_addition) %>% pivot_wider(names_from = OTU_ID, values_from = abund)                      

otu_PP_div <- otu_PP[1:6]

otu_PP_div$rich <- hill_taxa(otu_PP[7:48], q=0)                                                                                               
otu_PP_div$shan <- hill_taxa(otu_PP[7:48], q=1)                                                                                               
otu_PP_div$simp <- hill_taxa(otu_PP[7:48], q=2)
## RA
otu_RA <- otu_long %>% filter(Trophic=='RA') %>% select(ID, Site, Plot, Plant_removal, N_addition, Trophic, OTU_ID, abund) %>%                      
  group_by(ID, Site, Plot, Plant_removal, N_addition) %>% pivot_wider(names_from = OTU_ID, values_from = abund)                      

otu_RA_div <- otu_RA[1:6]

otu_RA_div$rich <- hill_taxa(otu_RA[7:38], q=0)                                                                                               
otu_RA_div$shan <- hill_taxa(otu_RA[7:38], q=1)                                                                                               
otu_RA_div$simp <- hill_taxa(otu_RA[7:38], q=2)
## combine into one sheet
div1 <- rbind(otu_BF_div,otu_FF_div)
div2 <- rbind(div1,otu_OM_div)
div3 <- rbind(div2,otu_PP_div)
div <- rbind(div3,otu_RA_div)
div$simp[!is.finite(div$simp)] <- 0
head(div)

## Check if trophic groups are significant ##
troph_rich <- lmer(rich ~ Plant_removal * N_addition * Trophic + (1|Site), data=div)
anova(troph_rich)

troph_shan <- lmer(shan ~ Plant_removal * N_addition * Trophic + (1|Site), data=div)
anova(troph_shan)

troph_simp <- lmer(simp ~ Plant_removal * N_addition * Trophic + (1|Site), data=div)
anova(troph_simp)

## Significance of each metric on each trophic group, test each metric for each trophic group
## anova with all individual trophic groups SHANNON
lm1 <- lmer(shan ~ Plant_removal * N_addition + (1|Site), data=div %>% filter(Trophic == 'BF'))
anova(lm1)

lm2 <- lmer(shan ~ Plant_removal * N_addition + (1|Site), data=div %>% filter(Trophic == 'FF'))
anova(lm2)

lm3 <- lmer(shan ~ Plant_removal * N_addition + (1|Site), data=div %>% filter(Trophic == 'PP'))
anova(lm3)

lm4 <- lmer(shan ~ Plant_removal * N_addition + (1|Site), data=div %>% filter(Trophic == 'OM'))
anova(lm4)

lm5 <- lmer(shan ~ Plant_removal * N_addition + (1|Site), data=div %>% filter(Trophic == 'RA'))
anova(lm5)

## anova with all individual trophic groups SIMPSON
lm6 <- lmer(simp ~ Plant_removal * N_addition + (1|Site), data=div %>% filter(Trophic == 'BF'))
anova(lm6)

lm7 <- lmer(simp ~ Plant_removal * N_addition + (1|Site), data=div %>% filter(Trophic == 'FF'))
anova(lm7)

lm8 <- lmer(simp ~ Plant_removal * N_addition + (1|Site), data=div %>% filter(Trophic == 'PP'))
anova(lm8)

lm9 <- lmer(simp ~ Plant_removal * N_addition + (1|Site), data=div %>% filter(Trophic == 'OM'))
anova(lm9)

lm10 <- lmer(simp ~ Plant_removal * N_addition + (1|Site), data=div %>% filter(Trophic == 'RA'))
anova(lm10)

## anova with all individual trophic groups RICHNESS
lm11 <- lmer(rich ~ Plant_removal * N_addition + (1|Site), data=div %>% filter(Trophic == 'BF'))
anova(lm11)

lm12 <- lmer(rich ~ Plant_removal * N_addition + (1|Site), data=div %>% filter(Trophic == 'FF'))
anova(lm12)

lm13 <- lmer(rich ~ Plant_removal * N_addition + (1|Site), data=div %>% filter(Trophic == 'PP'))
anova(lm13)

lm14 <- lmer(rich ~ Plant_removal * N_addition + (1|Site), data=div %>% filter(Trophic == 'OM'))
anova(lm14)

lm15 <- lmer(rich ~ Plant_removal * N_addition + (1|Site), data=div %>% filter(Trophic == 'RA'))
anova(lm15)

## visualize alpha diversity metrics ##
## I like to re-order my variables first to make the graph easier to read
# Turn your 'treatment' column into a character vector
div$N_addition <- as.character(div$N_addition)
# Then turn it back into a factor with the levels in the correct order
div$N_addition <- factor(div$N_addition, levels=c("X","N"))

# Turn your 'treatment' column into a character vector
div$Plant_removal <- as.character(div$Plant_removal)
# Then turn it back into a factor with the levels in the correct order
div$Plant_removal <- factor(div$Plant_removal, levels=c("X","G","D"))

## I also want trophic groups in a certain order
# Turn your 'treatment' column into a character vector
div$Trophic <- as.character(div$Trophic)
# Then turn it back into a factor with the levels in the correct order
div$Trophic <- factor(div$Trophic, levels=c("BF","PP","OM","RA","FF"))

#### Plot out each trophic group by each measure
ggplot(div, aes(x=N_addition, y=shan, fill=Plant_removal)) + geom_boxplot()  + facet_wrap(~Trophic, nrow=1) + theme_bw()+ xlab("Nitrogen Addition")+ylab("Shannon Diversity Index")+guides(fill=guide_legend(title="Plant Removal"))+scale_fill_manual("legend", values = c("X" = "grey", "G" = "blue", "D" = "red")) +theme( axis.text = element_text( size = 18 ),
                                                                                                                                                                                                                                                                                                                              axis.text.x = element_text( size = 20 ),
                                                                                                                                                                                                                                                                                                                              axis.title = element_text( size = 16, face = "bold" ),
                                                                                                                                                                                                                                                                                                                              strip.text = element_text(size = 20))
ggplot(div, aes(x=N_addition, y=simp, fill=Plant_removal)) + geom_boxplot()  + facet_wrap(~Trophic, nrow=1) + theme_bw()+ xlab("Nitrogen Addition")+ylab("Simpson Index")+guides(fill=guide_legend(title="Plant Removal"))+scale_fill_manual("legend", values = c("X" = "grey", "G" = "blue", "D" = "red")) +theme( axis.text = element_text( size = 18 ),
                                                                                                                                                                                                                                                                                                                              axis.text.x = element_text( size = 20 ),
                                                                                                                                                                                                                                                                                                                              axis.title = element_text( size = 16, face = "bold" ),
                                                                                                                                                                                                                                                                                                                              strip.text = element_text(size = 20))

ggplot(div, aes(x=N_addition, y=rich, fill=Plant_removal)) + geom_boxplot()  + facet_wrap(~Trophic, nrow=1) + theme_bw()+ xlab("Nitrogen Addition")+ylab("Richness")+guides(fill=guide_legend(title="Plant Removal"))+scale_fill_manual("legend", values = c("X" = "grey", "G" = "blue", "D" = "red")) +theme( axis.text = element_text( size = 18 ),
                                                                                                                                                                                                                                                                                                                              axis.text.x = element_text( size = 20 ),
                                                                                                                                                                                                                                                                                                                              axis.title = element_text( size = 16, face = "bold" ),
                                                                                                                                                                                                                                                                                                                              strip.text = element_text(size = 20))



#### Plot out entire community by each metric
ggplot(div, aes(x=N_addition, y=shan, fill=Plant_removal)) + geom_boxplot()  +  theme_bw()+ xlab("Nitrogen Addition")+ylab("Shannon Diversity Index")+guides(fill=guide_legend(title="Plant Removal"))+scale_fill_manual("legend", values = c("X" = "grey", "G" = "blue", "D" = "red")) +theme( axis.text = element_text( size = 18 ),
                                                                                                                                                                                                                                                                                                                              axis.text.x = element_text( size = 20 ),
                                                                                                                                                                                                                                                                                                                              axis.title = element_text( size = 16, face = "bold" ),
                                                                                                                                                                                                                                                                                                                              strip.text = element_text(size = 20))

ggplot(div, aes(x=N_addition, y=simp, fill=Plant_removal)) + geom_boxplot()   + theme_bw()+ xlab("Nitrogen Addition")+ylab("Simpson Index")+guides(fill=guide_legend(title="Plant Removal"))+scale_fill_manual("legend", values = c("X" = "grey", "G" = "blue", "D" = "red")) +theme( axis.text = element_text( size = 18 ),
                                                                                                                                                                                                                                                                                                                    axis.text.x = element_text( size = 20 ),
                                                                                                                                                                                                                                                                                                                    axis.title = element_text( size = 16, face = "bold" ),
                                                                                                                                                                                                                                                                                                                    strip.text = element_text(size = 20))
ggplot(div, aes(x=N_addition, y=rich, fill=Plant_removal)) + geom_boxplot()  + theme_bw()+ xlab("Nitrogen Addition")+ylab("Richness")+guides(fill=guide_legend(title="Plant Removal"))+scale_fill_manual("legend", values = c("X" = "grey", "G" = "blue", "D" = "red")) +theme( axis.text = element_text( size = 18 ),
                                                                                                                                                                                                                                                                                                               axis.text.x = element_text( size = 20 ),
                                                                                                                                                                                                                                                                                                               axis.title = element_text( size = 16, face = "bold" ),
                                                                                                                                                                                                                                                                                                               strip.text = element_text(size = 20))




########### Abundance Analyses of Morphological counts ########
## Input morphological data

dat<-read_csv("AbundMorphCounts.csv")
head(dat)

### Abundance of each trophic group by treatments
BF<-lmer(BF ~ N_addition*Plant_removal +(1|Site),data=dat)
anova(BF)

FF<-lmer(FF ~ N_addition*Plant_removal +(1|Site),data=dat)
anova(FF)

PP<-lmer(PP ~ N_addition*Plant_removal +(1|Site),data=dat)
anova(PP)

RA<-lmer(RA ~ N_addition*Plant_removal +(1|Site),data=dat)
anova(RA)

OM<-lmer(OM ~ N_addition*Plant_removal +(1|Site),data=dat)
anova(OM)

## Total abundance by treatments
SUM<-lmer(Sum ~ N_addition*Plant_removal +(1|Site),data=dat)
anova(SUM)


#######IMPORT DATA INTO NINJA SOFTWARE#######
####### Follow instructions on website below #######
## https://shiny.wur.nl/ninja/
## Input NINJA values
data<- read_csv("NINJA_forR.csv")
head(data)


##### Anova with NINJA data #######
## set up linear models & run anovas

MI<-lmer(Maturity.Index~N_addition*Plant_removal+(1|Site), data=data)
anova(MI)

MI25<-lmer(Maturity.Index.2.5~N_addition*Plant_removal+(1|Site), data=data)
anova(MI25)

SMI<-lmer(Sigma.Maturity.Index~N_addition*Plant_removal+(1|Site), data=data)
anova(SMI)

PPI<-lmer(Plant.Parasitic.Index~N_addition*Plant_removal+(1|Site), data=data)
anova(PPI)

BI<-lmer(Basal.Index~N_addition*Plant_removal+(1|Site), data=data)
anova(BI)

EI<-lmer(Enrichment.Index~N_addition*Plant_removal+(1|Site), data=data)
anova(EI)

SI<-lmer(Structure.Index~N_addition*Plant_removal+(1|Site), data=data)
anova(SI)

#####Test nematode Family Stats###
## Change families as needed
tax <- read_csv("Family_Count.csv")
head(tax)
## ALL BF

Ceph<-lmer(data=tax, Cephalobidae ~ N_addition*Plant_removal + (1|Site))
anova(Ceph)

Pris<-lmer(data=tax, Prismatolaimidae ~ N_addition*Plant_removal + (1|Site))
anova(Pris)

Plec<-lmer(data=tax, Plectidae ~ N_addition*Plant_removal + (1|Site))
anova(Plec)

Tera<-lmer(data=tax, Teratocephalidae ~ N_addition*Plant_removal + (1|Site))
anova(Tera)

Rhab<-lmer(data=tax, Rhabditidae ~ N_addition*Plant_removal + (1|Site))
anova(Rhab)

Alai<-lmer(data=tax, Alaimidae ~ N_addition*Plant_removal + (1|Site))
anova(Alai)

Desm<-lmer(data=tax, Desmodoridae ~ N_addition*Plant_removal + (1|Site))
anova(Desm)

Rhablaim<-lmer(data=tax, Rhabdolaimidae ~ N_addition*Plant_removal + (1|Site))
anova(Rhablaim)

Monh<-lmer(data=tax, Monhysteridae ~ N_addition*Plant_removal + (1|Site))
anova(Monh)

Pang<-lmer(data=tax, Pangrolaimidae ~ N_addition*Plant_removal + (1|Site))
anova(Pang)

Xyal<-lmer(data=tax, Xyalidae ~ N_addition*Plant_removal + (1|Site))
anova(Xyal)

## ALL RA
Tyle<-lmer(data=tax, Tylenchidae ~ N_addition*Plant_removal + (1|Site))
anova(Tyle)

## ALL FF

Aphe<-lmer(data=tax, Aphelenchidae ~ N_addition*Plant_removal + (1|Site))
anova(Aphe)

Apho<-lmer(data=tax, Aphelenchoididae ~ N_addition*Plant_removal + (1|Site))
anova(Apho)

Lamell<-lmer(data=tax, Tylencholaimellidae ~ N_addition*Plant_removal + (1|Site))
anova(Lamell)

## ALL OM 
Apor<-lmer(data=tax, Aporcelaimidae ~ N_addition*Plant_removal + (1|Site))
anova(Apor)

Belo<-lmer(data=tax, Belondiridae ~ N_addition*Plant_removal + (1|Site))
anova(Belo)

Dory<-lmer(data=tax, Dorylaimidae ~ N_addition*Plant_removal + (1|Site))
anova(Dory)

Nord<-lmer(data=tax, Nordiidae ~ N_addition*Plant_removal + (1|Site))
anova(Nord)

Qud<-lmer(data=tax, Qudsianematidae ~ N_addition*Plant_removal + (1|Site))
anova(Qud)

## ALL PP 

Ang<-lmer(data=tax, Anguinidae ~ N_addition*Plant_removal + (1|Site))
anova(Ang)

Cric<-lmer(data=tax, Criconematidae ~ N_addition*Plant_removal + (1|Site))
anova(Cric)

Doli<-lmer(data=tax, Dolichodoridae ~ N_addition*Plant_removal + (1|Site))
anova(Doli)

Hete<-lmer(data=tax, Heteroderidae ~ N_addition*Plant_removal + (1|Site))
anova(Hete)

Hopl<-lmer(data=tax, Hoplolaimidae ~ N_addition*Plant_removal + (1|Site))
anova(Hopl)

Melo<-lmer(data=tax, Meloidogynidae ~ N_addition*Plant_removal + (1|Site))
anova(Melo)

Prat<-lmer(data=tax, Pratylenchidae ~ N_addition*Plant_removal + (1|Site))
anova(Prat)

Ulid<-lmer(data=tax, Tylenchulidae ~ N_addition*Plant_removal + (1|Site))
anova(Ulid)

########################### Beta Diversity #############################
# import ASV 
#### NOTE: This can be done for each trophic group level by importing a file with only select groups
otu_taxa <- read_csv("OTU_FinalFinal.csv")
otu_only <- otu_taxa[,2:41]##select for only the count data
otu_only_order <- otu_only[,order(colnames(otu_only))]

##import map file
mapping<-read.csv("META_Final.csv", row.names = )
map <- mapping[order(rownames(mapping)),]

#check that names match for rows and columns
rownames(map) == colnames(otu_only_order)

#aggregate to ASV  level
sample_no <- length(map[,1])
otu_Troph <- aggregate(otu_only_order[,1:sample_no], by=list(otu_taxa$OTU_ID), FUN=sum)
###change the column name to match trophic levels instead of ASVs
rownames(otu_Troph) <- otu_Troph[,1]

##get rid of wierd first group.2 labeled column, modify based on your data
#check your data to see if you have this, odd column, if not ignore this step and move on. 
otu <- otu_Troph[,2:41]


#order treatments in specific order you wan to see them on the plot. So the example has all of my lake sediments, then shoreline, then prairies
#I have the sample names in the mapping file under id, change the variable after the $s according to what you want on the x axis
map$Plot = factor(map$Plot, levels=c("XX","GX","DX","XN", "GN", "DN"))

sample_design <- sample_data(map)
otu_all = phyloseq(Nema18S_otu, sample_data(map))
otu_all

#Bray Curtis
## Can use "jaccard" here
bray.dist<-vegdist(t(otu), method="bray",na.rm = TRUE) 


#Calculate Distances
P.sd_bacs_can_czm_clr.dist <- vegdist(t(bray.dist), method = "euclidean", diag = FALSE, upper = TRUE)

##### Beta-diversity Plotting #####
## Define variables for plot code
plot.dm.bray<-as.matrix(P.sd_bacs_can_czm_clr.dist) # Change the source of the plot data here



# FUNCTION:
Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor)
} 
# END OF FUNCTION
# Now run the Principle Coordinate Analysis ("pcoa") in package ape:
pcoa.pts<-pcoa(plot.dm.bray,correction="cailliez")
x.minimum<-min(pcoa.pts$vectors[,1])
x.maximum<-max(pcoa.pts$vectors[,1])
y.minimum<-min(pcoa.pts$vectors[,2])
y.maximum<-max(pcoa.pts$vectors[,2])

# calculate the % variation explained by each axis
# getting percent explained by each PC vector
vars_percent <- (pcoa.pts$values$Eigenvalues / sum(pcoa.pts$values$Eigenvalues)) * 100
# Step 1: make very simple plot, no colors #main="Nematode Community", # CHANGE this to be the title of the plot


plot(pcoa.pts$vectors[,2] ~ pcoa.pts$vectors[,1],
     xlab=paste("PCo1: ", round(vars_percent[1], 2), "%", sep=""),
     ylab=paste("PCo2: ", round(vars_percent[2], 2), "%", sep=""),
     cex=1,col="white", 
     xlim=c(x.minimum-0.1, x.maximum+0.05), ylim=c(y.minimum-0.05, y.maximum+0.05))


# Step 2: now subset vectors for the positions in the alpha file that correspond to cohcryo
map.GX <- which(map$Plot =="GX") ## CHANGE these depending on the number of categories in your variable of interest
map.XN <- which(map$Plot =="XN") ## CHANGE these depending on the number of categories in your variable of interest
map.XX <- which(map$Plot =="XX") ## CHANGE these depending on the number of categories in your variable of interest
map.DX <- which(map$Plot =="DX") ## CHANGE these depending on the number of categories in your variable of interest
map.DN <- which(map$Plot =="DN") ## CHANGE these depending on the number of categories in your variable of interest
map.GN <- which(map$Plot =="GN") ## CHANGE these depending on the number of categories in your variable of interest

# plot the points for each group using points()
points(pcoa.pts$vectors[map.GN,2] ~ pcoa.pts$vectors[map.GN,1],pch=16, cex=2, col="red")
points(pcoa.pts$vectors[map.GX,2] ~ pcoa.pts$vectors[map.GX,1],pch=16, cex=2, col="blue")
points(pcoa.pts$vectors[map.XN,2] ~ pcoa.pts$vectors[map.XN,1], pch=16, cex=2, col="firebrick")
points(pcoa.pts$vectors[map.XX,2] ~ pcoa.pts$vectors[map.XX,1],pch=16, cex=2, col="darkblue")
points(pcoa.pts$vectors[map.DX,2] ~ pcoa.pts$vectors[map.DX,1], pch=16, cex=2, col="lightblue")
points(pcoa.pts$vectors[map.DN,2] ~ pcoa.pts$vectors[map.DN,1],pch=16, cex=2, col="palevioletred1")

Plot_ConvexHull(xcoord =pcoa.pts$vectors[map.GN,1], ycoord = pcoa.pts$vectors[map.GN,2], lcolor = "red")
Plot_ConvexHull(xcoord = pcoa.pts$vectors[map.GX,1], ycoord = pcoa.pts$vectors[map.GX,2], lcolor = "blue")
Plot_ConvexHull(xcoord = pcoa.pts$vectors[map.XN,1], ycoord = pcoa.pts$vectors[map.XN,2],lcolor = "firebrick")
Plot_ConvexHull(xcoord = pcoa.pts$vectors[map.XX,1], ycoord = pcoa.pts$vectors[map.XX,2],lcolor = "darkblue")
Plot_ConvexHull(xcoord = pcoa.pts$vectors[map.DX,1], ycoord = pcoa.pts$vectors[map.DX,2],lcolor = "lightblue")
Plot_ConvexHull(xcoord = pcoa.pts$vectors[map.DN,1], ycoord = pcoa.pts$vectors[map.DN,2],lcolor = "palevioletred1")


legend(-.8,-.1, legend=c("XX","GX","DX","XN","GN","DN"), col=c("darkblue","blue","lightblue","firebrick","red", "palevioletred1"), pch=c(16,16,16,16,16,16), cex= 1.2)
dev.off()


# PERMANOVA, change to variables you want to assess
euk.czm.perm <- adonis(P.sd_bacs_can_czm_clr.dist~map$Plant_removal*map$N_addition,strata=map$Site)


euk.czm.perm 


########## distance based redundancy analyses ##############
### first check for correlation of metadata

head(mapping)

cormap<-mapping[,7:23] ### select all variables of interest
M=cor(cormap)
corrplot(M, method = 'circle', type = 'upper')

## decide based on your data what to keep in the dbRDA and what is too correlated

## Bray Curtis
bray = vegdist(t(otu))

## set up dbRDA with distance matrix and variables
dbrda.m1<-dbrda(data=map,bray ~ Plant_removal + N_addition+ pH+Geum.Biomass+Deschampsia.Biomass+Litter+Other.Plants.Biomass+Bacterial.Richness+Fungal.Richness+NO3+TDN+H20+Plant.Richness+Root.Biomass)

#plot model results
map$N_addition=as.factor(map$N_addition)
with(map, levels(map$N_addition))
scl<-2
colvec<-c("red","blue")
plot(dbrda.m1, type="n",scaling=scl)
with(map, points(dbrda.m1, display = "sites", col = colvec[N_addition],
                 scaling = scl, pch = 21, bg = colvec[N_addition]))
text(dbrda.m1, col="blue",display="bp")
with(map, legend(-2,1, legend = levels(N_addition), bty = "n",
                 col = colvec, pch = 21, pt.bg = colvec))

#test significance of model terms
anova.cca(dbrda.m1, by='term')
