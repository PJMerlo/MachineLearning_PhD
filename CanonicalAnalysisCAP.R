rm(list=ls())
##### 
#load packages
library(ecodist)
library(ggplot2)
library(BiodiversityR)
library(ggThemeAssist)
library(vegan)
library(ggrepel)

##### 
mypath <- "c:/Users/merlo/Documents" #modify to yout route
setwd(paste(mypath, sep=""))

dat <- read.csv("abundance_df.csv")

#Get the columns with species abundance values 
spp_abs_dat <- as.data.frame(dat[,1:13])

#Get the relative abundance by  sample 
spp_rel_dat <- as.data.frame(sapply(dat[,1:13], function(x) x / rowSums(dat[,1:13]))) #columns 1 to 13 contain spp
spp_rel_dat[is.na(spp_rel_dat)] <- 0 #replaces NaNs with zeros

rowSums(spp_rel_dat) #expected only ONEs or ZEROs (in samples with zero species registered)

# Already have spp dataframe, now create dataframe containing only factors:
factors_data <- as.data.frame(dat[,16:20])
#colnames(factors_data)  <- "Metodo"

#####
## CAP (package BiodiversityR)
#i will use a GOWER DISTANCE MATRIX AND RELATIVE DATA!!!"""
# 1°: With "MÉTODO" as FACTOR
# 2°: With "HABITAT" as FACTOR

#DEFINE THE VARIABLES AS FACTORS
factors_data$Metodo <-as.factor(factors_data$Metodo)
factors_data$Habitat <-as.factor(factors_data$Habitat)

#####
##Performing CAP analysis using "METODO" as the CATEGORICAL VARIABLE:
Ordination_model1 <- CAPdiscrim(spp_rel_dat ~ Metodo, data=factors_data, dist="gower",  axes=2, m=0, add=FALSE)
Ordination_model1

#Overall classification success (m=6) : 80 percent !
#BRUVS (n=30) correct: 70 percent
#UVC (n=30) correct: 90 percent

#          Df  Pillai approx F num Df den   Df  Pr(>F)  
#y[, group]  1 0.3414   4.579      6     53   0.0008199 ***
#Residuals  58

#Ahora tenés que agregar species scores manually with add.spec.scores
Ordination_model1 <- add.spec.scores(Ordination_model1,comm=spp_rel_dat,method="cor.scores",scaling="1") #comm=df with sites: rows, species: columns and spp abundance as cell values.

plot1 <- ordiplot(Ordination_model1, choices=c(1,2) )#, type="none")  borrar parentesis before '#'
ordisymbol(plot1, factors_data, "Metodo", legend=TRUE, cex= 1.75) #"managament dice que distinto colores y formas segun ese grupo
add.spec.scores(plot1,spp_rel_dat,method="cor.scores",multi=1,Rscale=F,scaling="1")

#Info on the locations of sites (circles in the ordiplot) is obtained via function sites.long. Info on spp is extracted by function species.long
sites.long1 <- sites.long(plot1, env.data=factors_data) #acá extraes los 2 primeros axis de la ordenacion para plot en ggplot

head(sites.long1)

species.long1 <- species.long(plot1)
species.long1

axis.long <- axis.long(Ordination_model1, choices=c(1, 2)) # info of labeling of the axes
axis.long

#Now we can generate the plot:
#set theme:
BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12, family="Arial"),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 13),
)#legend.key = element_blank())

#it is better to reduce the n of species to show in the plot to those that explain more of the variation:
#This info is obtained with the envfit function of vegan, using r (coef. of correlation)
spec.envfit <- envfit(plot1, env=spp_rel_dat) #env= transformed matrix of disimilarity
spec.data.envfit <- data.frame(r=spec.envfit$vectors$r, p=spec.envfit$vectors$pvals)
species.long2 <- species.long(plot1, spec.data=spec.data.envfit)
species.long2   # check r values of the species (Anderson et al. 2008)

species.long3 <- species.long2[species.long2$r >= 0.22, ] #filter spp with r>0.22 (arbitrary cut-off, Anderson et al. 2008)
species.long3

#now replot using that info:
plot_gower10_Metodo <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab("CAP1") +
  ylab("CAP2") +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
  geom_point(data=sites.long1, 
             aes(x=axis1, y=axis2, colour=Metodo, shape=Habitat), 
             size=3) +
  geom_segment(data=species.long3, 
               aes(x=0, y=0, xend=axis1*4, yend=axis2*4), 
               colour="red", size=0.5, arrow=arrow()) +
  geom_text_repel(data=species.long3, 
                  aes(x=axis1*4.99, y=axis2*4, label=labels),
                  colour="red",cex=4.5) +
  BioR.theme +
  ggsci::scale_colour_npg() +
  coord_fixed(ratio=1)
plot_gower10_Metodo
#save figure
#dev.copy(tiff,'CAP_gower_Metodo_relativedata.tiff', res = 300, width = 1920, height = 1440, pointsize = 10)
#dev.off()

#####
##Performing CAP analysis using "HABITAT" as the CATEGORICAL VARIABLE:
##CON HABITAT instead of Metodo as FACTOR
Ordination_model2 <- CAPdiscrim(spp_rel_dat ~ Habitat, dist="gower", data=factors_data,  axes=2, m=0, add=FALSE)
Ordination_model2

#Overall classification success (m=2) : 90 percent  
#Reef (n=30) correct: 86.66 percent
#Soft-bottom (n=30) correct: 93.333 percent

#'HÁBITAT' allows for a higher classification score (90% samples classified correctly) than 'METODO' (method)(80%)
Ordination_model2 <- add.spec.scores(Ordination_model2,comm=spp_rel_dat,method="cor.scores",scaling="1") #comm=df with sites: rows, species: columns and spp abundance as cell values.

plot2 <- ordiplot(Ordination_model2, choices=c(1,2) )#, type="none")  borrar parentesis before '#'
ordisymbol(plot2, factors_data, "Habitat", legend=TRUE, cex= 1.5) #"managament dice que distinto colores y formas segun ese grupo
add.spec.scores(plot2,spp_rel_dat ,method="cor.scores",multi=1,Rscale=F,scaling="1")

#Info on the locations of sites (circles in the ordiplot) is obtained via function sites.long. Info on spp is extracted by function species.long
sites.long2 <- sites.long(plot2, env.data=factors_data) #acá extraes los 2 primeros axis de la ordenacion para plot en ggplot

head(sites.long2)

species.long_b <- species.long(plot2)
species.long_b

axis.long2 <- axis.long(Ordination_model2, choices=c(1, 2)) # info of labeling of the axes
axis.long2

spec.envfit2 <- envfit(plot2, env=spp_rel_dat) #env= transformed matrix of disimilarity
spec.data.envfit2 <- data.frame(r=spec.envfit2$vectors$r, p=spec.envfit2$vectors$pvals)
species.long_b <- species.long(plot2, spec.data=spec.data.envfit2)
species.long_b   # ver que spp tienen r > 0.5 u otro valor de interes (anderson et al. 2008)

species.long_b2 <- species.long_b[species.long_b$r >= 0.22, ] #dejamos solo spp de r>0.22
species.long_b2

plot_gower10_Habitat <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab("CAP1") +
  ylab("CAP2") +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
  geom_point(data=sites.long2, 
             aes(x=axis1, y=axis2, colour=Habitat, shape=Metodo), 
             size=3) +
  geom_segment(data=species.long_b2, 
               aes(x=0, y=0, xend=axis1*3.8, yend=axis2*3.8), 
               colour="grey20", size=0.7, arrow=arrow(length = unit(0.12,"inches"))) +
  geom_text_repel(data=species.long_b2, 
                  aes(x=axis1*4.99, y=axis2*4, label=labels),
                  colour="red",cex=4.5) +
  BioR.theme +
  ggsci::scale_colour_npg() +
  coord_fixed(ratio=1)

plot_gower10_Habitat
#save figure
#dev.copy(tiff,'CAP_gower_Habitat_relativedata.tiff', res = 300, width = 2200, height = 2000, pointsize = 10)
#dev.off()

#References:
#Anderson, M. (2008). PERMANOVA+ for PRIMER: guide to software and statistical methods. Primer-E Limited.