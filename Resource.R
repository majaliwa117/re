###################### Eland Resource selection Analysis########
## Loading packages
rm(list=ls())#Clearing memory in the working environment 
options(warn =-1)#Supressing warning messages
require(ggplot2)# For data visualization
require(survival)# For conditional logistic regression
require(AICcmodavg)# For model selection
require(pROC)
require(visreg)#Visualizing regression models
require(dplyr) #Data manipulation
require(latex2exp)#Exponentiating expression

#Setting working directory
setwd("C:\\Users\\majal\\Documents\\MSc analysis")

#Listing files
list.files()

#Importing files to R
eland<- read.csv("elandG.csv", header=T)

attach(eland)

#Checking correlation of explanatory variables
explvars<-cbind(D_forbs,Vegetation.cover,Vegetation.biomass,Heap.dung, 
                Escape.impediment,protection_gradient,Proximity_road,
                Proximity_river,dist_Pasedges,Tourism_footprint,Elevation,slope,TWI)
cor.explvars<-cor(explvars, use="complete.obs")# correlation matrix
cor.explvars<- data.frame(cor.explvars)#putting in a dataframe
#write.csv(cor.explvars,"correlation.csv")#save it to CsV

####################### Checking for normality#####################
#Using Shapiro-Wilk normality test
shapiro.test(D_forbs);shapiro.test(Vegetation.cover);shapiro.test(Vegetation.biomass)
shapiro.test(Heap.dung);shapiro.test(Escape.impediment)
shapiro.test(protection_gradient);shapiro.test(Proximity_road);shapiro.test(Proximity_river)
shapiro.test(dist_Pasedges);shapiro.test(Tourism_footprint);shapiro.test(Elevation)
shapiro.test(slope);shapiro.test(TWI)#All variables are not normally distributed

############################## Transforming the variables ####################
# for continous data logarithmic transformation is used
eland$D_forbs<-log(D_forbs + 1)#Density of forbs
eland$Vegetation.cover<-log(Vegetation.cover + 1)#Vegetation cover
eland$Vegetation.biomass<-log10(Vegetation.biomass)#biomass
eland$protection_gradient<-log10(protection_gradient)#protection gradient
eland$Proximity_road<-log10(Proximity_road)#proximity to road
eland$Proximity_river<-log10(Proximity_river)#proximity to river
eland$dist_Pasedges<-log10(dist_Pasedges)#Distance to protected areas edge
eland$Elevation<-log10(Elevation) #Elevation
eland$slope<-log10(slope)#Slope
eland$TWI<-log10(TWI)#Topographic Wetness Index
eland$Tourism_footprint<- log10(Tourism_footprint)#Tourism footprint

# For count data square root transformation is used
eland$Heap.dung<- sqrt(Heap.dung +0.5)# Dung heap
eland$Escape.impediment<-sqrt(Escape.impediment + 0.5)#Escape impediment
############################ Standardizing the data ##############
eland$D_forbs<-(D_forbs-mean(D_forbs))/sd(D_forbs)
eland$Vegetation.cover<-(Vegetation.cover-mean(Vegetation.cover))/sd(Vegetation.cover)
eland$Vegetation.biomass<-(Vegetation.biomass-mean(Vegetation.biomass))/sd(Vegetation.biomass)
eland$protection_gradient<-(protection_gradient-mean(protection_gradient))/sd(protection_gradient)
eland$Proximity_road<-(Proximity_road-mean(Proximity_road))/sd(Proximity_road)
eland$Proximity_river<-(Proximity_river-mean(Proximity_river))/sd(Proximity_river)
eland$dist_Pasedges<-(dist_Pasedges-mean(dist_Pasedges))/sd(dist_Pasedges)
eland$Elevation<-(Elevation-mean(Elevation))/sd(Elevation)
eland$slope<-(slope-mean(slope))/sd(slope)
eland$TWI<-(TWI-mean(TWI))/sd(TWI)
eland$Tourism_footprint<- (Tourism_footprint-mean(Tourism_footprint))/sd(Tourism_footprint)
eland$Heap.dung<-(Heap.dung-mean(Heap.dung))/sd(Heap.dung)
eland$Escape.impediment<-(Escape.impediment-mean(Escape.impediment))/sd(Escape.impediment)

#####################Model selection ###############################

#Global model
M1<-clogit(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
             Heap.dung + Escape.impediment+ protection_gradient + Proximity_road + 
             Proximity_river + dist_Pasedges + Tourism_footprint + Micro.habitat.type + 
             slope + TWI + I(Tourism_footprint^2)+ I(Vegetation.biomass^2)+ cluster(Animal_ID)+ 
             strata(Season), method = "approximate")

#M2
M2<-clogit(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
             Heap.dung + Escape.impediment+ protection_gradient + Proximity_road + 
             Proximity_river + Vegetation.biomass*Proximity_river + dist_Pasedges + 
             Tourism_footprint + Micro.habitat.type + slope + TWI+ cluster(Animal_ID)+
             strata(Season), method = "approximate")
#M3
M3<-clogit(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
             Heap.dung + Escape.impediment+ protection_gradient + Proximity_road + 
             Proximity_river + dist_Pasedges + Tourism_footprint + Micro.habitat.type + 
             Vegetation.cover*Vegetation.biomass + slope + TWI+ cluster(Animal_ID)+
             strata(Season), method = "approximate")
#M4
M4<-clogit(Use~Vegetation.cover + Vegetation.biomass +
             Heap.dung + Escape.impediment+ protection_gradient + Proximity_road + 
             Proximity_river +  dist_Pasedges + Tourism_footprint + Micro.habitat.type + 
             slope + TWI+ cluster(Animal_ID)+ strata(Season), method = "approximate")
#M5
M5<-clogit(Use~ Vegetation.biomass + Heap.dung + Escape.impediment+ 
             protection_gradient + Proximity_road + Proximity_river + 
             dist_Pasedges + Tourism_footprint + Micro.habitat.type +
             slope + TWI+ cluster(Animal_ID)+ strata(Season), method = "approximate")
#M6
M6<-clogit(Use~ Heap.dung + Escape.impediment+ protection_gradient + Proximity_road + 
          Proximity_river +  dist_Pasedges + Tourism_footprint + Micro.habitat.type + 
         slope + TWI+ cluster(Animal_ID)+ strata(Season), 
           method = "approximate")
#M7
M7<-clogit(Use~ Heap.dung + Escape.impediment+ protection_gradient + Proximity_road + 
             Proximity_river + dist_Pasedges + Tourism_footprint + Micro.habitat.type + 
             slope + TWI+ cluster(Animal_ID)+ strata(Season), method = "approximate")
#M8
M8<-clogit(Use~Escape.impediment+ protection_gradient + Proximity_road + 
             Proximity_river + dist_Pasedges + Tourism_footprint + Micro.habitat.type + 
             slope + TWI+ cluster(Animal_ID)+ strata(Season), method = "approximate")

#M9
M9<-clogit(Use~ protection_gradient + Proximity_road + Proximity_river + 
              dist_Pasedges +Tourism_footprint + Micro.habitat.type + slope + TWI+ 
             cluster(Animal_ID)+ strata(Season), method = "approximate")

#M10
M10<-clogit(Use~ Proximity_road + Proximity_river + dist_Pasedges + 
             Tourism_footprint + Micro.habitat.type + slope + TWI+ cluster(Animal_ID)+ 
              strata(Season), method = "approximate")
#M11
M11<-clogit(Use~Proximity_river + dist_Pasedges + Tourism_footprint + Micro.habitat.type + 
             slope + TWI+ cluster(Animal_ID)+ strata(Season), method = "approximate")

#M12
M12<-clogit(Use~ dist_Pasedges + Tourism_footprint + Micro.habitat.type + 
             slope + TWI+ cluster(Animal_ID)+ strata(Season), method = "approximate")

#M13
M13<-clogit(Use~ Tourism_footprint + Micro.habitat.type + 
            slope + TWI+ cluster(Animal_ID)+ strata(Season), method = "approximate")

#M14
M14<-clogit(Use~ Micro.habitat.type + slope + TWI+ 
            cluster(Animal_ID)+ strata(Season), method = "approximate")

#M15
M15<-clogit(Use~ slope + TWI+ cluster(Animal_ID)+ strata(Season), method = "approximate")

#M16
M16<-clogit(Use~TWI+ cluster(Animal_ID)+ strata(Season), method = "approximate")

#Null
M17<-clogit(Use~ 1 + cluster(Animal_ID)+ strata(Season), method = "approximate")

best_model<-list(M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11,M12,M13,M14,
                 M15,M16,M17)
print(as.data.frame(aictab(best_model,sort=T)))# 
##The model selection seems not to work because it turns out that
## null model is the best model,therefore am proceeding with global model
summary(M1)

######################## Plotting ##############################
x11()
par(mfrow=c(2,3),mar=c(5,5,3,3),omi=c(0.2,0.2,0,0))#Arrange the graph in three columns and two rows

#Density of forbs
visreg(M1, "D_forbs",line=list(col="black"),
       xlab=expression('Density of forbs/M'^2),ylab="Selection coefficient")
abline(h=0, lty="dotted");text(2,0.5, "p<0.05", font=12)

#Vegetation biomass
visreg(M1, "Vegetation.biomass",xlab="Vegetation biomass (Kg/ha)",
       ylab="Selection coefficient",line=list(col="black"))
abline(h=0, lty="dotted");text(2,0.5, "p<0.05", font=12)

#Interspecies interaction
visreg(M1, "Heap.dung",xlab="Interspecies interaction",ylab="Selection coefficient",line=list(col="black"))
abline(h=0, lty="dotted");text(5,0.5, "p>0.05", font=12)

#Cover for lions
visreg(M1, "Vegetation.cover",xlab="Cover for lions(%)",
       ylab="Selection coefficient",line=list(col="black"))
abline(h=0, lty="dotted");text(40,0.5, "p<0.05", font=12)

#Escape impediment
visreg(M1, "Escape.impediment",xlab="Escape impediment",ylab="Selection coefficient",line=list(col="black"))
abline(h=0, lty="dotted");text(1,0.5, "p>0.05", font=12)

#Protected area edges
visreg(M1, "dist_Pasedges",xlab="Proximity to PAs edge",ylab="Selection coefficient",line=list(col="black"))
abline(h=0, lty="dotted");text(5,0.5, "p>0.05", font=12)

x11()
par(mfrow=c(2,3))#Arrange the graph in three columns and two rows

#Protection gradient
visreg(M1, "protection_gradient",xlab="Protection gradient(Km)",ylab="Selection coefficient",line=list(col="black"))
abline(h=0, lty="dotted");text(5,0.5, "p<0.05", font=12)

#Proximity to road
visreg(M1, "Proximity_road",xlab="Proximity to road(Km)",ylab="Selection coefficient",line=list(col="black"))
abline(h=0, lty="dotted");text(5,0.5, "p>0.05", font=12)

#Tourism footprint
visreg(M1, "Tourism_footprint",xlab="Tourism footprint",ylab="Selection coefficient",line=list(col="black"))
abline(h=0, lty="dotted");text(0.04,0.5, "p>0.05", font=12)

#Proximity to river
visreg(M1, "Proximity_river",xlab="Proximity to river(Km)",ylab="Selection coefficient",line=list(col="black"))
abline(h=0, lty="dotted");text(5,0.5, "p<0.05", font=12)

#Slope
visreg(M1, "slope",xlab="Slope",ylab="Selection coefficient",line=list(col="black"))
abline(h=0, lty="dotted");text(0.5,0.5, "p<0.05", font=12)

#Topographic Wetness Index
visreg(M1, "TWI",xlab="Topographic Wetness Index",ylab="Selection coefficient",line=list(col="black"))
abline(h=0, lty="dotted");text(12,0.5, "p>0.05", font=12)

##############Resource selection in dry and wet season##########
#Subsetting dry season only
eland_dry<-eland[eland$Season=="Dry",]

#Subsetting wet season only
eland_wet<-eland[eland$Season=="Wet",]

############################## Modelling for dry season ###########################
#Modelling for dry season
detach(eland)#detaching eland dataframe
attach(eland_dry) #Attaching eland dataframe for dry season only

Md<-clogit(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
             Heap.dung + Escape.impediment+ protection_gradient + Proximity_road + 
             Proximity_river + dist_Pasedges + Tourism_footprint + Micro.habitat.type + 
             slope + TWI + I(Tourism_footprint^2)+ I(Vegetation.biomass^2)+  Proximity_river*Vegetation.cover +
             cluster(Animal_ID)+ strata(Time), method = "approximate")

summary(Md)

#Plotting
x11()
par(mfrow=c(2,3))#Arrange the graph in three columns and two rows

#Density of forbs
visreg(Md, "D_forbs",line=list(col="black"),rug=2,
       xlab=expression('Density of forbs/M'^2),ylab="Selection coefficient")
abline(h=0, lty="dotted");text(2,0.5, "p<0.05", font=12)

#Vegetation biomass
visreg(Md, "Vegetation.biomass",xlab="Vegetation biomass (Kg/ha)",rug=2,
       ylab="Selection coefficient",line=list(col="black"))
abline(h=0, lty="dotted");text(2,1, "p<0.05", font=12)

#Interspecies interaction
visreg(Md, "Heap.dung",xlab="Interspecies interaction",ylab="Selection coefficient",line=list(col="black"),rug=2)
abline(h=0, lty="dotted");text(5,0.5, "p>0.05", font=12)

#Cover for lions
visreg(Md, "Vegetation.cover",xlab="Cover for lions(%)",rug=2,
       ylab="Selection coefficient",line=list(col="black"))
abline(h=0, lty="dotted");text(1,0.5, "p>0.05", font=12)

#Escape impediment
visreg(Md, "Escape.impediment",xlab="Escape impediment",ylab="Selection coefficient",line=list(col="black"),rug=2)
abline(h=0, lty="dotted");text(0.5,1, "p>0.05", font=12)

#Protected area edges
visreg(Md, "dist_Pasedges",xlab="Proximity to PAs edge",ylab="Selection coefficient",line=list(col="black"),rug=2)
abline(h=0, lty="dotted");text(1,1, "p>0.05", font=12)

x11()
par(mfrow=c(2,3))#Arrange the graph in three columns and two rows

#Protection gradient
visreg(Md, "protection_gradient",xlab="Protection gradient(Km)",ylab="Selection coefficient",line=list(col="black"),rug=2)
abline(h=0, lty="dotted");text(1,1, "p>0.05", font=12)

#Proximity to road
visreg(Md, "Proximity_road",xlab="Proximity to road(Km)",ylab="Selection coefficient",line=list(col="black"),rug=2)
abline(h=0, lty="dotted");text(1,1, "p>0.05", font=12)

#Tourism footprint
visreg(Md, "Tourism_footprint",xlab="Tourism footprint",ylab="Selection coefficient",line=list(col="black"),rug=2)
abline(h=0, lty="dotted");text(0.04,0.5, "p>0.05", font=12)

#Proximity to river
visreg(Md, "Proximity_river",xlab="Proximity to river(Km)",ylab="Selection coefficient",line=list(col="black"),rug=2)
abline(h=0, lty="dotted");text(1,1, "p>0.05", font=12)

#Slope
visreg(Md, "slope",xlab="Slope",ylab="Selection coefficient",line=list(col="black"),rug=2)
abline(h=0, lty="dotted");text(0.5,0.5, "p<0.05", font=12)

#Topographic Wetness Index
visreg(Md, "TWI",xlab="Topographic Wetness Index",ylab="Selection coefficient",line=list(col="black"),rug=2)
abline(h=0, lty="dotted");text(1,1, "p>0.05", font=12)

############################## Modelling for Wet season ###########################
#Modelling for dry season
detach(eland_dry)#detaching eland dry season dataframe
attach(eland_wet) #Attaching eland wet season dataframe

Mw<-clogit(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
             Heap.dung + Escape.impediment+ protection_gradient + Proximity_road + 
             Proximity_river + dist_Pasedges + Tourism_footprint + Micro.habitat.type + 
             slope + TWI + I(Tourism_footprint^2)+ I(Vegetation.biomass^2)+ 
             cluster(Animal_ID)+ strata(Time), method = "approximate")

summary(Mw)

x11()
par(mfrow=c(2,3))#Arrange the graph in three columns and two rows

#Density of forbs
visreg(Mw, "D_forbs",line=list(col="black"),rug=2,
       xlab=expression('Density of forbs/M'^2),ylab="Selection coefficient")
abline(h=0, lty="dotted");text(2,0.5, "p>0.05", font=12)

#Vegetation biomass
visreg(Mw, "Vegetation.biomass",xlab="Vegetation biomass (Kg/ha)",rug=2,
       ylab="Selection coefficient",line=list(col="black"))
abline(h=0, lty="dotted");text(1,1, "p<0.05", font=12)

#Interspecies interaction
visreg(Mw, "Heap.dung",xlab="Interspecies interaction",ylab="Selection coefficient",line=list(col="black"),rug=2)
abline(h=0, lty="dotted");text(1,1, "p>0.05", font=12)

#Cover for lions
visreg(Mw, "Vegetation.cover",xlab="Cover for lions(%)",rug=2,
       ylab="Selection coefficient",line=list(col="black"))
abline(h=0, lty="dotted");text(1,0.5, "p<0.05", font=12)

#Escape impediment
visreg(Mw, "Escape.impediment",xlab="Escape impediment",ylab="Selection coefficient",line=list(col="black"),rug=2)
abline(h=0, lty="dotted");text(4,1, "p>0.05", font=12)

#Protected area edges
visreg(Mw, "dist_Pasedges",xlab="Proximity to PAs edge(Km)",ylab="Selection coefficient",line=list(col="black"),rug=2)
abline(h=0, lty="dotted");text(1,1, "p>0.05", font=12)

x11()
par(mfrow=c(2,3))#Arrange the graph in three columns and two rows

#Protection gradient
visreg(Mw, "protection_gradient",xlab="Protection gradient(Km)",ylab="Selection coefficient",line=list(col="black"),rug=2)
abline(h=0, lty="dotted");text(1,0.5, "p>0.05", font=12)

#Proximity to road
visreg(Mw, "Proximity_road",xlab="Proximity to road(Km)",ylab="Selection coefficient",line=list(col="black"),rug=2)
abline(h=0, lty="dotted");text(1,0.5, "p>0.05", font=12)

#Tourism footprint
visreg(Mw, "Tourism_footprint",xlab="Tourism footprint",ylab="Selection coefficient",line=list(col="black"),rug=2)
abline(h=0, lty="dotted");text(2,0.5, "p>0.05", font=12)

#Proximity to river
visreg(Mw, "Proximity_river",xlab="Proximity to river(Km)",ylab="Selection coefficient",line=list(col="black"),rug=2)
abline(h=0, lty="dotted");text(1,0.5, "p<0.05", font=12)

#Slope
visreg(Mw, "slope",xlab="Slope",ylab="Selection coefficient",line=list(col="black"),rug=2)
abline(h=0, lty="dotted");text(0.5,0.5, "p<0.05", font=12)

#Topographic Wetness Index
visreg(Mw, "TWI",xlab="Topographic Wetness Index",ylab="Selection coefficient",line=list(col="black"),rug=2)
abline(h=0, lty="dotted");text(1,0.5, "p>0.05", font=12)
