######################Loading packages #############################
options(warn=-1)
require(lme4)
require(MASS)
require(ggeffects)
#require(cowplot)
require(AICcmodavg)
require(ggplot2)
require(car)#for multicollinearity
require(sjPlot) #Data visualization
require(sjlabelled)# labelling
require(sjmisc)#Data variable and transformation function

#Setting working directory
setwd("C:\\Users\\majal\\Documents\\MSc analysis")

#Listing files
list.files()

#Importing files to R
eland<- read.csv("elandG.csv", header=T)

attach(eland)#Attach eland

names(eland)#Gives the names of the column

#Checking correlation of explanatory variables
explvars<-cbind(D_forbs,Vegetation.cover,Vegetation.biomass,Heap.dung, 
                Escape.impediment,daily_rain,NDVI,dNDVI,protection_gradient,
                proximity_confluence,Proximity_road,Proximity_drainage,Proximity_river,Proximity_multiplePAs,Proximity_CorePAs_edge,
                Elevation,slope,TWI,Curvature,Proximity_fence,Tourism_footprint,Human_footprint,woody_density)
  
cor.explvars<-cor(explvars, use="complete.obs")# correlation matrix
cor.explvars<- data.frame(cor.explvars)#putting in a dataframe
#write.csv(cor.explvars,"correlation.csv")#save it to CsV
#Elevation and Proximity to multiple use area
#are left out because of strong correlation with other variables
############################## Transforming the variables ####################
# for continous data logarithmic transformation is used
eland$D_forbs<-log(D_forbs + 1)#Density of forbs
eland$Vegetation.cover<-log(Vegetation.cover + 1)#Vegetation cover
eland$Vegetation.biomass<-log10(Vegetation.biomass)#biomass
eland$protection_gradient<-log10(protection_gradient)#protection gradient
eland$Proximity_road<-log10(Proximity_road)#proximity to road
eland$Proximity_river<-log10(Proximity_river)#proximity to river
eland$Proximity_CorePAs_edge<-log10(Proximity_CorePAs_edge)#Distance to protected areas edge
eland$slope<-log10(slope)#Slope
eland$TWI<-log10(TWI)#Topographic Wetness Index
eland$Tourism_footprint<- log10(Tourism_footprint)#Tourism footprint
eland$Human_footprint<-log10(Human_footprint)#Distance to human
eland$woody_density<-log(woody_density +1)#woody density
eland$daily_rain<-log(daily_rain + 1) ##daily rain
eland$NDVI<-log(NDVI +1) #NDVI
eland$dNDVI<-log(dNDVI +1)#changes in NDVI from the previous
eland$proximity_confluence<-log10(proximity_confluence)
eland$Proximity_drainage<-log10(Proximity_drainage)
eland$Proximity_fence<-log10(Proximity_fence)

######### For count data square root transformation is used
eland$Heap.dung<- sqrt(Heap.dung +0.5)# Dung heap
eland$Escape.impediment<-sqrt(Escape.impediment + 0.5)#Escape impediment
detach(eland)
############################ Standardizing the data ##############
attach(eland)
eland$D_forbs<-(D_forbs-mean(D_forbs))/sd(D_forbs)
eland$Vegetation.cover<-(Vegetation.cover-mean(Vegetation.cover))/sd(Vegetation.cover)
eland$Vegetation.biomass<-(Vegetation.biomass-mean(Vegetation.biomass))/sd(Vegetation.biomass)
eland$protection_gradient<-(protection_gradient-mean(protection_gradient))/sd(protection_gradient)
eland$Proximity_road<-(Proximity_road-mean(Proximity_road))/sd(Proximity_road)
eland$Proximity_river<-(Proximity_river-mean(Proximity_river))/sd(Proximity_river)
eland$Proximity_CorePAs_edge<-(Proximity_CorePAs_edge-mean(Proximity_CorePAs_edge))/sd(Proximity_CorePAs_edge)
eland$slope<-(slope-mean(slope))/sd(slope)
eland$Curvature<-(Curvature-mean(Curvature))/sd(Curvature)
eland$TWI<-(TWI-mean(TWI))/sd(TWI)
eland$Tourism_footprint<- (Tourism_footprint-mean(Tourism_footprint))/sd(Tourism_footprint)
eland$Heap.dung<-(Heap.dung-mean(Heap.dung))/sd(Heap.dung)
eland$Escape.impediment<-(Escape.impediment-mean(Escape.impediment))/sd(Escape.impediment)
eland$Human_footprint<-(Human_footprint-mean(Human_footprint))/sd(Human_footprint)
eland$woody_density<-(woody_density-mean(woody_density,na.rm=T))/sd(woody_density,na.rm=T)
eland$daily_rain<-(daily_rain-mean(daily_rain,na.rm=T))/sd(daily_rain,na.rm=T)
eland$NDVI<-(NDVI-mean(NDVI,na.rm=T))/sd(NDVI,na.rm=T)
eland$dNDVI<-(dNDVI-mean(dNDVI,na.rm=T))/sd(dNDVI,na.rm=T)
eland$proximity_confluence<-(proximity_confluence-mean(proximity_confluence,na.rm=T))/sd(proximity_confluence,na.rm=T)
eland$Proximity_drainage<-(Proximity_drainage-mean(Proximity_drainage,na.rm=T))/sd(Proximity_drainage,na.rm=T)
eland$Proximity_fence<-(Proximity_fence-mean(Proximity_fence,na.rm=T))/sd(Proximity_fence,na.rm=T)

detach(eland)
#################### Generalized mixed model selection #########
attach(eland)
##Global model
M1<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + Heap.dung + Curvature +
            Escape.impediment+ protection_gradient + Proximity_road + Proximity_CorePAs_edge +
            Proximity_river + Tourism_footprint + Micro.habitat.type + proximity_confluence + Proximity_drainage +
            Proximity_fence + slope + TWI + Human_footprint+ woody_density + daily_rain +NDVI +dNDVI + I(daily_rain^2)+
            I(Vegetation.cover^2)+ I(Vegetation.biomass^2)+ I(NDVI^2)+ I(dNDVI^2)+ I(proximity_confluence^2)+ 
            I(Proximity_drainage^2)+ I(Heap.dung^2)+ I(Escape.impediment^2)+ I(protection_gradient^2) + I(TWI^2)+I(slope^2)+
            I(Proximity_CorePAs_edge^2)+ I(Proximity_fence^2)+I(Proximity_river^2)+ I(Human_footprint^2)+ I(Curvature^2)+
            I(woody_density^2)+ Vegetation.cover*Proximity_river + Vegetation.biomass*Human_footprint + 
           daily_rain*(dNDVI +NDVI + D_forbs + Vegetation.cover + Vegetation.biomass) +(1|Animal_ID), family=binomial)

##M2
M2<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + Heap.dung + protection_gradient + woody_density+ 
            Proximity_road + daily_rain + NDVI +dNDVI +Proximity_CorePAs_edge +  Human_footprint +
            Tourism_footprint + Micro.habitat.type + proximity_confluence + Proximity_drainage +slope+ Curvature+
            I(NDVI^2)+ I(Proximity_CorePAs_edge^2)+ I(daily_rain^2)+I(dNDVI^2)+I(proximity_confluence^2)+
            I(Vegetation.cover^2)+ I(Vegetation.biomass^2)+ I(Proximity_drainage^2)+ I(Heap.dung^2)+ 
            I(protection_gradient^2) + I(woody_density^2) + I(Human_footprint^2)+ I(Curvature^2)+
            (1|Animal_ID), family=binomial)
#M3
M3<-glmer(Use~D_forbs + Vegetation.biomass + Heap.dung + woody_density +
            NDVI +dNDVI + slope + daily_rain +Vegetation.cover +proximity_confluence+
            Proximity_drainage +Proximity_river + Escape.impediment+Curvature +
            Proximity_CorePAs_edge +protection_gradient + Proximity_road + 
            Tourism_footprint + Human_footprint+ Micro.habitat.type + 
            Vegetation.cover*(Proximity_river +proximity_confluence)+
            (1|Animal_ID), family=binomial)

#Remove proximity to confluence and TWI
M4<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + Heap.dung + Curvature +
            Escape.impediment+ protection_gradient + Proximity_road + Proximity_CorePAs_edge +
            Proximity_river + Tourism_footprint + Micro.habitat.type + Proximity_drainage +
            Proximity_fence + slope + Human_footprint+ woody_density + daily_rain +NDVI +dNDVI +
            daily_rain*(D_forbs + Vegetation.cover + Vegetation.biomass + NDVI + dNDVI)+
            Vegetation.cover*Proximity_river + Vegetation.biomass*Human_footprint +
            (1|Animal_ID), family=binomial)

#Remove proximity to river, TWI and vegetation cover
M5<- glmer(Use~D_forbs + Vegetation.biomass + Heap.dung + Curvature+
             Escape.impediment+ protection_gradient + Proximity_road + Proximity_CorePAs_edge +
             Tourism_footprint + Micro.habitat.type + proximity_confluence + Proximity_drainage +
             Proximity_fence + slope + Human_footprint+ woody_density + daily_rain +NDVI +dNDVI +
             Human_footprint*Vegetation.biomass +(1|Animal_ID), family=binomial)

#Remove proximity to confluence,TWI, and vegetation biomass
M6<-glmer(Use~D_forbs +  Heap.dung + woody_density + Proximity_river + Vegetation.cover + Proximity_drainage +
            Escape.impediment+ protection_gradient + Proximity_road + Proximity_CorePAs_edge +
            Tourism_footprint + Proximity_fence + slope + daily_rain + NDVI +dNDVI +Curvature+
            Micro.habitat.type +(1|Animal_ID), family=binomial)

#M7(Without interaction)
M7<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + Heap.dung + Curvature +
            Escape.impediment+ protection_gradient + Proximity_road + Proximity_CorePAs_edge +
            Proximity_river + Tourism_footprint + Micro.habitat.type + proximity_confluence + Proximity_drainage +
            Proximity_fence + slope + TWI + Human_footprint+ woody_density + daily_rain +NDVI +dNDVI +
            (1|Animal_ID), family=binomial)

best_model<-list(M1,M2,M3,M4,M5,M6,M7)# selecting the best model

print(as.data.frame(aictab(best_model)),sort=T)
summary(M6)
############################ Plotting Odd ratios############################
x11()
set_theme(theme_sjplot())
plot_model(M6,show.values = T,transform = "plogis",
           value.offset =.4,value.size = 4,
           dot.size =3,
           title="",
           line.size =.9,
           width =.3,colors = "black",
           vline.color = "Gray72",
           axis.labels = c("Habitat type(Woodland)","Habitat type(Short grassland)",
                           "Habitat type(Open woodland)","Curvature","dNDVI","NDVI","Rainfall","Slope",
                           "Proximity to fence","Tourism footprint","Proximity to CorePAs edge",
                           "Proximity to road","Protection gradient","Escape impediment",
                           "Proximity to drainage","Cover for lions","Proximity to river",
                           "Woody density","Interspecies interaction","Density of forbs")) +
  theme(text= element_text(size= 12, 
                           colour = "black", face="bold", family= "Times"),
        panel.border = element_rect(colour = "black", 
                                    fill=NA, size=2),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+ 
  ylab("Relative probability of selection")
###################################################################
###################### Dividing the data seasonally################
elandD<-subset(eland,Season=="Dry")#Dry season
elandW<-subset(eland, Season=="Wet")#Wet season

detach(eland)#Detaching eland data
attach(elandD)#Attaching dry season

##################### Modelling dry season ######################
##Global model
Md1<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
             Heap.dung + protection_gradient + woody_density+ 
             Proximity_road + daily_rain + NDVI +dNDVI +Proximity_CorePAs_edge + 
             Human_footprint + Tourism_footprint + Micro.habitat.type + 
             proximity_confluence + Proximity_drainage +slope +Proximity_fence +
             Curvature+ slope +TWI+ Proximity_river+ Escape.impediment+
             (1|Animal_ID), family=binomial)

##Removing proximity to drainage
Md2<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
             Heap.dung + protection_gradient + woody_density+ 
             Proximity_road + daily_rain + NDVI +dNDVI +Proximity_CorePAs_edge + 
             Human_footprint + Tourism_footprint + Micro.habitat.type + 
             proximity_confluence + slope +Proximity_fence + Curvature+ slope +
             TWI+ Proximity_river+ Escape.impediment+
             (1|Animal_ID), family=binomial)

##Remove proximity to confluence
Md3<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
             Heap.dung + protection_gradient + woody_density+ 
             Proximity_road + daily_rain + NDVI +dNDVI +Proximity_CorePAs_edge + 
             Human_footprint + Tourism_footprint + Micro.habitat.type + 
             slope +Proximity_fence + Curvature+ slope +TWI+ Proximity_river+ 
             Escape.impediment+(1|Animal_ID), family=binomial)

##Remove proximity to river
Md4<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
             Heap.dung + protection_gradient + woody_density+ 
             Proximity_road + daily_rain + NDVI +dNDVI +Proximity_CorePAs_edge + 
             Human_footprint + Tourism_footprint + Micro.habitat.type + 
             slope +Proximity_fence + Curvature+ slope +TWI+ 
             Escape.impediment+(1|Animal_ID), family=binomial)

#Remove TWI
Md5<- glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
              Heap.dung + protection_gradient + woody_density+ 
              Proximity_road + daily_rain + NDVI +dNDVI +Proximity_CorePAs_edge + 
              Human_footprint + Tourism_footprint + Micro.habitat.type + 
              slope +Proximity_fence +Curvature+ slope + Escape.impediment+
              (1|Animal_ID), family=binomial)

#Remove slope
Md6<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
             Heap.dung + protection_gradient + woody_density+ 
             Proximity_road + daily_rain + NDVI +dNDVI +Proximity_CorePAs_edge + 
             Human_footprint + Tourism_footprint + Micro.habitat.type + 
             Proximity_fence + Curvature+ Escape.impediment+
             (1|Animal_ID), family=binomial)

#Remove escape impediment
Md7<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
             Heap.dung + protection_gradient + woody_density+ 
             Proximity_road + daily_rain + NDVI + dNDVI +Proximity_CorePAs_edge + 
             Human_footprint + Tourism_footprint + Micro.habitat.type + 
             Proximity_fence + Curvature+ (1|Animal_ID), family=binomial)

#Remove Heap dung
Md8<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
             protection_gradient + woody_density + Proximity_road +
             daily_rain + NDVI + dNDVI +Proximity_CorePAs_edge + 
             Human_footprint + Tourism_footprint + Micro.habitat.type + 
             Proximity_fence + Curvature+ (1|Animal_ID), family=binomial)

#Based on randomly selected variables from the data frame
Md9<-glmer(Use ~ D_forbs +Vegetation.biomass + woody_density + Heap.dung +
             Proximity_river + Curvature + TWI+ slope + daily_rain + 
             NDVI  + Vegetation.cover + proximity_confluence + 
             Proximity_drainage + Proximity_road + Human_footprint +
             protection_gradient +(1 | Animal_ID),family = binomial)

best_modelD<-list(Md1,Md2,Md3,Md4,Md5,Md6,Md7,Md8,Md9)# selecting the best model

print(as.data.frame(aictab(best_modelD)),sort=T)# 
##Seems Md9 has lowest AIC
summary(Md9)
#############Plotting Odd ratios
x11()
set_theme(theme_sjplot())
plot_model(Md9,show.values = T,
           value.offset =.4,value.size = 4,
           dot.size =3,
           title="",
           line.size =.9,
           width =.3,colors = "black",
           vline.color = "white",
           axis.labels = c("Protection gradient","Human footprint","Proximity to road","Proximity to drainage",
                           "Proximity to confluence","Cover for lions","NDVI","Rainfall",
                           "Slope","Topographic wetness index","Landscape curvature","Proximity to river","Interspecies interaction",
                           "Woody density","Vegetation biomass","Density of forbs")) +
  theme(text= element_text(size= 12, 
                           colour = "black", face="bold", family= "Times"),
        panel.border = element_rect(colour = "black", 
                                    fill=NA, size=2),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+ 
  geom_hline(yintercept = 1,lty="dashed",size=1,colour="red")+
  ylab("Odd ratios")
###################################################################
############################# Wet season Modelling#################
detach(elandD)
attach(elandW)
##Global model
Mw1<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
             Heap.dung + protection_gradient + woody_density+ 
             Proximity_road + daily_rain + NDVI +dNDVI +Proximity_CorePAs_edge + 
             Human_footprint + Tourism_footprint + Micro.habitat.type + 
             proximity_confluence + Proximity_drainage +slope +Proximity_fence +
             Curvature+ slope + TWI + Proximity_river+ Escape.impediment+
             (1|Animal_ID), family=binomial)

#Remove escape impediment,woody density,proximity to drainage
Mw2<-glmer(Use~D_forbs +Vegetation.biomass + Heap.dung + 
             Proximity_river + TWI+ slope + daily_rain + Curvature + NDVI + dNDVI+
             Vegetation.cover+ proximity_confluence + Proximity_CorePAs_edge + 
             Tourism_footprint + Proximity_fence + protection_gradient + 
             Human_footprint + Proximity_road +Micro.habitat.type + (1|Animal_ID),family = binomial)

##Removing proximity to drainage
Mw3<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
             Heap.dung + protection_gradient + woody_density+ 
             Proximity_road + daily_rain + NDVI +dNDVI +Proximity_CorePAs_edge + 
             Human_footprint + Tourism_footprint + Micro.habitat.type + 
             proximity_confluence + slope +Proximity_fence +
             Curvature+ slope + TWI + Proximity_river+ Escape.impediment+
             (1|Animal_ID), family=binomial)

##Remove proximity to confluence
Mw4<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
             Heap.dung + protection_gradient + woody_density+ 
             Proximity_road + daily_rain + NDVI +dNDVI +Proximity_CorePAs_edge + 
             Human_footprint + Tourism_footprint + Micro.habitat.type + 
             Proximity_drainage +slope +Proximity_fence + Curvature+ slope + 
             TWI + Proximity_river+ Escape.impediment+
             (1|Animal_ID), family=binomial)

##Remove proximity to river
Mw5<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
             Heap.dung + protection_gradient + woody_density+ 
             Proximity_road + daily_rain + NDVI +dNDVI +Proximity_CorePAs_edge + 
             Human_footprint + Tourism_footprint + Micro.habitat.type + 
             proximity_confluence + Proximity_drainage +slope +Proximity_fence +
             Curvature+ slope + TWI + Escape.impediment+
             (1|Animal_ID), family=binomial)

#Remove TWI
Mw6<- glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
              Heap.dung + protection_gradient + woody_density+ 
              Proximity_road + daily_rain + NDVI +dNDVI +Proximity_CorePAs_edge + 
              Human_footprint + Tourism_footprint + Micro.habitat.type + 
              proximity_confluence + Proximity_drainage +slope +Proximity_fence +
              Curvature+ slope + Proximity_river+ Escape.impediment+
              (1|Animal_ID), family=binomial)

#Remove slope
Mw7<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
             Heap.dung + protection_gradient + woody_density+ 
             Proximity_road + daily_rain + NDVI +dNDVI +Proximity_CorePAs_edge + 
             Human_footprint + Tourism_footprint + Micro.habitat.type + 
             proximity_confluence + Proximity_drainage +slope +Proximity_fence +
             Curvature+ TWI + Proximity_river+ Escape.impediment+
             (1|Animal_ID), family=binomial)

#Remove escape impediment
Mw8<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
             Heap.dung + protection_gradient + woody_density+ 
             Proximity_road + daily_rain + NDVI +dNDVI +Proximity_CorePAs_edge + 
             Human_footprint + Tourism_footprint + Micro.habitat.type + 
             proximity_confluence + Proximity_drainage +slope +Proximity_fence +
             Curvature+ slope + TWI + Proximity_river+ (1|Animal_ID), family=binomial)

#Remove Heap dung
Mw9<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
             protection_gradient + woody_density+Proximity_road + 
             daily_rain + NDVI +dNDVI +Proximity_CorePAs_edge + 
             Human_footprint + Tourism_footprint + Micro.habitat.type + 
             proximity_confluence + Proximity_drainage +slope +Proximity_fence +
             Curvature+ slope + TWI + Proximity_river+ Escape.impediment+
             (1|Animal_ID), family=binomial)

best_modelW<-list(Mw1,Mw2,Mw3,Mw4,Mw5,Mw6,Mw7,Mw8,Mw9)# selecting the best model

print(as.data.frame(aictab(best_modelW)),sort=T)# 
#Mw2 has lowest AIC
summary(Mw2)
############################ Plotting Odd ratios############################
x11()
set_theme(theme_sjplot())
plot_model(Mw2,show.values = T,
    terms = c("D_forbs","Vegetation.biomass","Heap.dung","Proximity_river",
              "TWI","slope","daily_rain","Curvature","NDVI","dNDVI",
              "Vegetation.cover","proximity_confluence","Proximity_CorePAs_edge", 
              "Tourism_footprint","Proximity_fence","protection_gradient", 
              "Human_footprint","Proximity_road"),
           value.offset =.4,value.size = 4,
           dot.size =3,
           title="",
           line.size =.9,
           width =.3,colors = "black",
           vline.color = "white",
           axis.labels = c("Proximity to road","Human footprint","Protection gradient",
                           "Proximity to fence","Tourism footprint","Proximity to CorePAs edge",
                           "Proximity to confluence","Cover for lions","dNDVI","NDVI","Landscape curvature",
                           "Rainfall","Slope","Topographic wetness index","Proximity to river",
                           "Interspecies interaction","Vegetation biomass","Density of forbs")) +
  theme(text= element_text(size= 12, 
                           colour = "black", face="bold", family= "Times"),
        panel.border = element_rect(colour = "black", 
                                    fill=NA, size=2),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+ 
  geom_hline(yintercept = 1,lty="dashed",size=1,colour="blue")+
  ylab("Odd ratios")
###################################################################