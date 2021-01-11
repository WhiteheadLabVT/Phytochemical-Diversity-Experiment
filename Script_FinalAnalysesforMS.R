#load libraries
library(dplyr)
library(lme4)
library(multcomp)
library(emmeans)
library(MuMIn)
library(stringr)
library(ggplot2)
library(tidyr)
library(MASS) #for negative binomial GLMs in Pred 2a/2b and stepwise model selection in Pred 3b
#library(sjstats) #for ICC (intraclass correlation coefficient) calcs


##Data standardization must be run first and workspace saved. After that, all other
#sections of code can be run independently

# Data standardization and calculations-----
#--------------------

#load data
d.rich <- read.csv("Whitehead_et_al_Insects.csv")  
d.even <- read.csv("Whitehead_et_al_Insects_Evenness.csv")
d.fungi <- read.csv("Whitehead_et_al_Fungi.csv")

d.rich[c(1:3, 5:6, 8, 11)] <- lapply(d.rich[c(1:3, 5:6, 8, 11)], factor)
d.even[c(1:3, 5:6, 8, 11)] <- lapply(d.even[c(1:3, 5:6, 8, 11)], factor)
d.fungi[2:6] <- lapply(d.fungi[2:6], factor)

#creating standardized variables for pupal weight and development time that will be a ratio of the 
#performance on the treatment relative to the performance of relevant controls

#Also taking an inverse value of days to pupation, so that, for some analyses and plots,
##we can visualize all performance metrics 
#with the same direction (increase in value = better performance). This can be thought of as 
# "development speed"


for (i in 1:length(d.rich$Pupal.weight)){
  Sp <- d.rich$Species[i]
  Sex <- d.rich$Sex[i]
  exp <- d.rich$exp[i]
  Cval <- mean(d.rich$Pupal.weight[which(d.rich$Treatment=="C" & d.rich$Species==Sp &  
                                           d.rich$Sex==Sex & d.rich$exp==exp)], na.omit=TRUE)
  d.rich$Pupal.weight.ST[i] <- d.rich$Pupal.weight[i]/Cval
}


for (i in 1:length(d.rich$Days.to.pupation)){
  Sp <- d.rich$Species[i]
  Sex <-d.rich$Sex[i]
  exp <- d.rich$exp[i]
  Cval <- mean(d.rich$Days.to.pupation[which(d.rich$Treatment=="C" & d.rich$Species==Sp &
                                               d.rich$Sex==Sex & d.rich$exp==exp)], na.omit=TRUE)
  d.rich$Days.to.pupation.ST[i] <- d.rich$Days.to.pupation[i]/Cval
  #also calculating inverse value so that we can discuss all performance metrics with
  #positive numbers as better performance. This can be thought of as "development speed"
  d.rich$Days.to.pupation.ST.inv[i] <- Cval/d.rich$Days.to.pupation[i]
}
d.rich$Days.to.pupation.ST.inv <- d.rich$Days.to.pupation.ST.inv %>% replace(.,is.nan(.),  NA)


#same for evenness dataset
for (i in 1:length(d.even$Pupal.weight)){
  Sp <- d.even$Species[i]
  Sex <- d.even$Sex[i]
  exp <- d.even$exp[i]
  Cval <- mean(d.even$Pupal.weight[which(d.even$Treatment=="C" & d.even$Species==Sp &  
                                           d.even$Sex==Sex & d.even$exp==exp)], na.omit=TRUE)
  d.even$Pupal.weight.ST[i] <- d.even$Pupal.weight[i]/Cval
}


for (i in 1:length(d.even$Days.to.pupation)){
  Sp <- d.even$Species[i]
  Sex <-d.even$Sex[i]
  exp <- d.even$exp[i]
  Cval <- mean(d.even$Days.to.pupation[which(d.even$Treatment=="C" & d.even$Species==Sp &
                                               d.even$Sex==Sex & d.even$exp==exp)], na.omit=TRUE)
  d.even$Days.to.pupation.ST[i] <- d.even$Days.to.pupation[i]/Cval
  #also calculating inverse value so that we can discuss all performance metrics with
  #positive numbers as better performance. This can be thought of as "development speed"
  d.even$Days.to.pupation.ST.inv[i] <- Cval/d.even$Days.to.pupation[i]
}
d.even$Days.to.pupation.ST.inv <- d.even$Days.to.pupation.ST.inv %>% replace(.,is.nan(.),  NA)



#Creating summarized datasets for proportion survival (PS) that will group by diet treatment
#This was necessary for the GLMMs; models were unstable with the binary 0/1 variable for 
#survival and each herbivore as a replicate. With this summarized data we treat the
#specific diet treatment as the replicate and use a simplified random effects structure
d.rich.PS <- d.rich %>%  group_by(Species, Richness, SD, Treatment, exp) %>%
  summarise(dead=length(which(Surv==0)), alive=length(which(Surv==1)), 
            PropSurv=length(which(Surv==1))/length(Surv)) 
d.even.PS <- d.even %>%  group_by(Species, Evenness, SD, Treatment, exp) %>%
  summarise(dead=length(which(Surv==0)), alive=length(which(Surv==1)), 
            PropSurv=length(which(Surv==1))/length(Surv)) 

d.rich.PS$PropSurv.ST <- NA
#Creating standardized variable relative to controls
for (i in 1:length(d.rich.PS$Species)){
  exp <- d.rich.PS$exp[i]
  Sp <- d.rich.PS$Species[i]
  Cval <- mean(d.rich.PS$PropSurv[which(d.rich.PS$Treatment=="C" & d.rich.PS$exp==exp & 
                                          d.rich.PS$Species==Sp)])
  d.rich.PS$PropSurv.ST[i] <- d.rich.PS$PropSurv[i]/Cval
}

d.even.PS$PropSurv.ST <- NA
#Creating standardized variable relative to controls
for (i in 1:length(d.even.PS$Species)){
  exp <- d.even.PS$exp[i]
  Sp <- d.even.PS$Species[i]
  Cval <- mean(d.even.PS$PropSurv[which(d.even.PS$Treatment=="C" & d.even.PS$exp==exp & 
                                          d.even.PS$Species==Sp)])
  d.even.PS$PropSurv.ST[i] <- d.even.PS$PropSurv[i]/Cval
}


#Checking out sample sizes for richness dataset 
counts <- d.rich %>%
  group_by(Species, Treatment) %>%
  summarise(n.surv = sum(!is.na(Surv)), n.pw = sum(!is.na(Pupal.weight)))

hist(counts$n.surv[which(counts$Treatment != "C")])
range(counts$n.surv[which(counts$Treatment != "C")])
sum(counts$n.surv)
sum(counts$n.pw)

counts2 <- pivot_wider(counts[1:3], names_from = Species, values_from = n.surv) #total insects
counts3 <- pivot_wider(counts[c(1:2,4)], names_from = Species, values_from = n.pw)  #survivors

#Table S2 in supplement--gives sample sizes for each treatment/insect species combination
#as total number/number surviving. We only have data on pupal weights and days to pupation
#for insects that survived, so the sample size after the slash gives sample sizes for
#analyses with those variables
counts4 <- data.frame(Treatment=counts3$Treatment, Cp=paste(counts2$Cp, counts3$Cp, sep=" / "),
                      Hz=paste(counts2$Hz, counts3$Hz, sep=" / "), Px=paste(counts2$Px, counts3$Px, sep=" / "),
                      Sf=paste(counts2$Sf, counts3$Sf, sep=" / "), 
                      Total=paste(rowSums(counts2[2:5]), rowSums(counts3[2:5]), sep=" / "))

write.csv(counts4, "./Outputs/Tables/TableS2_InsectSampleSize.csv")

#Also checking sample sizes for evenness dataset 
counts <- d.even %>%
  group_by(Species, Treatment) %>%
  summarise(n.surv = sum(!is.na(Surv)), n.pw = sum(!is.na(Pupal.weight)))

hist(counts$n.surv[which(counts$Treatment != "C")])
range(counts$n.surv[which(counts$Treatment != "C")])
sum(counts$n.surv)
sum(counts$n.pw)
sum(counts$n.surv[which(counts$Treatment != "C")])
sum(counts$n.pw[which(counts$Treatment != "C")])

counts2 <- pivot_wider(counts[1:3], names_from = Species, values_from = n.surv) #total insects
counts3 <- pivot_wider(counts[c(1:2,4)], names_from = Species, values_from = n.pw)  #survivors



#For fungi data---creating standardized variable relative to controls

#taking mean of two duplicate readings
d.fungi <- d.fungi %>% rowwise() %>% mutate(MeanAbs=mean(c(Reading1, Reading2)))

#First checking out growth patterns of controls for each species
d.temp <- filter(d.fungi, Treatment=="DMSO" & Fungi=="Botrys")
plot(MeanAbs~Time, data=d.temp)

d.temp <- filter(d.fungi, Treatment=="DMSO" & Fungi=="Penicillium")
plot(MeanAbs~Time, data=d.temp)

d.temp <- filter(d.fungi, Treatment=="DMSO" & Fungi=="Collet")
plot(MeanAbs~Time, data=d.temp)

d.temp <- filter(d.fungi, Treatment=="DMSO" & Fungi=="Sclerotinia")
plot(MeanAbs~Time, data=d.temp)

#Also can switch out treatment here to look at patterns for any given treatment
d.temp <- filter(d.fungi, Treatment=="L10C")
plot(MeanAbs~Time, data=d.temp, col=as.factor(Fungi))

#some weird stuff is going on with quercetin, absorbance was high initially and then dropped
d.temp <- filter(d.fungi, Time==1)
plot(MeanAbs ~ as.factor(Treatment), data=d.temp)
#yes, quercetin at time 1 seems to be a major outlier, probably because it naturally has
#a yellow color

#Spreading data across time points so we can get delta Absorbance or slope

d.temp <- pivot_wider(data = d.fungi, id_cols=c(CellID, Fungi, Plate, Treatment),
            names_from = Time, 
            values_from = MeanAbs)
colnames(d.temp)[5:8] <- c("T1", "T2", "T3", "T4")

#calculating the slope of the line from T1-T4 for each well. For Quercetin,
#we will use the slope from T2-T4 because the first reading was so much higher
d.temp$slope <- NA
for(i in 1:length(d.temp$slope)){
  if(d.temp$Treatment[i]=="Q"){
    x <- c(2,3,4)
    y <- c(d.temp$T2[i], d.temp$T3[i], d.temp$T4[i])
    d.temp$slope[i] <- coef(lm(y~x))[2]
  } else{
    x <- c(1,2,3,4)
    y <- c(d.temp$T1[i], d.temp$T2[i], d.temp$T3[i], d.temp$T4[i])
    d.temp$slope[i] <- coef(lm(y~x))[2]
  }
}

d.fungi2 <- d.temp

#creating standardized variable--slope of treatment divided by slope of controls
for (i in 1:length(d.fungi2$slope)){
  Sp <- d.fungi2$Fungi[i]
  Cval <- mean(d.fungi2$slope[which(d.fungi2$Treatment=="DMSO" & 
                                      d.fungi2$Fungi==Sp)], na.omit=TRUE)   
  d.fungi2$dAbs.ST[i] <- (d.fungi2$slope[i]+1)/(Cval+1)
}

hist(d.fungi2$dAbs.ST)

#Adding richness and SD to fungi dataframe
d.rich.u <- distinct(d.rich[2:4], Treatment, .keep_all=T)  
d.fungi2$Treatment <- as.character(d.fungi2$Treatment)
d.rich.u$Treatment <- as.character(d.rich.u$Treatment)

d.fungi2 <- left_join(d.fungi2, d.rich.u, by="Treatment")

#taking out Captan (negative control using captan fungicide) and broth only control (no DMSO)
d.fungi2 <- filter(d.fungi2, Treatment !="Captan" & Treatment !="Broth")


#Use to get sample sizes for fungi 
counts <- d.fungi2 %>%
  group_by(Fungi, Treatment) %>%
  summarise(n= length(!is.na(dAbs.ST)))

sum(counts$n)

counts2 <- pivot_wider(counts, names_from = Fungi, values_from = n)
#four wells for everything!



#There was some concern during review that repetition of some treatments
#for insects (in just a few cases where we had low survival rates)
#could affect the outcome of the results. Thus, we also tried this analysis
#but dropping cases where there was repetition of treatments

#Checking out cases of repetition 
expreps <- d.rich %>%
  group_by(Species, Treatment, exp) %>%
  summarise(n.surv = sum(!is.na(Surv)), n.pw = sum(!is.na(Pupal.weight)))

table(expreps$Treatment)
#H8B for Px is in exp 6+7 (12 insects in 6, 11 died)
#H8C is in 6+7 for Sf, Hz, and Cp
#L2A is in 6+7 for Sf  (only one insect in 6)
#M2B is in 6+7 for Cp (only 2 insects in 6)
#M4B is in 6+7 for Px (should drop 7, only one insect)
#M8B is in 6+8, fully repeated in 8 for all insects
#R is in 6+8, fully repeated in 8 for all insects
#some of these with only one insect could be typos? 

d.rich2 <- d.rich
d.rich2 <- d.rich2[-which(d.rich2$Treatment=="H8B" & d.rich2$exp==6),] 
d.rich2 <- d.rich2[-which(d.rich2$Treatment=="H8C" & d.rich2$exp==6),] 
d.rich2 <- d.rich2[-which(d.rich2$Treatment=="L2A" & d.rich2$exp==6),] 
d.rich2 <- d.rich2[-which(d.rich2$Treatment=="M2B" & d.rich2$exp==6),] 
d.rich2 <- d.rich2[-which(d.rich2$Treatment=="M4B" & d.rich2$exp==7),] 
d.rich2 <- d.rich2[-which(d.rich2$Treatment=="M8B" & d.rich2$exp==6),] 
d.rich2 <- d.rich2[-which(d.rich2$Treatment=="R" & d.rich2$exp==6),] 

d.rich.PS2 <- d.rich2 %>%  group_by(Species, Richness, SD, Treatment, exp) %>%
  summarise(dead=length(which(Surv==0)), alive=length(which(Surv==1)), 
            PropSurv=length(which(Surv==1))/length(Surv)) 

d.rich.PS2$PropSurv.ST <- NA
#Creating standardized variable relative to controls
for (i in 1:length(d.rich.PS2$Species)){
  exp <- d.rich.PS2$exp[i]
  Sp <- d.rich.PS2$Species[i]
  Cval <- mean(d.rich.PS2$PropSurv[which(d.rich.PS2$Treatment=="C" & d.rich.PS2$exp==exp & 
                                          d.rich.PS2$Species==Sp)])
  d.rich.PS2$PropSurv.ST[i] <- d.rich.PS2$PropSurv[i]/Cval
}


#can try running analyses with d.rich2 to see if this affects outcomes


#Saving this workspace. Many other analyses below will depend on these standardized
#variables. 

save.image("./Outputs/Workspaces/StandardizedData")

#**************************************************************



# Effects of richness and SD on performance (Predictions 1a, 1b, 1c)-----
#--------------------

load("./Outputs/Workspaces/StandardizedData")

#------For Insects

#---all species combined

##Big picture first: looking at overall effects of richness and SD on standardized
#performance metrics, all species combined

##Now looking at overall effects of diversity using standardized variables

#Looking at effects of Richness, SD, Species, and Sex on pupal weights
d.temp <- filter (d.rich, Treatment != "C", !is.na(Sex), !is.na(Pupal.weight.ST))
m1.all.PW <- lmer(Pupal.weight.ST ~ SD*Richness*Species*Sex + (1|Treatment) + (1|exp/TrayID/CupID), 
                  data = d.temp, na.action=na.fail, REML="FALSE")
d1.all.PW <- dredge(m1.all.PW)
s <- subset(d1.all.PW, delta < 4)
s  #use for table S3
d1.all.PW.avg <- model.avg(d1.all.PW, subset = delta < 4, fit=TRUE)
summary(d1.all.PW.avg)  #use for table S4
m1.top.PW <- lmer(Pupal.weight.ST ~ Richness + Sex*Species + (1|Treatment) + (1|exp/TrayID/CupID), 
                  data = d.temp, na.action=na.fail) 
summary(m1.top.PW)
drop1(m1.top.PW, test="Chisq") #significant positive effect of richness in top model

plot(d.temp$Pupal.weight.ST ~ d.temp$Richness)
abline(lm(d.temp$Pupal.weight.ST ~ d.temp$Richness))

#Strong effects of species and sex and species*sex
###Marginal POSITIVE effect of richness on weights. 
#Effect size is small. but richness is in all top models and always positive
#In the top model, there is a significant positive effect of richness

#Effects on days to pupation
d.temp <- filter (d.rich, Treatment != "C", !is.na(Sex), !is.na(Days.to.pupation.ST.inv))
m1.all.dTp <- lmer(Days.to.pupation.ST.inv ~ SD*Richness*Species*Sex + (1|Treatment) + (1|exp/TrayID/CupID), 
                   data = d.temp, na.action=na.fail, REML="FALSE")
d1.all.dTp <- dredge(m1.all.dTp)
s <- subset(d1.all.dTp, delta < 4)
s  #use for table S3
d1.all.dTp.avg <- model.avg(d1.all.dTp, subset = delta < 4, fit=TRUE)
summary(d1.all.dTp.avg)  #use for table S4

m1.top.dTp <- lmer(Days.to.pupation.ST.inv ~ Richness + Sex*Species + Richness*Species + (1|Treatment) + (1|exp/TrayID/CupID), 
                   data = d.temp, na.action=na.fail) 
summary(m1.top.dTp)
drop1(m1.top.dTp, test="Chisq")

plot(d.temp$Days.to.pupation.ST.inv ~ d.temp$Richness)
abline(lm(d.temp$Days.to.pupation.ST.inv ~ d.temp$Richness))

#Strong effects of species
#Strong Richness *Species interactions
##Sex*Species interactions
#SD*Species SD*Sex interactions
#Richness is in every model, but effect is sometimes neg and sometimes pos


##Effects on survival

d.temp <- filter (d.rich.PS, Treatment != "C", !is.na(PropSurv.ST))
m1.all.S <- lmer(PropSurv.ST ~ SD*Richness*Species + (1|Treatment) + (1|exp), 
                 data = d.temp, na.action=na.fail, REML="FALSE")
d1.all.S <- dredge(m1.all.S)
s <- subset(d1.all.S, delta < 4)
s #use for table S3
d1.all.S.avg <- model.avg(d1.all.S, subset =  delta < 4, fit=TRUE)
summary(d1.all.S.avg)  #use for table S4

#OVerall only species is coming out as significant. Richness is in 3/5 top models, sometimes
#positive sometimes negative



#Considering the many interactions with species and sex, will also conduct
#these analyses separately for each herbivore species and performance metric for m and f

#----  Hz   

#making an empty table to store model results
m.sum <- data.frame(rxsd_chi=NA, rxsd_p=NA, r_chi=NA, r_p=NA, sd_chi=NA, sd_p=NA)

#female pupal weights
d.temp <- filter (d.rich, Species=="Hz" & Treatment != "C" & Sex=="f")
m1.Hz.PW.f <- lmer(Pupal.weight.ST ~ SD*Richness + (1|Treatment) + (1|exp), data = d.temp)
m2.Hz.PW.f <- lmer(Pupal.weight.ST ~ SD + Richness + (1|Treatment)+ (1|exp), data = d.temp)
d1 <- drop1(m1.Hz.PW.f, test="Chisq")
d2 <- drop1(m2.Hz.PW.f, test="Chisq")  #no effects
d1
d2
m.sum[1,] <- c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4])
vc <- data.frame(VarCorr(m1.Hz.PW.f))
vc <- mutate(vc, percVar=100*vcov/sum(vcov))
vc

#male pupal weights
d.temp <- filter (d.rich, Species=="Hz" & Treatment != "C" & Sex=="m")
m1.Hz.PW.m <- lmer(Pupal.weight.ST ~ SD*Richness + (1|Treatment) + (1|exp), data = d.temp)
m2.Hz.PW.m <- lmer(Pupal.weight.ST ~ SD + Richness + (1|Treatment) + (1|exp), data = d.temp)
d1 <- drop1(m1.Hz.PW.m, test="Chisq")
d2 <- drop1(m2.Hz.PW.m, test="Chisq")  #no effects
d1
d2
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))
vc <- data.frame(VarCorr(m1.Hz.PW.m))
vc <- mutate(vc, percVar=100*vcov/sum(vcov))
vc

#female days to pupation
d.temp <- filter (d.rich, Species=="Hz" & Treatment != "C" & Sex=="f")
m1.Hz.dTp.f <- lmer(Days.to.pupation.ST.inv ~ SD*Richness + (1|Treatment) + (1|exp), data = d.temp)
m2.Hz.dTp.f <- lmer(Days.to.pupation.ST.inv ~ SD + Richness + (1|Treatment) + (1|exp), data = d.temp)
d1 <- drop1(m1.Hz.dTp.f, test="Chisq")
d2 <- drop1(m2.Hz.dTp.f, test="Chisq")  #no effects
d1
d2
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))  
vc <- data.frame(VarCorr(m1.Hz.dTp.f))
vc <- mutate(vc, percVar=100*vcov/sum(vcov))
vc

#male days to pupation
d.temp <- filter (d.rich, Species=="Hz" & Treatment != "C" & Sex=="m")
m1.Hz.dTp.m <- lmer(Days.to.pupation.ST.inv ~ SD*Richness + (1|Treatment) + (1|exp), data = d.temp)
m2.Hz.dTp.m <- lmer(Days.to.pupation.ST.inv ~ SD + Richness + (1|Treatment) + (1|exp), data = d.temp)
d1 <- drop1(m1.Hz.dTp.m, test="Chisq")
d2 <- drop1(m2.Hz.dTp.m, test="Chisq") #no effects
d1
d2
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))
vc <- data.frame(VarCorr(m1.Hz.dTp.m))
vc <- mutate(vc, percVar=100*vcov/sum(vcov))
vc

#survival
d.temp <- filter (d.rich.PS, Species=="Hz" & Treatment != "C")

#First tried this as a glmer binomial, but models were not stable
#PS <- cbind(d.temp$alive, d.temp$dead)
#m1.Hz.Surv <- glmer(PS ~ SD*Richness + (1|exp), data = d.temp, family=binomial)
#m2.Hz.Surv <- glmer(PS ~ SD+Richness + (1|exp), data = d.temp, family=binomial)

m1.Hz.Surv <- lmer(PropSurv.ST ~ SD*Richness + (1|exp), data = d.temp)
m2.Hz.Surv <- lmer(PropSurv.ST ~ SD+Richness + (1|exp), data = d.temp)
d1 <- drop1(m1.Hz.Surv, test="Chisq")
d2 <- drop1(m2.Hz.Surv, test="Chisq")  #no effects; with glmer there is a significant effect of SD
d1
d2
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))

summary(m2.Hz.Surv)

#Tukey post-hoc for SD levels
m2.T <- glht(m2.Hz.Surv, linfct=mcp(SD="Tukey"))
summary(m2.T)
cld(m2.T, level = 0.05) #lower survival on Med SD
plot(d.temp$PropSurv ~ d.temp$SD)





#----  Sf 

#female pupal weights
d.temp <- filter (d.rich, Species=="Sf" & Treatment != "C" & Sex=="f")
m1.Sf.PW.f <- lmer(Pupal.weight.ST ~ SD*Richness + (1|Treatment) + (1|exp), data = d.temp)
m2.Sf.PW.f <- lmer(Pupal.weight.ST ~ SD + Richness + (1|Treatment)+ (1|exp), data = d.temp)
d1 <- drop1(m1.Sf.PW.f, test="Chisq")
d2 <- drop1(m2.Sf.PW.f, test="Chisq")  
d1
d2
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))
summary(m2.Sf.PW.f)  
plot(d.temp$Pupal.weight.ST ~d.temp$Richness)
abline(lm(d.temp$Pupal.weight.ST ~d.temp$Richness))
#marginal positive effect of richness on f pupal weights

#male pupal weights
d.temp <- filter (d.rich, Species=="Sf" & Treatment != "C" & Sex=="m")
m1.Sf.PW.m <- lmer(Pupal.weight.ST ~ SD*Richness + (1|Treatment) + (1|exp), data = d.temp)
m2.Sf.PW.m <- lmer(Pupal.weight.ST ~ SD + Richness + (1|Treatment) + (1|exp), data = d.temp)
d1 <- drop1(m1.Sf.PW.m, test="Chisq")
d2 <- drop1(m2.Sf.PW.m, test="Chisq")  #no effects
d1
d2
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))

#female days to pupation
d.temp <- filter (d.rich, Species=="Sf" & Treatment != "C" & Sex=="f")
m1.Sf.dTp.f <- lmer(Days.to.pupation.ST.inv ~ SD*Richness + (1|Treatment) + (1|exp), data = d.temp)
m2.Sf.dTp.f <- lmer(Days.to.pupation.ST.inv ~ SD  + Richness + (1|Treatment) + (1|exp), data = d.temp)
d1 <- drop1(m1.Sf.dTp.f, test="Chisq")
d2 <- drop1(m2.Sf.dTp.f, test="Chisq")  #marginal effect of Richness, slower development at high richness
d1
d2
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))

summary(m2.Sf.dTp.f)
plot(Days.to.pupation.ST.inv ~ Richness, data=d.temp)
abline(lm(Days.to.pupation.ST.inv ~ Richness, data=d.temp))

#male days to pupation
d.temp <- filter (d.rich, Species=="Sf" & Treatment != "C" & Sex=="m")
m1.Sf.dTp.m <- lmer(Days.to.pupation.ST.inv ~ SD*Richness + (1|Treatment) + (1|exp), data = d.temp)
m2.Sf.dTp.m <- lmer(Days.to.pupation.ST.inv ~ SD + Richness + (1|Treatment) + (1|exp), data = d.temp)
d1 <- drop1(m1.Sf.dTp.m, test="Chisq")
d2 <- drop1(m2.Sf.dTp.m, test="Chisq") #no effects
d1
d2
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))

#survival
d.temp <- filter (d.rich.PS, Species=="Sf" & Treatment != "C")
m1.Sf.Surv <- lmer(PropSurv.ST ~ SD*Richness + (1|exp), data = d.temp)
m2.Sf.Surv <- lmer(PropSurv.ST ~ SD+Richness + (1|exp), data = d.temp)
d1 <- drop1(m1.Sf.Surv, test="Chisq")
d2 <- drop1(m2.Sf.Surv, test="Chisq")  #positive effect of richness on survival
d1
d2
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))

summary(m2.Sf.Surv) 
plot(d.temp$PropSurv.ST ~ d.temp$Richness)
abline(lm(d.temp$PropSurv.ST ~ d.temp$Richness))
#when run with glmer binomial, there was a much
#stronger positive effect of richness on survival
#marginal effect of SD, with lower survival on Med SD


#----  Cp   

#female pupal weights
d.temp <- filter (d.rich, Species=="Cp" & Treatment != "C" & Sex=="f")
m1.Cp.PW.f <- lmer(Pupal.weight.ST ~ SD*Richness + (1|Treatment) + (1|exp), data = d.temp)
m2.Cp.PW.f <- lmer(Pupal.weight.ST ~ SD + Richness + (1|Treatment)+ (1|exp), data = d.temp)
d1 <- drop1(m1.Cp.PW.f, test="Chisq")
d2 <- drop1(m2.Cp.PW.f, test="Chisq")  #no effects
d1
d2
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))

#male pupal weights
d.temp <- filter (d.rich, Species=="Cp" & Treatment != "C" & Sex=="m")
m1.Cp.PW.m <- lmer(Pupal.weight.ST ~ SD*Richness + (1|Treatment) + (1|exp), data = d.temp)
m2.Cp.PW.m <- lmer(Pupal.weight.ST ~ SD + Richness + (1|Treatment) + (1|exp), data = d.temp)
d1 <- drop1(m1.Cp.PW.m, test="Chisq")
d2 <- drop1(m2.Cp.PW.m, test="Chisq")  #marginal positive effect of richness
d1
d2
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))
summary(m2.Cp.PW.m)
plot(d.temp$Pupal.weight.ST ~ d.temp$Richness)
abline(lm(d.temp$Pupal.weight.ST ~ d.temp$Richness))


#female days to pupation
d.temp <- filter (d.rich, Species=="Cp" & Treatment != "C" & Sex=="f")
m1.Cp.dTp.f <- lmer(Days.to.pupation.ST.inv ~ SD*Richness + (1|Treatment) + (1|exp), data = d.temp)
m2.Cp.dTp.f <- lmer(Days.to.pupation.ST.inv ~ SD + Richness + (1|Treatment) + (1|exp), data = d.temp)
d1 <- drop1(m1.Cp.dTp.f, test="Chisq")
d2 <- drop1(m2.Cp.dTp.f, test="Chisq")  #no effects
d1
d2
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))

#male days to pupation
d.temp <- filter (d.rich, Species=="Cp" & Treatment != "C" & Sex=="m")
m1.Cp.dTp.m <- lmer(Days.to.pupation.ST.inv ~ SD*Richness + (1|Treatment) + (1|exp), data = d.temp)
m2.Cp.dTp.m <- lmer(Days.to.pupation.ST.inv ~ SD + Richness + (1|Treatment) + (1|exp), data = d.temp)
d1 <- drop1(m1.Cp.dTp.m, test="Chisq")
d2 <- drop1(m2.Cp.dTp.m, test="Chisq") #no effects
d1
d2
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))

#survival
d.temp <- filter (d.rich.PS, Species=="Cp" & Treatment != "C")
m1.Cp.Surv <- lmer(PropSurv.ST ~ SD*Richness + (1|exp), data = d.temp)
m2.Cp.Surv <- lmer(PropSurv.ST ~ SD+Richness + (1|exp), data = d.temp)
d1 <- drop1(m1.Cp.Surv, test="Chisq")
d2 <- drop1(m2.Cp.Surv, test="Chisq")  
d1
d2
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))
summary(m1.Cp.Surv)   



#----  Px   

#female pupal weights
d.temp <- filter (d.rich, Species=="Px" & Treatment != "C" & Sex=="f")
m1.Px.PW.f <- lmer(Pupal.weight.ST ~ SD*Richness + (1|Treatment) + (1|exp), data = d.temp)
m2.Px.PW.f <- lmer(Pupal.weight.ST ~ SD + Richness + (1|Treatment)+ (1|exp), data = d.temp)
d1 <- drop1(m1.Px.PW.f, test="Chisq")
d2 <- drop1(m2.Px.PW.f, test="Chisq")  #no effects
d1
d2
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))

#male pupal weights
d.temp <- filter (d.rich, Species=="Px" & Treatment != "C" & Sex=="m")
m1.Px.PW.m <- lmer(Pupal.weight.ST ~ SD*Richness + (1|Treatment) + (1|exp), data = d.temp)
m2.Px.PW.m <- lmer(Pupal.weight.ST ~ SD + Richness + (1|Treatment) + (1|exp), data = d.temp)
d1 <- drop1(m1.Px.PW.m, test="Chisq")
d2 <- drop1(m2.Px.PW.m, test="Chisq")  #no effects
d1
d2
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))

#female days to pupation
d.temp <- filter (d.rich, Species=="Px" & Treatment != "C" & Sex=="f")
m1.Px.dTp.f <- lmer(Days.to.pupation.ST.inv ~ SD*Richness + (1|Treatment) + (1|exp), data = d.temp)
m2.Px.dTp.f <- lmer(Days.to.pupation.ST.inv ~ SD + Richness + (1|Treatment) + (1|exp), data = d.temp)
d1 <- drop1(m1.Px.dTp.f, test="Chisq")
d2 <- drop1(m2.Px.dTp.f, test="Chisq")  #no effects
d1
d2
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))


#male days to pupation
d.temp <- filter (d.rich, Species=="Px" & Treatment != "C" & Sex=="m")
m1.Px.dTp.m <- lmer(Days.to.pupation.ST.inv ~ SD*Richness + (1|Treatment) + (1|exp), data = d.temp)
m2.Px.dTp.m <- lmer(Days.to.pupation.ST.inv ~ SD + Richness + (1|Treatment) + (1|exp), data = d.temp)
d1 <- drop1(m1.Px.dTp.m, test="Chisq")
d2 <- drop1(m2.Px.dTp.m, test="Chisq") 
d1
d2
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))
summary(m1.Px.dTp.m)

#Tukey post-hoc for SD levels
m.T <- glht(m2.Px.dTp.m, linfct=mcp(SD="Tukey"))
summary(m.T)
cld(m.T, level = 0.05) #medium lower than low
plot(Days.to.pupation.ST.inv ~ SD, data=d.temp)




#significant SD: richness interaction and significant effect of SD
interaction.plot(d.temp$SD, d.temp$Richness, d.temp$Days.to.pupation.ST.inv)
interaction.plot(d.temp$Richness, d.temp$SD, d.temp$Days.to.pupation.ST.inv)

#separating by SD level
d.temp <- filter (d.rich, Species=="Px" & Treatment != "C" & SD=="L")
m.Px.dTp.SDL <- lmer(Days.to.pupation.ST.inv ~ Richness + (1|Treatment) + (1|exp), data = d.temp)
drop1(m.Px.dTp.SDL, test="Chisq")

d.temp <- filter (d.rich, Species=="Px" & Treatment != "C" & SD=="M")
m.Px.dTp.SDM <- lmer(Days.to.pupation.ST.inv ~ Richness + (1|Treatment) + (1|exp), data = d.temp)
drop1(m.Px.dTp.SDM, test="Chisq")

d.temp <- filter (d.rich, Species=="Px" & Treatment != "C" & SD=="H")
m.Px.dTp.SDH <- lmer(Days.to.pupation.ST.inv ~ Richness + (1|Treatment) + (1|exp), data = d.temp)
drop1(m.Px.dTp.SDH, test="Chisq")


#separating by richness level
d.temp <- filter (d.rich, Species=="Px" & Treatment != "C" & Richness==2)
m.Px.dTp.R2 <- lmer(Days.to.pupation.ST.inv ~ SD + (1|Treatment) + (1|exp), data = d.temp)
drop1(m.Px.dTp.R2, test="Chisq")
summary(m.Px.dTp.R2)  #marginal effect of SD
#Tukey post-hoc for SD levels
m.T <- glht(m.Px.dTp.R2, linfct=mcp(SD="Tukey"))
summary(m.T)
cld(m.T, level = 0.05) #marginal difference with medium lower than high
plot(Days.to.pupation.ST.inv ~ SD, data=d.temp)


d.temp <- filter (d.rich, Species=="Px" & Treatment != "C" & Richness==4)
m.Px.dTp.R4 <- lmer(Days.to.pupation.ST.inv ~ SD + (1|Treatment) + (1|exp), data = d.temp)
drop1(m.Px.dTp.R4, test="Chisq")


d.temp <- filter (d.rich, Species=="Px" & Treatment != "C" & Richness==6)
m.Px.dTp.R6 <- lmer(Days.to.pupation.ST.inv ~ SD + (1|Treatment) + (1|exp), data = d.temp)
drop1(m.Px.dTp.R6, test="Chisq")


d.temp <- filter (d.rich, Species=="Px" & Treatment != "C" & Richness==8)
m.Px.dTp.R8 <- lmer(Days.to.pupation.ST.inv ~ SD + (1|Treatment) + (1|exp), data = d.temp)
drop1(m.Px.dTp.R8, test="Chisq")


d.temp <- filter (d.rich, Species=="Px" & Treatment != "C" & Richness==10)
m.Px.dTp.R10 <- lmer(Days.to.pupation.ST.inv ~ SD + (1|Treatment) + (1|exp), data = d.temp)
drop1(m.Px.dTp.R10, test="Chisq")


#survival
d.temp <- filter (d.rich.PS, Species=="Px" & Treatment != "C")
m1.Px.Surv <- lmer(PropSurv.ST ~ SD*Richness + (1|exp), data = d.temp)
m2.Px.Surv <- lmer(PropSurv.ST ~ SD+Richness + (1|exp), data = d.temp)
d1 <- drop1(m1.Px.Surv, test="Chisq")
d2 <- drop1(m2.Px.Surv, test="Chisq")
d1
d2
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))


#---For Fungi

#---- all fungi combined

hist(d.fungi2$dAbs.ST)

#ALL fungal species together
d.temp <- filter (d.fungi2, Treatment != "DMSO")
m1.all.f <- lmer(dAbs.ST ~ SD*Richness*Fungi + (1|Treatment) + (1|Plate), 
                   data = d.temp, na.action=na.fail, REML="FALSE")
d1.all.f <- dredge(m1.all.f)
d1.all.f  #top model contains all factors and dAIC to next best is 36.12, retain all
summary(d1.all.f)

d1.all.avg <- model.avg(d1.all.f, fit=TRUE)

summary(d1.all.avg) #use for S4

m1.top.f <- lmer(dAbs.ST ~  SD*Richness*Fungi + (1|Treatment) + (1|Plate), 
                   data = d.temp, na.action=na.fail) 
summary(m1.top.f)
drop1(m1.top.f, test="Chisq")  #highly significant 3 way interaction

plot(d.temp$dAbs.ST ~ d.temp$Richness)
abline(lm(d.temp$dAbs.ST ~ d.temp$Richness))

#need to look at fungal species separately


#---- Botrys

d.temp <- filter (d.fungi2, Fungi=="Botrys" & Treatment != "DMSO")
m1.B <- lmer(dAbs.ST ~ SD*Richness + (1|Treatment) + (1|Plate), data = d.temp)
m2.B <- lmer(dAbs.ST ~ SD + Richness + (1|Treatment) + (1|Plate), data = d.temp)
d1 <- drop1(m1.B, test="Chisq")
d2 <- drop1(m2.B, test="Chisq")  
d1   #significant richness*SD interaction
d2   
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4])) 

hist(d.temp$dAbs.ST)
plot(m1.B)
hist(resid(m1.B))

plot(dAbs.ST ~ Richness, data=d.temp) 
interaction.plot(d.temp$Richness, d.temp$SD, d.temp$dAbs.ST)
#seems like a negative effect of richness on growth, but only at 
#med levels, and maybe a positive effect at low levels

#separate by SD

#--- low
d.temp <- filter (d.fungi2, Fungi=="Botrys" & Treatment != "DMSO" & SD=="L")
m1.B.L <- lmer(dAbs.ST ~ Richness + (1|Treatment) + (1|Plate), data = d.temp)
summary(m1.B.L)
d1 <- drop1(m1.B.L, test="Chisq")
d1   #significant positive effect of richness at low SD

#--- med
d.temp <- filter (d.fungi2, Fungi=="Botrys" & Treatment != "DMSO" & SD=="M")
m1.B.M <- lmer(dAbs.ST ~ Richness + (1|Treatment), data = d.temp)
summary(m1.B.M)
d1 <- drop1(m1.B.M, test="Chisq")
d1   #marginal negative effect of richness at med SD

#--- high
d.temp <- filter (d.fungi2, Fungi=="Botrys" & Treatment != "DMSO" & SD=="H")
m1.B.H <- lmer(dAbs.ST ~ Richness + (1|Treatment), data = d.temp)
summary(m1.B.H)
d1 <- drop1(m1.B.H, test="Chisq")
d1   #no effect of richness at high SD


#---- Colletotrichum 

d.temp <- filter (d.fungi2, Fungi=="Collet" & Treatment != "DMSO")
m1.B <- lmer(dAbs.ST ~ SD*Richness + (1|Treatment) + (1|Plate), data = d.temp)  
m2.B <- lmer(dAbs.ST ~ SD + Richness + (1|Treatment)+ (1|Plate), data = d.temp)
d1 <- drop1(m1.B, test="Chisq")
d2 <- drop1(m2.B, test="Chisq")  #no effects
d1
d2
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4])) 

hist(d.temp$dAbs.ST)
plot(m1.B)
hist(resid(m1.B))

plot(dAbs.ST ~ Richness, data=d.temp)  
abline(lm(dAbs.ST ~ Richness, data=d.temp))

#---- Penicillium  

d.temp <- filter (d.fungi2, Fungi=="Penicillium" & Treatment != "DMSO")
m1.B <- lmer(dAbs.ST ~ SD*Richness + (1|Treatment)+ (1|Plate), data = d.temp)
m2.B <- lmer(dAbs.ST ~ SD + Richness + (1|Treatment)+ (1|Plate), data = d.temp)
d1 <- drop1(m1.B, test="Chisq")
d2 <- drop1(m2.B, test="Chisq")  
d1
d2   #no effects
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4])) 

hist(d.temp$dAbs.ST)   ##left skewed
plot(m1.B)  
hist(resid(m1.B))

plot(dAbs.ST ~ Richness, data=d.temp) 
abline(lm(dAbs.ST ~ Richness, data=d.temp))

#---- Sclerotinia 
d.temp <- filter (d.fungi2, Fungi=="Sclerotinia" & Treatment != "DMSO")
m1.B <- lmer(dAbs.ST ~ SD*Richness + (1|Treatment)+ (1|Plate), data = d.temp)
m2.B <- lmer(dAbs.ST ~ SD + Richness + (1|Treatment)+ (1|Plate), data = d.temp)
d1 <- drop1(m1.B, test="Chisq")
d2 <- drop1(m2.B, test="Chisq")  
d1  #marginal interaction
d2   #no effects
m.sum <- rbind(m.sum, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4])) 

hist(d.temp$dAbs.ST)
plot(m1.B)
hist(resid(m1.B))

plot(dAbs.ST ~ Richness, data=d.temp)
abline(lm(dAbs.ST ~ Richness, data=d.temp))  #slightly positive
interaction.plot(d.temp$Richness, d.temp$SD, d.temp$dAbs.ST)
#again, interaction seems like it is all over...slightly negative at med, 
#slightly positive at low


#Save results table
m.sum
write.csv(m.sum, "./Outputs/Tables/TableS5_Pred1a-1c_EffectsBySpecies.csv")

#Save workspace
save.image("./Outputs/Workspaces/FinalAnalyses_Pred1a-1c")


#**************************************************************




#---Prediction 1a/1b/1c null models----
#----

load("./Outputs/Workspaces/StandardizedData")

#------For Insects

#---all species combined

#make table to store stats across all iterations
m.PW.null <- data.frame(SDL=numeric())

for (i in 1:100){
  #Looking at effects of Richness, SD, Species, and Sex on pupal weights
  d.temp <- filter (d.rich, Treatment != "C", !is.na(Sex), !is.na(Pupal.weight.ST))
  
  #randomizing data
  d.temp <- d.temp %>% 
    ungroup() %>%
    mutate(Pupal.weight.ST = Pupal.weight.ST[sample(row_number())]) 
  
  m1.all.PW <- lmer(Pupal.weight.ST ~ SD*Richness*Species*Sex + (1|Treatment) + (1|exp/TrayID/CupID), 
                    data = d.temp, na.action=na.fail, REML="FALSE")
  d1.all.PW <- dredge(m1.all.PW)
  s <- subset(d1.all.PW, delta < 4)
  s  #use for table S3
  d1.all.PW.avg <- model.avg(d1.all.PW, subset = delta < 4, fit=TRUE)
  sum <- summary(d1.all.PW.avg)  #use for table S4
  newrow <- data.frame(coef=rownames(sum$coefmat.full), P=sum$coefmat.full[,5])
  newrow <- pivot_wider(newrow, names_from = coef, values_from = P)
  m.PW.null <- full_join(m.PW.null, newrow)
}

T1 <- function(x){
  l <- length(x[which(x < 0.05)])
  na <- length(x[which(!is.na(x))])
  ifelse(na >= 10, return(l/na), return (NA))
}

T1rates <- m.PW.null %>%
  summarise_all(list(T1rate=T1)) 



#Effects on days to pupation
m.DtP.null <- data.frame(SDL=numeric())

for (i in 1:100){
  #Looking at effects of Richness, SD, Species, and Sex on DtP
  d.temp <- filter (d.rich, Treatment != "C", !is.na(Sex), !is.na(Days.to.pupation.ST.inv))
  
  #randomizing data
  d.temp <- d.temp %>% 
    ungroup() %>%
    mutate(Days.to.pupation.ST.inv = Days.to.pupation.ST.inv[sample(row_number())]) 
  
  m1.all.DtP <- lmer(Days.to.pupation.ST.inv ~ SD*Richness*Species*Sex + (1|Treatment) + (1|exp/TrayID/CupID), 
                    data = d.temp, na.action=na.fail, REML="FALSE")
  d1.all.DtP <- dredge(m1.all.DtP)
  s <- subset(d1.all.DtP, delta < 4)
  s  #use for table S3
  d1.all.DtP.avg <- model.avg(d1.all.DtP, subset = delta < 4, fit=TRUE)
  sum <- summary(d1.all.DtP.avg)  #use for table S4
  newrow <- data.frame(coef=rownames(sum$coefmat.full), P=sum$coefmat.full[,5])
  newrow <- pivot_wider(newrow, names_from = coef, values_from = P)
  m.DtP.null <- full_join(m.DtP.null, newrow)
}

T1 <- function(x){
  l <- length(x[which(x < 0.05)])
  na <- length(x[which(!is.na(x))])
  ifelse(na >= 10, return(l/na), return (NA))
}

T1rates.DtP <- m.DtP.null %>%
  summarise_all(list(T1rate=T1)) 



##Effects on survival

m.S.null <- data.frame(SDL=numeric())

for (i in 1:100){

  d.temp <- filter (d.rich.PS, Treatment != "C", !is.na(PropSurv.ST))
  
  #randomizing data
  d.temp <- d.temp %>% 
    ungroup() %>%
    mutate(PropSurv.ST = PropSurv.ST[sample(row_number())]) 
  
  
  m1.all.S <- lmer(PropSurv.ST ~ SD*Richness*Species + (1|Treatment) + (1|exp), 
                   data = d.temp, na.action=na.fail, REML="FALSE")
  d1.all.S <- dredge(m1.all.S)
  s <- subset(d1.all.S, delta < 4)
  s #use for table S3
  d1.all.S.avg <- model.avg(d1.all.S, subset =  delta < 4, fit=TRUE)
  sum <- summary(d1.all.S.avg)  #use for table S4
  newrow <- data.frame(coef=rownames(sum$coefmat.full), P=sum$coefmat.full[,5])
  newrow <- pivot_wider(newrow, names_from = coef, values_from = P)
  m.S.null <- full_join(m.S.null, newrow)

}  

T1rates.S <- m.S.null %>%
  summarise_all(list(T1rate=T1))   


#---For Fungi

#---- all fungi combined

m.f.null <- data.frame(SDL=numeric())

for (i in 1:100){

  #ALL fungal species together
  d.temp <- filter (d.fungi2, Treatment != "DMSO")
  
  #randomizing data
  d.temp <- d.temp %>% 
    ungroup() %>%
    mutate(dAbs.ST = dAbs.ST[sample(row_number())]) 
  
  m1.all.f <- lmer(dAbs.ST ~ SD*Richness*Fungi + (1|Treatment) + (1|Plate), 
                   data = d.temp, na.action=na.fail, REML="FALSE")
  d1.all.f <- dredge(m1.all.f)
  d1.all.f  #top model contains all factors and dAIC to next best is 36.12, retain all
  summary(d1.all.f)
  
  d1.all.avg <- model.avg(d1.all.f, fit=TRUE)
  
  sum <- summary(d1.all.avg) #use for S4
  newrow <- data.frame(coef=rownames(sum$coefmat.full), P=sum$coefmat.full[,5])
  newrow <- pivot_wider(newrow, names_from = coef, values_from = P)
  m.f.null <- full_join(m.f.null, newrow)
}

T1rates.f <- m.f.null %>%
  summarise_all(list(T1rate=T1)) 

save.image("./Outputs/Workspaces/FinalAnalyses_Pred1a-1c_nullmodels")


# Testing whether synergy increases with richness (Prediction 1d)----
#--------------------

load("./Outputs/Workspaces/StandardizedData")

#The basic idea here is to compare the observed effect to the expected based on 
#additive interactions among compounds. Methods are based on Tallarida 2000

#First need values for the "observed effect", aka the "specified level of effect" or the 
#  "Z" in Tallarida
#Since we do not have dose response curves and ED50 values, this will just be the 
#effect size estimate in comparisons between each treatment and the controls 
#from the same experiment

##Separate analyses for pupal weight and for days to pupation

#make table to store stats for richness vs synergy tests
m.sum.RvsS <- data.frame(Metric=NA, Sp=NA, Slope=NA, Z=NA, P=NA)

#for pupal weight

d.temp <- d.rich

Zs.PW <- data.frame(Species=factor(), Treatment=factor(), Zt=numeric(), Zt.SE=numeric())
for(i in 1:length(levels(d.temp$Treatment))){
  if(levels(d.temp$Treatment)[i] != "C"){
    Tx <- levels(d.temp$Treatment)[i] 
    d <- d.temp[which(d.temp$Treatment==Tx),]
    e <- as.character(unique(d$exp))
    d2 <- d.temp[which(d.temp$Treatment==Tx | d.temp$Treatment=="C"),]
    d2 <- filter(d2, exp %in% e)
    
    for(j in 1:length(levels(d.temp$Species))){
      Sp <- levels(d.temp$Species)[j]
      d3 <- d2[which(d2$Species==Sp),]
      d3$Treatment <- factor(d3$Treatment, levels=c("C", Tx))
      m.PW <- lm(Pupal.weight ~ Treatment, data=d3)
      b <- coef(summary(m.PW))[2,1]
      b.SE <- coef(summary(m.PW))[2,2]
      newrow = data.frame(Species=Sp, Treatment=Tx, Zt=b, Zt.SE=b.SE)
      Zs.PW <- rbind(Zs.PW, newrow)
    }
  }
}


##Now alculate Zadd and variance for each mixture
#from Tallarida, the variance of Zadd is
#var_Zadd= f^2 * V_A + (1-f)^2 * V_B
#where f=the fraction of each compound in mix


Zs.PW$Zadd <- NA
Zs.PW$Zadd.SE <- NA

#To calculate Zadd, need to know the conc of each compound in each treatment.
#Zs.PW <- left_join(Zs.PW, distinct(d.temp[,c(2, 12:25)], Treatment, .keep_all=TRUE), by="Treatment")

Tx.Comp <- distinct(d.temp[,c(2, 12:25)], Treatment, .keep_all=TRUE)
Tx.Comp$Treatment <- as.character(Tx.Comp$Treatment)

for(i in 1:length(Zs.PW$Zadd)){
  Zs.PW$Zadd[i] <- 0
  Zs.PW$Zadd.SE[i] <- 0
  H <- Zs.PW$Species[i]
  Tx <- Zs.PW$Treatment[i]
  for(j in 2:15){
    C.ID <- colnames(Tx.Comp)[j]
    e <- Zs.PW$Zt[which(Zs.PW$Treatment==C.ID & Zs.PW$Species==H)]  #effect of indiv compound
    var <- Zs.PW$Zt.SE[which(Zs.PW$Treatment==C.ID & Zs.PW$Species==H)]^2  ##var of slope from model coef (se^2)  see: https://stackoverflow.com/questions/14960868/how-do-i-print-the-variance-of-an-lm-in-r-without-computing-from-the-standard-er
    a <- e*(Tx.Comp[which(Tx.Comp$Treatment==Tx),j]/0.5)  ##effect * relative abundance (note total compounds=0.5mg, so dividing by that)
    #r <- d.temp$Richness[i]
    b <- var*((Tx.Comp[which(Tx.Comp$Treatment==Tx),j]/0.5))  ##in Tallarida equation the proportions are squared, 
    #but that gives really weird results with multiple compounds, where the variance estimate gets much lower as you increase in richness
    ##here I just made the variance a straight additive proportional function of the composition
    Zs.PW$Zadd[i] <- Zs.PW$Zadd[i] + a
    Zs.PW$Zadd.SE[i] <- Zs.PW$Zadd.SE[i] + b  #this is a actually generating the variance
  }
  Zs.PW$Zadd.SE [i] <- sqrt(Zs.PW$Zadd.SE[i]) #square root to get SE
}


# Define Richness: 0,1,2,4,6,8,10
Zs.PW$Richness <- NA
for (i in 1:nrow(Zs.PW)){
  Zs.PW$Richness[i] <- str_extract(Zs.PW$Treatment[i],'[0-9]+')}
Zs.PW$Richness[which(is.na(Zs.PW$Richness))] <- 1
Zs.PW$Richness <- as.numeric(Zs.PW$Richness)


#Calculate CIs for Zt and Zadd
#From Tallarida 2002, total sample size for calculating CI is based on total
#number of datapoints in drug1, drug 2, and combination all combined
#So in our case we have N=15 for each compound, plus ~N=5-8 for each mixture
#so use t-value 1.960 (for large sample size, anything over df=120)
#So 95% CI is Zadd +/- int

Zs.PW$Zt_upper <- Zs.PW$Zt + (Zs.PW$Zt.SE*1.96)
Zs.PW$Zt_lower <- Zs.PW$Zt - (Zs.PW$Zt.SE*1.96)
Zs.PW$Zadd_upper <- Zs.PW$Zadd + (Zs.PW$Zadd.SE*1.96)
Zs.PW$Zadd_lower <- Zs.PW$Zadd - (Zs.PW$Zadd.SE*1.96)



#Plots

#divide data
z1 <- Zs.PW[,c(1:4,7:9)]
z2 <- Zs.PW[,c(1:2, 5:7, 10:11)]
z1$Tx <- "obs"
z2$Tx <- "exp"
colnames(z1) <- c("Species", "Treatment", "Z", "SE", "Richness", "upper",  "lower", "Tx")
colnames(z2) <- c("Species", "Treatment", "Z", "SE", "Richness", "upper",  "lower", "Tx")

z3 <- rbind(z1, z2)

z4 <- z3[which(z3$Richness!=1),]
z4 <- droplevels(z4)
ggplot(z4[which(z4$Species=="Hz"),], aes(x=Treatment, y = Z, group = Tx)) +
  geom_point( aes(color=Tx)) +
  geom_errorbar(aes(ymin=lower, ymax=upper, color=Tx), width=.5)
ggplot(z4[which(z4$Species=="Sf"),], aes(x=Treatment, y = Z, group = Tx)) +
  geom_point( aes(color=Tx)) +
  geom_errorbar(aes(ymin=lower, ymax=upper, color=Tx), width=.5)
ggplot(z4[which(z4$Species=="Cp"),], aes(x=Treatment, y = Z, group = Tx)) +
  geom_point( aes(color=Tx)) +
  geom_errorbar(aes(ymin=lower, ymax=upper, color=Tx), width=.5)
ggplot(z4[which(z4$Species=="Px"),], aes(x=Treatment, y = Z, group = Tx)) +
  geom_point( aes(color=Tx)) +
  geom_errorbar(aes(ymin=lower, ymax=upper, color=Tx), width=.5)

#Note confidence intervals overlap for every single one!!


#or plot diff exp-obs
#0= no non-additive effects
#>0 indicates synergy
#<0 indicates antagonism
#positive effect of richness on Zdiff would indicate increasing synergy with more complex mixtures
##non overlapping CI would indicate significant non-additive effects

Zs.PW$diff <- Zs.PW$Zadd-Zs.PW$Zt   
plot(Zs.PW$diff~ Zs.PW$Richness)
abline(lm(Zs.PW$diff~ Zs.PW$Richness))


##Can test for increase in synergy with increasing richness

m1.Hz <- summary(lm(diff~ Richness, data=Zs.PW[which(Zs.PW$Species=="Hz" & Zs.PW$Richness !=1),]))  #no synergy
m1.Sf <- summary(lm(diff~ Richness, data=Zs.PW[which(Zs.PW$Species=="Sf" & Zs.PW$Richness !=1),]))  #no synergy
m1.Px <- summary(lm(diff~ Richness, data=Zs.PW[which(Zs.PW$Species=="Px" & Zs.PW$Richness !=1),]))  #marginal antagonism
m1.Cp <- summary(lm(diff~ Richness, data=Zs.PW[which(Zs.PW$Species=="Cp" & Zs.PW$Richness !=1),]))  #no synergy
m1.Hz
m1.Sf
m1.Cp
m1.Px


#marginally significant effect for Px; checking plot
plot(diff~ Richness, data=Zs.PW[which(Zs.PW$Species=="Px" & Zs.PW$Richness !=1),])
abline(lm(diff~ Richness, data=Zs.PW[which(Zs.PW$Species=="Px" & Zs.PW$Richness !=1),]))


#Add stats to summary table
m.sum.RvsS[1,] <- c("PW", "Hz", m1.Hz$coefficients[2,1], m1.Hz$coefficients[2,3], m1.Hz$coefficients[2,4])
m.sum.RvsS <- rbind(m.sum.RvsS,
                    c("PW", "Sf", m1.Sf$coefficients[2,1], m1.Sf$coefficients[2,3], m1.Sf$coefficients[2,4]))
m.sum.RvsS <- rbind(m.sum.RvsS,
                    c("PW", "Cp", m1.Cp$coefficients[2,1], m1.Cp$coefficients[2,3], m1.Cp$coefficients[2,4]))
m.sum.RvsS <- rbind(m.sum.RvsS,
                    c("PW", "Px", m1.Px$coefficients[2,1], m1.Px$coefficients[2,3], m1.Px$coefficients[2,4]))





####Now the same for DtP

d.temp <- d.rich

#taking the inverse of days to pupation--i.e. "development speed" so that the 
#directionality of effects is the same for all performance metrics
#d.temp$Days.to.pupation.inv <- 1/d.temp$Days.to.pupation
#hist(d.temp$Days.to.pupation.inv)

Zs.DtP <- data.frame(Species=factor(), Treatment=factor(), Zt=numeric(), Zt.SE=numeric())
for(i in 1:length(levels(d.temp$Treatment))){
  if(levels(d.temp$Treatment)[i] != "C"){
    Tx <- levels(d.temp$Treatment)[i] 
    d <- d.temp[which(d.temp$Treatment==Tx),]
    d2 <- d.temp[which(d.temp$Treatment==Tx | d.temp$Treatment=="C"),]
    
    for(j in 1:length(levels(d.temp$Species))){
      Sp <- levels(d.temp$Species)[j]
      d3 <- d2[which(d2$Species==Sp),]
      d3$Treatment <- factor(d3$Treatment, levels=c(Tx, "C"))
        #switched order of Tx and C here so that the directionality of the effect will be the same
        #across performance metrics. This will give the slope of the line from Tx to C
      m.DtP <- lm(Days.to.pupation ~ Treatment, data=d3)
      b <- coef(summary(m.DtP))[2,1]
      b.SE <- coef(summary(m.DtP))[2,2]
      newrow = data.frame(Species=Sp, Treatment=Tx, Zt=b, Zt.SE=b.SE)
      Zs.DtP <- rbind(Zs.DtP, newrow)
    }
  }
}

#**Tried doing these analyses comparing data for each treatment only to the controls from the 
#same experiment, but the models were getting "essentially perfect fit, may be unreliable" errors
#because in several cases all the controls for an experiment (7-9 caterpillars) pupated on the same day
#and so there was no variance. Thus, I compared to all the controls for that species for DtP, which
#increased the C sample size and I think provides a more clear picture 



##Now calculate Zadd and variance for each mixture
#from Tallarida, the variance of Zadd is
#var_Zadd= f^2 * V_A + (1-f)^2 * V_B
#where f=the fraction of each compound in mix


Zs.DtP$Zadd <- NA
Zs.DtP$Zadd.SE <- NA

Tx.Comp <- distinct(d.temp[,c(2, 12:25)], Treatment, .keep_all=TRUE)
Tx.Comp$Treatment <- as.character(Tx.Comp$Treatment)

for(i in 1:length(Zs.DtP$Zadd)){
  Zs.DtP$Zadd[i] <- 0
  Zs.DtP$Zadd.SE[i] <- 0
  H <- Zs.DtP$Species[i]
  Tx <- Zs.DtP$Treatment[i]
  for(j in 2:15){
    C.ID <- colnames(Tx.Comp)[j]
    e <- Zs.DtP$Zt[which(Zs.DtP$Treatment==C.ID & Zs.DtP$Species==H)]  #effect of indiv compound
    var <- Zs.DtP$Zt.SE[which(Zs.DtP$Treatment==C.ID & Zs.DtP$Species==H)]^2  ##var of slope from model coef (se^2)  see: https://stackoverflow.com/questions/14960868/how-do-i-print-the-variance-of-an-lm-in-r-without-computing-from-the-standard-er
    a <- e*(Tx.Comp[which(Tx.Comp$Treatment==Tx),j]/0.5)  ##effect * relative abundance (note total compounds=0.5mg, so dividing by that)
    #r <- d.temp$Richness[i]
    b <- var*((Tx.Comp[which(Tx.Comp$Treatment==Tx),j]/0.5))  ##in Tallarida equation the proportions are squared, 
    #but that gives really weird results with multiple compounds, where the variance estimate gets much lower as you increase in richness
    ##here I just made the variance a straight additive proportional function of the composition
    Zs.DtP$Zadd[i] <- Zs.DtP$Zadd[i] + a
    Zs.DtP$Zadd.SE[i] <- Zs.DtP$Zadd.SE[i] + b  #this is a actually generating the variance
  }
  Zs.DtP$Zadd.SE [i] <- sqrt(Zs.DtP$Zadd.SE[i]) #square root to get SE
}


# Define Richness: 0,1,2,4,6,8,10
Zs.DtP$Richness <- NA
for (i in 1:nrow(Zs.DtP)){
  Zs.DtP$Richness[i] <- str_extract(Zs.DtP$Treatment[i],'[0-9]+')}
Zs.DtP$Richness[which(is.na(Zs.DtP$Richness))] <- 1
Zs.DtP$Richness <- as.numeric(Zs.DtP$Richness)


#Calculate CIs for Zt and Zadd
#From Tallarida 2002, total sample size for calculating CI is based on total
#number of datapoints in drug1, drug 2, and combination all combined
#So in our case we have N=15 for each compound, plus ~N=5-8 for each mixture
#so use t-value 1.960 (for large sample size, anything over df=120)
#So 95% CI is Zadd +/- int

Zs.DtP$Zt_upper <- Zs.DtP$Zt + (Zs.DtP$Zt.SE*1.96)
Zs.DtP$Zt_lower <- Zs.DtP$Zt - (Zs.DtP$Zt.SE*1.96)
Zs.DtP$Zadd_upper <- Zs.DtP$Zadd + (Zs.DtP$Zadd.SE*1.96)
Zs.DtP$Zadd_lower <- Zs.DtP$Zadd - (Zs.DtP$Zadd.SE*1.96)




#Plots

#divide data
z1 <- Zs.DtP[,c(1:4,7:9)]
z2 <- Zs.DtP[,c(1:2, 5:7, 10:11)]
z1$Tx <- "obs"
z2$Tx <- "exp"
colnames(z1) <- c("Species", "Treatment", "Z", "SE", "Richness", "upper",  "lower", "Tx")
colnames(z2) <- c("Species", "Treatment", "Z", "SE", "Richness", "upper",  "lower", "Tx")

z3 <- rbind(z1, z2)

z4 <- z3[which(z3$Richness!=1),]
z4 <- droplevels(z4)
ggplot(z4[which(z4$Species=="Hz"),], aes(x=Treatment, y = Z, group = Tx)) +
  geom_point( aes(color=Tx)) +
  geom_errorbar(aes(ymin=lower, ymax=upper, color=Tx), width=.5)
ggplot(z4[which(z4$Species=="Sf"),], aes(x=Treatment, y = Z, group = Tx)) +
  geom_point( aes(color=Tx)) +
  geom_errorbar(aes(ymin=lower, ymax=upper, color=Tx), width=.5)
ggplot(z4[which(z4$Species=="Cp"),], aes(x=Treatment, y = Z, group = Tx)) +
  geom_point( aes(color=Tx)) +
  geom_errorbar(aes(ymin=lower, ymax=upper, color=Tx), width=.5)
ggplot(z4[which(z4$Species=="Px"),], aes(x=Treatment, y = Z, group = Tx)) +
  geom_point( aes(color=Tx)) +
  geom_errorbar(aes(ymin=lower, ymax=upper, color=Tx), width=.5)

#Some effects of Hz, Sf, and Px; in most cases values for observed are higher


#or plot diff exp-obs

#0= no non-additive effects
#>0 indicates synergy
#<0 indicates antagonism
#positive effect of richness on Zdiff would indicate increasing synergy with more complex mixtures
##non overlapping CI would indicate significant non-additive effects

Zs.DtP$diff <- Zs.DtP$Zadd-Zs.DtP$Zt   
plot(Zs.DtP$diff~ Zs.DtP$Richness)
abline(lm(Zs.DtP$diff~ Zs.DtP$Richness))


##Can test for increase in synergy with increasing richness

m1.Hz <- summary(lm(diff~ Richness, data=Zs.DtP[which(Zs.DtP$Species=="Hz" & Zs.DtP$Richness !=1),]))  #no synergy
m1.Sf <- summary(lm(diff~ Richness, data=Zs.DtP[which(Zs.DtP$Species=="Sf"& Zs.DtP$Richness !=1),]))  #no synergy
m1.Px <- summary(lm(diff~ Richness, data=Zs.DtP[which(Zs.DtP$Species=="Px"& Zs.DtP$Richness !=1),]))  #no synergy
m1.Cp <- summary(lm(diff~ Richness, data=Zs.DtP[which(Zs.DtP$Species=="Cp"& Zs.DtP$Richness !=1),]))  #no synergy
m1.Hz
m1.Sf
m1.Px
m1.Cp

#Add stats to summary table
m.sum.RvsS <- rbind(m.sum.RvsS, 
                    c("DtP", "Hz", m1.Hz$coefficients[2,1], m1.Hz$coefficients[2,3], m1.Hz$coefficients[2,4]))
m.sum.RvsS <- rbind(m.sum.RvsS,
                    c("DtP", "Sf", m1.Sf$coefficients[2,1], m1.Sf$coefficients[2,3], m1.Sf$coefficients[2,4]))
m.sum.RvsS <- rbind(m.sum.RvsS,
                    c("DtP", "Cp", m1.Cp$coefficients[2,1], m1.Cp$coefficients[2,3], m1.Cp$coefficients[2,4]))
m.sum.RvsS <- rbind(m.sum.RvsS,
                    c("DtP", "Px", m1.Px$coefficients[2,1], m1.Px$coefficients[2,3], m1.Px$coefficients[2,4]))



#same for fungi

d.temp <- d.fungi2
d.temp$Treatment <- as.factor(d.temp$Treatment)

Zs.F <- data.frame(Species=factor(), Treatment=factor(), Zt=numeric(), Zt.SE=numeric())
for(i in 1:length(levels(d.temp$Treatment))){
  if(levels(d.temp$Treatment)[i] != "DMSO"){
    Tx <- levels(d.temp$Treatment)[i] 
    d <- d.temp[which(d.temp$Treatment==Tx),]
    #e <- as.character(unique(d$exp)) #**see note below
    d2 <- d.temp[which(d.temp$Treatment==Tx | d.temp$Treatment=="DMSO"),]
    #d2 <- filter(d2, exp %in% e)
    
    for(j in 1:length(levels(d.temp$Fungi))){
      Sp <- levels(d.temp$Fungi)[j]
      d3 <- d2[which(d2$Fungi==Sp),]
      d3$Treatment <- factor(d3$Treatment, levels=c("DMSO", Tx))
      m.F <- lm(slope ~ Treatment, data=d3)
      b <- coef(summary(m.F))[2,1]
      b.SE <- coef(summary(m.F))[2,2]
      newrow = data.frame(Species=Sp, Treatment=Tx, Zt=b, Zt.SE=b.SE)
      Zs.F <- rbind(Zs.F, newrow)
    }
  }
}


##Now calculate Zadd and variance for each mixture
#from Tallarida, the variance of Zadd is
#var_Zadd= f^2 * V_A + (1-f)^2 * V_B
#where f=the fraction of each compound in mix


Zs.F$Zadd <- NA
Zs.F$Zadd.SE <- NA

for(i in 1:length(Zs.F$Zadd)){
  Zs.F$Zadd[i] <- 0
  Zs.F$Zadd.SE[i] <- 0
  H <- Zs.F$Species[i]
  Tx <- Zs.F$Treatment[i]
  for(j in 2:15){
    C.ID <- colnames(Tx.Comp)[j]
    e <- Zs.F$Zt[which(Zs.F$Treatment==C.ID & Zs.F$Species==H)]  #effect of indiv compound
    var <- Zs.F$Zt.SE[which(Zs.F$Treatment==C.ID & Zs.F$Species==H)]^2  ##var of slope from model coef (se^2)  see: https://stackoverflow.com/questions/14960868/how-do-i-print-the-variance-of-an-lm-in-r-without-computing-from-the-standard-er
    a <- e*(Tx.Comp[which(Tx.Comp$Treatment==Tx),j]/0.5)  ##effect * relative abundance (note total compounds=0.5mg, so dividing by that)
    #r <- d.temp$Richness[i]
    b <- var*((Tx.Comp[which(Tx.Comp$Treatment==Tx),j]/0.5))  ##in Tallarida equation the proportions are squared, 
    #but that gives really weird results with multiple compounds, where the variance estimate gets much lower as you increase in richness
    ##here I just made the variance a straight additive proportional function of the composition
    Zs.F$Zadd[i] <- Zs.F$Zadd[i] + a
    Zs.F$Zadd.SE[i] <- Zs.F$Zadd.SE[i] + b  #this is actually generating the variance
  }
  Zs.F$Zadd.SE [i] <- sqrt(Zs.F$Zadd.SE[i]) #square root to get SE
}


# Define Richness: 0,1,2,4,6,8,10
Zs.F$Richness <- NA
for (i in 1:nrow(Zs.F)){
  Zs.F$Richness[i] <- str_extract(Zs.F$Treatment[i],'[0-9]+')}
Zs.F$Richness[which(is.na(Zs.F$Richness))] <- 1
Zs.F$Richness <- as.numeric(Zs.F$Richness)


#Calculate CIs for Zt and Zadd
#From Tallarida 2002, total sample size for calculating CI is based on total
#number of datapoints in drug1, drug 2, and combination all combined
#So in our case we have N=15 for each compound, plus ~N=5-8 for each mixture
#so use t-value 1.960 (for large sample size, anything over df=120)
#So 95% CI is Zadd +/- int

Zs.F$Zt_upper <- Zs.F$Zt + (Zs.F$Zt.SE*1.96)
Zs.F$Zt_lower <- Zs.F$Zt - (Zs.F$Zt.SE*1.96)
Zs.F$Zadd_upper <- Zs.F$Zadd + (Zs.F$Zadd.SE*1.96)
Zs.F$Zadd_lower <- Zs.F$Zadd - (Zs.F$Zadd.SE*1.96)




#Plots

#divide data
z1 <- Zs.F[,c(1:4,7:9)]
z2 <- Zs.F[,c(1:2, 5:7, 10:11)]
z1$Tx <- "obs"
z2$Tx <- "exp"
colnames(z1) <- c("Species", "Treatment", "Z", "SE", "Richness", "upper",  "lower", "Tx")
colnames(z2) <- c("Species", "Treatment", "Z", "SE", "Richness", "upper",  "lower", "Tx")

z3 <- rbind(z1, z2)

z4 <- z3[which(z3$Richness!=1),]
z4 <- droplevels(z4)
ggplot(z4[which(z4$Species=="Botrys"),], aes(x=Treatment, y = Z, group = Tx)) +
  geom_point( aes(color=Tx)) +
  geom_errorbar(aes(ymin=lower, ymax=upper, color=Tx), width=.5)
ggplot(z4[which(z4$Species=="Collet"),], aes(x=Treatment, y = Z, group = Tx)) +
  geom_point( aes(color=Tx)) +
  geom_errorbar(aes(ymin=lower, ymax=upper, color=Tx), width=.5)
ggplot(z4[which(z4$Species=="Penicillium"),], aes(x=Treatment, y = Z, group = Tx)) +
  geom_point( aes(color=Tx)) +
  geom_errorbar(aes(ymin=lower, ymax=upper, color=Tx), width=.5)
ggplot(z4[which(z4$Species=="Sclerotinia"),], aes(x=Treatment, y = Z, group = Tx)) +
  geom_point( aes(color=Tx)) +
  geom_errorbar(aes(ymin=lower, ymax=upper, color=Tx), width=.5)

#Zs indicate the effect of treatment relative to control. Negative numbers indicate herbivores 
#did worse on treatment. 
#If observed slope is more negative than the expected additive slope, then we have synergy


#or plot diff exp-obs

#0= no non-additive effects
#>0 indicates synergy
#<0 indicates antagonism
#positive effect of richness on Zdiff would indicate increasing synergy with more complex mixtures
##non overlapping CI would indicate significant non-additive effects

Zs.F$diff <- Zs.F$Zadd-Zs.F$Zt   
plot(Zs.F$diff~ Zs.F$Richness)
abline(lm(Zs.F$diff~ Zs.F$Richness))


##Can test for increase in synergy with increasing richness

m1.Bo <- summary(lm(diff~ Richness, data=Zs.F[which(Zs.F$Species=="Botrys" & Zs.F$Richness !=1),]))  #no synergy
m1.Co <- summary(lm(diff~ Richness, data=Zs.F[which(Zs.F$Species=="Collet"& Zs.F$Richness !=1),]))  #no synergy
m1.Pe <- summary(lm(diff~ Richness, data=Zs.F[which(Zs.F$Species=="Penicillium"& Zs.F$Richness !=1),]))  #no synergy
m1.Sc <- summary(lm(diff~ Richness, data=Zs.F[which(Zs.F$Species=="Sclerotinia"& Zs.F$Richness !=1),]))  #no synergy
m1.Bo
m1.Co
m1.Pe
m1.Sc

#Add stats to summary table
m.sum.RvsS <- rbind(m.sum.RvsS, 
                    c("F", "Botrys", m1.Hz$coefficients[2,1], m1.Hz$coefficients[2,3], m1.Hz$coefficients[2,4]))
m.sum.RvsS <- rbind(m.sum.RvsS,
                    c("F", "Collet", m1.Sf$coefficients[2,1], m1.Sf$coefficients[2,3], m1.Sf$coefficients[2,4]))
m.sum.RvsS <- rbind(m.sum.RvsS,
                    c("F", "Penicillium", m1.Cp$coefficients[2,1], m1.Cp$coefficients[2,3], m1.Cp$coefficients[2,4]))
m.sum.RvsS <- rbind(m.sum.RvsS,
                    c("F", "Sclerotinia", m1.Px$coefficients[2,1], m1.Px$coefficients[2,3], m1.Px$coefficients[2,4]))

#Add column to Z tables showing whether or not the CIs overlapped
Zs.PW$CI.overlap <- as.factor(ifelse(Zs.PW$Zt_lower>Zs.PW$Zadd_upper | Zs.PW$Zadd_lower>Zs.PW$Zt_upper,1,0 ))
Zs.DtP$CI.overlap <- as.factor(ifelse(Zs.DtP$Zt_lower>Zs.DtP$Zadd_upper | Zs.DtP$Zadd_lower>Zs.DtP$Zt_upper,1,0 ))
Zs.F$CI.overlap <- as.factor(ifelse(Zs.F$Zt_lower>Zs.F$Zadd_upper | Zs.F$Zadd_lower>Zs.F$Zt_upper,1,0 ))

#How many synergies?
Zs.PW$syn <- ifelse(Zs.PW$Zadd_lower>Zs.PW$Zt_upper,1,0)
Zs.DtP$syn <- ifelse(Zs.DtP$Zadd_lower>Zs.DtP$Zt_upper,1,0)
Zs.F$syn <- ifelse(Zs.F$Zadd_lower>Zs.F$Zt_upper,1,0)
Zs.PW$ant <- ifelse(Zs.PW$Zt_lower>Zs.PW$Zadd_upper,1,0)
Zs.DtP$ant <- ifelse(Zs.DtP$Zt_lower>Zs.DtP$Zadd_upper,1,0)
Zs.F$ant <- ifelse(Zs.F$Zt_lower>Zs.F$Zadd_upper,1,0)

sum(Zs.PW$syn)
sum(Zs.DtP$syn)
sum(Zs.F$syn)
sum(Zs.PW$ant)
sum(Zs.DtP$ant)
sum(Zs.F$ant)
#these numbers are mentioned in the discussion

#Save workspace
#note for this save I just ran the data cleaning and standardization section and the synergy section
save.image("./Outputs/Workspaces/FinalAnalyses_Pred1d_synergy")

#export data for table 2
write.csv(m.sum.RvsS, "./Outputs/Tables/TableS6_Pred1d_synergy.csv")

#---Prediction 1d null models-----
#----

#During the review process, there was some concern that these models
#rely on using the output from one analysis to inform the next, and that
#it is not clear how error and variance propagate through the models. 
#It was suggested that we create a set of null models that use a 
#randomized dataset to compare the output to those models 

#This section repeats the analysis above for Prediction 1d 1000 times with 
#different randomizations of the dataset. Each run is a new randomization
#of the three response variables in the dataset used in the analysis: pupal weight
#days to pupation, and fungal growth

#Data randomization is by species,so values are randomized for all individuals of
#a species across all treatments, including controls 


#make table to store stats across all iterations
m.sum.RvsS.null <- data.frame(Metric=NA, Sp=NA, Slope=NA, Z=NA, P=NA, Run=NA)

for(k in 1:1000){

#table to store stats for models testing effects of richness on interaction index    
m.sum.RvsS <- data.frame(Metric=NA, Sp=NA, Slope=NA, Z=NA, P=NA)

#for pupal weight

d.temp <- d.rich

#randomizing data
d.temp <- d.temp %>% 
  group_by(Species) %>% 
  mutate(Pupal.weight = Pupal.weight[sample(row_number())]) 

Zs.PW <- data.frame(Species=factor(), Treatment=factor(), Zt=numeric(), Zt.SE=numeric())
for(i in 1:length(levels(d.temp$Treatment))){
  if(levels(d.temp$Treatment)[i] != "C"){
    Tx <- levels(d.temp$Treatment)[i] 
    d <- d.temp[which(d.temp$Treatment==Tx),]
    e <- as.character(unique(d$exp))
    d2 <- d.temp[which(d.temp$Treatment==Tx | d.temp$Treatment=="C"),]
    d2 <- filter(d2, exp %in% e)
    
    for(j in 1:length(levels(d.temp$Species))){
      Sp <- levels(d.temp$Species)[j]
      d3 <- d2[which(d2$Species==Sp),]
      d3$Treatment <- factor(d3$Treatment, levels=c("C", Tx))
      m.PW <- lm(Pupal.weight ~ Treatment, data=d3)
      b <- coef(summary(m.PW))[2,1]
      b.SE <- coef(summary(m.PW))[2,2]
      newrow = data.frame(Species=Sp, Treatment=Tx, Zt=b, Zt.SE=b.SE)
      Zs.PW <- rbind(Zs.PW, newrow)
    }
  }
}


##Now calculate Zadd and variance for each mixture
#from Tallarida, the variance of Zadd is
#var_Zadd= f^2 * V_A + (1-f)^2 * V_B
#where f=the fraction of each compound in mix


Zs.PW$Zadd <- NA
Zs.PW$Zadd.SE <- NA

#To calculate Zadd, need to know the conc of each compound in each treatment.
#Zs.PW <- left_join(Zs.PW, distinct(d.temp[,c(2, 12:25)], Treatment, .keep_all=TRUE), by="Treatment")

Tx.Comp <- distinct(d.temp[,c(2, 12:25)], Treatment, .keep_all=TRUE)
Tx.Comp$Treatment <- as.character(Tx.Comp$Treatment)

for(i in 1:length(Zs.PW$Zadd)){
  Zs.PW$Zadd[i] <- 0
  Zs.PW$Zadd.SE[i] <- 0
  H <- Zs.PW$Species[i]
  Tx <- Zs.PW$Treatment[i]
  for(j in 2:15){
    C.ID <- colnames(Tx.Comp)[j]
    e <- Zs.PW$Zt[which(Zs.PW$Treatment==C.ID & Zs.PW$Species==H)]  #effect of indiv compound
    var <- Zs.PW$Zt.SE[which(Zs.PW$Treatment==C.ID & Zs.PW$Species==H)]^2  ##var of slope from model coef (se^2)  see: https://stackoverflow.com/questions/14960868/how-do-i-print-the-variance-of-an-lm-in-r-without-computing-from-the-standard-er
    a <- e*(Tx.Comp[which(Tx.Comp$Treatment==Tx),j]/0.5)  ##effect * relative abundance (note total compounds=0.5mg, so dividing by that)
    #r <- d.temp$Richness[i]
    b <- var*((Tx.Comp[which(Tx.Comp$Treatment==Tx),j]/0.5))  ##in Tallarida equation the proportions are squared, 
    #but that gives really weird results with multiple compounds, where the variance estimate gets much lower as you increase in richness
    ##here I just made the variance a straight additive proportional function of the composition
    Zs.PW$Zadd[i] <- Zs.PW$Zadd[i] + a
    Zs.PW$Zadd.SE[i] <- Zs.PW$Zadd.SE[i] + b  #this is a actually generating the variance
  }
    Zs.PW$Zadd <- as.numeric(Zs.PW$Zadd)  
    Zs.PW$Zadd.SE <- as.numeric(Zs.PW$Zadd.SE) 
    Zs.PW$Zadd.SE [i] <- sqrt(Zs.PW$Zadd.SE[i]) #square root to get SE
}


# Define Richness: 0,1,2,4,6,8,10
Zs.PW$Richness <- NA
for (i in 1:nrow(Zs.PW)){
  Zs.PW$Richness[i] <- str_extract(Zs.PW$Treatment[i],'[0-9]+')}
Zs.PW$Richness[which(is.na(Zs.PW$Richness))] <- 1
Zs.PW$Richness <- as.numeric(Zs.PW$Richness)


#Calculate CIs for Zt and Zadd
#From Tallarida 2002, total sample size for calculating CI is based on total
#number of datapoints in drug1, drug 2, and combination all combined
#So in our case we have N=15 for each compound, plus ~N=5-8 for each mixture
#so use t-value 1.960 (for large sample size, anything over df=120)
#So 95% CI is Zadd +/- int

Zs.PW$Zt_upper <- Zs.PW$Zt + (Zs.PW$Zt.SE*1.96)
Zs.PW$Zt_lower <- Zs.PW$Zt - (Zs.PW$Zt.SE*1.96)
Zs.PW$Zadd_upper <- Zs.PW$Zadd + (Zs.PW$Zadd.SE*1.96)
Zs.PW$Zadd_lower <- Zs.PW$Zadd - (Zs.PW$Zadd.SE*1.96)



Zs.PW$diff <- Zs.PW$Zadd-Zs.PW$Zt   



##Can test for increase in synergy with increasing richness

m1.Hz <- summary(lm(diff~ Richness, data=Zs.PW[which(Zs.PW$Species=="Hz" & Zs.PW$Richness !=1),]))  #no synergy
m1.Sf <- summary(lm(diff~ Richness, data=Zs.PW[which(Zs.PW$Species=="Sf" & Zs.PW$Richness !=1),]))  #no synergy
m1.Px <- summary(lm(diff~ Richness, data=Zs.PW[which(Zs.PW$Species=="Px" & Zs.PW$Richness !=1),]))  #marginal antagonism
m1.Cp <- summary(lm(diff~ Richness, data=Zs.PW[which(Zs.PW$Species=="Cp" & Zs.PW$Richness !=1),]))  #no synergy
m1.Hz
m1.Sf
m1.Cp
m1.Px




#Add stats to summary table
m.sum.RvsS[1,] <- c("PW", "Hz", m1.Hz$coefficients[2,1], m1.Hz$coefficients[2,3], m1.Hz$coefficients[2,4])
m.sum.RvsS <- rbind(m.sum.RvsS,
                    c("PW", "Sf", m1.Sf$coefficients[2,1], m1.Sf$coefficients[2,3], m1.Sf$coefficients[2,4]))
m.sum.RvsS <- rbind(m.sum.RvsS,
                    c("PW", "Cp", m1.Cp$coefficients[2,1], m1.Cp$coefficients[2,3], m1.Cp$coefficients[2,4]))
m.sum.RvsS <- rbind(m.sum.RvsS,
                    c("PW", "Px", m1.Px$coefficients[2,1], m1.Px$coefficients[2,3], m1.Px$coefficients[2,4]))





####Now the same for DtP

d.temp <- d.rich

#randomizing data
d.temp <- d.temp %>% 
  group_by(Species) %>% 
  mutate(Days.to.pupation = Days.to.pupation[sample(row_number())]) 


#taking the inverse of days to pupation--i.e. "development speed" so that the 
#directionality of effects is the same for all performance metrics
#d.temp$Days.to.pupation.inv <- 1/d.temp$Days.to.pupation
#hist(d.temp$Days.to.pupation.inv)

Zs.DtP <- data.frame(Species=factor(), Treatment=factor(), Zt=numeric(), Zt.SE=numeric())
for(i in 1:length(levels(d.temp$Treatment))){
  if(levels(d.temp$Treatment)[i] != "C"){
    Tx <- levels(d.temp$Treatment)[i] 
    d <- d.temp[which(d.temp$Treatment==Tx),]
    d2 <- d.temp[which(d.temp$Treatment==Tx | d.temp$Treatment=="C"),]
    
    for(j in 1:length(levels(d.temp$Species))){
      Sp <- levels(d.temp$Species)[j]
      d3 <- d2[which(d2$Species==Sp),]
      d3$Treatment <- factor(d3$Treatment, levels=c(Tx, "C"))
      #switched order of Tx and C here so that the directionality of the effect will be the same
      #across performance metrics. This will give the slope of the line from Tx to C
      m.DtP <- lm(Days.to.pupation ~ Treatment, data=d3)
      b <- coef(summary(m.DtP))[2,1]
      b.SE <- coef(summary(m.DtP))[2,2]
      newrow = data.frame(Species=Sp, Treatment=Tx, Zt=b, Zt.SE=b.SE)
      Zs.DtP <- rbind(Zs.DtP, newrow)
    }
  }
}

#**Tried doing these analyses comparing data for each treatment only to the controls from the 
#same experiment, but the models were getting "essentially perfect fit, may be unreliable" errors
#because in several cases all the controls for an experiment (7-9 caterpillars) pupated on the same day
#and so there was no variance. Thus, I compared to all the controls for that species for DtP, which
#increased the C sample size and I think provides a more clear picture 



##Now calculate Zadd and variance for each mixture
#from Tallarida, the variance of Zadd is
#var_Zadd= f^2 * V_A + (1-f)^2 * V_B
#where f=the fraction of each compound in mix


Zs.DtP$Zadd <- NA
Zs.DtP$Zadd.SE <- NA

Tx.Comp <- distinct(d.temp[,c(2, 12:25)], Treatment, .keep_all=TRUE)
Tx.Comp$Treatment <- as.character(Tx.Comp$Treatment)

for(i in 1:length(Zs.DtP$Zadd)){
  Zs.DtP$Zadd[i] <- 0
  Zs.DtP$Zadd.SE[i] <- 0
  H <- Zs.DtP$Species[i]
  Tx <- Zs.DtP$Treatment[i]
  for(j in 2:15){
    C.ID <- colnames(Tx.Comp)[j]
    e <- Zs.DtP$Zt[which(Zs.DtP$Treatment==C.ID & Zs.DtP$Species==H)]  #effect of indiv compound
    var <- Zs.DtP$Zt.SE[which(Zs.DtP$Treatment==C.ID & Zs.DtP$Species==H)]^2  ##var of slope from model coef (se^2)  see: https://stackoverflow.com/questions/14960868/how-do-i-print-the-variance-of-an-lm-in-r-without-computing-from-the-standard-er
    a <- e*(Tx.Comp[which(Tx.Comp$Treatment==Tx),j]/0.5)  ##effect * relative abundance (note total compounds=0.5mg, so dividing by that)
    #r <- d.temp$Richness[i]
    b <- var*((Tx.Comp[which(Tx.Comp$Treatment==Tx),j]/0.5))  ##in Tallarida equation the proportions are squared, 
    #but that gives really weird results with multiple compounds, where the variance estimate gets much lower as you increase in richness
    ##here I just made the variance a straight additive proportional function of the composition
    Zs.DtP$Zadd[i] <- Zs.DtP$Zadd[i] + a
    Zs.DtP$Zadd.SE[i] <- Zs.DtP$Zadd.SE[i] + b  #this is a actually generating the variance
  }
  Zs.DtP$Zadd <- as.numeric(Zs.DtP$Zadd)
  Zs.DtP$Zadd.SE <- as.numeric(Zs.DtP$Zadd.SE)
  Zs.DtP$Zadd.SE [i] <- sqrt(Zs.DtP$Zadd.SE[i]) #square root to get SE
}


# Define Richness: 0,1,2,4,6,8,10
Zs.DtP$Richness <- NA
for (i in 1:nrow(Zs.DtP)){
  Zs.DtP$Richness[i] <- str_extract(Zs.DtP$Treatment[i],'[0-9]+')}
Zs.DtP$Richness[which(is.na(Zs.DtP$Richness))] <- 1
Zs.DtP$Richness <- as.numeric(Zs.DtP$Richness)


#Calculate CIs for Zt and Zadd
#From Tallarida 2002, total sample size for calculating CI is based on total
#number of datapoints in drug1, drug 2, and combination all combined
#So in our case we have N=15 for each compound, plus ~N=5-8 for each mixture
#so use t-value 1.960 (for large sample size, anything over df=120)
#So 95% CI is Zadd +/- int

Zs.DtP$Zt_upper <- Zs.DtP$Zt + (Zs.DtP$Zt.SE*1.96)
Zs.DtP$Zt_lower <- Zs.DtP$Zt - (Zs.DtP$Zt.SE*1.96)
Zs.DtP$Zadd_upper <- Zs.DtP$Zadd + (Zs.DtP$Zadd.SE*1.96)
Zs.DtP$Zadd_lower <- Zs.DtP$Zadd - (Zs.DtP$Zadd.SE*1.96)


Zs.DtP$diff <- Zs.DtP$Zadd-Zs.DtP$Zt   


##Can test for increase in synergy with increasing richness

m1.Hz <- summary(lm(diff~ Richness, data=Zs.DtP[which(Zs.DtP$Species=="Hz" & Zs.DtP$Richness !=1),]))  #no synergy
m1.Sf <- summary(lm(diff~ Richness, data=Zs.DtP[which(Zs.DtP$Species=="Sf"& Zs.DtP$Richness !=1),]))  #no synergy
m1.Px <- summary(lm(diff~ Richness, data=Zs.DtP[which(Zs.DtP$Species=="Px"& Zs.DtP$Richness !=1),]))  #no synergy
m1.Cp <- summary(lm(diff~ Richness, data=Zs.DtP[which(Zs.DtP$Species=="Cp"& Zs.DtP$Richness !=1),]))  #no synergy
m1.Hz
m1.Sf
m1.Px
m1.Cp


#Add stats to summary table
m.sum.RvsS <- rbind(m.sum.RvsS, 
                    c("DtP", "Hz", m1.Hz$coefficients[2,1], m1.Hz$coefficients[2,3], m1.Hz$coefficients[2,4]))
m.sum.RvsS <- rbind(m.sum.RvsS,
                    c("DtP", "Sf", m1.Sf$coefficients[2,1], m1.Sf$coefficients[2,3], m1.Sf$coefficients[2,4]))
m.sum.RvsS <- rbind(m.sum.RvsS,
                    c("DtP", "Cp", m1.Cp$coefficients[2,1], m1.Cp$coefficients[2,3], m1.Cp$coefficients[2,4]))
m.sum.RvsS <- rbind(m.sum.RvsS,
                    c("DtP", "Px", m1.Px$coefficients[2,1], m1.Px$coefficients[2,3], m1.Px$coefficients[2,4]))


#same for fungi

d.temp <- d.fungi2
d.temp$Treatment <- as.factor(d.temp$Treatment)

#randomizing data
d.temp <- d.temp %>% 
  group_by(Fungi) %>% 
  mutate(slope = slope[sample(row_number())]) 


Zs.F <- data.frame(Species=factor(), Treatment=factor(), Zt=numeric(), Zt.SE=numeric())
for(i in 1:length(levels(d.temp$Treatment))){
  if(levels(d.temp$Treatment)[i] != "DMSO"){
    Tx <- levels(d.temp$Treatment)[i] 
    d <- d.temp[which(d.temp$Treatment==Tx),]
    #e <- as.character(unique(d$exp)) #**see note below
    d2 <- d.temp[which(d.temp$Treatment==Tx | d.temp$Treatment=="DMSO"),]
    #d2 <- filter(d2, exp %in% e)
    
    for(j in 1:length(levels(d.temp$Fungi))){
      Sp <- levels(d.temp$Fungi)[j]
      d3 <- d2[which(d2$Fungi==Sp),]
      d3$Treatment <- factor(d3$Treatment, levels=c("DMSO", Tx))
      m.F <- lm(slope ~ Treatment, data=d3)
      b <- coef(summary(m.F))[2,1]
      b.SE <- coef(summary(m.F))[2,2]
      newrow = data.frame(Species=Sp, Treatment=Tx, Zt=b, Zt.SE=b.SE)
      Zs.F <- rbind(Zs.F, newrow)
    }
  }
}


##Now calculate Zadd and variance for each mixture
#from Tallarida, the variance of Zadd is
#var_Zadd= f^2 * V_A + (1-f)^2 * V_B
#where f=the fraction of each compound in mix


Zs.F$Zadd <- NA
Zs.F$Zadd.SE <- NA

for(i in 1:length(Zs.F$Zadd)){
  Zs.F$Zadd[i] <- 0
  Zs.F$Zadd.SE[i] <- 0
  H <- Zs.F$Species[i]
  Tx <- Zs.F$Treatment[i]
  for(j in 2:15){
    C.ID <- colnames(Tx.Comp)[j]
    e <- Zs.F$Zt[which(Zs.F$Treatment==C.ID & Zs.F$Species==H)]  #effect of indiv compound
    var <- Zs.F$Zt.SE[which(Zs.F$Treatment==C.ID & Zs.F$Species==H)]^2  ##var of slope from model coef (se^2)  see: https://stackoverflow.com/questions/14960868/how-do-i-print-the-variance-of-an-lm-in-r-without-computing-from-the-standard-er
    a <- e*(Tx.Comp[which(Tx.Comp$Treatment==Tx),j]/0.5)  ##effect * relative abundance (note total compounds=0.5mg, so dividing by that)
    #r <- d.temp$Richness[i]
    b <- var*((Tx.Comp[which(Tx.Comp$Treatment==Tx),j]/0.5))  ##in Tallarida equation the proportions are squared, 
    #but that gives really weird results with multiple compounds, where the variance estimate gets much lower as you increase in richness
    ##here I just made the variance a straight additive proportional function of the composition
    Zs.F$Zadd[i] <- Zs.F$Zadd[i] + a
    Zs.F$Zadd.SE[i] <- Zs.F$Zadd.SE[i] + b  #this is actually generating the variance
  }
  Zs.F$Zadd <- as.numeric(Zs.F$Zadd)
  Zs.F$Zadd.SE <- as.numeric(Zs.F$Zadd.SE)
  Zs.F$Zadd.SE [i] <- sqrt(Zs.F$Zadd.SE[i]) #square root to get SE
}


# Define Richness: 0,1,2,4,6,8,10
Zs.F$Richness <- NA
for (i in 1:nrow(Zs.F)){
  Zs.F$Richness[i] <- str_extract(Zs.F$Treatment[i],'[0-9]+')}
Zs.F$Richness[which(is.na(Zs.F$Richness))] <- 1
Zs.F$Richness <- as.numeric(Zs.F$Richness)


#Calculate CIs for Zt and Zadd
#From Tallarida 2002, total sample size for calculating CI is based on total
#number of datapoints in drug1, drug 2, and combination all combined
#So in our case we have N=15 for each compound, plus ~N=5-8 for each mixture
#so use t-value 1.960 (for large sample size, anything over df=120)
#So 95% CI is Zadd +/- int

Zs.F$Zt_upper <- Zs.F$Zt + (Zs.F$Zt.SE*1.96)
Zs.F$Zt_lower <- Zs.F$Zt - (Zs.F$Zt.SE*1.96)
Zs.F$Zadd_upper <- Zs.F$Zadd + (Zs.F$Zadd.SE*1.96)
Zs.F$Zadd_lower <- Zs.F$Zadd - (Zs.F$Zadd.SE*1.96)


Zs.F$diff <- Zs.F$Zadd-Zs.F$Zt   


##Can test for increase in synergy with increasing richness

m1.Bo <- summary(lm(diff~ Richness, data=Zs.F[which(Zs.F$Species=="Botrys" & Zs.F$Richness !=1),]))  #no synergy
m1.Co <- summary(lm(diff~ Richness, data=Zs.F[which(Zs.F$Species=="Collet"& Zs.F$Richness !=1),]))  #no synergy
m1.Pe <- summary(lm(diff~ Richness, data=Zs.F[which(Zs.F$Species=="Penicillium"& Zs.F$Richness !=1),]))  #no synergy
m1.Sc <- summary(lm(diff~ Richness, data=Zs.F[which(Zs.F$Species=="Sclerotinia"& Zs.F$Richness !=1),]))  #no synergy
m1.Bo
m1.Co
m1.Pe
m1.Sc

#Add stats to summary table
m.sum.RvsS <- rbind(m.sum.RvsS, 
                    c("F", "Botrys", m1.Hz$coefficients[2,1], m1.Hz$coefficients[2,3], m1.Hz$coefficients[2,4]))
m.sum.RvsS <- rbind(m.sum.RvsS,
                    c("F", "Collet", m1.Sf$coefficients[2,1], m1.Sf$coefficients[2,3], m1.Sf$coefficients[2,4]))
m.sum.RvsS <- rbind(m.sum.RvsS,
                    c("F", "Penicillium", m1.Cp$coefficients[2,1], m1.Cp$coefficients[2,3], m1.Cp$coefficients[2,4]))
m.sum.RvsS <- rbind(m.sum.RvsS,
                    c("F", "Sclerotinia", m1.Px$coefficients[2,1], m1.Px$coefficients[2,3], m1.Px$coefficients[2,4]))


m.sum.RvsS$Run <- k
m.sum.RvsS.null <- rbind(m.sum.RvsS.null, m.sum.RvsS)

}

m.sum.RvsS.null$P <- as.numeric(m.sum.RvsS.null$P)

hist(m.sum.RvsS.null$P)

length(m.sum.RvsS.null$P[which(m.sum.RvsS.null$P < 0.05)])

##the p-value was < 0.05 in 474 of 12000 tests (=0.0395 proportionally) 
#this suggests the multi-step analysis approach is not biasing the 
#results toward rejecting/supporting hypotheses


#check to see if these are equally distributed across species, we would expect
#that we should have ~100 for each insect and ~50 for each fungi, 
#assuming we had 600 total false positives out of 1200 (alpha=0.05)
length(m.sum.RvsS.null$P[which(m.sum.RvsS.null$P < 0.05 & m.sum.RvsS.null$Sp=="Hz")])
length(m.sum.RvsS.null$P[which(m.sum.RvsS.null$P < 0.05 & m.sum.RvsS.null$Sp=="Sf")])
length(m.sum.RvsS.null$P[which(m.sum.RvsS.null$P < 0.05 & m.sum.RvsS.null$Sp=="Px")])
length(m.sum.RvsS.null$P[which(m.sum.RvsS.null$P < 0.05 & m.sum.RvsS.null$Sp=="Cp")])
length(m.sum.RvsS.null$P[which(m.sum.RvsS.null$P < 0.05 & m.sum.RvsS.null$Sp=="Hz")])
length(m.sum.RvsS.null$P[which(m.sum.RvsS.null$P < 0.05 & m.sum.RvsS.null$Sp=="Botrys")])
length(m.sum.RvsS.null$P[which(m.sum.RvsS.null$P < 0.05 & m.sum.RvsS.null$Sp=="Collet")])
length(m.sum.RvsS.null$P[which(m.sum.RvsS.null$P < 0.05 & m.sum.RvsS.null$Sp=="Penicillium")])
length(m.sum.RvsS.null$P[which(m.sum.RvsS.null$P < 0.05 & m.sum.RvsS.null$Sp=="Sclerotinia")])
#these are variable, but I ran this several times and there is no consistent
#pattern across species


#and there should be ~200 per response variable
length(m.sum.RvsS.null$P[which(m.sum.RvsS.null$P < 0.05 & m.sum.RvsS.null$Metric=="PW")])
length(m.sum.RvsS.null$P[which(m.sum.RvsS.null$P < 0.05 & m.sum.RvsS.null$Metric=="DtP")])
length(m.sum.RvsS.null$P[which(m.sum.RvsS.null$P < 0.05 & m.sum.RvsS.null$Metric=="F")])

save.image("./Outputs/Workspaces/FinalAnalyses_Pred1d_nullmodels")


#********************************************************************





# Comparing mixtures to average effectiveness of singletons (Prediction 1e)----
#--------------------

load("./Outputs/Workspaces/StandardizedData")

##Using data from the richness experiment. 

#make table to store stats for singleton vs mixtures
m.sum.SvsM <- data.frame(Metric=NA, Sp=NA, Slope=NA, Z=NA, P=NA)

#----------------For Pupal weights

#first looking at individual compounds to see which is most effective

d.temp <- d.rich[which(d.rich$Richness==1),]
m1 <- lm(Pupal.weight.ST ~ Treatment, data=d.temp)
summary(m1)

plot(Pupal.weight.ST ~ droplevels(Treatment), data=d.temp[which(d.temp$Species=="Hz"),])
plot(Pupal.weight.ST ~ droplevels(Treatment), data=d.temp[which(d.temp$Species=="Sf"),])
plot(Pupal.weight.ST ~ droplevels(Treatment), data=d.temp[which(d.temp$Species=="Cp"),])
plot(Pupal.weight.ST ~ droplevels(Treatment), data=d.temp[which(d.temp$Species=="Px"),])

means <- d.temp %>%
  group_by(Species, Treatment) %>%
  summarise(avg=mean(Pupal.weight.ST, na.rm=TRUE))

#For Cp and Px, most effective compound overall is hyperin
#For Sf, most effective overall is pCA
#For Hz, most effective overall is GeA

#Now getting means for each treatment

d.temp <- d.rich
means.PW <- d.temp %>%
  group_by(Species, Treatment, Richness) %>%
  summarise(avg=mean(Pupal.weight.ST, na.rm=TRUE))

#Need to create a 0/1 variable that indicates whether the mixture is more effective than the most
#effective singleton in that mixture.
means.PW$exceeds <- NA

for(i in 1:length(means.PW$Species)){
  sp <- means.PW$Species[i]
  if(means.PW$Richness[i] > 1){
    t <- means.PW$Treatment[i]
    mix <- colnames(d.rich[which(d.rich[which(d.rich$Treatment==t)[1], 12:25]>0) + 11]) 
    singleton.mean <- mean(means.PW$avg[which(means.PW$Species==sp & means.PW$Treatment %in% mix)])
    means.PW$exceeds[i] <- ifelse(means.PW$avg[i] < singleton.mean, 1, 0)
  } else {
    means.PW$exceeds[i] <- NA
  }
}


##Can now do a binomial glm to ask whether the likelihood of a treatment being more effective
#than the most effective single compound increases with increasing richness

#  For Hz
d.temp <- means.PW[which(means.PW$Species=="Hz" & means.PW$Richness >1),]
m1.Hz <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
m1.Hz

ggplot(d.temp, aes(Richness, exceeds)) +
  stat_smooth(method="loess", formula=y~x,
              alpha=0.2, size=2) +
  geom_point(position=position_jitter(height=0.03, width=0.2)) +
  xlab("Richness") + ylab("Probability More Effective \n Than Best Singleton")

#For Sf
d.temp <- means.PW[which(means.PW$Species=="Sf" & means.PW$Richness>1),]
m1.Sf <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
m1.Sf

ggplot(d.temp, aes(Richness, exceeds)) +
  stat_smooth(method="loess", formula=y~x,
              alpha=0.2, size=2) +
  geom_point(position=position_jitter(height=0.03, width=0.2)) +
  xlab("Richness") + ylab("Probability More Effective \n Than Best Singleton")


#For Px
d.temp <- means.PW[which(means.PW$Species=="Px" & means.PW$Richness >1),]
m1.Px <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
m1.Px

ggplot(d.temp, aes(Richness, exceeds)) +
  stat_smooth(method="loess", formula=y~x,
              alpha=0.2, size=2) +
  geom_point(position=position_jitter(height=0.03, width=0.2)) +
  xlab("Richness") + ylab("Probability More Effective \n Than Best Singleton")


#for Cp
d.temp <- means.PW[which(means.PW$Species=="Cp" & means.PW$Richness >1),]
m1.Cp <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
m1.Cp

ggplot(d.temp, aes(Richness, exceeds)) +
  stat_smooth(method="loess", formula=y~x,
              alpha=0.2, size=2) +
  geom_point(position=position_jitter(height=0.03, width=0.2)) +
  xlab("Richness") + ylab("Probability More Effective \n Than Best Singleton")


#Add stats to summary table
m.sum.SvsM[1,] <- c("PW", "Hz", m1.Hz$coefficients[2,1], m1.Hz$coefficients[2,3], m1.Hz$coefficients[2,4])
m.sum.SvsM <- rbind(m.sum.SvsM,
                    c("PW", "Sf", m1.Sf$coefficients[2,1], m1.Sf$coefficients[2,3], m1.Sf$coefficients[2,4]))
m.sum.SvsM <- rbind(m.sum.SvsM,
                    c("PW", "Cp", m1.Cp$coefficients[2,1], m1.Cp$coefficients[2,3], m1.Cp$coefficients[2,4]))
m.sum.SvsM <- rbind(m.sum.SvsM,
                    c("PW", "Px", m1.Px$coefficients[2,1], m1.Px$coefficients[2,3], m1.Px$coefficients[2,4]))



#----------------For Days to Pupation

#first looking at individual compounds to see which is most effective

d.temp <- d.rich[which(d.rich$Richness==1),]
m1 <- lm(Days.to.pupation.ST.inv ~ Treatment, data=d.temp)
summary(m1)

plot(Days.to.pupation.ST.inv ~ droplevels(Treatment), data=d.temp[which(d.temp$Species=="Hz"),])
plot(Days.to.pupation.ST.inv ~ droplevels(Treatment), data=d.temp[which(d.temp$Species=="Sf"),])
plot(Days.to.pupation.ST.inv ~ droplevels(Treatment), data=d.temp[which(d.temp$Species=="Cp"),])
plot(Days.to.pupation.ST.inv ~ droplevels(Treatment), data=d.temp[which(d.temp$Species=="Px"),])

means <- d.temp %>%
  group_by(Species, Treatment) %>%
  summarise(avg=mean(Days.to.pupation.ST.inv, na.rm=TRUE))

#For Cp, most effective compound is SA, but effect is miniscule
#For Sf, most effective is Q
#For Hz, most effective is ChA
#For Px, most effective is H

#Now getting means for each treatment

d.temp <- d.rich
means.DtP <- d.temp %>%
  group_by(Species, Treatment, Richness) %>%
  summarise(avg=mean(Days.to.pupation.ST.inv, na.rm=TRUE))

#Need to create a 0/1 variable that indicates whether the mixture is more effective than the most
#effective singleton in that mixture.
means.DtP$exceeds <- NA

for(i in 1:length(means.DtP$Species)){
  sp <- means.DtP$Species[i]
  if(means.DtP$Richness[i] > 1){
    t <- means.DtP$Treatment[i]
    mix <- colnames(d.rich[which(d.rich[which(d.rich$Treatment==t)[1], 12:25]>0) + 11]) 
    singleton.mean <- mean(means.DtP$avg[which(means.DtP$Species==sp & means.DtP$Treatment %in% mix)])
    means.DtP$exceeds[i] <- ifelse(means.DtP$avg[i] < singleton.mean, 1, 0)
  } else {
    means.DtP$exceeds[i] <- NA
  }
}


##Can now do a binomial glm to ask whether the likelihood of a treatment being more effective
#than the most effective single compound increases with increasing richness

#  For Hz
d.temp <- means.DtP[which(means.DtP$Species=="Hz" & means.DtP$Richness >1),]
m1.Hz <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
m1.Hz

ggplot(d.temp, aes(Richness, exceeds)) +
  stat_smooth(method="loess", formula=y~x,
              alpha=0.2, size=2) +
  geom_point(position=position_jitter(height=0.03, width=0.2)) +
  xlab("Richness") + ylab("Probability More Effective \n Than Best Singleton")


#For Sf
d.temp <- means.DtP[which(means.DtP$Species=="Sf" & means.DtP$Richness>1),]
m1.Sf <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
m1.Sf

ggplot(d.temp, aes(Richness, exceeds)) +
  stat_smooth(method="loess", formula=y~x,
              alpha=0.2, size=2) +
  geom_point(position=position_jitter(height=0.03, width=0.2)) +
  xlab("Richness") + ylab("Probability More Effective \n Than Best Singleton")



#For Px
d.temp <- means.DtP[which(means.DtP$Species=="Px" & means.DtP$Richness >1),]
m1.Px <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
m1.Px

ggplot(d.temp, aes(Richness, exceeds)) +
  stat_smooth(method="loess", formula=y~x,
              alpha=0.2, size=2) +
  geom_point(position=position_jitter(height=0.03, width=0.2)) +
  xlab("Richness") + ylab("Probability More Effective \n Than Best Singleton")


#for Cp
d.temp <- means.DtP[which(means.DtP$Species=="Cp" & means.DtP$Richness >1),]
m1.Cp <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
m1.Cp

ggplot(d.temp, aes(Richness, exceeds)) +
  stat_smooth(method="loess", formula=y~x,
              alpha=0.2, size=2) +
  geom_point(position=position_jitter(height=0.03, width=0.2)) +
  xlab("Richness") + ylab("Probability More Effective \n Than Best Singleton")


#Add stats to summary table
m.sum.SvsM <- rbind(m.sum.SvsM,
                    c("DtP", "Hz", m1.Hz$coefficients[2,1], m1.Hz$coefficients[2,3], m1.Hz$coefficients[2,4]))
m.sum.SvsM <- rbind(m.sum.SvsM,
                    c("DtP", "Sf", m1.Sf$coefficients[2,1], m1.Sf$coefficients[2,3], m1.Sf$coefficients[2,4]))
m.sum.SvsM <- rbind(m.sum.SvsM,
                    c("DtP", "Cp", m1.Cp$coefficients[2,1], m1.Cp$coefficients[2,3], m1.Cp$coefficients[2,4]))
m.sum.SvsM <- rbind(m.sum.SvsM,
                    c("DtP", "Px", m1.Px$coefficients[2,1], m1.Px$coefficients[2,3], m1.Px$coefficients[2,4]))






#----------------For Survival

#first looking at individual compounds to see which is most effective

d.temp <- d.rich.PS[which(d.rich.PS$Richness==1),]
m1 <- lm(PropSurv.ST ~ Treatment, data=d.temp)
summary(m1)

plot(PropSurv.ST ~ droplevels(Treatment), data=d.temp[which(d.temp$Species=="Hz"),])
plot(PropSurv.ST ~ droplevels(Treatment), data=d.temp[which(d.temp$Species=="Sf"),])
plot(PropSurv.ST ~ droplevels(Treatment), data=d.temp[which(d.temp$Species=="Cp"),])
plot(PropSurv.ST ~ droplevels(Treatment), data=d.temp[which(d.temp$Species=="Px"),])

#For Cp, most effective compound is GeA
#For Sf, most effective is Q
#For Hz, most effective is ChA or R (same)
#For Px, most effective is GeA

#Now getting means for each treatment

S <- d.rich.PS

#Need to create a 0/1 variable that indicates whether the mixture is more effective than the most
#effective singleton in that mixture.
S$exceeds <- NA

for(i in 1:length(S$Species)){
  sp <- S$Species[i]
  if(S$Richness[i] > 1){
    t <- S$Treatment[i]
    mix <- colnames(d.rich[which(d.rich[which(d.rich$Treatment==t)[1], 12:25]>0) + 11]) 
    singleton.mean <- mean(S$PropSurv.ST[which(S$Species==sp & S$Treatment %in% mix)])
    S$exceeds[i] <- ifelse(S$PropSurv.ST[i] < singleton.mean, 1, 0)
  } else {
    S$exceeds[i] <- NA
  }
}



##Can now do a binomial glm to ask whether the likelihood of a treatment being more effective
#than the most effective single compound increases with increasing richness

#  For Hz
d.temp <- S[which(S$Species=="Hz" & S$Richness >1),]
m1.Hz <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
m1.Hz

ggplot(d.temp, aes(Richness, exceeds)) +
  stat_smooth(method="loess", formula=y~x,
              alpha=0.2, size=2) +
  geom_point(position=position_jitter(height=0.03, width=0.2)) +
  xlab("Richness") + ylab("Probability More Effective \n Than Best Singleton")


#For Sf
d.temp <- S[which(S$Species=="Sf" & S$Richness>1),]
m1.Sf <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
m1.Sf

ggplot(d.temp, aes(Richness, exceeds)) +
  stat_smooth(method="loess", formula=y~x,
              alpha=0.2, size=2) +
  geom_point(position=position_jitter(height=0.03, width=0.2)) +
  xlab("Richness") + ylab("Probability More Effective \n Than Best Singleton")




#For Px
d.temp <- S[which(S$Species=="Px" & S$Richness >1),]
m1.Px <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
m1.Px

ggplot(d.temp, aes(Richness, exceeds)) +
  stat_smooth(method="loess", formula=y~x,
              alpha=0.2, size=2) +
  geom_point(position=position_jitter(height=0.03, width=0.2)) +
  xlab("Richness") + ylab("Probability More Effective \n Than Best Singleton")


#for Cp
d.temp <- S[which(S$Species=="Cp" & S$Richness >1),]
m1.Cp <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
m1.Cp

ggplot(d.temp, aes(Richness, exceeds)) +
  stat_smooth(method="loess", formula=y~x,
              alpha=0.2, size=2) +
  geom_point(position=position_jitter(height=0.03, width=0.2)) +
  xlab("Richness") + ylab("Probability More Effective \n Than Best Singleton")


#Add stats to summary table
m.sum.SvsM <- rbind(m.sum.SvsM,
                    c("S", "Hz", m1.Hz$coefficients[2,1], m1.Hz$coefficients[2,3], m1.Hz$coefficients[2,4]))
m.sum.SvsM <- rbind(m.sum.SvsM,
                    c("S", "Sf", m1.Sf$coefficients[2,1], m1.Sf$coefficients[2,3], m1.Sf$coefficients[2,4]))
m.sum.SvsM <- rbind(m.sum.SvsM,
                    c("S", "Cp", m1.Cp$coefficients[2,1], m1.Cp$coefficients[2,3], m1.Cp$coefficients[2,4]))
m.sum.SvsM <- rbind(m.sum.SvsM,
                    c("S", "Px", m1.Px$coefficients[2,1], m1.Px$coefficients[2,3], m1.Px$coefficients[2,4]))




#----------------For Fungal Growth

#first looking at individual compounds to see which is most effective

d.temp <- d.fungi2[which(d.fungi2$Richness==1),]
d.temp$Treatment <- as.factor(d.temp$Treatment)
m1 <- lm(dAbs.ST ~ Treatment, data=d.temp)
summary(m1)

plot(dAbs.ST ~ droplevels(Treatment), data=d.temp[which(d.temp$Fungi=="Botrys"),])
plot(dAbs.ST ~ droplevels(Treatment), data=d.temp[which(d.temp$Fungi=="Collet"),])
plot(dAbs.ST ~ droplevels(Treatment), data=d.temp[which(d.temp$Fungi=="Penicillium"),])
plot(dAbs.ST ~ droplevels(Treatment), data=d.temp[which(d.temp$Fungi=="Sclerotinia"),])

means <- d.temp %>%
  group_by(Fungi, Treatment) %>%
  summarise(avg=mean(dAbs.ST, na.rm=TRUE))

#For Penicillium, most effective is Pht
#For Sclerotinia, most effective is Q
#For Collet, most effective is CA
#For Botrys, most effective is pCA

#Now getting means for each treatment

d.temp <- d.fungi2
means.F <- d.temp %>%
  group_by(Fungi, Treatment, Richness) %>%
  summarise(avg=mean(dAbs.ST, na.rm=TRUE))

#Need to create a 0/1 variable that indicates whether the mixture is more effective than the most
#effective singleton in that mixture.
means.F$exceeds <- NA
means.F$Richness[which(means.F$Treatment=="DMSO")] <- 0

for(i in 1:length(means.F$Fungi)){
  sp <- means.F$Fungi[i]
  if(means.F$Richness[i] > 1){
    t <- means.F$Treatment[i]
    mix <- colnames(d.rich[which(d.rich[which(d.rich$Treatment==t)[1], 12:25]>0) + 11]) 
    singleton.mean <- mean(means.F$avg[which(means.F$Fungi==sp & means.F$Treatment %in% mix)])
    means.F$exceeds[i] <- ifelse(means.F$avg[i] < singleton.mean, 1, 0)
  } else {
    means.F$exceeds[i] <- NA
  }
}



##Can now do a binomial glm to ask whether the likelihood of a treatment being more effective
#than the most effective single compound increases with increasing richness

#  For Botrys
d.temp <- means.F[which(means.F$Fungi=="Botrys" & means.F$Richness >1),]
m1.Botrys <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
m1.Botrys

ggplot(d.temp, aes(Richness, exceeds)) +
  stat_smooth(method="loess", formula=y~x,
              alpha=0.2, size=2) +
  geom_point(position=position_jitter(height=0.03, width=0.2)) +
  xlab("Richness") + ylab("Probability More Effective \n Than Best Singleton")

#For Collet
d.temp <- means.F[which(means.F$Fungi=="Collet" & means.F$Richness>1),]
m1.Collet <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
m1.Collet

ggplot(d.temp, aes(Richness, exceeds)) +
  stat_smooth(method="loess", formula=y~x,
              alpha=0.2, size=2) +
  geom_point(position=position_jitter(height=0.03, width=0.2)) +
  xlab("Richness") + ylab("Probability More Effective \n Than Best Singleton")


#For Sclerotinia
d.temp <- means.F[which(means.F$Fungi=="Sclerotinia" & means.F$Richness >1),]
m1.Sclerotinia <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
m1.Sclerotinia

ggplot(d.temp, aes(Richness, exceeds)) +
  stat_smooth(method="loess", formula=y~x,
              alpha=0.2, size=2) +
  geom_point(position=position_jitter(height=0.03, width=0.2)) +
  xlab("Richness") + ylab("Probability More Effective \n Than Best Singleton")


#for Penicillium
d.temp <- means.F[which(means.F$Fungi=="Penicillium" & means.F$Richness >1),]
m1.Penicillium <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
m1.Penicillium

ggplot(d.temp, aes(Richness, exceeds)) +
  stat_smooth(method="loess", formula=y~x,
              alpha=0.2, size=2) +
  geom_point(position=position_jitter(height=0.03, width=0.2)) +
  xlab("Richness") + ylab("Probability More Effective \n Than Best Singleton")

#significant negative effect of richness on probability

#Add stats to summary table
m.sum.SvsM <- rbind(m.sum.SvsM,
                    c("F", "Botrys", m1.Botrys$coefficients[2,1], m1.Botrys$coefficients[2,3], m1.Botrys$coefficients[2,4]))
m.sum.SvsM <- rbind(m.sum.SvsM,
                    c("F", "Collet", m1.Collet$coefficients[2,1], m1.Collet$coefficients[2,3], m1.Collet$coefficients[2,4]))
m.sum.SvsM <- rbind(m.sum.SvsM,
                    c("F", "Penicillium", m1.Penicillium$coefficients[2,1], m1.Penicillium$coefficients[2,3], m1.Penicillium$coefficients[2,4]))
m.sum.SvsM <- rbind(m.sum.SvsM,
                    c("F", "Sclerotinia", m1.Sclerotinia$coefficients[2,1], m1.Sclerotinia$coefficients[2,3], m1.Sclerotinia$coefficients[2,4]))


#export data for table S7
write.csv(m.sum.SvsM, "./Outputs/Tables/TableS7_Pred1e_singletonsVSmix.csv")


#Save workspace
#note for this save I just ran the data cleaning and standardization section and the singleton vs mixture section
save.image("./Outputs/Workspaces/FinalAnalyses_Pred1e_singletonsVSmix")

#**************************************************************


#----Prediction 1e null models----
#----

#During the review process, there was some concern that these models
#rely on using the output from one analysis to inform the next, and that
#it is not clear how error and variance propagate through the models. 
#It was suggested that we create a set of null models that use a 
#randomized dataset to compare the output to those models 

#This section repeats the analysis above for Prediction 1e 100 times with 
#different randomizations of the dataset. 


#make table to store stats across all iterations
m.sum.SvsM.null <- data.frame(Metric=NA, Sp=NA, Slope=NA, Z=NA, P=NA, Run=NA)

for(k in 1:1000){

  ##Using data from the richness experiment. 
  
  #make table to store stats for singleton vs mixtures
  m.sum.SvsM <- data.frame(Metric=NA, Sp=NA, Slope=NA, Z=NA, P=NA)
  
  #----------------For Pupal weights
  
  #Getting means for each treatment
  
  d.temp <- d.rich
  
  #randomize data
  d.temp <- d.temp %>% 
    group_by(Species) %>% 
    mutate(Pupal.weight.ST = Pupal.weight.ST[sample(row_number())]) 
  
  means.PW <- d.temp %>%
    group_by(Species, Treatment, Richness) %>%
    summarise(avg=mean(Pupal.weight.ST, na.rm=TRUE))
  
  #Need to create a 0/1 variable that indicates whether the mixture is more effective than the most
  #effective singleton in that mixture.
  means.PW$exceeds <- NA
  
  for(i in 1:length(means.PW$Species)){
    sp <- means.PW$Species[i]
    if(means.PW$Richness[i] > 1){
      t <- means.PW$Treatment[i]
      mix <- colnames(d.rich[which(d.rich[which(d.rich$Treatment==t)[1], 12:25]>0) + 11]) 
      singleton.mean <- mean(means.PW$avg[which(means.PW$Species==sp & means.PW$Treatment %in% mix)])
      means.PW$exceeds[i] <- ifelse(means.PW$avg[i] < singleton.mean, 1, 0)
    } else {
      means.PW$exceeds[i] <- NA
    }
  }
  
  
  ##Can now do a binomial glm to ask whether the likelihood of a treatment being more effective
  #than the most effective single compound increases with increasing richness
  
  #  For Hz
  d.temp <- means.PW[which(means.PW$Species=="Hz" & means.PW$Richness >1),]
  m1.Hz <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))

  #For Sf
  d.temp <- means.PW[which(means.PW$Species=="Sf" & means.PW$Richness>1),]
  m1.Sf <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
  
  #For Px
  d.temp <- means.PW[which(means.PW$Species=="Px" & means.PW$Richness >1),]
  m1.Px <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
 
  #for Cp
  d.temp <- means.PW[which(means.PW$Species=="Cp" & means.PW$Richness >1),]
  m1.Cp <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))

  
  #Add stats to summary table
  m.sum.SvsM[1,] <- c("PW", "Hz", m1.Hz$coefficients[2,1], m1.Hz$coefficients[2,3], m1.Hz$coefficients[2,4])
  m.sum.SvsM <- rbind(m.sum.SvsM,
                      c("PW", "Sf", m1.Sf$coefficients[2,1], m1.Sf$coefficients[2,3], m1.Sf$coefficients[2,4]))
  m.sum.SvsM <- rbind(m.sum.SvsM,
                      c("PW", "Cp", m1.Cp$coefficients[2,1], m1.Cp$coefficients[2,3], m1.Cp$coefficients[2,4]))
  m.sum.SvsM <- rbind(m.sum.SvsM,
                      c("PW", "Px", m1.Px$coefficients[2,1], m1.Px$coefficients[2,3], m1.Px$coefficients[2,4]))
  
  
  
  #----------------For Days to Pupation
  
  #Getting means for each treatment
  
  d.temp <- d.rich
  
  #randomize data
  d.temp <- d.temp %>% 
    group_by(Species) %>% 
    mutate(Days.to.pupation.ST.inv = Days.to.pupation.ST.inv[sample(row_number())]) 
  
  means.DtP <- d.temp %>%
    group_by(Species, Treatment, Richness) %>%
    summarise(avg=mean(Days.to.pupation.ST.inv, na.rm=TRUE))

  #Need to create a 0/1 variable that indicates whether the mixture is more effective than the most
  #effective singleton in that mixture.
  means.DtP$exceeds <- NA
  
  for(i in 1:length(means.DtP$Species)){
    sp <- means.DtP$Species[i]
    if(means.DtP$Richness[i] > 1){
      t <- means.DtP$Treatment[i]
      mix <- colnames(d.rich[which(d.rich[which(d.rich$Treatment==t)[1], 12:25]>0) + 11]) 
      singleton.mean <- mean(means.DtP$avg[which(means.DtP$Species==sp & means.DtP$Treatment %in% mix)])
      means.DtP$exceeds[i] <- ifelse(means.DtP$avg[i] < singleton.mean, 1, 0)
    } else {
      means.DtP$exceeds[i] <- NA
    }
  }
  
  
  ##Can now do a binomial glm to ask whether the likelihood of a treatment being more effective
  #than the most effective single compound increases with increasing richness
  
  #  For Hz
  d.temp <- means.DtP[which(means.DtP$Species=="Hz" & means.DtP$Richness >1),]
  m1.Hz <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
 
  #For Sf
  d.temp <- means.DtP[which(means.DtP$Species=="Sf" & means.DtP$Richness>1),]
  m1.Sf <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
  
  #For Px
  d.temp <- means.DtP[which(means.DtP$Species=="Px" & means.DtP$Richness >1),]
  m1.Px <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
  
  #for Cp
  d.temp <- means.DtP[which(means.DtP$Species=="Cp" & means.DtP$Richness >1),]
  m1.Cp <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
  
  #Add stats to summary table
  m.sum.SvsM <- rbind(m.sum.SvsM,
                      c("DtP", "Hz", m1.Hz$coefficients[2,1], m1.Hz$coefficients[2,3], m1.Hz$coefficients[2,4]))
  m.sum.SvsM <- rbind(m.sum.SvsM,
                      c("DtP", "Sf", m1.Sf$coefficients[2,1], m1.Sf$coefficients[2,3], m1.Sf$coefficients[2,4]))
  m.sum.SvsM <- rbind(m.sum.SvsM,
                      c("DtP", "Cp", m1.Cp$coefficients[2,1], m1.Cp$coefficients[2,3], m1.Cp$coefficients[2,4]))
  m.sum.SvsM <- rbind(m.sum.SvsM,
                      c("DtP", "Px", m1.Px$coefficients[2,1], m1.Px$coefficients[2,3], m1.Px$coefficients[2,4]))
  
  
  
  #----------------For Survival
  
  
  #Getting means for each treatment
  
  d.temp <- d.rich.PS
  
  #randomize data
  d.temp <- d.temp %>% 
    group_by(Species) %>% 
    mutate(PropSurv.ST = PropSurv.ST[sample(row_number())]) 
  
  
  #Need to create a 0/1 variable that indicates whether the mixture is more effective than the most
  #effective singleton in that mixture.
  d.temp$exceeds <- NA
  
  for(i in 1:length(d.temp$Species)){
    sp <- d.temp$Species[i]
    if(d.temp$Richness[i] > 1){
      t <- d.temp$Treatment[i]
      mix <- colnames(d.rich[which(d.rich[which(d.rich$Treatment==t)[1], 12:25]>0) + 11]) 
      singleton.mean <- mean(d.temp$PropSurv.ST[which(d.temp$Species==sp & d.temp$Treatment %in% mix)])
      d.temp$exceeds[i] <- ifelse(d.temp$PropSurv.ST[i] < singleton.mean, 1, 0)
    } else {
      d.temp$exceeds[i] <- NA
    }
  }
  
  
  
  ##Can now do a binomial glm to ask whether the likelihood of a treatment being more effective
  #than the most effective single compound increases with increasing richness
  
  #  For Hz
  d.temp2 <- d.temp[which(d.temp$Species=="Hz" & d.temp$Richness >1),]
  m1.Hz <- summary(glm(exceeds ~ Richness, data=d.temp2, family=binomial))
  
  #For Sf
  d.temp2 <- d.temp[which(d.temp$Species=="Sf" & d.temp$Richness>1),]
  m1.Sf <- summary(glm(exceeds ~ Richness, data=d.temp2, family=binomial))
  
  #For Px
  d.temp2 <- d.temp[which(d.temp$Species=="Px" & d.temp$Richness >1),]
  m1.Px <- summary(glm(exceeds ~ Richness, data=d.temp2, family=binomial))

  #for Cp
  d.temp2 <- d.temp[which(d.temp$Species=="Cp" & d.temp$Richness >1),]
  m1.Cp <- summary(glm(exceeds ~ Richness, data=d.temp2, family=binomial))

  
  #Add stats to summary table
  m.sum.SvsM <- rbind(m.sum.SvsM,
                      c("S", "Hz", m1.Hz$coefficients[2,1], m1.Hz$coefficients[2,3], m1.Hz$coefficients[2,4]))
  m.sum.SvsM <- rbind(m.sum.SvsM,
                      c("S", "Sf", m1.Sf$coefficients[2,1], m1.Sf$coefficients[2,3], m1.Sf$coefficients[2,4]))
  m.sum.SvsM <- rbind(m.sum.SvsM,
                      c("S", "Cp", m1.Cp$coefficients[2,1], m1.Cp$coefficients[2,3], m1.Cp$coefficients[2,4]))
  m.sum.SvsM <- rbind(m.sum.SvsM,
                      c("S", "Px", m1.Px$coefficients[2,1], m1.Px$coefficients[2,3], m1.Px$coefficients[2,4]))
  
  
  
  
  #----------------For Fungal Growth
  
  #Getting means for each treatment
  
  d.temp <- d.fungi2
  
  #randomize data
  d.temp <- d.temp %>% 
    group_by(Fungi) %>% 
    mutate(dAbs.ST = dAbs.ST[sample(row_number())]) 
  
  means.F <- d.temp %>%
    group_by(Fungi, Treatment, Richness) %>%
    summarise(avg=mean(dAbs.ST, na.rm=TRUE))
  
  #Need to create a 0/1 variable that indicates whether the mixture is more effective than the most
  #effective singleton in that mixture.
  means.F$exceeds <- NA
  means.F$Richness[which(means.F$Treatment=="DMSO")] <- 0
  
  for(i in 1:length(means.F$Fungi)){
    sp <- means.F$Fungi[i]
    if(means.F$Richness[i] > 1){
      t <- means.F$Treatment[i]
      mix <- colnames(d.rich[which(d.rich[which(d.rich$Treatment==t)[1], 12:25]>0) + 11]) 
      singleton.mean <- mean(means.F$avg[which(means.F$Fungi==sp & means.F$Treatment %in% mix)])
      #singleton.mean <- mean(means.F$avg[which(means.F$Fungi==sp & means.F$Richness==1)])
      means.F$exceeds[i] <- ifelse(means.F$avg[i] < singleton.mean, 1, 0)
    } else {
      means.F$exceeds[i] <- NA
    }
  }
  
  
  
  ##Can now do a binomial glm to ask whether the likelihood of a treatment being more effective
  #than the most effective single compound increases with increasing richness
  
  #  For Botrys
  d.temp <- means.F[which(means.F$Fungi=="Botrys" & means.F$Richness >1),]
  m1.Botrys <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
 
  #For Collet
  d.temp <- means.F[which(means.F$Fungi=="Collet" & means.F$Richness>1),]
  m1.Collet <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
  
  #For Sclerotinia
  d.temp <- means.F[which(means.F$Fungi=="Sclerotinia" & means.F$Richness >1),]
  m1.Sclerotinia <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
  
  #for Penicillium
  d.temp <- means.F[which(means.F$Fungi=="Penicillium" & means.F$Richness >1),]
  m1.Penicillium <- summary(glm(exceeds ~ Richness, data=d.temp, family=binomial))
 
  
  #Add stats to summary table
  m.sum.SvsM <- rbind(m.sum.SvsM,
                      c("F", "Botrys", m1.Botrys$coefficients[2,1], m1.Botrys$coefficients[2,3], m1.Botrys$coefficients[2,4]))
  m.sum.SvsM <- rbind(m.sum.SvsM,
                      c("F", "Collet", m1.Collet$coefficients[2,1], m1.Collet$coefficients[2,3], m1.Collet$coefficients[2,4]))
  m.sum.SvsM <- rbind(m.sum.SvsM,
                      c("F", "Penicillium", m1.Penicillium$coefficients[2,1], m1.Penicillium$coefficients[2,3], m1.Penicillium$coefficients[2,4]))
  m.sum.SvsM <- rbind(m.sum.SvsM,
                      c("F", "Sclerotinia", m1.Sclerotinia$coefficients[2,1], m1.Sclerotinia$coefficients[2,3], m1.Sclerotinia$coefficients[2,4]))

  
  
  m.sum.SvsM$Run <- k
  m.sum.SvsM.null <- rbind(m.sum.SvsM.null, m.sum.SvsM)
  
}

m.sum.SvsM.null$P <- as.numeric(m.sum.SvsM.null$P)

hist(m.sum.SvsM.null$P)

length(m.sum.SvsM.null$P[which(m.sum.SvsM.null$P < 0.05)])
length(m.sum.SvsM.null$Z[which(m.sum.SvsM.null$Z < 0 )])
length(m.sum.SvsM.null$Z[which(m.sum.SvsM.null$Z < 0 & m.sum.SvsM.null$P < 0.05)])
length(m.sum.SvsM.null$Z[which(m.sum.SvsM.null$Z > 0 & m.sum.SvsM.null$P < 0.05)])

##the p-value was < 0.05 in 711 of 16000 tests (=0.044 proportionally) 
#about half are negative and half are positive
#this suggests the multi-step analysis approach is not biasing the 
#results toward rejecting/supporting hypotheses


#check to see if these are equally distributed across species, we would expect
#that we should have ~150 for each insect and ~50 for each fungi, 
#assuming we had 800 total false positives out of 1600 (alpha=0.05)
length(m.sum.SvsM.null$P[which(m.sum.SvsM.null$P < 0.05 & m.sum.SvsM.null$Sp=="Hz")])
length(m.sum.SvsM.null$P[which(m.sum.SvsM.null$P < 0.05 & m.sum.SvsM.null$Sp=="Sf")])
length(m.sum.SvsM.null$P[which(m.sum.SvsM.null$P < 0.05 & m.sum.SvsM.null$Sp=="Px")])
length(m.sum.SvsM.null$P[which(m.sum.SvsM.null$P < 0.05 & m.sum.SvsM.null$Sp=="Cp")])
length(m.sum.SvsM.null$P[which(m.sum.SvsM.null$P < 0.05 & m.sum.SvsM.null$Sp=="Botrys")])
length(m.sum.SvsM.null$P[which(m.sum.SvsM.null$P < 0.05 & m.sum.SvsM.null$Sp=="Collet")])
length(m.sum.SvsM.null$P[which(m.sum.SvsM.null$P < 0.05 & m.sum.SvsM.null$Sp=="Penicillium")])
length(m.sum.SvsM.null$P[which(m.sum.SvsM.null$P < 0.05 & m.sum.SvsM.null$Sp=="Sclerotinia")])
#these are variable, but I ran this several times and there is no consistent
#pattern across species


#and there should be ~200 per response variable
length(m.sum.SvsM.null$P[which(m.sum.SvsM.null$P < 0.05 & m.sum.SvsM.null$Metric=="PW")])
length(m.sum.SvsM.null$P[which(m.sum.SvsM.null$P < 0.05 & m.sum.SvsM.null$Metric=="DtP")])
length(m.sum.SvsM.null$P[which(m.sum.SvsM.null$P < 0.05 & m.sum.SvsM.null$Metric=="F")])
length(m.sum.SvsM.null$P[which(m.sum.SvsM.null$P < 0.05 & m.sum.SvsM.null$Metric=="S")])


  
    
#*******************************************************************************
  
  
    



# Effects of Richness and SD on Number of Effects (Prediction 2a and 2b)----
#--------------------

load("./Outputs/Workspaces/StandardizedData")

d.temp <- d.rich[which(d.rich$Treatment != "C"),]
d.temp.PS <- d.rich.PS

d.sum <- d.temp %>% 
  group_by(Richness, SD, Treatment, Species) %>%
  summarise(PW.avg = mean(Pupal.weight.ST, na.rm=TRUE), PW.SD=sd(Pupal.weight.ST, na.rm=TRUE),
            PW.n = sum(!is.na(Pupal.weight.ST)), DtP.avg= mean(Days.to.pupation.ST.inv, na.rm=TRUE), 
            DtP.SD=sd(Days.to.pupation.ST.inv, na.rm=TRUE), DtP.n = sum(!is.na(Days.to.pupation.ST.inv)))

d.sum$PW.CI.low <- d.sum$PW.avg - (1.96*d.sum$PW.SD/sqrt(d.sum$PW.n))
d.sum$PW.CI.high <- d.sum$PW.avg + (1.96*d.sum$PW.SD/sqrt(d.sum$PW.n))

d.sum$DtP.CI.low <- d.sum$DtP.avg - (1.96*d.sum$DtP.SD/sqrt(d.sum$DtP.n))
d.sum$DtP.CI.high <- d.sum$DtP.avg + (1.96*d.sum$DtP.SD/sqrt(d.sum$DtP.n))
#Note we are using 95% CIs here as the cutoff for considering an organism
#"affected" but we get similar results with different thresholds, e.g. 
#using 90% CIs (change to 1.645)


##Add binary variable that indicates whether there is a significant effect (i.e. CI
#crosses 1) for each treatment

d.sum$TxMoreRes.PW <- NA
d.sum$TxMoreRes.DtP <- NA
for(i in 1:length(d.sum$TxMoreRes.PW)){ 
  d.sum$TxMoreRes.PW[i] <-ifelse(d.sum$PW.CI.high[i] < 1, 1, 0)
  d.sum$TxMoreRes.DtP[i] <-ifelse(d.sum$DtP.CI.high[i] < 1, 1, 0)
}   


#For survival, the data are already summarized by treatment

d.temp.PS <- d.rich.PS[which(d.rich.PS$Treatment !="C"),]
d.temp.PS.C <- d.rich.PS[which(d.rich.PS$Treatment=="C"),]

d.sum$S.Tx.alive <- NA
d.sum$S.Tx.dead <- NA
d.sum$S.C.alive <- NA
d.sum$S.C.dead <- NA
d.sum$S.pval <- NA
d.sum$TxMoreRes.S <- NA

for(i in 1:length(d.temp.PS$Treatment)){
  Sp <- d.temp.PS$Species[i]
  Tx <- d.temp.PS$Treatment[i]
  exp <- d.temp.PS$exp[i]
  Cs <- d.temp.PS.C[which(d.temp.PS.C$exp==exp & d.temp.PS.C$Species==Sp),]
  
  #Is the prob of survival different in Tx vs C?
  Chi.tbl <- rbind(c(d.temp.PS$dead[i], d.temp.PS$alive[i]),
                   c(Cs$dead, Cs$alive))
  p <- fisher.test(Chi.tbl)$p.value
  
  #Is prob of survival higher in C? 1=yes, 0=no
  TxHigh <- ifelse(d.temp.PS$dead[i]/(d.temp.PS$dead[i] + d.temp.PS$alive[i]) >
                     Cs$dead/(Cs$dead + Cs$alive), 1, 0)
  TxMoreRes <- ifelse(p < 0.05 & TxHigh==1, 1, 0)
  
  #Add values to d.sum table
  d.sum$S.Tx.alive[which(d.sum$Species==Sp & d.sum$Treatment==Tx)] <- d.temp.PS$alive[i]
  d.sum$S.Tx.dead[which(d.sum$Species==Sp & d.sum$Treatment==Tx)]  <- d.temp.PS$dead[i]
  d.sum$S.C.alive[which(d.sum$Species==Sp & d.sum$Treatment==Tx)]  <- Cs$alive
  d.sum$S.C.dead[which(d.sum$Species==Sp & d.sum$Treatment==Tx)]  <- Cs$dead
  d.sum$S.pval[which(d.sum$Species==Sp & d.sum$Treatment==Tx)]  <- p
  d.sum$TxMoreRes.S[which(d.sum$Species==Sp & d.sum$Treatment==Tx)]  <- TxMoreRes
}


d.sum$OrgAffected <- ifelse(d.sum$TxMoreRes.PW==1 | d.sum$TxMoreRes.DtP==1 |d.sum$TxMoreRes.S==1, 1, 0)

#Note that this counts an insect as affected if ANY of 
#the three performance metrics are affected. There was some concern during the
#review process that this "any will do" approach could mean that we are more
#likely to count insects as affected overall than fungi, which could affect
#the outcome of the results. Thus, we also tried the same approach, but only 
#using DtP for insects, which was the most commonly affected, by running same
#analysis with the below calculation for OrgAffected. Using that, we 
#still see an effect of richness and a marginal effect of SD on the number
#of organisms affected

#d.sum$OrgAffected <- ifelse(d.sum$TxMoreRes.DtP==1, 1, 0)



d.sum <- mutate(d.sum, NumEffects=TxMoreRes.PW + TxMoreRes.DtP + TxMoreRes.S)

#Now adding in fungi

d.temp.f <- filter(d.fungi2, Treatment != "DMSO" & Treatment !="Broth" & Treatment != "Captan")

d.sum.f <- d.temp.f %>% 
  group_by(Richness, SD, Treatment, Fungi) %>%
  summarise(dAbs.avg = mean(dAbs.ST, na.rm=TRUE), dAbs.SD=sd(dAbs.ST, na.rm=TRUE),
            dAbs.n = sum(!is.na(dAbs.ST)))

d.sum.f$dAbs.CI.low <- d.sum.f$dAbs.avg - (1.96*d.sum.f$dAbs.SD/sqrt(d.sum.f$dAbs.n))
d.sum.f$dAbs.CI.high <- d.sum.f$dAbs.avg + (1.96*d.sum.f$dAbs.SD/sqrt(d.sum.f$dAbs.n))

d.sum.f <- mutate(d.sum.f, OrgAffected=ifelse(dAbs.CI.high < 1, 1, 0))
d.sum.f <- mutate(d.sum.f, NumEffects=ifelse(dAbs.CI.high < 1, 1, 0))

colnames(d.sum.f)[4] <- "Species"

#Just checking more what this looks like:

d.sum.f.tot <- d.sum.f %>%
  group_by (Species, Richness, SD) %>%
  summarise(tot=sum(OrgAffected))



d.sum.hf <- rbind(d.sum[,c(1:4, 23:24)], d.sum.f[,c(1:4, 10:11)])



d.sum2 <- d.sum.hf %>%
  group_by(Richness, SD, Treatment) %>%
  summarise(NumOrgAffected=sum(OrgAffected, na.rm=TRUE), 
            TotEffects=sum(NumEffects, na.rm=TRUE))

#Some plots
plot(jitter(NumOrgAffected) ~ jitter(Richness), data=d.sum2)
abline(lm(NumOrgAffected ~ Richness, data=d.sum2))

plot(jitter(TotEffects) ~ jitter(Richness), data=d.sum2)
abline(lm(TotEffects ~ Richness, data=d.sum2))



#checking distribution of response variable
hist(d.sum2$NumOrgAffected)
shapiro.test(d.sum2$NumOrgAffected)
#not normal, but also does not approximate typical count distributions
#like poisson or negative binomial


#checking different models

m1 <- glm(NumOrgAffected ~ Richness, data=d.sum2)
summary(m1)
m2 <- glm(NumOrgAffected ~ Richness, data=d.sum2, family=poisson)
summary(m2)
m3 <- glm(NumOrgAffected ~ Richness, data=d.sum2, family=quasipoisson)
summary(m3)  #no AIC score
m4 <- glm.nb(NumOrgAffected ~ Richness, data=d.sum2) #negative binomial
summary(m4)

AIC(m1,m2,m3,m4)


#plot(m1)
#plot(m2)
#plot(m3)
#these all look fairly similar
#but AIC lower with gaussian than poisson

hist(resid(m1))  #these actually look pretty good
shapiro.test(resid(m1))

#Seems the best test is the gaussian, but all seem to be giving qualitatively similar results
#will use results from m1

summary(m1)

#checking means
d.sum2.means <- d.sum2 %>%
  group_by(Richness) %>%
  summarise(mean=mean(NumOrgAffected, na.rm=TRUE))
d.sum2.means  #use this to report effect size in paper. Increased by 37%  (5-3.64)/3.64



#---Effects of SD on number of effects

#Some plots
plot(jitter(NumOrgAffected) ~ SD, data=d.sum2)

m1 <- glm(NumOrgAffected ~ SD, data=d.sum2)
summary(m1)
summary(glht(m1, linfct=mcp(SD="Tukey")))
drop1(m1, test="Chisq")

hist(resid(m1))

#low has fewer effects than medium or high

#checking means
d.sum2.means2 <- d.sum2 %>%
  group_by(SD) %>%
  summarise(mean=mean(NumOrgAffected, na.rm=TRUE))
d.sum2.means2  #use this to report effect size in paper. 


#save values for fig
NumEffects.rich <- d.sum2
save.image("./Outputs/Workspaces/FinalAnalyses_Pred2a-2b_NumEffects")

#**************************************************************

#--------Predictions 2a/2b null models----
#----


#During the review process, there was some concern that these models
#rely on using the output from one analysis to inform the next, and that
#it is not clear how error and variance propagate through the models. 
#It was suggested that we create a set of null models that use a 
#randomized dataset to compare the output to those models 

#This section repeats the analysis above for Prediction 2a/2b 100 times with 
#different randomizations of the dataset. 

load("./Outputs/Workspaces/StandardizedData")


#make table to store stats across all iterations
rich.sum.null <- data.frame(run=NA, Estimate=NA, t=NA, p=NA)
SD.sum.null <- data.frame(run=NA, SDL=NA, SDM=NA, SDH=NA, Fstat=NA, p=NA)

for(k in 1:10000){

d.temp <- d.rich

#randomizing data, in this case first splitting by species, then randomizing
#such that the three performance metrics (PW, DtP, S) all stay together
#for a single individual, but are randomized across treatment

d.temp.Hz <- d.temp[which(d.temp$Species=="Hz"),]
s <- sample(row_number(d.temp.Hz$Pupal.weight.ST))
d.temp.Hz <- d.temp.Hz %>% 
  mutate(Pupal.weight.ST = Pupal.weight.ST[s],
         Days.to.pupation.ST.inv = Days.to.pupation.ST.inv[s],
         Surv = Surv[s]) 

d.temp.Sf <- d.temp[which(d.temp$Species=="Sf"),]
s <- sample(row_number(d.temp.Sf$Pupal.weight.ST))
d.temp.Sf <- d.temp.Sf %>% 
  mutate(Pupal.weight.ST = Pupal.weight.ST[s],
         Days.to.pupation.ST.inv = Days.to.pupation.ST.inv[s],
         Surv = Surv[s]) 

d.temp.Cp <- d.temp[which(d.temp$Species=="Cp"),]
s <- sample(row_number(d.temp.Cp$Pupal.weight.ST))
d.temp.Cp <- d.temp.Cp %>% 
  mutate(Pupal.weight.ST = Pupal.weight.ST[s],
         Days.to.pupation.ST.inv = Days.to.pupation.ST.inv[s],
         Surv = Surv[s]) 

d.temp.Px <- d.temp[which(d.temp$Species=="Px"),]
s <- sample(row_number(d.temp.Px$Pupal.weight.ST))
d.temp.Px <- d.temp.Px %>% 
  mutate(Pupal.weight.ST = Pupal.weight.ST[s],
         Days.to.pupation.ST.inv = Days.to.pupation.ST.inv[s],
         Surv = Surv[s]) 

#re-assembling data for four species
d.temp <- rbind(d.temp.Hz, d.temp.Sf, d.temp.Cp, d.temp.Px)

#summarizing by treatment
d.sum <- d.temp %>% 
  group_by(Richness, SD, Treatment, Species) %>%
  summarise(PW.avg = mean(Pupal.weight.ST, na.rm=TRUE), PW.SD=sd(Pupal.weight.ST, na.rm=TRUE),
            PW.n = sum(!is.na(Pupal.weight.ST)), DtP.avg= mean(Days.to.pupation.ST.inv, na.rm=TRUE), 
            DtP.SD=sd(Days.to.pupation.ST.inv, na.rm=TRUE), DtP.n = sum(!is.na(Days.to.pupation.ST.inv)))%>%
  ungroup()

d.sum$PW.CI.low <- d.sum$PW.avg - (1.96*d.sum$PW.SD/sqrt(d.sum$PW.n))
d.sum$PW.CI.high <- d.sum$PW.avg + (1.96*d.sum$PW.SD/sqrt(d.sum$PW.n))

d.sum$DtP.CI.low <- d.sum$DtP.avg - (1.96*d.sum$DtP.SD/sqrt(d.sum$DtP.n))
d.sum$DtP.CI.high <- d.sum$DtP.avg + (1.96*d.sum$DtP.SD/sqrt(d.sum$DtP.n))


##Add binary variable that indicates whether there is a significant effect (i.e. CI
#crosses 1) for each treatment

d.sum$TxMoreRes.PW <- NA
d.sum$TxMoreRes.DtP <- NA
for(i in 1:length(d.sum$TxMoreRes.PW)){ 
  d.sum$TxMoreRes.PW[i] <-ifelse(d.sum$PW.CI.high[i] < 1, 1, 0)
  d.sum$TxMoreRes.DtP[i] <-ifelse(d.sum$DtP.CI.high[i] < 1, 1, 0)
}            



#For survival
#note the survival data from d.temp was already
#randomized together with the PW and DtP data above

#summarizing by treatment
d.temp2 <- d.temp %>%  group_by(Species, Richness, SD, Treatment, exp) %>%
  summarise(dead=length(which(Surv==0)), alive=length(which(Surv==1)), 
            PropSurv=length(which(Surv==1))/length(Surv)) 

d.temp.PS <- d.temp2[which(d.temp2$Treatment !="C"),]
d.temp.PS.C <- d.temp2[which(d.temp2$Treatment=="C"),]


d.sum$S.Tx.alive <- NA
d.sum$S.Tx.dead <- NA
d.sum$S.C.alive <- NA
d.sum$S.C.dead <- NA
d.sum$S.pval <- NA
d.sum$TxMoreRes.S <- NA

for(i in 1:length(d.temp.PS$Treatment)){
  Sp <- d.temp.PS$Species[i]
  Tx <- d.temp.PS$Treatment[i]
  exp <- d.temp.PS$exp[i]
  Cs <- d.temp.PS.C[which(d.temp.PS.C$exp==exp & d.temp.PS.C$Species==Sp),]
  
  #Is the prob of survival different in Tx vs C?
  Chi.tbl <- rbind(c(d.temp.PS$dead[i], d.temp.PS$alive[i]),
                   c(Cs$dead, Cs$alive))
  p <- fisher.test(Chi.tbl)$p.value
  
  #Is prob of survival higher in C? 1=yes, 0=no
  TxHigh <- ifelse(d.temp.PS$dead[i]/(d.temp.PS$dead[i] + d.temp.PS$alive[i]) >
                     Cs$dead/(Cs$dead + Cs$alive), 1, 0)
  TxMoreRes <- ifelse(p < 0.05 & TxHigh==1, 1, 0)
 
  
  #Add values to d.sum table
  d.sum$S.Tx.alive[which(d.sum$Species==Sp & d.sum$Treatment==Tx)] <- d.temp.PS$alive[i]
  d.sum$S.Tx.dead[which(d.sum$Species==Sp & d.sum$Treatment==Tx)]  <- d.temp.PS$dead[i]
  d.sum$S.C.alive[which(d.sum$Species==Sp & d.sum$Treatment==Tx)]  <- Cs$alive
  d.sum$S.C.dead[which(d.sum$Species==Sp & d.sum$Treatment==Tx)]  <- Cs$dead
  d.sum$S.pval[which(d.sum$Species==Sp & d.sum$Treatment==Tx)]  <- p
  d.sum$TxMoreRes.S[which(d.sum$Species==Sp & d.sum$Treatment==Tx)]  <- TxMoreRes
}


d.sum$OrgAffected <- ifelse(d.sum$TxMoreRes.PW==1 | d.sum$TxMoreRes.DtP==1 | d.sum$TxMoreRes.S==1, 1, 0)
d.sum <- mutate(d.sum, NumEffects=TxMoreRes.PW + TxMoreRes.DtP + TxMoreRes.S)


#Now adding in fungi
d.temp.f <- filter(d.fungi2, Treatment != "DMSO" & Treatment !="Broth" & Treatment != "Captan")

#randomize data
d.temp.f <- d.temp.f %>% 
  group_by(Fungi) %>% 
  mutate(dAbs.ST = dAbs.ST[sample(row_number())]) %>%
  ungroup()

d.sum.f <- d.temp.f %>% 
  group_by(Richness, SD, Treatment, Fungi) %>%
  summarise(dAbs.avg = mean(dAbs.ST, na.rm=TRUE), dAbs.SD=sd(dAbs.ST, na.rm=TRUE),
            dAbs.n = sum(!is.na(dAbs.ST))) %>%
  ungroup()

d.sum.f$dAbs.CI.low <- d.sum.f$dAbs.avg - (1.96*d.sum.f$dAbs.SD/sqrt(d.sum.f$dAbs.n))
d.sum.f$dAbs.CI.high <- d.sum.f$dAbs.avg + (1.96*d.sum.f$dAbs.SD/sqrt(d.sum.f$dAbs.n))

d.sum.f <- mutate(d.sum.f, OrgAffected=ifelse(dAbs.CI.high < 1, 1, 0))
d.sum.f <- mutate(d.sum.f, NumEffects=ifelse(dAbs.CI.high < 1, 1, 0))
                  

colnames(d.sum.f)[4] <- "Species"


d.sum.hf <- rbind(d.sum[,c(1:4, 23:24)], d.sum.f[,c(1:4, 10:11)])


d.sum2 <- d.sum.hf %>%
  group_by(Richness, SD, Treatment) %>%
  summarise(NumOrgAffected=sum(OrgAffected, na.rm=TRUE), 
            TotEffects=sum(NumEffects, na.rm=TRUE)) %>%
  ungroup()

m1 <- summary(glm(NumOrgAffected ~ Richness, data=d.sum2))

newrow <- data.frame(run=k, Estimate=coef(m1)[2,1], t=coef(m1)[2,3], p=coef(m1)[2,4])
rich.sum.null <- rbind(rich.sum.null, newrow)


#---Effects of SD on number of effects

m1 <- anova(glm(NumOrgAffected ~ SD, data=d.sum2), test="F")
m1b <- summary(glm(NumOrgAffected ~ SD -1, data=d.sum2))


newrow <- data.frame(run=k, SDL=coef(m1b)[2,1], SDM=coef(m1b)[3,1],
                     SDH=coef(m1b)[1,1], Fstat=m1[2,5], p=m1[2,6])

SD.sum.null <- rbind(SD.sum.null, newrow)

}


hist(rich.sum.null$p)

length(rich.sum.null$p[which(rich.sum.null$p < 0.05)])
length(rich.sum.null$Estimate[which(rich.sum.null$Estimate < 0 )])
length(rich.sum.null$Estimate[which(rich.sum.null$Estimate < 0 & rich.sum.null$p < 0.05)])
length(rich.sum.null$Estimate[which(rich.sum.null$Estimate > 0 & rich.sum.null$p < 0.05)])


hist(SD.sum.null$p)
hist(SD.sum.null$SDL)
hist(SD.sum.null$SDM)
hist(SD.sum.null$SDH)

length(SD.sum.null$p[which(SD.sum.null$p < 0.05)])
length(SD.sum.null$p[which(SD.sum.null$p < 0.05 & SD.sum.null$SDL > SD.sum.null$SDM
                     & SD.sum.null$SDL > SD.sum.null$SDH)])
length(SD.sum.null$p[which(SD.sum.null$p < 0.05 & SD.sum.null$SDM > SD.sum.null$SDH
                           & SD.sum.null$SDM > SD.sum.null$SDL)])
length(SD.sum.null$p[which(SD.sum.null$p < 0.05 & SD.sum.null$SDH > SD.sum.null$SDM
                           & SD.sum.null$SDH > SD.sum.null$SDL)])


mean(SD.sum.null$SDL, na.rm=TRUE)
mean(SD.sum.null$SDM, na.rm=TRUE)
mean(SD.sum.null$SDH, na.rm=TRUE)

#With real data, effect size for L-H is -0.90575
ef <- SD.sum.null$SDL-SD.sum.null$SDH
min(ef, na.rm=TRUE)
length(ef[which(ef < -0.90)]) #3 cases of this magnitude

#with real data, effect size for M-L is 0.97
ef <- SD.sum.null$SDM-SD.sum.null$SDL
max(ef, na.rm=TRUE)
length(ef[which(ef > 0.97)]) #one case of this magnitude


save.image("./Outputs/Workspaces/Pred2a-2b_null-models")




#**************************************************************





# Most compounds effective against at least one consumer (Prediction 2d) OR NOT (Prediction 3a)----
#--------------------

load("./Outputs/Workspaces/StandardizedData")

d.temp <- d.rich[which(d.rich$Richness==1),]
d.temp.PS <- d.rich.PS

#We will use 95% confidence intervals for the standardized variables
#to test if there is support for an effect of each individual compound
#95% CIs that do not cross one would indicate support for that effect
#Because sample sizes are sometimes small, we should calculate 95% CIs
#based on the t-distribution rather than the Z-distribution, so can't just 
#calculate CIs with the standard 1.96* formula
#We can use output from one-sample t-tests for CIs and p

t.test.p <- function(x, mu){
  if(sd(x, na.rm=TRUE)>0){
    p <- t.test(x, mu=mu)[3]
  }else{
    p <- NA
  }
  return(as.numeric(p))
}

CI.low <- function(x, mu){
  if(sd(x, na.rm=TRUE)>0){
    ci <- t.test(x, mu=mu)$conf.int[1]
  }else{
    ci <- NA
  }
  return(as.numeric(ci))
}

CI.high <- function(x, mu){
  if(sd(x, na.rm=TRUE)>0){
    ci <- t.test(x, mu=mu)$conf.int[2]
  }else{
    ci <- NA
  }
  return(as.numeric(ci))
}

d.sum <- d.temp %>% 
  group_by(Treatment, Species) %>%
  summarise(PW.avg = mean(Pupal.weight.ST, na.rm=TRUE), 
            PW.SD=sd(Pupal.weight.ST, na.rm=TRUE),
            PW.n = sum(!is.na(Pupal.weight.ST)), 
            PW.p = t.test.p(Pupal.weight.ST, mu=1),
            PW.CI.low = CI.low(Pupal.weight.ST, mu=1),
            PW.CI.high = CI.high(Pupal.weight.ST, mu=1),
            DtP.avg= mean(Days.to.pupation.ST.inv, na.rm=TRUE), 
            DtP.SD=sd(Days.to.pupation.ST.inv, na.rm=TRUE), 
            DtP.n = sum(!is.na(Days.to.pupation.ST.inv)), 
            DtP.p = t.test.p(Days.to.pupation.ST.inv, mu=1),
            DtP.CI.low = CI.low(Days.to.pupation.ST.inv, mu=1),
            DtP.CI.high = CI.high(Days.to.pupation.ST.inv, mu=1))



##Add binary variable that indicates whether there is a significant effect (i.e. CI
#crosses 1) for each treatment

d.sum$TxMoreRes.PW <- NA
d.sum$TxMoreRes.DtP <- NA
for(i in 1:length(d.sum$TxMoreRes.PW)){ 
  d.sum$TxMoreRes.PW[i] <-ifelse(d.sum$PW.CI.high[i] < 1, 1, 0)
  d.sum$TxMoreRes.DtP[i] <-ifelse(d.sum$DtP.CI.high[i] < 1, 1, 0)
}   
d.sum$TxMoreRes.DtP[which(is.na(d.sum$TxMoreRes.DtP))] <- 0

#For survival, the data are already summarized by treatment. We will use chi-squared
#tests of count data to assess if there is a significant effect

d.temp.PS <- d.rich.PS[which(d.rich.PS$Richness==1),]
d.temp.PS.C <- d.rich.PS[which(d.rich.PS$Treatment=="C"),]

d.sum$S.Tx.alive <- NA
d.sum$S.Tx.dead <- NA
d.sum$S.C.alive <- NA
d.sum$S.C.dead <- NA
d.sum$S.pval <- NA
d.sum$TxMoreRes.S <- NA
d.sum$PropSurv.ST <- NA

for(i in 1:length(d.temp.PS$Treatment)){
  Sp <- d.temp.PS$Species[i]
  Tx <- d.temp.PS$Treatment[i]
  exp <- d.temp.PS$exp[i]
  Cs <- d.temp.PS.C[which(d.temp.PS.C$exp==exp & d.temp.PS.C$Species==Sp),]
  
  #Is the prob of survival different in Tx vs C?
  Chi.tbl <- rbind(c(d.temp.PS$dead[i], d.temp.PS$alive[i]),
                   c(Cs$dead, Cs$alive))
  p <- fisher.test(Chi.tbl)$p.value
  
  #Is prob of survival higher in C? 1=yes, 0=no
  TxHigh <- ifelse(d.temp.PS$dead[i]/(d.temp.PS$dead[i] + d.temp.PS$alive[i]) >
                     Cs$dead/(Cs$dead + Cs$alive), 1, 0)
  TxMoreRes <- ifelse(p < 0.05 & TxHigh==1, 1, 0)
  
  #Add values to d.sum table
  d.sum$S.Tx.alive[which(d.sum$Species==Sp & d.sum$Treatment==Tx)] <- d.temp.PS$alive[i]
  d.sum$S.Tx.dead[which(d.sum$Species==Sp & d.sum$Treatment==Tx)]  <- d.temp.PS$dead[i]
  d.sum$S.C.alive[which(d.sum$Species==Sp & d.sum$Treatment==Tx)]  <- Cs$alive
  d.sum$S.C.dead[which(d.sum$Species==Sp & d.sum$Treatment==Tx)]  <- Cs$dead
  d.sum$S.pval[which(d.sum$Species==Sp & d.sum$Treatment==Tx)]  <- p
  d.sum$TxMoreRes.S[which(d.sum$Species==Sp & d.sum$Treatment==Tx)]  <- TxMoreRes
  d.sum$PropSurv.ST[which(d.sum$Species==Sp & d.sum$Treatment==Tx)]  <- d.temp.PS$PropSurv.ST[i]
}


d.sum$HerbAffected <- ifelse(d.sum$TxMoreRes.PW==1 | d.sum$TxMoreRes.DtP==1 |d.sum$TxMoreRes.S==1, 1, 0)
d.sum$HerbAffected[which(is.na(d.sum$HerbAffected))] <- 0
d.sum <- mutate(d.sum, NumEffects=TxMoreRes.PW + TxMoreRes.DtP + TxMoreRes.S)


d.sum2 <- d.sum %>%
  group_by(Treatment) %>%
  summarise(NumHerbsAffected=sum(HerbAffected, na.rm=TRUE), NumEffectsTot=sum(NumEffects, na.rm=TRUE))



##Can also add to this the number of significant POSITIVE effects of phenolics on performance


d.sum$TxLessRes.PW <- NA
d.sum$TxLessRes.DtP <- NA
for(i in 1:length(d.sum$TxLessRes.PW)){ 
  d.sum$TxLessRes.PW[i] <-ifelse(d.sum$PW.CI.low[i] > 1, 1, 0)
  d.sum$TxLessRes.DtP[i] <-ifelse(d.sum$DtP.CI.low[i] > 1, 1, 0)
}   
d.sum$TxLessRes.DtP[which(is.na(d.sum$TxLessRes.DtP))] <- 0

d.sum$TxLessRes.S <- NA
for(i in 1:length(d.sum$TxLessRes.S)){ 
  d.sum$TxLessRes.S[i] <-ifelse(d.sum$S.pval[i] < 0.05 & (d.sum$S.Tx.alive[i]/(d.sum$S.Tx.alive[i]+d.sum$S.Tx.dead[i])) > 
                                  (d.sum$S.C.alive[i]/(d.sum$S.C.alive[i]+d.sum$S.C.dead[i])), 1, 0)
} 


d.sum$HerbPosAffected <- ifelse(d.sum$TxLessRes.PW==1 | d.sum$TxLessRes.DtP==1 |d.sum$TxLessRes.S==1, 1, 0)
d.sum$HerbPosAffected[which(is.na(d.sum$HerbPosAffected))] <- 0
d.sum <- mutate(d.sum, NumPosEffects=TxLessRes.PW + TxLessRes.DtP + TxLessRes.S)


d.sum2 <- d.sum %>%
  group_by(Treatment) %>%
  summarise(NumHerbsAffected=sum(HerbAffected, na.rm=TRUE), NumEffectsTot=sum(NumEffects, na.rm=TRUE),
            NumHerbsPosAffected=sum(HerbPosAffected, na.rm=TRUE), NumPosEffectsTot=sum(NumPosEffects, na.rm=TRUE))


####Now adding in fungi

d.temp <- d.fungi2[which(d.fungi2$Richness==1),]

d.sum.F <- d.temp %>% 
  group_by(Treatment, Fungi) %>%
  summarise(F.avg = mean(dAbs.ST, na.rm=TRUE), F.SD=sd(dAbs.ST, na.rm=TRUE),
            F.n = sum(!is.na(dAbs.ST)), F.p=t.test.p(dAbs.ST, mu=1),
            F.CI.low = CI.low(dAbs.ST, mu=1),
            F.CI.high = CI.high(dAbs.ST, mu=1))


##Add binary variable that indicates whether there is a significant effect (i.e. CI
#crosses 1) for each treatment

d.sum.F$TxMoreRes.F <- NA
for(i in 1:length(d.sum.F$TxMoreRes.F)){ 
  d.sum.F$TxMoreRes.F[i] <-ifelse(d.sum.F$F.CI.high[i] < 1, 1, 0)
}   

#And another variable to indicate whether there are significant POSITIVE effects
d.sum.F$TxLessRes.F <- NA
for(i in 1:length(d.sum.F$TxLessRes.F)){ 
  d.sum.F$TxLessRes.F[i] <-ifelse(d.sum.F$F.CI.low[i] > 1, 1, 0)
}  


#total numbers of effects for paper
sum(d.sum$HerbAffected)
sum(d.sum$HerbPosAffected)
length(d.sum$HerbAffected)

sum(sum(d.sum$TxMoreRes.PW) + sum(d.sum$TxMoreRes.DtP) + sum(d.sum$TxMoreRes.S))
sum(sum(d.sum$TxLessRes.PW) + sum(d.sum$TxLessRes.DtP) + sum(d.sum$TxLessRes.S))
length(d.sum$TxLessRes.PW) + length(d.sum$TxLessRes.DtP) + length(d.sum$TxLessRes.S)

sum(d.sum.F$TxMoreRes.F)
sum(d.sum.F$TxLessRes.F)
length(d.sum.F$TxMoreRes.F)


#Save these tables for figs

save.image("./Outputs/Workspaces/FinalAnalyses_Pred2d&3a_IndivComps")

#**************************************************************





# Effect of compound depends on herbivore identity (Prediction 2e)----
#--------------------

load("./Outputs/Workspaces/StandardizedData")

d.temp <- d.rich [which(d.rich$Richness !=0),] 
hist(d.temp$Pupal.weight.ST)
hist(d.temp$Days.to.pupation.ST)

#convergence issues due to the random effects structure they go away when I drop
#included exp/TrayID, but results are not substantially different either way
m.PW.1 <- lmer(Pupal.weight.ST ~ Species*Treatment + (1|exp/TrayID), data=d.temp)
m.PW.2 <- lmer(Pupal.weight.ST ~ Species+Treatment + (1|exp/TrayID), data=d.temp)
summary(m.PW.1)
plot(m.PW.1)
drop1(m.PW.1, test="Chisq")  #use for Table S8
drop1(m.PW.2, test="Chisq")

m.DtP.1 <- lmer(Days.to.pupation.ST ~ Species*Treatment + (1|exp), data=d.temp)
m.DtP.2 <- lmer(Days.to.pupation.ST ~ Species+Treatment + (1|exp), data=d.temp)
summary(m.DtP.1)
plot(m.DtP.1)
drop1(m.DtP.1, test="Chisq") #use for Table S8
drop1(m.DtP.2, test="Chisq")

d.temp <- d.rich.PS
hist(d.temp$PropSurv.ST)
m.S.1 <- lmer(PropSurv.ST ~ Species*Treatment + (1|exp), data=d.temp)
m.S.2 <- lmer(PropSurv.ST ~ Species+Treatment + (1|exp), data=d.temp)
summary(m.S.1)
drop1(m.S.1, test="Chisq") #Use for Table S8
drop1(m.S.2, test="Chisq")


d.temp <- d.fungi2[which(d.fungi2$Richness !=0),] 
hist(d.temp$dAbs.ST)
m.F.1 <- lmer(dAbs.ST ~ Fungi*Treatment + (1|Plate), data=d.temp)
m.F.2 <- lmer(dAbs.ST ~ Fungi+Treatment + (1|Plate), data=d.temp)
summary(m.F.1)
plot(m.F.1)
drop1(m.F.1, test="Chisq") #Use for Table S8
drop1(m.F.2, test="Chisq")

#**************************************************************




# One or few compounds sufficient to explain performance (Prediction 3b)----   
#--------------------

load("./Outputs/Workspaces/StandardizedData")

#for pupal weights
d.temp <- filter (d.rich, Treatment != "C", !is.na(Pupal.weight.ST))
d.temp <- mutate_at(d.temp, 12:25, list("pa" = ~ifelse(.==0, "a", "p")))
m1.PW <- lm(Pupal.weight.ST ~ ChA_pa + CA_pa + pCA_pa + FA_pa + GA_pa + SA_pa + GeA_pa + Ct_pa + eCt_pa + R_pa
            + H_pa + Q_pa + Phz_pa + Pht_pa + Species, data = d.temp, na.action= "na.fail")
summary(m1.PW)
s1.PW <- stepAIC(m1.PW, direction = "both", trace = FALSE)
summary(s1.PW)
#plot(s1.PW)
tbl.S9.PW <- data.frame(summary(s1.PW)$coefficients)


#for days to pupation
d.temp <- filter (d.rich, Treatment != "C", !is.na(Days.to.pupation.ST.inv))
d.temp <- mutate_at(d.temp, 12:25, funs("pa" = ifelse(.==0, "a", "p")))
m1.DtP <- lm(Days.to.pupation.ST.inv ~ ChA_pa + CA_pa + pCA_pa + FA_pa + GA_pa + SA_pa + GeA_pa + Ct_pa + eCt_pa + R_pa
             + H_pa + Q_pa + Phz_pa + Pht_pa + Species, data = d.temp, na.action= "na.fail")
s1.DtP <- stepAIC(m1.DtP, direction = "both", trace = FALSE)
summary(s1.DtP)
#plot(s1.DtP)
tbl.S9.DtP <- data.frame(summary(s1.DtP)$coefficients)

#for survival
d.temp <- filter (d.rich.PS, Treatment != "C", !is.na(PropSurv.ST))
d.temp <- left_join(d.temp, d.rich[,c(2, 12:25)], by="Treatment")
d.temp <- mutate_at(d.temp, 10:23, funs("pa" = ifelse(.==0, "a", "p")))
m1.S <- lm(PropSurv.ST ~ ChA_pa + CA_pa + pCA_pa + FA_pa + GA_pa + SA_pa + GeA_pa + Ct_pa + eCt_pa + R_pa
           + H_pa + Q_pa + Phz_pa + Pht_pa + Species, data = d.temp, na.action= "na.fail")
s1.S <- stepAIC(m1.S, direction = "both", trace = FALSE)
summary(s1.S)
#plot(s1.S)
tbl.S9.S <- data.frame(summary(s1.S)$coefficients)


#for fungal growth rates
d.temp <- filter (d.fungi2, Treatment != "DMSO", !is.na(dAbs.ST))
d.temp <- left_join(d.temp, d.rich[,c(2, 12:25)], by="Treatment")
d.temp <- mutate_at(d.temp, 13:26, funs("pa" = ifelse(.==0, "a", "p")))
m1.F <- lm(dAbs.ST ~ ChA_pa + CA_pa + pCA_pa + FA_pa + GA_pa + SA_pa + GeA_pa + Ct_pa + eCt_pa + R_pa +
             H_pa + Q_pa + Phz_pa + Pht_pa + Fungi, data = d.temp,na.action= "na.fail")
s1.F <- stepAIC(m1.F, direction = "both", trace = FALSE)
summary(s1.F)
#plot(s1.DtP)
tbl.S9.F <- data.frame(summary(s1.F)$coefficients)


#Use tbl.S9.-  PW, DtP, S, and F for Table S9

#**************************************************************




# Little or no variation is explained by specific mixture (Prediction 3c)----
#--------------------

load("./Outputs/Workspaces/StandardizedData")


#---- Hz; treatment as fixed effect

#for pupal weights
d.temp <- filter (d.rich, Treatment != "C", Species=="Hz", !is.na(Pupal.weight.ST))
m1 <- lmer(Pupal.weight.ST ~ Treatment + (1|exp/TrayID), data=d.temp)
summary(m1)  #note singular fit disappears with removal of TrayID from model but result is the same
drop1(m1, test="Chisq")

#for DtP
d.temp <- filter (d.rich, Treatment != "C", Species=="Hz", !is.na(Days.to.pupation.ST))
m1 <- lmer(Days.to.pupation.ST.inv ~ Treatment + (1|exp/TrayID) , data=d.temp)
summary(m1)
drop1(m1, test="Chisq")

#for S
d.temp <- filter (d.rich, Treatment != "C")
m1 <- glm(Surv ~ Treatment, data=d.temp, family=binomial)
summary(m1)  #convergence errors with glmm and exp as random effect, but this seems to work
drop1(m1, test="Chisq")

#strong effects of treatment in all cases

#---- Hz; Treatment as random effect

#for pupal weights
d.temp <- filter (d.rich, Treatment != "C", Species=="Hz", !is.na(Pupal.weight.ST))
m1.Hz.PW <- lmer(Pupal.weight.ST ~ (1|Treatment) + (1|exp), data=d.temp)
summary(m1.Hz.PW)
vc <- data.frame(VarCorr(m1.Hz.PW))
vc.Hz.PW <- mutate(vc, percVar=vcov/sum(vcov))
vc.Hz.PW
m2.Hz.PW <- lmer(Pupal.weight.ST ~ (1|exp), data=d.temp)
summary(m2.Hz.PW)
anova(m1.Hz.PW, m2.Hz.PW)  #model is better with treatment


#for DtP
d.temp <- filter (d.rich, Treatment != "C", Species=="Hz", !is.na(Days.to.pupation.ST))
m1.Hz.DtP <- lmer(Days.to.pupation.ST.inv ~ (1|Treatment) + (1|exp) , data=d.temp)
summary(m1.Hz.DtP)
vc <- data.frame(VarCorr(m1.Hz.DtP))
vc.Hz.DtP <- mutate(vc, percVar=vcov/sum(vcov))
vc.Hz.DtP   #43% of variance explained by treatment
m2.Hz.DtP <- lmer(Days.to.pupation.ST.inv ~ (1|exp), data=d.temp)
summary(m2.Hz.DtP)
anova(m1.Hz.DtP, m2.Hz.DtP)  #model is better with treatment

#for S
d.temp <- filter (d.rich, Treatment != "C", Species=="Hz")
m1.Hz.S <- glmer(Surv ~ (1|Treatment) + (1|exp), data=d.temp, family=binomial)
summary(m1.Hz.S)
vc <- data.frame(VarCorr(m1.Hz.S)) 
icc <- as.numeric(icc(m1.Hz.S)[1]) #estimates % of variance explained by all random effects combined
vc.Hz.S <- mutate(vc, percVar=vcov*icc/sum(vcov))
vc.Hz.S   #26.6% of variance explained by treatment
m2.Hz.S <- glmer(Surv ~ (1|exp), data=d.temp, family=binomial)
summary(m2.Hz.S)
anova(m1.Hz.S, m2.Hz.S)  #model is better with treatment

tbl.S10.Hz <- data.frame(Sp="Hz", metric=c("PW", "DtP", "S"), 
                       var=c(vc.Hz.PW$vcov[1], vc.Hz.DtP$vcov[1], vc.Hz.S$vcov[1]), 
                       percVar=c(vc.Hz.PW$percVar[1], vc.Hz.DtP$percVar[1], vc.Hz.S$percVar[1]), 
                       Chisq=c(anova(m1.Hz.PW, m2.Hz.PW)$Chisq[2],anova(m1.Hz.DtP, m2.Hz.DtP)$Chisq[2],
                               anova(m1.Hz.S, m2.Hz.S)$Chisq[2]), 
                       P=c(anova(m1.Hz.PW, m2.Hz.PW)$Pr[2], anova(m1.Hz.DtP, m2.Hz.DtP)$Pr[2],
                           anova(m1.Hz.S, m2.Hz.S)$Pr[2]))
               

#---- Sf; Treatment as random effect

#for pupal weights
d.temp <- filter (d.rich, Treatment != "C", Species=="Sf", !is.na(Pupal.weight.ST))
m1.Sf.PW <- lmer(Pupal.weight.ST ~ (1|Treatment) + (1|exp), data=d.temp)
summary(m1.Sf.PW)
vc <- data.frame(VarCorr(m1.Sf.PW))
vc.Sf.PW <- mutate(vc, percVar=vcov/sum(vcov))
vc.Sf.PW
m2.Sf.PW <- lmer(Pupal.weight.ST ~ (1|exp), data=d.temp)
summary(m2.Sf.PW)
anova(m1.Sf.PW, m2.Sf.PW)  #model is better with treatment


#for DtP
d.temp <- filter (d.rich, Treatment != "C", Species=="Sf", !is.na(Days.to.pupation.ST))
m1.Sf.DtP <- lmer(Days.to.pupation.ST.inv ~ (1|Treatment) + (1|exp) , data=d.temp)
summary(m1.Sf.DtP)
vc <- data.frame(VarCorr(m1.Sf.DtP))
vc.Sf.DtP <- mutate(vc, percVar=vcov/sum(vcov))
vc.Sf.DtP   #73% of variance explained by treatment
m2.Sf.DtP <- lmer(Days.to.pupation.ST.inv ~ (1|exp), data=d.temp)
summary(m2.Sf.DtP)
anova(m1.Sf.DtP, m2.Sf.DtP)  #model is better with treatment

#for S
d.temp <- filter (d.rich, Treatment != "C", Species=="Sf")
m1.Sf.S <- glmer(Surv ~ (1|Treatment) + (1|exp), data=d.temp, family=binomial)
summary(m1.Sf.S)
vc <- data.frame(VarCorr(m1.Sf.S))
icc <- as.numeric(icc(m1.Sf.S)[1]) #estimates % of variance explained by all random effects combined
vc.Sf.S <- mutate(vc, percVar=vcov*icc/sum(vcov))
vc.Sf.S   #40.89% of variance explained by treatment
m2.Sf.S <- glmer(Surv ~ (1|exp), data=d.temp, family=binomial)
summary(m2.Sf.S)
anova(m1.Sf.S, m2.Sf.S)  #model is better with treatment

tbl.S10.Sf <- data.frame(Sp="Sf", metric=c("PW", "DtP", "S"), 
                       var=c(vc.Sf.PW$vcov[1], vc.Sf.DtP$vcov[1], vc.Sf.S$vcov[1]), 
                       percVar=c(vc.Sf.PW$percVar[1], vc.Sf.DtP$percVar[1], vc.Sf.S$percVar[1]), 
                       Chisq=c(anova(m1.Sf.PW, m2.Sf.PW)$Chisq[2],anova(m1.Sf.DtP, m2.Sf.DtP)$Chisq[2],
                               anova(m1.Sf.S, m2.Sf.S)$Chisq[2]), 
                       P=c(anova(m1.Sf.PW, m2.Sf.PW)$Pr[2], anova(m1.Sf.DtP, m2.Sf.DtP)$Pr[2],
                           anova(m1.Sf.S, m2.Sf.S)$Pr[2]))

#---- Cp; Treatment as random effect

#for pupal weights
d.temp <- filter (d.rich, Treatment != "C", Species=="Cp", !is.na(Pupal.weight.ST))
m1.Cp.PW <- lmer(Pupal.weight.ST ~ (1|Treatment) + (1|exp), data=d.temp)
summary(m1.Cp.PW)
vc <- data.frame(VarCorr(m1.Cp.PW))
vc.Cp.PW <- mutate(vc, percVar=vcov/sum(vcov))
vc.Cp.PW
m2.Cp.PW <- lmer(Pupal.weight.ST ~ (1|exp), data=d.temp)
summary(m2.Cp.PW)
anova(m1.Cp.PW, m2.Cp.PW)  #model is better with treatment


#for DtP
d.temp <- filter (d.rich, Treatment != "C", Species=="Cp", !is.na(Days.to.pupation.ST))
m1.Cp.DtP <- lmer(Days.to.pupation.ST.inv ~ (1|Treatment) + (1|exp) , data=d.temp)
summary(m1.Cp.DtP)
vc <- data.frame(VarCorr(m1.Cp.DtP))
vc.Cp.DtP <- mutate(vc, percVar=vcov/sum(vcov))
vc.Cp.DtP   #6.7% of variance explained by treatment
m2.Cp.DtP <- lmer(Days.to.pupation.ST.inv ~ (1|exp), data=d.temp)
summary(m2.Cp.DtP)
anova(m1.Cp.DtP, m2.Cp.DtP)  #model is better with treatment

#for S
d.temp <- filter (d.rich, Treatment != "C", Species=="Cp")
m1.Cp.S <- glmer(Surv ~ (1|Treatment) + (1|exp), data=d.temp, family=binomial)
summary(m1.Cp.S)
vc <- data.frame(VarCorr(m1.Cp.S))
icc <- as.numeric(icc(m1.Cp.S)[1]) #estimates % of variance explained by all random effects combined
vc.Cp.S <- mutate(vc, percVar=vcov*icc/sum(vcov))
vc.Cp.S   #3.28% of variance explained by treatment
m2.Cp.S <- glmer(Surv ~ (1|exp), data=d.temp, family=binomial)
summary(m2.Cp.S)
anova(m1.Cp.S, m2.Cp.S)  #model is NOT better with treatment

tbl.S10.Cp <- data.frame(Sp="Cp", metric=c("PW", "DtP", "S"), 
                       var=c(vc.Cp.PW$vcov[1], vc.Cp.DtP$vcov[1], vc.Cp.S$vcov[1]), 
                       percVar=c(vc.Cp.PW$percVar[1], vc.Cp.DtP$percVar[1], vc.Cp.S$percVar[1]), 
                       Chisq=c(anova(m1.Cp.PW, m2.Cp.PW)$Chisq[2],anova(m1.Cp.DtP, m2.Cp.DtP)$Chisq[2],
                               anova(m1.Cp.S, m2.Cp.S)$Chisq[2]), 
                       P=c(anova(m1.Cp.PW, m2.Cp.PW)$Pr[2], anova(m1.Cp.DtP, m2.Cp.DtP)$Pr[2],
                           anova(m1.Cp.S, m2.Cp.S)$Pr[2]))

#---- Px; Treatment as random effect

#for pupal weights
d.temp <- filter (d.rich, Treatment != "C", Species=="Px", !is.na(Pupal.weight.ST))
m1.Px.PW <- lmer(Pupal.weight.ST ~ (1|Treatment) + (1|exp), data=d.temp)
summary(m1.Px.PW)
vc <- data.frame(VarCorr(m1.Px.PW))
vc.Px.PW <- mutate(vc, percVar=vcov/sum(vcov))
vc.Px.PW
m2.Px.PW <- lmer(Pupal.weight.ST ~ (1|exp), data=d.temp)
summary(m2.Px.PW)
anova(m1.Px.PW, m2.Px.PW)  #model is better with treatment


#for DtP
d.temp <- filter (d.rich, Treatment != "C", Species=="Px", !is.na(Days.to.pupation.ST))
m1.Px.DtP <- lmer(Days.to.pupation.ST.inv ~ (1|Treatment) + (1|exp) , data=d.temp)
summary(m1.Px.DtP)
vc <- data.frame(VarCorr(m1.Px.DtP))
vc.Px.DtP <- mutate(vc, percVar=vcov/sum(vcov))
vc.Px.DtP   #2.8% of variance explained by treatment
m2.Px.DtP <- lmer(Days.to.pupation.ST.inv ~ (1|exp), data=d.temp)
summary(m2.Px.DtP)
anova(m1.Px.DtP, m2.Px.DtP)  #model is better with treatment

#for S
d.temp <- filter (d.rich, Treatment != "C", Species=="Px")
m1.Px.S <- glmer(Surv ~ (1|Treatment) + (1|exp), data=d.temp, family=binomial)
summary(m1.Px.S)
vc <- data.frame(VarCorr(m1.Px.S))
icc <- as.numeric(icc(m1.Px.S)[1]) #estimates % of variance explained by all random effects combined
vc.Px.S <- mutate(vc, percVar=vcov*icc/sum(vcov))
vc.Px.S   #29% of variance explained by treatment
m2.Px.S <- glmer(Surv ~ (1|exp), data=d.temp, family=binomial)
summary(m2.Px.S)
anova(m1.Px.S, m2.Px.S)  #model is better with treatment

tbl.S10.Px <- data.frame(Sp="Px", metric=c("PW", "DtP", "S"), 
                       var=c(vc.Px.PW$vcov[1], vc.Px.DtP$vcov[1], vc.Px.S$vcov[1]), 
                       percVar=c(vc.Px.PW$percVar[1], vc.Px.DtP$percVar[1], vc.Px.S$percVar[1]), 
                       Chisq=c(anova(m1.Px.PW, m2.Px.PW)$Chisq[2],anova(m1.Px.DtP, m2.Px.DtP)$Chisq[2],
                               anova(m1.Px.S, m2.Px.S)$Chisq[2]), 
                       P=c(anova(m1.Px.PW, m2.Px.PW)$Pr[2], anova(m1.Px.DtP, m2.Px.DtP)$Pr[2],
                           anova(m1.Px.S, m2.Px.S)$Pr[2]))



#Now adding the fungi

#---- Botrys; Treatment as random effect

d.temp <- filter (d.fungi2, Treatment != "DMSO", Fungi=="Botrys", !is.na(dAbs.ST))
m1.B <- lmer(dAbs.ST ~ (1|Treatment) + (1|Plate), data=d.temp)
summary(m1.B)
vc <- data.frame(VarCorr(m1.B))
vc.B <- mutate(vc, percVar=vcov/sum(vcov))
vc.B
m2.B <- lmer(dAbs.ST ~ (1|Plate), data=d.temp)
summary(m2.B)
anova(m1.B, m2.B)  #model is better with treatment


#---- Collet; Treatment as random effect

d.temp <- filter (d.fungi2, Treatment != "DMSO", Fungi=="Collet", !is.na(dAbs.ST))
m1.C <- lmer(dAbs.ST ~ (1|Treatment) + (1|Plate), data=d.temp)
summary(m1.C)
vc <- data.frame(VarCorr(m1.C))
vc.C <- mutate(vc, percVar=vcov/sum(vcov))
vc.C
m2.C <- lmer(dAbs.ST ~ (1|Plate), data=d.temp)
summary(m2.C)
anova(m1.C, m2.C)  #model is better with treatment

#---- Penicillium; Treatment as random effect

d.temp <- filter (d.fungi2, Treatment != "DMSO", Fungi=="Penicillium", !is.na(dAbs.ST))
m1.P <- lmer(dAbs.ST ~ (1|Treatment) + (1|Plate), data=d.temp)
summary(m1.P)
vc <- data.frame(VarCorr(m1.P))
vc.P <- mutate(vc, percVar=vcov/sum(vcov))
vc.P
m2.P <- lmer(dAbs.ST ~ (1|Plate), data=d.temp)
summary(m2.P)
anova(m1.P, m2.P)  #model is better with treatment

#---- Sclerotinia; Treatment as random effect

d.temp <- filter (d.fungi2, Treatment != "DMSO", Fungi=="Sclerotinia", !is.na(dAbs.ST))
m1.S <- lmer(dAbs.ST ~ (1|Treatment) + (1|Plate), data=d.temp)
summary(m1.S)
vc <- data.frame(VarCorr(m1.S))
vc.S <- mutate(vc, percVar=vcov/sum(vcov))
vc.S
m2.S <- lmer(dAbs.ST ~ (1|Plate), data=d.temp)
summary(m2.S)
anova(m1.S, m2.S)  #model is better with treatment


tbl.S10.F <- data.frame(Sp=c("B", "C", "P", "S"), metric="GR", 
                       var=c(vc.B$vcov[1], vc.C$vcov[1], vc.P$vcov[1], vc.S$vcov[1]), 
                       percVar=c(vc.B$percVar[1], vc.C$percVar[1], vc.P$percVar[1],  vc.S$percVar[1]), 
                       Chisq=c(anova(m1.B, m2.B)$Chisq[2],anova(m1.C, m2.C)$Chisq[2],
                               anova(m1.P, m2.P)$Chisq[2], anova(m1.S, m2.S)$Chisq[2]), 
                       P=c(anova(m1.B, m2.B)$Pr[2], anova(m1.C, m2.C)$Pr[2],
                           anova(m1.P, m2.P)$Pr[2], anova(m1.S, m2.S)$Pr[2]))

tbl.S10 <- rbind(tbl.S10.Hz, tbl.S10.Sf, tbl.S10.Cp, tbl.S10.Px, tbl.S10.F)
tbl.S10 <- tbl.S10[c(1,4,7,10,2,5,8,11,3,6,9,12, 13:16),]
tbl.S10$percVar <- tbl.S10$percVar*100

write.csv(tbl.S10, "./Outputs/Tables/TableS10_Pred3c_TreatmentEffects.csv")

#**************************************************************





# Supplemental: Effects of evenness and SD on performance----
#--------------------

load("./Outputs/Workspaces/StandardizedData")

#---- all species combined

##Big picture first: looking at overall effects of evenness and SD on standardized
#performance metrics, all species combined

##Now looking at overall effects of diversity using standardized variables

#Looking at effects of Evenness, SD, Species, and Sex on pupal weights
d.temp <- filter (d.even, Treatment != "C", !is.na(Sex), !is.na(Pupal.weight.ST))
m1.E.all.PW <- lmer(Pupal.weight.ST ~ SD*Evenness*Species*Sex + (1|Treatment) + (1|exp/TrayID/CupID), 
                    data = d.temp, na.action=na.fail, REML="FALSE")
d1.E.all.PW <- dredge(m1.E.all.PW)
s <- subset(d1.E.all.PW, delta < 4)
s  #use for Table S12
d1.E.all.PW.avg <- model.avg(d1.E.all.PW, subset = delta < 4, fit=TRUE)
summary(d1.E.all.PW.avg)  #use for Table S13
m1.E.top.PW <- lmer(Pupal.weight.ST ~ Evenness + SD + Sex + Species + Evenness*Sex + Evenness*Species +
                      Sex*Species + Evenness*Sex*Species + (1|Treatment) + (1|exp/TrayID/CupID), 
                    data = d.temp, na.action=na.fail) 
summary(m1.E.top.PW)
drop1(m1.E.top.PW, test="Chisq")

plot(d.temp$Pupal.weight.ST ~ d.temp$Evenness)
abline(lm(d.temp$Pupal.weight.ST ~ d.temp$Evenness))

#Significant 3-way interaction evenness*sex*species, also model averaging suggests
#effects of sex*species, species, and sex. No overall effects of evenness


#Effects on days to pupation
d.temp <- filter (d.even, Treatment != "C", !is.na(Sex), !is.na(Days.to.pupation.ST.inv))
m1.E.all.dTp <- lmer(Days.to.pupation.ST.inv ~ SD*Evenness*Species*Sex + (1|Treatment) + (1|exp/TrayID/CupID), 
                     data = d.temp, na.action=na.fail, REML="FALSE")
d1.E.all.dTp <- dredge(m1.E.all.dTp)
s <- subset(d1.E.all.dTp, delta < 4)
s  #use for Table S12
d1.E.all.dTp.avg <- model.avg(d1.E.all.dTp, subset = delta < 4, fit=TRUE)
summary(d1.E.all.dTp.avg) #Use for Table S13

m1.E.top.dTp <- lmer(Days.to.pupation.ST.inv ~ Evenness*Species + SD*Species + Sex*Species + (1|Treatment) + (1|exp/TrayID/CupID), 
                     data = d.temp, na.action=na.fail) 
summary(m1.E.top.dTp)
drop1(m1.E.top.dTp, test="Chisq")

plot(d.temp$Days.to.pupation.ST.inv ~ d.temp$Evenness)
abline(lm(d.temp$Days.to.pupation.ST.inv ~ d.temp$Evenness))

#some effects of species
#Strong Evenness*Species and Evenness*SD interactions 
##Evenness*SD*Species interaction
#Evenness is in every model, but effect is sometimes neg and sometimes pos


##Effects on survival

d.temp <- filter (d.even.PS, Treatment != "C", !is.na(PropSurv.ST))
m1.E.all.S <- lmer(PropSurv.ST ~ SD*Evenness*Species + (1|Treatment) + (1|exp), 
                   data = d.temp, na.action=na.fail, REML="FALSE")
d1.E.all.S <- dredge(m1.E.all.S)
s <- subset(d1.E.all.S, delta < 4)
s   #Use for Table S12
d1.E.all.S.avg <- model.avg(d1.E.all.S, subset = delta < 4, fit=TRUE)
summary(d1.E.all.S.avg)  #Use for Table S13

#OVerall only species is coming out as significant. Evenness is in 3/5 top models, sometimes
#positive sometimes negative


#Considering the many interactions with species and sex, will also conduct
#these analyses separately for each herbivore species and performance metric for m and f

#----  Hz  

#making an empty table to store model results
m.sum.E <- data.frame(exsd_chi=NA, exsd_p=NA, e_chi=NA, e_p=NA, sd_chi=NA, sd_p=NA)

#female pupal weights
d.temp <- filter (d.even, Species=="Hz" & Treatment != "C" & Sex=="f")
m1.E.Hz.PW.f <- lmer(Pupal.weight.ST ~ SD*Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
m2.E.Hz.PW.f <- lmer(Pupal.weight.ST ~ SD + Evenness + (1|Treatment)+ (1|exp/TrayID), data = d.temp)
d1 <- drop1(m1.E.Hz.PW.f, test="Chisq")
d2 <- drop1(m2.E.Hz.PW.f, test="Chisq")  #no effects
d1
d2
m.sum.E[1,] <- c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4])

#male pupal weights
d.temp <- filter (d.even, Species=="Hz" & Treatment != "C" & Sex=="m")
m1.E.Hz.PW.m <- lmer(Pupal.weight.ST ~ SD*Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
m2.E.Hz.PW.m <- lmer(Pupal.weight.ST ~ SD + Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
d1 <- drop1(m1.E.Hz.PW.m, test="Chisq")
d2 <- drop1(m2.E.Hz.PW.m, test="Chisq")  #marginal SD x Even interaction
d1
d2
m.sum.E <- rbind(m.sum.E, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))

#female days to pupation
d.temp <- filter (d.even, Species=="Hz" & Treatment != "C" & Sex=="f")
m1.E.Hz.dTp.f <- lmer(Days.to.pupation.ST.inv ~ SD*Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
m2.E.Hz.dTp.f <- lmer(Days.to.pupation.ST.inv ~ SD + Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
d1 <- drop1(m1.E.Hz.dTp.f, test="Chisq")
d2 <- drop1(m2.E.Hz.dTp.f, test="Chisq")  #significant effect of SD
d1
d2
m.sum.E <- rbind(m.sum.E, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))  

#Tukey post-hoc for SD levels
m2.E.T <- glht(m2.E.Hz.dTp.f, linfct=mcp(SD="Tukey"))
summary(m2.E.T)
cld(m2.E.T, level = 0.05) #development faster on high SD compared to low or med
plot(d.temp$Days.to.pupation.ST.inv ~ d.temp$SD)

#male days to pupation
d.temp <- filter (d.even, Species=="Hz" & Treatment != "C" & Sex=="m")
m1.E.Hz.dTp.m <- lmer(Days.to.pupation.ST.inv ~ SD*Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
m2.E.Hz.dTp.m <- lmer(Days.to.pupation.ST.inv ~ SD + Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
d1 <- drop1(m1.E.Hz.dTp.m, test="Chisq")
d2 <- drop1(m2.E.Hz.dTp.m, test="Chisq") #significant effect of SD
d1
d2
m.sum.E <- rbind(m.sum.E, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))

#Tukey post-hoc for SD levels
m2.E.T <- glht(m2.E.Hz.dTp.m, linfct=mcp(SD="Tukey"))
summary(m2.E.T)
cld(m2.E.T, level = 0.05) #development faster on high SD compared to low or med
plot(d.temp$Days.to.pupation.ST.inv ~ d.temp$SD)

#survival
d.temp <- filter (d.even.PS, Species=="Hz" & Treatment != "C")
m1.E.Hz.Surv <- lmer(PropSurv.ST ~ SD*Evenness + (1|exp), data = d.temp)
m2.E.Hz.Surv <- lmer(PropSurv.ST ~ SD+Evenness + (1|exp), data = d.temp)
d1 <- drop1(m1.E.Hz.Surv, test="Chisq")
d2 <- drop1(m2.E.Hz.Surv, test="Chisq")  #no effects
d1
d2
m.sum.E <- rbind(m.sum.E, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))





#----  Sf 

#female pupal weights
d.temp <- filter (d.even, Species=="Sf" & Treatment != "C" & Sex=="f")
m1.E.Sf.PW.f <- lmer(Pupal.weight.ST ~ SD*Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
m2.E.Sf.PW.f <- lmer(Pupal.weight.ST ~ SD + Evenness + (1|Treatment)+ (1|exp/TrayID), data = d.temp)
d1 <- drop1(m1.E.Sf.PW.f, test="Chisq")
d2 <- drop1(m2.E.Sf.PW.f, test="Chisq")  
d1
d2  #no effects
m.sum.E <- rbind(m.sum.E, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))


#male pupal weights
d.temp <- filter (d.even, Species=="Sf" & Treatment != "C" & Sex=="m")
m1.E.Sf.PW.m <- lmer(Pupal.weight.ST ~ SD*Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
m2.E.Sf.PW.m <- lmer(Pupal.weight.ST ~ SD + Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
d1 <- drop1(m1.E.Sf.PW.m, test="Chisq")
d2 <- drop1(m2.E.Sf.PW.m, test="Chisq")  #significant positive effect of evenness
d1
d2
m.sum.E <- rbind(m.sum.E, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))

summary(m2.E.Sf.PW.m)


#female days to pupation
d.temp <- filter (d.even, Species=="Sf" & Treatment != "C" & Sex=="f")
m1.E.Sf.dTp.f <- lmer(Days.to.pupation.ST.inv ~ SD*Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
m2.E.Sf.dTp.f <- lmer(Days.to.pupation.ST.inv ~ SD  + Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
d1 <- drop1(m1.E.Sf.dTp.f, test="Chisq")
d2 <- drop1(m2.E.Sf.dTp.f, test="Chisq")  #significant effect of SD
d1
d2
m.sum.E <- rbind(m.sum.E, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))

#Tukey post-hoc for SD levels
m2.E.T <- glht(m2.E.Sf.dTp.f, linfct=mcp(SD="Tukey"))
summary(m2.E.T)
cld(m2.E.T, level = 0.05) #development faster on high SD compared to low or med
plot(d.temp$Days.to.pupation.ST.inv ~ d.temp$SD)


#male days to pupation
d.temp <- filter (d.even, Species=="Sf" & Treatment != "C" & Sex=="m")
m1.E.Sf.dTp.m <- lmer(Days.to.pupation.ST.inv ~ SD*Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
m2.E.Sf.dTp.m <- lmer(Days.to.pupation.ST.inv ~ SD + Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
d1 <- drop1(m1.E.Sf.dTp.m, test="Chisq")
d2 <- drop1(m2.E.Sf.dTp.m, test="Chisq") #marginal effect of SD
d1
d2
m.sum.E <- rbind(m.sum.E, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))

#Tukey post-hoc for SD levels
m2.E.T <- glht(m2.E.Sf.dTp.m, linfct=mcp(SD="Tukey"))
summary(m2.E.T)
cld(m2.E.T, level = 0.05) #development marginally faster on high SD compared to low; med is intermediate
plot(d.temp$Days.to.pupation.ST.inv ~ d.temp$SD)

#survival
d.temp <- filter (d.even.PS, Species=="Sf" & Treatment != "C")
m1.E.Sf.Surv <- lmer(PropSurv.ST ~ SD*Evenness + (1|exp), data = d.temp)
m2.E.Sf.Surv <- lmer(PropSurv.ST ~ SD+Evenness + (1|exp), data = d.temp)
d1 <- drop1(m1.E.Sf.Surv, test="Chisq")
d2 <- drop1(m2.E.Sf.Surv, test="Chisq")  #no effects
d1
d2
m.sum.E <- rbind(m.sum.E, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))




#----  Cp 

#female pupal weights
d.temp <- filter (d.even, Species=="Cp" & Treatment != "C" & Sex=="f")
m1.E.Cp.PW.f <- lmer(Pupal.weight.ST ~ SD*Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
m2.E.Cp.PW.f <- lmer(Pupal.weight.ST ~ SD + Evenness + (1|Treatment)+ (1|exp/TrayID), data = d.temp)
d1 <- drop1(m1.E.Cp.PW.f, test="Chisq")
d2 <- drop1(m2.E.Cp.PW.f, test="Chisq")  #no effects
d1
d2
m.sum.E <- rbind(m.sum.E, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))

#male pupal weights
d.temp <- filter (d.even, Species=="Cp" & Treatment != "C" & Sex=="m")
m1.E.Cp.PW.m <- lmer(Pupal.weight.ST ~ SD*Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
m2.E.Cp.PW.m <- lmer(Pupal.weight.ST ~ SD + Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
d1 <- drop1(m1.E.Cp.PW.m, test="Chisq")
d2 <- drop1(m2.E.Cp.PW.m, test="Chisq")  #no effects
d1
d2
m.sum.E <- rbind(m.sum.E, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))



#female days to pupation
d.temp <- filter (d.even, Species=="Cp" & Treatment != "C" & Sex=="f")
m1.E.Cp.dTp.f <- lmer(Days.to.pupation.ST.inv ~ SD*Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
m2.E.Cp.dTp.f <- lmer(Days.to.pupation.ST.inv ~ SD + Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
d1 <- drop1(m1.E.Cp.dTp.f, test="Chisq")
d2 <- drop1(m2.E.Cp.dTp.f, test="Chisq")  #marginal interaction
d1
d2
m.sum.E <- rbind(m.sum.E, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))


#male days to pupation
d.temp <- filter (d.even, Species=="Cp" & Treatment != "C" & Sex=="m")
m1.E.Cp.dTp.m <- lmer(Days.to.pupation.ST.inv ~ SD*Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
m2.E.Cp.dTp.m <- lmer(Days.to.pupation.ST.inv ~ SD + Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
d1 <- drop1(m1.E.Cp.dTp.m, test="Chisq")
d2 <- drop1(m2.E.Cp.dTp.m, test="Chisq") #marginal negative effect of evenness
d1
d2
m.sum.E <- rbind(m.sum.E, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))

summary(m2.E.Cp.dTp.m)


#survival
d.temp <- filter (d.even.PS, Species=="Cp" & Treatment != "C")
m1.E.Cp.Surv <- lmer(PropSurv.ST ~ SD*Evenness + (1|exp), data = d.temp)
m2.E.Cp.Surv <- lmer(PropSurv.ST ~ SD+Evenness + (1|exp), data = d.temp)
d1 <- drop1(m1.E.Cp.Surv, test="Chisq")
d2 <- drop1(m2.E.Cp.Surv, test="Chisq")  
d1
d2   #no effects

m.sum.E <- rbind(m.sum.E, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))


#----------Px


#female pupal weights
d.temp <- filter (d.even, Species=="Px" & Treatment != "C" & Sex=="f")
m1.E.Px.PW.f <- lmer(Pupal.weight.ST ~ SD*Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
m2.E.Px.PW.f <- lmer(Pupal.weight.ST ~ SD + Evenness + (1|Treatment)+ (1|exp/TrayID), data = d.temp)
d1 <- drop1(m1.E.Px.PW.f, test="Chisq")
d2 <- drop1(m2.E.Px.PW.f, test="Chisq")  #no effects
d1
d2
m.sum.E <- rbind(m.sum.E, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))

#male pupal weights
d.temp <- filter (d.even, Species=="Px" & Treatment != "C" & Sex=="m")
m1.E.Px.PW.m <- lmer(Pupal.weight.ST ~ SD*Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
m2.E.Px.PW.m <- lmer(Pupal.weight.ST ~ SD + Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
d1 <- drop1(m1.E.Px.PW.m, test="Chisq")
d2 <- drop1(m2.E.Px.PW.m, test="Chisq")  #significant negative effect of evenness
d1
d2
m.sum.E <- rbind(m.sum.E, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))

summary(m2.E.Px.PW.m)

#female days to pupation
d.temp <- filter (d.even, Species=="Px" & Treatment != "C" & Sex=="f")
m1.E.Px.dTp.f <- lmer(Days.to.pupation.ST.inv ~ SD*Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
m2.E.Px.dTp.f <- lmer(Days.to.pupation.ST.inv ~ SD + Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
d1 <- drop1(m1.E.Px.dTp.f, test="Chisq")
d2 <- drop1(m2.E.Px.dTp.f, test="Chisq")  #no effects
d1
d2
m.sum.E <- rbind(m.sum.E, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))


#male days to pupation
d.temp <- filter (d.even, Species=="Px" & Treatment != "C" & Sex=="m")
m1.E.Px.dTp.m <- lmer(Days.to.pupation.ST.inv ~ SD*Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
m2.E.Px.dTp.m <- lmer(Days.to.pupation.ST.inv ~ SD + Evenness + (1|Treatment) + (1|exp/TrayID), data = d.temp)
d1 <- drop1(m1.E.Px.dTp.m, test="Chisq")
d2 <- drop1(m2.E.Px.dTp.m, test="Chisq") #significant negative effect of evenness
d1
d2
summary(m2.E.Px.dTp.m)
m.sum.E <- rbind(m.sum.E, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))



#survival
d.temp <- filter (d.even.PS, Species=="Px" & Treatment != "C")
m1.E.Px.Surv <- lmer(PropSurv.ST ~ SD*Evenness + (1|exp), data = d.temp)
m2.E.Px.Surv <- lmer(PropSurv.ST ~ SD+Evenness + (1|exp), data = d.temp)
d1 <- drop1(m1.E.Px.Surv, test="Chisq")
d2 <- drop1(m2.E.Px.Surv, test="Chisq") #no effects
d1
d2
m.sum.E <- rbind(m.sum.E, c(d1[2,3], d1[2,4], d2[3,3], d2[3,4], d2[2,3], d2[2,4]))

#Save results table
write.csv(m.sum.E, file="./Outputs/Tables/TableS14_Evenness_Performance.csv")

#Save workspace
save.image("./Outputs/Workspaces/FinalAnalyses_Evenness_Performance")

#**************************************************************





# Supplemental: Effects of evenness on number of effects----
#--------------------

load("./Outputs/Workspaces/StandardizedData")


d.temp <- d.even[which(d.even$Treatment != "C"),]

d.sum <- d.temp %>% 
  group_by(Evenness, SD, Treatment, Species) %>%
  summarise(PW.avg = mean(Pupal.weight.ST, na.rm=TRUE), PW.SD=sd(Pupal.weight.ST, na.rm=TRUE),
            PW.n = sum(!is.na(Pupal.weight.ST)), DtP.avg= mean(Days.to.pupation.ST.inv, na.rm=TRUE), 
            DtP.SD=sd(Days.to.pupation.ST.inv, na.rm=TRUE), DtP.n = sum(!is.na(Days.to.pupation.ST.inv)))

d.sum$PW.CI.low <- d.sum$PW.avg - (1.96*d.sum$PW.SD/sqrt(d.sum$PW.n))
d.sum$PW.CI.high <- d.sum$PW.avg + (1.96*d.sum$PW.SD/sqrt(d.sum$PW.n))

d.sum$DtP.CI.low <- d.sum$DtP.avg - (1.96*d.sum$DtP.SD/sqrt(d.sum$DtP.n))
d.sum$DtP.CI.high <- d.sum$DtP.avg + (1.96*d.sum$DtP.SD/sqrt(d.sum$DtP.n))


##Add binary variable that indicates whether there is a significant effect (i.e. CI
#crosses 1) for each treatment

d.sum$TxMoreRes.PW <- NA
d.sum$TxMoreRes.DtP <- NA
for(i in 1:length(d.sum$TxMoreRes.PW)){ 
  d.sum$TxMoreRes.PW[i] <-ifelse(d.sum$PW.CI.high[i] < 1, 1, 0)
  d.sum$TxMoreRes.DtP[i] <-ifelse(d.sum$DtP.CI.high[i] < 1, 1, 0)
}   


#For survival, the data are already summarized by treatment

d.temp.PS <- d.even.PS[which(d.even.PS$Treatment !="C"),]
d.temp.PS.C <- d.even.PS[which(d.even.PS$Treatment=="C"),]

d.sum$S.Tx.alive <- NA
d.sum$S.Tx.dead <- NA
d.sum$S.C.alive <- NA
d.sum$S.C.dead <- NA
d.sum$S.pval <- NA
d.sum$TxMoreRes.S <- NA

for(i in 1:length(d.temp.PS$Treatment)){
  Sp <- d.temp.PS$Species[i]
  Tx <- d.temp.PS$Treatment[i]
  exp <- d.temp.PS$exp[i]
  Cs <- d.temp.PS.C[which(d.temp.PS.C$exp==exp & d.temp.PS.C$Species==Sp),]
  
  #Is the prob of survival different in Tx vs C?
  Chi.tbl <- rbind(c(d.temp.PS$dead[i], d.temp.PS$alive[i]),
                   c(Cs$dead, Cs$alive))
  p <- fisher.test(Chi.tbl)$p.value
  
  #Is prob of survival higher in C? 1=yes, 0=no
  TxHigh <- ifelse(d.temp.PS$dead[i]/(d.temp.PS$dead[i] + d.temp.PS$alive[i]) >
                     Cs$dead/(Cs$dead + Cs$alive), 1, 0)
  TxMoreRes <- ifelse(p < 0.05 & TxHigh==1, 1, 0)
  
  #Add values to d.sum table
  d.sum$S.Tx.alive[which(d.sum$Species==Sp & d.sum$Treatment==Tx)] <- d.temp.PS$alive[i]
  d.sum$S.Tx.dead[which(d.sum$Species==Sp & d.sum$Treatment==Tx)]  <- d.temp.PS$dead[i]
  d.sum$S.C.alive[which(d.sum$Species==Sp & d.sum$Treatment==Tx)]  <- Cs$alive
  d.sum$S.C.dead[which(d.sum$Species==Sp & d.sum$Treatment==Tx)]  <- Cs$dead
  d.sum$S.pval[which(d.sum$Species==Sp & d.sum$Treatment==Tx)]  <- p
  d.sum$TxMoreRes.S[which(d.sum$Species==Sp & d.sum$Treatment==Tx)]  <- TxMoreRes
}


d.sum$HerbAffected <- ifelse(d.sum$TxMoreRes.PW==1 | d.sum$TxMoreRes.DtP==1 |d.sum$TxMoreRes.S==1, 1, 0)
d.sum <- mutate(d.sum, NumEffects=TxMoreRes.PW + TxMoreRes.DtP + TxMoreRes.S)


d.sum2 <- d.sum %>%
  group_by(Evenness, SD, Treatment) %>%
  summarise(NumHerbsAffected=sum(HerbAffected, na.rm=TRUE), NumEffectsTot=sum(NumEffects, na.rm=TRUE))


#Some plots
plot(jitter(NumHerbsAffected) ~ jitter(Evenness), data=d.sum2)
abline(lm(NumHerbsAffected ~ Evenness, data=d.sum2))

plot(jitter(NumEffectsTot) ~ jitter(Evenness), data=d.sum2)
abline(lm(NumEffectsTot ~ Evenness, data=d.sum2))


#checking distribution of response variable
hist(d.sum2$NumHerbsAffected)
shapiro.test(d.sum2$NumHerbsAffected)
#not normal, but also does not approximate typical count distributions
#like poisson or negative binomial

#checking different models
m1 <- lm(NumHerbsAffected ~ Evenness, data=d.sum2)
summary(m1)
m2 <- glm(NumHerbsAffected ~ Evenness, data=d.sum2, family=poisson)
summary(m2)
m3 <- glm(NumHerbsAffected ~ Evenness, data=d.sum2, family=quasipoisson)
summary(m3)  #no AIC score

plot(m1)
plot(m2)
plot(m3)

#again the normal distribution looks the best
#but AIC lower with gaussian than poisson

hist(resid(m1))  #these actually look pretty good

AIC(m1, m2, m3)

#Seems the best test is the gaussian, but all seem to be giving qualitatively similar results
#Use gaussian for results
summary(m1)

#---Effects of SD on number of effects

#Some plots
plot(NumHerbsAffected ~ SD, data=d.sum2)


m1 <- lm(NumHerbsAffected ~ SD, data=d.sum2)
summary(m1)
#no differences

#----checking means
d.sum2.means2 <- d.sum2 %>%
  group_by(Evenness) %>%
  summarise(mean=mean(NumHerbsAffected, na.rm=TRUE))
d.sum2.means2  


#Save workspace for fig
NumEffects.even <- d.sum2 
save.image("./Outputs/Workspaces/FinalAnalyses_Evenness_NumEffects")

#**************************************************************




# Supplemental: Molecular weight and bioactivity----
#--------------------

##in response to reviews and thinking about whether we can compare 
#compound bioactivities because they were present in different molar concentrations
#we are curious whether larger, more complex compounds might have
#broader overall potential for bioactivity
#Checking whether total number of effects depends on size
#for the individual compounds

load("./Outputs/Workspaces/FinalAnalyses_Pred2a-2b_NumEffects")
MW <- read.csv("Whitehead_et_al_CIDs.csv")

#checking relationship between MW and total number of effects
d.temp <- d.sum2[which(d.sum2$Richness==1),]
colnames(MW)[3] <- "Treatment"
d.temp <- left_join(d.temp, MW, by="Treatment")

plot(TotEffects ~ MW, data=d.temp)
summary(lm(TotEffects ~ MW, data=d.temp))

#can also look at relationship for particular performance metrics

#Sf_PW
d.temp <- d.sum[which(d.sum$Richness==1 & d.sum$Species=="Sf"),]
d.temp <- left_join(d.temp, MW, by="Treatment")
plot(PW.avg ~ MW, data=d.temp)
summary(lm(PW.avg ~ MW, data=d.temp))

#Hz_PW
d.temp <- d.sum[which(d.sum$Richness==1 & d.sum$Species=="Hz"),]
d.temp <- left_join(d.temp, MW, by="Treatment")
plot(PW.avg ~ MW, data=d.temp)
summary(lm(PW.avg ~ MW, data=d.temp))

#Cp_PW
d.temp <- d.sum[which(d.sum$Richness==1 & d.sum$Species=="Cp"),]
d.temp <- left_join(d.temp, MW, by="Treatment")
plot(PW.avg ~ MW, data=d.temp)
summary(lm(PW.avg ~ MW, data=d.temp))

#Px_PW
d.temp <- d.sum[which(d.sum$Richness==1 & d.sum$Species=="Px"),]
d.temp <- left_join(d.temp, MW, by="Treatment")
plot(PW.avg ~ MW, data=d.temp)
summary(lm(PW.avg ~ MW, data=d.temp))


#Sf_DtP
d.temp <- d.sum[which(d.sum$Richness==1 & d.sum$Species=="Sf"),]
d.temp <- left_join(d.temp, MW, by="Treatment")
plot(DtP.avg ~ MW, data=d.temp)
summary(lm(DtP.avg ~ MW, data=d.temp))

#Hz_DtP
d.temp <- d.sum[which(d.sum$Richness==1 & d.sum$Species=="Hz"),]
d.temp <- left_join(d.temp, MW, by="Treatment")
plot(DtP.avg ~ MW, data=d.temp)
summary(lm(DtP.avg ~ MW, data=d.temp))

#Cp_DtP
d.temp <- d.sum[which(d.sum$Richness==1 & d.sum$Species=="Cp"),]
d.temp <- left_join(d.temp, MW, by="Treatment")
plot(DtP.avg ~ MW, data=d.temp)
summary(lm(DtP.avg ~ MW, data=d.temp))

#Px_DtP
d.temp <- d.sum[which(d.sum$Richness==1 & d.sum$Species=="Px"),]
d.temp <- left_join(d.temp, MW, by="Treatment")
plot(DtP.avg ~ MW, data=d.temp)
summary(lm(DtP.avg ~ MW, data=d.temp))



#Sf_S
d.temp <- d.rich.PS[which(d.rich.PS$Richness==1 & d.rich.PS$Species=="Sf"),]
d.temp <- left_join(d.temp, MW, by="Treatment")
plot(PropSurv.ST ~ MW, data=d.temp)
summary(lm(PropSurv.ST ~ MW, data=d.temp))
#here is one where there is a negative relationships...survival lower
#on larger compounds, even though larger compounds have
#fewer moles

#Hz_S
d.temp <- d.rich.PS[which(d.rich.PS$Richness==1 & d.rich.PS$Species=="Hz"),]
d.temp <- left_join(d.temp, MW, by="Treatment")
plot(PropSurv.ST ~ MW, data=d.temp)
summary(lm(PropSurv.ST ~ MW, data=d.temp))

#Cp_S
d.temp <- d.rich.PS[which(d.rich.PS$Richness==1 & d.rich.PS$Species=="Cp"),]
d.temp <- left_join(d.temp, MW, by="Treatment")
plot(PropSurv.ST ~ MW, data=d.temp)
summary(lm(PropSurv.ST ~ MW, data=d.temp))

#Px_S
d.temp <- d.rich.PS[which(d.rich.PS$Richness==1 & d.rich.PS$Species=="Px"),]
d.temp <- left_join(d.temp, MW, by="Treatment")
plot(PropSurv.ST ~ MW, data=d.temp)
summary(lm(PropSurv.ST ~ MW, data=d.temp))


#Botrys
d.temp <- d.sum.f[which(d.sum.f$Richness==1 & d.sum.f$Species=="Botrys"),]
d.temp <- left_join(d.temp, MW, by="Treatment")
plot(dAbs.avg ~ MW, data=d.temp)
summary(lm(dAbs.avg ~ MW, data=d.temp))

#Collet
d.temp <- d.sum.f[which(d.sum.f$Richness==1 & d.sum.f$Species=="Collet"),]
d.temp <- left_join(d.temp, MW, by="Treatment")
plot(dAbs.avg ~ MW, data=d.temp)
summary(lm(dAbs.avg ~ MW, data=d.temp))

#Penicillium
d.temp <- d.sum.f[which(d.sum.f$Richness==1 & d.sum.f$Species=="Penicillium"),]
d.temp <- left_join(d.temp, MW, by="Treatment")
plot(dAbs.avg ~ MW, data=d.temp)
summary(lm(dAbs.avg ~ MW, data=d.temp))

#Sclerotinia
d.temp <- d.sum.f[which(d.sum.f$Richness==1 & d.sum.f$Species=="Sclerotinia"),]
d.temp <- left_join(d.temp, MW, by="Treatment")
plot(dAbs.avg ~ MW, data=d.temp)
summary(lm(dAbs.avg ~ MW, data=d.temp))

#****************************************************************
