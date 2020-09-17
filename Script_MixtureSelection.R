#This script uses ChemmineR to calculate structural similarity among
#compounds, and then semi-randomly selects mixtures for the experiment
#at each level of richness/evenness and structural diversity

#ChemmineR and fmcsR are bioconductor packages. They need to be sourced
#from there--they are not available on the regular CRAN mirrors
# See https://bioconductor.org/install/#install-bioconductor-packages

#to first install, run...
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()

#BiocManager::install(c("ChemmineR", "fmcsR"))


library(ChemmineR) 
library(fmcsR) #add on package for Maximum Common Substructure (MCS) searching
library(abind)  #For nice pasting together of matrices
vignette("ChemmineR") #Opens this PDF manual from R 
library(dplyr)


#The compound structure data is stored as an SDF file from PubChem
#I got this from: https://pubchem.ncbi.nlm.nih.gov/
#had to upload a list of identifiers for the compounds (e.g. CIDs) as csv file 
#(I used the list in "Whitehead_et_al_CIDs" but took out header and names column)
#then choose "push to entrez" then you can download the structure data as an sdf 
#file with gz compression
sdfset <- read.SDFset("Whitehead_et_al_structures.sdf.gz") 
names <- read.csv("Whitehead_et_al_CIDs.csv")


#Prune compound set to things we can get commercially
sdfset <- sdfset[c(1:9, 12, 14, 18:20)]
names <- names[c(1:9, 12, 14, 18:20),]


##General ways to look at data
plot(sdfset[1:9]) #plots structures in R, can't do them all at once
#sdf.visualize(sdfset) #opens window in browser to view structures in ChemMine online 
sdfid(sdfset) #see list of CIDs in sdfset
propma <- data.frame(ID=names$Compound, MF=MF(sdfset), MW=MW(sdfset), atomcountMA(sdfset))
datablock(sdfset) <- propma  #Assign matrix data to data block:
datablock(sdfset[1:14]) 

grepSDFset("rutin", sdfset, field="datablock", mode="subset")  #to search dataset for any compound, can also use CID or MW, etc
grepSDFset("rutin", sdfset, field="datablock", mode="index") #to get index number of compound




##Different ways to calculate chemical similarity among compounds

##1) MCS 

#pairwise comparison
MCStest <- fmcs(sdfset[1], sdfset[2], au=2, bu=1) # Searches for MCS with mismatches, au specifies number of atom mismatches allowed, bu specifies the number of bond mismatches allowed 
plotMCS(MCStest) # Plots both query compounds with MCS in color 

#batch comparison
fmcsBatch(sdfset[1], sdfset, au=0, bu=0) 

#Compute a similarity matrix for all pairwise comparisons
d.mcs <- sapply(cid(sdfset), function(x) fmcsBatch(sdfset[x], sdfset, au=0, bu=0)[,"Tanimoto_Coefficient"]) 
d.mcs 

#Plot as dendrogram

par(mfrow=c(1,3), mar=c(5.1, 4.1, 4.1, 5.1)) 
hc <- hclust(as.dist(1-d.mcs), method="complete")
plot(as.dendrogram(hc), edgePar=list(col=4, lwd=2), horiz=TRUE, main="MCS") 


###2) Atom Pair Descriptors
apset <- sdf2ap(sdfset) # Generate atom pair descriptor database for searching 
cmp.search(apset, apset[1], type=3, cutoff = 0, quiet=TRUE) # Search apset database to see which other compounds are most similar to 1. 
cmp.cluster(db=apset, cutoff = c(0.65, 0.5), quiet=TRUE) # Binning clustering using variable similarity cutoffs. 
cmp.similarity(apset[1], apset[2])

db.explain(apset[1])

#Compute a similarity matrix for all pairwise comparisons

d.ap <- NULL
for (i in 1:length(apset)) {
  x <- cmp.search(apset, apset[i], type=3, cutoff = 0, quiet=TRUE)
  y <- x[order(x$index),]
  d.ap <- rbind(d.ap, y$scores)
}
colnames(d.ap) <- cid(sdfset)
rownames(d.ap) <- cid(sdfset)
d.ap

#Plot as dendrogram
hc <- hclust(as.dist(1-d.ap), method="complete")
plot(as.dendrogram(hc), edgePar=list(col=4, lwd=2), horiz=TRUE, main="Atom Pair Similarity") 


###3) Fingerprints

fpset <- desc2fp(apset)
view(fpset[1:2]) 
fpSim(fpset[1], fpset, method="Tanimoto", cutoff=0) 



d.fp <- NULL
for (i in 1:length(fpset)) {
  x <- fpSim(fpset[i], fpset, method="Tanimoto", cutoff=0)
  ids <- as.numeric(gsub("CMP","", names(x)))
  names(x) <- ids
  y <- x[order(as.numeric(names(x)))]
  d.fp <- rbind(d.fp, y)
}
colnames(d.fp) <- cid(sdfset)
rownames(d.fp) <- cid(sdfset)
d.fp

#Plot as dendrogram
hc <- hclust(as.dist(1-d.fp), method="complete")
plot(as.dendrogram(hc), edgePar=list(col=4, lwd=2), horiz=TRUE, main="Fingerprints") 


cor(as.vector(d.fp), as.vector(d.ap))
cor(as.vector(d.fp), as.vector(d.mcs))
cor(as.vector(d.ap), as.vector(d.mcs))

d.avg <- (d.fp + d.ap + d.mcs)/3

#d.sub <- d.avg[c(1:9,12, 14, 18:20),c(1:9,12, 14, 18:20)]
d.sub <- d.avg


#Plot as dendrogram
hc <- hclust(as.dist(1-d.avg), method="complete")
plot(as.dendrogram(hc), edgePar=list(col=4, lwd=2), horiz=TRUE, main="Average") 


#save distance matrix for dendrogram for figures
colnames(d.avg) <- names$Compound
rownames(d.avg) <- names$Compound

write.csv(d.avg, "./Outputs/Tables/CompoundDistances.csv")



#############################################################
####For richness experiments-----------selecting mixtures
#########################################

###Set parameters for selection of sets
levels <- c(2,4,6,8,10) #The levels of richness to use; we would also have 1 and 14 but there is only one way to make those
cut <- .10  #the percentage cutoff for what is considered a high, mid, and low structural diversity group
            #candidate sets of 3 will be randomly drawn from this pool, then checked to see if they meet criteria
max.shared <- c(0,2,4,5,7) #set maximum number of compounds that can be shared per group for each richness level
                  #must be the same length as levels
                  ##So at richness=4, only two compounds can be shared per group
max.occur <- c(1,2,2,2,3) #set number of times an individual compound can occur in a set of three for each richness level

IDs <- colnames(d.sub)
all.combos <- list()
sel.combos <- list()
for (i in 1:length(levels)){  #For each richness level
  nr <- max.occur[i]
  g <- max.shared[i] 
  d <- t(combn(IDs, levels[i]))    #create a list of all possible compound combinations
  d <- cbind(d, NA)
  
    for (j in 1:nrow(d)){     #For each compound combo
      #Get the average similarity index across all pw contrasts in combo
      list <- as.character(d[j,1:levels[i]])   
      pw.combos <- combn(list, 2)
      pw.combos <- rbind(pw.combos, NA)
        for (k in 1:ncol(pw.combos)){
          x <- d.sub[as.character(pw.combos[1,k]), as.character(pw.combos[2,k])]  
          pw.combos[3,k] <- x
        }
      z <- mean(as.numeric(pw.combos[3,]))
      d[j,levels[i]+1] <- z
     }
    
  d <- as.data.frame(d) 
  names(d)[levels[i]+1] <- "MeanDist" 
  d2 <- d[order(d$MeanDist),]
  #d2 is a list of all possible combinations ordered by their mean Similarity index for the current richness level

  #Separate out groups for high, mid, and low functional diversity based on specified cutoff
  high <- d2[1:round((nrow(d2)*cut)),]
  
    #for mid group
    up <- round(nrow(d2)/2 + nrow(high)/2)
    down <- round(nrow(d2)/2 - nrow(high)/2)
  mid <- d2[up:down,]
  
  low <- d2[(nrow(d2)-nrow(high)):nrow(d2),]
  
  #Take a random sample of three sets from each group
  high.set <- high[sample(nrow(high), 3, replace=FALSE),]
  mid.set <- mid[sample(nrow(mid), 3, replace=FALSE),]
  low.set <- low[sample(nrow(low), 3, replace=FALSE),] 
    

  ###check whether the random samples of three groups meet the specified criteria

  #First for high
  high.set.final <- 0
  while(high.set.final==0){

    t <- c(as.character(unname(unlist(high.set[1,1:levels[i]]))),
           as.character(unname(unlist(high.set[2,1:levels[i]]))),
           as.character(unname(unlist(high.set[3,1:levels[i]]))))
            
    #make a list of overlaps between each pairwise comparison
    one.two <- intersect(as.character(unname(unlist(high.set[1,1:levels[i]]))),
                         as.character(unname(unlist(high.set[2,1:levels[i]]))))
    two.three <- intersect(as.character(unname(unlist(high.set[2,1:levels[i]]))),
                           as.character(unname(unlist(high.set[3,1:levels[i]]))))
    one.three <- intersect(as.character(unname(unlist(high.set[1,1:levels[i]]))),
                           as.character(unname(unlist(high.set[3,1:levels[i]]))))

     ##check if combos meet criteria for g and nr
      if(all(table(t)<=nr) &&
         length(one.two) <= g  &&
         length(one.three) <= g  &&
         length(two.three) <= g  ){ 
        high.set.final <- high.set  
        } else {
        high.set <- high[sample(nrow(high), 3, replace=FALSE),]
      }
  }
      
  #For mid
  mid.set.final <- 0
  while(mid.set.final==0){
    
    t <- c(as.character(unname(unlist(mid.set[1,1:levels[i]]))),
           as.character(unname(unlist(mid.set[2,1:levels[i]]))),
           as.character(unname(unlist(mid.set[3,1:levels[i]]))))
    
    #make a list of overlaps between each pairwise comparison
    one.two <- intersect(as.character(unname(unlist(mid.set[1,1:levels[i]]))),
                         as.character(unname(unlist(mid.set[2,1:levels[i]]))))
    two.three <- intersect(as.character(unname(unlist(mid.set[2,1:levels[i]]))),
                           as.character(unname(unlist(mid.set[3,1:levels[i]]))))
    one.three <- intersect(as.character(unname(unlist(mid.set[1,1:levels[i]]))),
                           as.character(unname(unlist(mid.set[3,1:levels[i]]))))
    
    ##check if combos meet criteria for g and nr
    if(all(table(t)<=nr) &&
       length(one.two) <= g  &&
       length(one.three) <= g  &&
       length(two.three) <= g  ){ 
      mid.set.final <- mid.set  
    } else {
      mid.set <- mid[sample(nrow(mid), 3, replace=FALSE),]
    }
  }  
  
  #For low
  low.set.final <- 0
  while(low.set.final==0){
    
    t <- c(as.character(unname(unlist(low.set[1,1:levels[i]]))),
           as.character(unname(unlist(low.set[2,1:levels[i]]))),
           as.character(unname(unlist(low.set[3,1:levels[i]]))))
    
    #make a list of overlaps between each pairwise comparison
    one.two <- intersect(as.character(unname(unlist(low.set[1,1:levels[i]]))),
                         as.character(unname(unlist(low.set[2,1:levels[i]]))))
    two.three <- intersect(as.character(unname(unlist(low.set[2,1:levels[i]]))),
                           as.character(unname(unlist(low.set[3,1:levels[i]]))))
    one.three <- intersect(as.character(unname(unlist(low.set[1,1:levels[i]]))),
                           as.character(unname(unlist(low.set[3,1:levels[i]]))))
    
    ##check if combos meet criteria for g and nr
    if(all(table(t)<=nr) &&
       length(one.two) <= g  &&
       length(one.three) <= g  &&
       length(two.three) <= g  ){ 
      low.set.final <- low.set  
    } else {
      low.set <- low[sample(nrow(low), 3, replace=FALSE),]
    }
  }
    
  #Add the  complete sets from high, mid, and low to a list so we can see distributions from which the random
  #samples were drawn
  names.all <- c(names(all.combos), paste("high", levels[i], sep=""),paste("mid", levels[i], sep=""),
                 paste("low", levels[i], sep=""))
  all.combos <- c(all.combos, list(high, mid, low))
  names(all.combos) <- names.all
            
  ##Add the final choices for high, mid, low to final list
  names.sel <- c(names(sel.combos), paste("high", levels[i], sep=""),paste("mid", levels[i], sep=""),
                    paste("low", levels[i], sep=""))    
  sel.combos <- c(sel.combos, list(high.set.final, mid.set.final, low.set.final))
  names(sel.combos) <- names.sel
    
}
sel.combos

low2 <- mutate(sel.combos$low2, Tx=c("L2A", "L2B", "L2C"))
low4 <- mutate(sel.combos$low4, Tx=c("L4A", "L4B", "L4C"))
low6 <- mutate(sel.combos$low6, Tx=c("L6A", "L6B", "L6C"))
low8 <- mutate(sel.combos$low8, Tx=c("L8A", "L8B", "L8C"))
low10 <- mutate(sel.combos$low10, Tx=c("L10A", "L10B", "L10C"))
mid2 <- mutate(sel.combos$mid2, Tx=c("M2A", "M2B", "M2C"))
mid4 <- mutate(sel.combos$mid4, Tx=c("M4A", "M4B", "M4C"))
mid6 <- mutate(sel.combos$mid6, Tx=c("M6A", "M6B", "M6C"))
mid8 <- mutate(sel.combos$mid8, Tx=c("M8A", "M8B", "M8C"))
mid10 <- mutate(sel.combos$mid10, Tx=c("M10A", "M10B", "M10C"))
high2 <- mutate(sel.combos$high2, Tx=c("H2A", "H2B", "H2C"))
high4 <- mutate(sel.combos$high4, Tx=c("H4A", "H4B", "H4C"))
high6 <- mutate(sel.combos$high6, Tx=c("H6A", "H6B", "H6C"))
high8 <- mutate(sel.combos$high8, Tx=c("H8A", "H8B", "H8C"))
high10 <- mutate(sel.combos$high10, Tx=c("H10A", "H10B", "H10C"))


sel.combos.table <- full_join(low2, low4)
sel.combos.table <- full_join(sel.combos.table, low6)
sel.combos.table <- full_join(sel.combos.table, low8)
sel.combos.table <- full_join(sel.combos.table, low10)
sel.combos.table <- full_join(sel.combos.table, mid2)
sel.combos.table <- full_join(sel.combos.table, mid4)
sel.combos.table <- full_join(sel.combos.table, mid6)
sel.combos.table <- full_join(sel.combos.table, mid8)
sel.combos.table <- full_join(sel.combos.table, mid10)
sel.combos.table <- full_join(sel.combos.table, high2)
sel.combos.table <- full_join(sel.combos.table, high4)
sel.combos.table <- full_join(sel.combos.table, high6)
sel.combos.table <- full_join(sel.combos.table, high8)
sel.combos.table <- full_join(sel.combos.table, high10)

sel.combos.table <- sel.combos.table[c(4,1,2,5:12,3)]


#rename compounds with meaningful names
sel.combos.table$V1 <- recode(sel.combos.table$V1, CMP1="ChA", CMP2="CA", CMP3="pCA",
                                     CMP4="FA", CMP5="GA", CMP6="SA",
                                     CMP7="GeA", CMP8="Ct", CMP9="eCt",
                                     CMP12="R",CMP14="H", CMP18="Q",
                                     CMP19="Phz", CMP20="Pht")
sel.combos.table$V2 <- recode(sel.combos.table$V2, CMP1="ChA", CMP2="CA", CMP3="pCA",
                                     CMP4="FA", CMP5="GA", CMP6="SA",
                                     CMP7="GeA", CMP8="Ct", CMP9="eCt",
                                     CMP12="R",CMP14="H", CMP18="Q",
                                     CMP19="Phz", CMP20="Pht")
sel.combos.table$V3 <- recode(sel.combos.table$V3, CMP1="ChA", CMP2="CA", CMP3="pCA",
                                     CMP4="FA", CMP5="GA", CMP6="SA",
                                     CMP7="GeA", CMP8="Ct", CMP9="eCt",
                                     CMP12="R",CMP14="H", CMP18="Q",
                                     CMP19="Phz", CMP20="Pht")
sel.combos.table$V4 <- recode(sel.combos.table$V4, CMP1="ChA", CMP2="CA", CMP3="pCA",
                                     CMP4="FA", CMP5="GA", CMP6="SA",
                                     CMP7="GeA", CMP8="Ct", CMP9="eCt",
                                     CMP12="R",CMP14="H", CMP18="Q",
                                     CMP19="Phz", CMP20="Pht")
sel.combos.table$V5 <- recode(sel.combos.table$V5, CMP1="ChA", CMP2="CA", CMP3="pCA",
                                     CMP4="FA", CMP5="GA", CMP6="SA",
                                     CMP7="GeA", CMP8="Ct", CMP9="eCt",
                                     CMP12="R",CMP14="H", CMP18="Q",
                                     CMP19="Phz", CMP20="Pht")
sel.combos.table$V6 <- recode(sel.combos.table$V6, CMP1="ChA", CMP2="CA", CMP3="pCA",
                                     CMP4="FA", CMP5="GA", CMP6="SA",
                                     CMP7="GeA", CMP8="Ct", CMP9="eCt",
                                     CMP12="R",CMP14="H", CMP18="Q",
                                     CMP19="Phz", CMP20="Pht")
sel.combos.table$V7 <- recode(sel.combos.table$V7, CMP1="ChA", CMP2="CA", CMP3="pCA",
                              CMP4="FA", CMP5="GA", CMP6="SA",
                              CMP7="GeA", CMP8="Ct", CMP9="eCt",
                              CMP12="R",CMP14="H", CMP18="Q",
                              CMP19="Phz", CMP20="Pht")
sel.combos.table$V8 <- recode(sel.combos.table$V8, CMP1="ChA", CMP2="CA", CMP3="pCA",
                              CMP4="FA", CMP5="GA", CMP6="SA",
                              CMP7="GeA", CMP8="Ct", CMP9="eCt",
                              CMP12="R",CMP14="H", CMP18="Q",
                              CMP19="Phz", CMP20="Pht")
sel.combos.table$V9 <- recode(sel.combos.table$V9, CMP1="ChA", CMP2="CA", CMP3="pCA",
                              CMP4="FA", CMP5="GA", CMP6="SA",
                              CMP7="GeA", CMP8="Ct", CMP9="eCt",
                              CMP12="R",CMP14="H", CMP18="Q",
                              CMP19="Phz", CMP20="Pht")
sel.combos.table$V10 <- recode(sel.combos.table$V10, CMP1="ChA", CMP2="CA", CMP3="pCA",
                              CMP4="FA", CMP5="GA", CMP6="SA",
                              CMP7="GeA", CMP8="Ct", CMP9="eCt",
                              CMP12="R",CMP14="H", CMP18="Q",
                              CMP19="Phz", CMP20="Pht")


colnames(sel.combos.table)[2:11] <- c("CMP1", "CMP2", "CMP3", "CMP4", "CMP5", "CMP6",
                                      "CMP7", "CMP8", "CMP9", "CMP10")

#Adding single compounds to the table as treatments
sel.combos.table$CMP1 <- as.character(sel.combos.table$CMP1)

sel.combos.table <- rbind(sel.combos.table, c("ChA", "ChA", rep(NA, 10)))
sel.combos.table <- rbind(sel.combos.table, c("CA", "CA", rep(NA, 10)))
sel.combos.table <- rbind(sel.combos.table, c("pCA", "pCA", rep(NA, 10)))
sel.combos.table <- rbind(sel.combos.table, c("FA", "FA", rep(NA, 10)))
sel.combos.table <- rbind(sel.combos.table, c("GA", "GA", rep(NA, 10)))
sel.combos.table <- rbind(sel.combos.table, c("SA", "SA", rep(NA, 10)))
sel.combos.table <- rbind(sel.combos.table, c("GeA", "GeA", rep(NA, 10)))
sel.combos.table <- rbind(sel.combos.table, c("Ct", "Ct", rep(NA, 10)))
sel.combos.table <- rbind(sel.combos.table, c("eCt", "eCt", rep(NA, 10)))
sel.combos.table <- rbind(sel.combos.table, c("R", "R", rep(NA, 10)))
sel.combos.table <- rbind(sel.combos.table, c("H", "H", rep(NA, 10)))
sel.combos.table <- rbind(sel.combos.table, c("Q", "Q", rep(NA, 10)))
sel.combos.table <- rbind(sel.combos.table, c("Phz", "Phz", rep(NA, 10)))
sel.combos.table <- rbind(sel.combos.table, c("Pht", "Pht", rep(NA, 10)))

#randomizing the order of the mixtures for the experiments
sel.combos.table$rand <- sample(nrow(sel.combos.table))
sel.combos.table <- sel.combos.table[order(sel.combos.table$rand),]

sel.combos.table <- sel.combos.table[-13]

#Saving table. Note that due to the random selection, 
#these tables would be different each time we re-run the code

write.csv(sel.combos.table, file="./Outputs/Tables/CompoundMixtures_Richness2.csv", row.names=FALSE)
  



#################################################################
###For evenness experiments----selecting mixtures
#############################################################

num.comp <- 6  ##Set the number of compounds we want in the mixtures
min.conc <- 1  #Set the minimum percent out of total conc that we can include per compound

###Set parameters for selection of sets
levels <- c(0.2, 0.4, 0.6, 0.8, 1.0) #The levels of evenness to use 
var <- 0.01 #set allowable level of variation in evenness (e.g. 0.1 values can vary from 0.11 to 0.09) 
cut <- .10  #the percentage cutoff for what is considered a high, mid, and low structural diversity group
#candidate sets of 3 will be randomly drawn from this pool, then checked to see if they meet criteria
max.shared <- 4 #set maximum number of compounds that can be shared per group of 6 compounds 
max.occur <- 2 #set number of times an individual compound can occur in the set of three

IDs <- colnames(d.sub)
sel.combos.even <- list()
full.table <- data.frame(FD=factor(), Elevel=factor(),Comps=factor(), Percent=numeric()) #list of all compounds and amounts for calculation of needs
for (i in 1:length(levels)){  #For each evenness level
  nr <- max.occur
  g <- max.shared 
  d <- t(combn(IDs, num.comp))    #create a list of all possible compound combinations
  d <- cbind(d, NA)
  
    for (j in 1:nrow(d)){     #For each compound combo
      #Get the average similarity index across all pw contrasts in combo
      list <- as.character(d[j,1:num.comp])   
      pw.combos <- combn(list, 2)
      pw.combos <- rbind(pw.combos, NA)
      for (k in 1:ncol(pw.combos)){
        x <- d.sub[as.character(pw.combos[1,k]), as.character(pw.combos[2,k])]  
        pw.combos[3,k] <- x
      }
      z <- mean(as.numeric(pw.combos[3,]))
      d[j,num.comp+1] <- z
    }
  
  d <- as.data.frame(d) 
  names(d)[num.comp+1] <- "MeanDist" 
  d2 <- d[order(d$MeanDist),]
  #d2 is a list of all possible combinations ordered by their mean Similarity index for the current richness level
  
  #Separate out groups for high, mid, and low functional diversity based on specified cutoff
  high <- d2[1:round((nrow(d2)*cut)),]
  
  #for mid group
  up <- round(nrow(d2)/2 + nrow(high)/2)
  down <- round(nrow(d2)/2 - nrow(high)/2)
  mid <- d2[up:down,]
  
  low <- d2[(nrow(d2)-nrow(high)):nrow(d2),]
  
  #Take a random sample of three sets from each group
  high.set <- high[sample(nrow(high), 3, replace=FALSE),]
  mid.set <- mid[sample(nrow(mid), 3, replace=FALSE),]
  low.set <- low[sample(nrow(low), 3, replace=FALSE),] 
  
  
  ###check whether the random samples of three groups meet the specified criteria
  
  #First for high
  high.set.final <- 0
  while(high.set.final==0){
    
    t <- c(as.character(unname(unlist(high.set[1,1:num.comp]))),
           as.character(unname(unlist(high.set[2,1:num.comp]))),
           as.character(unname(unlist(high.set[3,1:num.comp]))))
    
    #make a list of overlaps between each pairwise comparison
    one.two <- intersect(as.character(unname(unlist(high.set[1,1:num.comp]))),
                         as.character(unname(unlist(high.set[2,1:num.comp]))))
    two.three <- intersect(as.character(unname(unlist(high.set[2,1:num.comp]))),
                           as.character(unname(unlist(high.set[3,1:num.comp]))))
    one.three <- intersect(as.character(unname(unlist(high.set[1,1:num.comp]))),
                           as.character(unname(unlist(high.set[3,1:num.comp]))))
    
      ##check if combos meet criteria for g and nr
      if(all(table(t)<=nr) &&
        length(one.two) <= g  &&
        length(one.three) <= g  &&
        length(two.three) <= g  ){ 
        high.set.final <- high.set  
        } else {
        high.set <- high[sample(nrow(high), 3, replace=FALSE),]
      }
    }
  
  #For mid
  mid.set.final <- 0
  while(mid.set.final==0){
    
    t <- c(as.character(unname(unlist(mid.set[1,1:num.comp]))),
           as.character(unname(unlist(mid.set[2,1:num.comp]))),
           as.character(unname(unlist(mid.set[3,1:num.comp]))))
    
    #make a list of overlaps between each pairwise comparison
    one.two <- intersect(as.character(unname(unlist(mid.set[1,1:num.comp]))),
                         as.character(unname(unlist(mid.set[2,1:num.comp]))))
    two.three <- intersect(as.character(unname(unlist(mid.set[2,1:num.comp]))),
                           as.character(unname(unlist(mid.set[3,1:num.comp]))))
    one.three <- intersect(as.character(unname(unlist(mid.set[1,1:num.comp]))),
                           as.character(unname(unlist(mid.set[3,1:num.comp]))))
    
    ##check if combos meet criteria for g and nr
      if(all(table(t)<=nr) &&
        length(one.two) <= g  &&
        length(one.three) <= g  &&
        length(two.three) <= g  ){ 
        mid.set.final <- mid.set  
        } else {
        mid.set <- mid[sample(nrow(mid), 3, replace=FALSE),]
      }
    }  
  
  #For low
  low.set.final <- 0
  while(low.set.final==0){
    
    t <- c(as.character(unname(unlist(low.set[1,1:num.comp]))),
           as.character(unname(unlist(low.set[2,1:num.comp]))),
           as.character(unname(unlist(low.set[3,1:num.comp]))))
    
    #make a list of overlaps between each pairwise comparison
    one.two <- intersect(as.character(unname(unlist(low.set[1,1:num.comp]))),
                         as.character(unname(unlist(low.set[2,1:num.comp]))))
    two.three <- intersect(as.character(unname(unlist(low.set[2,1:num.comp]))),
                           as.character(unname(unlist(low.set[3,1:num.comp]))))
    one.three <- intersect(as.character(unname(unlist(low.set[1,1:num.comp]))),
                           as.character(unname(unlist(low.set[3,1:num.comp]))))
    
    ##check if combos meet criteria for g and nr
      if(all(table(t)<=nr) &&
         length(one.two) <= g  &&
         length(one.three) <= g  &&
         length(two.three) <= g  ){ 
        low.set.final <- low.set  
        } else {
        low.set <- low[sample(nrow(low), 3, replace=FALSE),]
      }
    }

  
  #Now randomly assign evenness values to each set
  done <- 0
  sel.mix <- matrix(ncol=num.comp+1)
  Es <- vector()
    while(done==0){
      if(levels[i] != 1.0){
        #Create a vector of potential concentrations
        y <- vector(mode="numeric", length=num.comp)
          for (j in 1:(num.comp-1)){
            y[j] <- sample(min.conc:(100-sum(y)-((num.comp-j)*min.conc)), 1)
          }
        y[num.comp] <- 100-sum(y)
        y <- sample(y) #Randomly mix vector to make sure first value isn't always biggest
        } else {  #For when evenness == 1
          y <- rep(100/num.comp, num.comp)
        }
      #Check what the evenness value would be for those conc values
      pi2 <- (y/sum(y))^2
      D <- 1/sum(pi2)
      E <- D/num.comp  #Evenness
      Es <- c(Es, E)  ##Put all evenness values in a vector so we can see distribution
        if(E <= levels[i]+var && E >= levels[i]-var){
          y <- c(y, E)
          sel.mix <- rbind(sel.mix, y)
          colnames(sel.mix)[(num.comp+1)] <- "Evenness"
        }
      if (nrow(sel.mix)==10){done=1}
    }
  
  ##Create final list of compound mixes and ratios
  names.sel <- c(names(sel.combos.even), paste("high", levels[i], sep=""),paste("mid", levels[i], sep=""),
                 paste("low", levels[i], sep=""))    
  sel.combos.even <- c(sel.combos.even, list(abind(data.frame(high.set.final), data.frame(sel.mix[2:4,]), along=2),
                                             abind(data.frame(mid.set.final), data.frame(sel.mix[5:7,]), along=2),
                                             abind(data.frame(low.set.final), data.frame(sel.mix[8:10,]), along=2)))
  names(sel.combos.even) <- names.sel

    for (j in 1:num.comp){
      full.table <- rbind(full.table, cbind(rep("high", 3),rep(levels[i], 3), as.character(high.set.final[,j]), as.numeric(sel.mix[2:4,j])))
      full.table <- rbind(full.table, cbind(rep("mid", 3),rep(levels[i], 3), as.character(mid.set.final[,j]), as.numeric(sel.mix[5:7,j])))
      full.table <- rbind(full.table, cbind(rep("low", 3),rep(levels[i], 3), as.character(low.set.final[,j]), as.numeric(sel.mix[8:10,j])))
    }
}



sel.combos.even
names(sel.combos.even)

sel.combos.even.table <- data.frame(Tx=paste(rep(names(sel.combos.even), each=3), rep(c("a","b","c"), 6), sep=""), 
                                    CMP1=c(sel.combos.even$high0.2[,1], sel.combos.even$mid0.2[,1], sel.combos.even$low0.2[,1],
                                         sel.combos.even$high0.4[,1], sel.combos.even$mid0.4[,1], sel.combos.even$low0.4[,1],
                                         sel.combos.even$high0.6[,1], sel.combos.even$mid0.6[,1], sel.combos.even$low0.6[,1],
                                         sel.combos.even$high0.8[,1], sel.combos.even$mid0.8[,1], sel.combos.even$low0.8[,1],
                                         sel.combos.even$high1[,1], sel.combos.even$mid1[,1], sel.combos.even$low1[,1]),
                                    prop_CMP1=round(as.numeric(c(sel.combos.even$high0.2[,8], sel.combos.even$mid0.2[,8], sel.combos.even$low0.2[,8],
                                         sel.combos.even$high0.4[,8], sel.combos.even$mid0.4[,8], sel.combos.even$low0.4[,8],
                                         sel.combos.even$high0.6[,8], sel.combos.even$mid0.6[,8], sel.combos.even$low0.6[,8],
                                         sel.combos.even$high0.8[,8], sel.combos.even$mid0.8[,8], sel.combos.even$low0.8[,8],
                                         sel.combos.even$high1[,8], sel.combos.even$mid1[,8], sel.combos.even$low1[,8]))/100, digits=3),
                                    CMP2=c(sel.combos.even$high0.2[,2], sel.combos.even$mid0.2[,2], sel.combos.even$low0.2[,2],
                                         sel.combos.even$high0.4[,2], sel.combos.even$mid0.4[,2], sel.combos.even$low0.4[,2],
                                         sel.combos.even$high0.6[,2], sel.combos.even$mid0.6[,2], sel.combos.even$low0.6[,2],
                                         sel.combos.even$high0.8[,2], sel.combos.even$mid0.8[,2], sel.combos.even$low0.8[,2],
                                         sel.combos.even$high1[,2], sel.combos.even$mid1[,2], sel.combos.even$low1[,2]),
                                    prop_CMP2=round(as.numeric(c(sel.combos.even$high0.2[,9], sel.combos.even$mid0.2[,9], sel.combos.even$low0.2[,9],
                                         sel.combos.even$high0.4[,9], sel.combos.even$mid0.4[,9], sel.combos.even$low0.4[,9],
                                         sel.combos.even$high0.6[,9], sel.combos.even$mid0.6[,9], sel.combos.even$low0.6[,9],
                                         sel.combos.even$high0.8[,9], sel.combos.even$mid0.8[,9], sel.combos.even$low0.8[,9],
                                         sel.combos.even$high1[,9], sel.combos.even$mid1[,9], sel.combos.even$low1[,9]))/100, digits=3),
                                    CMP3=c(sel.combos.even$high0.2[,3], sel.combos.even$mid0.2[,3], sel.combos.even$low0.2[,3],
                                         sel.combos.even$high0.4[,3], sel.combos.even$mid0.4[,3], sel.combos.even$low0.4[,3],
                                         sel.combos.even$high0.6[,3], sel.combos.even$mid0.6[,3], sel.combos.even$low0.6[,3],
                                         sel.combos.even$high0.8[,3], sel.combos.even$mid0.8[,3], sel.combos.even$low0.8[,3],
                                         sel.combos.even$high1[,3], sel.combos.even$mid1[,3], sel.combos.even$low1[,3]),
                                    prop_CMP3=round(as.numeric(c(sel.combos.even$high0.2[,10], sel.combos.even$mid0.2[,10], sel.combos.even$low0.2[,10],
                                         sel.combos.even$high0.4[,10], sel.combos.even$mid0.4[,10], sel.combos.even$low0.4[,10],
                                         sel.combos.even$high0.6[,10], sel.combos.even$mid0.6[,10], sel.combos.even$low0.6[,10],
                                         sel.combos.even$high0.8[,10], sel.combos.even$mid0.8[,10], sel.combos.even$low0.8[,10],
                                         sel.combos.even$high1[,10], sel.combos.even$mid1[,10], sel.combos.even$low1[,10]))/100, digits=3),
                                    CMP4=c(sel.combos.even$high0.2[,4], sel.combos.even$mid0.2[,4], sel.combos.even$low0.2[,4],
                                         sel.combos.even$high0.4[,4], sel.combos.even$mid0.4[,4], sel.combos.even$low0.4[,4],
                                         sel.combos.even$high0.6[,4], sel.combos.even$mid0.6[,4], sel.combos.even$low0.6[,4],
                                         sel.combos.even$high0.8[,4], sel.combos.even$mid0.8[,4], sel.combos.even$low0.8[,4],
                                         sel.combos.even$high1[,4], sel.combos.even$mid1[,4], sel.combos.even$low1[,4]),
                                    prop_CMP4=round(as.numeric(c(sel.combos.even$high0.2[,11], sel.combos.even$mid0.2[,11], sel.combos.even$low0.2[,11],
                                         sel.combos.even$high0.4[,11], sel.combos.even$mid0.4[,11], sel.combos.even$low0.4[,11],
                                         sel.combos.even$high0.6[,11], sel.combos.even$mid0.6[,11], sel.combos.even$low0.6[,11],
                                         sel.combos.even$high0.8[,11], sel.combos.even$mid0.8[,11], sel.combos.even$low0.8[,11],
                                         sel.combos.even$high1[,11], sel.combos.even$mid1[,11], sel.combos.even$low1[,11]))/100, digits=3),
                                    CMP5=c(sel.combos.even$high0.2[,5], sel.combos.even$mid0.2[,5], sel.combos.even$low0.2[,5],
                                         sel.combos.even$high0.4[,5], sel.combos.even$mid0.4[,5], sel.combos.even$low0.4[,5],
                                         sel.combos.even$high0.6[,5], sel.combos.even$mid0.6[,5], sel.combos.even$low0.6[,5],
                                         sel.combos.even$high0.8[,5], sel.combos.even$mid0.8[,5], sel.combos.even$low0.8[,5],
                                         sel.combos.even$high1[,5], sel.combos.even$mid1[,5], sel.combos.even$low1[,5]),
                                    prop_CMP5=round(as.numeric(c(sel.combos.even$high0.2[,12], sel.combos.even$mid0.2[,12], sel.combos.even$low0.2[,12],
                                         sel.combos.even$high0.4[,12], sel.combos.even$mid0.4[,12], sel.combos.even$low0.4[,12],
                                         sel.combos.even$high0.6[,12], sel.combos.even$mid0.6[,12], sel.combos.even$low0.6[,12],
                                         sel.combos.even$high0.8[,12], sel.combos.even$mid0.8[,12], sel.combos.even$low0.8[,12],
                                         sel.combos.even$high1[,12], sel.combos.even$mid1[,12], sel.combos.even$low1[,12]))/100, digits=3),
                                    CMP6=c(sel.combos.even$high0.2[,6], sel.combos.even$mid0.2[,6], sel.combos.even$low0.2[,6],
                                         sel.combos.even$high0.4[,6], sel.combos.even$mid0.4[,6], sel.combos.even$low0.4[,6],
                                         sel.combos.even$high0.6[,6], sel.combos.even$mid0.6[,6], sel.combos.even$low0.6[,6],
                                         sel.combos.even$high0.8[,6], sel.combos.even$mid0.8[,6], sel.combos.even$low0.8[,6],
                                         sel.combos.even$high1[,6], sel.combos.even$mid1[,6], sel.combos.even$low1[,6]),
                                    prop_CMP6=round(as.numeric(c(sel.combos.even$high0.2[,13], sel.combos.even$mid0.2[,13], sel.combos.even$low0.2[,13],
                                         sel.combos.even$high0.4[,13], sel.combos.even$mid0.4[,13], sel.combos.even$low0.4[,13],
                                         sel.combos.even$high0.6[,13], sel.combos.even$mid0.6[,13], sel.combos.even$low0.6[,13],
                                         sel.combos.even$high0.8[,13], sel.combos.even$mid0.8[,13], sel.combos.even$low0.8[,13],
                                         sel.combos.even$high1[,13], sel.combos.even$mid1[,13], sel.combos.even$low1[,13]))/100, digits=3),
                                    MeanDist=round(as.numeric(c(sel.combos.even$high0.2[,7], sel.combos.even$mid0.2[,7], sel.combos.even$low0.2[,7],
                                         sel.combos.even$high0.4[,7], sel.combos.even$mid0.4[,7], sel.combos.even$low0.4[,7],
                                         sel.combos.even$high0.6[,7], sel.combos.even$mid0.6[,7], sel.combos.even$low0.6[,7],
                                         sel.combos.even$high0.8[,7], sel.combos.even$mid0.8[,7], sel.combos.even$low0.8[,7],
                                         sel.combos.even$high1[,7], sel.combos.even$mid1[,7], sel.combos.even$low1[,7])), digits=3),
                                    Evenness=round(as.numeric(c(sel.combos.even$high0.2[,14], sel.combos.even$mid0.2[,14], sel.combos.even$low0.2[,14],
                                         sel.combos.even$high0.4[,14], sel.combos.even$mid0.4[,14], sel.combos.even$low0.4[,14],
                                         sel.combos.even$high0.6[,14], sel.combos.even$mid0.6[,14], sel.combos.even$low0.6[,14],
                                         sel.combos.even$high0.8[,14], sel.combos.even$mid0.8[,14], sel.combos.even$low0.8[,14],
                                         sel.combos.even$high1[,14], sel.combos.even$mid1[,14], sel.combos.even$low1[,14])), digits=3))


#rename compounds with meaningful names
sel.combos.even.table$CMP1 <- recode(sel.combos.even.table$CMP1, CMP1="ChA", CMP2="CA", CMP3="pCA",
                                     CMP4="FA", CMP5="GA", CMP6="SA",
                                     CMP7="GeA", CMP8="Ct", CMP9="eCt",
                                     CMP12="R",CMP14="H", CMP18="Q",
                                     CMP19="Phz", CMP20="Pht")
sel.combos.even.table$CMP2 <- recode(sel.combos.even.table$CMP2, CMP1="ChA", CMP2="CA", CMP3="pCA",
                                     CMP4="FA", CMP5="GA", CMP6="SA",
                                     CMP7="GeA", CMP8="Ct", CMP9="eCt",
                                     CMP12="R",CMP14="H", CMP18="Q",
                                     CMP19="Phz", CMP20="Pht")
sel.combos.even.table$CMP3 <- recode(sel.combos.even.table$CMP3, CMP1="ChA", CMP2="CA", CMP3="pCA",
                                     CMP4="FA", CMP5="GA", CMP6="SA",
                                     CMP7="GeA", CMP8="Ct", CMP9="eCt",
                                     CMP12="R",CMP14="H", CMP18="Q",
                                     CMP19="Phz", CMP20="Pht")
sel.combos.even.table$CMP4 <- recode(sel.combos.even.table$CMP4, CMP1="ChA", CMP2="CA", CMP3="pCA",
                                     CMP4="FA", CMP5="GA", CMP6="SA",
                                     CMP7="GeA", CMP8="Ct", CMP9="eCt",
                                     CMP12="R",CMP14="H", CMP18="Q",
                                     CMP19="Phz", CMP20="Pht")
sel.combos.even.table$CMP5 <- recode(sel.combos.even.table$CMP5, CMP1="ChA", CMP2="CA", CMP3="pCA",
                                     CMP4="FA", CMP5="GA", CMP6="SA",
                                     CMP7="GeA", CMP8="Ct", CMP9="eCt",
                                     CMP12="R",CMP14="H", CMP18="Q",
                                     CMP19="Phz", CMP20="Pht")
sel.combos.even.table$CMP6 <- recode(sel.combos.even.table$CMP6, CMP1="ChA", CMP2="CA", CMP3="pCA",
                                     CMP4="FA", CMP5="GA", CMP6="SA",
                                     CMP7="GeA", CMP8="Ct", CMP9="eCt",
                                     CMP12="R",CMP14="H", CMP18="Q",
                                     CMP19="Phz", CMP20="Pht")



#randomize order for experiments
sel.combos.even.table$rand <- sample(nrow(sel.combos.even.table))
sel.combos.even.table <- sel.combos.even.table[order(sel.combos.even.table$rand),]

sel.combos.even.table <- sel.combos.even.table[-16]

#Saving table. Note that due to the random selection, 
#these tables would be different each time we re-run the code

write.csv(sel.combos.even.table, "./Outputs/Tables/CompoundMixtures_Evenness2.csv", row.names=FALSE)




#######################################################################
##Calculate amounts of compounds needed for experiments
###############################################################

#Set parameters

#Enter total concentration to be used in mG/100G
conc <- 200  #200 mg/100g is average total phenolic conc reported from apples on Phenol Explorer
diet.vol <- 5 #enter volume of diet in mL to be added to each cup
cups.per.Tx <- 22 #20 for N=5 reps per species, could also be 22--3 each for Cp and Px (so N=9), and 8 each for Sf and Hz  
cups.per.single <- 44 #60 for N=15 per species, could be 44--6 for each Cp and Px (so N=18), and 16 each for Sf and Hz   


diet.weight <- diet.vol*1092/1000
cmp.per.cup <- 200*diet.weight/100 #gives total mg of compound per cup




#For richness
needs.rich <- data.frame(IDs=IDs, NeedsPerRep=0)

for (i in 1:length(IDs)){   #For each compound
  for (j in 1:length(sel.combos)){  #For each matrix
    for (k in 1:nrow(sel.combos[[j]])){  #For each row
      for(l in 1:(ncol(sel.combos[[j]])-1)){  #For each element in each row
        if (sel.combos[[j]][k,l]==IDs[i]){
          to.add <- cmp.per.cup/(ncol(sel.combos[[j]])-1)
          needs.rich[i,2] <- needs.rich[i,2] + to.add    
        }
      }
    } 
  }
}

needs.rich$TotalNeeds.Mix <- needs.rich$NeedsPerRep*cups.per.Tx
needs.rich$TotalNeeds <- needs.rich$TotalNeeds.Mix + cups.per.single*cmp.per.cup + cups.per.Tx*cmp.per.cup/length(IDs)

needs.rich
sum(needs.rich$TotalNeeds)


##For evenness

needs.even <- data.frame(IDs=IDs, NeedsPerRep=0)

for (i in 1:length(IDs)){   #For each compound
  for (j in 1:length(sel.combos.even)){  #For each matrix
    for (k in 1:nrow(sel.combos.even[[j]])){  #For each row
      for(l in 1:num.comp){  #For each compound in each row
        if (sel.combos.even[[j]][k,l]==IDs[i]){
          to.add <- as.numeric(sel.combos.even[[j]][k,(l+num.comp+1)])*cmp.per.cup/100
          needs.even[i,2] <- needs.even[i,2] + to.add    
        }
      }
    } 
  }
}

needs.even$TotalNeeds <- needs.even$NeedsPerRep*cups.per.Tx


needs.even
sum(needs.rich$TotalNeeds)


needs.total <- cbind(as.character(needs.rich$IDs), (needs.even$TotalNeeds + needs.rich$TotalNeeds))








