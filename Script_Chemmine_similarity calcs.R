#load("ChemDiversity.RData")  #saved outputs
#load("ChemDiversity_V2.RData")  #saved outputs


#ChemmineR and fmcsR are bioconductor packages. They need to be sourced
#from there--they are not available on the regular CRAN mirrors
# See https://bioconductor.org/install/#install-bioconductor-packages

#to first install, run...
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c("ChemmineR", "fmcsR"))


library(ChemmineR) # Loads the package
library(fmcsR) #add on package for Maximum Common Substructure (MCS) searching
library(abind)  #For nice pasting together of matrices
vignette("ChemmineR") #Opens this PDF manual from R 

#The compound structure data is stored as an SDF file from PubChem
#Go to: https://pubchem.ncbi.nlm.nih.gov/
#upload a list of identifiers for the compounds (e.g. CIDs) as csv file 
#(I used the list in "Whitehead_et_al_CIDs" but took out header and names column)
#then choose "push to entrez" then you can download the structure data as an sdf 
#file with gz compression
sdfset <- read.SDFset("Whitehead_et_al_structures.sdf.gz") 
names <- read.csv("Whitehead_et_al_CIDs.csv")


#Prune compound set to things we can get
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
####For richness experiments-----------
#########################################

###Set parameters for selection of sets
levels <- c(2,4,6,8,10) #The levels of richness to use; we would also have 1 and 14 but there is only one way to make those
cut <- .10  #the percentage cutoff for what is considered a high, mid, and low functional diversity group
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



save.image("ChemDiversity.RData") 


write.csv(sel.combos, file="Sel_combos.csv")
  



#################################################################
###For evenness experiments
#############################################################

#conc <- 200  ##set the total concentration that we want to use in experiments
num.comp <- 6  ##Set the number of compounds we want in the mixtures
min.conc <- 1  #Set the minimum percent out of total conc that we can include per compound

###Set parameters for selection of sets
levels <- c(0.2, 0.4, 0.6, 0.8, 1.0) #The levels of evenness to use 
var <- 0.01 #set allowable level of variation in evenness (e.g. 0.1 values can vary from 0.11 to 0.09) 
cut <- .10  #the percentage cutoff for what is considered a high, mid, and low functional diversity group
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
                                    V1=c(sel.combos.even$high0.2[,1], sel.combos.even$mid0.2[,1], sel.combos.even$low0.2[,1],
                                         sel.combos.even$high0.4[,1], sel.combos.even$mid0.4[,1], sel.combos.even$low0.4[,1],
                                         sel.combos.even$high0.6[,1], sel.combos.even$mid0.6[,1], sel.combos.even$low0.6[,1],
                                         sel.combos.even$high0.8[,1], sel.combos.even$mid0.8[,1], sel.combos.even$low0.8[,1],
                                         sel.combos.even$high1[,1], sel.combos.even$mid1[,1], sel.combos.even$low1[,1]),
                                    E1=round(as.numeric(c(sel.combos.even$high0.2[,8], sel.combos.even$mid0.2[,8], sel.combos.even$low0.2[,8],
                                         sel.combos.even$high0.4[,8], sel.combos.even$mid0.4[,8], sel.combos.even$low0.4[,8],
                                         sel.combos.even$high0.6[,8], sel.combos.even$mid0.6[,8], sel.combos.even$low0.6[,8],
                                         sel.combos.even$high0.8[,8], sel.combos.even$mid0.8[,8], sel.combos.even$low0.8[,8],
                                         sel.combos.even$high1[,8], sel.combos.even$mid1[,8], sel.combos.even$low1[,8]))/100, digits=3),
                                    V2=c(sel.combos.even$high0.2[,2], sel.combos.even$mid0.2[,2], sel.combos.even$low0.2[,2],
                                         sel.combos.even$high0.4[,2], sel.combos.even$mid0.4[,2], sel.combos.even$low0.4[,2],
                                         sel.combos.even$high0.6[,2], sel.combos.even$mid0.6[,2], sel.combos.even$low0.6[,2],
                                         sel.combos.even$high0.8[,2], sel.combos.even$mid0.8[,2], sel.combos.even$low0.8[,2],
                                         sel.combos.even$high1[,2], sel.combos.even$mid1[,2], sel.combos.even$low1[,2]),
                                    E2=round(as.numeric(c(sel.combos.even$high0.2[,9], sel.combos.even$mid0.2[,9], sel.combos.even$low0.2[,9],
                                         sel.combos.even$high0.4[,9], sel.combos.even$mid0.4[,9], sel.combos.even$low0.4[,9],
                                         sel.combos.even$high0.6[,9], sel.combos.even$mid0.6[,9], sel.combos.even$low0.6[,9],
                                         sel.combos.even$high0.8[,9], sel.combos.even$mid0.8[,9], sel.combos.even$low0.8[,9],
                                         sel.combos.even$high1[,9], sel.combos.even$mid1[,9], sel.combos.even$low1[,9]))/100, digits=3),
                                    V3=c(sel.combos.even$high0.2[,3], sel.combos.even$mid0.2[,3], sel.combos.even$low0.2[,3],
                                         sel.combos.even$high0.4[,3], sel.combos.even$mid0.4[,3], sel.combos.even$low0.4[,3],
                                         sel.combos.even$high0.6[,3], sel.combos.even$mid0.6[,3], sel.combos.even$low0.6[,3],
                                         sel.combos.even$high0.8[,3], sel.combos.even$mid0.8[,3], sel.combos.even$low0.8[,3],
                                         sel.combos.even$high1[,3], sel.combos.even$mid1[,3], sel.combos.even$low1[,3]),
                                    E3=round(as.numeric(c(sel.combos.even$high0.2[,10], sel.combos.even$mid0.2[,10], sel.combos.even$low0.2[,10],
                                         sel.combos.even$high0.4[,10], sel.combos.even$mid0.4[,10], sel.combos.even$low0.4[,10],
                                         sel.combos.even$high0.6[,10], sel.combos.even$mid0.6[,10], sel.combos.even$low0.6[,10],
                                         sel.combos.even$high0.8[,10], sel.combos.even$mid0.8[,10], sel.combos.even$low0.8[,10],
                                         sel.combos.even$high1[,10], sel.combos.even$mid1[,10], sel.combos.even$low1[,10]))/100, digits=3),
                                    V4=c(sel.combos.even$high0.2[,4], sel.combos.even$mid0.2[,4], sel.combos.even$low0.2[,4],
                                         sel.combos.even$high0.4[,4], sel.combos.even$mid0.4[,4], sel.combos.even$low0.4[,4],
                                         sel.combos.even$high0.6[,4], sel.combos.even$mid0.6[,4], sel.combos.even$low0.6[,4],
                                         sel.combos.even$high0.8[,4], sel.combos.even$mid0.8[,4], sel.combos.even$low0.8[,4],
                                         sel.combos.even$high1[,4], sel.combos.even$mid1[,4], sel.combos.even$low1[,4]),
                                    E4=round(as.numeric(c(sel.combos.even$high0.2[,11], sel.combos.even$mid0.2[,11], sel.combos.even$low0.2[,11],
                                         sel.combos.even$high0.4[,11], sel.combos.even$mid0.4[,11], sel.combos.even$low0.4[,11],
                                         sel.combos.even$high0.6[,11], sel.combos.even$mid0.6[,11], sel.combos.even$low0.6[,11],
                                         sel.combos.even$high0.8[,11], sel.combos.even$mid0.8[,11], sel.combos.even$low0.8[,11],
                                         sel.combos.even$high1[,11], sel.combos.even$mid1[,11], sel.combos.even$low1[,11]))/100, digits=3),
                                    V5=c(sel.combos.even$high0.2[,5], sel.combos.even$mid0.2[,5], sel.combos.even$low0.2[,5],
                                         sel.combos.even$high0.4[,5], sel.combos.even$mid0.4[,5], sel.combos.even$low0.4[,5],
                                         sel.combos.even$high0.6[,5], sel.combos.even$mid0.6[,5], sel.combos.even$low0.6[,5],
                                         sel.combos.even$high0.8[,5], sel.combos.even$mid0.8[,5], sel.combos.even$low0.8[,5],
                                         sel.combos.even$high1[,5], sel.combos.even$mid1[,5], sel.combos.even$low1[,5]),
                                    E5=round(as.numeric(c(sel.combos.even$high0.2[,12], sel.combos.even$mid0.2[,12], sel.combos.even$low0.2[,12],
                                         sel.combos.even$high0.4[,12], sel.combos.even$mid0.4[,12], sel.combos.even$low0.4[,12],
                                         sel.combos.even$high0.6[,12], sel.combos.even$mid0.6[,12], sel.combos.even$low0.6[,12],
                                         sel.combos.even$high0.8[,12], sel.combos.even$mid0.8[,12], sel.combos.even$low0.8[,12],
                                         sel.combos.even$high1[,12], sel.combos.even$mid1[,12], sel.combos.even$low1[,12]))/100, digits=3),
                                    V6=c(sel.combos.even$high0.2[,6], sel.combos.even$mid0.2[,6], sel.combos.even$low0.2[,6],
                                         sel.combos.even$high0.4[,6], sel.combos.even$mid0.4[,6], sel.combos.even$low0.4[,6],
                                         sel.combos.even$high0.6[,6], sel.combos.even$mid0.6[,6], sel.combos.even$low0.6[,6],
                                         sel.combos.even$high0.8[,6], sel.combos.even$mid0.8[,6], sel.combos.even$low0.8[,6],
                                         sel.combos.even$high1[,6], sel.combos.even$mid1[,6], sel.combos.even$low1[,6]),
                                    E6=round(as.numeric(c(sel.combos.even$high0.2[,13], sel.combos.even$mid0.2[,13], sel.combos.even$low0.2[,13],
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


sel.combos.even.table$rand <- sample(nrow(sel.combos.even.table))
sel.combos.even.table <- sel.combos.even.table[order(sel.combos.even.table$rand),]



write.csv(sel.combos.even.table, file="Sel_combos_even.csv")

save.image("ChemDiversity_V2.RData") 

full.table_0.2 <- data.frame(full.table)
#Es_0.2_ShanWein <- Es  #stored objects to compare distribution of Es with SW index vs Simpson
#Es_0.2_Simpson <- Es




##############################
##Richness vs Evenness
#############################################

###Set parameters for selection of sets
Rlevels <- c(2,4,6,8,10) #The levels of richness to use; we would also have 1 and 14 but there is only one way to make those
Elevels <- c(0.2, 0.4, 0.6, 0.8, 1)
cut <- .10  #the percentage cutoff for what is considered a high, mid, and low functional diversity group
#candidate sets of 3 will be randomly drawn from this pool, then checked to see if they meet criteria
max.shared <- c(0,2,4,5,7) #set maximum number of compounds that can be shared per group for each richness level
#must be the same length as levels
##So at richness=4, only two compounds can be shared per group
max.occur <- c(1,2,2,2,3) #set number of times an individual compound can occur in a set of three for each richness level

min.conc <- 1  #Set the minimum percent out of total conc that we can include per compound
var <- 0.01 #set allowable level of variation in evenness (e.g. 0.1 values can vary from 0.11 to 0.09) 


IDs <- colnames(d.sub)
sel.combos.RvE <- list()
for (i in 1:length(Rlevels)){  #For each richness level
  nr <- max.occur[i]
  g <- max.shared[i] 
  d <- t(combn(IDs, Rlevels[i]))    #create a list of all possible compound combinations
  d <- cbind(d, NA)
  
  for (j in 1:nrow(d)){     #For each compound combo
    #Get the average similarity index across all pw contrasts in combo
    list <- as.character(d[j,1:Rlevels[i]])   
    pw.combos <- combn(list, 2)
    pw.combos <- rbind(pw.combos, NA)
    for (k in 1:ncol(pw.combos)){
      x <- d.sub[as.character(pw.combos[1,k]), as.character(pw.combos[2,k])]  
      pw.combos[3,k] <- x
    }
    z <- mean(as.numeric(pw.combos[3,]))
    d[j,Rlevels[i]+1] <- z
  }
  
  d <- as.data.frame(d) 
  names(d)[Rlevels[i]+1] <- "MeanDist" 
  d2 <- d[order(d$MeanDist),]
  #d2 is a list of all possible combinations ordered by their mean Similarity index for the current richness level
  
  #Separate out groups for high, mid, and low functional diversity based on specified cutoff
  high <- d2[1:round((nrow(d2)*cut)),]
  
  #for mid group
  up <- round(nrow(d2)/2 + nrow(high)/2)
  down <- round(nrow(d2)/2 - nrow(high)/2)
  mid <- d2[up:down,]
  
  low <- d2[(nrow(d2)-nrow(high)):nrow(d2),]
  
  for (j in 1:length(Elevels)){
    set <- mid[sample(nrow(mid), 3, replace=FALSE),]
    set.final <- 0
    while(set.final==0){  ##Checking if set of 3 meets criteria
      
      t <- c(as.character(unname(unlist(set[1,1:levels[i]]))),
             as.character(unname(unlist(set[2,1:levels[i]]))),
             as.character(unname(unlist(set[3,1:levels[i]]))))
      
      #make a list of overlaps between each pairwise comparison
      one.two <- intersect(as.character(unname(unlist(set[1,1:levels[i]]))),
                           as.character(unname(unlist(set[2,1:levels[i]]))))
      two.three <- intersect(as.character(unname(unlist(set[2,1:levels[i]]))),
                             as.character(unname(unlist(set[3,1:levels[i]]))))
      one.three <- intersect(as.character(unname(unlist(set[1,1:levels[i]]))),
                             as.character(unname(unlist(set[3,1:levels[i]]))))
      
      ##check if combos meet criteria for g and nr
      if(all(table(t)<=nr) &&
         length(one.two) <= g  &&
         length(one.three) <= g  &&
         length(two.three) <= g  ){ 
        set.final <- set  
        } else {
        set <- mid[sample(nrow(high), 3, replace=FALSE),]
      }
    }
    
    #Now randomly assign evenness values to final set
    done <- 0
    sel.mix <- matrix(ncol=Rlevels[i]+1)
    Es <- vector()
    while(done==0){
      if(Elevels[j] != 1.0){
        #Create a vector of potential concentrations
        y <- vector(mode="numeric", length=Rlevels[i])
        for (k in 1:(Rlevels[i]-1)){
          y[k] <- sample(min.conc:(100-sum(y)-((Rlevels[i]-k)*min.conc)), 1)
        }
        y[Rlevels[i]] <- 100-sum(y)
        y <- sample(y) #Randomly mix vector to make sure first value isn't always biggest
      } else {  #For when evenness == 1
        y <- rep(100/Rlevels[i], Rlevels[i])
      }
      #Check what the evenness value would be for those conc values
      pi2 <- (y/sum(y))^2
      D <- 1/sum(pi2)
      E <- D/Rlevels[i]  #Evenness
      Es <- c(Es, E)  ##Put all evenness values in a vector so we can see distribution
      if(E <= Elevels[i]+var && E >= Elevels[i]-var){
        y <- c(y, E)
        sel.mix <- rbind(sel.mix, y)
        colnames(sel.mix)[(Rlevels[i]+1)] <- "Evenness"
      }
      if (nrow(sel.mix)==10){done=1}
    }
    
    ##Create final list of compound mixes and ratios
    names.sel <- c(names(sel.combos.RvE), paste("Rich",Rlevels[i], ";Even", Elevels[j], sep=""))    
    sel.combos.RvE <- c(sel.combos.RvE, list(abind(data.frame(set.final), data.frame(sel.mix[2:4,]), along=2)))
    names(sel.combos.RvE) <- names.sel
  }
}  
  


save.image("ChemDiversity_V3.RData") 


#######################################################################
##Calculate amounts of compounds needed for experiments
###############################################################

#Set parameters

#Enter total concentration to be used in mG/100G
conc <- 200  #200 mg/100g is average total phenolc conc reported from apples on Phenol Explorer
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








######################------------------------------------------------------------------------------------
####OLD CODE

full.table <- data.frame(full.table)
needs.even <- aggregate(as.numeric(as.character(full.table$V4)), by=list(full.table$V3), FUN=sum)
aggregate(as.numeric(as.character(full.table$V4)), by=list(full.table$V2, full.table$V3), FUN=sum)
aggregate(as.numeric(as.character(full.table$V4)), by=list(full.table$V1, full.table$V3), FUN=sum)




##########Calculate evenness


done <- 0
sel.mix <- list()
IDs <- 1:length(num.comp)
amounts <- 1:length(num.comp)
mix <- cbind(IDs, amounts)
Es <- vector()
while(done==0){
  y <- vector(mode="numeric", length=num.comp)
  for (i in 1:(num.comp-1)){
    y[i] <- sample(min.conc:(conc-sum(y)-((num.comp-i)*min.conc)), 1)
  }
  y[num.comp] <- conc-sum(y)
  

  
  pi <- y/sum(y)
  pi2 <- pi * log(pi)
  D <- sum(pi2)*-1 
  E <- D/log(num.comp)
  Es <- c(Es, E)
  if(E<= 0.22 && E >= 0.18){
    done <- 0
    sel.mix <- c(sel.mix, as.list(y))
  }
}  


#Option 3 for generating y, need to add more y's if we increase richness    
y1 <- sample(min.conc:(conc-((num.comp-1)*min.conc)), 1)
y2 <- sample(min.conc:(conc-y1-((num.comp-2)*min.conc)), 1); y <- c(y1, y2)
y3 <- sample(min.conc:(conc-sum(y)-((num.comp-3)*min.conc)), 1); y <- c(y, y3)
y4 <- sample(min.conc:(conc-sum(y)-((num.comp-4)*min.conc)), 1); y <- c(y, y4)
y5 <- sample(min.conc:(conc-sum(y)-((num.comp-5)*min.conc)), 1); y <- c(y, y5)
y6 <- conc-sum(y)



y <- NA
y1 <- sample(1:(conc-(num.comp-1)), 1)
y2 <- sample(1:(conc-y1-(num.comp-2)), 1); y <- c(y1, y2)
y3 <- sample(1:(conc-sum(y)-(num.comp-3)), 1); y <- c(y, y3)
y4 <- sample(1:(conc-sum(y)-(num.comp-4)), 1); y <- c(y, y4)
y5 <- sample(1:(conc-sum(y)-(num.comp-5)), 1); y <- c(y, y5)
y6 <- conc-sum(y)




sum(y)


y3 <- runif(1, min=0, max=1-sum(y)); y <- c(y,y3)
y4 <- runif(1, min=0, max=1-sum(y)); y <- c(y,y4)
y5 <- runif(1, min=0, max=1-sum(y)); y <- c(y,y5)
y6 <- runif(1, min=0, max=1-sum(y)); y <- c(y,y6)
y7 <- runif(1, min=0, max=1-sum(y)); y <- c(y,y7)
y8 <- runif(1, min=0, max=1-sum(y)); y <- c(y,y8)


rand_vect <- function(N, M, sd = 1, pos.only = TRUE) {
  vec <- rnorm(N, M/N, sd)
  if (abs(sum(vec)) < 0.01) vec <- vec + 1
  vec <- round(vec / sum(vec) * M)
  deviation <- M - sum(vec)
  for (. in seq_len(abs(deviation))) {
    vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
  }
  if (pos.only) while (any(vec < 0)) {
    negs <- vec < 0
    pos  <- vec > 0
    vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
    vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
  }
  vec
}


t <- runif(num.comp-1, min=0, max=conc)
xy <- c(0, sort(t), conc)
di <- vector()
for (i in 1:num.comp){
  di[i] <- xy[i+1]-xy[i]
}


di <- xy{}
hist(xy(:,1),100)

#Option 1 for generating y
#x <- sample(1:100, num.comp)
#y <- x/sum(x) * conc
#Option 2 for generating y  
#x <- runif(num.comp-1, min=0, max=conc)
#x2 <- c(0, sort(x), conc)
#y <- vector()
#for (i in 1:num.comp){
#  y[i] <- x2[i+1]-x2[i]
#}


Generate X???1 random numbers in [0,Y] and insert 0 and 100 to produce a list of X+1X+1 numbers (xi)Xi=0(xi)i=0X. Sort the list, and then calculate their successive differences

di=xi+1???xi


##Add the final choices for high, mid, low to final list
names.sel <- c(names(sel.mix), paste("high", levels[i], sep=""),paste("mid", levels[i], sep=""),
               paste("low", levels[i], sep=""))    
sel.combos <- c(sel.combos, list(high.set.final, mid.set.final, low.set.final))
names(sel.combos) <- names.sel




x<-runif(3)
> x/sum(x)
[1] 0.1130642 0.4098608 0.4770750



randomValue = rand;
value = CurrentVector + randomValue;
if 0.005 <= value && value <= 0.03
fprintf('Criteria satisfied with a random value of %f\n', randomValue);
break; % out of the while loop
else
  fprintf('Criteria NOT satisfied with a random value of %f\n', randomValue);
end
end





###For richness vs evenness experiments





i=2




i=2
z=1



    
    
    
        
}
high.sel <- matrix(,ncol=levels[i], nrow=3)
while(z<3){
  high.sel[z,] <- as.character(unname(unlist(high[z,1:levels[i]])))
  if (
    sum(duplicated(c(as.character(unname(unlist(high[z,1:levels[i]]))), 
                     as.character(unname(unlist(high[z+1,1:levels[i]]))))))
    >2){
    z=z+1
  }
  else{
    
  }
}    








####For richness experiments-----------

###Set parameters for groupings
cut <- .10  #the proportion cutoff level for top, middle, and low
#N.combos <- 6 #number of random samples to pull from each group
levels <- c(2,4,6,8,10) #The levels of richness to use; we would also have 1 and 14 but there is only one way to make those

IDs <- colnames(d.sub)
all.combos <- list()
sel.combos <- list()
for (i in 1:length(levels)){
  d <- t(combn(IDs, levels[i]))
  d <- cbind(d, NA)
  
  for (j in 1:nrow(d)){
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
  names.all <- c(names(all.combos), paste("level",levels[i], sep=""))
  all.combos <- c(all.combos, list(d2))
  names(all.combos) <- names.all
  
  
  high <- d2[1:round((nrow(d2)*cut)),]
  
  #for mid group
  up <- round(nrow(d2)/2 + nrow(high)/2)
  down <- round(nrow(d2)/2 - nrow(high)/2)
  mid <- d2[up:down,]
  
  low <- d2[(nrow(d2)-nrow(high)):nrow(d2),]
  
  
  
  #high.sample <- high[sample(nrow(high), N.combos, replace=FALSE),]
  #mid.sample <- mid[sample(nrow(mid), N.combos, replace=FALSE),]
  #low.sample <- low[sample(nrow(low), N.combos, replace=FALSE),] 
  
  high.sample <- high
  mid.sample <- mid
  low.sample <- low
  
  names.sel <- c(names(sel.combos), paste("high", levels[i], sep=""),paste("mid", levels[i], sep=""),
                 paste("low", levels[i], sep=""))    
  sel.combos <- c(sel.combos, list(high.sample, mid.sample, low.sample))
  names(sel.combos) <- names.sel
  
}



sel.combos

file.names <- vector()
for(i in 1:length(sel.combos)){
  file.name <- paste("2016PD_", names(sel.combos[i]), ".csv", sep="")
  file.names <- c(file.names, file.name)
}

for(i in 1:length(file.names)){
  write.csv(sel.combos[[i]], file = file.names[i], row.names = FALSE)
}


if(nrow(d2) > cut*3){  #for those lists that have lots of possibilities, take a subset
  high <- d2[1:cut,]
  mid <- d2[(nrow(d2)/2 - cut/2):(nrow(d2)/2 + cut/2),]
  low <- d2[(nrow(d2)-cut):nrow(d2),]
} else {   ##for sets of compound combinations < 150 
  high <- d2[1:10,]
  mid <- d2[(nrow(d2)/2 -5):(nrow(d2)/2 +5),]
  low <- d2[(nrow(d2)-10):nrow(d2),]
}


high.sets <- combn(rownames(high), 3) #all possible sets of 3




mid.sets <- combn(rownames(mid), 3) #all possible sets of 3

#check which sets of 3 meet the criteria
ok.mid.sets <- list()
names.ok.mid.sets <- list()
c <- 0
for (k in 1:ncol(mid.sets)){
  #Make a combined list of all compounds in all three sets
  t <- c(as.character(unname(unlist(mid[mid.sets[1,k],1:levels[i]]))),
         as.character(unname(unlist(mid[mid.sets[2,k],1:levels[i]]))),
         as.character(unname(unlist(mid[mid.sets[3,k],1:levels[i]]))))
  
  #make a list of overlaps between each pairwise comparison
  one.two <- intersect(as.character(unname(unlist(mid[mid.sets[1,k],1:levels[i]]))),
                       as.character(unname(unlist(mid[mid.sets[2,k],1:levels[i]]))))
  two.three <- intersect(as.character(unname(unlist(mid[mid.sets[2,k],1:levels[i]]))),
                         as.character(unname(unlist(mid[mid.sets[3,k],1:levels[i]]))))
  one.three <- intersect(as.character(unname(unlist(mid[mid.sets[1,k],1:levels[i]]))),
                         as.character(unname(unlist(mid[mid.sets[3,k],1:levels[i]]))))
  if(   ##choose only combos that meet criteria for g and nr
    all(table(t)<=nr) &&
    length(one.two) < g  &&
    length(one.three) < g  &&
    length(two.three) < g
  ){  ##add those combos to a list, ordered by the total sum of MeanDist
    c <- c+1
    ok.mid.sets <- c(ok.mid.sets, list(as.matrix(rbind(mid[mid.sets[1,k],] ,
                                                       mid[mid.sets[2,k],] , 
                                                       mid[mid.sets[3,k],]))))
    b <- as.data.frame(ok.mid.sets[[c]])
    names.ok.mid.sets <- c(names.ok.mid.sets, sum(as.numeric(as.character(b$MeanDist))))
    names(ok.mid.sets) <- names.ok.mid.sets
    ok.mid.sets <- ok.mid.sets[order(as.numeric(names(ok.mid.sets)))] 
  }
}

##ok.mid.sets is now an ordered list of all possible sets that meet criteria
mid.set.final <- ok.mid.sets[[round(length(ok.mid.sets)/2)]]  ##use the one with mid functional diversity


low.sets <- combn(rownames(low), 3) #all possible sets of 3

#check which sets of 3 meet the criteria
ok.low.sets <- list()
names.ok.low.sets <- list()
c <- 0
for (k in 1:ncol(low.sets)){
  #Make a combined list of all compounds in all three sets
  t <- c(as.character(unname(unlist(low[low.sets[1,k],1:levels[i]]))),
         as.character(unname(unlist(low[low.sets[2,k],1:levels[i]]))),
         as.character(unname(unlist(low[low.sets[3,k],1:levels[i]]))))
  
  #make a list of overlaps between each pairwise comparison
  one.two <- intersect(as.character(unname(unlist(low[low.sets[1,k],1:levels[i]]))),
                       as.character(unname(unlist(low[low.sets[2,k],1:levels[i]]))))
  two.three <- intersect(as.character(unname(unlist(low[low.sets[2,k],1:levels[i]]))),
                         as.character(unname(unlist(low[low.sets[3,k],1:levels[i]]))))
  one.three <- intersect(as.character(unname(unlist(low[low.sets[1,k],1:levels[i]]))),
                         as.character(unname(unlist(low[low.sets[3,k],1:levels[i]]))))
  if(   ##choose only combos that meet criteria for g and nr
    all(table(t)<=nr) &&
    length(one.two) < g  &&
    length(one.three) < g  &&
    length(two.three) < g
  ){  ##add those combos to a list, ordered by the total sum of MeanDist
    c <- c+1
    ok.low.sets <- c(ok.low.sets, list(as.matrix(rbind(low[low.sets[1,k],] ,
                                                       low[low.sets[2,k],] , 
                                                       low[low.sets[3,k],]))))
    b <- as.data.frame(ok.low.sets[[c]])
    names.ok.low.sets <- c(names.ok.low.sets, sum(as.numeric(as.character(b$MeanDist))))
    names(ok.low.sets) <- names.ok.low.sets
    ok.low.sets <- ok.low.sets[order(as.numeric(names(ok.low.sets)))] 
  }
}

##ok.low.sets is now an ordered list of all possible sets that meet criteria
low.set.final <- ok.low.sets[[length(ok.low.sets)]]  ##use the one with highest functional diversity



x <- matrix(1:12,3,4)
y <- x+100
dim(abind(x,y,along=0))     # binds on new dimension before first
dim(abind(x,y,along=1))     # binds on first dimension
dim(abind(x,y,along=1.5))
dim(abind(x,y,along=2))
dim(abind(x,y,along=3))
dim(abind(x,y,rev.along=1)) # binds on last dimension
dim(abind(x,y,rev.along=0)) # binds on new dimension after last

