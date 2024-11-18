#Phylogenetic Alpha Diversity


####set_up####
#Install and load relevant pacakges 
library(tidyverse)
library(dplyr)
library(ape)
library(vegan)
library(picante)

#call in data
#the pitcher by ASV table (rowa are pitchers (n=515, columns are ASVs (n=3376))
asv <- read.table ("Data/final/pitcher_plant_ASV_table.csv", sep= ",", check.names = F,row=1, header =T)
dim(asv)

#meta is the covariate data associated with each pitcher (n=515)
meta <- read.table("Data/final/pitcher_plant_meta_table.csv", sep= ",", header=T)
dim(meta)

#phylogenetic tree for the ASVs
asv16s.tree <- read.tree("Data/final/pitcher_plant_tree.nwk")

#ses.pd run null model for each site separately
####subset for QC####
asv2<-as.data.frame(asv[grep('^QC', rownames(asv)),])%>%
      mutate_if(is.character, as.numeric)%>%
      select_if(negate(function(col) is.numeric(col) && sum(col) == 0))#remove asvs that are not present at site

#prune tree to site level asv
pdTree<-prune.sample(asv2, asv16s.tree)

#calculate PD and SESPD for site
pitcher_pd<-ses.pd(asv2, pdTree, null.model="taxa.labels", runs=999)

#write.csv(pitcher_pd, "Data/final/PD_QC.csv")

####subset for FL####
asvFL<-as.data.frame(asv[grep('^FL', rownames(asv)),])%>%
      mutate_if(is.character, as.numeric)%>%
      select_if(negate(function(col) is.numeric(col) && sum(col) == 0))#remove asvs that are not present at site

#prune tree to site level asv
pdTree<-prune.sample(asvFL, asv16s.tree)

#calculate PD and SESPD for site
pitcher_pd<-ses.pd(asvFL, pdTree, null.model="taxa.labels", runs=999)

#write.csv(pitcher_pd, "Data/final/PD_FL.csv")

####subset for GA####
asvGA<-as.data.frame(asv[grep('^GA', rownames(asv)),])%>%
      mutate_if(is.character, as.numeric)%>%
      select_if(negate(function(col) is.numeric(col) && sum(col) == 0))#remove asvs that are not present at site

#prune tree to site level asv
pdTree<-prune.sample(asvGA, asv16s.tree)

#calculate PD and SESPD for site
pitcher_pd<-ses.pd(asvGA, pdTree, null.model="taxa.labels", runs=999)

#write.csv(pitcher_pd, "Data/final/PD_GA.csv")

####subset for IN####
asvIN<-as.data.frame(asv[grep('^IN', rownames(asv)),])%>%
        mutate_if(is.character, as.numeric)%>%
        select_if(negate(function(col) is.numeric(col) && sum(col) == 0))#remove asvs that are not present at site

#prune tree to site level asv
pdTree<-prune.sample(asvIN, asv16s.tree)

#calculate PD and SESPD for site
pitcher_pd<-ses.pd(asvIN, pdTree, null.model="taxa.labels", runs=999)

#write.csv(pitcher_pd, "Data/final/PD_IN.csv")

####subset for MA####
asvMA<-as.data.frame(asv[grep('^MA', rownames(asv)),])%>%
      mutate_if(is.character, as.numeric)%>%
      select_if(negate(function(col) is.numeric(col) && sum(col) == 0))#remove asvs that are not present at site

#prune tree to site level asv
pdTree<-prune.sample(asvMA, asv16s.tree)

#calculate PD and SESPD for site
pitcher_pd<-ses.pd(asvMA, pdTree, null.model="taxa.labels", runs=999)

#write.csv(pitcher_pd, "Data/final/PD_MA.csv")

####subset for NC ####
asvNC<-as.data.frame(asv[grep('^NC', rownames(asv)),])%>%
      mutate_if(is.character, as.numeric)%>%
      select_if(negate(function(col) is.numeric(col) && sum(col) == 0))#remove asvs that are not present at site

#prune tree to site level asv
pdTree<-prune.sample(asvNC, asv16s.tree)

#calculate PD and SESPD for site
pitcher_pd<-ses.pd(asvNC, pdTree, null.model="taxa.labels", runs=999)

#write.csv(pitcher_pd, "Data/final/PD_NC.csv")

####subset for WI ####
asvWI<-as.data.frame(asv[grep('^WI', rownames(asv)),])%>%
      mutate_if(is.character, as.numeric)%>%
      select_if(negate(function(col) is.numeric(col) && sum(col) == 0))#remove asvs that are not present at site

#prune tree to site level asv
pdTree<-prune.sample(asvWI, asv16s.tree)

#calculate PD and SESPD for site
pitcher_pd<-ses.pd(asvWI, pdTree, null.model="taxa.labels", runs=999)

#write.csv(pitcher_pd, "Data/final/PD_WI.csv")