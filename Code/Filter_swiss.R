#Remove Swiss data and make new ASV table, Metadata table, Taxonomic list 
#and a Pruned phylogeny for use in all downstream analysis

####set up####
#install packages 
library(tidyverse)
library(dplyr)
library(phyloseq)
library(picante)


#load data
#these are the cleaned, filtered, and rarefied data for all sites
asv16s.physeq3 <- readRDS("Data/dMTM_asv16s.physeq3.RDS")

#Taxonomic list
tax <- as.data.frame(tax_table(asv16s.physeq3))

#Phylogenetic tree for the ASVs
asv16s.tree <-phy_tree(asv16s.physeq3) 
  





#the pitcher bu OTU table (rows are ASVs (n=3806), columns are pitchers (n=693))
asv <- t(data.matrix((otu_table(asv16s.physeq3))))
dim(asv)

#meta is the covariate data associated with each pitcher (n=693)

#function that allows for classification of numeric variables in meta below
is_all_numeric <- function(x) {
  !any(is.na(suppressWarnings(as.numeric(na.omit(x))))) & is.character(x)
}
meta<- HTSSIP::phyloseq2df(asv16s.physeq3, table_func=phyloseq::sample_data)%>%
        filter(!country=='CH')%>%# remove swiss sites
        select(volume_ml,midge_number,mosquito_number,#select columns we want
               fleshfly_number,intactprey_number,latitude,ph,longitude,Actual_Day,site, pitcher_number)%>%
        mutate_if(is_all_numeric, as.numeric)%>%# is_all_numeric from above
        mutate(longitude=replace(longitude, longitude==-1.740490, -71.74049))#fix typo on longitude for Quebec site
str(meta)

#add pitcher covariates to asv table data 
comb_dat<-merge(asv,meta, by=0) %>% #merge asv with meta
          column_to_rownames(var = "Row.names")


#separate our only pitcher-by-asv matrix and remove asvs with 0 col sum (likely swiss) (484 samples-by-3323 asvs)
filtered_asv<-comb_dat[,1:3806]%>%#these are the columns that are ASVs
              select_if(negate(function(col) is.numeric(col) && sum(col) == 0))#remove asvs with 0 col sum


#filter taxonomic table bsed on the asv
filtered_tax<-tax %>%  
              filter(ASVs %in% colnames(filtered_asv)) 



#Prune the Tree
#prune tree asv table with swiss filtered out
pdTree<-prune.sample(filtered_asv, asv16s.tree)




#check dimensions
dim(filtered_tax) #3376 ASVs
dim(meta) #515 samples
dim(filtered_asv) #3376 ASV x 515 samples


#save out processed data for use in all downstream analyses
#write.tree(pdTree,"Data/pitcher_plant_tree.nwk")
#write.csv( filtered_asv,"Data/pitcher_plant_ASV_table.csv")
#write.csv(meta, "Data/pitcher_plant_meta_table.csv")
#write.csv(filtered_tax,"Data/pitcher_plant_tax_table.csv")







