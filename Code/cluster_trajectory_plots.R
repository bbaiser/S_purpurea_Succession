
####set up####
#install packages  (ladcov from github)
R.version
library(devtools)
devtools::install_github("gilsonshimizu/ldacov")
library(ldacov)
library(tidyverse)
library(dplyr)
library(phyloseq)
library(textmineR)
library(HTSSIP)
library(reshape2)
library(ggridges)
library(bayestestR)
library(hutils)
library(ggpubr)
#load data
#these are the cleaned, filtered, and rarefied data from Jessica Berdardin in "MTM_Diversity_Analysis_2023.Rmd"
asv16s.physeq3 <- readRDS("Data/dMTM_asv16s.physeq3.RDS")

#the pitcher bu OTU table (rows are ASVs (n=3806), columns are pitchers (n=693))
asv <- t(data.matrix((otu_table(asv16s.physeq3))))


#meta is the covariate data associated with each pitcher (n=693)
meta<- HTSSIP::phyloseq2df(asv16s.physeq3, table_func=phyloseq::sample_data)%>%
        filter(!country=='CH')%>%# remove swiss sites
        select(volume_ml,midge_number,mosquito_number,#select columns we want
               fleshfly_number,intactprey_number,latitude,ph,longitude,day,Actual_Day)%>%
        mutate_if(is.character, as.numeric)%>%
        mutate(longitude=replace(longitude, longitude==-1.740490, -71.74049))#fix typo on longitude

#add pitcher covariates to asv table data and remove swiss sites
comb_dat<-merge(asv,meta, by=0) %>% #merge asv with meta
          column_to_rownames(var = "Row.names")%>%
          drop_na()#drop rows where the covariate has an NA


#separate our only pitcher-by-asv matrix (484 samples-by-3323 asvs)
comm<-comb_dat[,1:3806]%>%
  select_if(negate(function(col) is.numeric(col) && sum(col) == 0))#remove asvs with 0 col sum (likely swiss)

comm2<-as.matrix(comm)#make data fram a matrix

dim(comm2)#check dimensions

#separate our pitcher-by-predictor variable matrix
xvars<-apply(comb_dat[,3807:3816], 2, as.numeric)#rows 3807-3815 are the covariates

#scale covariates
#xvars<-scale(as.matrix(xvars))

#add intercept
#xvars<-cbind(int = 1, xvars)



#call in LDA OBJECT
lda_with_covariates2<-readRDS("Data/dlda_with_cov_full_6.RDS")

#theta matrix
seq1=100:1000
tmp=matrix(colMeans(lda_with_covariates2$nlk[seq1,]),nrow=nrow(comm2),ncol=6)
theta <- tmp/rowSums(tmp)
colnames(theta)=paste0('Cluster',1:6)
rownames(theta)=rownames(comm2)
head(round(theta,2))



#plot individual pitchers
plot_test<-as.data.frame(cbind(theta,xvars))%>% #merge asv with meta
          #select(-(int))%>%
          rownames_to_column(var = "pitcher")%>%
          mutate(pitcherID=str_sub(pitcher, end = -3))%>%
          select(-(pitcher))%>%
          select(Cluster2,Actual_Day,pitcherID,latitude)#replace cluster number to make plots for each cluster 

meadz<-plot_test%>%
       group_by(latitude, Actual_Day) %>% 
       summarise_at(vars(Cluster2), median)#replace cluster number to make plots for each cluster 


##before running each of these plots, replace cluster number (e.g.., Cluster2, Cluster3) above in the first 
one<-ggplot() + 
    geom_line(data=plot_test,aes(x=Actual_Day,y=Cluster1, group = interaction(pitcherID,latitude), color = latitude), alpha = 0.2)+
    scale_colour_gradient(low = "red",high = "blue")+
    geom_line(data=meadz,aes(x=Actual_Day,y=Cluster1,group = latitude, color = latitude),linewidth=1.5)+
    ylab("Proportion")+
    ggtitle("Group 1")+
    theme(plot.title = element_text(hjust = 0.5))+
    theme_bw()+
    theme(axis.title.x = element_blank())

two<-ggplot() + 
  geom_line(data=plot_test,aes(x=Actual_Day,y=Cluster2, group = interaction(pitcherID,latitude), color = latitude), alpha = 0.2)+
  scale_colour_gradient(low = "red",high = "blue")+
  geom_line(data=meadz,aes(x=Actual_Day,y=Cluster2,group = latitude, color = latitude),linewidth=1.5)+
  ylab("Proportion")+
  ggtitle("Group 2")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  theme(axis.title.x = element_blank())

three<-ggplot() + 
  geom_line(data=plot_test,aes(x=Actual_Day,y=Cluster3, group = interaction(pitcherID,latitude), color = latitude), alpha = 0.2)+
  scale_colour_gradient(low = "red",high = "blue")+
  geom_line(data=meadz,aes(x=Actual_Day,y=Cluster3,group = latitude, color = latitude),linewidth=1.5)+
  ylab("Proportion")+
  ggtitle("Group 3")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  theme(axis.title.x = element_blank())

four<-ggplot() + 
  geom_line(data=plot_test,aes(x=Actual_Day,y=Cluster4, group = interaction(pitcherID,latitude), color = latitude), alpha = 0.2)+
  scale_colour_gradient(low = "red",high = "blue")+
  geom_line(data=meadz,aes(x=Actual_Day,y=Cluster4,group = latitude, color = latitude),linewidth=1.5)+
  ylab("Proportion")+
  ggtitle("Group 4")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  theme(axis.title.x = element_blank())

five<-ggplot() + 
  geom_line(data=plot_test,aes(x=Actual_Day,y=Cluster5, group = interaction(pitcherID,latitude), color = latitude), alpha = 0.2)+
  scale_colour_gradient(low = "red",high = "blue")+
  geom_line(data=meadz,aes(x=Actual_Day,y=Cluster5,group = latitude, color = latitude),linewidth=1.5)+
  ylab("Proportion")+
  ggtitle("Group 5")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  theme(axis.title.x = element_blank())

six<-ggplot() + 
  geom_line(data=plot_test,aes(x=Actual_Day,y=Cluster6, group = interaction(pitcherID,latitude), color = latitude), alpha = 0.2)+
  scale_colour_gradient(low = "red",high = "blue")+
  geom_line(data=meadz,aes(x=Actual_Day,y=Cluster6,group = latitude, color = latitude),linewidth=1.5)+
  ylab("Proportion")+
  ggtitle("Group 6")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  theme(axis.title.x = element_blank())

#combine figures into panel
figure <- ggarrange(one, two, three,four,five,six,
                    labels = c("A)", "B)", "C)", "D)", "E)", "F)"),
                    ncol = 2, nrow = 3, common.legend = TRUE,legend = "right")
figure

annotate_figure(figure, bottom = text_grob("Time (days)", color = "black", size = 14))

