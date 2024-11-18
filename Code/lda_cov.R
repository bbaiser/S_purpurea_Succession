# LDAcov Models

####set up####
#install packages 
library(devtools)
devtools::install_github("gilsonshimizu/ldacov")# ladcov from github
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
library(arsenal)

#load data
#the pitcher bu OTU table (columns are ASVs (n=3376), rows are pitchers (n=515))
asv <- read.csv("Data/final/pitcher_plant_ASV_table.csv", row=1)
dim(asv)

#meta is the covariate data (n=11 varaibles/columns) associated with each pitcher (n=515 rows)
meta <- read.csv("Data/final/pitcher_plant_meta_table.csv", row=1)


#add pitcher covariates to asv table data 
comb_dat<-merge(asv,meta, by=0) %>% #merge asv with meta
          column_to_rownames(var = "Row.names")%>%
          drop_na()#drop rows where the covariate has an NA


#separate our only pitcher-by-asv matrix (484 samples-by-3323 asvs)
comm<-comb_dat[,1:3376]%>%
      select_if(negate(function(col) is.numeric(col) && sum(col) == 0))#remove asvs with 0 col sum 

comm2<-as.matrix(comm)#make data frame a matrix
dim(comm2)#check dimensions

#separate our pitcher-by-predictor variable matrix
xvars<-apply(comb_dat[,3377:3385], 2, as.numeric)#rows 3807-3815 are the covariates

#scale covariates
xvars<-scale(as.matrix(xvars))

#add intercept
xvars<-cbind(int = 1, xvars)

#save out covariates if desired
#write.csv(xvars, "Outputs/xvars2_matrix.csv")


####LDA with out data, no covariates to find group#### 
#TSB prior


#run lda with no covariates to assess group# (we did this with 10 goups and found that 6 groups contains >75% elements)
set.seed(1)
lda_no_covariates2=gibbs.LDA(y=comm2,
                             ncomm=6,#number of clusters to recover
                             ngibbs=5000,#number of iterations for the Gibbs sampler
                             nburn=2500,# number of iterations to be discarded as burn-in;
                             psi=0.01,#parameter for the prior Dirichlet distribution used for ϕk
                             gamma=0.1)#parameter for the TSB prior used for theta

#save out LDA object
#saveRDS(lda_no_covariates2, file = "Data/lda_no_cov_full_6.RDS")

#call in lda object if not running
#lda_no_covariates2<-readRDS("Data/lda_no_cov_full_6.RDS")

#check convergence (converges around 500 iterations)
plot(lda_no_covariates2$llk,type='l',xlab='Iterations',ylab='Log-likelihood')

#get theta matrix
array.lsk.init=lda_no_covariates2$array.lsk
nlk=apply(array.lsk.init,c(1,3),sum)
theta=nlk/apply(nlk,1,sum)
colnames(theta)=paste0('Cluster',1:6)
rownames(theta)=rownames(comm2)
head(round(theta,2))


#assess cluster with box plot and histogtram
boxplot(theta,ylab=expression(theta),xlab='Clusters',ylim=c(0,1))


#elements contained in clusters, used to determine how many clusters to keep
cumsum1=cumsum(colSums(theta,na.rm=TRUE)/nrow(comm2))
cumsum1[1:6]

#scree plot
plot(colSums(theta,na.rm=TRUE)/nrow(comm2))


#phi matrix
seq1=1500:2500#number of iterations after convergence to use
dim(lda_no_covariates2$phi)
phi <- matrix(colMeans(lda_no_covariates2$phi[seq1,]),nrow=10,ncol=ncol(comm2))
rownames(phi)=paste0('Cluster',1:10)
colnames(phi)=paste0(colnames(comm2))
head(round(phi[,1:10],2))

#plot phi matrix
data <- expand.grid(X=1:ncol(comm2), Y=1:10)
data$proportion <- as.vector(phi)
head(as.vector(phi))

ggplot(data, aes(X, Y, fill= proportion)) + geom_tile() + 
  scale_fill_gradient(low="green", high="darkblue",na.value="green") +
  labs(x = "Category distribution") + labs(y = "Cluster") +
  labs(fill = " ")+theme_minimal()+
  ggtitle("Estimated Phi matrix")

#identify main cluster categories
phi.max=matrix(NA,nrow=10,ncol=ncol(comm2))
for (i in 1:10){phi.max[i,]=apply(phi[-i,], 2, max)}
results=phi/phi.max
colnames(results)=colnames(phi)
rownames(results)=rownames(phi)
head(round(results[,1:10],2))

#order most importatn asv's per cluster'
max_categ_cluster=apply(phi/phi.max, 1, function(x) names(sort(x,decreasing=TRUE)))
head(max_categ_cluster)
topten<-max_categ_cluster[1:10,]

#write.csv(topten, "Outputs/toptenv2.csv")

####LDA WITH COVARIATES####
set.seed(1)
lda_with_covariates2 <- gibbs.LDA.cov(ncomm=6,#number of clusters above tsb
ngibbs=1000,#number of iterations for the Gibbs sampler
y=comm2,# pitcher by asv matrix
xmat=xvars,# pitcher by environment matrix
phi.prior=0.01,#parameter to be used in the Dirichlet prior for Φ if this matrix is estimated;
array.lsk.init=lda_no_covariates2$array.lsk,#initial values for the array.lsk array
var.betas=rep(dim(comm2)[1],ncol(xvars)),# variance parameters for the normal distribution priors for the regression coefficients;
phi.init=lda_no_covariates2$phi,#use posterior samples from no covarites if below is false
estimate.phi=TRUE)

#save out lda object
#saveRDS(lda_with_covariates2, file = "Data/dlda_with_cov_full_6.RDS")

#read in lda object if not running
lda_with_covariates2<-readRDS("Data/dlda_with_cov_full_6.RDS")

# check convergance
plot(lda_with_covariates2$llk,type='l',xlab="iterations",ylab='log-likelihood')

#theta matrix
seq1=100:1000
tmp=matrix(colMeans(lda_with_covariates2$nlk[seq1,]),nrow=nrow(comm2),ncol=6)
theta <- tmp/rowSums(tmp)
colnames(theta)=paste0('Cluster',1:6)
rownames(theta)=rownames(comm2)
head(round(theta,2))



#plot theta matrix (Appendix S3:FigS2)
data <- expand.grid(X=1:nrow(comm2), Y=1:6)
data$proportion <- as.vector(theta)

ggplot(data, aes(X, Y, fill= proportion)) + geom_tile() + 
  scale_fill_gradient(low="green", high="darkblue",na.value="green") +
  labs(x = "Pitcher Samples") + labs(y = "Group") +
  labs(fill = " ")+theme_minimal() +
  ggtitle("Estimated Theta matrix")+
  guides(fill=guide_legend(title="Proportion"))+
  scale_y_continuous(breaks= c(1,2,3,4,5,6),labels = c(1,2,3,4,5,6) )



#phi matrix
phi <- matrix(colMeans(lda_with_covariates2$phi[seq1,]),nrow=6,ncol=ncol(comm2))
rownames(phi)=paste0('Cluster',1:6)
colnames(phi)=paste0(colnames(comm2))
head(round(phi[,1:6],2))

#plot phi matrix (Appendix S3:FigS2)
data <- expand.grid(X=1:ncol(comm2), Y=1:6)
data$proportion <- as.vector(phi)

ggplot(data, aes(X, Y, fill= proportion)) + geom_tile() + 
  scale_fill_gradient(low="green", high="darkblue",na.value="green") +
  labs(x = "ASV") + labs(y = "Group") +
  labs(fill = " ")+theme_minimal()+
  ggtitle("Estimated Phi matrix")+
  guides(fill=guide_legend(title="Relative abundance"))+
  scale_y_continuous(breaks= c(1,2,3,4,5,6),labels = c(1,2,3,4,5,6))

#most important ASVs for each cluster
phi.max=matrix(NA,nrow=6,ncol=ncol(comm2))
for (i in 1:6){phi.max[i,]=apply(phi[-i,], 2, max)}
results=phi/phi.max
colnames(results)=colnames(phi)
rownames(results)=rownames(phi)
head(round(results[,1:6],2))

#order most important asv's per cluster'
max_categ_cluster=apply(phi/phi.max, 1, function(x) names(sort(x,decreasing=TRUE)))
head(max_categ_cluster)
topten<-max_categ_cluster[1:10,]

#save out list of top 10 asvs per cluster
#write.csv(topten, "Outputs/dtoptenv_with_cov_full_6.csv")


# Extract the posterior mean of the regression coefficients from the LDA (used to make table/Fig 3)
tmp=colMeans(lda_with_covariates2$betas[seq1,])
betas=matrix(tmp,ncol=6)
colnames(betas)=paste0('Cluster',1:6)
rownames(betas)=colnames(xvars)
coefs<-t(round(betas,2))
#write.csv(coefs, "Outputs/finals/lda_coefs.csv")


#get CIs for posterior distributions
post <- data.frame(lda_with_covariates2$betas[seq1,])
CIs<-as.data.frame(sapply(post,ci,ci=0.95))#manually went through to identify which didn't cross zero



#plot individual pitcher from Wisconsin for Manuscript Figure 6
plot_test<-as.data.frame(cbind(theta,xvars))%>% #merge asv with meta
          #select(-(int))%>%
          rownames_to_column(var = "pitcher")%>%
          filter(stringr::str_detect(pitcher, "WI_P02"))%>%
          select(Cluster1,Cluster2, Cluster3,Cluster4,Cluster5,Cluster6)%>%
          rename("Group 1"=Cluster1,"Group 2"=Cluster2,"Group 3"=Cluster3,
                 "Group 4"=Cluster4,"Group 5"=Cluster5,"Group 6"=Cluster6)%>%
          mutate(dayraw = c(3, 7, 14, 28, 42, 56, 70, 84))


#Figure 6 main text
#format for plot
dd_tidyr <- pivot_longer(plot_test, cols = -c(dayraw), names_to = "LDA Group")

ggplot(dd_tidyr) +
  geom_line(aes(x = dayraw, y = value, colour = `LDA Group`),size=1) +
  scale_colour_manual(values = c("black","red", "green", "blue","purple","pink","orange"))+
  theme_bw()+
  labs(y= "Proportion of Community", x= "Time (days)")



#### plot posterior distributions per cluster, run for each cluster and then make panel figure below####

#Format data for  cluster/group 1
number <- 1:length(seq1)
aux1 <- lda_with_covariates2$betas[seq1,2:10]#number reflect groups of 8 (i.e. 2-10, 12-20, etc)
colnames(aux1)=colnames(xvars[,2:10])
aux2<-data.frame(number,aux1)%>%
      rename( time = Actual_Day)
aux3<- melt(aux2, id=c("number"))


#plots individual groups. (you have to re)
one<-ggplot(aux3, aes(x = value, y = variable, fill = variable)) +
  geom_density_ridges(scale=0.9) +
  theme_ridges() + 
  theme(legend.position = "N" ) +
  xlab("") +
  ylab("")+
  ggtitle("Group 1")


#Format data for plotcluster/group 2
number <- 1:length(seq1)
aux1 <- lda_with_covariates2$betas[seq1,12:20]#number reflect groups of 8 (i.e. 2-10, 12-20, etc)
colnames(aux1)=colnames(xvars[,2:10])
aux2<-data.frame(number,aux1)%>%
  rename( time = Actual_Day)
aux3<- melt(aux2, id=c("number"))


#plots individual groups
two<-ggplot(aux3, aes(x = value, y = variable, fill = variable)) +
  geom_density_ridges(scale=0.9) +
  theme_ridges() + 
  theme(legend.position = "N" ) +
  xlab("") +
  ylab("")+
  ggtitle("Group 2")

#Format data for plotcluster/group 3
number <- 1:length(seq1)
aux1 <- lda_with_covariates2$betas[seq1,22:30]#number reflect groups of 8 (i.e. 2-10, 12-20, etc)
colnames(aux1)=colnames(xvars[,2:10])
aux2<-data.frame(number,aux1)%>%
  rename( time = Actual_Day)
aux3<- melt(aux2, id=c("number"))


#plots individual groups. 
three<-ggplot(aux3, aes(x = value, y = variable, fill = variable)) +
  geom_density_ridges(scale=0.9) +
  theme_ridges() + 
  theme(legend.position = "N" ) +
  xlab("") +
  ylab("")+
  ggtitle("Group 3")

#Format data for plot cluster/group 4
number <- 1:length(seq1)
aux1 <- lda_with_covariates2$betas[seq1,32:40]#number reflect groups of 8 (i.e. 2-10, 12-20, etc)
colnames(aux1)=colnames(xvars[,2:10])
aux2<-data.frame(number,aux1)%>%
  rename( time = Actual_Day)
aux3<- melt(aux2, id=c("number"))


#plots individual groups
four<-ggplot(aux3, aes(x = value, y = variable, fill = variable)) +
  geom_density_ridges(scale=0.9) +
  theme_ridges() + 
  theme(legend.position = "N" ) +
  xlab("") +
  ylab("")+
  ggtitle("Group 4")

#Format data for plot cluster/group 5
number <- 1:length(seq1)
aux1 <- lda_with_covariates2$betas[seq1,42:50]#number reflect groups of 8 (i.e. 2-10, 12-20, etc)
colnames(aux1)=colnames(xvars[,2:10])
aux2<-data.frame(number,aux1)%>%
  rename( time = Actual_Day)
aux3<- melt(aux2, id=c("number"))


#plots individual groups
five<-ggplot(aux3, aes(x = value, y = variable, fill = variable)) +
  geom_density_ridges(scale=0.9) +
  theme_ridges() + 
  theme(legend.position = "N" ) +
  xlab("") +
  ylab("")+
  ggtitle("Group 5")

#Format data for plot cluster/group 6
number <- 1:length(seq1)
aux1 <- lda_with_covariates2$betas[seq1,52:60]#number reflect groups of 8 (i.e. 2-10, 12-20, etc)
colnames(aux1)=colnames(xvars[,2:10])
aux2<-data.frame(number,aux1)%>%
  rename( time = Actual_Day)
aux3<- melt(aux2, id=c("number"))


#plots individual groups
six<-ggplot(aux3, aes(x = value, y = variable, fill = variable)) +
    geom_density_ridges(scale=0.9) +
    theme_ridges() + 
    theme(legend.position = "N" ) +
    xlab("") +
    ylab("")+
    ggtitle("Group 6")




#combine figures into panel Appendix S3:Fig S1
figure <- ggarrange(one, two, three,four,five,six,
                    labels = c("A)", "B)", "C)", "D)", "E)", "F)"),
                    ncol = 2, nrow = 3)

annotate_figure(figure, bottom = text_grob("Coefficent Estimate", color = "black", size = 14))


