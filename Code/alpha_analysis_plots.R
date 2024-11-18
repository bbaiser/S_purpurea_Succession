# Alpha diversity models and plots 

####set up####

#Install and load relevant packages 

#install.packages("hillR", repos = c('https://daijiang.r-universe.dev', 'https://cloud.r-project.org'))
#install.packages('TMB', type = 'source')
BiocManager::install("phyloseq")

library(ggfortify)
library(phyloseq)
library(tidyverse)
library(dplyr)
library(ape)
library(hillR)
library(vegan)
library(glmmTMB)
library(DHARMa)
library(ggplot2)
library(MuMIn)
library(performance)
library(ggpubr)

#load data
#the pitcher bu OTU table (columns are ASVs (n=3376), rows are pitchers (n=515))
asv <- read.csv("Data/final/pitcher_plant_ASV_table.csv", row=1)
dim(asv)

#meta is the covariate data (n=11 varaibles/columns) associated with each pitcher (n=515 rows)
meta <- read.csv("Data/final/pitcher_plant_meta_table.csv", row=1)
str(meta) 

#read in the individual PD files
data_dir <- "C:/Users/bbaiser/OneDrive - University of Florida/Documents/GitHub/Pplant_succession/Data"#or the name of your directory where the pd files are
csv_files <- fs::dir_ls(data_dir, regexp = "PD_")#grab site PD files
csv_files

#make a column of site level ses.pd
site_PD<-readr::read_csv(csv_files)%>%
        select(...1,pd.obs.z)%>%
        rename (.,Row.names=...1)%>%
        rename(.,site_SES_PD=pd.obs.z)


# calculate alpha diversity  (using hillR, q=0 for richness, q=1 for Shannon entropy, q=2 for inverse Simpson)####

q0<- hill_taxa(asv, q = 0, MARGIN = 1, base = exp(1))
q1<- hill_taxa(asv, q = 1, MARGIN = 1, base = exp(1))
q2<- hill_taxa(asv, q = 2, MARGIN = 1, base = exp(1)) 

as.data.frame(q0)
#make a data frame with all three metrics
alpha_div<-cbind(q0, q1, q2)     


#add pitcher covariates to the diversity data 

#this is the raw data for plotting below
full_data_raw<-merge(alpha_div,meta, by=0) %>% #merge with hill#'s
              left_join(.,site_PD, by="Row.names")%>%#merge with PD from site pool
              column_to_rownames(var = "Row.names")%>%
              rename("Time" = "Actual_Day")%>%
              drop_na()

#this is for scaling an use in model
full_data<-merge(alpha_div,meta, by=0) %>% #merge with hill#'s
          left_join(.,site_PD, by="Row.names")%>%#merge with PD from site pool
          column_to_rownames(var = "Row.names")%>%
          rename("Time" = "Actual_Day")%>%
          drop_na()


#z-scale model variables that are continiuos
for (i in 4:12) {
  full_data[, i] <- (full_data[, i] - mean(full_data[, i])) / sd(full_data[, i])
}         


####Taxonomic alpha diversity models
####q0####
#full model (no dispersion model)
Mfull <- glmmTMB(q0 ~ (Time) +I(Time^2)+
                   volume_ml+
                   intactprey_number+
                   midge_number+
                   fleshfly_number+
                   mosquito_number+
                   latitude+
                   longitude+
                   (1|site) + (1|site:pitcher_number),
                 family = "gaussian",
                 na.action = na.fail,
                 data = full_data)

#check out full model
summary(Mfull)
check_collinearity(Mfull)#no colinearity
r.squaredGLMM(Mfull)


#DHARMa residuals check full model
simulationOutput1 <- simulateResiduals(fittedModel = Mfull,n=10000, plot=F)
plot(simulationOutput1)  #Heteroskedastic


#look at each covariate to identify what is driving heteroskedasticity
testQuantiles(simulationOutput1, predictor =full_data$Time)#yes
testQuantiles(simulationOutput1, predictor =full_data$volume_ml)#yes
testQuantiles(simulationOutput1, predictor =full_data$midge_number)#no
testQuantiles(simulationOutput1, predictor =full_data$mosquito_number)#yes
testQuantiles(simulationOutput1, predictor =full_data$intactprey_number)#no
testCategorical(Mfull, catPred = full_data$latitude)#yes
testCategorical(Mfull, catPred = full_data$longitude)#yes
testCategorical(Mfull, catPred = full_data$fleshfly_number)#no



#full model (with dispersion model based on heteroskedasticity above)
Mfull_disp <- glmmTMB(q0 ~ (Time) +I(Time^2)+
                        volume_ml+
                        intactprey_number+
                        midge_number+
                        fleshfly_number+
                        mosquito_number+
                        latitude+
                        longitude+
                        (1|site) + (1|site:pitcher_number),
                      dispformula=~latitude+longitude+Time+volume_ml+mosquito_number,
                      family = "gaussian",
                      na.action = na.fail,
                      data = full_data)

#check out full model
summary(Mfull_disp)
check_collinearity(Mfull_disp)#no colinearity
r.squaredGLMM(Mfull_disp)#not available because of dispersion model


#DHARMa residuals check full model with dispersion
simulationOutput1 <- simulateResiduals(fittedModel = Mfull_disp,n=10000, plot=F)
plot(simulationOutput1) #residuals look good 


#model selection
stats::drop1(Mfull_disp, test="Chisq")

#best model q0 from drop1
Mbest<- glmmTMB(q0 ~ (Time) +I(Time^2)+
                  latitude+
                  longitude+
                  (1|site) + (1|site:pitcher_number),
                dispformula=~latitude+longitude+Time,#account for the heteroskedsicity in model
                family = "gaussian",
                na.action = na.fail,
                data = full_data)

#check out best-fit model
summary(Mbest)
check_collinearity(Mbest)#no colinearity
r.squaredGLMM(Mbest) #no R2 available because of dispersion model


#DHARMa residuals check best fit model
simulationOutput1 <- simulateResiduals(fittedModel = Mbest,n=10000, plot=F)
plot(simulationOutput1)  #residuals look good


#model sketch q0
range(full_data$Time) # -1.286801  2.215573

MyData <- expand.grid(latitude  = mean(full_data$latitude),
                      longitude = mean(full_data$longitude),
                      Time  = seq(-1.28 , 2.22, length = 16))

X <- model.matrix( ~ Time+I(Time^2)+latitude+longitude,
                   data = MyData)

#Calculate ETA and MU
eta <- X%*% fixef(Mbest)$cond #(use is the best fit model)
mu <- eta 

#Get the standard errors
SE   <- sqrt(diag(X %*%vcov(Mbest)$cond %*% t(X))  )#hand calculates standard error
seup <- eta + 1.96 * SE
selo <- eta - 1.96 * SE
mu
NewData<-cbind(MyData,eta,mu,SE,seup,selo)


dayz<-c(0,7,14,21,28,35,42,49,56,63,70,77,84,91,98,105)#set x axis for days
day2<-full_data_raw$Time#set actual raw time for data points in plot

#plot model sketch
rich<-ggplot(NewData, aes(x=dayz, y=mu))+
      geom_point(data =full_data , 
                 mapping = aes(x = day2, y = q0),color = "gray")+#add raw points
      geom_ribbon(aes(x=dayz,ymin=selo, ymax=seup), fill="#21918c")+
      geom_line(size=1.2, color = "black")+
      xlab("Time (days)")+
      ylab("ASV Richness")+
      ylim(0, 125)+
      theme_bw()+
      theme(text = element_text(size=15),
            plot.caption.position = "plot",
            plot.caption = element_text(hjust = 0))

rich
####q1####

#results in Appendix S2: Table S1
#full model
Mfullq1 <- glmmTMB(q1 ~ (Time) +I(Time^2)+
                     volume_ml+
                     intactprey_number+
                     midge_number+
                     fleshfly_number+
                     mosquito_number+
                     latitude+
                     longitude+
                     (1|site) + (1|site:pitcher_number),
                   family = "gaussian",
                   na.action = na.fail,
                   data = full_data)

#check out full model
summary(Mfullq1)
check_collinearity(Mfullq1)#no colinearity
r.squaredGLMM(Mfullq1)


#DHARMa residuals check
simulationOutput1 <- simulateResiduals(fittedModel = Mfullq1,n=10000, plot=F)
plot(simulationOutput1) #heteroskedastic

#look at each covariate to identify what is driving heteroskedasticity
testQuantiles(simulationOutput1, predictor =full_data$Time)#yes
testQuantiles(simulationOutput1, predictor =full_data$volume_ml)#yes
testQuantiles(simulationOutput1, predictor =full_data$midge_number)#yes
testQuantiles(simulationOutput1, predictor =full_data$mosquito_number)#yes
testQuantiles(simulationOutput1, predictor =full_data$intactprey_number)#yes
testCategorical(simulationOutput1, catPred = full_data$latitude)#yes
testCategorical(simulationOutput1, catPred = full_data$longitude)#yes
testCategorical(simulationOutput1, catPred = full_data$fleshfly_number)#no

#full model w dispersion
Mfullq1_disp <- glmmTMB(q1 ~ (Time) +I(Time^2)+
                          volume_ml+
                          intactprey_number+
                          midge_number+
                          fleshfly_number+
                          mosquito_number+
                          latitude+
                          longitude+
                          (1|site) + (1|site:pitcher_number),
                        dispformula=~latitude+longitude+Time+volume_ml,
                        family = "gaussian",
                        na.action = na.fail,
                        data = full_data)

#check out full model
summary(Mfullq1_disp)
check_collinearity(Mfullq1)#no colinearity
r.squaredGLMM(Mfullq1)


#DHARMa residuals check
simulationOutput1 <- simulateResiduals(fittedModel = Mfullq1_disp,n=10000, plot=F)
plot(simulationOutput1) #good

#model selection
drop1(Mfullq1_disp, test="Chisq")

#best model q1 from drop1
Mbestq1<- glmmTMB(q1 ~ Time+
                    longitude+
                    latitude+
                    volume_ml+
                    fleshfly_number+
                    (1|site) + (1|site:pitcher_number),
                  dispformula=~latitude+volume_ml+longitude+Time,#account for the heteroskedsicity in model
                  family = "gaussian",
                  na.action = na.fail,
                  data = full_data)

#check out full model
summary(Mbestq1)
check_collinearity(Mbestq1)#no colinearity
AIC(Mbest)



#DHARMa residuals check
simulationOutput1 <- simulateResiduals(fittedModel = Mbestq1,n=10000, plot=F)
plot(simulationOutput1) #ok 




####q2####
#results in Appendix S2: Table S2
#full model
Mfullq2 <- glmmTMB(q2 ~ (Time) +I(Time^2)+
                     volume_ml+
                     intactprey_number+
                     midge_number+
                     fleshfly_number+
                     mosquito_number+
                     latitude+
                     longitude+
                     (1|site) + (1|site:pitcher_number),
                   #dispformula=~latitude+longitude+day+volume_ml,
                   family = "gaussian",
                   na.action = na.fail,
                   data = full_data)

#check out full model
summary(Mfullq2)
check_collinearity(Mfullq2)#no colinearity
r.squaredGLMM(Mfullq2)


#DHARMa residuals check

simulationOutput1 <- simulateResiduals(fittedModel = Mfullq2,n=10000, plot=F)
plot(simulationOutput1) #heteroskedastic

#model w/dispersion
Mfullq2_dsip <- glmmTMB(q2 ~ (Time) +I(Time^2)+
                          volume_ml+
                          intactprey_number+
                          midge_number+
                          fleshfly_number+
                          mosquito_number+
                          latitude+
                          longitude+
                          (1|site) + (1|site:pitcher_number),
                        dispformula=~latitude+longitude+Time+volume_ml,
                        family = "gaussian",
                        na.action = na.fail,
                        data = full_data)

#check out full model
summary(Mfullq2_dsip)
check_collinearity(Mfullq2)#no colinearity



#DHARMa residuals check
simulationOutput1 <- simulateResiduals(fittedModel = Mfullq2_dsip,n=10000, plot=F)
plot(simulationOutput1) #better


#model selection
drop1(Mfullq2_dsip, test="Chisq")

#best model q2 from drop1
Mbestq2<- glmmTMB(q2 ~ Time+
                    longitude+
                    latitude+
                    volume_ml+
                    fleshfly_number+
                    (1|site) + (1|site:pitcher_number),
                  dispformula=~longitude+Time+volume_ml+latitude,#account for the heteroskedasicity in model
                  family = "gaussian",
                  na.action = na.fail,
                  data = full_data)

#check out full model
summary(Mbestq2)
check_collinearity(Mbestq2)#no colinearity



#DHARMa residuals check
simulationOutput1 <- simulateResiduals(fittedModel = Mbestq2,n=10000, plot=F)
plot(simulationOutput1) #ok 

####phylogenetic (SES_PD) alpha diversity model####

#full model
Mfullp <- glmmTMB( site_SES_PD~ (Time) +I(Time^2)+
                     volume_ml+
                     intactprey_number+
                     (midge_number)+
                     I(midge_number^2)+
                     fleshfly_number+
                     mosquito_number+
                     latitude+
                     longitude+
                     (1|site) + (1|site:pitcher_number),
                   family = "gaussian",
                   na.action = na.fail,
                   data = full_data)

#check out full model
summary(Mfullp)
check_collinearity(Mfullp)#low/moderate (midge) colinearity


#DHARMa residuals check

simulationOutput1 <- simulateResiduals(fittedModel = Mfullp,n=10000, plot=F)
plot(simulationOutput1)  #heteroskedastic


testQuantiles(simulationOutput1, predictor =full_data$Time)#no
testQuantiles(simulationOutput1, predictor =full_data$volume_ml)#no
testQuantiles(simulationOutput1, predictor =full_data$midge_number)#no
testQuantiles(simulationOutput1, predictor =full_data$mosquito_number)#no
testQuantiles(simulationOutput1, predictor =full_data$intactprey_number)#no
testCategorical(Mfull, catPred = full_data$latitude)#yes
testCategorical(Mfull, catPred = full_data$longitude)#yes
testCategorical(Mfull, catPred = full_data$fleshfly_number)#no


#full model with dispersion
Mfullp_disp <- glmmTMB( site_SES_PD~ (Time) +I(Time^2)+
                          volume_ml+
                          intactprey_number+
                          (midge_number)+
                          I(midge_number^2)+
                          fleshfly_number+
                          mosquito_number+
                          latitude+
                          longitude+
                          (1|site) + (1|site:pitcher_number),
                        dispformula=~latitude+longitude,
                        family = "gaussian",
                        na.action = na.fail,
                        data = full_data)

#check out full modelwith dispersion
summary(Mfullp_disp)
check_collinearity(Mfullp_disp)#low colinearity


simulationOutput1 <- simulateResiduals(fittedModel = Mfullp_disp,n=10000, plot=F)
plot(simulationOutput1)  #looks good


#model selection
mod_sel<-drop1(Mfullp_disp, test = "Chisq") 



#best model is actually fine without the dispersion model (see Dharma residuals below)
Mbestp<- glmmTMB(site_SES_PD~ Time +
                   intactprey_number+
                   midge_number+
                   longitude+(1|site) + 
                   (1|site:pitcher_number),
                 family = "gaussian",
                 na.action = na.fail,
                 data = full_data)

#check out full model
summary(Mbestp)
check_collinearity(Mbestp)#no colinearity
r.squaredGLMM(Mbestp) #R2m= 0.16, R2c=0.25




#DHARMa residuals check
simulationOutput1 <- simulateResiduals(fittedModel = Mbestp,n=10000, plot=F)
plot(simulationOutput1)  #looks fine




####model sketch ses.pd
range(full_data$Time) # -1.289882  2.22
range(full_data$midge_number)#-0.4023389 , 5.9456754
range(full_data$longitude)#-1.186735  1.435534

#Sketch for Time
MyData <- expand.grid(
          longitude = mean(full_data$longitude),
          midge_number=mean(full_data$midge),
          intactprey_number = mean(full_data$intactprey_number),
          Time  = seq(-1.28 , 2.22, length = 16))

X <- model.matrix( ~ Time + intactprey_number +  
                     midge_number+ longitude,
                   data = MyData)


#Calculate ETA and MU
eta <- X%*% fixef(Mbestp)$cond #(use is the best fit model)
mu <- eta 


#Get the standard errors
SE   <- sqrt(diag(X %*%vcov(Mbestp)$cond %*% t(X))  )#hand calculates standard error
seup <- eta + 1.96 * SE
selo <- eta - 1.96 * SE
NewData<-cbind(MyData,eta,mu,SE,seup,selo) 

dayz<-c(0,7,14,21,28,35,42,49,56,63,70,77,84,91,98,105)
day2<-full_data_raw$Time


#Figure 2B
pd_day<-ggplot(NewData, aes(x=dayz, y=mu))+
  geom_point(data =full_data , 
             mapping = aes(x = day2, y = site_SES_PD),color = "gray")+
  #geom_point(data =full_data,
  #mapping = aes(x = day2, y = site_SES_PD), color="gray")+#add raw points
  geom_ribbon(aes(x=dayz,ymin=selo, ymax=seup), fill="#21918c")+
  geom_line(size=1.2, color = "black")+
  geom_hline(yintercept=1.96, color="blue")+
  geom_hline(yintercept=-1.96, color="blue")+
  xlab("Time (days)")+
  ylab("SES PD")+
  ylim(-8, 3)+
  theme_bw()+
  theme(text = element_text(size=15),
        plot.caption.position = "plot",
        plot.caption = element_text(hjust = 0))



#Sketch for Midge
MyData <- expand.grid(
          longitude = mean(full_data$longitude),
          intactprey_number = mean(full_data$intactprey_number),
          volume_ml =mean(full_data$volume_ml),
          Time = mean(full_data$Time),
          midge_number  = seq(-0.40 , 5.95, length = 21))


X <- model.matrix( ~ Time+ intactprey_number +  
                     midge_number+ longitude,
                   data = MyData)


#Calculate ETA and MU
eta <- X%*% fixef(Mbestp)$cond #(use is the best fit model)
mu <- eta 


#Get the standard errors
SE   <- sqrt(diag(X %*%vcov(Mbestp)$cond %*% t(X))  )#hand calculates standard error
seup <- eta + 1.96 * SE
selo <- eta - 1.96 * SE

NewData<-cbind(MyData,eta,mu,SE,seup,selo)  

#raw midge data
midge<-0:20
mm<-full_data_raw$midge_number


#Fig 2C
pd_midge<-ggplot(NewData, aes(x=midge, y=mu))+
          geom_point(data =full_data , 
                     mapping = aes(x = mm, y = site_SES_PD), color="gray")+#add raw points
          geom_ribbon(aes(x=midge,ymin=selo, ymax=seup), fill="#21918c")+
          geom_line(size=1.2, color = "black")+
          geom_hline(yintercept=1.96, color="blue")+
          geom_hline(yintercept=-1.96, color="blue")+
          xlab("Midge Abundance")+
          ylab("SES PD")+
          ylim(-7.5, 3)+
          theme_bw()+
          theme(text = element_text(size=15),
                plot.caption.position = "plot",
                plot.caption = element_text(hjust = 0))


#make panel figure that includes class trajectories from 'class.plots.R'
figure <- ggarrange(rich, pd_day, pd_midge,class,
                    labels = c("A)", "B)", "C)", "D)"),
                    ncol = 2, nrow = 2, common.legend = F)
figure

annotate_figure(figure, bottom = text_grob("Sample Day", color = "black", size = 14))



###############################################

