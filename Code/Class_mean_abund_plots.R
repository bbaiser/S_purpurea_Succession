
#Calculate and plot mean relative abundance of the top ten most abundant bacterial classes
#Fig. 2D

####set up####

#Install and load relevant packages 
library(tidyverse)
library(dplyr)
library(ggplot2)



#data

#the pitcher bu ASV table (columns are ASVs (n=3376), rows are pitchers (n=515))
asv <- read.csv("Data/pitcher_plant_ASV_table.csv", row=1)
dim(asv)

#meta is the covariate data (n=11 varaibles/columns) associated with each pitcher (n=515 rows)
meta <- read.csv("Data/pitcher_plant_meta_table.csv", row=1)
dim(meta)
                         
#Taxonomic list of ASVs (n=3376)
tax <- read.csv("Data/pitcher_plant_tax_table.csv", row=1)

#change row names to a column for down stream         
meta<-meta%>%
      rownames_to_column(var = "pitcher_id") %>%
      rename("Time" = "Actual_Day")


# Plotting the relative abundance at class level
asv16.g <- data.frame(Class=tax$Class,t(asv))%>%
           mutate(Class = case_when(
            Class == "Acidobacteriae"~ "Acidobacteriia",
            Class == "Chthonomonadetes" ~ "Chthonomonadia",
            Class == "Planctomycetes" ~ "Planctomycetia",
            Class=="Verrucomicrobiae" ~ "Verrucomicrobiia",
            TRUE ~ Class ))
  

asv16.g$Class[is.na(asv16.g$Class)] <- "Unknown"
asv16.ga <- aggregate(. ~ asv16.g$Class, asv16.g[,2:ncol(asv16.g)], sum)
row.names(asv16.ga) <- asv16.ga[,1]
asv16.ga <-as.data.frame( asv16.ga[,2:ncol(asv16.ga)])


#TO MAKE RELATIVE ABUNDANCE TABLE
row_sums <- rowSums(asv16.ga)
asv16.ga <- asv16.ga[row_sums != 0, ]
asv16.gao <- as.matrix(asv16.ga[order(rowSums(asv16.ga),decreasing = T),])
asv16.gao1 <- asv16.gao
asv16.gao1_per <- apply(asv16.gao1, 2, function(x){x/sum(x,na.rm=T)})#calculate relative abundance in each sample



#format as time series for each class
asv_class <-as.data.frame( t(asv16.ga[,2:ncol(asv16.ga)]))%>% 
            rownames_to_column( var = "pitcher_id")%>% 
            left_join(meta, by="pitcher_id")%>% 
            select(!c(volume_ml, latitude,longitude, mosquito_number))%>% 
            group_by(Time)%>% 
            summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))%>% 
            mutate(Time = as.numeric(Time))


#10 most common classes (highest mean abundances)
top_ten_class<-as.data.frame(colMeans(asv_class))%>% 
               arrange(desc(.))%>%
               filter(row_number() %in% c(1,2,3,4,6,7,8,9,10,11))#this removes the means for time (i.e., row 5)


#make into long format to plot
long <- asv_class %>% 
       pivot_longer(
       cols = rownames(top_ten_class) , 
       names_to = "Class",
       values_to = "Abundance")
                              
#plot mean relative abundances of top ten ASVs (Figure 2D) get integrated into panel with alpha diversity plots
class <- ggplot(data=long, aes(x=Time, y=Abundance, color=Class)) + 
         geom_smooth( method = "loess", se = FALSE)+
         xlab("Time (days)")+
         ylab("Mean Relative Abundance")+
         theme_bw()+
         theme(text = element_text(size=15),
            plot.caption.position = "plot",
            plot.caption = element_text(hjust = 0))+
         guides(shape = guide_legend(override.aes = list(size = .15)))+
         guides(color = guide_legend(override.aes = list(size = .15)))+
         theme(legend.title = element_text(size = 10), 
                   legend.text = element_text(size = 8))
  

