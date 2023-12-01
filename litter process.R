## LEAF LITTER CLEANING ##

library(dplyr)

# Bring in raw data files from Jenn Follstad Shah and Carri LeRoy
# File can also be found here: https://raw.githubusercontent.com/andrew-hipp/decomposition-phylogeny-2019/master/data/k-values_with_select_env_variables_CJL_ALH_9-26-18.csv

leaf1 <- read.csv(file="LeRoy.ExpandedDataset.Kvalues.csv",na.strings = ".")
leaf_cond <- read.csv(file="LeafConditionKey2.csv")

leaf2 <- merge(leaf1,leaf_cond,by='Sorting.code')

# Fix the extra space in "green "
leaf2[leaf2$Leaf.condition=="green ","Leaf.condition"] <- "green"

# Simplify to just needed columns
litter_var <- c('Long_DD','Lat_DD',
                'Leaf.condition','mesh.size.category',
                'Genus','k..d.','mean..ºC.','Citation','Sorting.code')

# Rename columns
leaf3 <- leaf2[,litter_var]
colnames(leaf3) <- c('longitude','latitude','Leaf.condition','Mesh.size',
                     'Genus','kd','mean_temp','Citation','Sorting.code')

# Retain only senesced litter with coarse or fine bags
table(leaf3$Leaf.condition)
table(leaf3$Mesh.size)

# Air-dried was determined to be same as senesced
leaf3[leaf3$Leaf.condition=="air-dried",'Leaf.condition'] <- "senesced"

leaf4 <- leaf3[leaf3$Leaf.condition %in% c("senesced","green"),]
leaf5 <- leaf4[leaf4$Mesh.size %in% c("fine","coarse"),]

# Calculate mean kd and temp by location, litter type, and mesh size
# Retain sorting code to recapture citations
leaf6 <- leaf5 %>% group_by(Mesh.size,Leaf.condition,Genus,latitude,longitude) %>% 
  summarise(mean_kd=mean(kd),Sorting.code=min(Sorting.code),mean_mean_daily_temp=mean(mean_temp,na.rm=T),.groups = 'drop') %>%
  as.data.frame()

dim(leaf6) #1005 measurements of kd

# Merge back citations
leaf7 <- merge(leaf6,leaf5[,c("Sorting.code","Citation")],by='Sorting.code')

# Remove any litter genera where n < 3
gen_tab <- as.data.frame(table(leaf7$Genus))
dim(gen_tab) #105 total genera
gen_tab3 <- gen_tab[gen_tab$Freq>3,] #Only genera with n > 3
dim(gen_tab3) #35 genera

leaf8 <- leaf7[leaf7$Genus %in% gen_tab3$Var1,] 
dim(leaf8) #895 measurements

# Write the final cleaned file for BRT analysis
#write.csv(leaf8,"litter_processed.csv",row.names = F)


## TRAIT CLEANING ##

# Bring in raw data file
# From lit review done by Jenn Follstad Shah
litter_trait <- read.csv(file="Litter_trait_review.csv")
summary(litter_trait)

# How many genera
length(unique(litter_trait$Genus)) #172 genera

# Aggregate by genus-level means
litter_Mn <- litter_trait %>% group_by(Genus) %>% 
  summarise(N_Litter_Mn=mean(perN,na.rm=T),
            P_Litter_Mn=mean(perP,na.rm=T),
            C_Litter_Mn=mean(perC,na.rm=T),
            Lignin_Litter_Mn=mean(perLignin,na.rm=T),
            Cellulose_Litter_Mn=mean(perCellulose,na.rm=T),
            CtoN_Litter_Mn=mean(CtoN,na.rm=T),
            NtoP_Litter_Mn=mean(NtoP,na.rm=T),
            .groups = 'drop') %>%
  as.data.frame()

# Bring in leaf traits from TRY database
leaf_TRY <- read.csv("TRY_traits.csv")

#Look for mismatches
leaf_TRY[!(leaf_TRY$Genus %in% litter_Mn$Genus),'Genus'] 
  #Rubus and Carex are in TRY but not the litter review

gen_tab3[!(gen_tab3$Var1 %in% leaf_TRY$Genus),'Var1']
  #8 genera with >3 decomp rates but no leaf traits

gen_tab3[!(gen_tab3$Var1 %in% litter_Mn$Genus),'Var1']
  #All genera have at least 1 litter trait


# Merge leaf and litter traits
traits <- merge(litter_Mn,leaf_TRY,by='Genus',all.x=T,all.y=T)

# Write database of traits
#write.csv(traits,file="traits.csv")
