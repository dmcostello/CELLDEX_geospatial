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
litter_var <- c('Longitude.2','Latitude.2',
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
dim(gen_tab[gen_tab$Freq>3,]) #35 genera with n > 3

leaf8 <- leaf7[leaf7$Genus %in% gen_tab[gen_tab$Freq>3,'Var1'],] 
dim(leaf8) #895 measurements

# Write the final cleaned file for BRT analysis
write.csv(leaf8,"litter_processed.csv",row.names = F)
