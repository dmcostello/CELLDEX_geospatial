# Packages required
library(sp)
library(raster)
library(geodata)
library(fields)
library(sf)
library(caret)
library(gbm)
library(pdp)
library(gridExtra)
library(grid)

# Set random start to reproduce published model iteration
set.seed(15)

###########################################
#### Retrieve CELLDEX data from GitHub ####
###########################################

# SOURCE: Tiegs et al. 2019 https://doi.org/10.1126/sciadv.aav0486
# DATA REPOSITORY: https://github.com/dmcostello/CELLDEX2018 
# DATA DOI: https://doi.org/10.5281/zenodo.2591572

# Site characteristics
CELLDEXsites <- read.csv("https://raw.githubusercontent.com/dmcostello/CELLDEX2018/edddf5a9477585cbbcb8280e686fc8a6000a29e8/CELLDEX_SITE_DATA.csv")
CELLDEXsites$part.str <- paste(CELLDEXsites$partnerid,CELLDEXsites$stream)

# Temperature during deployment
tempdata <- read.csv("https://raw.githubusercontent.com/dmcostello/CELLDEX2018/edddf5a9477585cbbcb8280e686fc8a6000a29e8/CELLDEX_TEMPERATURE.csv")
tempdata$part.str <- paste(tempdata$partnerid,tempdata$stream)
tempdata2 <- aggregate(tempdata[,'mean_mean_daily_temp'],list(tempdata$part.str,tempdata$habitat),mean)
colnames(tempdata2)[1:3] <- c("part.str","habitat","mean_mean_daily_temp")
strtemp <- subset(tempdata2,habitat=="STR")

# Cotton decay rate
CELLDEXkd <- read.csv("https://raw.githubusercontent.com/dmcostello/CELLDEX2018/edddf5a9477585cbbcb8280e686fc8a6000a29e8/str_k.csv")
CELLDEXkd$logk <- log(CELLDEXkd$k)

# bring in deployment dates by CELLDEX location
dates<-read.csv("CELLDEX_deploy_dates.csv")
dates$part.str <- paste(dates$partnerid,dates$stream)
dates$deploy_date <- as.Date(dates$deploy_date,format="%m/%d/%y")
dates$month <- format(dates$deploy_date,"%m")

#Merge CELLDEX datasets
mer1 <- merge(CELLDEXkd,strtemp[,c('part.str','mean_mean_daily_temp')],by='part.str')
mer2 <- merge(mer1,dates[,c('deploy_date','part.str','month')],by='part.str',all.y=F)
CELLDEX <- merge(mer2,CELLDEXsites[,c('latitude','longitude','part.str')],
                 by='part.str',all.y=F)

# Turn CELLDEX stream locations and data into a spatial points data frame
Cpts<-SpatialPointsDataFrame(CELLDEX[,c('longitude','latitude')],CELLDEX,
                             proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))


######################################################
#### Rasters of nutrients (stream and deposition) ####
######################################################

# Use WorldClim data to set raster extent and resolution
clim <- worldclim_global("tmin",res=10,path=tempdir())
r <- raster::raster(clim$wc2.1_10m_tmin_12)

# Bring in stream concentration rasters
# SOURCE: McDowell et al. (2021) https://doi.org/10.1002/gdj3.111
# DATA REPOSITORY: https://doi.org/10.25400/lincolnuninz.11894697 (rasters provided by corr. author)
# UNITS: kg (NO3-N or DRP-P)/ha/yr

NO3<-raster::raster("NO3_Conc.tif")
DRP<-raster::raster("DRP_Conc.tif")

# Bring in N deposition and create raster
# SOURCE: Ackerman et al. (2019) https://doi.org/10.1029/2018GB005990
# DATA REPOSITORY: https://doi.org/10.13020/D6KX2R
# UNITS: kg N/km2/yr

Ndep <- read.csv("https://conservancy.umn.edu/bitstream/handle/11299/197613/inorganic_N_deposition.csv?sequence=28&isAllowed=y")

# turn N deposition locations and data into a spatial points data frame
Npts<-SpatialPointsDataFrame(Ndep[,c('longitude','latitude')],Ndep)
crs(Npts)<-crs(Cpts) # set coordinate system

# create raster of total N deposition with WorldClim raster as template
r2<-raster(extent(r),nrows=91,ncols=144,crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
TNdep<-raster::rasterize(Npts,r2,field="tot_2016")

# Bring in P deposition and raster
# SOURCE 1: Brahney et al. (2015) https://doi.org/10.1002/2015GB005137
# DATA REPOSITORY 1: No repository, downloaded from supporting information (Table S1)
# UNITS 1: kg P/km2/yr
# SOURCE 2: Mahowald et al. (2008)  https://doi.org/10.1029/2008GB003240
# DATA REPOSITORY 2: No repository, downloaded from supporting information and converted units (Table S2b)
# UNITS 2: kg P/km2/yr

dat1<-read.csv("gbc20320-sup-0002-supinfo.csv")
dat2<-read.csv("gbc1549-sup-0009-ts02.csv")
colnames(dat2)[1:2]<-c("Latitude","Longitude") #Rename columns to match dat1 and dat2
TPdepdat<-rbind(dat1[,c(3,4,12)],dat2[,c(1,2,4)])

# Interpolate a raster of total P deposition from TPdepdat using Worldclim raster as template
TPdep <- Tps(TPdepdat[,c("Longitude","Latitude")], TPdepdat[,"TP.mg.m.2.yr.1"],lon.lat=TRUE)
TPdep_interp <- raster::interpolate(r, TPdep)


###############################
#### Load HydroBASINS data ####
###############################

# The HydroBASINS shape file is large (>1 GB) and not uploaded in this repository
# These data can be downloaded at  https://www.hydrosheds.org/products/hydroatlas

# SOURCE: Linke et al. (2019) https://doi.org/10.1038/s41597-019-0300-6
# DATA REPOSITORY: https://www.hydrosheds.org/products/hydroatlas
# RESOLUTION: 12-digit subwatersheds

basins<-st_read("BasinATLAS_v10_lev12.shp")

# Subset subwatersheds to get only those with >1 area in ha
basins<-subset(basins,basins$ria_ha_ssu>1)

# Transform subwatershed polygons to equal distance Mollenweide projection
bas_moll <- st_transform(basins, "+proj=moll")
rm(basins)


#################################
#### Join and clean datasets ####
#################################

# Join nutrient data at CELLDEX locations
nut_data <- list(NO3c=NO3,DRPc=DRP,TNdep=TNdep,TPdep=TPdep_interp)
merge_temp <- data.frame(Cpts@data$part.str)
for(i in 1:length(nut_data)){
  tempdat <- raster::extract(nut_data[[i]],Cpts,df=T)
  merge_temp <- cbind(merge_temp,tempdat[,2])
  colnames(merge_temp)[1+i] <- names(nut_data)[i]
}

Cpts@data<-cbind(Cpts@data,merge_temp[,2:5])

# Convert CELLDEX points to sf object
C_st<-as(Cpts,"sf")

# Transform CELLDEX spatial points to Mollweide
C_moll<-st_transform(C_st,"+proj=moll")

# Spatially join CELLDEX stream points to attributes of nearest subwatershed
C_stb<-st_join(C_moll, bas_moll, join = st_nearest_feature)
C_stb$geometry<-NULL

# Turn -999 values (no data as represented in Hydrobasins) to NA
C_stb[C_stb==-999]<-NA

# Pull out climate data from Hydrobasins subwatershed for the month of cotton deployment
C_stb$tempmonth<-NULL
C_stb$precipmonth<-NULL
C_stb$PETmonth<-NULL
C_stb$AETmonth<-NULL
C_stb$moist_indexmonth<-NULL
C_stb$snowmonth<-NULL
C_stb$soilwatermonth<-NULL
for(i in 1:nrow(C_stb)) {
  month<-C_stb$month[i]
  if (!is.na(month)){
    C_stb$tempmonth[i]<-C_stb[i,paste0("tmp_dc_s",month)]
    C_stb$precipmonth[i]<-C_stb[i,paste0("pre_mm_s",month)]
    C_stb$PETmonth[i]<-C_stb[i,paste0("pet_mm_s",month)]
    C_stb$AETmonth[i]<-C_stb[i,paste0("aet_mm_s",month)]
    C_stb$moist_indexmonth[i]<-C_stb[i,paste0("cmi_ix_s",month)]
    C_stb$snowmonth[i]<-C_stb[i,paste0("snw_pc_s",month)]
    C_stb$soilwatermonth[i]<-C_stb[i,paste0("swc_pc_s",month)]
  }
}

# Build subset including only variables for the model
# This excludes (1) all monthly data, (2) all class variables,(3) variables with 
# no variance, and (4) composite indices (i.e., human development and biome).
# Simultaneously, log transform variables that have large skew and back-transform 
# 10x and 100x corrections done in HydroBASINS
# Select variables needing transformation or back-transformation are identified
# in the variable names datasheet.
mod_vars <- read.csv(file="var_names.csv")

#Setup empty dataframe
Cdat <- data.frame(logk=C_stb$logk)

#For loop to select variables and complete transformations 
for(i in 2:length(mod_vars$Variables)){
  if(mod_vars$Transform[i]=="none"){
    tempname <- mod_vars$Variables[i]
    Cdat <- cbind(Cdat,C_stb[,mod_vars$Variables[i]])
    colnames(Cdat)[i] = tempname
  } else
    if(mod_vars$Transform[i]=="log"){
      tempname <- paste0("log10",mod_vars$Variables[i])
      Cdat <- cbind(Cdat,log10(C_stb[,mod_vars$Variables[i]]))
      colnames(Cdat)[i] = tempname
    } else
      if(mod_vars$Transform[i]=="log1"){
        tempname <- paste0("log10",mod_vars$Variables[i])
        Cdat <- cbind(Cdat,log10(C_stb[,mod_vars$Variables[i]]+1))
        colnames(Cdat)[i] = tempname
      } else
        if(mod_vars$Transform[i]=="xten"){
          tempname <- mod_vars$Variables[i]
          Cdat <- cbind(Cdat,(C_stb[,mod_vars$Variables[i]]/10))
          colnames(Cdat)[i] = tempname
        } else
          if(mod_vars$Transform[i]=="xhund"){
            tempname <- mod_vars$Variables[i]
            Cdat <- cbind(Cdat,(C_stb[,mod_vars$Variables[i]]/100))
            colnames(Cdat)[i] = tempname
          } 
}

dim(Cdat) #102 total variables


######################################
#### BRT to predict cotton decomp ####
######################################

# Run boosted regression tree model with Gaussian error (e.g. linear regression)
Cgbm<- gbm(logk~., 
           data=Cdat, 
           distribution="gaussian",
           n.trees=20000,
           shrinkage=0.001,
           interaction.depth=5,
           cv.folds=20)

# Check performance
(best.iter <- gbm.perf(Cgbm,method="cv"))

# Plot variable influence based on the estimated best number of trees
sum<-summary(Cgbm,n.trees=best.iter,method=permutation.test.gbm) 

head(sum,20)

# Calculate pseudo-R2
print(1-sum((Cdat$logk - predict(Cgbm, newdata=Cdat, n.trees=best.iter))^2)/
            sum((Cdat$logk - mean(Cdat$logk))^2))

plot(Cgbm,i.var='tmp_dc_smx')


###########################################
#### Figure 1 partial dependence plots ####
###########################################

# Extract partial dependence values for top importance variables
drp<-plot(Cgbm,i.var='log10DRPc',return.grid=TRUE)
pet<-plot(Cgbm,i.var='pet_mm_uyr',return.grid=TRUE)
meantemp<-plot(Cgbm,i.var='mean_mean_daily_temp',return.grid=TRUE)
limn<-plot(Cgbm,i.var='log10lka_pc_sse',return.grid=TRUE)
NO3<-plot(Cgbm,i.var='log10NO3c',return.grid=TRUE)
mntmp<-plot(Cgbm,i.var='tmp_dc_smn',return.grid=TRUE)

# Create individual panels
ps1<-ggplot(data = drp, aes(log10DRPc, y)) +
  geom_line(color = "steelblue", size = 1) +
  ylab("") + 
  xlab(expression(~log[10]~"Stream dissolved reactive P kg ha"^"-1"~"yr"^"-1")) +
  geom_rug(sides="b",position = "jitter",color="cornflowerblue") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

ps2<-ggplot(data = pet, aes(pet_mm_uyr, y)) +
  geom_line(color = "steelblue", size = 1) +
  ylab("") + 
  labs(x=bquote('Upstream mean PET mm yr'^'-1')) +
  geom_rug(sides="b",position = "jitter",color="cornflowerblue") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

ps3<-ggplot(data = meantemp, aes(mean_mean_daily_temp, y)) +
  geom_line(color = "steelblue", size = 1) +
  ylab("") + 
  labs(x="Mean daily stream temp C during deployment") +
  geom_rug(sides="b",position = "jitter",color="cornflowerblue") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

ps4<-ggplot(data = limn, aes(log10lka_pc_sse, y)) +
  geom_line(color = "steelblue", size = 1) +
  ylab("") + labs(x=~log[10]~"Subwatershed lake area %") +
  geom_rug(sides="b",position = "jitter",color="cornflowerblue") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

ps5<-ggplot(data = NO3, aes(log10NO3c, y)) +
  geom_line(color = "steelblue", size = 1) +
  ylab("") + 
  xlab(expression(~log[10]~"Stream NO"["2"]~"- NO"["3"]~"kg ha"^"-1"~"yr"^"-1")) +
  geom_rug(sides="b",position = "jitter",color="cornflowerblue") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

ps6<-ggplot(data = mntmp, aes(tmp_dc_smn, y)) +
  geom_line(color = "steelblue", size = 1) +
  ylab("") + 
  xlab("Subwatershed minimum temperature C") +
  geom_rug(sides="b",position = "jitter",color="cornflowerblue") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Make trellis of partial dependence plots for top 6 important variables
pdf(file = "stream_pdps_6.pdf",width=4,height=10)
grid.arrange(ps1,ps2,ps3,ps4,ps5,ps6, ncol = 1,left=textGrob(bquote('ln' ~K[d]),rot=90))
dev.off()

########################################
#### Load and join leaf litter data ####
########################################

# SOURCE: Follstad Shah et al. (2017) https://doi.org/10.1111/gcb.13609
# DATA REPOSITORY: None. Received updated files from authors
# DATA PROCESSING: (1) Calculated mean decay rates for sites with duplicate measures.
# (2) Removed genera with <3 measures of decomposition rate.

fs <- read.csv("FS2017_expanded.csv")
fs$logk <- log(fs$kd)
fs$mean_mean_daily_temp <- NA #### CHANGE WHEN MEAN DAILY TEMP IS ADDED TO FS DATA

# Turn litter data into spatial points data frame
fs_pts<-SpatialPointsDataFrame(fs[,c("Longitude","Latitude")],fs,
                               proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# Join nutrient data at litter data locations
merge_temp2 <- data.frame(fs_pts@data$Genus)
for(i in 1:length(nut_data)){
  tempdat <- raster::extract(nut_data[[i]],fs_pts,df=T)
  merge_temp2 <- cbind(merge_temp2,tempdat[,2])
  colnames(merge_temp2)[1+i] <- names(nut_data)[i]
}

fs_pts@data<-cbind(fs_pts@data,merge_temp2[,2:5])

# Convert FS points to sf object
fs_st<-as(fs_pts,"sf")

# Transform FS spatial points to Mollweide
fs_moll<-st_transform(fs_st,"+proj=moll")

# Spatially join FS stream points to attributes of nearest subwatershed
fs_stb<-st_join(fs_moll, bas_moll, join = st_nearest_feature)
fs_stb$geometry<-NULL

# Turn -999 values (no data as represented in Hydrobasins) to NA
fs_stb[fs_stb==-999]<-NA

# Litter database does not include time of deployment
# All month of deployment variables are NA
fs_stb$tempmonth<-NA
fs_stb$precipmonth<-NA
fs_stb$PETmonth<-NA
fs_stb$AETmonth<-NA
fs_stb$moist_indexmonth<-NA
fs_stb$snowmonth<-NA
fs_stb$soilwatermonth<-NA

# Build subset including predictor variables used for CELLDEX model
# Simultaneously, log transform variables that have large skew and back-transform 
# 10x and 100x corrections done in HydroBASINS

#Setup empty dataframe
Ldat <- data.frame(logk=fs_stb$logk)

#For loop to select variables and complete transformations 
for(i in 2:length(mod_vars$Variables)){
  if(mod_vars$Transform[i]=="none"){
    tempname <- mod_vars$Variables[i]
    Ldat <- cbind(Ldat,fs_stb[,mod_vars$Variables[i]])
    colnames(Ldat)[i] = tempname
  } else
    if(mod_vars$Transform[i]=="log"){
      tempname <- paste0("log10",mod_vars$Variables[i])
      Ldat <- cbind(Ldat,log10(fs_stb[,mod_vars$Variables[i]]))
      colnames(Ldat)[i] = tempname
    } else
      if(mod_vars$Transform[i]=="log1"){
        tempname <- paste0("log10",mod_vars$Variables[i])
        Ldat <- cbind(Ldat,log10(fs_stb[,mod_vars$Variables[i]]+1))
        colnames(Ldat)[i] = tempname
      } else
        if(mod_vars$Transform[i]=="xten"){
          tempname <- mod_vars$Variables[i]
          Ldat <- cbind(Ldat,(fs_stb[,mod_vars$Variables[i]]/10))
          colnames(Ldat)[i] = tempname
        } else
          if(mod_vars$Transform[i]=="xhund"){
            tempname <- mod_vars$Variables[i]
            Ldat <- cbind(Ldat,(fs_stb[,mod_vars$Variables[i]]/100))
            colnames(Ldat)[i] = tempname
          } 
}

dim(Ldat) #102 total variables

# Predict ln(k) for FS sites using stream ln(k) model (Cgbm) above
fs$ln_pred_k<-predict(Cgbm, newdata=Ldat, n.trees=best.iter)

# Combine ln(k) predictions with litter condition and and substrate quality
# in the form of genus level traits (mean and median) from TRY (Kattge et al., 2011)
traits <- read.csv("TRY_traits.csv")

fs2 <- merge(fs,traits,by='Genus')

# Clean up variables and make factors
fs2$mesh.size.category<-factor(fs2$mesh.size.category)
fs2$Leaf.condition[which(fs2$Leaf.condition=="unknown")]<-NA
fs2$Leaf.condition[which(fs2$Leaf.condition=="air-dried")]<-NA
fs2$Leaf.condition[which(fs2$Leaf.condition=="mixed")]<-NA
fs2$Leaf.condition<-factor(fs2$Leaf.condition)

# Build the dataframe containing only variables for the model
val_var <- read.csv(file="validation_variables.csv")
Fdat <- fs2[,colnames(fs2) %in% val_var$Variables]

# Run the BRT validation model
Fgbm<- gbm(logk~., 
           data=Fdat, 
           distribution="gaussian",
           n.trees=50000,
           shrinkage=0.001,
           interaction.depth=5,
           cv.folds=20)

# Check performance
(best.iter2 <- gbm.perf(Fgbm,method="cv"))

# Plot variable influence based on the estimated best number of trees
Fsum<-summary(Fgbm,n.trees=best.iter2,method=permutation.test.gbm) 
Fsum

# Calculate pseudo-R2
print(1-sum((Fdat$logk - predict(Fgbm, newdata=Fdat, n.trees=best.iter2))^2)/
        sum((Fdat$logk - mean(Fdat$logk))^2))

plot(Fgbm,i.var='ln_pred_k')
plot(Fgbm,i.var='mesh.size.category')
plot(Fgbm,i.var='Leaf.condition')

toc()