# packages required
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
library(fasterize)
library(glwdr)
library(parallel)

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
#clim <- worldclim_global("tmin",res=10,path=tempdir())
#r <- raster::raster(clim$wc2.1_10m_tmin_12)

r<-raster::raster("wc2.1_10m_tmin_12.tif") # from https://www.worldclim.org/


# Bring in stream concentration rasters
# SOURCE: McDowell et al. (2021) https://doi.org/10.1002/gdj3.111
# DATA REPOSITORY: https://doi.org/10.25400/lincolnuninz.11894697 (rasters provided by corr. author)
# UNITS: kg (NO3-N or DRP-P)/ha/yr

NO3<-raster::raster("NO3_Conc.tif")
DRP<-raster::raster("DRP_Conc.tif")

# Bring in N deposition and create raster
# SOURCE: Ackerman et al. (2019) https://doi.org/10.1029/2018GB005990
# DATA REPOSITORY: https://hdl.handle.net/11299/197613
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

# Subset watersheds to get only those with >1 ha of upstream river area
basins1<-subset(basins,basins$ria_ha_ssu>1)

# Transform subwatershed polygons to equal distance Mollenweide projection
bas_moll <- st_transform(basins1, "+proj=moll")


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

#######################################
#### Figure 1 - Cotton BRT results ####
#######################################

# Get grids of x variable values and predictions to make partial dependence plots
drp<-plot(Cgbm,i.var='log10DRPc',return.grid=TRUE)
pet<-plot(Cgbm,i.var='pet_mm_uyr',return.grid=TRUE)
meantemp<-plot(Cgbm,i.var='mean_mean_daily_temp',return.grid=TRUE)
limn<-plot(Cgbm,i.var='log10lka_pc_sse',return.grid=TRUE)
no3<-plot(Cgbm,i.var='log10NO3c',return.grid=TRUE)
mntmp<-plot(Cgbm,i.var='tmp_dc_smn',return.grid=TRUE)

# Get quartiles of the x variable to mark in plots
drprug<-as.data.frame(quantile(Cdat$log10DRPc,probs = seq(0, 1, 0.25),na.rm=T))
colnames(drprug)[1]<-"rug"
petrug<-as.data.frame(quantile(Cdat$pet_mm_uyr),probs = seq(0, 1, 0.25))
colnames(petrug)[1]<-"rug"
mtrug<-as.data.frame(quantile(Cdat$mean_mean_daily_temp),probs = seq(0, 1, 0.25))
colnames(mtrug)[1]<-"rug"
limnrug<-as.data.frame(quantile(Cdat$log10lka_pc_sse,probs=seq(0, 1, 0.25)))
colnames(limnrug)[1]<-"rug"
limnrug[2,1]<-0.01
no3rug<-as.data.frame(quantile(Cdat$log10NO3c,probs=seq(0, 1, 0.25),na.rm=T))
colnames(no3rug)[1]<-"rug"
mntrug<-as.data.frame(quantile(Cdat$tmp_dc_smn,probs=seq(0, 1, 0.25),na.rm=T))
colnames(mntrug)[1]<-"rug"

# Make partial dependence plots in ggplot
ps1<-ggplot(data = drp, aes(log10DRPc, y)) +
  geom_point(aes(x=rug,y=min(drp$y)),drprug,color="gray",shape = 15) +
  geom_line(color = "steelblue", size = 1) +
  ylab("") + 
  xlab(expression(~log[10]~"Stream dissolved reactive P kg ha"^"-1"~"yr"^"-1")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

ps2<-ggplot(data = pet, aes(pet_mm_uyr, y)) +
  geom_point(aes(x=rug,y=min(pet$y)),petrug,color="gray",shape = 15) +
  geom_line(color = "steelblue", size = 1) +
  ylab("") + 
  labs(x=bquote('Upstream mean PET mm yr'^'-1')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

ps3<-ggplot(data = meantemp, aes(mean_mean_daily_temp, y)) +
  geom_point(aes(x=rug,y=min(meantemp$y)),mtrug,color="gray",shape = 15) +
  geom_line(color = "steelblue", size = 1) +
  ylab("") + 
  xlab(expression(x="Mean daily stream temp during deployment "*degree*C)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

ps4<-ggplot(data = limn, aes(log10lka_pc_sse, y)) +
  geom_point(aes(x=rug,y=min(limn$y)),limnrug,color="gray",shape = 15) +
  geom_line(color = "steelblue", size = 1) +
  ylab("") + labs(x=~log[10]~"Subwatershed lake area %") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

ps5<-ggplot(data = no3, aes(log10NO3c, y)) +
  geom_point(aes(x=rug,y=min(no3$y)),no3rug,color="gray",shape = 15) +
  geom_line(color = "steelblue", size = 1) +
  ylab("") + 
  xlab(expression(~log[10]~"Stream NO"["2"]~"- NO"["3"]~"kg ha"^"-1"~"yr"^"-1")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

ps6<-ggplot(data = mntmp, aes(tmp_dc_smn, y)) +
  geom_point(aes(x=rug,y=min(mntmp$y)),mntrug,color="gray",shape = 15) +
  geom_line(color = "steelblue", size = 1) +
  ylab("") + 
  xlab(expression("Subwatershed minimum temperature "*degree*C)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Make trellis of the top 6 partial dependence plots
pdf(file = "stream_pdps_6.pdf",width=4,height=10)
grid.arrange(ps1,ps2,ps3,ps4,ps5,ps6, ncol = 1,left=textGrob(bquote('ln' ~K[d]),rot=90))
dev.off()


##################################################
#### Generate global predictions of cotton kd ####
##################################################

# Extract values from Hydrobasins
# Rasterize basins with WorldClim raster as template using fasterize package
gr<-fasterize(basins1,r,field = "HYBAS_ID",fun="max")

# Load large lakes from glwdr to mask out lake pixels
lakes <- glwd_load(level = 1)
crs(lakes)<-crs(Cpts)
gr<-mask(gr,lakes,inverse=T)
rm(lakes)

# Convert raster object to dataframe 
grv<-as.data.frame(gr,xy=T)
colnames(grv)<-c("longitude","latitude","HYBAS_ID")

# Resample nutrient rasters to correct resolution
NO3res <- raster::resample(NO3,r)
DRPres <- raster::resample(DRP,r)
TNdepres <- raster::resample(TNdep,r)

# Stack and merge nutrient rasters
nstack<-raster::stack(NO3res,DRPres,TNdepres,TPdep_interp)
v<-as.data.frame(nstack)
colnames(v)<-c("NO3c","DRPc","TNdep","TPdep")

grv2 <- cbind(grv,v)
grv2$Sort_code <- seq(1,dim(grv2)[1])

# Remove monthly variables for global prediction
# Variables were not among most important in BRT
mod_vars_nomo <- mod_vars[!mod_vars$Variables %in% c("tempmonth",
                                                     "precipmonth",
                                                     "PETmonth",
                                                     "AETmonth",
                                                     "moist_indexmonth",
                                                     "snowmonth",
                                                     "soilwatermonth")&
                            mod_vars$Source=="HydroBASINS",]

# Build dataframe from HydroBASINS
basin_dat <- as.data.frame(basins1)

# Build HydroBASINS data and complete transformations 
for(i in 1:length(mod_vars_nomo$Variables)){
  if(mod_vars_nomo$Transform[i]=="none"){
    tempname <- mod_vars_nomo$Variables[i]
    tempdat <- cbind(basin_dat$HYBAS_ID,basin_dat[,mod_vars_nomo$Variables[i]])
    grv2 <- merge(grv2,tempdat,all.x=T,by.x="HYBAS_ID",by.y=1)
    colnames(grv2)[8+i] = tempname
  } else
    if(mod_vars_nomo$Transform[i]=="log"){
      tempname <- paste0("log10",mod_vars_nomo$Variables[i])
      tempdat <- cbind(basin_dat$HYBAS_ID,log10(basin_dat[,mod_vars_nomo$Variables[i]]))
      grv2 <- merge(grv2,tempdat,all.x=T,by.x="HYBAS_ID",by.y=1)
      colnames(grv2)[8+i] = tempname
    } else
      if(mod_vars_nomo$Transform[i]=="log1"){
        tempname <- paste0("log10",mod_vars_nomo$Variables[i])
        tempdat <- cbind(basin_dat$HYBAS_ID,log10(basin_dat[,mod_vars_nomo$Variables[i]]+1))
        grv2 <- merge(grv2,tempdat,all.x=T,by.x="HYBAS_ID",by.y=1)
        colnames(grv2)[8+i] = tempname
      } else
        if(mod_vars_nomo$Transform[i]=="xten"){
          tempname <- mod_vars_nomo$Variables[i]
          tempdat <- cbind(basin_dat$HYBAS_ID,basin_dat[,mod_vars_nomo$Variables[i]]/10)
          grv2 <- merge(grv2,tempdat,all.x=T,by.x="HYBAS_ID",by.y=1)
          colnames(grv2)[8+i] = tempname
        } else
          if(mod_vars_nomo$Transform[i]=="xhund"){
          tempname <- mod_vars_nomo$Variables[i]
        tempdat <- cbind(basin_dat$HYBAS_ID,basin_dat[,mod_vars_nomo$Variables[i]]/100)
        grv2 <- merge(grv2,tempdat,all.x=T,by.x="HYBAS_ID",by.y=1)
        colnames(grv2)[8+i] = tempname
        } 
}

# Resort to original position
grv2 <- grv2[order(grv2$Sort_code),]

# Monthly averages and mean_mean_daily_temp as NA
grv2$tempmonth <- NA
grv2$precipmonth <- NA
grv2$PETmonth <- NA
grv2$AETmonth <- NA
grv2$moist_indexmonth <- NA
grv2$snowmonth <- NA
grv2$soilwatermonth <- NA
grv2$mean_mean_daily_temp <- NA

# Log transform nutrients
grv2$log10NO3c <- log10(grv2$NO3c)
grv2$log10DRPc <- log10(grv2$DRPc)
grv2$log10TNdep <- log10(grv2$TNdep)
grv2$log10TPdep <- log10(grv2$TPdep)

# Make predictions of kd
grv2$gl_kd <- predict(Cgbm,newdata=grv2,n.trees=best.iter)

# Produce the global raster image if cotton ln(kd)
global_kd <- setValues(gr,grv2$gl_kd)
global_kd<-mask(global_kd,gr)
plot(global_kd)

writeRaster(global_kd,"global_stream_lnkd.tif",overwrite=T)

###############################################
#### Figure 2 - Global raster of cotton kd ####
###############################################

# Project raster in Mollenweide
glkd_moll<-projectRaster(global_kd,crs="+proj=moll")

mask_glkd_moll <- as.data.frame(glkd_moll,xy=T)

mask_glkd <- as.data.frame(global_kd,xy=T)

# make map figures in ggplot
map<-ggplot() + borders(fill="lightgray") + geom_raster(data = mask_glkd , aes(x = x, y = y,fill = layer)) + 
  scale_fill_gradientn(colors=rev(c("lightgray","darkred", "red", "orange", "yellow","darkgreen","green", "lightgreen", "lightblue","blue","violet","lightgray")),na.value=NA,name=bquote('ln Mean Annual Stream' ~K[d]))+
  xlab("") + ylab("") + theme(legend.position = "none",legend.box.background = element_blank(),legend.background = element_blank()) + theme(panel.background = element_rect(fill = "white",colour = "white",size = 1, linetype = "solid")) +
  theme(axis.text.x=element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.ticks.x=element_blank()) + theme(axis.ticks.y=element_blank()) + theme(axis.text.y=element_blank())  
map


########################################
#### Load and join leaf litter data ####
########################################

# SOURCE: LeRoy et al. (2020) https://doi.org/10.1111/1365-2745.13262
# DATA REPOSITORY: https://github.com/andrew-hipp/decomposition-phylogeny-2019
# DATA PROCESSING: See 'litter process.R'

litter <- read.csv("litter_processed.csv")
litter$logk <- log(litter$mean_kd)

# Turn litter data into spatial points data frame
litter_pts<-SpatialPointsDataFrame(litter[,c("longitude","latitude")],litter,
                               proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# Join nutrient data at litter data locations
merge_temp2 <- data.frame(litter_pts@data$Genus)
for(i in 1:length(nut_data)){
  tempdat <- raster::extract(nut_data[[i]],litter_pts,df=T)
  merge_temp2 <- cbind(merge_temp2,tempdat[,2])
  colnames(merge_temp2)[1+i] <- names(nut_data)[i]
}

litter_pts@data<-cbind(litter_pts@data,merge_temp2[,2:5])

# Convert litter points to sf object
litter_st<-as(litter_pts,"sf")

# Transform litter spatial points to Mollweide
litter_moll<-st_transform(litter_st,"+proj=moll")

# Spatially join litter stream points to attributes of nearest subwatershed
litter_stb<-st_join(litter_moll, bas_moll, join = st_nearest_feature)
litter_stb$geometry<-NULL

# Turn -999 values (no data as represented in Hydrobasins) to NA
litter_stb[litter_stb==-999]<-NA

# Litter database does not include time of deployment
# All month of deployment variables are NA
litter_stb$tempmonth<-NA
litter_stb$precipmonth<-NA
litter_stb$PETmonth<-NA
litter_stb$AETmonth<-NA
litter_stb$moist_indexmonth<-NA
litter_stb$snowmonth<-NA
litter_stb$soilwatermonth<-NA

# Build subset including predictor variables used for CELLDEX model
# Simultaneously, log transform variables that have large skew and back-transform 
# 10x and 100x corrections done in HydroBASINS

#Setup empty dataframe
Ldat <- data.frame(logk=litter_stb$logk)

#For loop to select variables and complete transformations 
for(i in 2:length(mod_vars$Variables)){
  if(mod_vars$Transform[i]=="none"){
    tempname <- mod_vars$Variables[i]
    Ldat <- cbind(Ldat,litter_stb[,mod_vars$Variables[i]])
    colnames(Ldat)[i] = tempname
  } else
    if(mod_vars$Transform[i]=="log"){
      tempname <- paste0("log10",mod_vars$Variables[i])
      Ldat <- cbind(Ldat,log10(litter_stb[,mod_vars$Variables[i]]))
      colnames(Ldat)[i] = tempname
    } else
      if(mod_vars$Transform[i]=="log1"){
        tempname <- paste0("log10",mod_vars$Variables[i])
        Ldat <- cbind(Ldat,log10(litter_stb[,mod_vars$Variables[i]]+1))
        colnames(Ldat)[i] = tempname
      } else
        if(mod_vars$Transform[i]=="xten"){
          tempname <- mod_vars$Variables[i]
          Ldat <- cbind(Ldat,(litter_stb[,mod_vars$Variables[i]]/10))
          colnames(Ldat)[i] = tempname
        } else
          if(mod_vars$Transform[i]=="xhund"){
            tempname <- mod_vars$Variables[i]
            Ldat <- cbind(Ldat,(litter_stb[,mod_vars$Variables[i]]/100))
            colnames(Ldat)[i] = tempname
          } 
}

dim(Ldat) #102 total variables, same as cotton dataframe

# Predict ln(k) for litter sites using stream ln(k) model (Cgbm) above
litter$ln_pred_k<-predict(Cgbm, newdata=Ldat, n.trees=best.iter)

# Combine ln(k) predictions with litter condition and and substrate quality
# in the form of genus-level leaf traits (mean) from TRY (Kattge et al., 2011)
# also genus-level litter traits (mean) from literature review.
# See 'litter process.R' for details.

traits <- read.csv("traits.csv")
litter2 <- merge(litter,traits,by='Genus')

length(unique(litter2$Genus)) #33 genera with leaf OR litter traits

# Clean up variables and make factors
litter2$Mesh.size<-factor(litter2$Mesh.size)
litter2$Leaf.condition<-factor(litter2$Leaf.condition)

# Build the dataframe containing only variables for the model
# Variables for the validation model include cotton kd, experimental conditions,
# chemistry of fresh leaves, and chemistry of senesced leaves 
val_var <- read.csv(file="validation_variables.csv")
Fdat <- litter2[,colnames(litter2) %in% val_var$Variables]

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


##########################################
#### Figure 3 - Litter validation BRT ####
##########################################

# Get grids of x variable values and predictions to make partial dependence plots
lnpredk<-plot(Fgbm,i.var=c('ln_pred_k','Mesh.size'),return.grid=TRUE)
c2n<-plot(Fgbm,i.var='CtoN_Litter_Mn',return.grid=TRUE)
cel<-plot(Fgbm,i.var='Cellulose_Litter_Mn',return.grid=TRUE)
lig<-plot(Fgbm,i.var='Lignin_Litter_Mn',return.grid=TRUE)
nmn<-plot(Fgbm,i.var='N_Litter_Mn',return.grid=TRUE)
cmn<-plot(Fgbm,i.var='C_Litter_Mn',return.grid=TRUE)

# Get quartiles of the x variable to mark in plots
kdrug<-as.data.frame(quantile(Fdat$ln_pred_k,probs = seq(0, 1, 0.25),na.rm=T))
colnames(kdrug)[1]<-"rug"
c2nrug<-as.data.frame(quantile(Fdat$CtoN_Litter_Mn,probs = seq(0, 1, 0.25),na.rm=T))
colnames(c2nrug)[1]<-"rug"
celrug<-as.data.frame(quantile(Fdat$Cellulose_Litter_Mn,probs = seq(0, 1, 0.25),na.rm=T))
colnames(celrug)[1]<-"rug"
ligrug<-as.data.frame(quantile(Fdat$Lignin_Litter_Mn,probs=seq(0, 1, 0.25),na.rm=T))
colnames(ligrug)[1]<-"rug"
nmnrug<-as.data.frame(quantile(Fdat$N_Litter_Mn,probs=seq(0, 1, 0.25),na.rm=T))
colnames(nmnrug)[1]<-"rug"
cmnrug<-as.data.frame(quantile(Fdat$C_Litter_Mn,probs=seq(0, 1, 0.25),na.rm=T))
colnames(cmnrug)[1]<-"rug"

# Make partial dependence plots in ggplot
cols<-c("brown","forestgreen")
pf1<-ggplot(data = lnpredk, aes(ln_pred_k, y)) +
  geom_smooth(aes(color=Mesh.size),method="gam",se=T) +
  geom_line(aes(color=Mesh.size),linetype=1,alpha=0.25,linewidth=0.25) +
  geom_point(aes(x=rug,y=min(lnpredk$y)),kdrug,color="gray",shape = 15) +
  scale_color_manual(values=cols) +
  ylab("") + xlab(bquote('ln Predicted' ~K[d])) + theme(legend.position=c(0.25,0.95)) + theme(legend.title = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  guides(linetype = guide_legend(ncol = 1)) + guides(color=guide_legend(ncol=2))

pf2<-ggplot(data = lig, aes(Lignin_Litter_Mn, y)) +
  geom_point(aes(x=rug,y=min(lig$y)),ligrug,color="gray",shape = 15) +
  geom_line(color = "steelblue", size = 1) +
  ylab("") + 
  xlab("Litter % lignin") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

pf3<-ggplot(data = c2n, aes(CtoN_Litter_Mn, y)) +
  geom_point(aes(x=rug,y=min(c2n$y)),c2nrug,color="gray",shape = 15) +
  geom_line(color = "steelblue", size = 1) +
  ylab("") + 
  xlab("Litter C:N molar ratio") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

pf4<-ggplot(data = nmn, aes(N_Litter_Mn, y)) +
  geom_point(aes(x=rug,y=min(nmn$y)),nmnrug,color="gray",shape = 15) +
  geom_line(color = "steelblue", size = 1) +
  ylab("") + 
  xlab("Litter % N content") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

pf5<-ggplot(data = cel, aes(Cellulose_Litter_Mn, y)) +
  geom_point(aes(x=rug,y=min(cel$y)),celrug,color="gray",shape = 15) +
  geom_line(color = "steelblue", size = 1) +
  ylab("") + xlab("Litter % cellulose") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

pf6<-ggplot(data = cmn, aes(C_Litter_Mn, y)) +
  geom_point(aes(x=rug,y=min(cmn$y)),cmnrug,color="gray",shape = 15) +
  geom_line(color = "steelblue", size = 1) +
  ylab("") + 
  xlab("Litter % C content") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Make trellis of the top 6 partial dependence plots
pdf(file = "validation_pdps_6_meshsize2.pdf",width=4,height=12)
grid.arrange(pf1,pf2,pf3,pf4,pf5,pf6, ncol = 1,left=textGrob(bquote('ln' ~K[d]),rot=90))
dev.off()

