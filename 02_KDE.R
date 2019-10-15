#####################################################################################
# 02_KDE.R
# script to run KDEs
# following from 00_TelemDataPrep.R (and O1_MCP.R)
# adapted from script written by genevieve perkins (genevieve.perkins@gov.bc.ca)
# modified by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 06-Oct-2019
#####################################################################################
.libPaths("C:/Program Files/R/R-3.6.0/library")# to ensure reading/writing libraries from C drive (H drive too slow)

# overall process: 
#- Define the area of interest per animal; availability/ vs use (MCP)
#- Export KDe shapefiles

# help files: 
#https://cran.r-project.org/web/packages/adehabitatHS/adehabitatHS.pdf

# run libraries
library(bcmaps)
library(dplyr)
library(ggplot2)
library(Cairo)
library(sf)
library(sp)
library(adehabitatHR)

# set up working directories on H drive
InputDir <- c("H:/R/Analysis/Generic_HomeRange/Input")
OutputDir <- c("H:/R/Analysis/Generic_HomeRange/Output")
GISDir <- c("H:/R/Analysis/Generic_HomeRange/GISDir")

##############################################################
#### LOAD and REVIEW DATA (BEGINNING)
#############################################################
## load data into R
setwd(InputDir)
load("HR_InputData.RData")

# review spatial object data - check to see if loaded properly
st_geometry(HR.sf) # currently in 4326 CRS, lat/long
HR.sf$AnimalID <- as.factor(HR.sf$AnimalID)
names(HR.sf)
summary(HR.sf)

# plot to check
# check the spread of animals with maped locations 
bc <- bc_bound()
SC <- nr_districts() %>% filter(ORG_UNIT %in% c("DCK", "DSQ", "DSC")) # revise as appropriate to study area

# Plot by AnimalID
unique(HR.sf$AnimalID) #29 individuals
ggplot() +
  geom_sf(data=SC, fill="white", col="gray") +
  geom_sf(data=HR.sf, aes(fill=AnimalID, col=AnimalID))+
  coord_sf() +
  theme_minimal() +
  ggtitle("Animal GPS locations")


# Plot by Animal_Season
unique(HR.sf.AS$Animal_Season) # 29 animal_seasons
ggplot() +
  geom_sf(data=SC, fill="white", col="gray") +
  geom_sf(data=HR.sf.AS, aes(fill=Animal_Season, col=Animal_Season))+
  coord_sf() +
  theme_minimal() +
  ggtitle("Animal_Season GPS locations")


# Plot by Animal_Year
unique(HR.sf.AY$Animal_Year) # 24 individuals
ggplot() +
  geom_sf(data=SC, fill="white", col="gray") +
  geom_sf(data=HR.sf.AY, aes(fill=Animal_Year, col=Animal_Year))+
  coord_sf() +
  theme_minimal() +
  ggtitle("Animal_Year GPS locations")

##############################################################
#### LOAD and REVIEW DATA (END)
#############################################################


###################################################################
#### RUN HREF KDE MODEL - HREF (BEGINNING)
###################################################################

# smoothing parameter h controls the "width" of the kernel functions placed over each point
# the larger the value of h the larger the UD
# the reference bandwith option supposes that the UD is a bivariate normal distribution and will overestimate if the animal has multiple activity centres
# first try with "href", not specifying grid (need to specify extent as default not working with 95%)

###--- ANIMAL_SEASON 
st_geometry(HR.sf.AS)
HR.sf.AS_utm <- st_transform(HR.sf.AS, crs=26910) # utm and m units

HR.sf.AS_utm$geometry
HR.sp.AS <- as(HR.sf.AS_utm, "Spatial")
class(HR.sp.AS)

kde1.AS  <- kernelUD(HR.sp.AS[c("Animal_Season")], h = "href", kern = c("bivnorm"), extent = 2) # default grid and extent

length(kde1.AS) # 29 Animal_Season

# create vector listing the h-value for each KDE
kde1.AS_href <- rep(NA, length(kde1.AS))
for (i in 1:length(kde1.AS_href )){ 
  i.href <- kde1.AS[[i]]@h$h
  kde1.AS_href[i] <- i.href
}

min(kde1.AS_href); max(kde1.AS_href); mean(kde1.AS_href)
# [1] 105.0504
# [1] 9484.046
# [1] 1546.877
# huge variation in href depending on Animal_Season

# create KDEs
ver95 <- getverticeshr(kde1.AS,95) ;   ver95.sf<- st_as_sf(ver95)
ver50 <- getverticeshr(kde1.AS,50) ;   ver50.sf<- st_as_sf(ver50)

plot(st_geometry(ver95.sf),col = "red")
plot(st_geometry(ver50.sf),col = "yellow", add = TRUE)
plot(HR.sp.AS, pch = 1, size = 0.5, add = TRUE)     # Add points 

# add in meta-data and then export shapefiles
colnames(ver95.sf)[1] <- "Animal_Season"
ver95.sf <- dplyr::left_join(ver95.sf, 
                             unique(HR.df %>% dplyr::select("Animal_Season","AnimalID","Group.New","Season","Species","Sex", "Age_Class")), 
                             by = "Animal_Season")
ver95.sf$HR_Type <- "KDE_95"

colnames(ver50.sf)[1] <- "Animal_Season"
ver50.sf <- dplyr::left_join(ver50.sf, 
                             unique(HR.df %>% dplyr::select("Animal_Season","AnimalID","Group.New","Season","Species","Sex", "Age_Class")), 
                             by = "Animal_Season") 
ver50.sf$HR_Type <- "KDE_50"

ver1AS.sf <- rbind(ver50.sf,ver95.sf)

setwd(GISDir)
st_write(ver1AS.sf,"KDE_Animal_Season_href.shp")

###--- ANIMAL_YEAR 
st_geometry(HR.sf.AY)
HR.sf.AY_utm <- st_transform(HR.sf.AY, crs=26910) # utm and m units

HR.sf.AY_utm$geometry
HR.sp.AY <- as(HR.sf.AY_utm, "Spatial")
class(HR.sp.AY)

kde1.AY  <- kernelUD(HR.sp.AY[c("Animal_Year")], h = "href", kern = c("bivnorm"), extent = 2) # default grid and extent

length(kde1.AY) # 24 Animal_Year

# create vector listing the h-value for each KDE
kde1.AY_href <- rep(NA, length(kde1.AY))
for (i in 1:length(kde1.AY_href )){ 
  i.href <- kde1.AY[[i]]@h$h
  kde1.AY_href[i] <- i.href
}

min(kde1.AY_href); max(kde1.AY_href); mean(kde1.AY_href)
# [1] 105.0504
# [1] 13872.24
# [1] 1835.504
# huge variation in href depending on Animal_Year

# create KDEs
ver95 <- getverticeshr(kde1.AY,95) ;   ver95.sf<- st_as_sf(ver95)
ver50 <- getverticeshr(kde1.AY,50) ;   ver50.sf<- st_as_sf(ver50)

plot(st_geometry(ver95.sf),col = "red")
plot(st_geometry(ver50.sf),col = "yellow", add = TRUE)
plot(HR.sp.AY, pch = 1, size = 0.5, add = TRUE)     # Add points 

# add in meta-data and then export shapefiles
colnames(ver95.sf)[1] <- "Animal_Year"
ver95.sf <- dplyr::left_join(ver95.sf, 
                             unique(HR.df %>% dplyr::select("Animal_Year","AnimalID","Group.New","Year","Species","Sex", "Age_Class")), 
                             by = "Animal_Year") 
ver95.sf$HR_Type <- "KDE_95"

colnames(ver50.sf)[1] <- "Animal_Year"
ver50.sf <- dplyr::left_join(ver50.sf, 
                             unique(HR.df %>% dplyr::select("Animal_Year","AnimalID","Group.New","Year","Species","Sex", "Age_Class")), 
                             by = "Animal_Year") 
ver50.sf$HR_Type <- "KDE_50"

ver1AY.sf <- rbind(ver50.sf,ver95.sf)

setwd(GISDir)
st_write(ver1AY.sf,"KDE_Animal_Year_href.shp")

# for housekeeping and to reduce space, remove kde objects once output has been saved
rm(kde1.AS, kde1.AY)

###################################################################
#### RUN HREF KDE MODEL (END)
###################################################################


###################################################################
#### RUN USER DEFINED KDE MODELS (BEGINNING)
###################################################################

###--- create KDE function to run user defined models 

# data =  the sp object; e.g., HR.sp.AS
# hvalue = user defined h value; e.g., 1000
# HR_group =  HR grouping from 00_TelemDataPrep.R; e.g., "Animal_Season"
# grid = user defined grid value; e.g., 500
# extent = user defined extent value; e.g., 2

KDE.function <- function(data, HR_group, hvalue, grid, extent){
  kde <- kernelUD(data[c(HR_group)], h = hvalue, kern = c("bivnorm"), grid = grid, extent = extent)
  ver95.2 <- getverticeshr(kde,95) ;   ver95.2.sf<- st_as_sf(ver95.2)
  ver50.2 <- getverticeshr(kde,50) ;   ver50.2.sf<- st_as_sf(ver50.2)
  
  colnames(ver95.2.sf)[1] <- HR_group
  ver95.2.sf <- dplyr::left_join(ver95.2.sf, 
                                 unique(HR.df %>% dplyr::select(HR_group,"AnimalID","Group.New","Season","Year","Species","Sex", "Age_Class")), 
                                 by = HR_group) 
  ver95.sf$HR_Type <- "KDE_95"
  
  colnames(ver50.2.sf)[1] <- HR_group
  ver50.2.sf <- dplyr::left_join(ver50.2.sf, 
                                 unique(HR.df %>% dplyr::select(HR_group,"AnimalID","Group.New","Season","Year","Species","Sex", "Age_Class")), 
                                 by = HR_group) 
  ver50.sf$HR_Type <- "KDE_50"
  
  ver.sf <- rbind(ver50.2.sf,ver95.2.sf)
  
  st_write(ver.sf, paste("KDE_",HR_group,"_h", hvalue,".shp", sep=""))}

###--- end of function

###################################################################

###--- Run kernel density estimates for HR_groups (Animal_Season and Animal_Year) 
# with user selected h ref (1000, 2000)
# try for grid = 500 and extent = 5 but may need to adjust

setwd(GISDir)

AS_h1000 <- KDE.function(data = HR.sp.AS, HR_group = "Animal_Season", hvalue = 1000, grid =500, extent = 5)

AY_h1000 <- KDE.function(data = HR.sp.AY, HR_group = "Animal_Year", hvalue = 1000, grid =500, extent = 5)

AS_h2000 <- KDE.function(data = HR.sp.AS, HR_group = "Animal_Season", hvalue = 2000, grid =500, extent = 5)

AY_h2000 <- KDE.function(data = HR.sp.AY, HR_group = "Animal_Year", hvalue = 2000, grid =500, extent = 5)

###################################################################
#### RUN USER DEFINED KDE MODELS (END)
###################################################################


##############################################################
#### HOUSEKEEPING (BEGINNING)
#############################################################

###-- Save workspace
setwd(InputDir)
save.image("02_KDE.RData")
#load("02_KDE.RData")

##############################################################
#### HOUSEKEEPING (END)
#############################################################
