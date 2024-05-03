#####################################################################################
# 02_KDE.R
# script to run KDEs
# following from 00_TelemDataPrep.R (and O1_MCP.R)
# adapted from script written by genevieve perkins (genevieve.perkins@gov.bc.ca)
# modified by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 06-Oct-2019
#####################################################################################
#####################################################################################
# Copyright 2021 Province of British Columbia
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and limitations under the License.
#####################################################################################

R_version <- paste0("R-",version$major,".",version$minor)

.libPaths(paste0("C:/Program Files/R/",R_version,"/library")) # to ensure reading/writing libraries from C drive
tz = Sys.timezone() # specify timezone in BC

GISDir <- "//spatialfiles.bcgov/work/wlap/sry/Workarea/jburgar"

# overall process: 
#- Define the area of interest per animal; availability/ vs use (MCP)
#- Export KDe shapefiles

# help files: 
#https://cran.r-project.org/web/packages/adehabitatHS/adehabitatHS.pdf

# run libraries
library(bcmaps)
library(tidyverse)
library(Cairo)
library(sf)
library(sp)
library(adehabitatHR)

##############################################################
#### LOAD and REVIEW DATA (BEGINNING)
#############################################################
## load data into R
load("HR_InputData.RData")

# review spatial object data - check to see if loaded properly
st_geometry(HR.sf) # currently in 4326 CRS, lat/long
HR.sf$AnimalID <- as.factor(HR.sf$AnimalID)
glimpse(HR.sf)

# plot to check
# check the spread of animals with maped locations 
bc <- bc_bound()
NRD <- nr_districts() %>% filter(ORG_UNIT %in% c("DCC", "DMH"))

# Plot by AnimalID
unique(HR.sf$AnimalID) #2 individuals
ggplot() +
  geom_sf(data=NRD, fill="white", col="gray") +
  geom_sf(data=HR.sf, aes(fill=AnimalID, col=AnimalID))+
  coord_sf() +
  theme_minimal() +
  ggtitle("Animal GPS locations")

# Plot by Animal_Year
unique(HR.sf.AY$Animal_Year) # 2 individuals
ggplot() +
  geom_sf(data=NRD, fill="white", col="gray") +
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

###--- ANIMAL_YEAR 
st_geometry(HR.sf.AY)
HR.sf.AY_utm <- st_transform(HR.sf.AY, crs=26910) # utm and m units

HR.sf.AY_utm$geometry
HR.sp.AY <- as(HR.sf.AY_utm, "Spatial")
class(HR.sp.AY)

kde1.AY  <- kernelUD(HR.sp.AY[c("Animal_Year")], h = "href", kern = c("bivnorm"), extent = 2) # default grid and extent

length(kde1.AY) # 2 Animal_Year

# create vector listing the h-value for each KDE
kde1.AY_href <- rep(NA, length(kde1.AY))
for (i in 1:length(kde1.AY_href )){ 
  i.href <- kde1.AY[[i]]@h$h
  kde1.AY_href[i] <- i.href
}

min(kde1.AY_href); max(kde1.AY_href); mean(kde1.AY_href)
# [1] 603.8
# [1] 789.1566
# [1] 696.4783
#  not much variation in href depending on Animal_Year

# create KDEs
ver95 <- getverticeshr(kde1.AY,95) ;   ver95.sf<- st_as_sf(ver95)
ver50 <- getverticeshr(kde1.AY,50) ;   ver50.sf<- st_as_sf(ver50)

plot(st_geometry(ver95.sf),col = "red")
plot(st_geometry(ver50.sf),col = "yellow", add = TRUE)
plot(HR.sp.AY, pch = 1, size = 0.5, add = TRUE)     # Add points 

# add in meta-data and then export shapefiles
anml$Animal_Year <- paste(anml$AnimalID, anml$Cptr_Year, sep="_")

colnames(ver95.sf)[1] <- "Animal_Year"
ver95.sf <- left_join(ver95.sf, anml %>% dplyr::select("Animal_Year","AnimalID","Sex", "Species","Age_Class"))
ver95.sf$HR_Type <- "KDE_95"

colnames(ver50.sf)[1] <- "Animal_Year"
ver50.sf <- dplyr::left_join(ver50.sf, anml %>% dplyr::select("Animal_Year","AnimalID","Sex", "Species","Age_Class")) 
ver50.sf$HR_Type <- "KDE_50"

ver1AY.sf <- rbind(ver50.sf,ver95.sf)

st_write(ver1AY.sf,"KDE_Animal_Year_href.shp")


st_write(ver1AY.sf %>% filter(AnimalID==35060),"KDE_Fisher_35060.shp")

ggplot()+
  geom_sf(data=ver1AY.sf %>% filter(AnimalID==35060))

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

# for each animal 
KDE.function <- function(data, HR_group, hvalue, grid, extent){
  kde <- kernelUD(data[c(HR_group)], h = hvalue, kern = c("bivnorm"), grid = grid, extent = extent)
  ver95 <- getverticeshr(kde,95) ;   ver95.sf<- st_as_sf(ver95)
  ver50 <- getverticeshr(kde,50) ;   ver50.sf<- st_as_sf(ver50)
  
  colnames(ver95.sf)[1] <- HR_group
  ver95.sf <- dplyr::left_join(ver95.sf, 
                               unique(HR.df %>% dplyr::select(HR_group,"Collar","Sex", "Age_Class")), 
                               by = HR_group) 
  ver95.sf$HR_Type <- "KDE_95"
  
  colnames(ver50.sf)[1] <- HR_group
  ver50.sf <- dplyr::left_join(ver50.sf, 
                               unique(HR.df %>% dplyr::select(HR_group,"Collar","Sex", "Age_Class")), 
                               by = HR_group) 
  ver50.sf$HR_Type <- "KDE_50"
  
  ver.sf <- rbind(ver50.sf,ver95.sf)
  
  st_write(ver.sf, paste("KDE_",HR_group,"_h", hvalue,".shp", sep=""))}

# for each animal-year
KDE.function.Yr <- function(data, HR_group, hvalue, grid, extent){
  kde <- kernelUD(data[c(HR_group)], h = hvalue, kern = c("bivnorm"), grid = grid, extent = extent)
  ver95 <- getverticeshr(kde,95) ;   ver95.sf<- st_as_sf(ver95)
  ver50 <- getverticeshr(kde,50) ;   ver50.sf<- st_as_sf(ver50)
  
  colnames(ver95.sf)[1] <- HR_group
  ver95.sf <- dplyr::left_join(ver95.sf, 
                                 unique(HR.df %>% dplyr::select(HR_group,"AnimalID","Collar","Year","Sex", "Age_Class")), 
                                 by = HR_group) 
  ver95.sf$HR_Type <- "KDE_95"
  
  colnames(ver50.sf)[1] <- HR_group
  ver50.sf <- dplyr::left_join(ver50.sf, 
                                 unique(HR.df %>% dplyr::select(HR_group,"AnimalID","Collar","Year","Sex", "Age_Class")), 
                                 by = HR_group) 
  ver50.sf$HR_Type <- "KDE_50"
  
  ver.sf <- rbind(ver50.sf,ver95.sf)
  
  st_write(ver.sf, paste("KDE_",HR_group,"_h", hvalue,".shp", sep=""))}

# for each animal-season
KDE.function.Sn <- function(data, HR_group, hvalue, grid, extent){
  kde <- kernelUD(data[c(HR_group)], h = hvalue, kern = c("bivnorm"), grid = grid, extent = extent)
  ver95 <- getverticeshr(kde,95) ;   ver95.sf<- st_as_sf(ver95)
  ver50 <- getverticeshr(kde,50) ;   ver50.sf<- st_as_sf(ver50)
  
  colnames(ver95.sf)[1] <- HR_group
  ver95.sf <- dplyr::left_join(ver95.sf, 
                               unique(HR.df %>% dplyr::select(HR_group,"AnimalID","Collar","Season","Sex", "Age_Class")), 
                               by = HR_group) 
  ver95.sf$HR_Type <- "KDE_95"
  
  colnames(ver50.sf)[1] <- HR_group
  ver50.sf <- dplyr::left_join(ver50.sf, 
                               unique(HR.df %>% dplyr::select(HR_group,"AnimalID","Collar","Season","Sex", "Age_Class")), 
                               by = HR_group) 
  ver50.sf$HR_Type <- "KDE_50"
  
  ver.sf <- rbind(ver50.sf,ver95.sf)
  
  st_write(ver.sf, paste("KDE_",HR_group,"_h", hvalue,".shp", sep=""))}

###--- end of functions



###################################################################

###--- Run kernel density estimates for HR_groups (Animal_Season and Animal_Year) 
# with user selected h ref (1000, 2000)
# try for grid = 500 and extent = 5 but may need to adjust

setwd(GISDir)

AS_h1000 <- KDE.function.Sn(data = HR.sp.AS, HR_group = "Animal_Season", hvalue = 1000, grid =500, extent = 5)

AY_h1000 <- KDE.function.Yr(data = HR.sp.AY, HR_group = "Animal_Year", hvalue = 1000, grid =500, extent = 5)

AS_h2000 <- KDE.function.Sn(data = HR.sp.AS, HR_group = "Animal_Season", hvalue = 2000, grid =500, extent = 5)

AY_h2000 <- KDE.function.Yr(data = HR.sp.AY, HR_group = "Animal_Year", hvalue = 2000, grid =500, extent = 5)

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
