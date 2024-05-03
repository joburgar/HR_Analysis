#####################################################################################
# 01_MCP.R
# script to run MCPs
# following from 00_TelemDataPrep.R
# adapted from script written by genevieve perkins (genevieve.perkins@gov.bc.ca)
# modified by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 06-Oct-2019
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
#- Export MCP shapefiles

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
# check the spread of fisher with mapped locations 
bc <- bc_bound()
NRD <- nr_districts() %>% filter(ORG_UNIT %in% c("DCC", "DMH"))

# Plot by AnimalID (2 individuals in 2 Districts)
unique(HR.sf$AnimalID)
ggplot() +
  geom_sf(data=NRD, fill="white", col="gray") +
  geom_sf(data=HR.sf, aes(fill=AnimalID, col=AnimalID))+
  coord_sf() +
  theme_minimal() +
  ggtitle("Animal GPS locations")


# Plot by Animal_Year (2 individuals)
unique(HR.sf.AY$Animal_Year)
ggplot() +
  geom_sf(data=NRD, fill="white", col="gray") +
  geom_sf(data=HR.sf.AY, aes(fill=Animal_Year, col=Animal_Year))+
  coord_sf() +
  theme_minimal() +
  ggtitle("Animal_Year GPS locations")

##############################################################
#### LOAD and REVIEW DATA (END)
#############################################################

##############################################################
#### RUN MINIMUM CONVEX POLYGON - ANIMAL_SEASON (BEGINNING)
#############################################################

###--- Change to a SpatialPointsDataFrame and set to utm (m units)
HR.sf.AY_utm <- st_transform(HR.sf.AY, crs=26910) # utm and m units
HR.sf.AY_utm %>% count(Animal_Year)

HR.sf.AY_utm$geometry
HR.sp.AY <- as(HR.sf.AY_utm, "Spatial"); class(HR.sp.AY)

# Calculate MCPs for each animal season
names(HR.sp.AY)
AYmcp.95 <- mcp(HR.sp.AY[,c("Animal_Year")], percent = 95)
AYmcp.50 <- mcp(HR.sp.AY[,c("Animal_Year")], percent = 50)

# Plot
plot(HR.sp.AY) # looks similar to before
plot(AYmcp.95, col = alpha(1:73, 0.5), add = TRUE)
plot(AYmcp.50, col = alpha(1:73, 0.5), add = TRUE)

# create shapefiles
# add in meta data covariates
glimpse(anml)
glimpse(AYmcp.95.sf)
anml$Animal_Year <- paste(anml$AnimalID, anml$Cptr_Year, sep="_")

AYmcp.95.sf <- st_as_sf(AYmcp.95)
colnames(AYmcp.95.sf)[1] <- "Animal_Year" # change id back to Animal_Year
AYmcp.95.sf <- left_join(AYmcp.95.sf, anml %>% dplyr::select("Animal_Year","AnimalID","Sex", "Species","Age_Class"))
AYmcp.95.sf$HR_Type <- "MCP_95"

AYmcp.50.sf <- st_as_sf(AYmcp.50)
colnames(AYmcp.50.sf)[1] <- "Animal_Year" # change id back to Animal_Year
AYmcp.50.sf <- left_join(AYmcp.50.sf, anml %>% dplyr::select("Animal_Year","AnimalID","Sex", "Species","Age_Class"))
AYmcp.50.sf$HR_Type <- "MCP_50"


# combine and write as one shapefile
AYmcp.sf <- rbind(AYmcp.50.sf,AYmcp.95.sf)
st_write(AYmcp.sf, "MCP_Fisher_Year.shp")

ggplot()+
  geom_sf(data=AYmcp.sf %>% filter(HR_Type=="MCP_95"), aes(col=AnimalID))+
  geom_sf(data=AYmcp.sf %>% filter(HR_Type=="MCP_50"), aes(col=AnimalID))


###---
# Calculate the MCP by including 50 to 100 percent of points
par(mar=c(1,1,1,1)) # to fit in window
AY.hrs <- mcp.area(HR.sp.AY[c("Animal_Year")], percent = seq(50, 100, by = 10))
# visual inspection shows much variation in animals use of home range between 50-100% of points
AY.hrs # examine dataframe
AY.hrs.df <- as.data.frame(AY.hrs)
write.csv(AY.hrs.df, "MCP_HRS_Animal_Year.csv")

##############################################################
#### RUN MINIMUM CONVEX POLYGON - ANIMAL_Year (END)
#############################################################


##############################################################
#### SUMMARISE MINIMUM CONVEX POLYGON ESTIMATES (BEGINNING)
#############################################################

###--- summarise data
AYmcp.sf$Year <- "2024"
# Animal Year by Group Type
AYmcp.sf %>% group_by(HR_Type) %>% 
  summarise(mean = mean(area), se = sd(area)/sqrt(n())) %>% st_drop_geometry()

# HR_Type  Year   mean     se
# 1 MCP_50  2024  1981. 1398.
# 2 MCP_95  2024  6532. 3952.
# plot out the MCP area sensitivity by year and group type
setwd(OutputDir)
AY.mcp.sens <- ggplot(AYmcp.sf, aes(x=as.factor(Sex), y=area, color=HR_Type)) +
  geom_boxplot() +
  scale_color_brewer(palette="Dark2") +
  labs(title="MCP Area Sensitivity", y = "MCP Area (ha)")+
  theme(axis.title.x=element_blank()) +
  facet_wrap(.~Year, scales="free_y")

Cairo(800, 500, pointsize = 36,
      file="MCP_Year_Sex.png", type="png", bg="white")
AY.mcp.sens
dev.off()


##############################################################
#### SUMMARISE MINIMUM CONVEX POLYGON ESTIMATES (END)
#############################################################


##############################################################
#### HOUSEKEEPING (BEGINNING)
#############################################################

###-- Save workspace and move to KDE home range analysis - 02_KDE.R
save.image("01_MCP.RData")

##############################################################
#### HOUSEKEEPING (END)
#############################################################


