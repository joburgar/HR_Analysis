#####################################################################################
# 01_MCP.R
# script to run MCPs
# following from 00_TelemDataPrep.R
# adapted from script written by genevieve perkins (genevieve.perkins@gov.bc.ca)
# modified by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 06-Oct-2019
#####################################################################################
.libPaths("C:/Program Files/R/R-3.6.0/library")# to ensure reading/writing libraries from C drive (H drive too slow)

# overall process: 
#- Define the area of interest per animal; availability/ vs use (MCP)
#- Export MCP shapefiles

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
# check the spread of elk with maped locations 
bc <- bc_bound()
SC <- nr_districts() %>% filter(ORG_UNIT %in% c("DCK", "DSQ", "DSC"))

# Plot by AnimalID (29 individuals)
unique(HR.sf$AnimalID)
ggplot() +
  geom_sf(data=SC, fill="white", col="gray") +
  geom_sf(data=HR.sf, aes(fill=AnimalID, col=AnimalID))+
  coord_sf() +
  theme_minimal() +
  ggtitle("Animal GPS locations")


# Plot by Animal_Season (29 animal_seasons)
unique(HR.sf.AS$Animal_Season)
ggplot() +
  geom_sf(data=SC, fill="white", col="gray") +
  geom_sf(data=HR.sf.AS, aes(fill=Animal_Season, col=Animal_Season))+
  coord_sf() +
  theme_minimal() +
  ggtitle("Animal_Season GPS locations")


# Plot by Animal_Year (24 individuals)
unique(HR.sf.AY$Animal_Year)
ggplot() +
  geom_sf(data=SC, fill="white", col="gray") +
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
HR.sf.AS_utm <- st_transform(HR.sf.AS, crs=26910) # utm and m units

HR.sf.AS_utm$geometry
HR.sp.AS <- as(HR.sf.AS_utm, "Spatial")
class(HR.sp.AS)

# Calculate MCPs for each animal season
names(HR.sp.AS)
ASmcp.95 <- mcp(HR.sp.AS[,c("Animal_Season")], percent = 95)
ASmcp.50 <- mcp(HR.sp.AS[,c("Animal_Season")], percent = 50)

# Plot
plot(HR.sp.AS) # looks similar to before
plot(ASmcp.95, col = alpha(1:73, 0.5), add = TRUE)
plot(ASmcp.50, col = alpha(1:73, 0.5), add = TRUE)

# create shapefiles
# add in meta data covariates
head(HR.df)

ASmcp.95.sf <- st_as_sf(ASmcp.95)
colnames(ASmcp.95.sf)[1] <- "Animal_Season" # change id back to Animal_Season
ASmcp.95.sf <- dplyr::left_join(ASmcp.95.sf, 
                                HR.df[c("Animal_Season","AnimalID","Group.New","Season","Species","Sex", "Age_Class")], 
                                by = "Animal_Season") 
ASmcp.95.sf$HR_Type <- "MCP_95"

ASmcp.50.sf <- st_as_sf(ASmcp.50)
colnames(ASmcp.50.sf)[1] <- "Animal_Season" # change id back to Animal_Season
ASmcp.50.sf <- dplyr::left_join(ASmcp.50.sf, 
                                HR.df[c("Animal_Season","AnimalID","Group.New","Season","Species","Sex", "Age_Class")], 
                                by = "Animal_Season") 
ASmcp.50.sf$HR_Type <- "MCP_50"


# combine and write as one shapefile
setwd(GISDir)

ASmcp.sf <- rbind(ASmcp.50.sf,ASmcp.95.sf)
st_write(ASmcp.sf, "MCP_Animal_Season.shp")

###---
# Calculate the MCP by including 50 to 100 percent of points
par(mar=c(1,1,1,1)) # to fit in window
AS.hrs <- mcp.area(HR.sp.AS[c("Animal_Season")], percent = seq(50, 100, by = 10))
# visual inspection shows much variation in animals use of home range between 50-100% of points
AS.hrs # examine dataframe
AS.hrs.df <- as.data.frame(AS.hrs)
setwd(OutputDir)
write.csv(AS.hrs.df, "MCP_HRS_Animal_Season.csv")

##############################################################
#### RUN MINIMUM CONVEX POLYGON - ANIMAL_SEASON (END)
#############################################################


##############################################################
#### RUN MINIMUM CONVEX POLYGON - ANIMAL_YEAR (BEGINNING)
#############################################################

###--- Change to a SpatialPointsDataFrame and set to utm (m units)
HR.sf.AY_utm <- st_transform(HR.sf.AY, crs=26910) # utm and m units

HR.sf.AY_utm$geometry
HR.sp.AY <- as(HR.sf.AY_utm, "Spatial")
class(HR.sp.AY)

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
head(HR.df)

AYmcp.95.sf <- st_as_sf(AYmcp.95)
colnames(AYmcp.95.sf)[1] <- "Animal_Year" # change id back to Animal_Year
AYmcp.95.sf <- dplyr::left_join(AYmcp.95.sf, 
                                HR.df[c("Animal_Year","AnimalID","Group.New","Year","Species","Sex", "Age_Class")], 
                                by = "Animal_Year") 
AYmcp.95.sf$HR_Type <- "MCP_95"

AYmcp.50.sf <- st_as_sf(AYmcp.50)
colnames(AYmcp.50.sf)[1] <- "Animal_Year" # change id back to Animal_Year
AYmcp.50.sf <- dplyr::left_join(AYmcp.50.sf, 
                                HR.df[c("Animal_Year","AnimalID","Group.New","Year","Species","Sex", "Age_Class")], 
                                by = "Animal_Year") 
AYmcp.50.sf$HR_Type <- "MCP_50"


# combine and write as one shapefile
setwd(GISDir)
AYmcp.sf <- rbind(AYmcp.50.sf,AYmcp.95.sf)
st_write(AYmcp.sf, "MCP_Animal_Year.shp")

###---
# Calculate the MCP by including 50 to 100 percent of points
par(mar=c(1,1,1,1)) # to fit in window
AY.hrs <- mcp.area(HR.sp.AY[c("Animal_Year")], percent = seq(50, 100, by = 10))
# visual inspection shows much variation in animals use of home range between 50-100% of points
AY.hrs # examine dataframe
AY.hrs.df <- as.data.frame(AY.hrs)
setwd(OutputDir)
write.csv(AY.hrs.df, "MCP_HRS_Animal_Year.csv")

##############################################################
#### RUN MINIMUM CONVEX POLYGON - ANIMAL_YEAR (END)
#############################################################


##############################################################
#### SUMMARISE MINIMUM CONVEX POLYGON ESTIMATES (BEGINNING)
#############################################################

###--- summarise data

# Animal Year by Group Type
AYmcp.sf %>% group_by(HR_Type, Year) %>% 
  summarise(mean = mean(area), se = sd(area)/sqrt(n())) %>% st_drop_geometry()

# HR_Type  Year   mean     se
# 1 MCP_50   2017 14425. 1128. 
# 2 MCP_50   2018   276.   21.3
# 3 MCP_95   2017 46098. 3147. 
# 4 MCP_95   2018  1108.  108. 

# plot out the MCP area sensitivty by year and group type
setwd(OutputDir)
AY.mcp.sens <- ggplot(AYmcp.sf, aes(x=as.factor(Year), y=area, color=HR_Type)) +
  geom_boxplot() +
  scale_color_brewer(palette="Dark2") +
  labs(title="MCP Area Sensitivity", y = "MCP Area (ha)")+
  theme(axis.title.x=element_blank()) +
  facet_wrap(.~Group.New, scales="free_y")

Cairo(800, 500, pointsize = 36,
      file="MCP_Year_Group.png", type="png", bg="white")
AY.mcp.sens
dev.off()

# Animal season by Group Type
ASmcp.sf %>% group_by(HR_Type, Season) %>% 
  summarise(mean = mean(area), se = sd(area)/sqrt(n())) %>% st_drop_geometry()

# HR_Type Season         mean    se
# 1 MCP_50  Breeding      3272.  261.
# 2 MCP_50  Non-Breeding  2408.  197.
# 3 MCP_95  Breeding     12972.  963.
# 4 MCP_95  Non-Breeding  8347.  667.

# plot out the MCP area sensitivty by year and group type
setwd(OutputDir)
AS.mcp.sens <- ggplot(ASmcp.sf, aes(x=Season, y=area, color=HR_Type)) +
  geom_boxplot() +
  scale_color_brewer(palette="Dark2") +
  labs(title="MCP Area Sensitivity", y = "MCP Area (ha)")+
  theme(axis.title.x=element_blank()) +
  facet_wrap(.~Group.New, scales="free_y")

Cairo(800, 500, pointsize = 36,
      file="MCP_Season_Group.png", type="png", bg="white")
AS.mcp.sens
dev.off()

##############################################################
#### SUMMARISE MINIMUM CONVEX POLYGON ESTIMATES (END)
#############################################################


##############################################################
#### HOUSEKEEPING (BEGINNING)
#############################################################

###-- Save workspace and move to KDE home range analysis - 02_KDE.R
setwd(InputDir)
save.image("01_MCP.RData")

##############################################################
#### HOUSEKEEPING (END)
#############################################################


