#####################################################################################
# 00_TelemDataPrep.R
# script to prepare telemetry data for running home range analyses
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - Updated 2-May-2024
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

GISDir <- "//spatialfiles.bcgov/work/wlap/sry/Workarea/jburgar/MMP"

# overall process: 
#- Upload Animal metadata and collar metadata (if applicable)
#- Upload shapefiles or csv of telemetry data
#- modified for fisher GPS collar data
#- Output R object ready for home range analyses (MCP & KDE) 

# run libraries
library(tidyverse)  # for plotting, data manipulation, formatting character data
library(lubridate)  # for date-time conversions
library(sf)         # for uploading shapefiles and working with sf objects
library(bcdata); library(bcmaps)


##############################################################
#### METADATA EXPLORATION & FORMATTING (BEGINNING) 
#############################################################
list.files()

###--- upload animal and collar metadatda files
anml.full <- read.csv("Input/PEPE_Metadata.csv", # point to appropriate metadata file
                     header = TRUE, stringsAsFactors = TRUE, na.strings=c("", "NA"),)

glimpse(anml.full) # view data
head(anml.full)

# # subset to smaller dataframe, and rename columns for consistency
names(anml.full)
anml <- anml.full[c(1:8)] # if not in the suggested order, re-order to reflect order of colnames below
colnames(anml) <-  c("AnimalID","Species","Sex", "Age_Class","Cptr_Northing","Cptr_Easting", "Cptr_Date","Rls_Date")
head(anml)  # check

# format dates for R
anml$Cptr_Datep <- as.POSIXct(strptime(anml$Cptr_Date, format = "%d-%b-%y"))
anml$Rls_Datep <- as.POSIXct(strptime(anml$Rls_Date, format = "%d-%b-%y"))

# add in Year
anml$Cptr_Year <- year(anml$Cptr_Datep)
glimpse(anml) # check dates

###--- Check  data for number of captures per year by species, sex and min/max annual capture dates
anml %>% group_by(Cptr_Year, Species) %>% count(Sex)
anml %>% group_by(Cptr_Year, Species) %>% summarise(min(Cptr_Datep), max(Cptr_Datep))

###--- to turn the capture points into spatial object, using the CRS for NAD83 UTM Zone 10
# may need to clean some data during the import
Cpt.sf <- st_as_sf(anml, coords=c("Cptr_Easting","Cptr_Northing"), crs=26910)

# plot to check
bc <- bc_bound()

ggplot() +
  geom_sf(data = bc) +
  geom_sf(data = Cpt.sf, aes(fill=Species, col=Species)) +
  coord_sf() +
  theme_minimal()

ggplot() +
  geom_sf(data = bc) +
  geom_sf(data = Cpt.sf, aes(fill=as.factor(Cptr_Year), col=as.factor(Cptr_Year))) +
  coord_sf() +
  theme_minimal()

ggplot() +
  geom_sf(data = bc) +
  geom_sf(data = Cpt.sf, aes(fill=Sex, col=Sex)) +
  coord_sf() +
  theme_minimal()




MMP_BCGov <- st_read(dsn=GISDir, layer="MMP_BCGov_grid_March2024") 

##############################################################
#### METADATA EXPLORATION & FORMATTING (END) 
#############################################################

##############################################################
#### TELEMETRY DATA EXPLORATION & FORMATTING (BEGINNING) 
#############################################################

#############################################################

###--- function to retrieve geodata from BCGW

retrieve_geodata_aoi <- function (ID=ID){
  aoi.geodata <- bcdc_query_geodata(ID) %>%
    filter(BBOX(st_bbox(aoi))) %>%
    collect()
  aoi.geodata <- aoi.geodata %>% st_intersection(aoi)
  aoi.geodata$Area_km2 <- st_area(aoi.geodata)*1e-6
  aoi.geodata <- drop_units(aoi.geodata)
  return(aoi.geodata)
}


MMP_BCGov_grid_March2024
# biogeoclimatic zones
# bcdc_search("Biogeoclimatic zone", res_format = "wms")
# 3: BEC Map (other, wms, kml)
# ID: f358a53b-ffde-4830-a325-a5a03ff672c3
# Name: bec-map
aoi.BEC <- retrieve_geodata_aoi(ID = "f358a53b-ffde-4830-a325-a5a03ff672c3")


###--- For csv upload
telem <- read.csv("BDOW_TelemData.csv", header=TRUE, 
                  na.strings = c("NA",""), stringsAsFactors = TRUE)
glimpse(telem) # check columns are coming in as appropriate class
head(telem) # check that data is reading correctly
summary(telem) # check if spelling inconsistencies or NAs
#############################################################

###--- check data quality and remove any objectional rows
# latitude coordinate off - delete that row
telem <- telem[order(telem$Latitude),]
tail(telem) # only 1 latitiude coordinate off (probably a typo where 9 should be 5 but delete to be safe) 
# otherwise data was cleaned prior to uploading and is good to go

# use the CRS for lat/long and WGS84
telem.sf <- st_as_sf(telem[telem$Latitude<90,], coords=c("Longitude","Latitude"), crs=4326) %>% 
  select(USFW.Band, Date, Time, Group) %>% rename(AnimalID = USFW.Band) # selects only specific columns and chagnes the name to match column name in anml object

summary(telem.sf)
# create new Group names - revise as appropriate
telem.sf$Group.New <- as.factor(ifelse(grepl("Captive", telem.sf$Group),"Captive",
                                       ifelse(grepl("Relocated", telem.sf$Group), "Relocated",
                                              ifelse(grepl("Control", telem.sf$Group), "Control",
                                                     ifelse(grepl("Resident", telem.sf$Group), "Resident", NA)))))

table(telem.sf$Group, telem.sf$Group.New) # check that grouping worked

# format dates for R
telem.sf$Date.Time <- paste(telem.sf$Date, telem.sf$Time, sep=" ")
telem.sf$Date.Timep <- as.POSIXct(strptime(telem.sf$Date.Time, format = "%Y%m%d %H:%M:%S"))
telem.sf$Year <- year(telem.sf$Date.Timep)
telem.sf$Month <- month(telem.sf$Date.Timep)
telem.sf$Day.j <- julian(telem.sf$Date.Timep) # Julian day

glimpse(telem.sf)
summary(telem.sf)

# plot to check
ggplot() +
  geom_sf(data = telem.sf, aes(fill=Group.New, col=Group.New)) +
  coord_sf() +
  theme_minimal()

ggplot() +
  geom_sf(data = telem.sf, aes(fill=as.factor(Year), col=as.factor(Year))) +
  coord_sf() +
  theme_minimal()

ggplot() +
  geom_sf(data = telem.sf, aes(fill=AnimalID, col=AnimalID)) +
  facet_wrap(~Year) +
  theme(legend.position = "bottom") +
  coord_sf() 

# looks good - all plotting in general area

###--- Create Season (breeding) dates
# check date range by year for telemetry data
telem.sf %>% group_by(Year) %>% summarise(min(Date.Timep), max(Date.Timep)) %>% st_drop_geometry()
# Year `min(Date.Timep)`   `max(Date.Timep)`  
# 1  2017 2017-04-21 09:59:28 2017-10-31 10:01:36
# 2  2018 2018-03-11 18:40:00 2018-10-25 02:01:36

# Breeding = Apr 1 to Sep 30; Non-Breeding = Oct 1 - Mar 30
telem.sf$Season <-as.factor(ifelse(telem.sf$Month < 4 | telem.sf$Month > 9, "Non-Breeding",
                                 ifelse(telem.sf$Month > 3 | telem.sf$Month < 10, "Breeding", NA)))

table(telem.sf$Month, telem.sf$Season) # check to make sure Seasons are pulling correct dates
# pulling dates correctly, note that telemetry data is minimal in non-breeding season

###--- check number of fixes per animal / day / etc
# Check the multiple counts of animals per day
counts.per.day <- telem.sf %>% group_by(AnimalID,Day.j) %>% 
  summarise(total = n(),unique = unique(Day.j)) %>% 
  group_by(AnimalID,total) %>% 
  summarise(total.j = n()) 

pp = ggplot(data = counts.per.day, aes(total,total.j,col=AnimalID)) + 
  geom_point() + 
  ggtitle('Number of fixes per julien day per animal'); pp 
# only 1 fix per day


## Q : How many fixes do we have per animal?
p1 <- ggplot(telem.sf) + 
  geom_bar(aes(AnimalID)) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  facet_wrap(~Year, ncol=1)  ; p1

unique(telem.sf$AnimalID) # 29 animals
table(telem.sf$AnimalID, telem.sf$Year) # 4 with data for both years

## Q:  What about temporal variability - number of years? time of year? 
id.date <- telem.sf %>% group_by(AnimalID, Year, Month) %>% summarise(count = n()) ; id.date
id.season.yr <- telem.sf %>% group_by(AnimalID, Year, Season) %>% summarise(count = n()) ; id.season.yr
id.season <- telem.sf %>% group_by(AnimalID, Season) %>% summarise(count = n()) ; id.season
id.jdate <- telem.sf %>% group_by(AnimalID, Year) %>% summarise(count = length(unique(Day.j))) ; id.jdate 


p2 <- ggplot(id.date,aes(x = as.factor(Year), count)) + 
  geom_bar(stat ="identity") + 
  facet_wrap(~AnimalID) + 
  ggtitle("Total telemetry fixes per animal per year (all months)"); p2

p2.1 <- ggplot(id.jdate, aes(x = as.factor(Year), count)) + 
  geom_bar(stat ="identity") +  
  geom_hline(yintercept=50, linetype="dashed", color = "red") +
  facet_wrap(~AnimalID) + 
  ggtitle("Total telemetry days per animal per year (all months)"); p2.1

p3 <- ggplot(id.date,aes(x = as.factor(Month), count)) + 
  geom_bar(stat ="identity") + 
  facet_wrap(~AnimalID)+ 
  ggtitle("Total telemetry fixes per animal per month (all years)"); p3

p4 <- ggplot(id.season,aes(x = Season, count)) + 
  geom_bar(stat ="identity") +  
  geom_hline(yintercept=50, linetype="dashed", color = "red") +
  facet_wrap(~AnimalID) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  ggtitle("Total telemetry fixes per animal per season "); p4 # red line to highlight above/below 50 fixes


# Q how much data do we have when we combine all years and all animals? 
p5 <- ggplot(id.date,aes(x = Month, count)) + 
  geom_bar(stat ="identity"); p5
p6 <- ggplot(id.season,aes(x = Season,count)) + 
  geom_bar(stat ="identity"); p6
p7 <- ggplot(id.jdate,aes(x = as.factor(Year), count)) + 
  geom_bar(stat ="identity"); p7

##############################################################
#### TELEMETRY DATA EXPLORATION & FORMATTING (END) 
#############################################################

##############################################################
#### MERGE INTO ONE (FULL) SF OBJECT (BEGINNING)
#############################################################
head(anml)
head(telem.sf)

HR.sf <- left_join(telem.sf %>% select(AnimalID, Group.New, Date.Timep, Year, Month, Day.j, Season),
                   anml %>% select(-Cptr_Date, -Rls_Date), by = "AnimalID")

glimpse(HR.sf)
summary(HR.sf) # not all animals in both databases, 412 NAs

# drop levels if no entries
HR.sf$Species <- droplevels(HR.sf$Species)
HR.sf$Sex <- droplevels(HR.sf$Sex)
HR.sf$Age_Class <- droplevels(HR.sf$Age_Class)

summary(HR.sf) # check drop levels - looks good

# plot to check
ggplot() +
  geom_sf(data = HR.sf, aes(fill=Sex, col=Sex))
  
ggplot() +
  geom_sf(data = HR.sf, aes(fill=as.factor(Cptr_Year), col=as.factor(Cptr_Year)))


HR.df <- HR.sf %>% st_drop_geometry() # create non-spatial attribute table for joining
##############################################################
#### MERGE INTO ONE (FULL) SF OBJECT (END)
#############################################################

######################################################################
#### SUBSET DATA TO ANIMAL_SEASON AND ANIMAL_YEAR MIN OBS (BEGINNING) 
#####################################################################
# need to drop animals with less than minimum number of observations
summary(HR.sf)
table(HR.sf$Year, HR.sf$Season)
# Breeding Non-Breeding
# 2017     1071          218
# 2018     1156          191

###--- Calculate MCPs for each animal, annually and seasonally (i.e.,combining years but separating seasons)
# group into seasons for each animal 
HR.sf$Animal_Season <- as.factor(paste(HR.sf$AnimalID, HR.sf$Season, sep="_"))
HR.sf$Animal_Year <- as.factor(paste(HR.sf$AnimalID, as.factor(HR.sf$Year), sep="_"))

as.data.frame(HR.sf %>% group_by(AnimalID) %>% count(Season, sort=TRUE)) 
# 43 unique animal-seasons but only 24 with >= 50 obs, 29 with >= 25 obs
as.data.frame(HR.sf %>% group_by(AnimalID) %>% count(Year, sort=TRUE)) 
# 34 unique animal-years but only 20 with >= 50 obs, 24 with >= 25 obs

###--- DECISION POINT - how many obs minimum for HR estimate?
# remove animals with < 25 points per Animal_Season
head(HR.sf)
HR.sf.AS <- HR.sf
HR.sf.AS <- HR.sf.AS[HR.sf.AS$Animal_Season %in% names(table(HR.sf.AS$Animal_Season)) [table(HR.sf.AS$Animal_Season) >= 25], ]
HR.sf.AS$Animal_Season <- droplevels(HR.sf.AS$Animal_Season)
nrow(HR.sf) - nrow(HR.sf.AS) # dropped 145 points

# remove animals with < 25 points per Animal_Year
HR.sf.AY <- HR.sf
HR.sf.AY <- HR.sf.AY[HR.sf.AY$Animal_Year %in% names(table(HR.sf.AY$Animal_Year)) [table(HR.sf.AY$Animal_Year) >= 25], ]
HR.sf.AY$Animal_Year <- droplevels(HR.sf.AY$Animal_Year)
nrow(HR.sf) - nrow(HR.sf.AY) # dropped 108 points

as.data.frame(HR.sf.AS %>% group_by(AnimalID) %>% count(Season, sort=TRUE)) # 0 animal-seasons below 25 obs; 29 unique animal-seasons
as.data.frame(HR.sf.AY %>% group_by(AnimalID) %>% count(Year, sort=TRUE)) # 0 animal-seasons below 25 obs; 24 unique animal-seasons

######################################################################
#### SUBSET DATA TO ANIMAL_SEASON AND ANIMAL_YEAR MIN OBS (END) 
#####################################################################

##############################################################
#### HOUSEKEEPING (BEGINNING)
#############################################################

###-- Save workspace and move to home range analysis - 01_MCP.R
setwd(InputDir)
save.image("00_TelemDataPrep.RData")
# load("00_TelemDataPrep.RData")

###--- For housekeeping, remove all but necessary objects to use in next script
rm(list=setdiff(ls(), c("telem.sf","HR.sf", "HR.df","anml", "Cpt.sf", "HR.sf.AS", "HR.sf.AY")))
save.image("HR_InputData.RData")

##############################################################
#### HOUSEKEEPING (END)
#############################################################

