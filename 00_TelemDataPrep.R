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

GISDir <- "//spatialfiles.bcgov/work/wlap/sry/Workarea/jburgar"

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
library(units)

# source function
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

##############################################################
#### METADATA EXPLORATION & FORMATTING (BEGINNING) 
#############################################################

###--- upload animal and collar metadatda files
anml.full <- read.csv("Input/PEPE_Metadata.csv", # point to appropriate metadata file
                     header = TRUE, stringsAsFactors = TRUE, na.strings=c("", "NA"),)

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

# create a 100 km2 buffer around capture locations as study area
Cpt.sf.buff <- Cpt.sf %>% summarise( geometry = st_combine( geometry ) ) %>% 
  st_convex_hull() %>% st_buffer(dist=100000)
aoi <- Cpt.sf.buff

# plot to check
bc <- bc_bound()

ggplot() +
  geom_sf(data = bc) +
  geom_sf(data = Cpt.sf.buff)+
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


MMP_BCGov <- st_read(dsn=paste0(GISDir,"/MMP"), layer="MMP_BCGov_grid_March2024") 

ggplot() +
  geom_sf(data = bc) +
  geom_sf(data=MMP_BCGov) +
  geom_sf(data=Cpt.sf, aes(col=Species))

# FHE <- st_read(dsn=paste0(GISDir,"/Fisher"), layer="FHE_two_pops")
# FHE <- st_simplify(FHE)
# FCol <- FHE %>% filter(Hab_zone!="Boreal")
# Created a list of traplines that intersect the Columbian fisher range
# bcdc_search("trapline", res_format = "wms")
# traplines <- bcdc_query_geodata("f8e27889-fb07-4f9d-b2ba-591578274b7c") %>% collect()
# Col.trapline <- traplines %>% st_intersection(FCol)
# write.csv(Col.trapline %>% select(TRAPLINE_AREA_IDENTIFIER) %>% st_drop_geometry, file="./Traplines_ColRange.csv", row.names = FALSE)

# using the bc data warehouse option to clip to aoi
aoi <- aoi %>% st_transform(3005)

# biogeoclimatic zones
# bcdc_search("Biogeoclimatic zone", res_format = "wms")
aoi.BEC <- retrieve_geodata_aoi(ID = "f358a53b-ffde-4830-a325-a5a03ff672c3")
aoi.BEC %>% group_by(ZONE) %>% summarise(sum(Area_km2)) %>% st_drop_geometry()

# approved WHAs & UWRs (GAR Orders)
# bcdc_search("WHA", res_format = "wms")
aoi.WHA <- retrieve_geodata_aoi(ID = "b19ff409-ef71-4476-924e-b3bcf26a0127")
aoi.WHA %>% group_by(COMMON_SPECIES_NAME) %>% summarise(sum(Area_km2)) %>% st_drop_geometry()

# bcdc_search("UWR", res_format = "wms")
aoi.UWR <- retrieve_geodata_aoi(ID = "712bd887-7763-4ed3-be46-cdaca5640cc1")
aoi.UWR %>% group_by(SPECIES_1) %>% summarise(sum(Area_km2)) %>% st_drop_geometry()

# bcdc_search("Forest District", res_format = "wms")
# 2: Natural Resource (NR) Districts (multiple, wms, oracle_sde)
# ID: 0bc73892-e41f-41d0-8d8e-828c16139337
# Name: natural-resource-nr-district
aoi.NRD <- retrieve_geodata_aoi(ID = "0bc73892-e41f-41d0-8d8e-828c16139337")
aoi.NRD %>% group_by(DISTRICT_NAME) %>% summarise(sum(Area_km2)) %>% st_drop_geometry()

##############################################################
#### METADATA EXPLORATION & FORMATTING (END) 
#############################################################

##############################################################
#### TELEMETRY DATA EXPLORATION & FORMATTING (BEGINNING) 
#############################################################

###--- For csv upload
# list.files("Input")

telem <- fs::dir_ls(path="Input",regexp = "^Input/GPS", recurse = TRUE) %>%
  map_dfr(read_csv, col_types = cols(.default = 'c'), .id = "source") %>% type.convert()

telem$AnimalID <- as.factor(word(telem$source, 2, sep="\\_"))
telem$Date_PST <- ymd(telem$`RTC-date`, tz=tz) # RTC is real-time clock that we programmed
telem$Time_PST <- hms(telem$`RTC-time`)
telem$DateTime_PST <- ymd_hms(paste(telem$Date_PST, telem$Time_PST, tz=tz))


glimpse(telem) # check columns are coming in as appropriate class
head(telem) # check that data is reading correctly

#############################################################
###--- clean data
# remove telem locations outside of range when on animals
anml %>% group_by(AnimalID) %>% summarise(min(Rls_Datep))


telem <- telem %>% mutate(time.keep = case_when(AnimalID == "35060" & Date_PST < '2024-02-19 PST' ~ "remove",
                                       AnimalID == "35060" & Date_PST > '2024-04-01 PST' ~ "remove",
                                       AnimalID == "94092" & Date_PST < '2024-02-17 PST' ~ "remove",
                                       AnimalID == "94092" & Date_PST > '2024-04-01 PST' ~ "remove",
                                       TRUE ~ "keep"))

telem %>% group_by(time.keep, AnimalID) %>% summarise(min(Date_PST), max(Date_PST))

telem <- telem %>% filter(time.keep=="keep")

validGPS <- telem %>% group_by(AnimalID) %>% count(Status)
2915/(1213+2915); 3754/(495+3754)
# 35060 had 71% of locations as valid, 94092 had 88% valid locations (within the time period in the field)

telem <- telem %>% filter(Status=="Valid")

###--- check data quality and remove any objectionable rows

# latitude coordinate off - delete that row
tail(telem[order(telem$Latitude),])
tail(telem) # only 2 latitude coordinates off (one coord for each animal) 
# otherwise data was cleaned prior to uploading and is good to go
# use the CRS for lat/long and WGS84
telem %>% summarise(min(Longitude))
telem.sf <- st_as_sf(telem %>% filter(Latitude < 58) %>% filter(Latitude > 50), coords=c("Longitude","Latitude"), crs=4326) 

ggplot()+
  geom_sf(data=telem.sf, aes(col=AnimalID))

# a few points are off
glimpse(telem)
prelim <- telem %>% filter(Sats %in% c("04-Apr","05-May","06-Jun","07-Jul","08-Aug")) %>%
  filter(Latitude < 51.9) %>% filter(Longitude < -121.5)
prelim.sf <- st_as_sf(prelim, coords=c("Longitude","Latitude"), crs=4326) 


ggplot()+
  geom_sf(data=prelim.sf, aes(col=AnimalID))

# format dates for R
telem.sf <- prelim.sf
telem.sf$Year <- year(telem.sf$DateTime_PST)
telem.sf$Month <- month(telem.sf$DateTime_PST)
telem.sf$Day.j <- julian(telem.sf$DateTime_PST) # Julian day

telem.sf <- telem.sf %>% filter(!is.na(Month))

ggplot() +
  geom_sf(data = telem.sf, aes(fill=AnimalID, col=AnimalID)) +
  facet_wrap(~Year) +
  theme(legend.position = "bottom") +
  coord_sf() 

# looks good - all plotting in general area


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
glimpse(anml)
anml$AnimalID <- as.factor(anml$AnimalID)
head(telem.sf)

HR.sf <- left_join(telem.sf %>% select(AnimalID, DateTime_PST, Year, Month, Day.j),
                   anml %>% select(-Cptr_Date, -Rls_Date), by = "AnimalID")

glimpse(HR.sf)
summary(HR.sf)

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
#### SUBSET DATA TO ANIMAL_YEAR MIN OBS (BEGINNING) 
#####################################################################
# need to drop animals with less than minimum number of observations
summary(HR.sf)
table(HR.sf$Year, HR.sf$Sex)

###--- Calculate MCPs for each animal, annually and seasonally (i.e.,combining years but separating seasons)
# group into seasons for each animal 
HR.sf$Animal_Year <- as.factor(paste(HR.sf$AnimalID, as.factor(HR.sf$Year), sep="_"))

as.data.frame(HR.sf %>% group_by(AnimalID) %>% count(Year, sort=TRUE)) 
# 2 unique animal-years and each with sufficient n

###--- DECISION POINT - how many obs minimum for HR estimate?
# remove animals with < 25 points per Animal_Year
HR.sf.AY <- HR.sf
HR.sf.AY <- HR.sf.AY[HR.sf.AY$Animal_Year %in% names(table(HR.sf.AY$Animal_Year)) [table(HR.sf.AY$Animal_Year) >= 25], ]
HR.sf.AY$Animal_Year <- droplevels(HR.sf.AY$Animal_Year)
nrow(HR.sf) - nrow(HR.sf.AY) # did not remove any  points

######################################################################
#### SUBSET DATA TO ANIMAL_YEAR MIN OBS (END) 
#####################################################################

##############################################################
#### HOUSEKEEPING (BEGINNING)
#############################################################

###-- Save workspace and move to home range analysis - 01_MCP.R
save.image("00_TelemDataPrep.RData")
# load("00_TelemDataPrep.RData")

###--- For housekeeping, remove all but necessary objects to use in next script
rm(list=setdiff(ls(), c("telem.sf","HR.sf", "HR.df","anml", "Cpt.sf", "HR.sf.AS", "HR.sf.AY")))
save.image("HR_InputData.RData")

##############################################################
#### HOUSEKEEPING (END)
#############################################################

