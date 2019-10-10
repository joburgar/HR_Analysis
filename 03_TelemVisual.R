#############################################################
# 03_TelemVisual.R
# script for visualizing telemetry of collared animals
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 08-Oct-2019
#############################################################

.libPaths("C:/Program Files/R/R-3.6.0/library") # to ensure reading/writing libraries from C drive (H drive too slow)

###--- Load packages
library(dplyr)    # for viewing and manipulating data
library(anipaths) # for making animated graphics
library(sf)       # working with sf objects
library(OpenStreetMap)  # if using Open Street Maps (OSM) data as background map
library(ggmap)   # for plotting background map that uses osmdata


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


###--- ANIMAL_SEASON 
st_geometry(HR.sf.AS)
HR.sf_utm <- st_transform(HR.sf.AS, crs=26910) # utm and m units

head(HR.sf_utm)
st_coordinates(HR.sf_utm)
HR.sf_latlon <- st_transform(HR.sf_utm, 4326) # can be either utm or lat/lon for antipaths, but easier to download as latlon for OSM
summary(HR.sf_latlon)

coords <- as.array(st_coordinates(HR.sf_latlon))
glimpse(coords)
class(coords)
head(coords)
colnames(coords) <- c("Longitude", "Latitude")

HR.ani <- as(HR.sf_latlon, "Spatial")

###--- view OSM data and download appropriate section for study area
st_bbox(HR.sf_latlon)
st_bbox(HR.sf_latlon)[4]

LAT1 = st_bbox(HR.sf_latlon)[2] ; LAT2 = st_bbox(HR.sf_latlon)[4]
LON1 = st_bbox(HR.sf_latlon)[3] ; LON2 = st_bbox(HR.sf_latlon)[1]

#our background map
map <- openmap(c(LAT2,LON1), c(LAT1,LON2), zoom = NULL,
               type = c("osm", "stamen-toner", "stamen-terrain","stamen-watercolor", "esri","esri-topo")[6],
               mergeTiles = TRUE)


## OSM CRS :: "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs"
map.latlon <- openproj(map, projection = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
coord <- st_coordinates(HR.sf_latlon)


#######################################################################################
# for all animals, with background map, on a weekly interval, grouped by "Group.New"
setwd(paste(OutputDir,"/Animations/AS_All", sep=""))
animate_paths(paths = HR.ani, # SpatialPointsDataFrame object
              delta.t = "week",
              coord = coord(c("Longitude","Latitude")), # set up as Longitude, then Latitude
              Time.name = "Date.Timep",
              ID.name = "Animal_Season",
              covariate = "Group.New",
              covariate.colors = RColorBrewer::brewer.pal(n = 4, "RdYlGn"),
              legend.loc = NA, #removes legend - gets too messy with this many animals
              background = map.latlon) # using OpenStreetMaps to avoid google maps API key issues
              
# for all Relocated individuals, with background map, on a daily interval
setwd(paste(OutputDir,"/Animations/AS_Relocated", sep=""))
coord.Relocated <- st_coordinates(HR.sf_latlon[HR.sf_latlon$Group.New=="Relocated",])

animate_paths(paths = HR.ani[HR.ani$Group.New=="Relocated",], # SpatialPointsDataFrame object
              delta.t = "day",
              coord = coord.Relocated(c("Longitude","Latitude")), # set up as Longitude, then Latitude
              Time.name = "Date.Timep",
              ID.name = "Animal_Season",
              #legend.loc = NA, #removes legend - gets too messy with this many animals
              background = map.latlon) # using OpenStreetMaps to avoid google maps API key issues


setwd(InputDir)
save.image("03_TelemVisual.RData")

