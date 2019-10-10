###################################################################################################
# 02_KDE_advanced.R
# script to run advanced KDEs (barriers and Brownian Bridge)
# written by Joanna Burgar (Joanna.Burgar@gov.bc.ca) - 06-Aug-2019
# still in beta phase and needs tweaking as data dependent
# also need to add in better code for determining distance to nearest point within same sf object
###################################################################################################
.libPaths("C:/Program Files/R/R-3.6.0/library")

# overall process: 
#- Explore data to see what we have 
#- Define the area of interest per animal; availability vs use (Kernel Density Estimate)
#- Export KDE shapefiles

# help files: 
#https://cran.r-project.org/web/packages/adehabitatHS/adehabitatHS.pdf

# instal packages and run libraries

# run libraries
library(bcmaps)
library(bcdata)
library(lattice)  
library(ggplot2)
library(scales)
library(dplyr)
library(GGally)
library(lubridate)
library(raster)
library(smoothr)
library(sf)
library(sp)
library(adehabitatHR)
library(ggmap)
library(scales)
library(ggmap)
library(tidyr)
library(Cairo)

# set up working directories on H drive
InputDir <- c("H:/R/Analysis/Generic_HomeRange/Input")
OutputDir <- c("H:/R/Analysis/Generic_HomeRange/Output")
GISDir <- c("H:/R/Analysis/Generic_HomeRange/GISDir")

###################################################################################################
#### LOAD, REVIEW DATA AND FORMATE FOR KDE (BEGINNING)
###################################################################################################
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
unique(HR.sf$AnimalID) # individuals
ggplot() +
  geom_sf(data=SC, fill="white", col="gray") +
  geom_sf(data=HR.sf, aes(fill=AnimalID, col=AnimalID))+
  coord_sf() +
  theme_minimal() +
  ggtitle("Animal GPS locations")

###--- Convert to sp for all data
st_geometry(HR.sf)
HR.sf_utm <- st_transform(HR.sf, crs=26910) # utm and m units

HR.sf_utm$geometry
HR.sp <- as(HR.sf_utm, "Spatial")
class(HR.sp)

###--- ANIMAL_SEASON 
st_geometry(HR.sf.AS)
HR.sf.AS_utm <- st_transform(HR.sf.AS, crs=26910) # utm and m units

HR.sf.AS_utm$geometry
HR.sp.AS <- as(HR.sf.AS_utm, "Spatial")
class(HR.sp.AS)

###--- ANIMAL_YEAR 
st_geometry(HR.sf.AY)
HR.sf.AY_utm <- st_transform(HR.sf.AY, crs=26910) # utm and m units

HR.sf.AY_utm$geometry
HR.sp.AY <- as(HR.sf.AY_utm, "Spatial")
class(HR.sp.AY)

###################################################################################################
#### LOAD, REVIEW DATA AND FORMATE FOR KDE (END)
###################################################################################################


###################################################################################################
#### BRING IN SPATIAL BARRIER DATA (BEGINNING)
##################################################################################################
###--- considering constraints in adehabitatHR, best to go with simple, non-tortuous barrier boundaries. 
# For this example, bringing in elevation data and considering animals don't go above 1500 m elev
# bring in layers from internal folders
# study area comprises mapsheets 92 F,G,H,J,K
# \\imagefiles.bcgov\imagery\dem\elevation\trim_25m\utm\esri_ascii_grid
# raster file names variation of "92h-utm-elevation.asc"
DEMDir <- c("//imagefiles.bcgov/imagery/dem/") # can also use bcdata to pull in barrier data, depending on barrier
setwd(DEMDir)

list.files(DEMDir, pattern='\\.asc$', recursive = TRUE)

###--- First attempt - create boundary at 1500 m using contours converted from DEM
# Load in elevation asc (raster) files for relevant mapsheets
elev.92fr <- raster("./elevation/trim_25m/utm/esri_ascii_grid/92f-utm-elevation.asc")
elev.92gr <- raster("./elevation/trim_25m/utm/esri_ascii_grid/92g-utm-elevation.asc")
elev.92hr <- raster("./elevation/trim_25m/utm/esri_ascii_grid/92h-utm-elevation.asc")
elev.92jr <- raster("./elevation/trim_25m/utm/esri_ascii_grid/92j-utm-elevation.asc")
elev.92kr <- raster("./elevation/trim_25m/utm/esri_ascii_grid/92k-utm-elevation.asc")

# plot to make sure loading correctly
plot(elev.92hr,elev.92gr)
str(elev.92hr)

# convert raster to contours and merge sp objects to create one file
elev.levels <- as.vector(c(500,1000,1500,2000,2500))

elev.92fr_cntrs <- rasterToContour(elev.92fr, levels=elev.levels)
elev.92gr_cntrs <- rasterToContour(elev.92gr, levels=elev.levels)
elev.92hr_cntrs <- rasterToContour(elev.92hr, levels=elev.levels)
elev.92jr_cntrs <- rasterToContour(elev.92jr, levels=elev.levels)
elev.92kr_cntrs <- rasterToContour(elev.92kr, levels=elev.levels)

elev.92_cntrs <- rbind(elev.92fr_cntrs,elev.92gr_cntrs,elev.92hr_cntrs,elev.92jr_cntrs,elev.92kr_cntrs)

# plot to check 
plot(elev.92_cntrs[elev.92_cntrs$level==1000,])
plot(elev.92_cntrs[elev.92_cntrs$level==1500,], add = TRUE, col = "red")
plot(elev.92_cntrs[elev.92_cntrs$level==2000,], add = TRUE, col = "blue")
# looks good

# let's try out 1500 m as the barrier
plot(elev.92_cntrs[elev.92_cntrs$level==1500,])
plot(HR.sp, pch = 19, size = 0.5, col = "red", add = TRUE) 

elev.barrier <- elev.92_cntrs[elev.92_cntrs$level==1500,]
class(elev.barrier)

# will need to smooth the lines as currently don't meet kernelUD() constraints
sf_elev.barrier <- st_as_sf(elev.barrier) # first convert to sf object

# (https://cran.r-project.org/web/packages/smoothr/vignettes/smoothr.html)
elev.barrier_Chaikin <- smooth(sf_elev.barrier, method = "chaikin") # try out Chaikin's corner cutting algorithm 
elev.barrier_Kernel <- smooth(sf_elev.barrier, method = "ksmooth", smoothness = 50) # try out kernel method, increase smoothing paramter
elev.barrier_spline <- smooth(sf_elev.barrier, method = "spline") # try out spline method 

elev.barrier_Kernel <- st_set_crs(elev.barrier_Kernel, 26910)
# problem might be small line lengths - drop all lines below 1st quartile
# on reviewing plots with various lenths omitted it looks like lines are not contiguous
# so removing line lengths below threshold doesn't always remove the appropriate lines
elev.barrier_Kernel$length <- st_length(elev.barrier_Kernel)
summary(elev.barrier_Kernel$length)
elev.barrier_Kernel$length <- as.numeric(elev.barrier_Kernel$length)
# filter out lines that are < 500 km in length
elev.barrier_Kernel2 <- elev.barrier_Kernel %>% filter(length > 500000)

setwd(OutputDir)
Cairo(1000, 800, pointsize = 10,
      file="kde2.1BK_smooth50_sub.png", type="png", bg="white")
ggplot() +
  geom_sf(data = elev.barrier_Kernel2, col = "black") +
  geom_sf(data = HR.sf_utm, col = "red", size = 0.5)
dev.off()

# # plot out contours and animal data points
# for (i in 1:nrow(elev.barrier_Kernel)) {
#   plot(st_geometry(elev.barrier_Kernel[i, ]), col = "grey20", lwd = 2)
# }
# plot(HR.sp, pch = 19, size = 0.5, col = "red", add = TRUE) 


# convert sf objects back to sp for use in kernelUD()
sp_elev.barrier_Chaikin <- as(elev.barrier_Chaikin, "Spatial")
sp_elev.barrier_Kernel <- as(elev.barrier_Kernel, "Spatial")
sp_elev.barrier_spline <- as(elev.barrier_spline, "Spatial")

###--- create UD with barrier - use same paramters as user defined optimal in 02_KDE.R 
# for this example using h = 500, grid = 500, extent = 2 and using sp object Animal_Season as data
kde2.1B  <- kernelUD(HR.sp.AS[c("Animal_Season")],
                     h = 500, kern = c("bivnorm"), grid = 500,extent = 2, 
                     boundary=elev.barrier) # too tortuous, didn't work (turning angles > pi/2)

kde2.1BC  <- kernelUD(HR.sp.AS[c("Animal_Season")],
                     h = 500, kern = c("bivnorm"), grid = 500,extent = 2, 
                     boundary=sp_elev.barrier_Chaikin) # didn't work (turning angles > pi/2) 

kde2.1BK  <- kernelUD(HR.sp.AS[c("Animal_Season")],
                     h = 500, kern = c("bivnorm"), grid = 250,extent = 2, 
                     boundary=sp_elev.barrier_Kernel) # didn't work, even when smoothing parameter increased to 50, small line lengths dropped (those <500 km; although looks like lines might not be contiguous), for either h=250 or h=500

kde2.1BS  <- kernelUD(HR.sp.AS[c("Animal_Season")],
                     h = 500, kern = c("bivnorm"), grid = 500,extent = 2, 
                     boundary=sp_elev.barrier_spline) # didn't work, turning angles > pi/2

###--- housekeeping - remove large elevation files
rm(sp_elev.barrier,sp_elev.barrier_Chaikin, sp_elev.barrier_Kernel, sp_elev.barrier_spline,
   elev.barrier_Chaikin, elev.barrier_Kernel, elev.barrier_spline)

#############################################################################
###--- Second attempt - create boundary at 150 % slope using data converted from DEM
# Load in slope asc (raster) files for relevant mapsheets
slp.92fr <- raster("./slope/trim_25m/percent/utm/esri_ascii_grid/92f-utm-slope-pct.asc")
slp.92gr <- raster("./slope/trim_25m/percent/utm/esri_ascii_grid/92g-utm-slope-pct.asc")
slp.92hr <- raster("./slope/trim_25m/percent/utm/esri_ascii_grid/92h-utm-slope-pct.asc")
slp.92jr <- raster("./slope/trim_25m/percent/utm/esri_ascii_grid/92j-utm-slope-pct.asc")
slp.92kr <- raster("./slope/trim_25m/percent/utm/esri_ascii_grid/92k-utm-slope-pct.asc")

plot(slp.92fr)

# convert raster to contours/lines and merge sp objects to create one file
slope.levels <- as.vector(c(100, 200, 300, 400, 500))

slp.92fr_cntrs <- rasterToContour(slp.92fr, levels=slope.levels) 
plot(slp.92fr_cntrs) # too tortuous, won't work in kernelUD() function
# for this example, recommend omitting area post-hoc

###################################################################################################
#### BRING IN SPATIAL BARRIER DATA (END)
##################################################################################################

###################################################################################################
#### CONSIDERING TIME DEPENDENCE BETWEEN RELOCATIONS (BEGINNING)
##################################################################################################
###--- the Brownian bridge kernel HR, instead of classic kernel HR
###--- two smoothing paramters
# sig1: This parameter controls the width of the \bridge" connecting successive relocations
# The larger it is and the larger the bridge is, related to the speed of the animal (stress on related)
# sig2: This parameter controls the width of the \bumps" added over the relocations
# It is similar to the smoothing parameter h of the classical kernel method (related to the imprecision of the relocations)

# figure out mean distance SD - use sf object
class(HR.sf_utm)
anml.dist <- HR.sf_utm %>% group_by(Animal_Season) %>% st_distance() # return matrix of all possible distances
# not the best method, but not sure how to calculate summary statistics on distances for each animal (need to think through matrix algebra)
dim(anml.dist)

### wouldn't it be the mean between anml.dist[1,] and anml.dist[,1]??? need to rethink this

test <- as.data.frame(anml.dist[,1])
test$Animal_Season <- HR.sf_utm$Animal_Season
dist.AS.smrysts <- as.data.frame(test %>% group_by(Animal_Season) %>% summarise(mean = mean(`anml.dist[, 1]`), sd = sd(`anml.dist[, 1]`)))
mean(dist.AS.smrysts$sd) # mean of all Animal_Season
mean(dist.AS.smrysts$sd/sqrt(nrow(dist.AS.smrysts))) # mean se
# for housekeeping, remove anml.dist
rm(anml.dist)

# The class ltraj is intended to store trajectories of animals. 
# Trajectories of type II correspond to trajectories for which the time is available for each relocation (mainly GPS and radio-tracking). 

# try for one animal first (note - using animal_seasons in test)
# create a trajectory for one animal
# use dist measure of all animal seasons for sig2
# ?as.ltraj # to check help page for details

HR.sp.AS$Datep <- as.POSIXct(strptime(HR.sp.AS$Date.Timep, format = "%Y-%m-%d"))

geom.anml <- geometry(HR.sp.AS)
anml.traj <- as.ltraj(geom.anml@coords, date=HR.sp.AS$Datep, 
                     id=as.character(HR.sp.AS$Animal_Season), burst=HR.sp.AS$Animal_Season, 
                     typeII = TRUE, slsp = c("remove", "missing"))

setwd(OutputDir)
Cairo(1000, 800, pointsize = 10,
      file="anml_trajectory.png", type="png", bg="white")
plot(anml.traj)
dev.off()

# using SE of mean dist for all anml as sig2, just for anml would be 165
# doesn't seem too sensitive and as just for illustrative purposes, staying with 242
#anml_lik <- liker(anml.traj, sig2 = 242, rangesig1 = c(1,10)) # optimum value is 2 and 2.99 (winter, non-winter)
anml_lik <- liker(anml.traj, sig2 = 165, rangesig1 = c(1,10)) # optimum value is 2.10 and 3.04 (winter, non-winter)
anml_lik
#Maximization of the log-likelihood for parameter
#sig1 of brownian bridge
#anml_Non-winter : Sig1 = 2.991 Sig2 = 242 
#anml_Winter : Sig1 = 2 Sig2 = 242 
mean(c(2.991, 2)) # 2.50
mean(c(3.036, 2.0901)) # 2.57

#anml.tata <- kernelbb(anml.traj, sig1 = 2.50, sig2 = 242, grid = 500, extent = 2) # keeping grid and extent consistent with kernel UD
anml.tata <- kernelbb(anml.traj, sig1 = 2.50, sig2 = 165, grid = 500, extent = 2) # keeping grid and extent consistent with kernel UD
#changing sig2 didn't seem to change the UD, perhaps more sensitive to sig1?

image(anml.tata) 

veranml.tata <- getverticeshr(anml.tata,95) ;   veranml.tata.sf<- st_as_sf(veranml.tata)
veranml.tata.sf <- st_set_crs(veranml.tata.sf, 26910)

setwd(OutputDir)
Cairo(1000, 1000, pointsize = 10,
      file="KDBB_anml_Winter_meansig2.png", type="png", bg="white")
plot(st_geometry(ver95.2.1.sf[ver95.2.1.sf$id=="anml_Winter",]),col = "green", 
     xlim = c(436296, 459255), ylim = c(5470160, 5494927))
plot(st_geometry(veranml.tata.sf[veranml.tata.sf$id=="anml_Winter",], 95), add = TRUE, lwd=2)
plot(HR.sp[HR.sp$Animal_Season=="anml_Winter",], pch = 19, size = 0.5, add = TRUE) # too difficult to see UDs if plotted
dev.off()


#try for all animals by animal season
geom.AS <- geometry(HR.sp.AS)
head(geom.AS)
glimpse(geom.AS)

AS.traj <- as.ltraj(geom.AS@coords, date=HR.sp.AS$Datep, 
                   id=as.character(HR.sp.AS$Animal_Season), burst=HR.sp.AS$Animal_Season, 
                   typeII = TRUE, slsp = c("remove", "missing"))
plot(AS.traj)
# sd seems too large, trying with se
# plot not working with multiple animal seasons
lik <- liker(anml.traj, sig2 = 242, rangesig1 = c(0.1,10), plotit = FALSE) 
head(lik)

# create vector listing the sig1 value for each UD
anml_sig1 <- rep(NA, length(lik))
for (i in 1:length(anml_sig1)){ 
  i.sig1 <- lik[[i]]$sig1
  anml_sig1[i] <- i.sig1
}
min(anml_sig1); max(anml_sig1); mean(anml_sig1)

###--- We can now estimate the kernel Brownian bridge home range with the parameters
# Note - running kernelbb on 66 animal_seasons takes 15-30 min #
# attempt 1: sig1=2.58 and sig2=242 # HR estimates too large in many cases
# attempt 2: sig1=0.87 and sig2=242
# neither really worked - if this is something to pursue reccomend finding optimal sig1 and sig2 for each HR and looping through
# rather than using one value for all animals (considering variaibility in movement, doesn't really make sense to have same values for all animals)
anml.tata <- kernelbb(anml.traj, sig1 = 2.58, sig2 = 242, grid = 500, extent = 2) # keeping grid and extent consistent with kernel UD
#anml.tata <- kernelbb(anml.traj, sig1 = 0.87, sig2 = 242, grid = 500, extent = 2) # keeping grid and extent consistent with kernel UD
str(anml.tata$`anml_Breeding`)
#veranml.tata <- getverticeshr(anml.tata,95) ;   veranml.tata.sf<- st_as_sf(veranml.tata, proj4string())

#image(anml.tata) only  useful if plotting 1, otherwise hard to fit in graphics window (66 individual plots)
image(anml.tata$`anml_Breeding`) 

# for all animal_seasons
setwd(OutputDir)

Cairo(1000, 1000, pointsize = 10,
      file="KDBB_sig1mean.png", type="png", bg="white")
plot(st_geometry(ver95.2.1.sf),col = "green")
plot(getverticeshr(anml.tata, 95), add = TRUE, lwd=2)
#plot(HR.sp, pch = 19, size = 0.5, add = TRUE) # too difficult to see UDs if plotted
dev.off()

###################################################################################################
#### CONSIDERING TIME DEPENDENCE BETWEEN RELOCATIONS (END)
##################################################################################################

##################################################################################################
#### HOUSEKEEPING (BEGINNING)
##################################################################################################
###-- Save workspace
setwd(InputDir)
save.image("02_KDE_advanced.RData") # may be too large as contains tata

###--- For housekeeping, save only tata output and then delete from environment
setwd(InputDir)
rm(list=setdiff(ls(), c("anml.tata")))
save.image("anml.tata.RData")
##################################################################################################
#### HOUSEKEEPING (END)
##################################################################################################

