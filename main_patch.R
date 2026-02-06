############################################################################################
## 
##  Queru et al. Better modelling habitat patches // Project ERC SCALED - 949812
##  Main script to perform cutting procedure on habitat contiguity-based patches - assessing initial patches, 
## assessing cutting parameters, cutting patches, creating the final patch layer  
##
############################################################################################
############################################################################################

### Load packages
library(terra)
library(raster)
library(deldir)
library(RColorBrewer)
library("rstudioapi")


### import functions
source("fct_dilatation_erosion.R")
source("fct_tesselation_and_cutting.R")
source("fct_situation_report_patches.R")
source("fct_patchs_analysis.R")
source("fct_patch_cut_proc.R")

######################################################################################################################
#### PREPARE DATA AND SET UP PARAMETERS
######################################################################################################################

### Prepare habitat map // IF NEEDED // Example with a thresholded layer of forest cover density, to map habitat for the European red squirrel
habitat_threshold = 30
treecoverdensity = rast("inputs/forest_density100.tif")
treecoverdensity_th = treecoverdensity>habitat_threshold
study_area = vect("inputs/study_area.shp")
habitat = crop(treecoverdensity_th,study_area)
writeRaster(habitat,file=("inputs/habitat_Alpilles.tif"),overwrite=TRUE)
  
### OR Load habitat map ## Test data on Alpilles landscape ## 
# Habitat Layer (must be a binary raster or a text file, containing only 0s and 1s)
# 0 = matrix pixel ; 1 = Habitat pixel
HabitatLayer = rast("inputs/habitat_Alpilles.tif")
plot(HabitatLayer)

### If necessary add fragmenting elements on top, e.g. roads and rivers, as shapefile
# Has to be shapefile, cropped to fit the study area limits
# If you don't have / need fragmenting elements, please do not set up this parameter
Road = vect("inputs/roads.shp")
River = vect("inputs/rivers.shp")
#crs(river, warn=FALSE)<-crs(roads)
ListFragmentingElements = list(Road, River)

### OPTIONAL STEP - Perform a DILATATION EROSION on habitat map, to remove the small gaps and spurs
# Do you want to perform a dilatation-erosion step? If yes, set TRUE. This additional step allows to fill small gaps in habitat
# Default value is FALSE
DilatationErosionChoice = F # TRUE
# Set here the threshold value for dilatation - erosion, ie the maximum distance from habitat to be dilated (m)
# e.g : threshold = 100, resolution = 100, dilatation will be performed on one pixel)
# DilatationThreshold has to be a strictly positive integer
DilatationThreshold = 100 #m

### Set up DESIRED PATCH SIZE RANGE   
# Set here the minimum and maximum boundaries for patch definition, ie the minimum surface
# for a habitat zone to be considered as a patch (smaller areas can be considered as 
# stepping stones) - example here given for the European red squirrel
MiniArea = 22500 #m²
MaxiArea = 125000 #m²

# AESTHETIC PREFERENCES
# Would you like to display the intermediate figures : set "TRUE" if you do
# Default value is FALSE
Display = T

# How many colors do you want for displaying final habitat patches in the map plots : Must be between 3 and 12
# Default value is 12
NbColorsWanted = 12

######################################################################################################################
#### EVALUATE THE HABITAT LAYER - define how many patches (based on habitat contiguity definition) fall within the desired patch size range
#### AND DISTRIBUTE THE INITIAL PATCHES INTO THREE SUB-LAYERS : too small, correct size and too large patches
######################################################################################################################

## USE THE FUNCTION
situation_report_patches (habitat_layer = HabitatLayer,
                          list_fragmenting_elements = ListFragmentingElements,
                          dilatation_erosion_choice = DilatationErosionChoice,
                          dilatation_threshold = DilatationThreshold,
                          mini_area = MiniArea,
                          maxi_area = MaxiArea,
                          nb_colors_wanted = NbColorsWanted,
                          display = Display, 
                          log_scale = T)


# too-small patches are saved at (0,1 layer):  
TooSmallPatches = rast("patches_too_small.tif")
# correct-size patches are saved at (stored with a patch number):  
CorrectPatches = rast("patches_correct_size.tif")
# too-large patches are saved at (0,1 layer):  
TooLargePatches = rast("patches_too_large.tif")

######################################################################################################################
#### DEFINE MAIN PARAMETERS TO CUT THE TOO LARGE PATCHES AND 
#### CHOOSE THE NUMBER OF POINTS TO BE DRAWN FOR WITHIN-PATCH POLYGON DELIMITATION (higher number of points mean more numerous and smaller patches)
######################################################################################################################

### TESSELATION PARAMETERS
# Do you want to draw random ("random") or regular ("regular") points for tesselation within too oversized contiguity-based patches?
Type = "random"

# Regarding the optimal points to draw, the default option is to maximize the percentage of patches (once cut) that fall within
# the desired path size range (option "range")  
# An alternative option is to minimize the difference between the mean patch size (once cut) and a mean desired patch size (option "goal")
OptimalPointNumberMeth = "range"
# If "range" was chosen, please define the accuracy of the search with parameter tol, min is 1, default is 10, 
# tol should be quite low to have an accurate value for the optimal number of points
Tol = 5
# If "goal" has been chosen, please define how close you want to be to the goal
# ie. if you accept a deviation of 20%, set 0.2 ..., default = 0.1
DistanceToGoal = 0.1

### FIND THE RIGHT NB OF POINTS TO DRAW
# Percent of patches within the desired range will form a bell-shape with increasing number of points, first increasing, then decreasing after having reached a maximum
# bounds_min and bounds_max are to be played with to ensure the full range of option is tested around what could be a first approximation: mean_patch_size_of_too_large_patches/mean_wanted_patch_size
# if percent of correctly ranked patches only increases, increase bounds_max; if it only decreases, decrease bounds_min; otherwise test a large range: 0.1 - 4 for instance
# bounds_min cannot be too small, otherwise the number of points to be drawn will be 0 (the function corrects for that)
# bounds_min cannot be too large, otherwise the range explored will be too large and the optimum might be missed, a warning is returned if the number of points to be drawn exceeds 
# the number of available pixels
# this choice may depend on Type, bounds_max can be larger for "regular" than "random" 
 
# length.out: number of points numbers for which the desired size range is calculated to search for the optimum - numbers of points are along a regular sequence between 
# bounds_min*theoretical_number and bounds_max*theoretical_number, default is length.out = 10, it should be an integer. Note that increasing length.out will increase the 
# accuracy of the search, but will also increase calculation time

# rep defines how many times the procedure should be repeated (drawing N points, cutting, and assessing percent of patches that fall within
# the desired patch size range) for each number of points. This process includes some stochasticity, so each repetition will be a bit different, 
# and increasing rep will help infer this potential variability. Yet, increasing rep will also lead to higher calculation time; default is 1

## USE THE FUNCTION
NbOpt = find_random_nb_opt_based_on_range(too_large_p = TooLargePatches,
                                  mini_area = MiniArea, 
                                  maxi_area = MaxiArea, 
                                  type = Type,
                                  tol = Tol, 
                                  bounds_min=0.1, 
                                  bounds_max=4, 
                                  length.out=10, 
                                  rep = 1)


######################################################################################################################
#### CREATE THE FINAL PATCH LAYER
######################################################################################################################

patch_cut_proc(too_large_p = TooLargePatches,
          correct_p = CorrectPatches,
          nb_opt = NbOpt,
          type = Type,
          optimal = OptimalPointNumberMeth,
          distance_to_goal= DistanceToGoal,
          tol = Tol,
          nb_colors_wanted = NbColorsWanted,
          display = Display) 

# final patches are saved at   
FinalPatches = rast("final_patches.tif")

COL = brewer.pal(n = NbColorsWanted, name = "Set3")[rep(sample(1:NbColorsWanted,NbColorsWanted,replace=F),ceiling(minmax(FinalPatches, compute=TRUE)[2]/NbColorsWanted))] 
plot(FinalPatches,col=COL)
plot(Road,add=T)

# check characteristics of new patch layer
assess_patch_size(FinalPatches,mini_area=MiniArea,maxi_area=MaxiArea)
plot_histo(FinalPatches,log_scale=T); abline(v=log(MiniArea),col="red");abline(v=log(MaxiArea),col="red") 

### Create a resistance layer from final habitat patches and info on surrounding matrix  
# reclassify final patches as habitat
#FinalPatchesRec = classify(FinalPatches, matrix(c(0,99999,1,NA,NA,100), ncol=3, byrow=TRUE), include.lowest=TRUE)
#TooSmallPatchesRec = classify(TooSmallPatches, matrix(c(0,99999,10,NA,NA,1), ncol=3, byrow=TRUE), include.lowest=TRUE) 
# give a breadth to roads
#RoadBuffer = buffer(Road,40)
#RoadBufferRast = rasterize(RoadBuffer,FinalPatches,field="LONGUEUR")
#RoadBufferRastRec = classify(RoadBufferRast, matrix(c(0, 99999, 900, NA, NA, 0), ncol=3, byrow=TRUE), include.lowest=TRUE)
#plot(RoadBufferRastRec)
# combine
#CostLayer = FinalPatchesRec+RoadBufferRastRec
#CostLayer = CostLayer/TooSmallPatchesRec 
#writeRaster(CostLayer, "cost_layer.tif", overwrite=TRUE)

 


