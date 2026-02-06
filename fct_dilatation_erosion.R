############################################################################################

# Queru et al. Better modelling habitat patches // Project ERC SCALED - 949812

############################################################################################

# FUNCTION TO PERFORM THE DILATATION - EROSION STEP, MEANING FILLING SMALL GAPS WITHIN HABITAT AREAS

## REQUIRED PARAMETERS : 
# r = raster of habitat layer in binary format 0/1 (0=matrix ; 1=habitat)
# threshold = size (in m) of the maximum gap size allowed to perform the dilation-erosion

dilatation_erosion <- function(r,threshold){

    ### Replace 0 (ie matrix) with NA 
    r[r[[1]]== 0] = NA
    
    ### Calculate distance to habitat areas
    distances = terra::distance(r)
    
    ### distance-to-habitat thresholding (cell nb)
    distances_hab = distances>threshold
    
    ### Replace 0 (ie distance-to-habitat < threshold ) by NAs
    distances_hab[distances_hab[[1]]== 0] <- NA
    
    ### Calculating distance to non-habitat
    distances_non_habitat = terra::distance(distances_hab)
    
    ### distance-to-habitat thresholding (cell nb)
    final_rast = as.numeric(distances_non_habitat>threshold)
    
    return(final_rast)

}


