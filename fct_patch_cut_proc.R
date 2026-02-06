############################################################################################

# Queru et al. Better modelling habitat patches // Project ERC SCALED - 949812

############################################################################################
##  This function allows to cut oversized patches using tessellation principle, and mosaic them with patches of the correct size


#PARAMETERS
# too_large_p = too-large-patch raster (layer of the patches that are oversized and need to be cut), in which the tessellation will be performed 
# correct_p = correct-size-patch raster  to be combined with the patches that will result from the cutting procedure
## WARNING, in these rasters areas that are not habitat are coded as NA,  others as 1 or an integer >0 coding for patch id
# nb_opt = number of points to be drawn, can come from a search of the optimal number of points with find_random_nb_opt_based_on_range
# type = c("random","regular"), random or regular selection to perform the tessellation
# optimal = type of search, the default is "range"
# tol = the accuracy wanted when searching for the optimal number of points to be drawn
# distance_to_goal = define how close the optimal needs to be from the goal if optimal == "goal" 
# nb_colors_wanted = number of colours to deplay maps, default is set to 12
# display = TRUE

patch_cut_proc = function(too_large_p, correct_p, nb_opt, type=c("random","regular"),
                      optimal="range", tol=10, distance_to_goal=0.1, 
                      nb_colors_wanted=12, display=TRUE){
                     
  # Assign new id to correct patches
  temp = terra::freq(correct_p)
  rcl = data.frame(from = temp[,2], to = 1:nrow(temp))
  correct_patch = classify(correct_p,rcl)
  id_max = minmax(correct_patch, compute=TRUE)[2]
  COL = brewer.pal(n = nb_colors_wanted, name = "Set3")[rep(sample(1:nb_colors_wanted,nb_colors_wanted,replace=F),ceiling(id_max/nb_colors_wanted))] 
  
  # Cut too large patches into smaller patches
  out_cut = patch_cutting_all(r=too_large_p,n=nb_opt,id=id_max,type="random")

  # Add back habitat that could have been lost during tessellation
  #id_max2 = minmax(out_cut, compute=TRUE)[2]
  #large_patch_no_id = mask(too_large_p,out_cut, inverse=TRUE)
  #large_patch_no_id = patches(large_patch_no_id,directions = 8)#, zeroAsNA=TRUE, allowGaps=FALSE)
  #large_patch_no_id = large_patch_no_id + id_max2
  #plot(large_patch_no_id)

  # Group cut and correct patches into the same layer
  #out_mos = mosaic(out_cut, large_patch_no_id)
  out_mos2 = mosaic(out_cut, correct_patch)
  plot(out_mos2,col = COL)
  
  ##Save the new file
  writeRaster(x = out_mos2, filename = "final_patches.tif",overwrite=TRUE)
}




