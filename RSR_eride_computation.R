# Load required packages
library(raster)
library(landscapemetrics)
require(tidyverse)
library(doParallel)

setwd("D:/PROJECTS/OneHealth/Asia_eride_roads")

# Load Rhinolophid bat diversity x forest raster layer, should be UTM projected (m units)
gc <- raster("D:/PROJECTS/OneHealth/Final results/Current/hide/Rhinolophids_SDM_sum_all_win_forfrags_UTM.tif") #1km version, TEAsia

# Convert to binary
gc.bin <- gc 
gc.bin[gc>0] <- 1

# Load roads layer
roads <- raster("D:/PROJECTS/OneHealth/Infrastructure/Final_Roads_All_UTM_resampled.tif")

# plot to check
windows()
par(mfrow=c(1,1))
plot(gc)
plot(gc.bin)
plot(roads)

# Crop forest layer according to extent of roads layer
gc.bin <- crop(gc.bin, roads)


### Fragmenting habitat layer with roads layer ####

# Convert roads to binary first
roads2 <- roads
roads2[roads>=1] <- 1
roads2[is.na(roads)] <- 0
plot(roads2)

writeRaster(roads2, filename = "D:/PROJECTS/OneHealth/Infrastructure/Final_Roads_All_UTM_resampled_binary.tif", overwrite = T)

# subtract roads layer from habitat layer
gc.bin2 <- gc.bin - roads2
gc.bin2[gc.bin<=0] <- 0
plot(gc.bin2)

writeRaster(gc.bin2, filename = "Rhinolophids_forfrags_roadfragd_utm.tif", overwrite = T)
#gc.bin2 <- raster("Rhinolophids_forfrags_roadfragd_utm.tif")

#check if can be analysed in landscape metrics
check_landscape(gc.bin2)

# define forest patches with 4-pixel connectivity
gc_frags <- get_patches(gc.bin2, class = "1", directions = 4, return_raster = T)

# Plot forest patches
plot(gc_frags$layer_1$class_1)

# write it out
writeRaster(gc_frags$layer_1$class_1, filename = "TEA_ch25m_1km_frags_roaded.tif", overwrite = T)

# convert non-forest classes to 0 for patch stats computation
gc3 <- gc.bin2
gc3[gc.bin2==0] <- NA
plot(gc3)

# check if can be analysed in landscape metrics
check_landscape(gc3)

# get patch areas 
patch_areas <- lsm_p_area(gc3, directions = 4)

# get patch perimeters
patch_perims <- lsm_p_perim(gc3, directions = 4)

# bind area and perimeter metric with the patch IDs
pids <- gc_frags$layer_1$class_1$layer
pids.df <- raster::as.data.frame(pids, xy=T, df=T, na.rm=T)
head(pids.df)
colnames(pids.df)[3] <- "id" #1km


pids.df2 <- left_join(pids.df, patch_areas, by = "id") 
pids.df2 <- left_join(pids.df2, patch_perims, by = "id")
head(pids.df2)

pids.df2 <- dplyr::select(pids.df2, x, y, id, value.x, value.y)
colnames(pids.df2) <- c("x", "y", "patchid", "area", "perim") 


#compute biodiversity values
pids.df2$biodiv <- (pids.df2$area)^0.28 
pids.df2$biodivlog <- log(pids.df2$biodiv)
head(pids.df2)
write.csv(pids.df2, "TEA_ch25m_1km_fragstats_roaded.csv", row.names = F)

summary(pids.df2) 


# rasterize 

pids.biodiv <- rasterize(x = pids.df2[,1:2], 
                         y = pids, 
                         field = pids.df2$biodiv, 
                         fun = mean, 
                         background = 0)

pids.biodivlog <- log(pids.biodiv)

par(mfrow=c(1,2))
plot(pids)
plot(pids.biodiv, main = "species div")
plot(pids.biodivlog, main = "species div log")


#save
writeRaster(pids, filename = "TEA_ch25m_1km_frags_nums_roaded.tif", overwrite = T)
writeRaster(pids.biodiv, filename = "TEA_ch25m_1km_frags_spdiv_roaded.tif", overwrite = T)
writeRaster(pids.biodivlog, filename = "TEA_ch25m_1km_frags_spdivlog_roaded.tif", overwrite = T)




#### determining patch edge pixels for each patch ###
# create a focal function to check for edge cells
edge_fun <- function(x) {
  if (any(x == 0 | is.na(x))) {
    return(1)
  } else {
    return(0)
  }
}

# load binary raster
r <- pids.biodiv
r[pids.biodiv>0] <- 1


# apply the focal function to the input raster
edge_raster <- focal(r, w = matrix(1, nrow = 3, ncol = 3), fun = edge_fun, pad = TRUE, silent = F, na.policy = "omit")

# remove NA values from the edge raster
edge_raster <- mask(edge_raster, r)

#raster calculator to find the edges
edge_final <- raster(edge_raster)
values(edge_final) <- 0
edge_final <- edge_raster + r

edge_final[edge_final == 1] <- 0 
edge_final[edge_final == 2] <- 1

edge_final[is.na(edge_final)] <- 0
#edge_final <- mask(edge_final, r)

par(mfrow=c(1,2))
plot(r); plot(edge_final)

#write it out
writeRaster(edge_final, filename = "TEA_ch25m_1km_frags_edges2_roaded.tif", overwrite = T )

### count the number of edge pixels for each patch ###

# calculate the sum of values in the other raster within the zones defined by the boundary raster
patchedge_sum_frag <- zonal(edge_final, pids, fun = "sum", na.rm = TRUE)
head(patchedge_sum_frag)
colnames(patchedge_sum_frag) <- c("patchid", "edgepxsum")

pids.df2 <- left_join(pids.df2, patchedge_sum_frag, by = "patchid", `copy` = TRUE)
head(pids.df2)

#convert to raster

pids.edgepxsum <- rasterize(x = pids.df2[,1:2], 
                            y = pids, 
                            field = pids.df2$edgepxsum, 
                            fun = mean)
plot(pids.edgepxsum)

writeRaster(pids.edgepxsum, filename = "TEA_ch25m_1km_frags_edgepxsum_roaded.tif")


### eRIDE index: sum of the product of the number of edge px from each patch x Rhinolophids patch biodiversity ###
### counting the number of patch pixels within the 20x20 window (nrow,ncol = 21)

patchedge_px_sum <- focal(edge_final, w = matrix(1, nrow = 21, ncol = 21), fun = sum, na.policy = "omit")
plot(patchedge_px_sum)

writeRaster(patchedge_px_sum, filename = "TEA_1km_frags_patchedgepxsumfocal_roaded.tif", overwrite=T)


#eride1 = direct multiplication edge pixel focal sums and biodiv layer 
gc.x <- gc.bin2 * gc
plot(gc.x)
writeRaster(gc.x, filename = "Rhinolophids_SDM_sum_all_win_forfrags_roaded_UTM.tif", overwrite = T)

eride_px <- patchedge_px_sum*gc.x
eride_px_log <- log(eride_px)


par(mfrow=c(1,2))
plot(eride_px, main = "eRIDE1")
plot(eride_px_log, main = "eRIDE1 log")



writeRaster(eride_px, filename = "TEA_ch25m_1km_frags_eride1_roaded.tif", overwrite = T)
writeRaster(eride_px_log, filename = "TEA_ch25m_1km_frags_eridelog1_roaded.tif", overwrite = T)



#read back in..
#eride_px <- raster("TEA_ch25m_1km_frags_eride1.tif")


## Population at Risk (PAR) ####

#read in population density raster
afrpop <- raster("D:/PROJECTS/OneHealth/Asia_eride/WORLDPOP_ASIA_PD_2020_TEA_UTM.tif")
plot(afrpop)

afrpop2 <- crop(afrpop, eride_px)
afrpop2 <- resample(afrpop2, eride_px, method = "bilinear")
plot(afrpop2)


poprisk1 <- eride_px*afrpop2
poprisklog1 <- log10(poprisk1)
poprisklog1[poprisklog1 < 0] <- 0

par(mfrow=c(1,2))
plot(log(afrpop2), main = "Population Density (2020)")
plot(poprisklog1, col = heat.colors(30), main = "log10(population at risk)")

writeRaster(poprisk1, filename = "TEA_1km_frags_poprisk1_roaded.tif", overwrite = T)
writeRaster(poprisklog1, filename = "TEA_1km_frags_poprisklog1_roaded.tif", overwrite = T)


##### END OF CODE #####

