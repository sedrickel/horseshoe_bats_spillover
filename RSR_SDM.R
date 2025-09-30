require(tidyverse)
library(raster)
library(rasterVis)
library(rJava)
library(usdm)
library(ENMeval)
library(maps)
library(mapdata)
library(dismo)  
library(maptools)
library(jsonlite)
library(rworldmap)
library(RColorBrewer)
library(doParallel)
library(readxl)
library(tmap)
library(ENMTools)
library(forcats)
library(hrbrthemes)
library(viridis)
library(multcompView)
library(visreg)
library(jtools)
library(caret)
library(CoordinateCleaner)
library(sf)

registerDoParallel()
getDoParWorkers() 

setwd("D:/PROJECTS/OneHealth/Maxent")

### TRUE SKILL STATISTIC TEST FUNCTION #####
#to be able to run this script you need to have told the Maxent model to produce background predictions. If you are running MaxEnt in R this means putting the argument (after "args") "writebackgroundpredictions=true" as true not false. 

#FUNCTION CODE
TSS_calculations <- function(background_preds, presence_preds, threshold) {
  # Count predictions above/below threshold for background
  maior_bb <- sum(background_preds > threshold)
  menor_bb <- sum(background_preds <= threshold)
  
  # Count predictions above/below threshold for presence
  maior_test <- sum(presence_preds > threshold)
  menor_test <- sum(presence_preds <= threshold)
  
  # Calculate sensitivity and specificity
  sensitivity <- maior_test / (maior_test + menor_test)
  specificity <- menor_bb / (menor_bb + maior_bb)
  
  # Calculate TSS
  tss <- sensitivity + specificity - 1
  
  # Return as a list
  return(list(
    TSS = tss,
    Sensitivity = sensitivity,
    Specificity = specificity
  ))
}

#### READ IN AND PROCESS INPUT FILES ####

windows()
#wrld_simpl <- read_sf("D:/DATA/GIS_Layers/World_Shpfiles/WORLD_Dissolved.shp")

#species points
batspp <- read_xlsx("Rhinolophids_SE_Asia_compiled_correct_names v4 02112023.xlsx")
head(batspp)

batspp <- filter(batspp, Latitude <= 36 & Latitude >= -16 & Longitude >= 66 & Longitude <= 156)
summary(batspp)

spnames <- data.frame(unique(batspp$Species))
spcounts <- batspp %>% count(Species)

## Spatial thinning

gridsize <- 0.008333  # value for half a degree = 0.5; value for 10 minutes: 0.167; 5 mins is 0.083; 30 secs is 0.0083

z <- batspp$Longitude
batspp$longmid <- ifelse(z >= 0.0,
                         (z%/%1)+(((z-z%/%1)/gridsize)%/%1)*gridsize,
                         (((z*(-1.0))%/%1)+((((z*(-1.0))-(z*(-1.0))%/%1)/gridsize)%/%1)*gridsize)*(-1.0))


z <- batspp$Latitude
batspp$latmid <- ifelse(z >= 0.0, (z%/%1)+(((z-z%/%1)/gridsize)%/%1)*gridsize,
                        (((z*(-1.0))%/%1)+((((z*(-1.0))-(z*(-1.0))%/%1)/gridsize)%/%1)*gridsize)*(-1.0))

#now remove repeats
batspp2 <- distinct(batspp, Species, longmid, latmid, .keep_all = TRUE)
dim(batspp); dim(batspp2)

#perform clean coordinates
batspp2 <- clean_coordinates(batspp2, lon = "longmid", lat = "latmid",
                             species = "Species", tests = c("capitals","centroids", "equal", "gbif", 
                                                            "institutions", "outliers", "seas","zeros"), 
                             capitals_rad = 10000, 
                             centroids_rad = 1000, centroids_detail = "both", inst_rad = 100,
                             outliers_method = "quantile", outliers_mtp = 5, 
                             outliers_size = 7, range_rad = 0, zeros_rad = 0.5,
                             seas_scale = 50, value = "spatialvalid",
                             verbose = TRUE, report = FALSE)
summary(batspp2)
batspp3 <- filter(batspp2, .summary=="TRUE") #remove false coordinates based on CCtest
summary(batspp3)

#select relevant columns
batspp3 <- batspp3[,1:5]
write.csv(batspp3, file = "input/Rhinolophids_SE_Asia_compiled_correct_names_clean.csv", row.names = F)

# get final counts
spcounts2 <- batspp3 %>% count(Species)
spcounts2 <- left_join(spcounts, spcounts2, by = "Species")

# Get species with enough coordinates for modeling
spcounts2 <- batspp3 %>% count(Species)
spnames2 <- filter(spcounts2, n >= 15) #44 species


# Convert to sf object with lon-lat coordinates
species_sf <- st_as_sf(batspp3, coords = c("longmid", "latmid"), crs = 4326)

# Split species by name
species_split <- species_sf %>%
  group_by(Species) %>%
  summarise()

#save to disk as separate shapefiles for plotting
for (sp in unique(species_sf$Species)) {
  st_write(species_split[species_split$Species == sp, ], paste0("Rhino_shapefiles2/", sp, ".shp"))
}


### Environmental data layers ###
#bioclim, aridity, soils
#read back in
predictors <- brick("D:/PROJECTS/TEAMAP/DATA/TEA_climsoil_predictors_lean_30s_20102023.tif")
names(predictors) <- c("bio01"   ,      "bio03"    ,     "bio04"    ,     "bio07"    ,     "bio13"   ,      "bio15"      ,   "bio19"  ,       "aridity", "bedrock_depth" ,"cec_30cm"   ,   "clay_30cm"   ,  "soilpH_30cm" ,  "sand_30cm")

predictors <- subset(predictors, 1:9)

#ndvi
ndvi <- raster("input/NDVI_2020.tif")
ndvi <- crop(ndvi, predictors)
ndvi <- resample(ndvi, predictors$bio01, method = "bilinear")
ndvi <- ndvi*0.0001
plot(ndvi)
writeRaster(ndvi, filename = "input/NDVI_2020_TEA_resampled.tif", datatype = "FLT4S", overwrite = T)
ndvi <- raster("input/NDVI_2020_TEA_resampled.tif")

#distance
wbdist <- raster("D:/DATA/GIS_Layers/Water Bodies/TEA_WatBod_dist.tif")

wbdist <- crop(wbdist, predictors$bio01)
wbdist <- resample(wbdist, predictors$bio01, method = "bilinear")
wbdist <- mask(wbdist, predictors$bio01)
plot(wbdist)
writeRaster(wbdist, filename = "input/TEA_WatBod_dist_resampled.tif", datatype = "FLT4S")
wbdist <- raster("input/TEA_WatBod_dist_resampled.tif")

# tree height
treeht <- raster("D:/PROJECTS/OneHealth/Asia_eride/TEA_canopy_ht.tif")

treeht <- crop(treeht, predictors$bio01)
treeht <- resample(treeht, predictors$bio01, method = "bilinear")
treeht <- mask(treeht, predictors$bio01)
plot(treeht)
writeRaster(treeht, filename = "input/TEA_Tree_Height_resampled.tif", datatype = "FLT4S")
treeht <- raster("input/TEA_Tree_Height_resampled.tif")

treeht2 <- mask(treeht, predictors_model$bio01)
predictors_model <- stack(predictors_model, treeht2)

#stack them

predictors <- stack(predictors, ndvi, wbdist, treeht)
names(predictors)
names(predictors)[10:12] <- c("NDVI", "dist_wb", "tree_ht")
writeRaster(predictors, filename = "input/TEA_predictors_all_25102023.tif")

writeRaster(predictors_model, filename = "input/TEA_predictors_model_lean_30s_masked_08112023.tif")

# create predictors_model layer for maxent modeling 

# Convert the occurrence data to a SpatialPoints object
coordinates(batspp3) <- c("longmid", "latmid")
proj4string(batspp3) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# Create a raster object with a specified extent and resolution
raster_template <- raster(extent(65.99986, 156.0082, -16.00014, 36.00819), resolution = 0.008333333, crs = "+proj=longlat +datum=WGS84 +no_defs")

# Convert the occurrence data to a raster using the template
species_raster <- rasterize(batspp3, raster_template, field = 1, fun = "count")

# Plot to check
plot(species_raster)

# Save the species raster as a TIFF file
writeRaster(species_raster, "input/Rhinolophids_TEA_allspecies_raster.tif", format = "GTiff", datatype = "INT2S")

## sample 500k points from the predictors extent across TEA

tea500 <- sampleRandom(predictors$bio01, size=500000, na.rm=T, asRaster=T)

#save this layer for the model training
writeRaster(tea500, filename = "input/TEA_500K_samplepts2.tif", format = "GTiff", datatype = "INT2S")

## Combine species_raster with the sampling raster for modeling
tea500[!is.na(tea500)] <- 1
tea500[is.na(tea500)] <- 0
plot(tea500)

species_raster[!is.na(species_raster)] <- 1
species_raster[is.na(species_raster)] <- 0
plot(species_raster)

predictors_mask <- species_raster + tea500
predictors_mask[predictors_mask==2] <- 1
plot(predictors_mask)
predictors_mask <- mask(predictors_mask, predictors$bio01)
writeRaster(predictors_mask, filename = "input/TEA500k_allspecies_sample_mask_08112023.tif", overwrite=T, datatype="INT2S")
predictors_mask[predictors_mask==0] <- NA

# Prepare layers for modeling
#mask out predictors for modeling only
predictors_model <- mask(predictors, predictors_mask)
predictors_model
plot(predictors_model$bio01)
writeRaster(predictors_model, filename = "input/TEA_predictors_model_lean_30s_masked_25102023.tif")


#### READ back all inputs (to clear memory before maxent modeling!) ####

# perform correlation on soils using the same sampling points 
x3 <- as.data.frame(predictors_model, na.rm = T)
x3.cor <- cor(x3, method="spearman", use = "pairwise.complete.obs")

write.csv(x3.cor, file = "input/x3_cor_masked_w_canopyht.csv", row.names = T)

# drop layers with correlation r>0.7

predictors <- dropLayer(predictors, c("bio03", "bio04", "bio19", "NDVI"))
#predictors_model <- dropLayer(predictors_model, c("bio03", "bio04", "bio19", "NDVI"))

predictors[[1]] <- predictors[[1]]*0.1 - 273.15
#predictors_model[[1]] <- predictors_model[[1]]*0.1 - 273.15

# align NAs across all layers
predictors <- mask(predictors, calc(predictors,fun = sum))
par(mfrow=c(1,1))
plot(predictors$dist_wb)

writeRaster(predictors, filename = "input/TEA_predictors_lean_30s_masked_25102023.tif")


# FINAL read back in before modeling 
predictors <- brick("input/TEA_predictors_lean_30s_masked_08112023.tif")
predictors_model <- brick("input/TEA_predictors_model_lean_30s_masked_08112023.tif")
batspp3 <- read.csv("input/Rhinolophids_SE_Asia_compiled_correct_names_clean.csv") 

names(predictors) <- c("bio01", "bio07", "bio13", "bio15", "aridity", "bedrock_depth",  "dist_wb", "tree_ht")
names(predictors_model) <- names(predictors)

### MAXENT MODELING #### 
set.seed(123)
sp.name <- "Rhinolophus rouxii"
sp.occ <- filter(batspp3, Species == sp.name) 

sp.occ <- as.data.frame(sp.occ)
sp.occ <- dplyr::select(sp.occ, longmid, latmid)
head(sp.occ)


# Read the polygon convex hull shapefile
sp.ch <- st_read(paste0("MCP/CH_",sp.name, ".shp"))

# Project the point data to Eckert IV projection
eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
sp.ch <- st_transform(sp.ch, crs = eckertIV)

# Buffer the polygons by 500 km
sp.ch.buf <- st_buffer(sp.ch, dist = 500000)

# Reproject the buffered polygons back to WGS84
sp.ch.buf <- st_transform(sp.ch.buf, crs = raster::crs(predictors))

# Plot the raster and points
plot(predictors$bio01, main = names(predictors)[1])
points(sp.occ)

# Plot the buffered polygons
plot(sp.ch.buf, border = "blue", col = NA, add = TRUE)


# Crop environmental rasters to match the study extent
predictors.bg <- raster::crop(predictors, sp.ch.buf)
# Next, mask the rasters to the shape of the buffers
predictors.bg <- raster::mask(predictors.bg, sp.ch.buf)
plot(predictors.bg[[1]], main = names(predictors)[1])
points(sp.occ)
plot(sp.ch.buf, border = "blue", col = NA, add = T)


# Randomly sample 10,000 background points from one background extent raster 
# (only one per cell without replacement).
bg_points <- sampleRandom(predictors.bg[[1]], 10000, sp = TRUE)

# Convert the sampled points to a data frame and set column names
bg <- as.data.frame(bg_points[,-1])
colnames(bg) <- colnames(sp.occ)

# Notice how we have pretty good coverage (every cell).
plot(predictors.bg[[1]])
points(bg, pch = 20, cex = 0.2)

# ENMEval ####
set.seed(123)
rm <- seq(from =1, to = 5, by = 1)
tune.args <- list(fc=c("L", "LQ", "LQH", "H"), rm=rm)

# 5. No raster data (a.k.a, samples with data, or SWD): no full model raster predictions 
# created, so will run faster; also, both cbi.train and cbi.val will be calculated on the 
# point data (training and validation background) instead of on the "envs" rasters (default).
# For this implementation, assigning the categorical variable to factor with the argument 
# "categoricals" is easier, as ENMevaluate() internally assigns the levels based on both 
# occs and bg, avoiding any errors associated with different factor levels when combining data.
occs.z <- cbind(sp.occ, raster::extract(predictors, sp.occ))
bg.z <- cbind(bg, raster::extract(predictors, bg))

e.swd <- ENMevaluate(occs.z, bg = bg.z, algorithm = "maxent.jar", tune.args = tune.args, parallel = T, numCores = 3, taxon.name = sp.name,
                     partitions = "block")
write.csv(e.swd@results, paste0("outputs/enmeval/", sp.name,"_enmeval_results_30s_swd1.csv"), row.names = F)

# Overall results
res <- eval.results(e.swd)
# Select the model with delta AICc equal to 0, or the one with the lowest AICc score.
# In practice, models with delta AICc scores less than 2 are usually considered 
# statistically equivalent.
opt.aicc <- res %>% filter(delta.AICc == 0)
opt.aicc

# This dplyr operation executes the sequential criteria explained above.
opt.seq <- res %>% 
  filter(or.10p.avg == min(or.10p.avg)) %>% 
  filter(auc.val.avg == max(auc.val.avg))
opt.seq

#Let’s now choose the optimal model settings based on the sequential criteria and examine it.
# We can select a single model from the ENMevaluation object using the tune.args of our
# optimal model.

# lowest delta AIC
mod.seq <- eval.models(e.swd)[[opt.aicc$tune.args]]
# Here are the non-zero coefficients in our model.
mod.seq@lambdas

# sequential criteria
mod.seq2 <- eval.models(e.swd)[[opt.seq$tune.args]]
# Here are the non-zero coefficients in our model.
mod.seq2@lambdas

# And these are the marginal response curves for the predictor variables wit non-zero 
# coefficients in our model. We define the y-axis to be the cloglog transformation, which
# is an approximation of occurrence probability (with assumptions) bounded by 0 and 1
# (Phillips et al. 2017).
plot(mod.seq, type = "cloglog")
#plot(mod.seq2, type = "cloglog")

# maxent.jar models use the dismo::response() function for this
dismo::response(eval.models(e.swd)[[opt.aicc$tune.args]])


# predict the model to space
dir.create(paste0("outputs/", sp.name))

#optimal aicc
pred.seq <- predict(object=mod.seq,
                    x=predictors,
                    filename = paste0("outputs/", sp.name,"/",sp.name,"_R_enm_",opt.aicc$tune.args),
                    na.rm=TRUE,
                    format='GTiff',
                    overwrite=TRUE,
                    doclamp= TRUE)

#optimal sequential
pred.seq2 <- predict(object=mod.seq2,
                    x=predictors,
                    filename = paste0("outputs/", sp.name,"/",sp.name,"_R_enm_",opt.seq$tune.args),
                    na.rm=TRUE,
                    format='GTiff',
                    overwrite=TRUE,
                    doclamp= TRUE)

par(mfrow = c(1,2))
plot(pred.seq); plot(pred.seq)
points(sp.occ, pch=3, cex = 0.5)


# TSS 
# extract the cloglog/logistic results of the background pts
backgroundclog <- predict(mod.seq, mod.seq@absence, response)

# extract the cloglog/logistic results of the species pts
sampleclog <- predict(mod.seq, mod.seq@presence, response)

#set n the number of pseudoabsences used for background predictions by MaxEnt
n <- nrow(e.swd@bg)

### Set the threshold rule here
maxres <- as.data.frame(t(mod.seq@results))
th <- maxres$X10.percentile.training.presence.Cloglog.threshold

#run the function, the input values are the sampleclog values, then the background clog values, the sample number for the pseudo absences and then threshold value
ts <- TSS_calculations(backgroundclog, sampleclog, th)
ts
maxres$Training.AUC


#omission rates
OR <- opt.aicc[,12:15]

# accuracy results
df <- data.frame(sp.name, opt.aicc[,1:2], maxres$Training.AUC, ts, OR) #optimal aicc


# save to disk
file_path <- "outputs/Maxent_accuracy_test_results.csv"
# If the file doesn't exist, create a new file and save the data frame
#write.csv(df, file_path, row.names = FALSE)

# If the file already exists, append the new rows to the existing file
existing_df <- read.csv(file_path)
updated_df <- rbind(existing_df, df)
write.csv(updated_df, file_path, row.names = FALSE)


######## END OF CODE ##########


