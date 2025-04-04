library(raster)
library(tidyverse)
library(tidygraph)
library(igraph)
library(doParallel)
library(progress)
library(rgdal)
library(parallel)
library(parabar)

setwd("D:/PROJECTS/OneHealth/Maxent")

#create network from the raster file

#load population density raster
r <- raster("eride/WORLDPOP_ASIA_PD_2020_TEA_FINAL_UTM.tif")

#load population at risk (PAR) raster
risk_raster <- raster("eride/TEA_1km_frags_poprisk1.tif")

shpfl <- readOGR("Ecoregions/Asia_mainland_diss_UTM.shp")
shpfl <- readOGR("Ecoregions/Asia_islands_UTM.shp")

windows()
plot(r)
plot(risk_raster)
plot(shpfl, add = T)

r <- mask(r, shpfl)
risk_raster <- mask(risk_raster, shpfl)



# convert risk_raster to data frame
risk_raster[r==0] <- NaN  #make sure the 0 in the population density raster is also NA
risk_raster[is.na(r)] <- NaN
risk_df <- as.data.frame(risk_raster, xy = TRUE, na.rm = T)
colnames(risk_df)[3] <- "risk" 

#write.csv(risk_df, file = "eride/TEA_ISLANDS_1km_frags_poprisk1_10km_riskdf.csv", row.names = F)

#identify nodes
nodes = tibble(value = getValues(r)) %>%
  bind_cols(as_tibble(coordinates(r))) 
nodes[nodes==0] <- NaN #convert 0's to NA

#write.csv(nodes, file = "eride/WORLDPOP_ASIA_ISLANDS_PD_2020_TEA_FINAL_UTM_10km_nodes.csv", row.names = F)

#identify connected pixels, and assign edge weights based on the population density
edges = as_tibble(adjacent(r, cell = 1:ncell(r), directions = 4, sorted = TRUE)) %>%
  mutate(weight = 1/(nodes$value[from] * nodes$value[to]))

edges[edges < 0] <- 0

#write.csv(edges, file = "eride/WORLDPOP_ASIA_ISLANDS_PD_2020_TEA_FINAL_UTM_10km_edges.csv", row.names = F)


#create graph
g = tbl_graph(nodes, edges, directed = FALSE) %>%
  filter(!is.na(value)) %>%
  delete_vertices(V(.)[degree(.) == 0]) #delete vertices with degree 0


#plot to check
windows()
ly <- layout_with_fr(g)
plot(g, layout = ly)

# Update edge weights to reflect inverse
E(g)$weight <- 1 / E(g)$weight

# keep 2 decimal places for edge weights to reduce file sizes
t <- E(g)$weight
t2 <- format(round(t, 3), nsmall = 3)
E(g)$weight <- t2

##### PANDEMIC SPREAD POTENTIAL ####
### Note: This code was run using an HPC and the vertices were split into multiple runs for efficiency.


# loop over all vertices in the graph in parallel
cl <- makeCluster(detectCores())
clusterExport(cl, c("g", "risk_df"))

# Calculate the potential spread of pandemic for each vertex in parallel
pandemic_spread_list <- parLapply(cl, 1:vcount(g), function(start_vertex) {
  # Define a function to calculate the potential spread of pandemic for a start vertex
  calculate_pandemic_spread <- function(start_vertex, g, risk_df) {
    # Obtain the shortest path distances and ending vertices
    shortest_paths <- igraph::shortest_paths(g, from = start_vertex, to = igraph::V(g), output = "both", weights = igraph::E(g)$weight)$epath
    
    # Initialize empty vector to store shortest path edge weights
    shortest_path_weights <- c()
    
    # Loop through each shortest path and extract the sum of edge weights
    for (j in 1:length(shortest_paths)) {
      path <- shortest_paths[[j]]
      path_weights <- as.numeric(igraph::E(g)[path]$weight)
      shortest_path_weights[j] <- sum(path_weights)
    }
    
    
    #Calculate the potential spread of pandemic from each pixel
    # Extract values from the risk data frame
    values <- data.frame(risk_df$risk, shortest_path_weights)
    colnames(values) <- c("risk", "dist")
    values$ps <- values$risk * values$dist 
    
    # Sum the values
    pandemic_spread <- sum(values$ps, na.rm = TRUE)
    
    # Return the potential spread of pandemic for each pixel
    pandemic_spread
    
    
  }
  
  calculate_pandemic_spread(start_vertex, g, risk_df)
})


stopCluster(cl)

# Summarize pandemic_spread_list by summing values per [[x]], ignoring NAs
pandemic_spread_df <- unlist(pandemic_spread_list)

# Separate the two lists into separate vectors
potential_spread_vec <- unlist(lapply(pandemic_spread_list, "[[", 1))

# bind the xy coordinates
pandemic_spread_df <- data.frame(risk_df[151:300,1:2], potential_spread_vec)

# save to disk 
write.csv(pandemic_spread_df, file = "eride/Asia_mainland_test/r003.csv", row.names = F)

# Create a raster of the values
pandemic_spr_ras <- raster::rasterize(pandemic_spread_df[,1:2], r, field = pandemic_spread_df[,3], fun = mean)

# plot
par(mfrow = c(1,2))

#plot(log(r), main = "pop density (log)")
plot(risk_raster, main = "pop at risk")
plot(pandemic_spr_ras, main = "potential spread")

# save to disk
writeRaster(pandemic_spr_ras, filename = "eride/psy.tif", overwrite = T)
writeRaster(log(pandemic_spr_ras), filename = "eride/psylog.tif" , overwrite = T)


##### END OF CODE ######

