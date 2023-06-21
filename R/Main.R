library(terra)
library(sf)
library(leastcostpath)
library(foreach)
library(tmap)
library(dplyr)
library(ggplot2)
library(scales)

#### SET-UP ####

set.seed(1, kind = "L'Ecuyer-CMRG")
nsims <- 1000
sard_nsims <- 50
neigh <- 4
ncores <- 50

source("./R/Functions.R")

#### CASE STUDY ONE: TACTICAL SIMULATION ####

#### PREPROCESSS DIGITAL ELEVATION MODEL ####

dem_220_scenario_1 <- terra::rast("./Data/Tactical_sim/FS_220.tif")
dem_220_scenario_1 <- normalise_raster(dem_220_scenario_1) * 5
names(dem_220_scenario_1) <- "Fractal Dimension: 2.2"

dem_220_scenario_2 <- terra::rast("./Data/Tactical_sim/FS_220.tif")
dem_220_scenario_2 <- normalise_raster(dem_220_scenario_2) * 10
names(dem_220_scenario_2) <- "Fractal Dimension: 2.2"

dem_230_scenario_1 <- terra::rast("./Data/Tactical_sim/FS_230.tif")
dem_230_scenario_1 <- normalise_raster(dem_230_scenario_1) * 5
names(dem_230_scenario_1) <- "Fractal Dimension: 2.3"

dem_230_scenario_2 <- terra::rast("./Data/Tactical_sim/FS_230.tif")
dem_230_scenario_2 <- normalise_raster(dem_230_scenario_2) * 10
names(dem_230_scenario_2) <- "Fractal Dimension: 2.3"

dem_240_scenario_1 <- terra::rast("./Data/Tactical_sim/FS_240.tif")
dem_240_scenario_1 <- normalise_raster(dem_240_scenario_1) * 5
names(dem_240_scenario_1) <- "Fractal Dimension: 2.4"

dem_240_scenario_2 <- terra::rast("./Data/Tactical_sim/FS_240.tif")
dem_240_scenario_2 <- normalise_raster(dem_240_scenario_2) * 10
names(dem_240_scenario_2) <- "Fractal Dimension: 2.4"

dem_250_scenario_1 <- terra::rast("./Data/Tactical_sim/FS_250.tif")
dem_250_scenario_1 <- normalise_raster(dem_250_scenario_1) * 5
names(dem_250_scenario_1) <- "Fractal Dimension: 2.5"

dem_250_scenario_2 <- terra::rast("./Data/Tactical_sim/FS_250.tif")
dem_250_scenario_2 <- normalise_raster(dem_250_scenario_2) * 10
names(dem_250_scenario_2) <- "Fractal Dimension: 2.5"

dem_260_scenario_1 <- terra::rast("./Data/Tactical_sim/FS_260.tif")
dem_260_scenario_1 <- normalise_raster(dem_260_scenario_1) * 5
names(dem_260_scenario_1) <- "Fractal Dimension: 2.6"

dem_260_scenario_2 <- terra::rast("./Data/Tactical_sim/FS_260.tif")
dem_260_scenario_2 <- normalise_raster(dem_260_scenario_2) * 10
names(dem_260_scenario_2) <- "Fractal Dimension: 2.6"

#### CREATE COST FUNCTION TABLE ####

cfs <- c("tobler", "modified tobler", "irmischer-clarke male", "irmischer-clarke offpath male", "irmischer-clarke female", "irmischer-clarke offpath female", 'rees', "davey", 'garmy', 'kondo-saino', 'naismith', "campbell 2019 50", "herzog", "llobera-sluckin")

index_table <- expand.grid(1:length(cfs), 1:length(cfs))
index_table <- index_table[!duplicated(apply(index_table,1,function(x) paste(sort(x),collapse=''))),]
index_table <- index_table[index_table$Var1 != index_table$Var2,]
rownames(index_table) <- 1:nrow(index_table)

cs_list_220_scenario_1 <- apply(X = data.frame(cfs), MARGIN = 1, FUN = function(x) {  
  cs <- leastcostpath::create_slope_cs(x = dem_220_scenario_1, cost_function = x, neighbours = neigh, crit_slope = NULL)
  return(cs)
})

cs_list_220_scenario_2 <- apply(X = data.frame(cfs), MARGIN = 1, FUN = function(x) {  
  cs <- leastcostpath::create_slope_cs(x = dem_220_scenario_2, cost_function = x, neighbours = neigh, crit_slope = NULL)
  return(cs)
})

cs_list_230_scenario_1 <- apply(X = data.frame(cfs), MARGIN = 1, FUN = function(x) {  
  cs <- leastcostpath::create_slope_cs(x = dem_230_scenario_1, cost_function = x, neighbours = neigh, crit_slope = NULL)
  return(cs)
})

cs_list_230_scenario_2 <- apply(X = data.frame(cfs), MARGIN = 1, FUN = function(x) {  
  cs <- leastcostpath::create_slope_cs(x = dem_230_scenario_2, cost_function = x, neighbours = neigh, crit_slope = NULL)
  return(cs)
})

cs_list_240_scenario_1 <- apply(X = data.frame(cfs), MARGIN = 1, FUN = function(x) {  
  cs <- leastcostpath::create_slope_cs(x = dem_240_scenario_1, cost_function = x, neighbours = neigh, crit_slope = NULL)
  return(cs)
})

cs_list_240_scenario_2 <- apply(X = data.frame(cfs), MARGIN = 1, FUN = function(x) {  
  cs <- leastcostpath::create_slope_cs(x = dem_240_scenario_2, cost_function = x, neighbours = neigh, crit_slope = NULL)
  return(cs)
})

cs_list_250_scenario_1 <- apply(X = data.frame(cfs), MARGIN = 1, FUN = function(x) {  
  cs <- leastcostpath::create_slope_cs(x = dem_250_scenario_1, cost_function = x, neighbours = neigh, crit_slope = NULL)
  return(cs)
})

cs_list_250_scenario_2 <- apply(X = data.frame(cfs), MARGIN = 1, FUN = function(x) {  
  cs <- leastcostpath::create_slope_cs(x = dem_250_scenario_2, cost_function = x, neighbours = neigh, crit_slope = NULL)
  return(cs)
})

cs_list_260_scenario_1 <- apply(X = data.frame(cfs), MARGIN = 1, FUN = function(x) {  
  cs <- leastcostpath::create_slope_cs(x = dem_260_scenario_1, cost_function = x, neighbours = neigh, crit_slope = NULL)
  return(cs)
})

cs_list_260_scenario_2 <- apply(X = data.frame(cfs), MARGIN = 1, FUN = function(x) {  
  cs <- leastcostpath::create_slope_cs(x = dem_260_scenario_2, cost_function = x, neighbours = neigh, crit_slope = NULL)
  return(cs)
})

#### CALCULATE AND COMPARE LCPs CREATED USING DIFFERENT COST FUNCTIONS ####

# Used to calculate the two random points within the bounding box of the DEM. Terra::wrap required to make DEM available when parallelisng the calculation
dem_packaged <- terra::wrap(dem_220_scenario_1)

#### LCP 220 Fractal Dimension, Scenario 1 ####

myCluster <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(myCluster)

lcps_220_scenario_1 <- foreach(i=1:nsims) %dopar% {
  
  pts <- sf::st_sample(sf::st_as_sfc(sf::st_bbox(terra::rast(dem_packaged))), 2, type = "random")
  pts <- sf::st_as_sf(pts)
  
  compare_lcp_220 <- apply(X = index_table, MARGIN = 1, FUN = function(x) { compare_lcps(cs_list = cs_list_220_scenario_1, cfs = cfs, cf1_index = x[1], cf2_index = x[2], pts = pts)})
  
  compare_lcp_220 <- do.call(rbind, compare_lcp_220)
  compare_lcp_220$sim_no <- i
  compare_lcp_220$fractal <- "2.20"
  
  return(compare_lcp_220)
  
}

lcps_220_scenario_1_2 <- do.call(rbind, lcps_220_scenario_1)
lcps_220_scenario_1_3 <- sf::st_drop_geometry(lcps_220_scenario_1_2)

parallel::stopCluster(myCluster)

write.csv(lcps_220_scenario_1_3, "./Output/Data/simulated_lcps_220_scenario_1.csv")

#### LCP 220 Fractal Dimension, Scenario 2 ####

myCluster <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(myCluster)

lcps_220_scenario_2 <- foreach(i=1:nsims) %dopar% {
  
  pts <- sf::st_sample(sf::st_as_sfc(sf::st_bbox(terra::rast(dem_packaged))), 2, type = "random")
  pts <- sf::st_as_sf(pts)
  
  compare_lcp_220 <- apply(X = index_table, MARGIN = 1, FUN = function(x) { compare_lcps(cs_list = cs_list_220_scenario_2, cfs = cfs, cf1_index = x[1], cf2_index = x[2], pts = pts)})
  
  compare_lcp_220 <- do.call(rbind, compare_lcp_220)
  compare_lcp_220$sim_no <- i
  compare_lcp_220$fractal <- "2.20"
  
  return(compare_lcp_220)
  
}

lcps_220_scenario_2_2 <- do.call(rbind, lcps_220_scenario_2)
lcps_220_scenario_2_3 <- sf::st_drop_geometry(lcps_220_scenario_2_2)

parallel::stopCluster(myCluster)

write.csv(lcps_220_scenario_2_3, "./Output/Data/simulated_lcps_220_scenario_2.csv")

#### LCP 230 Fractal Dimension, Scenario 1 ####

myCluster <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(myCluster)

lcps_230_scenario_1 <- foreach(i=1:nsims) %dopar% {
  
  pts <- sf::st_sample(sf::st_as_sfc(sf::st_bbox(terra::rast(dem_packaged))), 2, type = "random")
  pts <- sf::st_as_sf(pts)
  
  compare_lcp_230 <- apply(X = index_table, MARGIN = 1, FUN = function(x) { compare_lcps(cs_list = cs_list_230_scenario_1, cfs = cfs, cf1_index = x[1], cf2_index = x[2], pts = pts)})
  
  compare_lcp_230 <- do.call(rbind, compare_lcp_230)
  compare_lcp_230$sim_no <- i
  compare_lcp_230$fractal <- "2.30"
  
  return(compare_lcp_230)
  
}

lcps_230_scenario_1_2 <- do.call(rbind, lcps_230_scenario_1)
lcps_230_scenario_1_3 <- sf::st_drop_geometry(lcps_230_scenario_1_2)

parallel::stopCluster(myCluster)

write.csv(lcps_230_scenario_1_3, "./Output/Data/simulated_lcps_230_scenario_1.csv")

#### LCP 230 Fractal Dimension, Scenario 2 ####

myCluster <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(myCluster)

lcps_230_scenario_2 <- foreach(i=1:nsims) %dopar% {
  
  pts <- sf::st_sample(sf::st_as_sfc(sf::st_bbox(terra::rast(dem_packaged))), 2, type = "random")
  pts <- sf::st_as_sf(pts)
  
  compare_lcp_230 <- apply(X = index_table, MARGIN = 1, FUN = function(x) { compare_lcps(cs_list = cs_list_230_scenario_2, cfs = cfs, cf1_index = x[1], cf2_index = x[2], pts = pts)})
  
  compare_lcp_230 <- do.call(rbind, compare_lcp_230)
  compare_lcp_230$sim_no <- i
  compare_lcp_230$fractal <- "2.30"
  
  return(compare_lcp_230)
  
}

lcps_230_scenario_2_2 <- do.call(rbind, lcps_230_scenario_2)
lcps_230_scenario_2_3 <- sf::st_drop_geometry(lcps_230_scenario_2_2)

parallel::stopCluster(myCluster)

write.csv(lcps_230_scenario_2_3, "./Output/Data/simulated_lcps_230_scenario_2.csv")

#### LCP 240 Fractal Dimension, Scenario 1 ####

myCluster <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(myCluster)

lcps_240_scenario_1 <- foreach(i=1:nsims) %dopar% {
  
  pts <- sf::st_sample(sf::st_as_sfc(sf::st_bbox(terra::rast(dem_packaged))), 2, type = "random")
  pts <- sf::st_as_sf(pts)
  
  compare_lcp_240 <- apply(X = index_table, MARGIN = 1, FUN = function(x) { compare_lcps(cs_list = cs_list_240_scenario_1, cfs = cfs, cf1_index = x[1], cf2_index = x[2], pts = pts)})
  
  compare_lcp_240 <- do.call(rbind, compare_lcp_240)
  compare_lcp_240$sim_no <- i
  compare_lcp_240$fractal <- "2.40"
  
  return(compare_lcp_240)
  
}

lcps_240_scenario_1_2 <- do.call(rbind, lcps_240_scenario_1)
lcps_240_scenario_1_3 <- sf::st_drop_geometry(lcps_240_scenario_1_2)

parallel::stopCluster(myCluster)

write.csv(lcps_240_scenario_1_3, "./Output/Data/simulated_lcps_240_scenario_1.csv")

#### LCP 240 Fractal Dimension, Scenario 2 ####

myCluster <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(myCluster)

lcps_240_scenario_2 <- foreach(i=1:nsims) %dopar% {
  
  pts <- sf::st_sample(sf::st_as_sfc(sf::st_bbox(terra::rast(dem_packaged))), 2, type = "random")
  pts <- sf::st_as_sf(pts)
  
  compare_lcp_240 <- apply(X = index_table, MARGIN = 1, FUN = function(x) { compare_lcps(cs_list = cs_list_240_scenario_2, cfs = cfs, cf1_index = x[1], cf2_index = x[2], pts = pts)})
  
  compare_lcp_240 <- do.call(rbind, compare_lcp_240)
  compare_lcp_240$sim_no <- i
  compare_lcp_240$fractal <- "2.40"
  
  return(compare_lcp_240)
  
}

lcps_240_scenario_2_2 <- do.call(rbind, lcps_240_scenario_2)
lcps_240_scenario_2_3 <- sf::st_drop_geometry(lcps_240_scenario_2_2)

parallel::stopCluster(myCluster)

write.csv(lcps_240_scenario_2_3, "./Output/Data/simulated_lcps_240_scenario_2.csv")

#### LCP 250 Fractal Dimension, Scenario 1 ####

myCluster <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(myCluster)

lcps_250_scenario_1 <- foreach(i=1:nsims) %dopar% {
  
  pts <- sf::st_sample(sf::st_as_sfc(sf::st_bbox(terra::rast(dem_packaged))), 2, type = "random")
  pts <- sf::st_as_sf(pts)
  
  compare_lcp_250 <- apply(X = index_table, MARGIN = 1, FUN = function(x) { compare_lcps(cs_list = cs_list_250_scenario_1, cfs = cfs, cf1_index = x[1], cf2_index = x[2], pts = pts)})
  
  compare_lcp_250 <- do.call(rbind, compare_lcp_250)
  compare_lcp_250$sim_no <- i
  compare_lcp_250$fractal <- "2.50"
  
  return(compare_lcp_250)
  
}

lcps_250_scenario_1_2 <- do.call(rbind, lcps_250_scenario_1)
lcps_250_scenario_1_3 <- sf::st_drop_geometry(lcps_250_scenario_1_2)

parallel::stopCluster(myCluster)

write.csv(lcps_250_scenario_1_3, "./Output/Data/simulated_lcps_250_scenario_1.csv")

#### LCP 250 Fractal Dimension, Scenario 2 ####

myCluster <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(myCluster)

lcps_250_scenario_2 <- foreach(i=1:nsims) %dopar% {
  
  pts <- sf::st_sample(sf::st_as_sfc(sf::st_bbox(terra::rast(dem_packaged))), 2, type = "random")
  pts <- sf::st_as_sf(pts)
  
  compare_lcp_250 <- apply(X = index_table, MARGIN = 1, FUN = function(x) { compare_lcps(cs_list = cs_list_250_scenario_2, cfs = cfs, cf1_index = x[1], cf2_index = x[2], pts = pts)})
  
  compare_lcp_250 <- do.call(rbind, compare_lcp_250)
  compare_lcp_250$sim_no <- i
  compare_lcp_250$fractal <- "2.50"
  
  return(compare_lcp_250)
  
}

lcps_250_scenario_2_2 <- do.call(rbind, lcps_250_scenario_2)
lcps_250_scenario_2_3 <- sf::st_drop_geometry(lcps_250_scenario_2_2)

parallel::stopCluster(myCluster)

write.csv(lcps_250_scenario_2_3, "./Output/Data/simulated_lcps_250_scenario_2.csv")

#### LCP 260 Fractal Dimension, Scenario 1 ####

myCluster <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(myCluster)

lcps_260_scenario_1 <- foreach(i=1:nsims) %dopar% {
  
  pts <- sf::st_sample(sf::st_as_sfc(sf::st_bbox(terra::rast(dem_packaged))), 2, type = "random")
  pts <- sf::st_as_sf(pts)
  
  compare_lcp_260 <- apply(X = index_table, MARGIN = 1, FUN = function(x) { compare_lcps(cs_list = cs_list_260_scenario_1, cfs = cfs, cf1_index = x[1], cf2_index = x[2], pts = pts)})
  
  compare_lcp_260 <- do.call(rbind, compare_lcp_260)
  compare_lcp_260$sim_no <- i
  compare_lcp_260$fractal <- "2.60"
  
  return(compare_lcp_260)
  
}

lcps_260_scenario_1_2 <- do.call(rbind, lcps_260_scenario_1)
lcps_260_scenario_1_3 <- sf::st_drop_geometry(lcps_260_scenario_1_2)

parallel::stopCluster(myCluster)

write.csv(lcps_260_scenario_1_3, "./Output/Data/simulated_lcps_260_scenario_1.csv")

#### LCP 260 Fractal Dimension, Scenario 2 ####

myCluster <- parallel::makeCluster(ncores)
doParallel::registerDoParallel(myCluster)

lcps_260_scenario_2 <- foreach(i=1:nsims) %dopar% {
  
  pts <- sf::st_sample(sf::st_as_sfc(sf::st_bbox(terra::rast(dem_packaged))), 2, type = "random")
  pts <- sf::st_as_sf(pts)
  
  compare_lcp_260 <- apply(X = index_table, MARGIN = 1, FUN = function(x) { compare_lcps(cs_list = cs_list_260_scenario_2, cfs = cfs, cf1_index = x[1], cf2_index = x[2], pts = pts)})
  
  compare_lcp_260 <- do.call(rbind, compare_lcp_260)
  compare_lcp_260$sim_no <- i
  compare_lcp_260$fractal <- "2.60"
  
  return(compare_lcp_260)
  
}

lcps_260_scenario_2_2 <- do.call(rbind, lcps_260_scenario_2)
lcps_260_scenario_2_3 <- sf::st_drop_geometry(lcps_260_scenario_2_2)

parallel::stopCluster(myCluster)

write.csv(lcps_260_scenario_2_3, "./Output/Data/simulated_lcps_260_scenario_2.csv")

#### METHODS PLOT ONE ######

cols <- grey.colors(10)

methods_plot1a <- tmap::tm_shape(c(dem_220_scenario_1, 
                                    dem_230_scenario_1, 
                                    dem_240_scenario_1, 
                                    dem_250_scenario_1, 
                                    dem_260_scenario_1)) + 
  tmap::tm_raster(palette = cols[1:5], title = "Elevation (m)", legend.is.portrait = FALSE, breaks = seq(0, 5, 1)) + 
  tmap::tm_facets(nrow = 1) + 
  tmap::tm_layout(panel.label.size = 0.8,
                  legend.outside.position = "bottom",
                  legend.outside.size = 0.35,
                  legend.outside = TRUE,
                  main.title = "A",
                  main.title.position = "left")

methods_plot1b <- tmap::tm_shape(c(dem_220_scenario_2, 
                                    dem_230_scenario_2, 
                                    dem_240_scenario_2, 
                                    dem_250_scenario_2, 
                                    dem_260_scenario_2)) + 
  tmap::tm_raster(palette = cols[1:10], title = "Elevation (m)", legend.is.portrait = FALSE, breaks = seq(0, 10, 1)) + 
  tmap::tm_facets(nrow = 1) + 
  tmap::tm_layout(panel.label.size = 0.8,
                  legend.outside.position = "bottom",
                  legend.outside.size = 0.35,
                  legend.outside = TRUE,
                  main.title = "B",
                  main.title.position = "left")

methods_plot1c <- tmap::tmap_arrange(methods_plot1a, methods_plot1b, ncol = 1)

tmap::tmap_save(tm = methods_plot1c, filename = "./Output/Figures/methods_plot_01.png", dpi = 300)

#### METHODS PLOT TWO #####

dems_scenario_1 <- list(dem_220_scenario_1, 
                        dem_230_scenario_1,
                        dem_240_scenario_1,
                        dem_250_scenario_1,
                        dem_260_scenario_1)

dems_scenario_1_df <- lapply(dems_scenario_1, function(x) { 
  dem_slope_vals <- calculate_slope(x, 4)
  
  slope2deg <- function(slope) { (atan(slope) * 180/pi)}
  
  slope_vals <- slope2deg(dem_slope_vals)
  
  df <- data.frame(scenario = 1, slope_vals = slope_vals)
  
  return(df)
  
}

)

dems_scenario_1_df2 <- do.call(rbind, dems_scenario_1_df)
dems_scenario_1_df2$fractal_dimension <- rep(c(2.2, 2.3, 2.4, 2.5, 2.6), each = nrow(dems_scenario_1_df2) /5)
dems_scenario_1_df2$fractal_dimension <- factor(dems_scenario_1_df2$fractal_dimension, levels = rev(unique(dems_scenario_1_df2$fractal_dimension)))

dems_scenario_2 <- list(dem_220_scenario_2, 
                        dem_230_scenario_2,
                        dem_240_scenario_2,
                        dem_250_scenario_2,
                        dem_260_scenario_2)

dems_scenario_2_df <- lapply(dems_scenario_2, function(x) { 
  
  dem_slope_vals <- calculate_slope(x, 4)
  
  slope2deg <- function(slope) { (atan(slope) * 180/pi)}
  
  slope_vals <- slope2deg(dem_slope_vals)
  
  df <- data.frame(scenario = 2, slope_vals = slope_vals)
  
  return(df)
  
}

)

dems_scenario_2_df2 <- do.call(rbind, dems_scenario_2_df)
dems_scenario_2_df2$fractal_dimension <- rep(c(2.2, 2.3, 2.4, 2.5, 2.6), each = nrow(dems_scenario_2_df2) /5)
dems_scenario_2_df2$fractal_dimension <- factor(dems_scenario_2_df2$fractal_dimension, levels = rev(unique(dems_scenario_2_df2$fractal_dimension)))

methods_plot2a <- ggplot(dems_scenario_1_df2) +
  ggridges::geom_density_ridges(aes(x = slope_vals, y = fractal_dimension, group = fractal_dimension), scale = 3) + 
  xlim(c(-50, 50)) + 
  labs(y = "Fractal Dimension", x = "Slope Degrees (°)", title = "A") + 
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

methods_plot2b <- ggplot(dems_scenario_2_df2) +
  ggridges::geom_density_ridges(aes(x = slope_vals, y = fractal_dimension, group = fractal_dimension), scale = 3) + 
  xlim(c(-50, 50)) +
  labs(y = "Fractal Dimension", x = "Slope Degrees (°)", title = "B") + 
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

methods_plot2c <- gridExtra::grid.arrange(methods_plot2a, methods_plot2b, nrow = 1)

ggplot2::ggsave(filename = "methods_plot_02.png", plot = methods_plot2c, path = "./Output/Figures/", dpi = 300, width = 7, height = 4)

#### METHODS PLOT THREE #####

cf_df <- list()

for (i in 1:length(cfs)) {
  
  print(i)
  print(cfs[i])
  
  time_cfs <- c("tobler", "modified tobler", "irmischer-clarke male", "irmischer-clarke offpath male", "irmischer-clarke female", "irmischer-clarke offpath female", 'rees', "davey", 'garmy', 'kondo-saino', 'naismith', "campbell 2019 50")
  
  energy_cfs <- c("herzog", "llobera-sluckin")
  
  if(cfs[i] %in% time_cfs) { 
    type <- "Time-based"
  } else if (cfs[i] %in% energy_cfs) { 
    type <- "Energy-based"  
  }
  
  cf <- cost(cost_function = cfs[i])
  
  slope_scenario_1 <- seq(-0.5, 0.5, 0.01)
  vals_scenario_1 <- 1/cf(slope_scenario_1)
  vals_scenario_1 <- (vals_scenario_1 - min(vals_scenario_1)) / (max(vals_scenario_1) - min(vals_scenario_1))
  
  slope_scenario_2 <- seq(-1, 1, 0.01)
  vals_scenario_2 <- 1/cf(slope_scenario_2)
  vals_scenario_2 <- (vals_scenario_2 - min(vals_scenario_2)) / (max(vals_scenario_2) - min(vals_scenario_2))
  
  
  cf_df[[i]] <- rbind(
    data.frame(slope = slope_scenario_1, cf_vals = vals_scenario_1, cf = cfs[i], type = type, scenario = 1),
    data.frame(slope = slope_scenario_2, cf_vals = vals_scenario_2, cf = cfs[i], type = type, scenario = 2))
  
}

cf_df <- do.call(rbind, cf_df)
cf_df$type <- factor(cf_df$type, levels = c("Time-based", "Energy-based"))
cf_df$cf <- tools::toTitleCase(cf_df$cf)
cf_df$cf <- factor(cf_df$cf, levels = tools::toTitleCase(cfs))
cf_df$slope <- slope2deg(cf_df$slope)
cf_df$scenario <- factor(cf_df$scenario)

methods_plot3a <- ggplot(cf_df) + 
  geom_vline(xintercept = 0, lty = 2, colour = "grey80", size = 0.2) +
  geom_line(aes(x = slope, y = cf_vals, group = scenario, colour = type, lty = scenario)) +
  facet_wrap(~cf, nrow = 5, scales = "free") + 
  scale_colour_manual(values = c("#003f5c", "#bc5090")) + 
  labs(x = "Degrees Slope", y = "Normalised cost", colour = "Hypothesis type") + 
  scale_x_continuous(breaks = seq(-50, 50, 10)) + 
  theme_classic() + 
  theme(legend.position = "bottom", strip.text.x = element_text(size = 8))

ggplot2::ggsave(filename = "methods_plot_03.png", plot = methods_plot3a, path = "./Output/Figures/", dpi = 300, width = 7, height = 8)

#### METHODS PLOT FOUR #####

methods_ts_1 <- rbind(data.frame(dist = rexp(n = 1000, rate = 0.5) + runif(n = 1000, min = 0, 5),
                          type = "Time"),
               data.frame(dist = rexp(n = 1000, rate = 0.5) + runif(n = 1000, min = 0, 5),
                          type = "Energy"),
               data.frame(dist = rnorm(n = 1000, mean = 100, sd = 2) + runif(n = 1000, min = 0, max = 50),
                          type = "Difference"))

methods_ts_1$scenario <- "1"

methods_ts_1_arrows <-
  data.frame(
    x1 = c(50, 100, 100),
    x2 = c(110, 25, 25),
    y1 = c(0.8, 2.35, 2.65),
    y2 = c(1.2, 2, 3)
  )

methods_plot4a <- ggplot(methods_ts_1) + 
  geom_boxplot(aes(x = dist,y = type, fill = type)) + 
  xlim(c(0, 200)) +
  labs(x = "Euclidean distance between least-cost paths", y =  NULL, fill = "Hypothesis Type", title = "A") +
  theme_bw() + 
  theme(legend.position = "bottom") + 
  guides(fill = guide_legend(reverse=T)) + 
  annotate("text", x = 50, y = 0.65, label = "Not\nunderdetermined") + 
  annotate("text", x = 100, y = 2.5, label = "Hypotheses\nare robust") + 
  geom_curve(data = methods_ts_1_arrows[1,], aes(x = x1, y = y1, xend = x2, yend = y2),
             arrow = arrow(length = unit(0.08, "inch")), size = 1,
             color = c("#F8766D"), curvature = -0.4) + 
  geom_curve(data = methods_ts_1_arrows[2,], aes(x = x1, y = y1, xend = x2, yend = y2),
             arrow = arrow(length = unit(0.08, "inch")), size = 1,
             color = c("#00BA38"), curvature = -0.4) + 
  geom_curve(data = methods_ts_1_arrows[3,], aes(x = x1, y = y1, xend = x2, yend = y2),
             arrow = arrow(length = unit(0.08, "inch")), size = 1,
             color = c("#619CFF"), curvature = 0.4)

methods_ts_2 <- rbind(data.frame(dist = rnorm(n = 1000, mean = 100, sd = 2) + runif(n = 1000, min = 0, max = 50),
                          type = "Time"),
               data.frame(dist = rnorm(n = 1000, mean = 100, sd = 2) + runif(n = 1000, min = 0, max = 50),
                          type = "Energy"),
               data.frame(dist = rnorm(n = 1000, mean = 100, sd = 2) + runif(n = 1000, min = 0, max = 50),
                          type = "Different"))

methods_ts_2$scenario <- "2"

methods_ts_2_arrows <-
  data.frame(
    x1 = c(40, 60, 60),
    x2 = c(110, 110, 110),
    y1 = c(0.8, 2.55, 2.85),
    y2 = c(1.2, 2.2, 3.2)
  )

methods_plot4b <- ggplot(methods_ts_2) + 
  geom_boxplot(aes(x = dist,y = type, fill = type)) + 
  xlim(c(0, 200)) +
  labs(x = "Euclidean distance (m) between least-cost paths", y =  NULL, fill = "Hypothesis Type", title = "B") +
  theme_bw() + 
  theme(legend.position = "bottom") + 
  guides(fill = guide_legend(reverse=T)) + 
  annotate("text", x = 40, y = 0.6, label = "Not\nunderdetermined") + 
  annotate("text", x = 60, y = 2.7, label = "Hypotheses are\nnot robust") + 
  geom_curve(data = methods_ts_2_arrows[1,], aes(x = x1, y = y1, xend = x2, yend = y2),
             arrow = arrow(length = unit(0.08, "inch")), size = 1,
             color = c("#F8766D"), curvature = -0.8) + 
  geom_curve(data = methods_ts_2_arrows[2,], aes(x = x1, y = y1, xend = x2, yend = y2),
             arrow = arrow(length = unit(0.08, "inch")), size = 1,
             color = c("#00BA38"), curvature = 0.4) + 
  geom_curve(data = methods_ts_2_arrows[3,], aes(x = x1, y = y1, xend = x2, yend = y2),
             arrow = arrow(length = unit(0.08, "inch")), size = 1,
             color = c("#619CFF"), curvature = -0.4)

methods_ts_3 <- rbind(data.frame(dist = rexp(n = 1000, rate = 0.5) * 5,
                          type = "Difference"),
               data.frame(dist = rexp(n = 1000, rate = 0.5) + runif(n = 1000, min = 0, 5),
                          type = "Time"),
               data.frame(dist = rexp(n = 1000, rate = 0.5) + runif(n = 1000, min = 0, 5),
                          type = "Energy"))

methods_ts_3$scenario <- "3"

methods_ts_3_arrows <-
  data.frame(
    x1 = c(150, 120, 120),
    x2 = c(30, 25, 25),
    y1 = c(0.8, 2.55, 2.85),
    y2 = c(1.2, 2.2, 3.2)
  )

methods_plot4c <- ggplot(methods_ts_3) + 
  geom_boxplot(aes(x = dist,y = type, fill = type)) + 
  xlim(c(0, 200)) +
  labs(x = "Euclidean distance (m) between least-cost paths", y =  NULL, fill = "Hypothesis Type", title = "C") +
  theme_bw() + 
  theme(legend.position = "bottom") + 
  guides(fill = guide_legend(reverse=T)) + 
  annotate("text", x = 140, y = 0.75, label = "Underdetermined") + 
  annotate("text", x = 120, y = 2.7, label = "Hypotheses\nare robust") + 
  geom_curve(data = methods_ts_3_arrows[1,], aes(x = x1, y = y1, xend = x2, yend = y2),
             arrow = arrow(length = unit(0.08, "inch")), size = 1,
             color = c("#F8766D"), curvature = 0.8) + 
  geom_curve(data = methods_ts_3_arrows[2,], aes(x = x1, y = y1, xend = x2, yend = y2),
             arrow = arrow(length = unit(0.08, "inch")), size = 1,
             color = c("#00BA38"), curvature = -0.4) + 
  geom_curve(data = methods_ts_3_arrows[3,], aes(x = x1, y = y1, xend = x2, yend = y2),
             arrow = arrow(length = unit(0.08, "inch")), size = 1,
             color = c("#619CFF"))

methods_ts_4 <- rbind(data.frame(dist = rexp(n = 1000, rate = 0.5) * 5,
                          type = "Difference"),
               data.frame(dist = rnorm(n = 1000, mean = 100, sd = 2) + runif(n = 1000, min = 0, max = 50),
                          type = "Time"),
               data.frame(dist = rnorm(n = 1000, mean = 100, sd = 2) + runif(n = 1000, min = 0, max = 50),
                          type = "Energy"))

methods_ts_4$scenario <- "4"

methods_ts_4_arrows <-
  data.frame(
    x1 = c(80, 50, 50),
    x2 = c(20, 110, 110),
    y1 = c(1.45, 2.55, 2.85),
    y2 = c(1.3, 2.2, 3.2)
  )

methods_plot4d <- ggplot(methods_ts_4) + 
  geom_boxplot(aes(x = dist,y = type, fill = type)) + 
  xlim(c(0, 200)) +
  labs(x = "Euclidean distance (m) between least-cost paths", y =  NULL, fill = "Hypothesis Type", title = "D") +
  theme_bw() + 
  theme(legend.position = "bottom") + 
  guides(fill = guide_legend(reverse=T)) + 
  annotate("text", x = 80, y = 1.35, label = "Underdetermined") + 
  annotate("text", x = 50, y = 2.7, label = "Hypotheses are\nnot robust") + 
  geom_curve(data = methods_ts_4_arrows[1,], aes(x = x1, y = y1, xend = x2, yend = y2),
             arrow = arrow(length = unit(0.08, "inch")), size = 1,
             color = c("#F8766D"), curvature = 0.4) + 
  geom_curve(data = methods_ts_4_arrows[2,], aes(x = x1, y = y1, xend = x2, yend = y2),
             arrow = arrow(length = unit(0.08, "inch")), size = 1,
             color = c("#00BA38"), curvature = 0.4) + 
  geom_curve(data = methods_ts_4_arrows[3,], aes(x = x1, y = y1, xend = x2, yend = y2),
             arrow = arrow(length = unit(0.08, "inch")), size = 1,
             color = c("#619CFF"), curvature = -0.4)

methods_plot4e <- gridExtra::grid.arrange(methods_plot4a, methods_plot4b, methods_plot4c, methods_plot4d, ncol = 2)

ggsave("./Output/Figures/methods_plot_04.png", methods_plot4e, dpi = 300, width = 12, height = 10)

##### RESULTS PLOT ONE #####

sim_lcps_scenario_1 <- rbind(
  read.csv("C:/Users/Joe/Documents/28_06_2021/Cost_Functions_and_Robustness/Output/Data/simulated_lcps_220_scenario_1.csv"),
  read.csv("C:/Users/Joe/Documents/28_06_2021/Cost_Functions_and_Robustness/Output/Data/simulated_lcps_230_scenario_1.csv"),
  read.csv("C:/Users/Joe/Documents/28_06_2021/Cost_Functions_and_Robustness/Output/Data/simulated_lcps_240_scenario_1.csv"),
  read.csv("C:/Users/Joe/Documents/28_06_2021/Cost_Functions_and_Robustness/Output/Data/simulated_lcps_250_scenario_1.csv"),
  read.csv("C:/Users/Joe/Documents/28_06_2021/Cost_Functions_and_Robustness/Output/Data/simulated_lcps_260_scenario_1.csv"))

sim_lcps_scenario_1$type[sim_lcps_scenario_1$cf1_type == "Time-based"  & sim_lcps_scenario_1$cf2_type == "Time-based"] <- "Time"
sim_lcps_scenario_1$type[sim_lcps_scenario_1$cf1_type == "Energy-based"  & sim_lcps_scenario_1$cf2_type == "Energy-based"] <- "Energy"
sim_lcps_scenario_1$type[sim_lcps_scenario_1$cf1_type != sim_lcps_scenario_1$cf2_type] <- "Different"
sim_lcps_scenario_1$fractal <- factor(sim_lcps_scenario_1$fractal, levels = (c(2.2, 2.3, 2.4, 2.5, 2.6)))
sim_lcps_scenario_1$type <- factor(sim_lcps_scenario_1$type, levels = c("Different", "Energy", "Time"))
sim_lcps_scenario_1$scenario <- 1

sim_lcps_scenario_2 <- rbind(
  read.csv("C:/Users/Joe/Documents/28_06_2021/Cost_Functions_and_Robustness/Output/Data/simulated_lcps_220_scenario_2.csv"),
  read.csv("C:/Users/Joe/Documents/28_06_2021/Cost_Functions_and_Robustness/Output/Data/simulated_lcps_230_scenario_2.csv"),
  read.csv("C:/Users/Joe/Documents/28_06_2021/Cost_Functions_and_Robustness/Output/Data/simulated_lcps_240_scenario_2.csv"),
  read.csv("C:/Users/Joe/Documents/28_06_2021/Cost_Functions_and_Robustness/Output/Data/simulated_lcps_250_scenario_2.csv"),
  read.csv("C:/Users/Joe/Documents/28_06_2021/Cost_Functions_and_Robustness/Output/Data/simulated_lcps_260_scenario_2.csv"))

sim_lcps_scenario_2$type[sim_lcps_scenario_2$cf1_type == "Time-based"  & sim_lcps_scenario_2$cf2_type == "Time-based"] <- "Time"
sim_lcps_scenario_2$type[sim_lcps_scenario_2$cf1_type == "Energy-based"  & sim_lcps_scenario_2$cf2_type == "Energy-based"] <- "Energy"
sim_lcps_scenario_2$type[sim_lcps_scenario_2$cf1_type != sim_lcps_scenario_2$cf2_type] <- "Different"
sim_lcps_scenario_2$fractal <- factor(sim_lcps_scenario_2$fractal, levels = (c(2.2, 2.3, 2.4, 2.5, 2.6)))
sim_lcps_scenario_2$type <- factor(sim_lcps_scenario_2$type, levels = c("Different", "Energy", "Time"))
sim_lcps_scenario_2$scenario <- 2

sim_lcps_scenario_3 <- rbind(sim_lcps_scenario_1, sim_lcps_scenario_2)

x <- terra::xyFromCell(dem_220_scenario_1, sim_lcps_scenario_3$fromCell)
y <- terra::xyFromCell(dem_220_scenario_1, sim_lcps_scenario_3$toCell)

xy3 <- (x[,1] - y[,1])^2
xy4 <- (x[,2] - y[,2])^2

sim_lcps_scenario_3$euc_dist <- sqrt(xy3 + xy4)

results_plot1a <- ggplot(sim_lcps_scenario_3[sim_lcps_scenario_3$scenario == 1,]) + 
  geom_boxplot(aes(x = sqrt(dist),y = type, fill = type)) + 
  facet_wrap(~fractal, nrow = 5) + 
  labs(x = "Square Root Euclidean distance (m)\nbetween least-cost paths", y = "Hypothesis Type", fill = "Hypothesis Type", title = "A") +
  theme_bw() + 
  theme(legend.position = "bottom") + 
  guides(fill = guide_legend(reverse=T))

results_plot1b <- ggplot(sim_lcps_scenario_3[sim_lcps_scenario_3$scenario == 2,]) + 
  geom_boxplot(aes(x = sqrt(dist),y = type, fill = type)) + 
  facet_wrap(~fractal, nrow = 5) + 
  labs(x = "Square Root Euclidean distance (m)\nbetween least-cost paths", y = "Hypothesis Type", fill = "Hypothesis Type", title = "B") +
  theme_bw() + 
  theme(legend.position = "bottom") + 
  guides(fill = guide_legend(reverse=T))

results_plot1c <- gridExtra::grid.arrange(results_plot1a, results_plot1b, ncol = 2)

ggplot2::ggsave("./Output/Figures/results_plot_01.png", results_plot1c, dpi = 300, width = 10, height = 7)

##### RESULTS PLOT TWO #####

results_plot2a <- ggplot(sim_lcps_scenario_3[sim_lcps_scenario_3$scenario == 1,]) + 
  geom_boxplot(aes(x = fractal,y = sqrt(dist), fill = type)) + 
  facet_wrap(~factor(type, levels = c("Time", "Energy", "Different")), nrow = 3) + 
  labs(x = "Fractal Dimension", y = "Square Root Euclidean distance (m)\nbetween least-cost paths", fill = "Hypothesis Type", title = "A") +
  ylim(c(0, 20)) + 
  theme_bw() + 
  theme(legend.position = "bottom") + 
  guides(fill = guide_legend(reverse=T))

results_plot2b <- ggplot(sim_lcps_scenario_3[sim_lcps_scenario_3$scenario == 2,]) + 
  geom_boxplot(aes(x = fractal,y = sqrt(dist), fill = type)) + 
  facet_wrap(~factor(type, levels = c("Time", "Energy", "Different")), nrow = 3) +
  labs(x = "Fractal Dimension", y = "Square Root Euclidean distance (m)\nbetween least-cost paths", fill = "Hypothesis Type", title = "B") +
  ylim(c(0, 20)) + 
  theme_bw() + 
  theme(legend.position = "bottom") + 
  guides(fill = guide_legend(reverse=T))

results_plot2c <- gridExtra::grid.arrange(results_plot2a, results_plot2b, ncol = 2)

ggplot2::ggsave("./Output/Figures/results_plot_02.png", results_plot2c, dpi = 300, width = 10, height = 7)

##### RESULTS PLOT THREE #####

sim_lcps_scenario_3$type <- factor(sim_lcps_scenario_3$type, levels = c("Different", "Energy", "Time"))
sim_lcps_scenario_3$scenario <- factor(sim_lcps_scenario_3$scenario)

### need to alter colour scheme to be reversed... currently red for time etc. needs to be blue

sim_lcps_scenario_3$alpha <- NA
sim_lcps_scenario_3$alpha[sim_lcps_scenario_3$type == "Time"] <- 0.0025
sim_lcps_scenario_3$alpha[sim_lcps_scenario_3$type == "Energy"] <- 0.05
sim_lcps_scenario_3$alpha[sim_lcps_scenario_3$type == "Different"] <- 0.0025

results_plot3a <- ggplot(sim_lcps_scenario_3, aes(x = euc_dist, y = dist, colour = type)) + 
  geom_point(aes(colour = type), alpha = sim_lcps_scenario_3$alpha, show.legend = FALSE) +
  geom_smooth(method='lm', aes(linetype = scenario), se = FALSE) + 
  facet_grid(~factor(type, levels = c("Time", "Energy", "Different"))~fractal, scales = "free_y") + 
  labs(x = "Euclidean distance (m) from origin to destination", y = "Euclidean distance (m) between least-cost paths", colour = "Hypothesis Type", title = "A") + 
  theme_bw() + 
  theme(legend.position = "bottom",
        plot.margin = unit(c(1,1,1,1), "cm")) + 
  theme(legend.position = "bottom") + 
  guides(colour = guide_legend(order = 1, reverse=T), linetype = guide_legend(order = 2))

ggplot2::ggsave("./Output/Figures/results_plot_03.png", results_plot3a, dpi = 300, width = 12, height = 8)

#### CASE STUDY TWO: SARDINIA ####

dem <- terra::rast("./Data/Case_study/SARD_SUL_CAR_10m.tif")
dem[dem <= 0] <- NA

locs <- sf::st_read("./Data/Case_study/SW_sites.shp")

for (i in 1:length(cfs)) {
  
  print(paste0("i = ", i))
  
  lcps <- list()
  
  for(j in 1:sard_nsims) { 
    
    print(paste0("i = ", j))
    
    print("adding error")
    dem_error <- leastcostpath::add_dem_error(x = dem, rmse = 4.3, type = "u")
    print("calculating slope_cs")
    slope_cs <- leastcostpath::create_slope_cs(x = dem_error, cost_function = cfs[i], neighbours = 16)
    print("calculating lcp")
    lcps[[j]] <- leastcostpath::create_lcp(x = slope_cs, origin = locs[1,], destination = locs[2,])
    
  }
  
  lcps <- do.call(rbind, lcps)
  sf::write_sf(lcps, paste0("./Output/LCPs/", cfs[i], "_lcp.shp"))
}

time_based_cfs <- c("tobler", "tobler offpath", "modified tobler", "irmischer-clarke male", "irmischer-clarke offpath male", "irmischer-clarke female", "irmischer-clarke offpath female", 'rees', "davey", 'garmy', 'kondo-saino', 'naismith', "campbell 2019 50")

energy_based_cfs <- c("herzog", "llobera-sluckin")

lcp_files <- list.files(path = "./Output/LCPs/", pattern = "\\.shp$", full.names = TRUE, recursive = TRUE)

time_based_lcps <- do.call(rbind, lapply(lcp_files[apply(sapply(time_based_cfs, grepl, lcp_files), 1, sum) != 0], sf::st_read))
energy_based_lcps <- do.call(rbind, lapply(lcp_files[apply(sapply(energy_based_cfs, grepl, lcp_files), 1, sum) != 0], sf::st_read))

time_based_lcp_dens <- lcp_density(lcps = time_based_lcps, dem = dem)
energy_based_lcp_dens <- lcp_density(lcps = energy_based_lcps, dem = dem)

time_based_lcp_dens <- time_based_lcp_dens / nrow(time_based_lcps)
energy_based_lcp_dens <- energy_based_lcp_dens / nrow(energy_based_lcps)

time_based_lcp_dens[time_based_lcp_dens == 0] <- NA
energy_based_lcp_dens[energy_based_lcp_dens == 0] <- NA

terra::writeRaster(time_based_lcp_dens, "./Output/LCPs/time_based_lcp_density.tif")
terra::writeRaster(energy_based_lcp_dens, "./Output/LCPs/energy_based_lcp_density.tif")

energy_based_lcps_herzog <- energy_based_lcps[energy_based_lcps$cstFnct == "herzog",]
energy_based_lcps_llobera_sluckin <- energy_based_lcps[energy_based_lcps$cstFnct == "llobera-sluckin",]

energy_based_lcps_herzog_dens <- lcp_density(lcps = energy_based_lcps_herzog, dem = dem)
energy_based_lcps_llobera_sluckin_dens <- lcp_density(lcps = energy_based_lcps_llobera_sluckin, dem = dem)

energy_based_lcps_herzog_dens <- energy_based_lcps_herzog_dens / nrow(energy_based_lcps_herzog)
energy_based_lcps_llobera_sluckin_dens <- energy_based_lcps_llobera_sluckin_dens / nrow(energy_based_lcps_llobera_sluckin)

energy_based_lcps_herzog_dens[energy_based_lcps_herzog_dens == 0] <- NA
energy_based_lcps_llobera_sluckin_dens[energy_based_lcps_llobera_sluckin_dens == 0] <- NA

terra::writeRaster(energy_based_lcps_herzog_dens, "./Output/LCPs/energy_based_lcps_herzog_dens.tif")
terra::writeRaster(energy_based_lcps_llobera_sluckin_dens, "./Output/LCPs/energy_based_lcps_llobera_sluckin_dens.tif")
