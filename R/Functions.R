normalise_raster <- function(x){(x-terra::minmax(x)[1])/(terra::minmax(x)[2]-terra::minmax(x)[1])}

euclidean_distance <- function(lcp, comparison) {
  
  lcp_pts <- sf::st_cast(lcp, "POINT", warn = FALSE)
  
  distance <- max(sf::st_distance(x = lcp_pts, y = comparison))
  distance <- as.numeric(distance)
  
  return(distance)
  
}

compare_lcps <- function(cs_list, cfs, cf1_index, cf2_index, pts) {
  
  slope_cs1 <- cs_list[[cf1_index]]
  slope_cs2 <- cs_list[[cf2_index]]
  
  lcp1 <- leastcostpath::create_lcp(x = slope_cs1, origin = pts[1,], destination = pts[2,])
  lcp2 <- leastcostpath::create_lcp(x = slope_cs2, origin = pts[1,], destination = pts[2,])
  
  dist <- euclidean_distance(lcp = lcp1, comparison = lcp2)
  
  lcp1$dist <- dist
  lcp1$cf1 <- cfs[cf1_index]
  lcp1$cf2 <- cfs[cf2_index]
  lcp1$cf1_type <- NA
  lcp1$cf2_type <- NA

  time_cfs <- c("tobler", "rees", "davey", "irmischer-clarke male", "irmischer-clarke offpath male", "irmischer-clarke female", "irmischer-clarke offpath female", "modified tobler", "garmy", "kondo-saino", "naismith", "campbell 2019 50")
  
  energy_cfs <- c("herzog", "llobera-sluckin")
  
  if(lcp1$cf1 %in% time_cfs) {
    lcp1$cf1_type <- "Time-based"
  } else if (lcp1$cf1 %in% energy_cfs) {
    lcp1$cf1_type <- "Energy-based"
  }
  
  if(lcp1$cf2 %in% time_cfs) {
    lcp1$cf2_type <- "Time-based"
  } else if (lcp1$cf2 %in% energy_cfs) {
    lcp1$cf2_type <- "Energy-based"
  }
  
  return(lcp1)
  
}

cost <- function(cost_function, crit_slope, percentile) {
  
  cfs <- c("tobler", "tobler offpath", "davey", 'rees', "irmischer-clarke male", "irmischer-clarke offpath male", "irmischer-clarke female", "irmischer-clarke offpath female", "modified tobler", 'garmy', 'kondo-saino', "wheeled transport", "herzog", "llobera-sluckin", "naismith", "minetti", "campbell", "campbell 2019 50", "sullivan")
  
  if (inherits(cost_function, "character")) {
    if (!cost_function %in% cfs) {
      stop("cost_function argument is invalid. See details for accepted cost functions")
    }
    
    # mathematical slope to degrees
    slope2deg <- function(slope) { (atan(slope) * 180/pi)}
    
    # degrees to radians
    deg2rad <- function(deg) {(deg * pi) / (180)}
    
    # Tobler Hiking Function measured in km/h. Divide by 3.6 to turn into m/s
    if (cost_function == "tobler") {
      
      # 3.6 converts from km/h to m/s
      cf <- function(x) {
        (6 * exp(-3.5 * abs(x + 0.05))) / 3.6
      }
    }
    
    if (cost_function == "tobler offpath") {
      
      cf <- function(x) {
        ((6 * exp(-3.5 * abs(x + 0.05))) * 0.6) / 3.6
      }
      
    }
    
    # requires mathematical slope to be in radians. Necessary to convert mathematical slope to degrees and then to radians
    if (cost_function == "davey") { 
      
      cf <- function(x, y = 1.40, z = 2.8) { 
        x_deg <- slope2deg(x)
        x_rad <- deg2rad(x_deg)
        
        (y * exp(-z*abs(x_rad)))
      }
      
    }
    
    if (cost_function == "rees") { 
      
      cf <- function(x) { 
        (1 / (0.75 + 0.09 * abs(x) + 14.6 * (abs(x))^2))
      }
      
    }
    
    if (cost_function == "irmischer-clarke male") {
      
      cf <- function(x) {
        (0.11 + exp(-(abs(x) * 100 + 5)^2/(2 * 30^2)))
        
      }
      
    }
    
    if (cost_function == "irmischer-clarke offpath male") {
      
      cf <- function(x) {
        (0.11 + 0.67 * exp(-(abs(x) * 100 + 2)^2/(2 * 30^2)))
      }
      
    }
    
    if (cost_function == "irmischer-clarke female") {
      
      cf <- function(x) {
        (0.95 * (0.11 + exp(-(abs(x) * 100 + 5)^2/(2 * 30^2))))
      }
      
    }
    
    if (cost_function == "irmischer-clarke offpath female") {
      
      cf <- function(x) {
        0.95 * (0.11 + 0.67 * exp(-(abs(x) * 100 + 2)^2/(2 * 30^2)))
      }
      
    }
    
    if (cost_function == "modified tobler") {
      
      cf <- function(x) {
        (4.8 * exp(-5.3 * abs((x * 0.7) + 0.03)))
      }
      
    }
    
    if (cost_function == "garmy") {
      
      cf <- function(x) { 
        (4 * exp(-0.008 * (slope2deg(x)^2)))
      }
      
    }
    
    if (cost_function == "kondo-saino") {
      
      cf <- function(x) { 
        ifelse(abs(x) >= -0.07, (5.1 * exp(-2.25 * abs(x + 0.07))), (5.1 * exp(-1.5 * abs(x + 0.07))))
      }
      
    }
    
    if (cost_function == "wheeled transport") {
      
      cf <- function(x) {
        (1/(1 + ((abs(x) * 100)/crit_slope)^2))
      }
      
    }
    
    if (cost_function == "herzog") {
      
      cf <- function(x) {
        
        (1/((1337.8 * (x)^6) + (278.19 * (x)^5) - (517.39 * (x)^4) - (78.199 * (x)^3) + (93.419 * (x)^2) + (19.825 * (x)) + 1.64))
        
      }
      
    }
    
    if (cost_function == "llobera-sluckin") {
      
      cf <- function(x) {
        (1/(2.635 + (17.37 * (x)) + (42.37 * (x)^2) - (21.43 * (x)^3) + (14.93 * (x)^4)))
      }
      
    }
    
    if (cost_function == "naismith") {
      
      cf <- function(x) {
        
        x_deg <- slope2deg(x)
        x_rad <- deg2rad(x_deg)
        
        ifelse(x_deg > 0, (1 / (0.72 + 6 * tan(x_rad))), ifelse(x_deg <= -12, (1 / (0.72 - 2 * tan(x_rad))), ifelse(x_deg > -12 & x_deg <= -5, (1 / (0.72 + 2 * tan(x_rad))), 1.40)))
      }
    }
    
    if (cost_function == "minetti") { 
      
      cf <- function(x) {
        (1/((280.5 * abs(x)^5) - (58.7 * abs(x)^4) - (76.8 * abs(x)^3) + (51.9 * abs(x)^2) + (19.6 * abs(x)) + 2.5))
      }
    }
    
    if (cost_function == "campbell") { 
      
      cf <- function(x) { 
        x_deg <- slope2deg(x)
        
        1.662 - (5.191 * 10^-3) * x_deg - (1.127*10^-3) * x_deg^2
      }
      
    }
    
    if (cost_function == "campbell 2019 50") {
      
      cf <- function(x) { 
        (63.66*(1/(pi*10.064*(1+((slope2deg(x)--2.171)/10.064)^2)))+0.628)+-0.00463*slope2deg(x)
      }
    }
    
    if (cost_function == "sullivan") {
      
      percentile_choice <- c(0.167, 0.5, 0.833)
      
      if (!percentile %in% percentile_choice) {
        stop("percentile argument is invalid. Expecting percentile value of 0.167, 0.5, 0.833")
      }
      
      term_a <- c(-3.3717, -2.8292, -2.2893)
      term_b <- c(25.8255, 20.9482, 19.4024)
      term_c <- c(92.6594, 77.6346, 65.3577)
      term_d <- c(-0.1624, 0.2228, 0.6226)
      term_e <- c(0.0019, -0.0004, -0.0020)
      
      lorentz_function_terms <- data.frame(percentile = percentile_choice, term_a, term_b, term_c, term_d, term_e)
      terms <- lorentz_function_terms[lorentz_function_terms$percentile == percentile, ]
      
      cf <- function(x){
        (terms$term_c*(1/(pi*terms$term_b*(1+((slope2deg(x)-terms$term_a)/terms$term_b)^2)))+terms$term_d)+terms$term_e*slope2deg(x)
      }
      
    }
    
  }
  
  if(is.function(cost_function)) {
    
    cf <- cost_function
    
  }
  
  return(cf)
}

calculate_distance <- function(x, adj) { 
  
  xy1 <- terra::xyFromCell(x, adj[, 1])
  xy2 <- terra::xyFromCell(x,adj[, 2])
  
  xy3 <- (xy1[,1] - xy2[,1])^2
  xy4 <- (xy1[,2] - xy2[,2])^2
  
  dist <- sqrt(xy3 + xy4)
  
  return(dist)
  
}

calculate_slope <- function(x, neighbours) { 
  
  cells <- which(!is.na(terra::values(x)))
  na_cells <- which(is.na(terra::values(x)))
  
  adj <- terra::adjacent(x = x, cells = cells, directions = neighbours, pairs = TRUE)
  adj <- adj[!adj[,2] %in% na_cells,]
  
  elev_values <- terra::values(x)[,1]
  
  message("calculating slope...")
  
  rise <- (elev_values[adj[,2]] - elev_values[adj[,1]])
  run <- calculate_distance(x = x, adj = adj)
  
  mathematical_slope <- rise/run
  
  return(mathematical_slope)
  
}

slope2deg <- function(slope) { (atan(slope) * 180/pi)}

lcp_density <- function(lcps, dem) { 
  
  dem_base <- dem
  dem_base[] <- 0
  
  res <- max(terra::res(dem)) / 2
  
  for (i in 1:nrow(lcps)) { 
    
    print(i)
    
    pts <- sf::st_line_sample(x = lcps[i,], density = 1/res, type = "regular")
    pts_coord <- sf::st_coordinates(pts)[,1:2]
    
    cells <- terra::cellFromXY(object = dem, xy = pts_coord)
    cells <- unique(cells)
    
    dem_base[cells] <- dem_base[cells] + 1
    
  }
  
  return(dem_base)
  
}