rm(list = ls())

library(sf);library(tmap); library(dplyr); library(SpatialEpi); library(ggplot2); library(polyCub); library(raster)
library(sp); library(spdep); library(exactextractr); library(INLA);

#Source the functions
# SSDAM and SSDEM functions
file_list = c("functions/Krig_fit.R","functions/cov.spatial.R", "functions/help_functions.R")
for (file in file_list) {
  source(file)
}

###################### Underlying continuous surface ###########################
raster_true <- raster(extent(0, 100, 0, 100), res = 1) 

raster_true_poly <- rasterToPolygons(raster_true, dissolve = FALSE)
plot(raster_true_poly, border = "black", col = NA, lwd = 1)


raster_area2 <- raster(extent(0, 100, 0, 100), res = 10)
map_area2 = rasterToPolygons(raster_area2, dissolve = FALSE)
map_area2 = st_as_sf(map_area2)


# Raster for estimation
raster1 <- raster(extent(0, 100, 0, 100), res = 1) 
raster1_poly <- rasterToPolygons(raster1, dissolve = FALSE)
plot(raster1_poly, border = "black", col = NA, lwd = 1)


# Number of simulations
N <- 100
beta0 <- -4 #default baseline
#beta0 <- -6 #reduced baseline


# Simulate from spatial process
sigma2 = 0.7
range = 3
#range = 10
cov.method <- "matern"

spat_range <- seq(0,100, length = 500)

cov_true <- sigma2 * 1/range * spat_range * besselK(1/range * spat_range, 1)

plot(spat_range, cov_true, type = "l")



set.seed(2903)
sim_spat = geoR::grf(n = prod(dim(raster_true)), 
                     grid = xyFromCell(raster_true, 1:ncell(raster_true)), 
                     cov.pars = c(sigma2, range), nsim = N, cov.model = "matern",
                     kappa = 1)


# Plot simulated spatial effect from iteration 1
values(raster_true) = scale(sim_spat$data[,1], scale = F)

limits = c(-3.5,3.5)

ggplot() +
  geom_raster(data = xyFromCell(raster_true, 1:ncell(raster_true)), 
              aes(x = x, y = y, fill = values(raster_true))) +
  geom_sf(data = st_cast(st_union(map_area2), "POLYGON"), fill = NA, color = "black") +
  scale_fill_viridis_c(limits = limits) + 
  theme_minimal() +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA)) +
  ggtitle("True spatial effect")+
  labs(fill = "Spatial effect")


# Linear covariate
linear_cov <- function(x, y) {
  x = x/100*6-3
  y = y/100*6-3  
  return(0.2*(((x-y)^2)/15+ sin(x) * cos(y))+0.5)
}

# True effects
smooth_spat_poly_weighted <- function(linear_xy, spatial_xy) {
  return(exp(beta0 + spatial_xy-1.5*linear_xy))
}

# Covariance functions
cov_function = function(sigma, rho, spat_range){
  if (cov.method == "exponential"){
    return(sigma * exp(-rho *spat_range))
  }else if (cov.method == "spherical"){
    cov_return = NULL
    for (iter in 1:length(spat_range)){
      cov_return[iter] = sigma*ifelse(spat_range[iter]<= 1/rho,
                                      1-1.5*rho*spat_range[iter]+
                                        0.5*rho^3 * spat_range[iter]^3,0)
    }
    return(cov_return)
  }else if(cov.method == "circular"){
    cov_return = NULL  
    for (iter in 1:length(spat_range)){
      theta = rho*spat_range[iter]
      cov_return[iter] = sigma*ifelse(spat_range[iter] <= 1/rho,
                                      1-2/pi*(theta*sqrt(1-theta^2)+asin(theta)),0)
      
    }
    return(cov_return)
  }else if (cov.method == "matern"){
    return(exp(-rho * spat_range)*(1+rho * spat_range) * sigma)
  }
}

# Population size (same in every cell)
pop_den <- raster_true
values(pop_den) <- 1
rstr_centers = xyFromCell(raster_true, 1:ncell(raster_true))


# Find knots with space filling algorithm

mtime_knots2 <- proc.time()
knots = find_knots(map_area2, nd = 200)
time_knots2 <- (proc.time()-mtime_knots2)[3]

    
    
    # Save the results
    
    rho_raster <- sigma_raster <- sigma_error_raster <- time_raster <- numeric(N)
    cov_raster <- cor_raster <- matrix(0, nrow = N, ncol = 500)  
    CI_spat_raster <- Bias_spat_raster <- RMSE_spat_raster <- matrix(0, nrow = N, ncol = prod(dim(raster_true)))  
    CI_continuous_incidence <- Bias_continuous_incidence <- RMSE_continuous_incidence <- matrix(0, nrow = N, ncol = prod(dim(raster_true)))  
    CI_area <- Bias_area <- RMSE_area <- matrix(0, nrow = N, ncol = dim(map_area2)[1])  
    
    rho_disaggr <- sigma_disaggr <- time_disaggr <- numeric(N)
    cov_disaggr <- cor_disaggr <- matrix(0, nrow = N, ncol = 500)  
    CI_spat_disaggr <- Bias_spat_disaggr <- RMSE_spat_disaggr <- matrix(0, nrow = N, ncol = prod(dim(raster_true)))  
    CI_area_disaggr <- Bias_area_disaggr <- RMSE_area_disaggr <- matrix(0, nrow = N, ncol = dim(map_area2)[1])  
    CI_continuous_incidence_disaggr <- Bias_continuous_incidence_disaggr <- RMSE_continuous_incidence_disaggr <- matrix(0, nrow = N, ncol = prod(dim(raster_true)))  
    
    rho_exact <- sigma_exact <- sigma_error_exact <- time_exact <- numeric(N)
    cov_exact <- cor_exact <- matrix(0, nrow = N, ncol = 500)  
    CI_spat_exact <- Bias_spat_exact <- RMSE_spat_exact <- matrix(0, nrow = N, ncol = prod(dim(raster_true)))  
    CI_continuous_incidence_exact <- Bias_continuous_incidence_exact <- RMSE_continuous_incidence_exact <- matrix(0, nrow = N, ncol = prod(dim(raster_true)))  
    CI_area_exact <- Bias_area_exact <- RMSE_area_exact <- matrix(0, nrow = N, ncol = dim(map_area2)[1])  

    
    for (ii in 1:N){
      tryCatch({
      
        set.seed(ii)
        
        # Rasters at fine scale
        rstr_linear = raster_true
        values(rstr_linear) = 0 # No covariates (all values are 0)
        
        rstr_true_spat = raster_true
        values(rstr_true_spat) = scale(sim_spat$data[,ii], scale = F)
        
        rstr_true_incidence = raster_true
        values(rstr_true_incidence) = smooth_spat_poly_weighted(values(rstr_linear),
                                                                values(rstr_true_spat))
        
        rstr_true_cases = rstr_true_incidence
        values(rstr_true_cases) = rpois(n = length(values(rstr_true_cases)), 
                                        lambda = values(rstr_true_incidence))
        
        # Aggregate to largest area level
        true_inc <- exact_extract(rstr_true_incidence, map_area2, 'sum', progress = F)
        
        # Simulate the data
        y <- exact_extract(rstr_true_cases, map_area2, 'sum', progress = F)
        
        # True population size
        true_pop <- exact_extract(pop_den,map_area2, 'sum', progress = F)
        
        
        model <- glm(y ~ offset(log(true_pop)), family = "poisson")
        start_sigma = c(1/exp(1),min(mean(model$residuals^2),1))                                 
        
                                        
       ##############################   SSDAM  ###############################

      # Preparing the data (which cells intersect with which area)
      intersecting_cells_raster <- find_intersecting_cells(map = map_area2,
                                                           rstr_pop_density = pop_den,
                                                           covariate_missing = F)
   

      dat_raster <- prepare_data(y = y, map = map_area2, 
                                 intersecting_cells = intersecting_cells_raster, # output of find_intersecting_cells
                                 exact = F, weighted = F)
      
      
      # Fitting the model 
      mtime_raster<- proc.time()
      model_raster <- Krig.fit(data = dat_raster, # Output of prepare_data
                               cov.method = cov.method,
                               knots = knots, # Output of find_knots
                               offset = true_pop, # Area-level population size
                               exact = F, n_init = 25, parallel = F,
                               sigma_init = start_sigma)
      
      time_raster[ii] <- (proc.time()-mtime_raster)[3]
      
      
      
      
      rho_raster[ii] = exp(model_raster$v_mode[length(model_raster$v_mode)])
      sigma_raster[ii] = 1/exp(model_raster$v_mode[length(model_raster$v_mode)-1])
      sigma_error_raster[ii] = 1/exp(model_raster$v_mode[length(model_raster$v_mode)-2])
      
      cov_raster[ii,] = cov_function(sigma = sigma_raster[ii],
                                     rho = rho_raster[ii],
                                     spat_range = spat_range)
      
      cor_raster[ii,] = cov_function(sigma = 1,
                                     rho = rho_raster[ii],
                                     spat_range = spat_range)
      
      
      #################  spatial component ##########      
      spat_raster = predict_spatial(newdata = xyFromCell(rstr_true_spat, 1:ncell(rstr_true_spat)),
                                    model = model_raster, plot = T,
                                    limits = limits, map =  st_cast(st_union(map_area2), "POLYGON"))
      
      
      ysp.0 <- values(rstr_true_spat)
      
      # Coverage
      ci_spat <- numeric(length(ysp.0))
      for (i in 1:length(ysp.0)) {
        ci_spat[i] <- as.numeric(ifelse(ysp.0[i] >= spat_raster$lower[i] & ysp.0[i] <= spat_raster$upper[i], 1, 0))
      }
      CI_spat_raster[ii,] <- ci_spat
      Bias_spat_raster[ii,] <- ysp.0 - spat_raster$fit
      RMSE_spat_raster[ii,] <- (ysp.0 - spat_raster$fit)^2
      
      
      ################# incidence rate (continuous) ###########
      
      pred_continuous = predict_incidence_continuous(xyFromCell(rstr_true_spat, 1:ncell(rstr_true_spat)), 
                                                     model_raster, plot = F)
      
      # True effect
      true_effect = values(rstr_true_incidence)
      
      
      ci_y0 <- numeric(length(true_effect))
      for (i in 1:length(true_effect)) {
        ci_y0[i] <- as.numeric(ifelse(log(true_effect)[i] >= pred_continuous$lower[i] & log(true_effect)[i] <= pred_continuous$upper[i], 1, 0))
      }
      CI_continuous_incidence[ii,] <- ci_y0
      Bias_continuous_incidence[ii,] <- true_effect-exp(pred_continuous$fit)
      RMSE_continuous_incidence[ii,] <- (true_effect-exp(pred_continuous$fit))^2
      
      
      
      #################  discrete incidence ###########
      
      discrete_raster = predict_discrete_mean(map_area2, model_raster, plot = F)
      
      ci_y0 <- numeric(length(true_inc))
       for (i in 1:length(true_inc)) {
        ci_y0[i] <- as.numeric(ifelse(log(true_inc)[i] >= discrete_raster$lower[i] & log(true_inc)[i] <= discrete_raster$upper[i], 1, 0))
       }
      CI_area[ii,] <- ci_y0
      Bias_area[ii,] <- true_inc - exp(discrete_raster$fit)
      RMSE_area[ii,] <- (true_inc-exp(discrete_raster$fit))^2
      
      
        
      ##############################   SSDEM  ###############################

      dat_exact <- prepare_data(y = y, map = map_area2, 
                                 intersecting_cells = intersecting_cells_raster,
                                 exact = T)
      
      
      # Fit model
      mtime_exact<- proc.time()
      model_exact <- Krig.fit(data = dat_exact,
                               cov.method = cov.method,
                               knots = knots,
                               offset = dat_exact$rstr_centers$pop, # Population size in each intersecting grid cell
                               exact = T, n_init = 25, parallel = F,
                               sigma_init = start_sigma)
      
      time_exact[ii] <- (proc.time()-mtime_exact)[3]
      
      
      
      
      rho_exact[ii] = exp(model_exact$v_mode[length(model_exact$v_mode)])
      sigma_exact[ii] = 1/exp(model_exact$v_mode[length(model_exact$v_mode)-1])
      sigma_error_exact[ii] = 1/exp(model_exact$v_mode[length(model_exact$v_mode)-2])
      
      cov_exact[ii,] = cov_function(sigma = sigma_exact[ii],
                                    rho = rho_exact[ii],
                                    spat_range = spat_range)
      
      cor_exact[ii,] = cov_function(sigma = 1,
                                    rho = rho_exact[ii],
                                    spat_range = spat_range)
      
      
      ################# spatial component ###########

      spat_exact = predict_spatial(newdata = xyFromCell(rstr_true_spat, 1:ncell(rstr_true_spat)),
                                    model = model_exact, plot = T,
                                    limits = limits, map =  st_cast(st_union(map_area2), "POLYGON"))
      
      
      ysp.0 <- values(rstr_true_spat)
      
      # Coverage
      ci_spat <- numeric(length(ysp.0))
      for (i in 1:length(ysp.0)) {
        ci_spat[i] <- as.numeric(ifelse(ysp.0[i] >= spat_exact$lower[i] & ysp.0[i] <= spat_exact$upper[i], 1, 0))
      }
      CI_spat_exact[ii,] <- ci_spat
      Bias_spat_exact[ii,] <- ysp.0 - spat_exact$fit
      RMSE_spat_exact[ii,] <- (ysp.0 - spat_exact$fit)^2
      
      
      ################# incidence rate (continuous) ###########
      
      pred_continuous_exact = predict_incidence_continuous(xyFromCell(rstr_true_spat, 1:ncell(rstr_true_spat)), 
                                                     model_exact, plot = F)
      
      # True effect
      true_effect = values(rstr_true_incidence)
      
      
      ci_y0 <- numeric(length(true_effect))
      for (i in 1:length(true_effect)) {
        ci_y0[i] <- as.numeric(ifelse(log(true_effect)[i] >= pred_continuous_exact$lower[i] & log(true_effect)[i] <= pred_continuous_exact$upper[i], 1, 0))
      }
      CI_continuous_incidence_exact[ii,] <- ci_y0
      Bias_continuous_incidence_exact[ii,] <- true_effect-exp(pred_continuous_exact$fit)
      RMSE_continuous_incidence_exact[ii,] <- (true_effect-exp(pred_continuous_exact$fit))^2
      
      
      
      ################# discrete incidence ###########
      
      discrete_exact = predict_discrete_mean(map_area2, model_exact, plot = F)
      
      ci_y0 <- numeric(length(true_inc))
      for (i in 1:length(true_inc)) {
        ci_y0[i] <- as.numeric(ifelse(log(true_inc)[i] >= discrete_exact$lower[i] & log(true_inc)[i] <= discrete_exact$upper[i], 1, 0))
      }
      CI_area_exact[ii,] <- ci_y0
      Bias_area_exact[ii,] <- true_inc - exp(discrete_exact$fit)
      RMSE_area_exact[ii,] <- (true_inc - exp(discrete_exact$fit))^2
      
      
     },error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
      
      print(ii) 
      
    }

    
    result <- data.frame(
      Metric = c("time","Bias spatial","RMSE spatial", "Coverage spatial",
                 "Bias incidence", "RMSE incidence", "Coverage incidence",
                 "Bias area", "RMSE area", "Coverage area"),
      SSDAM = round(c(mean(time_raster[ind]),mean(Bias_spat_raster[ind,]),
                         sqrt(mean(RMSE_spat_raster[ind,], na.rm = T)), mean(CI_spat_raster[ind,]),
                         mean(Bias_continuous_incidence[ind,]),
                         sqrt(mean(RMSE_continuous_incidence[ind,], na.rm=T)),
                         mean(CI_continuous_incidence[ind,]),
                         mean(Bias_area[ind,]), 
                         sqrt(mean(RMSE_area[ind,], na.rm = T)),
                      mean(CI_area[ind,])),2),
      SSDEM = round(c(mean(time_exact[ind]),mean(Bias_spat_exact[ind,]),
                      sqrt(mean(RMSE_spat_exact[ind,], na.rm = T)), mean(CI_spat_exact[ind,]),
                      mean(Bias_continuous_incidence_exact[ind,]),
                      sqrt(mean(RMSE_continuous_incidence_exact[ind,], na.rm=T)),
                      mean(CI_continuous_incidence_exact[ind,]),
                      mean(Bias_area_exact[ind,]), 
                      sqrt(mean(RMSE_area_exact[ind,], na.rm = T)),
                      mean(CI_area_exact[ind,])),2)
      
    )  
    
    result


