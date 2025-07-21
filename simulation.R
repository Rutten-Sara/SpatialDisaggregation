rm(list = ls())

library(sf);library(tmap); library(dplyr); library(SpatialEpi); library(ggplot2); library(polyCub); library(raster)
library(sp); library(spdep); library(exactextractr); library(INLA);

#Source the functions
# SSDAM and SSDEM functions
file_list = c("functions/Krig_fit.R","functions/cov.spatial.R", "functions/help_functions.R")
for (file in file_list) {
  source(file)
}

wd <- getwd()
# Additional functions disaggregation package, SDALGCP and ATP
source("functions/Other functions/Disaggregation/disaggregation_help_functions.R")


# If SDALGCP not available for your version of R, reading dm_functions.R is sufficient
library(SDALGCP)
#source("functions/Other functions/SDALGCP/dm_functions.R")

setwd('functions/Other functions/ATA-Poisson-CoKriging-main/R')
source('ATA-ATP-Functions.R')

setwd(wd)

###################### Underlying continuous surface ###########################
raster_true <- raster(extent(0, 100, 0, 100), res = 1) 

raster_true_poly <- rasterToPolygons(raster_true, dissolve = FALSE)
plot(raster_true_poly, border = "black", col = NA, lwd = 1)

raster_area1 <- raster(extent(0, 100, 0, 100), res = 20)
map_area1 = rasterToPolygons(raster_area1, dissolve = FALSE)
map_area1 = st_as_sf(map_area1)

raster_area2 <- raster(extent(0, 100, 0, 100), res = 10)
map_area2 = rasterToPolygons(raster_area2, dissolve = FALSE)
map_area2 = st_as_sf(map_area2)


raster_area3 <- raster(extent(0, 100, 0, 100), res = 4)
map_area3 = rasterToPolygons(raster_area3, dissolve = FALSE)
map_area3 = st_as_sf(map_area3)

map_area_list = list(map_area1 = map_area1, map_area2 = map_area2,
                     map_area3 = map_area3)



pdf("Figures/raster_area3.pdf",height=6,width=6)
plot(map_area1, border = "black", col = NA, lwd = 1, main = "Area 3", cex.main = 2)
dev.off()

pdf("Figures/raster_area2.pdf",height=6,width=6)
plot(map_area2, border = "black", col = NA, lwd = 1, main = "Area 2", cex.main = 2)
dev.off()

pdf("Figures/raster_area1.pdf",height=6,width=6)
plot(map_area3, border = "black", col = NA, lwd = 1, main = "Area 1", cex.main = 2)
dev.off()


# Raster for estimation
raster1 <- raster(extent(0, 100, 0, 100), res = 1) 
raster1_poly <- rasterToPolygons(raster1, dissolve = FALSE)
plot(raster1_poly, border = "black", col = NA, lwd = 1)


# Number of simulations
N <- 100
beta0 <- -3


# Simulate from spatial process
sigma2 = 0.7
#range = 3
range = 10
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
  geom_sf(data = st_cast(st_union(map_area1), "POLYGON"), fill = NA, color = "black") +
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
mtime_knots1 <- proc.time()
knots_area1 = find_knots(map_area1, nd = 50)
time_knots1 <- (proc.time()-mtime_knots1)[3]


mtime_knots2 <- proc.time()
knots_area2 = find_knots(map_area2, nd = 200)
time_knots2 <- (proc.time()-mtime_knots2)[3]


mtime_knots3 <- proc.time()
knots_area3 = find_knots(map_area3, nd = 350)
time_knots3 <- (proc.time()-mtime_knots3)[3]


knots_list = list(knots_area1, knots_area2, knots_area3)

  for (area_i in 1:3){
    
    # Find knots with space filling algorithm
    knots = knots_list[[area_i]]
    
    # Graph of neighbours (for BYM)
    nb <- poly2nb(map_area_list[[area_i]], queen=T)
    adj <- unlist(nb)
    nb2INLA("adj",nb)
    g <- inla.read.graph(filename="adj")
    
    
    # Save the results
    
    rho_raster <- sigma_raster <- sigma_error_raster <- time_raster <- numeric(N)
    cov_raster <- cor_raster <- matrix(0, nrow = N, ncol = 500)  
    CI_spat_raster <- Bias_spat_raster <- RMSE_spat_raster <- matrix(0, nrow = N, ncol = prod(dim(raster_true)))  
    CI_continuous_incidence <- Bias_continuous_incidence <- RMSE_continuous_incidence <- matrix(0, nrow = N, ncol = prod(dim(raster_true)))  
    CI_area <- Bias_area <- RMSE_area <- matrix(0, nrow = N, ncol = dim(map_area_list[[area_i]])[1])  
    
    rho_disaggr <- sigma_disaggr <- time_disaggr <- numeric(N)
    cov_disaggr <- cor_disaggr <- matrix(0, nrow = N, ncol = 500)  
    CI_spat_disaggr <- Bias_spat_disaggr <- RMSE_spat_disaggr <- matrix(0, nrow = N, ncol = prod(dim(raster_true)))  
    CI_area_disaggr <- Bias_area_disaggr <- RMSE_area_disaggr <- matrix(0, nrow = N, ncol = dim(map_area_list[[area_i]])[1])  
    CI_continuous_incidence_disaggr <- Bias_continuous_incidence_disaggr <- RMSE_continuous_incidence_disaggr <- matrix(0, nrow = N, ncol = prod(dim(raster_true)))  
    
    rho_exact <- sigma_exact <- sigma_error_exact <- time_exact <- numeric(N)
    cov_exact <- cor_exact <- matrix(0, nrow = N, ncol = 500)  
    CI_spat_exact <- Bias_spat_exact <- RMSE_spat_exact <- matrix(0, nrow = N, ncol = prod(dim(raster_true)))  
    CI_continuous_incidence_exact <- Bias_continuous_incidence_exact <- RMSE_continuous_incidence_exact <- matrix(0, nrow = N, ncol = prod(dim(raster_true)))  
    CI_area_exact <- Bias_area_exact <- RMSE_area_exact <- matrix(0, nrow = N, ncol = dim(map_area_list[[area_i]])[1])  
    
    CI_continuous_incidence_atp <- Bias_continuous_incidence_atp <- RMSE_continuous_incidence_atp <-  matrix(0, nrow = N, ncol = prod(dim(raster_true)))  
    rho_atp <- sigma_atp <- time_atp <- numeric(N)
    cov_atp <- cor_atp <- matrix(0, nrow = N, ncol = 500)  
    
    
    rho_sdalgcp <- sigma_sdalgcp <- time_sdalgcp <- numeric(N)
    cov_sdalgcp <- cor_sdalgcp <- matrix(0, nrow = N, ncol = 500)  
    CI_spat_sdalgcp <- Bias_spat_sdalgcp <- RMSE_spat_sdalgcp <- matrix(0, nrow = N, ncol = prod(dim(raster_true)))  
    CI_area_sdalgcp <- Bias_area_sdalgcp <- RMSE_area_sdalgcp <- matrix(0, nrow = N, ncol = dim(map_area_list[[area_i]])[1])  
    
    
    CI_area_BYM <- Bias_area_BYM <- RMSE_area_BYM <- matrix(0, nrow = N, ncol = dim(map_area_list[[area_i]])[1])  
    time_BYM <- numeric(N)
    
    for (ii in 1:N){
      tryCatch({
      
      set.seed(ii)
      
      # Linear covariate raster
      rstr_linear = raster_true
      values(rstr_linear) = linear_cov(rstr_centers[,1],rstr_centers[,2]) + rnorm(prod(dim(rstr_linear)),0,0.1)
      
      # Spatial effect raster   
      rstr_true_spat = raster_true
      values(rstr_true_spat) = scale(sim_spat$data[,ii], scale = F)
      
      # True incidence raster
      rstr_true_incidence = raster_true
      values(rstr_true_incidence) = smooth_spat_poly_weighted(values(rstr_linear),
                                                              values(rstr_true_spat))
      
      # Simulate the true number of cases in every grid cell
      rstr_true_cases = rstr_true_incidence
      values(rstr_true_cases) = rpois(n = length(values(rstr_true_cases)), 
                                      lambda = values(rstr_true_incidence))
      
      # Aggregate mean incidence to area level
      true_inc <- exact_extract(rstr_true_incidence, map_area_list[[area_i]], 'sum', progress = F)
      
      # Aggregate the response data to area level
      y <- exact_extract(rstr_true_cases, map_area_list[[area_i]], 'sum', progress = F)
      
      # True population size at area level
      true_pop <- exact_extract(pop_den,map_area_list[[area_i]], 'sum', progress = F)
      
      # Linear covariate at area level (for BYM model) 
      xlin.0 = exact_extract(rstr_linear,map_area_list[[area_i]], 'mean', progress = F)
      
      # Initial values for SSDAM and SSDEM algorithm
      model <- glm(y ~xlin.0 + offset(log(true_pop)), family = "poisson")
      start_sigma = c(1/exp(1),min(mean(model$residuals^2),1))                                 
                                        
       ##############################   SSDAM  ###############################

      # Preparing the data (which cells intersect with which area)
      intersecting_cells_raster <- find_intersecting_cells(map = map_area_list[[area_i]],
                                                           rstr_covariates = list(linear = rstr_linear),
                                                           covariate_missing = F)
   

      dat_raster <- prepare_data(y = y, map = map_area_list[[area_i]], 
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
                                    limits = limits, map =  st_cast(st_union(map_area1), "POLYGON"))
      
      
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
                                                     model_raster, plot = F,
                                                     rstr_covariates = 
                                                       list(linear = rstr_linear))
      
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
      
      discrete_raster = predict_discrete_mean(map_area_list[[area_i]], model_raster, plot = F)
      
      ci_y0 <- numeric(length(true_inc))
       for (i in 1:length(true_inc)) {
        ci_y0[i] <- as.numeric(ifelse(log(true_inc)[i] >= discrete_raster$lower[i] & log(true_inc)[i] <= discrete_raster$upper[i], 1, 0))
       }
      CI_area[ii,] <- ci_y0
      Bias_area[ii,] <- true_inc - exp(discrete_raster$fit)
      RMSE_area[ii,] <- (true_inc-exp(discrete_raster$fit))^2
      
      
        
      ##############################   SSDEM  ###############################

      dat_exact <- prepare_data(y = y, map = map_area_list[[area_i]], 
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
                                    limits = limits, map =  st_cast(st_union(map_area1), "POLYGON"))
      
      
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
                                                     model_exact, plot = F,
                                                     rstr_covariates = 
                                                       list(linear = rstr_linear))
      
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
      
      discrete_exact = predict_discrete_mean(map_area_list[[area_i]], model_exact, plot = F)
      
      ci_y0 <- numeric(length(true_inc))
      for (i in 1:length(true_inc)) {
        ci_y0[i] <- as.numeric(ifelse(log(true_inc)[i] >= discrete_exact$lower[i] & log(true_inc)[i] <= discrete_exact$upper[i], 1, 0))
      }
      CI_area_exact[ii,] <- ci_y0
      Bias_area_exact[ii,] <- true_inc - exp(discrete_exact$fit)
      RMSE_area_exact[ii,] <- (true_inc - exp(discrete_exact$fit))^2
      
      
      ########################### Disaggregation package #########################
      map_disaggr = map_area_list[[area_i]]
      map_disaggr$y = y
      map_disaggr$ID = 1:dim(map_disaggr)[1]
      polygon_shapefile = map_disaggr[,c("ID", "y")]

      covariate_raster_full = c(terra::rast(rstr_linear))
      names(covariate_raster_full) = c("xlin")
      
   
      mtime_disaggr <- proc.time()
      dis_data = disaggregation::prepare_data(polygon_shapefile = polygon_shapefile,
                                              aggregation_raster = terra::rast(pop_den),
                                              covariate_rasters = covariate_raster_full,
                                              id_var = "ID", response_var = "y", na_action = T)
      
      
      
      
      fitted_disaggregation = disaggregation::disag_model(data = dis_data, iterations = 1000, family = "poisson", link = "log",
                                                           hess_control_ndeps = 1e-3)
      
      time_disaggr[ii] <- (proc.time()-mtime_disaggr)[3]
      
      
      range_disaggr = exp(fitted_disaggregation$opt$par["log_rho"])
      sigma_disaggr[ii] = exp(fitted_disaggregation$opt$par["log_sigma"])^2
      nu_disaggr = 1
      rho_disaggr[ii] = sqrt(8)/range_disaggr #similar to rho_approx
      
      cov_disaggr[ii,] = sigma_disaggr[ii] * rho_disaggr[ii] * spat_range * besselK(rho_disaggr[ii] * spat_range, nu_disaggr)
      cor_disaggr[ii,] = rho_disaggr[ii] * spat_range * besselK(rho_disaggr[ii] * spat_range, nu_disaggr)
      
      
      ######## Predict

      # Spatial effect
      spatial_prediction = disaggr_predict_spatial(fitted_disaggregation, 
                                                   xyFromCell(rstr_true_spat, 1:ncell(rstr_true_spat)))
      
      mean_spat = mean(spatial_prediction)
      spatial_prediction = spatial_prediction - mean_spat
      
      
      plot_estimated_disagg<- ggplot() +
        geom_raster(data = data.frame(x=xyFromCell(rstr_true_spat, 1:ncell(rstr_true_spat))[,1],y = xyFromCell(rstr_true_spat, 1:ncell(rstr_true_spat))[,2]),
                    aes(x = x, y = y, fill = spatial_prediction)) +
        geom_sf(data = st_cast(st_union(map_area1), "POLYGON"), fill = NA, color = "black") +
        scale_fill_viridis_c(limits = limits) + 
        theme_minimal() +
        theme(panel.background = element_rect(fill = "transparent", color = NA),
              plot.background = element_rect(fill = "transparent", color = NA)) +
        ggtitle("Estimated spatial effect disaggregation")+
        labs(fill = "Spatial effect")
      
      spatial_prediction_CI = disaggr_predict_spatial_uncertainty(fitted_disaggregation, xyFromCell(rstr_true_spat, 1:ncell(rstr_true_spat)), n = 1000)
      ci_spat_disaggr <- numeric(length(ysp.0))
      for (i in 1:length(ysp.0)) {
        ci_spat_disaggr[i] <- as.numeric(ifelse(ysp.0[i] >= spatial_prediction_CI[i, 1]- mean_spat & ysp.0[i] <= spatial_prediction_CI[i, 2]- mean_spat, 1, 0))
      }
      CI_spat_disaggr[ii,] <- ci_spat_disaggr
      Bias_spat_disaggr[ii,] <- ysp.0 - spatial_prediction
      RMSE_spat_disaggr[ii,] <- (ysp.0 - spatial_prediction)^2
      

      # Predict continuous incidence
      pred_continuous = disaggregation::predict_model(fitted_disaggregation, predict_iid = T)
      pred_continuous_CI = disaggregation::predict_uncertainty(fitted_disaggregation, predict_iid = T)
      pred_continuous_LL = pred_continuous_CI$predictions_ci[,,1]
      pred_continuous_UL = pred_continuous_CI$predictions_ci[,,2]
      
      ci_y0 <- numeric(length(true_effect))
      for (i in 1:length(true_effect)) {
        ci_y0[i] <- as.numeric(ifelse(true_effect[i] >= pred_continuous_LL[i] & 
                                        true_effect[i] <= pred_continuous_UL[i], 1, 0))
      }
      CI_continuous_incidence_disaggr[ii,] <- ci_y0
      Bias_continuous_incidence_disaggr[ii,] <- true_effect-values(pred_continuous$prediction)
      RMSE_continuous_incidence_disaggr[ii,] <- (true_effect-values(pred_continuous$prediction))^2
      
      
      ################## Area to point kriging ###################################
      
      map_atp = map_area_list[[area_i]]

      intersecting_cells_coordinates <- find_intersecting_cells(map = map_area_list[[area_i]],
                                                                rstr_covariates = list(linear = rstr_linear,
                                                                                       true_spat = rstr_true_spat,
                                                                                       true_inc = rstr_true_incidence),
                                                                covariate_missing = F)
      
      discretize.points = intersecting_cells_coordinates %>%
        rename(areaId = ID, ptx = x, pty = y, weight = tot_weight) %>%
        dplyr::select(areaId, ptx, pty, weight)

      xy = intersecting_cells_coordinates %>% group_by(ID) %>%
        summarize(weightedpop(as.vector(x),
                       as.vector(y),
                       as.vector(tot_weight)))
      
      data.test = data.frame(areaId = xy$ID, centx = xy$x, centy = xy$y,
                             counts = y, size = apply(matrix(true_pop), 1,int.pop)) 
      
      createDataset <- function(data.test, discretePoints){
        rslt <- list(areaValues = data.test, discretePoints = discretePoints)
        class(rslt) <- c("list", "discreteArea")
        return(rslt)
       }
      
      
      data_atp = createDataset(data.test, discretize.points)
      
      mtime_atp <- proc.time()
      point_variogram <- deconvPointVgm(data_atp, model="Exp", rd = 1, fig = T)
      
      mean_inc = mean(y/true_pop)
      sigma_atp[ii] <- log(point_variogram$pointVariogram$psill/mean_inc^2+1)
      range_atp <- point_variogram$pointVariogram$range
      
      rho_atp[ii] <- 1/range_atp
      
      cov_atp[ii,] <- sigma_atp[ii] * exp(-rho_atp[ii] * spat_range)
      cor_atp[ii,] <- exp(-rho_atp[ii] * spat_range)
      
      
      dat_pred <- data.frame( xyFromCell(rstr_true_spat, 1:ncell(rstr_true_spat)))
      
      nmax_all = c(2,5,10)
      predictions_atp <- atpKriging(data_atp, unknown = dat_pred, point_variogram,
                                    nmax = nmax_all[area_i], showProgress = T)
      
      time_atp[ii] <- (proc.time()-mtime_atp)[3]
      ci_y0 <- numeric(length(true_effect))
      for (i in 1:length(true_effect)) {
        ci_y0[i] <- as.numeric(ifelse(log(true_effect)[i] >= log(max(predictions_atp$pred[i], 1e-16))-
                                        qnorm(0.975)*sqrt(max(predictions_atp$var[i],0))/max(predictions_atp$pred[i], 1e-16) & 
                                        log(true_effect)[i] <= log(max(predictions_atp$pred[i], 1e-16))+
                                        qnorm(0.975)*sqrt(max(predictions_atp$var[i],0))/max(predictions_atp$pred[i], 1e-16), 1, 0))
      }
      CI_continuous_incidence_atp[ii,] <- ci_y0
      Bias_continuous_incidence_atp[ii,] <- true_effect-predictions_atp$pred
      RMSE_continuous_incidence_atp[ii,] <- (true_effect-predictions_atp$pred)^2
      
      
      ################################### SDALGCP ###################################
      
      xlin.0 = exact_extract(rstr_linear,map_area_list[[area_i]], 'mean', progress = F)
      
      data_lgcp <- data.frame(y = y, xlin = xlin.0, area_pop = true_pop)
      FORM <- y ~ xlin +  offset(log(area_pop))
      
      map_lgcp <-  as(map_area_list[[area_i]], "Spatial")
      
      # create the set of phi's to use (same as for SSDAM and SSDEM algorithm)
      r_ind = 2^seq(0,25*0.25, by=0.25)[6:25]
      phi <- max(model_raster$covs$Z.dist)/r_ind
      
      
      #set the mcmc parameters
      control.mcmc <- controlmcmcSDA(n.sim = 110000, burnin = 10000, thin= 10, 
                                     h=1.65/((dim(map_lgcp)[1])^(1/6)),
                                     c1.h = 0.01, c2.h = 1e-04)
      
      mtime_SDALGCP <- proc.time()
      
      my_est <- SDALGCPMCML(data=data_lgcp, formula=FORM, my_shp=map_lgcp,
                            delta=1, phi = phi, method=1, 
                            weighted=F,  plot=F, par0=NULL, control.mcmc=control.mcmc)
      
      
      range_sdalgcp = my_est$phi_opt
      sigma_sdalgcp[ii] = my_est$sigma2_opt
      rho_sdalgcp[ii] = 1/range_sdalgcp 
      cov_sdalgcp[ii,] = sigma_sdalgcp * exp(-rho_sdalgcp[ii]*spat_range)
      cor_sdalgcp[ii,] = exp(-rho_sdalgcp[ii]*spat_range)
      Con_pred <- SDALGCPPred(para_est=my_est,  cellsize=1,  continuous=TRUE, 
                              pred.loc = xyFromCell(rstr_true_spat, 1:ncell(rstr_true_spat)), 
                              control.mcmc = control.mcmc)
      
      time_sdalgcp[ii] <- (proc.time()-mtime_SDALGCP)[3]
      
      # Spatial level
      spatial_sdalgcp = log(as.numeric(Con_pred$pred))
      spatial_sdalgcp_CI = t(apply(Con_pred$pred.draw, 2, function(row) quantile(row, probs = c(0.025, 0.975))))
      
      mean_spat = mean(log(as.numeric(Con_pred$pred)))
      
      spatial_sdalgcp = spatial_sdalgcp - mean_spat
      
      ci_spat_sdalgcp <- numeric(length(ysp.0))
      for (i in 1:length(ysp.0)) {
        ci_spat_sdalgcp[i] <- as.numeric(ifelse(ysp.0[i] >= spatial_sdalgcp_CI[i, 1]- mean_spat &
                                                  ysp.0[i] <= spatial_sdalgcp_CI[i, 2]- mean_spat, 1, 0))
      }
      CI_spat_sdalgcp[ii,] <- ci_spat_sdalgcp
      Bias_spat_sdalgcp[ii,] <- ysp.0 - spatial_sdalgcp
      RMSE_spat_sdalgcp[ii,] <- (ysp.0 - spatial_sdalgcp)^2
      
 
      
      # Area level
      
      discrete_SDALGCP = Con_pred$my_shp$pMean_RR
      sd_SDALGCP = t(apply(Con_pred$S.draw, 2, function(row) quantile(row, probs = c(0.025, 0.975))))
      
      ci_y0 <- numeric(length(true_inc))
      for (i in 1:length(true_inc)) {
        ci_y0[i] <- as.numeric(ifelse(log(true_inc)[i] >= sd_SDALGCP[i, 1] + log(true_pop[i]) & 
                                        log(true_inc)[i] <= sd_SDALGCP[i, 2] + log(true_pop[i]) , 1, 0))
      }
      CI_area_sdalgcp[ii,] <- ci_y0
      Bias_area_sdalgcp[ii,] <- true_inc -discrete_SDALGCP*true_pop
      RMSE_area_sdalgcp[ii,] <- (true_inc-discrete_SDALGCP*true_pop)^2
      
      
      ################################### BYM ######################################

      region <- 1:dim(map_area_list[[area_i]])[1]
      formula <- y ~ xlin + f(region, model = "bym2", graph = g, param = c(2,0.5))
      my.data.inla = data.frame(y = y)
      my.data.inla$xlin = as.numeric(xlin.0)
      my.data.inla$area_pop = true_pop
      
      
      mtime_BYM <- proc.time()
      
      model.inla <- inla(formula, family = "Poisson", data = my.data.inla, offset = log(area_pop))
      
      time_BYM[ii] <- (proc.time()-mtime_BYM)[3]
      
      # Credible interval
      ci_spat_area_BYM <- numeric(length(true_inc))
      for (i in 1:length(true_inc)) {
        ci_spat_area_BYM[i] <- as.numeric(ifelse(true_inc[i] >= as.numeric(model.inla$summary.fitted.values[i,"0.025quant"]) & 
                                                   true_inc[i] <= as.numeric(model.inla$summary.fitted.values[i,"0.975quant"]), 1, 0))
      }
      CI_area_BYM[ii,] <- ci_spat_area_BYM
      Bias_area_BYM[ii,] <-true_inc - as.numeric(model.inla$summary.fitted.values$mean)
      RMSE_area_BYM[ii,] <- (true_inc - as.numeric(model.inla$summary.fitted.values$mean))^2
      
      },error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
      
      print(ii) 
      
    }
    
    ind <- !(time_BYM == 0) # All indices for which SDALGCP converged
    
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
                      mean(CI_area_exact[ind,])),2),
      disaggr = round(c(mean(time_disaggr[ind]),mean(Bias_spat_disaggr[ind,]),
                        sqrt(mean(RMSE_spat_disaggr[ind,])), mean(CI_spat_disaggr[ind,]),
                        mean(Bias_continuous_incidence_disaggr[ind,]),
                        sqrt(mean(RMSE_continuous_incidence_disaggr[ind,], na.rm=T)),
                        mean(CI_continuous_incidence_disaggr[ind,]), NA, NA, NA),2),
      SDALGCP = round(c(mean(time_sdalgcp[ind]), mean(Bias_spat_sdalgcp[ind,]),
                        sqrt(mean(RMSE_spat_sdalgcp[ind,], na.rm=T)),mean(CI_spat_sdalgcp[ind,]),
                        NA, NA, NA, mean(Bias_area_sdalgcp[ind,]), sqrt(mean(RMSE_area_sdalgcp[ind,], na.rm=T)),
                        mean(CI_area_sdalgcp[ind,])),2),
      atp = round(c(mean(time_atp[ind]),NA,NA,NA, mean(Bias_continuous_incidence_atp[Bias_continuous_incidence_atp<Inf & ind], na.rm=T),
                    sqrt(mean(RMSE_continuous_incidence_atp[RMSE_continuous_incidence_atp<Inf & ind], na.rm=T)),
                    mean(CI_continuous_incidence_atp[ind,], na.rm = T), NA, NA, NA),2),
      BYM = round(c(mean(time_BYM[ind]), NA, NA, NA, NA, NA, NA, mean(Bias_area_BYM[ind,]),
                    sqrt(mean(RMSE_area_BYM[ind,], na.rm = T)),
                    mean(CI_area_BYM[ind,])),2)
      
    )  
    
    result
    
    save.image(paste("Results/results_area",area_i,".RData", sep=""))

    print(area_i)
    
  }



# Compare areas 
load("Results/results_area1.RData")

data_compare_area3 = data.frame(RMSE_spat = as.numeric(result[3,-1]),
                                RMSE_cont_inc = as.numeric(result[6,-1]),
                                RMSE_area_inc = as.numeric(result[9,-1]),
                                Method = colnames(result)[-1],
                                area = "Area 3")

data_sigma_area3 = data.frame(sigma = c(sigma_raster[ind], sigma_exact[ind],
                                        sigma_disaggr[ind], sigma_sdalgcp[ind],
                                        sigma_atp[ind]),
                              Method = c(rep("SSDAM",sum(ind)), rep("SSDEM",sum(ind)), rep("disaggr",sum(ind)),
                                         rep("SDALGCP",sum(ind)), rep("atp",sum(ind))),
                              area = "Area 3")



data_cor_area3 = data.frame(cor = c(colMeans(cor_raster[ind,]), colMeans(cor_exact[ind,]),
                                    colMeans(cor_disaggr[ind,]), colMeans(cor_sdalgcp[ind,]),
                                    colMeans(cor_atp[ind,]),
                                    1/range * spat_range * besselK(1/range * spat_range, 1)),
                            Method = c(rep("SSDAM",length(spat_range)), 
                                       rep("SSDEM",length(spat_range)),
                                       rep("disaggr",length(spat_range)),
                                       rep("SDALGCP",length(spat_range)), 
                                       rep("atp",length(spat_range)),
                                       rep("true", length(spat_range))),
                            area = "Area 3",
                            x = rep(spat_range, 6))



result_area3 = result 

data_cor_area3$cor = ifelse(data_cor_area3$x==0,1,data_cor_area3$cor)


load("Results/results_area2.RData")

data_compare_area2 = data.frame(RMSE_spat = as.numeric(result[3,-1]),
                                RMSE_cont_inc = as.numeric(result[6,-1]),
                                RMSE_area_inc = as.numeric(result[9,-1]),
                                Method = colnames(result[-1]),
                                area = "Area 2")

data_sigma_area2 = data.frame(sigma = c(sigma_raster[ind], sigma_exact[ind],
                                        sigma_disaggr[ind], sigma_sdalgcp[ind],
                                        sigma_atp[ind]),
                              Method = c(rep("SSDAM",sum(ind)), rep("SSDEM",sum(ind)), rep("disaggr",sum(ind)),
                                         rep("SDALGCP",sum(ind)), rep("atp",sum(ind))),
                              area = "Area 2")



data_cor_area2 = data.frame(cor = c(colMeans(cor_raster[ind,]), colMeans(cor_exact[ind,]),
                                    colMeans(cor_disaggr[ind,]), colMeans(cor_sdalgcp[ind,]),
                                    colMeans(cor_atp[ind,]),
                                    1/range * spat_range * besselK(1/range * spat_range, 1)),
                            Method = c(rep("SSDAM",length(spat_range)), 
                                       rep("SSDEM",length(spat_range)),
                                       rep("disaggr",length(spat_range)),
                                       rep("SDALGCP",length(spat_range)), 
                                       rep("atp",length(spat_range)),
                                       rep("true", length(spat_range))),
                            area = "Area 2",
                            x = rep(spat_range, 6))



result_area2 = result 

data_cor_area2$cor = ifelse(data_cor_area2$x==0,1,data_cor_area2$cor)


load("Results/results_area3.RData")

data_compare_area1 = data.frame(RMSE_spat = as.numeric(result[3,-1]),
                                RMSE_cont_inc = as.numeric(result[6,-1]),
                                RMSE_area_inc = as.numeric(result[9,-1]),
                                Method = colnames(result[-1]),
                                area = "Area 1")

data_sigma_area1 = data.frame(sigma = c(sigma_raster[ind], sigma_exact[ind],
                                        sigma_disaggr[ind], sigma_sdalgcp[ind],
                                        sigma_atp[ind]),
                              Method = c(rep("SSDAM",sum(ind)), rep("SSDEM",sum(ind)), rep("disaggr",sum(ind)),
                                         rep("SDALGCP",sum(ind)), rep("atp",sum(ind))),
                              area = "Area 1")



data_cor_area1 = data.frame(cor = c(colMeans(cor_raster[ind,]), colMeans(cor_exact[ind,]),
                                    colMeans(cor_disaggr[ind,]), colMeans(cor_sdalgcp[ind,]),
                                    colMeans(cor_atp[ind,]),
                                    1/range * spat_range * besselK(1/range * spat_range, 1)),
                            Method = c(rep("SSDAM",length(spat_range)), 
                                       rep("SSDEM",length(spat_range)),
                                       rep("disaggr",length(spat_range)),
                                       rep("SDALGCP",length(spat_range)), 
                                       rep("atp",length(spat_range)),
                                       rep("true", length(spat_range))),
                            area = "Area 1",
                            x = rep(spat_range, 6))



result_area1 = result 

data_cor_area1$cor = ifelse(data_cor_area1$x==0,1,data_cor_area1$cor)


data_compare_areas = rbind(data_compare_area1, data_compare_area2, data_compare_area3)
data_sigma_areas <- rbind(data_sigma_area1, data_sigma_area2, data_sigma_area3)
data_cor_areas <- rbind(data_cor_area1, data_cor_area2, data_cor_area3)


# Plot 

pdf("Figures/RMSE_spat_simulation.pdf",height=4,width=6)
ggplot(data_compare_areas, aes(x = factor(area), y = RMSE_spat, group = Method, color = Method))+
  geom_line(position = position_dodge(width = 0.05))+
  geom_point(position = position_dodge(width = 0.05), size = 3,
             aes(shape = Method)) +
  labs(x = "Area", y = "RMSE spatial")+
  theme_minimal() +
  ylim(c(0,0.9))
dev.off()

pdf("Figures/RMSE_cont_simulation.pdf",height=4,width=6)
ggplot(data_compare_areas, aes(x = factor(area), y = RMSE_cont_inc, group = Method, color = Method,
                               linetype = Method))+
  geom_line(position = position_dodge(width = 0.05))+
  geom_point(position = position_dodge(width = 0.05), size = 3,
             aes(shape = Method)) +
  labs(x = "Area", y = "RMSE continuous incidence")+
  theme_minimal() +
  ylim(c(0,0.035))
dev.off()

pdf("Figures/RMSE_area_simulation.pdf",height=4,width=6)
ggplot(data_compare_areas, aes(x = factor(area), y = RMSE_area_inc, group = Method, color = Method))+
  geom_line(position = position_dodge(width = 0.05))+
  geom_point(position = position_dodge(width = 0.05), size = 3,
             aes(shape = Method)) +
  labs(x = "Area", y = "RMSE area counts")+
  theme_minimal() +
  ylim(c(0,4))
dev.off()


pdf("Figures/sigma.pdf",height=4,width=6)
ggplot(data_sigma_areas, aes(x = factor(area), y = sigma, fill = Method)) +
  geom_boxplot() +
  labs(x = "area", y = "variance") + theme_bw() +
  geom_hline(aes(yintercept = 0.7), linetype = 2) 
dev.off()

pdf("Figures/cor_mean.pdf",height=4,width=6)
ggplot(data_cor_areas, aes(x = x, y = cor, group = Method, color = Method,
                           linetype = Method))+
  geom_line()+
  facet_wrap(~area) +
  theme_bw() +
  ylab("Mean correlation")+
  theme(legend.position = "bottom")+
  scale_linetype_manual(values = c(6,5,4,3,2,1))
dev.off()
  





