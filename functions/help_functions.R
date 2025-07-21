################################################################################
#       Help functions to prepare data and summarize results                   #
################################################################################
# Different functions:
#   - find_knots: Calculate the knots with space filling algorithm
#   - find_intersecting_cells: Prepare input data: which grid cells intersect with which areas?
#   - prepare_data: Prepare data object that can be used as input to Krig.fit
#   - predict_spatial: Predict disaggregated spatial effect
#   - predict_smooth: Predict smooth effect
#   - predict_smooth_RR: Predict smooth effect RR compared to reference value
#   - predict_incidence_continuous: Predict disaggregated incidence
#   - predict_discret_mean: Predict area level mean count
#   - predict_excprob_continuous: Calculate exceedance probabilities at disaggregated scale


################################################################################
#   function to calculate knots
#   Input:
#       - map: sf map object
#       - nd: requested number of knots (default = 350)
#       - seed: Seed for random knot generation

find_knots <- function(map, seed = 150150, nd = 350){
  set.seed(seed)
  samples_map <- sf::st_sample(map, 5000)
  samples_matrix <- do.call(rbind, lapply(samples_map, sf::st_coordinates))
  knot.sp <- fields::cover.design(R = cbind(samples_matrix[,1], samples_matrix[,2]), 
                                  nd = nd)
  
  knots <- as.matrix(knot.sp$design)
  return(knots)
}

################################################################################
#   function to find intersecting cells
#   Input:
#       - rstr: Disaggregation raster object
#               If not provided, disaggregation raster is extracted from rstr_pop_density or rstr_covariates
#       - map: sf map object
#       - rstr_pop_density: Population density at disaggregation raster scale
#                           Assumed to have the same geometry as rstr and rstr_covariates, if provided
#       - rstr_covariates: List of covariate rasters
#                          Assumed to have the same geometry as rstr and rstr_pop_density, if provided
#                          List elements that represent linear covariates should contain "linear" in the name
#                          List elements that represent smooth covariates should contain "smooth" in the name
#       - covariate_missing: True if there are missing covariates that need to be estimated using IDW


find_intersecting_cells <- function(rstr = NULL, map, rstr_pop_density = NULL,
                                    rstr_covariates = NULL, covariate_missing = F){
  library(exactextractr)
  
  if(is.null(rstr) & is.null(rstr_pop_density) & is.null(rstr_covariates)){
    stop("Please provide a raster")
  }else if(!is.null(rstr_pop_density)){
    if (!is.null(rstr) | !is.null(rstr_covariates)){
      print("Warning: You provided multiple rasters: assuming all geometries are the same")
    }
    rstr_one = rstr_pop_density
    rstr_one@data@values = rep(1, length(rstr_one@data@values))
  }else if(!is.null(rstr_covariates) ){
    if (!is.null(rstr) | length(rstr_covariates)>1){
      print("Warning: You provided multiple rasters: assuming all geometries are the same")
    }
    rstr_one = rstr_covariates[[1]]
    rstr_one@data@values = rep(1, length(rstr_one@data@values))  
  }else{
    rstr_one = rasterize(as(map,"Spatial"), rstr, field = 1, 
                         fun = mean, background = NA)
    }
  
  intersecting_cells <- exact_extract(rstr_one,map, include_xy = T, progress = F)
  
  intersecting_cells_coordinates <- do.call(rbind, lapply(1:length(intersecting_cells), 
                                                          function(i) {
                                                            cbind(intersecting_cells[[i]],ID = rep(i, dim(intersecting_cells[[i]])[1]))
                                                          }))
  
  intersecting_cells_coordinates = intersecting_cells_coordinates[,-1]
  names(intersecting_cells_coordinates) = c("x","y","weight_intr","ID")
  
  intersecting_cells_coordinates$pop = rep(NA, dim(intersecting_cells_coordinates)[1])
  
  if(!is.null(rstr_pop_density)){
    true_pop_raster <- exact_extract(rstr_pop_density,map, progress = F) 
    intersecting_cells_coordinates$pop <- unlist(lapply(true_pop_raster, function(df) df[, 1]))
    
    intersecting_cells_coordinates$pop = ifelse(is.na(intersecting_cells_coordinates$pop),0,
                                                intersecting_cells_coordinates$pop)
  }
  
  if(! is.null(rstr_pop_density)){
  intersecting_cells_coordinates$tot_weight = intersecting_cells_coordinates$weight_intr*
    intersecting_cells_coordinates$pop+1e-16
  }else{
    intersecting_cells_coordinates$pop = 1
    intersecting_cells_coordinates$tot_weight = intersecting_cells_coordinates$weight_intr*
      intersecting_cells_coordinates$pop

  }
  
  if(!is.null(rstr_covariates)){
    name_cov = rep("",length(rstr_covariates))
    for (cov_i in 1:length(rstr_covariates)){
      true_raster <- exact_extract(rstr_covariates[[cov_i]],map, progress = F) 
      intersecting_cells_coordinates <- cbind(intersecting_cells_coordinates,
                                              unlist(lapply(true_raster, function(df) df[, 1])))
      name_cov[cov_i] = paste("linear_", cov_i)
    }
    if(!is.null(names(rstr_covariates))){
      names(intersecting_cells_coordinates)[((dim(intersecting_cells_coordinates)[2]-length(rstr_covariates)+1):dim(intersecting_cells_coordinates)[2])] = names(rstr_covariates)
      
    }else{
      names(intersecting_cells_coordinates)[((dim(intersecting_cells_coordinates)[2]-length(rstr_covariates)+1):dim(intersecting_cells_coordinates)[2])] = name_cov
    }
  }

  
  if(covariate_missing == T){
    for (cov_i in 1:length(rstr_covariates)){
      name = names(rstr_covariates)[cov_i]
      known_values <- intersecting_cells_coordinates[!is.na(intersecting_cells_coordinates[name]),]
      missing_values <- intersecting_cells_coordinates[is.na(intersecting_cells_coordinates[name]),]
      if (!is.null(missing_values)){
      coordinates(known_values) = ~x+y
      coordinates(missing_values) = ~x+y
      formula_str = paste(name, " ~ 1")
      idw_result<- gstat::idw(formula = as.formula(formula_str), known_values, missing_values)
      intersecting_cells_coordinates[is.na(intersecting_cells_coordinates[name]),name] <- idw_result$var1.pred
      }
    }
  }
  
  
  return(intersecting_cells_coordinates)
  
}


################################################################################
#   function to prepare data
#   Input:
#       - y: vector of response variable
#       - map: sf map object
#       - intersecting_cells: Output of find_intersecting_cells
#       - exact: True if you are preparing your data to use the exact method
#       - weighted: True if you are preparing your data to use the weighted version of the approximate method

prepare_data <- function(y, map, intersecting_cells, exact = F, weighted = F){
  if((exact == F & weighted == F) | exact == T){
    intersecting_cells = intersecting_cells %>% na.omit()
  return (list(y = y, map = map, rstr_centers = intersecting_cells, 
               weight = intersecting_cells$weight_intr, weighted = weighted))
  }else if (exact == F){
    intersecting_cells = intersecting_cells %>% na.omit()
    return (list(y = y, map = map, rstr_centers = intersecting_cells, 
                 weight = intersecting_cells$tot_weight, weighted = weighted))
  }
}

################################################################################
#   function to predict the disaggregated log spatial effect
#   Input:
#       - newdata: Matrix of coordinates at which predictions need to be made (e.g. center of grid cells)
#       - model: output model (result of Krig.fit)
#       - plot: True to provide a plot of the disaggregated surface
#       - limits: Vector with lower and upper limit of the log spatial effect in the plot
#       - map: sf map object, should only be provided to construct a plot


predict_spatial <- function(newdata, model, plot = T, limits = c(-2,2), map = NULL){
  xi_estim<- model$xi_estim
  Sigma <- model$Sigma
  covs <- model$covs
  v_mode <- model$v_mode
  K_spat <- model$covs$K

  smooth_cov = sum(grepl("smooth",names(model$data_points)))>0
  
  X_sd = model$covs$sd_X
  
  w1.0 <- newdata[,1]
  w2.0 <- newdata[,2]
  
  Z1.dist.0 <- outer(w1.0,covs$knots.s[,1],"-")
  Z2.dist.0 <- outer(w2.0,covs$knots.s[,2],"-")
  Z.dist.0 <- sqrt(Z1.dist.0^2 + Z2.dist.0^2)

  
  Zw.0 <- covs$Scov(Z.dist.0, v = v_mode)

  
  constspat <- c(covs$constX, covs$constZ(v_mode))
  Cspat.cent <- as.matrix(cbind(w1.0, w2.0, Zw.0) - matrix(rep(constspat, each = length(w1.0), byrow=T), nrow=length(w1.0)))
  Cspat.cent[,1] = Cspat.cent[,1]/X_sd[1]
  Cspat.cent[,2] = Cspat.cent[,2]/X_sd[2]
 
  xi_spat <- xi_estim[(length(xi_estim)-(K_spat +2+length(unique(model$data_points$ID)))+1):(length(xi_estim)-length(unique(model$data_points$ID)))]
  fit.w0 <-  as.numeric(Cspat.cent%*%xi_spat)

  fit.w0 <- fit.w0 - mean(fit.w0)
  pred_grid = data.frame(X = newdata[,1],
                         Y = newdata[,2],
                         fit.w0 = fit.w0)

  if (plot & is.null(map)){stop("Please provide a map for plotting")}
  if (plot){
    p <- ggplot() +
      geom_raster(data = as.data.frame(pred_grid), aes(x = X, y = Y, fill = fit.w0)) +
      geom_sf(data = map, fill = NA, color = "black") +
      scale_fill_viridis_c(limits = limits) + 
      theme_minimal() +
      theme(panel.background = element_rect(fill = "transparent", color = NA),
            plot.background = element_rect(fill = "transparent", color = NA)) +
      ggtitle("Estimated spatial effect")+
      labs(fill = "Spatial effect")

    print(p)
  }
  # Coverage
  
  if(smooth_cov){
    all_smooth_cov = sum(grepl("smooth",names(model$data_points)))
    sd_ID = 1/sqrt(exp(v_mode[all_smooth_cov+1]))
  }else{
    sd_ID = 1/sqrt(exp(v_mode[1]))
  }
  
  
  ind.spat <- (length(xi_estim)-(K_spat +2+length(unique(model$data_points$ID)))+1):(length(xi_estim)-length(unique(model$data_points$ID)))
  covarxi_spat <- Sigma[ind.spat, ind.spat]
  sd.spat <- sqrt(pmax(0,Matrix::rowSums(((Cspat.cent%*%covarxi_spat)*Cspat.cent))))
  
  set.seed(1)
  quantiles_spat <- mapply(function(mean_val, sd_val) {
    rsim = rnorm(n = 100, mean = mean_val, sd = sd_val)+rnorm(100,0,sd_ID)
    return (as.numeric(quantile(rsim, prob = c(0.025,0.975), na.rm=T)))
  }, fit.w0, sd.spat)
  quantiles_spat <- t(quantiles_spat)
  
  
  result = data.frame(fit = fit.w0, lower = quantiles_spat[,1], upper = quantiles_spat[,2])

  return(result)
  
}


################################################################################
#   function to predict the smooth covariate(s)
#   Input:
#       - xgrid: Grid of x values for which to predict the smooth effect
#       - model: output model (result of Krig.fit)
#       - var: Variable for which to predict the smooth effect (same name as provided in raster_covariates)


predict_smooth <- function(xgrid, model, var){
  if(is.null(model$sm)){
    stop("You did not specify any smooth effect")
  }
  xi_estim<- model$xi_estim
  Sigma <- model$Sigma
  covs <- model$covs
  v_mode <- model$v_mode
  K_spat <- model$covs$K
  nbasis<- dim(model$sm[[1]])[2]
  
  all_smooth_cov = names(model$data_points)[grepl("smooth",names(model$data_points))]
  ind_var = which(var == all_smooth_cov)
  if(length(ind_var)==0){stop("Variable not found in smooth effects")}
  
  xsmooth <- Matrix::Matrix(model$sm[[ind_var]])
  knots_sm <- Matrix::Matrix(attr(model$sm[[ind_var]], "knots"))
  X_sd = model$covs$sd_X
  
  ind =  (length(xi_estim)-(K_spat +2+length(unique(model$data_points$ID)))-nbasis*(length(all_smooth_cov)-ind_var+1)+1):(length(xi_estim)-(K_spat +2+length(unique(model$data_points$ID)))-nbasis*(length(all_smooth_cov)-ind_var))
  xi_smooth = xi_estim[ind]  
  Zsm.grid <- Matrix::Matrix(ps(xgrid, df = nbasis, knots = knots_sm))
  
  
  fshat.xgrid <- as.numeric(Zsm.grid %*% xi_smooth)
  fshat.xgrid <- fshat.xgrid - mean(fshat.xgrid)
  
  covarxi_smooth <- Sigma[ind, ind]
  sd <- sqrt(pmax(0,Matrix::rowSums(((Zsm.grid%*%covarxi_smooth)*Zsm.grid))))
  quantiles_xgrid <- mapply(function(mean_val, sd_val) {
    qnorm(p = c(0.025, 0.975), mean = mean_val, sd = sd_val)
  }, fshat.xgrid, sd)
  quantile_xsmooth <- t(quantiles_xgrid)

  result = data.frame(fit = fshat.xgrid, lower = quantile_xsmooth[,1], upper = quantile_xsmooth[,2])
  
  
}


################################################################################
#   function to predict the RR of the smooth covariate(s)
#   Input:
#       - xgrid: Grid of x values for which to predict the smooth effect
#       - model: output model (result of Krig.fit)
#       - cen: Reference/centering value (for which RR is 1)
#       - var: Variable for which to predict the smooth effect (same name as provided in raster_covariates)

predict_smooth_RR <- function(xgrid, model, cen = NULL, var){
  if(is.null(model$sm)){
    stop("You did not specify any smooth effect")
  }
  xi_estim<- model$xi_estim
  Sigma <- model$Sigma
  covs <- model$covs
  v_mode <- model$v_mode
  K_spat <- model$covs$K
  nbasis<- dim(model$sm[[1]])[2]
  
  all_smooth_cov = names(model$data_points)[grepl("smooth",names(model$data_points))]
  ind_var = which(var == all_smooth_cov)
  if(length(ind_var)==0){stop("Variable not found in smooth effects")}
  
  xsmooth <- Matrix::Matrix(model$sm[[ind_var]])
  knots_sm <- Matrix::Matrix(attr(model$sm[[ind_var]], "knots"))
  
  X_sd = model$covs$sd_X
  
  ind =  (length(xi_estim)-(K_spat +2+length(unique(model$data_points$ID)))-nbasis*(length(all_smooth_cov)-ind_var+1)+1):(length(xi_estim)-(K_spat +2+length(unique(model$data_points$ID)))-nbasis*(length(all_smooth_cov)-ind_var))
  xi_smooth = xi_estim[ind]  
  Zsm.grid <- Matrix::Matrix(ps(xgrid, df = nbasis, knots = knots_sm))
  
  if(is.null(cen)){
    cen = min(xgrid)
  }
  cenvec <- Matrix::Matrix(ps(cen, df = nbasis, knots = knots_sm))
  Zsm.grid <- scale(Zsm.grid, center=cenvec, scale=F)
  
  fshat.xgrid <- as.numeric(Zsm.grid %*% xi_smooth)
  
  covarxi_smooth <- Sigma[ind, ind]
  sd <- sqrt(pmax(0,Matrix::rowSums(((Zsm.grid%*%covarxi_smooth)*Zsm.grid))))
  quantiles_xgrid <- mapply(function(mean_val, sd_val) {
    qnorm(p = c(0.025, 0.975), mean = mean_val, sd = sd_val)
  }, fshat.xgrid, sd)
  quantile_xsmooth <- t(quantiles_xgrid)
  
  result = data.frame(fit = exp(fshat.xgrid), lower = exp(quantile_xsmooth[,1]), upper = exp(quantile_xsmooth[,2]))
  
  
}



################################################################################
#   function to predict the continuous incidence 
#   Input:
#       - newdata: grid of coordinates at which to predict (data frame with two columns)
#       - model: output model (result of Krig.fit)
#       - plot: True to plot the continuous incidence
#       - limits: Vector with lower and upper limit of the continuous incidence in the plot
#       - map: sf map object, should only be provided to construct a plot
#       - rstr_covariates: List of covariate rasters

predict_incidence_continuous <- function(newdata, model, plot = T, limits = c(-2,2), map = NULL,
                                         rstr_covariates = NULL){
  
  # Extract the covariate values in the newdata coordinates
  if(!is.null(rstr_covariates)){
    name_cov = rep("",length(rstr_covariates))
    for (cov_i in 1:length(rstr_covariates)){
      true_raster <- raster::extract(rstr_covariates[[cov_i]],newdata[,1:2], na.rm=T) 
      newdata <- cbind(newdata,true_raster)
      name_cov[cov_i] = paste("linear_", cov_i)
    }
    if(!is.null(names(rstr_covariates))){
      names(newdata)[((dim(newdata)[2]-length(rstr_covariates)+1):dim(newdata)[2])] = names(rstr_covariates)
      
    }else{
      names(newdata)[((dim(newdata)[2]-length(rstr_covariates)+1):dim(newdata)[2])] = name_cov
    }
  }
  
  
  xi_estim<- model$xi_estim
  Sigma <- model$Sigma
  covs <- model$covs
  v_mode <- model$v_mode
  K_spat <- model$covs$K
  if(!is.null(model$sm)){
    nbasis<- dim(model$sm[[1]])[2]
    knots_spat = list()
    
    for (s in 1:length(model$sm)){
    knots_spat[[s]] <- Matrix::Matrix(attr(model$sm[[s]], "knots"))
    }
  }
  X_sd = model$covs$sd_X

  
  w1.0 <-  newdata[,1]
  w2.0 <-  newdata[,2]
  
  Z1.dist.0 <- outer(w1.0,covs$knots.s[,1],"-")
  Z2.dist.0 <- outer(w2.0,covs$knots.s[,2],"-")
  Z.dist.0 <- sqrt(Z1.dist.0^2 + Z2.dist.0^2)
  
  Zw.0 <- covs$Scov(Z.dist.0, v = v_mode)
  constspat <- c(covs$constX, covs$constZ(v_mode))
  Cspat.cent <- as.matrix(cbind(w1.0, w2.0, Zw.0) - matrix(rep(constspat, each = length(w1.0), byrow=T), nrow=length(w1.0)))
  Cspat.cent[,1] = Cspat.cent[,1]/X_sd[1]
  Cspat.cent[,2] = Cspat.cent[,2]/X_sd[2]
  
  smooth_cov = sum(grepl("smooth",names(newdata)))>0
  
  xlin = as.matrix(data.frame(newdata[,grep("linear", names(newdata))]))

  
  if(smooth_cov){
  all_smooth_cov = newdata[grepl("smooth",names(newdata))]
  Csmooth = list()
  for (s in 1:dim(all_smooth_cov)[2]){  
    Csmooth[[s]] =  Matrix::Matrix(ps(all_smooth_cov[,s], df = nbasis, knots = knots_spat[[s]]))
    
  }
    
  Csmooth <- do.call(cbind, Csmooth)
  C0 <- Matrix::Matrix(cbind(1,xlin, Csmooth, as.matrix(Cspat.cent)), sparse = T)
  sd_ID = 1/sqrt(exp(v_mode[dim(all_smooth_cov)[2]+1]))
  }else{
    sd_ID = 1/sqrt(exp(v_mode[1]))
    C0 <- Matrix::Matrix(cbind(1,xlin, as.matrix(Cspat.cent)), sparse = T)
  }

  xi <- xi_estim[1:(length(xi_estim)-length(unique(model$data_points$ID)))]
  
  # Compute mean and variance for log(mu)
  logmu <- exp(as.numeric(C0%*%xi))
  vcov = Sigma[1:(length(xi_estim)-length(unique(model$data_points$ID))),1:(length(xi_estim)-length(unique(model$data_points$ID)))]
  logmu.var <- sqrt(pmax(0,Matrix::rowSums((C0%*%vcov)*C0)))
  
  # Quantiles
  set.seed(1)
  quantiles_continuous <- mapply(function(mean_val, sd_val) {
    rsim = rnorm(n = 100, mean = mean_val, sd = sd_val)+rnorm(100,0,sd_ID)
    return (as.numeric(quantile(rsim, prob = c(0.025,0.975), na.rm=T)))
  }, log(logmu), logmu.var)
  quantiles_continuous <- t(quantiles_continuous)
  
 
  result = data.frame(fit = log(logmu), lower = quantiles_continuous[,1], upper = quantiles_continuous[,2])
  
  if (plot){
    plot_continuous_incidence <- ggplot() +
      geom_raster(data = data.frame(x=w1.0,y =w2.0), aes(x = x, y = y, fill = logmu)) +
      geom_sf(data = map, fill = NA, color = "black") +
      scale_fill_viridis_c(limits = limits) + 
      theme_minimal() +
      theme(panel.background = element_rect(fill = "transparent", color = NA),
            plot.background = element_rect(fill = "transparent", color = NA)) +
      ggtitle("Estimated log incidence continuous")+
      labs(fill = "Incidence")
    
    print(plot_continuous_incidence)

  }

  
  return(result)
  
}

################################################################################
#   function to predict the discrete mean (area-level expected COUNT, not incidence)
#   Input:
#       - map: sf map object
#       - model: output model (result of Krig.fit)
#       - plot: True to plot the area-level incidence
#       - limits: Vector with lower and upper limit of the area-level incidence in the plot
#       - pop: Area level population, should only be provided to plot the area-level incidence

predict_discrete_mean <- function(map, model, plot = T, limits = c(0,6), pop = NULL){
  xi_estim<- model$xi_estim
  Sigma <- model$Sigma
  covs <- model$covs
  v_mode <- model$v_mode
  K_spat <- model$covs$K
  
  if(!is.null(model$sm)){
    nbasis<- dim(model$sm[[1]])[2]
    knots_spat = list()
    
    for (s in 1:length(model$sm)){
      knots_spat[[s]] <- Matrix::Matrix(attr(model$sm[[s]], "knots"))
    }
  }

  
  X_sd = model$covs$sd_X
  
  w1.0 <-  model$data_points$x
  w2.0 <-  model$data_points$y
  
  
  Z1.dist.0 <- outer(w1.0,covs$knots.s[,1],"-")
  Z2.dist.0 <- outer(w2.0,covs$knots.s[,2],"-")
  Z.dist.0 <- sqrt(Z1.dist.0^2 + Z2.dist.0^2)
  
  Zw.0 <- covs$Scov(Z.dist.0, v = v_mode)
  constspat <- c(covs$constX, covs$constZ(v_mode))
  Cspat.cent <- as.matrix(cbind(w1.0, w2.0, Zw.0) - matrix(rep(constspat, each = length(w1.0), byrow=T), nrow=length(w1.0)))
  Cspat.cent[,1] = Cspat.cent[,1]/X_sd[1]
  Cspat.cent[,2] = Cspat.cent[,2]/X_sd[2]
  
  ID = 1:dim(map)[1]
  
  if(model$exact == F ){
    if(model$weighted == F){
    weight = model$data_points$weight_intr
    
    
    Weight_matrix_data <- data.frame(ID =  model$data_points$ID, weight = weight) %>%
      group_by(ID)%>%
      mutate(tot_weight = sum(weight))%>%
      ungroup()%>%
      mutate(full_weight = weight/tot_weight)
    
    
    weight_matrix <- Matrix::sparseMatrix(
      i = 1:length(model$data_points$ID),               
      j = model$data_points$ID,         
      x = Weight_matrix_data$full_weight,
      dims = c(length(model$data_points$ID), dim(map)[1])
    )
    }else{
      weight = model$data_points$tot_weight
      
      
      Weight_matrix_data <- data.frame(ID =  model$data_points$ID, weight = weight) %>%
        group_by(ID)%>%
        mutate(tot_weight = sum(weight))%>%
        ungroup()%>%
        mutate(full_weight = weight/tot_weight)
      
      
      weight_matrix <- Matrix::sparseMatrix(
        i = 1:length(model$data_points$ID),               
        j = model$data_points$ID,         
        x = Weight_matrix_data$full_weight,
        dims = c(length(model$data_points$ID), dim(map)[1])
      )
    }
    
    xlin_cont = as.matrix(data.frame(model$data_points[,grep("linear", names(model$data_points))]))
    X_linear <- cbind(rep(1,dim(map)[1]), Matrix::t(weight_matrix)%*%xlin_cont)
    discrete_C <- Matrix::t(weight_matrix)%*%Cspat.cent
    
    smooth_cov = sum(grepl("smooth",names(model$data_points)))>0
    
    if(smooth_cov){
      all_smooth_cov = data.frame(model$data_points[grepl("smooth",names(model$data_points))])
      Csmooth = list()
      
      for (s in 1:dim(all_smooth_cov)[2]){
      Csmooth[[s]] <- Matrix::Matrix(ps(all_smooth_cov[,s], df = nbasis, knots = knots_spat[[s]]))
      
      }
      Csmooth = do.call(cbind, Csmooth)
      discrete_Csmooth <- Matrix::t(weight_matrix)%*%Csmooth

      Z.rand <- Matrix::sparse.model.matrix(~ as.factor(ID) + 0)
      
      C0 = Matrix::Matrix(cbind(as.matrix(X_linear),discrete_Csmooth, discrete_C, Z.rand), sparse = T)
    }else{
      Z.rand <- Matrix::sparse.model.matrix(~ as.factor(ID) + 0)
      
      C0 = Matrix::Matrix(cbind(as.matrix(X_linear), discrete_C, Z.rand), sparse = T)
    }
    
    # Compute mean and variance for log(mu)
    logmu <- as.numeric(C0%*%xi_estim+log(model$offset))
    logmu.var <- sqrt(pmax(0,diag(as.matrix(C0%*%Sigma%*%t(C0)))))
    
    # Quantiles
    quantiles_spat_area <- mapply(function(mean_val, sd_val) {
      qnorm(p = c(0.025, 0.975), mean = mean_val, sd = sd_val)
    }, logmu, logmu.var)
    quantiles_w0_area <- t(quantiles_spat_area)
    
    result = data.frame(fit = logmu, lower = quantiles_w0_area[,1], upper = quantiles_w0_area[,2])
    
    
    
    if (plot){
      map$pred <- logmu - log(model$offset)
      
      map_plot <- tm_shape(map) +
        tm_polygons("pred", palette = "-RdYlGn",title="Estimated log incidence", 
                    breaks = seq(limits[1],limits[2],length=6)) +
        tm_layout(legend.outside = T,legend.outside.position="bottom") 
      
      print(map_plot)
    }
    
  }else{
    if (!require("dlnm", character.only = TRUE)) {
      message(paste("Package dlnm", "is required but not installed."))
    }

    weight = model$data_points$weight_intr
    
    Weight_matrix_data <- data.frame(ID =  model$data_points$ID, weight = weight) #%>%

    weight_matrix <- Matrix::sparseMatrix(
      i = 1:length(model$data_points$ID),               
      j = model$data_points$ID,         
      x = Weight_matrix_data$weight,
      dims = c(length(model$data_points$ID), dim(map)[1])
    )
    
    xlin_cont = as.matrix(data.frame(model$data_points[,grep("linear", names(model$data_points))]))
    X_linear <- cbind(rep(1,length(model$data_points$ID)), as.matrix(xlin_cont))
    
    
    smooth_cov = sum(grepl("smooth",names(model$data_points)))>0

    if(smooth_cov){
      all_smooth_cov = data.frame(model$data_points[grepl("smooth",names(model$data_points))])
      Csmooth = list()
      
      for (s in 1:dim(all_smooth_cov)[2]){
        Csmooth[[s]] <- Matrix::Matrix(as.matrix(ps(all_smooth_cov[,s], df = nbasis, knots = knots_spat[[s]])))
        
      }
      
    Csmooth <- do.call(cbind, Csmooth)
    Z.rand <- Matrix::sparse.model.matrix(~ as.factor(model$data_points$ID) + 0)
    C0 <- Matrix::Matrix(cbind(as.matrix(X_linear),
                               Csmooth, Cspat.cent, as.matrix(Z.rand)), sparse = T)
    }else{
      Z.rand <- Matrix::sparse.model.matrix(~ as.factor(model$data_points$ID) + 0)
      C0 <- Matrix::Matrix(cbind(as.matrix(X_linear),
                                 Cspat.cent, as.matrix(Z.rand)), sparse = T)
    }
    
    pop = model$data_points$pop
    
    logmu = as.numeric(C0%*%xi_estim)
    mu_area = as.numeric(Matrix::t(weight_matrix)%*%Matrix::Matrix(pop*exp(logmu)))
    
    # Quantiles
    set.seed(1)
    k <- length(xi_estim)
    eigen <- eigen(model$Sigma)
    
    eigen$values <- ifelse(eigen$values<0,0, eigen$values)
    X <- matrix(rnorm(length(xi_estim)*2000),2000)
    coefsim <- xi_estim + eigen$vectors %*% diag(sqrt(eigen$values),k) %*% t(X)
    
    quantiles_spat_area <-apply(coefsim,2, function(coefi) {
      rsim_logmu = as.numeric(C0%*%coefi)
      rsim_mu = as.numeric(Matrix::t(weight_matrix)%*%Matrix::Matrix(pop*exp(rsim_logmu)))
      return (as.numeric(log(rsim_mu)))
    })

            
    quantiles_spat_area_level = apply(quantiles_spat_area, 1, function(row) quantile(row, probs = c(0.025, 0.975)))
    quantiles_w0_area <- t(quantiles_spat_area_level)

    
    medians = apply(quantiles_spat_area, 1, function(row) quantile(row, probs = c(0.5)))
    result = data.frame(fit = medians, lower = quantiles_w0_area[,1], upper = quantiles_w0_area[,2])

  
    if (plot){
      if (!is.null(pop)){
      map$pred <- log(logmu_area/pop)
      map_plot <- tm_shape(map) +
        tm_polygons("pred", palette = "-RdYlGn",title="Estimated log incidence", 
                    breaks = seq(limits[1],limits[2],length=6)) +
        tm_layout(legend.outside = T,legend.outside.position="bottom") 
      
      }else{
        map$pred <- log(logmu_area)
        map_plot <- tm_shape(map) +
          tm_polygons("pred", palette = "-RdYlGn",title="Estimated log mean", 
                      breaks = seq(limits[1],limits[2],length=6)) +
          tm_layout(legend.outside = T,legend.outside.position="bottom") 
        
      }
      print(map_plot)
      }

    
    
  }


  
  return(result)
  
}


################################################################################
#   function to calculate the exceedance probability of the continuous incidence (grid level)
#   Input:
#       - newdata: grid of coordinates at which to predict (data frame with two columns)
#       - model: output model (result of Krig.fit)
#       - plot: True to plot the exceedance probabilities
#       - limits: Vector with lower and upper limit of the exceedance probability in the plot
#       - map: sf map object, should only be provided to construct a plot
#       - rstr_covariates: List of covariate rasters
#       - Threshold: Threshold for which to calculate the exceedance probability

predict_excprob_continuous <- function(newdata, model, plot = T, limits = c(0,1), map = NULL,
                                         rstr_covariates = NULL, threshold = 0.5){
  
  # Extract the covariate values in the newdata coordinates
  if(!is.null(rstr_covariates)){
    name_cov = rep("",length(rstr_covariates))
    for (cov_i in 1:length(rstr_covariates)){
      true_raster <- raster::extract(rstr_covariates[[cov_i]],newdata[,1:2], na.rm=T) 
      newdata <- cbind(newdata,true_raster)
      name_cov[cov_i] = paste("linear_", cov_i)
    }
    if(!is.null(names(rstr_covariates))){
      names(newdata)[((dim(newdata)[2]-length(rstr_covariates)+1):dim(newdata)[2])] = names(rstr_covariates)
      
    }else{
      names(newdata)[((dim(newdata)[2]-length(rstr_covariates)+1):dim(newdata)[2])] = name_cov
    }
  }
  
  
  xi_estim<- model$xi_estim
  Sigma <- model$Sigma
  covs <- model$covs
  v_mode <- model$v_mode
  K_spat <- model$covs$K
  if(!is.null(model$sm)){
    nbasis<- dim(model$sm[[1]])[2]
    knots_spat = list()
    
    for (s in 1:length(model$sm)){
      knots_spat[[s]] <- Matrix::Matrix(attr(model$sm[[s]], "knots"))
    }
  }
  X_sd = model$covs$sd_X
  
  
  w1.0 <-  newdata[,1]
  w2.0 <-  newdata[,2]
  
  Z1.dist.0 <- outer(w1.0,covs$knots.s[,1],"-")
  Z2.dist.0 <- outer(w2.0,covs$knots.s[,2],"-")
  Z.dist.0 <- sqrt(Z1.dist.0^2 + Z2.dist.0^2)
  
  Zw.0 <- covs$Scov(Z.dist.0, v = v_mode)
  constspat <- c(covs$constX, covs$constZ(v_mode))
  Cspat.cent <- as.matrix(cbind(w1.0, w2.0, Zw.0) - matrix(rep(constspat, each = length(w1.0), byrow=T), nrow=length(w1.0)))
  Cspat.cent[,1] = Cspat.cent[,1]/X_sd[1]
  Cspat.cent[,2] = Cspat.cent[,2]/X_sd[2]
  
  smooth_cov = sum(grepl("smooth",names(newdata)))>0
  
  xlin = as.matrix(data.frame(newdata[,grep("linear", names(newdata))]))
  
  
  if(smooth_cov){
    all_smooth_cov = newdata[grepl("smooth",names(newdata))]
    Csmooth = list()
    for (s in 1:dim(all_smooth_cov)[2]){  
      Csmooth[[s]] =  Matrix::Matrix(ps(all_smooth_cov[,s], df = nbasis, knots = knots_spat[[s]]))
      
    }
    
    Csmooth <- do.call(cbind, Csmooth)
    C0 <- Matrix::Matrix(cbind(1,xlin, Csmooth, as.matrix(Cspat.cent)), sparse = T)
    sd_ID = 1/sqrt(exp(v_mode[dim(all_smooth_cov)[2]+1]))
  }else{
    sd_ID = 1/sqrt(exp(v_mode[1]))
    C0 <- Matrix::Matrix(cbind(1,xlin, as.matrix(Cspat.cent)), sparse = T)
  }
  
  xi <- xi_estim[1:(length(xi_estim)-length(unique(model$data_points$ID)))]
  
  # Compute mean and variance for log(mu)
  logmu <- exp(as.numeric(C0%*%xi))
  vcov = Sigma[1:(length(xi_estim)-length(unique(model$data_points$ID))),1:(length(xi_estim)-length(unique(model$data_points$ID)))]
  logmu.var <- sqrt(pmax(0,Matrix::rowSums((C0%*%vcov)*C0)))
  
  # Quantiles
  set.seed(1)
  exceedance_prob <- mapply(function(mean_val, sd_val) {
    rsim = rnorm(n = 1000, mean = mean_val, sd = sd_val)+rnorm(1000,0,sd_ID)
    return (exp(as.numeric(rsim)))
  }, log(logmu), logmu.var)
  exceedance_prob <- t(exceedance_prob)
  
  
  exc_prob <- apply(exceedance_prob, 1, function(row){
    sum(row>threshold)/1000
  })
  
  
  result = exc_prob
  
  if (plot){
    plot_continuous_incidence <- ggplot() +
      geom_raster(data = data.frame(x=w1.0,y =w2.0), aes(x = x, y = y, fill = exc_prob)) +
      geom_sf(data = map, fill = NA, color = "black") +
      scale_fill_viridis_c(limits = limits) + 
      theme_minimal() +
      theme(panel.background = element_rect(fill = "transparent", color = NA),
            plot.background = element_rect(fill = "transparent", color = NA)) +
      ggtitle("Exceedance probability")+
      labs(fill = paste("Prob(inc>",round(threshold,2),")"), x = "", y = "")
    
    print(plot_continuous_incidence)
    
  }
  
  
  return(result)
  
}
