
disaggr_predict_spatial = function(output, newdata){
    Amatrix = fmesher::fm_evaluator(output$data$mesh, loc = as.matrix(newdata))$proj$A
    pars <- output$obj$env$last.par.best
    pars <- split(pars, names(pars))
    
    # Spatial objects
    field = (Amatrix %*% pars$nodemean)[,1]
    
    return(field)
}

disaggr_predict_spatial_uncertainty = function(output, newdata, n = 100){
  parameters <- output$obj$env$last.par.best
  
  if(output$model_setup$iid | output$model_setup$field){
    ch <- Matrix::Cholesky(output$sd_out$jointPrecision)
    par_draws <- sparseMVN::rmvn.sparse(n, parameters, ch, prec = TRUE)
  } else {
    covariance_matrix <- Matrix::Matrix(output$sd_out$cov.fixed, sparse = TRUE)
    ch <- Matrix::Cholesky(covariance_matrix)
    par_draws <- sparseMVN::rmvn.sparse(n, parameters, ch, prec = FALSE)
  }
  
  predictions <- list()
  Amatrix = fmesher::fm_evaluator(output$data$mesh, loc = as.matrix(newdata))$proj$A
  
  for(r in seq_len(n)) {
    
    p <- split(par_draws[r, ], names(parameters))
    
    predictions[[r]] <- (Amatrix %*% p$nodemean)[,1]
  }
  
  sim_matrix = do.call(cbind, predictions)
  
  predictions_ci <- apply(sim_matrix, 1, quantile, probs = c(0.025, 0.975))
  predictions_ci <- t(predictions_ci)
  
  
}



disaggr_pred_continuous_incidence = function(output, newdata, predict_iid = F){
  
  grid_area = data.frame(x = newdata$x, y = newdata$y)
  Amatrix = fmesher::fm_evaluator(output$data$mesh, loc = as.matrix(grid_area))$proj$A
  
  pars <- output$obj$env$last.par.best
  pars <- split(pars, names(pars))
  
  
  # linear covariates
  covariates = output$data$covariate_rasters
  points_vec = terra::vect(grid_area, geom = names(grid_area))
  extracted_covariates = terra::extract(covariates, points_vec)
  prediction_data = cbind(grid_area, cov = extracted_covariates)
  linear_pred =  as.numeric(matrix(pars$slope, nrow=1) %*% t(as.matrix(prediction_data[,4:(dim(prediction_data)[2])], nrow = dim(prediction_data)[1]))) + pars$intercept
  
  
  # Spatial objects
  field = (Amatrix %*% pars$nodemean)[,1]

  if(predict_iid) {
    tmp_shp <- output$data$polygon_shapefile
    ID_polygon <- st_join(st_as_sf(grid_area, coords = c("x","y"), crs = st_crs(tmp_shp)),
                          tmp_shp)
    
    # IID effect
      iideffect_sd <- 1/sqrt(exp(pars$iideffect_log_tau))
      iideffect_est = pars$iideffect[ID_polygon$ID]
      iideffect_est = ifelse(is.na(iideffect_est), rnorm(1,0,iideffect_sd), iideffect_est)
    }else{
      iideffect_est = 0 
    }
  
  
  pred = linear_pred + field + iideffect_est
  
   return(pred)
}








