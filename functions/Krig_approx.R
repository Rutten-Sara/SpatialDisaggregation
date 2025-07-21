#######################################################################
#       Fit model based on first order Taylor approximation           #
#######################################################################
# This function is called by Krig.fit if the approximate method is chosen
# Input of the function
#   - data: result of prepare_data function
#   - cov.method: "exponential", "matern", "circular" or "spherical" 
#   - knots: Spatial knots (result of find_knots function or a matrix with x and y coordinates as columns)
#   - offset: Offset per area 
#   - multiple.init: True to restart the algorithm with different initial values for spatial range (recommended)
#   - nr_cores: Number of cores in case of parallel computation (only with multiple.init = T)
#   - parallel: True to allow for parallel computation
#   - spat_init: Initial values for spatial range parameter log(rho)
#                If NULL, n_init initial values are constructed based on the maximum spatial distance
#                If provided, this overwrites n_init
#   - n_init: Number of initial values spat_init to generate
#   - sigma_init: Initial value for spatial variance parameter (if not provided, estimated based on GAM)
#   - pen_init: Initial value for penalty parameter of smooth terms (vector of same length as number of smooth parameters)
#               If not provided, estimated based on GAM
#   - B: number of B-spline basis functions for smooth covariate effects


Krig.approx <- function(data, cov.method, knots = knots, offset = offset,
                      multiple.init = F, nr_cores = 9, parallel = F, spat_init,
                      n_init, sigma_init, pen_init, B){
  library(dlnm)
  if (!require("Matrix", character.only = TRUE)) {
    message(paste("Package Matrix", "is required but not installed."))
  }
  
  y <- data$y
  map <- data$map
  data_points = data$rstr_centers 
  weight = data$weight

  ID = 1:dim(map)[1]
  
  coord <- cbind(data_points$x,data_points$y)
  xlin_cont = as.matrix(data.frame(data_points[,grep("linear", names(data_points))]))

  
  smooth_cov = sum(grepl("smooth",names(data_points)))>0

  if (smooth_cov){
    
  all_smooth_cov = data.frame(data_points[grepl("smooth",names(data_points))])
  
  if(is.null(pen_init)){
    pen_init = rep(5,dim(all_smooth_cov)[2])
  }
  
  sm = list()
  D = list()
  Csmooth = list()
  
  nr_smooth = dim(all_smooth_cov)[2]
  
  for (s in 1:dim(all_smooth_cov)[2]){  
  xsmooth = all_smooth_cov[,s]
  ind = !is.na(xsmooth)
  
  data_points = data_points[ind,]
  coord = coord[ind,]
  xsmooth = xsmooth[ind]
  weight = weight[ind]
  xlin_cont = xlin_cont[ind,]
  
  sm[[s]] <- ps(xsmooth, df = B)
  nbasis <- B
  Csmooth[[s]] <- Matrix::Matrix(sm[[s]])
  
  D[[s]] = Matrix::t(Matrix::Matrix(attr(sm[[s]],"S")))%*%Matrix::Matrix(attr(sm[[s]],"S"))+ Matrix::Diagonal(n = B, x = 1e-12)
  
  }
  }
  
  
  n <- dim(data_points)[1]
  
  data_points = data_points %>% group_by(ID) %>%
    mutate(ID = cur_group_id()) %>% ungroup()

  # Spatial component
  knots <- knots
  K <- dim(knots)[1]
  covs <- cov.spatial(coord = coord,
                      K = K,
                      cov.method = cov.method,
                      knots = knots)

  
  #weighted_average <- function(col,w){
  #  return(sum(w*col)/sum(w))
  #}
  
  weighted_precision <- function(w){
    return(sum(w)^2/sum(w^2))
  }
 

  
  Z.rand <- Matrix::sparse.model.matrix(~ as.factor(ID) + 0)
  q.rand <- dim(Z.rand)[2]

  Weight_matrix_data <- data.frame(ID =  data_points$ID, weight = weight) %>%
    group_by(ID)%>%
    mutate(tot_weight = sum(weight))%>%
    ungroup()%>%
    mutate(full_weight = weight/tot_weight)
  
  
  weight_matrix <- Matrix::sparseMatrix(
    i = 1:length(data_points$ID),               
    j = data_points$ID,         
    x = Weight_matrix_data$full_weight,
    dims = c(length(data_points$ID), length(ID))
  )
  
  X_linear <- cbind(rep(1,length(ID)), Matrix::t(weight_matrix)%*%xlin_cont)
  
  C <- function(v) {
    continuous_C <- Matrix::Matrix(cbind(covs$X,covs$Zw(v)))
    discrete_C <- Matrix::t(weight_matrix)%*%continuous_C
    if (smooth_cov){
    discrete_Csmooth <- Matrix::t(weight_matrix)%*%do.call(cbind,Csmooth)
    result = Matrix::Matrix(cbind(as.matrix(X_linear),discrete_Csmooth, discrete_C, Z.rand), sparse = T)
    }else{
      result = Matrix::Matrix(cbind(as.matrix(X_linear),discrete_C, Z.rand), sparse = T)
    }
    return(result)
  }
  
  
  Gv <- function(v) {
    continuous_error = data.frame(ID = data_points$ID, weight = weight) 
    discrete_error = continuous_error %>% group_by(ID) %>%
      summarize(prec = weighted_precision(weight))
    return(Matrix::Diagonal(n = q.rand, x = exp(v)*discrete_error$prec))
  }
  
  # Starting values sigma and penalty
  if(smooth_cov){
    if(is.null(sigma_init) | is.null(pen_init)){
    discrete_Csmooth <- as.matrix(Matrix::t(weight_matrix)%*%do.call(cbind,Csmooth))
    model_smooth <- gam(y ~-1+as.matrix(X_linear) + discrete_Csmooth + offset(log(offset)), 
                        family = "poisson")
    if(is.null(sigma_init)){
    sigma_init = c(1/exp(1),min(mean(model_smooth$residuals^2),1))
    }
    }
  }else{
    if(is.null(sigma_init)){
      model_linear <- glm(y ~-1+as.matrix(X_linear) + offset(log(offset)), family = "poisson")
      sigma_init = c(1/exp(1),min(mean(model_linear$residuals^2),1))
    }
  }
  
  ################# Estimation ###################
  # Hyperparameters for Gamma prior of delta
  a <- b <- 1e-05
  # nu prior parameter for the penalty
  nu <- 3
  # Prior precision fixed effect
  zeta <-  1e-05 
  zeta_intercept <- 1e-05
  zeta_slope <- 1e-05
 
  
  Q <- function(v) {
    if (smooth_cov){
      result_matrix <- Matrix::bdiag(c(list(Matrix::Diagonal(n = dim(X_linear)[2],
                                                      x = c(zeta_intercept, rep(zeta_slope,dim(X_linear)[2]-1)))), 
                                    Map(function(di, vi) exp(vi) * Matrix::Matrix(di), D, v[1:nr_smooth]),
                                   list(Matrix::Diagonal(n = 2, x = zeta),
                                   exp(v[length(v) - 1]) * covs$Covw(v),
                                   Gv(v[length(v)-2]))))
    }else{
      result_matrix <- Matrix::bdiag(Matrix::Diagonal(n = dim(X_linear)[2],
                                                      x = c(zeta_intercept, rep(zeta_slope,dim(X_linear)[2]-1))),
                                     Matrix::Diagonal(n = 2, x = zeta),
                                     exp(v[length(v) - 1]) * covs$Covw(v),
                                     Gv(v[length(v)-2]))
    }
    
    return(result_matrix)
  }
  
  dimxi <- ifelse(smooth_cov, K + 2 + dim(X_linear)[2] + q.rand + nr_smooth*B,
                  K + 2 + dim(X_linear)[2] + q.rand)
  if(smooth_cov){
    rk = c(rep(B,nr_smooth), K)
  }else{
    rk = K
  }
  

  
  
  # Log conditional posterior of xi given v
  logpxi <- function(xi, Cv, Qv) {
    Cvxi <- as.numeric(Matrix::Matrix(Cv %*% xi)+log(offset))
    value <- sum(y * Cvxi - exp(Cvxi)) - .5 * Matrix::t(xi) %*% Qv %*% xi
    as.numeric(value)
  }
  
  # Gradient of xi
  grad.logpxi <- function(xi,Cv, Qv){
    value <- Matrix::t(Cv)%*%Matrix::Matrix((y - exp(as.numeric(Matrix::Matrix(Cv %*% xi)+log(offset))))) - Matrix::Matrix(Qv%*%xi)
    as.numeric(value)
  }
  
  Hess.logpxi <- function(xi,Cv, Qv){
    value <- - Matrix::t(Cv)%*%Matrix::Diagonal(x = exp(as.numeric(Matrix::Matrix(Cv %*% xi)+log(offset))))%*%Cv - Qv + Matrix::Diagonal(n = dimxi, x = 1e-05)
    value
  }

  
  # Laplace approximation to conditional posterior of xi
  Laplace <- function(xi0, Cv, Qv){
    
    epsilon <- 1e-03 # Stop criterion
    maxiter <- 100   # Maximum iterations
    iter <- 0        # Iteration counter
    
    for (k in 1:maxiter) {
      dxi <- as.numeric((-1) * Matrix::solve(Matrix::Matrix(Hess.logpxi(xi0, Cv, Qv)),
                        Matrix::Matrix(grad.logpxi(xi0, Cv, Qv)), tol = 1e-25))
      xi.new <- xi0 + dxi
      step <- 1
      iter.halving <- 1
      logpxi.current <- logpxi(xi0, Cv, Qv)
      while (logpxi(xi.new, Cv, Qv) <= logpxi.current) {
        step <- step * .5
        xi.new <- xi0 + (step * dxi)
        iter.halving <- iter.halving + 1
        if (iter.halving > 30) {
          break
        }
      }
      dist <- sqrt(sum((xi.new - xi0) ^ 2))
      iter <- iter + 1
      xi0 <- xi.new
      if(dist < epsilon) break
    }
    
    xistar <- xi0 
    return(xistar)
  }
  
  
  log_pv <- function(v){
    vs <- v[-c(length(v)-2, length(v))]
    vphi <- v[length(v)]
    vu <- v[length(v)-2]
    Cv <- C(v)
    Cvxi <- as.numeric(Matrix::Matrix(Cv %*% xi_hat)+log(offset))
    Qv <- Q(v)
    a1 <- sum(y * Cvxi - exp(Cvxi)) - .5 * t(xi_hat) %*% Qv %*% xi_hat
    a2=0.5*as.numeric(Matrix::determinant(covs$Covw(v), logarithm = TRUE)$modulus)
    a3=-0.5*as.numeric(Matrix::determinant(Matrix::Matrix(Matrix::t(Cv)%*%Matrix::Diagonal(x = exp(Cvxi))%*%Cv + Qv), logarithm = TRUE)$modulus)
    a4 <- sum(0.5 * (rk + nu) * vs)
    au <- 0.5 * (nu+q.rand) * vu-(0.5*nu + a)*log(b + 0.5*nu*exp(vu))
    a5 <- - sum(((0.5 * nu) + a) *  log(0.5*(nu * exp(vs)) + b))
    a6 <- (0.5 * nu) * vphi - ((0.5 * nu) + a) *  log(0.5*(nu * exp(vphi)) + b)
    value <- a1+a2+a3+a4+a5+au+a6
    return(as.numeric(value))
  }
  
  log_pv_partly <- function(v, xi_hat, Cvxi, vphi){
    vs <- v[-c(length(v)-1)] 
    vu <- v[length(v)-1]
    
    Qv <- Q(c(v, v_init_3_rho))
    
    value_hessian = -(value_hessian_init - Qv) 
    
    a3=-0.5*as.numeric(Matrix::determinant(value_hessian, logarithm = TRUE)$modulus)
    a1 <- sum(y * Cvxi - exp(Cvxi)) - .5 * t(xi_hat) %*% Qv %*% xi_hat
    a4 <- sum(0.5 * (rk + nu) * vs)
    au <- 0.5 * (nu+q.rand) * vu-(0.5*nu + a)*log(b + 0.5*nu*exp(vu))
    a5 <- - sum(((0.5 * nu) + a) *  log(0.5*(nu * exp(vs)) + b))
    value <- a1+a2+a3+a4+a5+au
    return(as.numeric(value))
  }
  
  # Multiple starting values for v_init?
  
  if (!multiple.init){
    if(smooth_cov){
      if(is.null(sigma_init)){
          v_init <- c(as.numeric(pen_init), 1, 0, log(1/max(covs$Z.dist)))

      }else{
         v_init <- c(as.numeric(pen_init), log(1/sigma_init[1]), log(1/sigma_init[2]), log(1/max(covs$Z.dist)))
        }
    }else{
      if(is.null(sigma_init)){
        v_init <- c(1, 0, log(1/max(covs$Z.dist)))
      }else{
        v_init <- c(log(1/sigma_init[1]), log(1/sigma_init[2]), log(1/max(covs$Z.dist)))
      }
    }
  

  Cv_init <- C(v_init)
  Qv_init <- Q(v_init)

  # Initial estimate for xi
  xi_hat <- Laplace(xi0 = rep(0,dimxi), 
                   Cv = Cv_init, 
                   Qv = Qv_init)

  
  v_mode <- optim(par = v_init, 
                  fn = log_pv, 
                  method = "Nelder-Mead", 
                  control = list(fnscale = -1, 
                                 reltol = 1e-8))$par

  Cv_mode <- C(v_mode)
  Qv_mode <- Q(v_mode)
  
  # Mode a posteriori estimate for xi
  xi_estim <- Laplace(xi_hat, 
                            Cv = Cv_mode, 
                            Qv = Qv_mode)
  Sigma <- nearPD(-Matrix::solve(Matrix::Matrix(Hess.logpxi(xi = xi_estim,
                                                            Cv = Cv_mode,
                                                            Qv = Qv_mode)), tol = 1e-25))$mat
  }else{
    if(!is.null(spat_init)){
      n_sim = length(spat_init)
      v_init_rho = spat_init
    }else if (!is.null(n_init)){
      v_init_rho_max = 1/max(covs$Z.dist)
      r_ind = 2^seq(0,n_init*0.25, by=0.25)
      v_init_rho = log(r_ind*v_init_rho_max)
      n_sim = length(v_init_rho)
    }else{
      v_init_rho_max = 1/max(covs$Z.dist)
      r_ind = 2^seq(0,(10-2)*0.25, by=0.25)
      v_init_rho = log(r_ind*v_init_rho_max)
      n_sim = length(v_init_rho)
    }
    
    log_xi_v = rep(0,n_sim)
    v_mode_final = list()
    xi_estim_final = list()
    
    if (parallel == F){
      
      for (v_init_3 in 1:n_sim){
      loop_try = tryCatch({
        v_init_3_rho = v_init_rho[v_init_3]
        if(smooth_cov){
          if(is.null(sigma_init)){
              v_init <- c(pen_init, 1, 0)
          }else{
              v_init <- c(pen_init, log(1/sigma_init[1]), log(1/sigma_init[2]))
          }
        }else{
          if(is.null(sigma_init)){
            v_init <- c(1, 0)
          }else{
            v_init <- c(log(1/sigma_init[1]), log(1/sigma_init[2]))
          }
        }
        
        Cv_init <- C(c(v_init, v_init_3_rho))
        Qv_init <- Q(c(v_init, v_init_3_rho))
        
        # Initial estimate for xi
        xi_hat <- Laplace(xi0 = rep(0,dimxi), 
                          Cv = Cv_init, 
                          Qv = Qv_init)
        
        Cvxi <- as.numeric(Matrix::Matrix(Cv_init %*% xi_hat)+log(offset))
        
        a2=0.5*as.numeric(Matrix::determinant(covs$Covw(c(v_init, v_init_3_rho)), logarithm = TRUE)$modulus)
        
        value_hessian_init = - Matrix::Matrix(Matrix::t(Cv_init)%*%Matrix::Diagonal(x = exp(Cvxi))%*%Cv_init)
        v_mode <- optim(par = v_init, 
                        fn = log_pv_partly, 
                        xi_hat = xi_hat,
                        Cvxi = Cvxi,
                        vphi = v_init_3_rho,
                        method = "Nelder-Mead", 
                        control = list(fnscale = -1, 
                                       reltol = 1e-8))$par
        
        
        Cv_mode <-  C(c(v_mode, v_init_3_rho))
        Qv_mode <- Q(c(v_mode, v_init_3_rho))
        
        # Mode a posteriori estimate for xi
        xi_estim <- Laplace(xi_hat, 
                            Cv = Cv_mode, 
                            Qv = Qv_mode)
        
        
        Cvxi_mode = as.numeric(Matrix::Matrix(Cv_mode %*% xi_hat)+log(offset))
        value_hessian_init = - Matrix::Matrix(Matrix::t(Cv_mode)%*%Matrix::Diagonal(x = exp(Cvxi_mode))%*%Cv_mode)
        log_xi_v[v_init_3] =log_pv_partly(v_mode, xi_hat, Cvxi_mode, v_init_3_rho) 
        v_mode_final[[v_init_3]] = c(v_mode, v_init_3_rho)
        xi_estim_final[[v_init_3]] = xi_estim
      }, error= function(err){
        log_xi_v[v_init_3] = -Inf
        v_mode_final[[v_init_3]] = NA
        xi_estim_final[[v_init_3]] = NA
      })
      }


    v_mode = v_mode_final[[which.max(log_xi_v)]]
    xi_estim = xi_estim_final[[which.max(log_xi_v)]]
    Cv_mode <- C(v_mode)
    Qv_mode <- Q(v_mode)
    
    Sigma <- nearPD(-Matrix::solve(Matrix::Matrix(Hess.logpxi(xi = xi_estim,
                                                              Cv = Cv_mode,
                                                              Qv = Qv_mode)), tol = 1e-25))$mat
    }else{
      library(foreach)
      library(doParallel)
      # Number of available cores
      detectCores()
      
      # Make sure that the number provided is less than the number of cores.
      clust <- makeCluster(nr_cores)
      registerDoParallel(clust)
      
      output_loop <- foreach(v_init_3 = 1:n_sim,
                             .packages = c("dlnm",  "Matrix", "dplyr"),
                             .inorder = F) %dopar%
        {
          
          source("functions/Krig_approx.R")
          loop_try = tryCatch({
            v_init_3_rho = v_init_rho[v_init_3]
            if(smooth_cov){
              if(is.null(sigma_init)){
                  v_init <- c(pen_init, 1, 0)
              }else{
                  v_init <- c(pen_init, log(1/sigma_init[1]), log(1/sigma_init[2]))
  
              }
            }else{
              if(is.null(sigma_init)){
                v_init <- c(1, 0)
              }else{
                v_init <- c(log(1/sigma_init[1]), log(1/sigma_init[2]))
              }
            }
            
            Cv_init <- C(c(v_init, v_init_3_rho))
            Qv_init <- Q(c(v_init, v_init_3_rho))
            
            # Initial estimate for xi
            xi_hat <- Laplace(xi0 = rep(0,dimxi), 
                              Cv = Cv_init, 
                              Qv = Qv_init)
            
            Cvxi <- as.numeric(Matrix::Matrix(Cv_init %*% xi_hat)+log(offset))
            
            a2=0.5*as.numeric(Matrix::determinant(covs$Covw(c(v_init, v_init_3_rho)), logarithm = TRUE)$modulus)
            
            value_hessian_init = - Matrix::Matrix(Matrix::t(Cv_init)%*%Matrix::Diagonal(x = exp(Cvxi))%*%Cv_init)
            v_mode <- optim(par = v_init, 
                            fn = log_pv_partly, 
                            xi_hat = xi_hat,
                            Cvxi = Cvxi,
                            vphi = v_init_3_rho,
                            method = "Nelder-Mead", 
                            control = list(fnscale = -1, 
                                           reltol = 1e-8))$par
            
            
            Cv_mode <-  C(c(v_mode, v_init_3_rho))
            Qv_mode <- Q(c(v_mode, v_init_3_rho))
            
            # Mode a posteriori estimate for xi
            xi_estim <- Laplace(xi_hat, 
                                Cv = Cv_mode, 
                                Qv = Qv_mode)
            
            
            Cvxi_mode = as.numeric(Matrix::Matrix(Cv_mode %*% xi_hat)+log(offset))
            value_hessian_init = - Matrix::Matrix(Matrix::t(Cv_mode)%*%Matrix::Diagonal(x = exp(Cvxi_mode))%*%Cv_mode)
            
            log_xi_v_par = log_pv_partly(v_mode, xi_hat, Cvxi_mode, v_init_3_rho)
            v_mode_final_par = c(v_mode, v_init_3_rho)
            xi_estim_final_par = xi_estim
            return (list(log_xi_v = log_xi_v_par, v_mode_final = v_mode_final_par, xi_estim_final = xi_estim_final_par))
          }, error= function(err){
            log_xi_v_par = -Inf
            v_mode_final_par = NA
            xi_estim_final_par = NA
            return (list(log_xi_v = log_xi_v_par, v_mode_final = v_mode_final_par, xi_estim_final = xi_estim_final_par))
          })
          
        }
      max_index = which.max(sapply(output_loop, function(x) x$log_xi_v))
      v_mode = (output_loop[[max_index]])$v_mode_final
      xi_estim = (output_loop[[max_index]])$xi_estim_final
      
      Cv_mode <- C(v_mode)
      Qv_mode <- Q(v_mode)
      
      Sigma <- nearPD(-Matrix::solve(Matrix::Matrix(Hess.logpxi(xi = xi_estim,
                                                         Cv = Cv_mode,
                                                         Qv = Qv_mode)), tol = 1e-25))$mat
      stopCluster(clust)
    }
    
  }
  
  
  if(!smooth_cov){
    sm = NULL
  }
  

  output <- list("xi_estim" = xi_estim,
                 "Sigma" = Sigma,
                 "covs" = covs, 
                 "v_mode" = v_mode, 
                 "K" = K,
                 "sm" = sm,
                 "data_points" = data_points,
                 "exact" = F,
                 "weighted" = data$weighted,
                 "offset" = offset)

}