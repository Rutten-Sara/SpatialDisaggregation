#######################################################################
#                      Fit model exact                                #
#######################################################################
# This function is called by Krig.fit if the exact method is chosen
# Input of the function
#   - data: result of prepare_data function
#   - cov.method: "exponential", "matern", "circular" or "spherical" 
#   - knots: Spatial knots (result of find_knots function or a matrix with x and y coordinates as columns)
#   - offset: Offset per grid cell
#   - multiple.init: True to restart the algorithm with different initial values for spatial range (recommended)
#   - nr_cores: Number of cores in case of parallel computation (only with multiple.init = T)
#   - parallel: True to allow for parallel computation
#   - exact: True for exact calculation, otherwise first order Taylor approximation is used
#   - spat_init: Initial values for spatial range parameter log(rho)
#                If NULL, n_init initial values are constructed based on the maximum spatial distance
#                If provided, this overwrites n_init
#   - n_init: Number of initial values spat_init to generate
#   - sigma_init: Initial value for spatial variance parameter (if not provided, estimated based on GAM)
#   - pen_init: Initial value for penalty parameter of smooth terms (vector of same length as number of smooth parameters)
#               If not provided, estimated based on GAM
#   - B: number of B-spline basis functions for smooth covariate effects

Krig.exact <- function(data, cov.method, knots = knots, offset = offset,
                       multiple.init = F, nr_cores = 9, parallel = F, spat_init,
                       n_init, sigma_init, pen_init, B){
  
  library(dlnm)
  if (!require("Matrix", character.only = TRUE)) {
    message(paste("Package Matrix", "is required but not installed."))
  }
  
  if (!is.null(data$areal_cov)){
    print("This is the exact method, areal covariates will be ignored.")
  }
  
  y <- data$y
  map <- data$map
  data_points = data$rstr_centers 
  weight = data$weight

  ID = 1:dim(map)[1]
  
  coord <- cbind(data_points$x,data_points$y)
  xlin = data.frame(data_points[,grep("linear", names(data_points))])
  
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


  weighted_precision <- function(w){
    return(sum(w)^2/sum(w^2))
  }
 
  X_linear <- cbind(rep(1,length(data_points$ID)), as.matrix(xlin))
  

  
  Z.rand <- Matrix::sparse.model.matrix(~ as.factor(data_points$ID) + 0)
  q.rand <- dim(Z.rand)[2]
  
  Weight_matrix_data <- data.frame(ID =  data_points$ID, weight = weight) #%>%
  
  weight_matrix <- Matrix::sparseMatrix(
    i = 1:length(data_points$ID),               
    j = data_points$ID,         
    x = Weight_matrix_data$weight,
    dims = c(length(data_points$ID), length(ID))
  )
  
  C <- function(v) {
    if (smooth_cov){
      result = Matrix::Matrix(cbind(X_linear,do.call(cbind,Csmooth), cbind(covs$X,covs$Zw(v)),Z.rand), sparse = T)
    } else{
      result = Matrix::Matrix(cbind(X_linear, cbind(covs$X,covs$Zw(v)),Z.rand), sparse = T)
    }
    return(result) #Z.rand
  }
  
  
  # Starting values sigma and penalty
  if(smooth_cov){
    if(is.null(sigma_init) | is.null(pen_init)){
      
      Weight_matrix_data_init <- data.frame(ID =  data_points$ID, weight = weight) %>%
        group_by(ID)%>%
        mutate(tot_weight = sum(weight))%>%
        ungroup()%>%
        mutate(full_weight = weight/tot_weight)
      
      
      weight_matrix_init <- Matrix::sparseMatrix(
        i = 1:length(data_points$ID),               
        j = data_points$ID,         
        x = Weight_matrix_data_init$full_weight,
        dims = c(length(data_points$ID), length(ID))
      )
      
      X_linear_init <- cbind(rep(1,length(ID)), Matrix::t(weight_matrix_init)%*%as.matrix(xlin))
      offset_agg <- as.numeric(Matrix::t(weight_matrix)%*%offset)
      
      discrete_Csmooth <- as.matrix(Matrix::t(weight_matrix_init)%*%do.call(cbind,Csmooth))
      model_smooth <- gam(y ~-1+as.matrix(X_linear_init) + discrete_Csmooth + offset(log(offset_agg)), 
                          family = "poisson")
      if(is.null(sigma_init)){
        sigma_init = c(1/exp(1),min(mean(model_smooth$residuals^2),1))
      }
    }
  }else{
    if(is.null(sigma_init)){
      Weight_matrix_data_init <- data.frame(ID =  data_points$ID, weight = weight) %>%
        group_by(ID)%>%
        mutate(tot_weight = sum(weight))%>%
        ungroup()%>%
        mutate(full_weight = weight/tot_weight)
      
      
      weight_matrix_init <- Matrix::sparseMatrix(
        i = 1:length(data_points$ID),               
        j = data_points$ID,         
        x = Weight_matrix_data_init$full_weight,
        dims = c(length(data_points$ID), length(ID))
      )
      
      X_linear_init <- cbind(rep(1,length(ID)), Matrix::t(weight_matrix_init)%*%as.matrix(xlin))
      offset_agg <- as.numeric(Matrix::t(weight_matrix)%*%offset)
      
      model_linear <- glm(y ~-1+as.matrix(X_linear_init) + offset(log(offset_agg)), family = "poisson")
      sigma_init = c(1/exp(1),min(mean(model_linear$residuals^2),1))
    }
  }
  
  
  
  
  # Prior precision fixed effect
  zeta <- 1e-05 
  zeta_intercept <- 1e-05
  zeta_slope <- 1e-05

  
   Gv <- function(v) {
    continuous_error = data.frame(ID = data_points$ID, weight = weight*(offset-1e-16)+1e-16) 
    discrete_error = continuous_error %>% group_by(ID) %>%
      summarize(prec = weighted_precision(weight))
    return(Matrix::Diagonal(n = q.rand, x = exp(v)*discrete_error$prec))
  }
 
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
  
  ################# Estimation ###################
  # Hyperparameters for Gamma prior of delta
  a <- b <- 1e-05
  # nu prior parameter for the penalty
  nu <- 3

  
  
  Cv_xi <- function(xi, Cv){
    Cvxi <- Matrix::t(weight_matrix)%*%exp(Matrix::Matrix(Cv %*% xi)+log(offset))
    return(Cvxi)
  }



  # Log conditional posterior of xi given v
  logpxi <- function(xi, Cv, Qv, Cvxi) {
    value <- sum(y * log(Cvxi) - Cvxi) - .5 * Matrix::t(xi) %*% Qv %*% xi
    as.numeric(value)
  }


  # Gradient of xi
  grad.logpxi <- function(xi,Cv, Qv, Cvxi){
    dlambda_beta <- Matrix::t(Cv)%*%(Matrix::Diagonal(x=exp(as.numeric(Matrix::Matrix(Cv %*% xi)+log(offset)))))
    dlambda_beta_weighted = dlambda_beta%*%weight_matrix
    value= dlambda_beta_weighted%*%Matrix::Matrix(y*(Cvxi)^(-1)-1)- Matrix::Matrix(Qv%*%xi)
    return(as.numeric(value))
  }
  
  num_area = data_points %>%group_by(ID) %>% summarize(tot_points = n()) %>% ungroup()
  Hess.logpxi <- function(xi,Cv, Qv, Cvxi){
    Diag_matrix = Matrix::Diagonal(x=exp(as.numeric(Matrix::Matrix(Cv %*% xi)+log(offset))))
    dlambda_beta <-  Matrix::t(Cv)%*%Diag_matrix
    dlambda_beta_weighted = dlambda_beta%*%weight_matrix
    value1 = dlambda_beta_weighted%*%Matrix::Diagonal(x = as.numeric(-y*(Cvxi)^(-2))) %*% Matrix::t(dlambda_beta_weighted)
    
    weight_hess = Matrix::Matrix(y*(Cvxi)^(-1)-1)
    nr_rep =rep(as.numeric(weight_hess), times = num_area$tot_points)
    first_part = Matrix::Diagonal(x=sqrt(exp(as.numeric(Matrix::Matrix(Cv %*% xi)+log(offset)))))%*%Cv
    first_part_weight = Matrix::Diagonal(x=sqrt(abs(as.numeric(Weight_matrix_data$weight))))%*%first_part
    first_part_weight_scaled <- first_part_weight * as.numeric(nr_rep)
    value2 <- Matrix::crossprod(first_part_weight_scaled, first_part_weight)
    
    value = (value1+value2) - Qv + Matrix::Diagonal(n = dimxi, x = 1e-05)
    return(value)
  }
  
  Hess.logpxi_stable <- function(xi,Cv, Qv, Cvxi){
    Diag_matrix = Matrix::Diagonal(x=exp(as.numeric(Matrix::Matrix(Cv %*% xi)+log(offset))))
    dlambda_beta <-  Matrix::t(Cv)%*%Diag_matrix
    dlambda_beta_weighted = dlambda_beta%*%weight_matrix
    value1 = dlambda_beta_weighted%*%Matrix::Diagonal(x = as.numeric(-y*(Cvxi)^(-2))) %*% Matrix::t(dlambda_beta_weighted)
    
    weight_hess = Matrix::Matrix(y*(Cvxi)^(-1)-1)
    nr_rep =rep(as.numeric(weight_hess), times = num_area$tot_points)
    first_part = Matrix::Diagonal(x=sqrt(exp(as.numeric(Matrix::Matrix(Cv %*% xi)+log(offset)))))%*%Cv
    first_part_weight = Matrix::Diagonal(x=sqrt(abs(as.numeric(Weight_matrix_data$weight))))%*%first_part
    first_part_weight_scaled <- first_part_weight * as.numeric(nr_rep)
    value2 <- Matrix::crossprod(first_part_weight_scaled, first_part_weight)
    
    value = (value1+value2) - Qv + Matrix::Diagonal(n = dimxi, x = 1e-02)
    return(value)
  }

  
  # Laplace approximation to conditional posterior of xi
  Laplace <- function(xi0, Cv, Qv){
    
    epsilon <- 1e-03 # Stop criterion
    maxiter <- 100   # Maximum iterations
    iter <- 0        # Iteration counter
    
    for (k in 1:maxiter) {
      Cvxi = Cv_xi(xi0, Cv)
      dxi <- as.numeric((-1) * Matrix::solve(Matrix::Matrix(Hess.logpxi(xi0, Cv, Qv,Cvxi)),
                                     Matrix::Matrix(grad.logpxi(xi0, Cv, Qv, Cvxi)), tol = 1e-25))

      
      xi.new <- xi0 + dxi
      step <- 1
      iter.halving <- 1
      logpxi.current <- logpxi(xi0, Cv, Qv, Cvxi)
      while (max(logpxi(xi.new, Cv, Qv, Cv_xi(xi.new,Cv)),logpxi.current,na.rm=T)==logpxi.current) {
        step <- step * .5
        xi.new <- xi0 + (step * dxi)
        iter.halving <- iter.halving + 1
        if (iter.halving > 1000) {
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
  
  Laplace_stable <- function(xi0, Cv, Qv){
    
    epsilon <- 1e-03 # Stop criterion
    maxiter <- 100   # Maximum iterations
    iter <- 0        # Iteration counter
    
    for (k in 1:maxiter) {
      Cvxi = Cv_xi(xi0, Cv)
      dxi <- as.numeric((-1) * Matrix::solve(Matrix::Matrix(Hess.logpxi_stable(xi0, Cv, Qv,Cvxi)),
                                             Matrix::Matrix(grad.logpxi(xi0, Cv, Qv, Cvxi)), tol = 1e-25))
      
      
      xi.new <- xi0 + dxi
      step <- 1
      iter.halving <- 1
      logpxi.current <- logpxi(xi0, Cv, Qv, Cvxi)
      while (max(logpxi(xi.new, Cv, Qv, Cv_xi(xi.new,Cv)),logpxi.current,na.rm=T)==logpxi.current) {
        step <- step * .5
        xi.new <- xi0 + (step * dxi)
        iter.halving <- iter.halving + 1
        if (iter.halving > 1000) {
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
    vs <- v[-c(length(v)-2,length(v))] 
    vphi <- v[length(v)]
    vu <- v[length(v)-2]
    Cv <- C(v)
    Cvxi <- Cv_xi(xi_hat, Cv)
    
    Qv <- Q(v)
    
    dlambda_beta <- Matrix::t(Cv)%*%(Matrix::Diagonal(x=exp(as.numeric(Matrix::Matrix(Cv %*% xi_hat)+log(offset)))))
    dlambda_beta_weighted = dlambda_beta%*%weight_matrix
    
    value1_hessian = dlambda_beta_weighted%*%Matrix::Diagonal(x = as.numeric(-y*(Cvxi)^(-2))) %*% Matrix::t(dlambda_beta_weighted)
    
    weight_hess = Matrix::Matrix(y*(Cvxi)^(-1)-1)
    nr_rep =rep(as.numeric(weight_hess), times = num_area$tot_points)
    first_part = Matrix::Diagonal(x=sqrt(exp(as.numeric(Matrix::Matrix(Cv %*% xi_hat)+log(offset)))))%*%Cv
    first_part_weight = Matrix::Diagonal(x=sqrt(abs(as.numeric(Weight_matrix_data$weight))))%*%first_part
    first_part_weight_scaled <- first_part_weight * as.numeric(nr_rep)
    value2_hessian <- Matrix::crossprod(first_part_weight_scaled, first_part_weight)

    
    value_hessian = -(value1_hessian+value2_hessian - Qv) 
    
    a2=0.5*as.numeric(Matrix::determinant(covs$Covw(v), logarithm = TRUE)$modulus)
    a3=-0.5*as.numeric(Matrix::determinant(value_hessian, logarithm = TRUE)$modulus)
    a1 <- sum(y * log(Cvxi) - Cvxi) - .5 * t(xi_hat) %*% Qv %*% xi_hat
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
    a1 <- sum(y * log(Cvxi) - Cvxi) - .5 * t(xi_hat) %*% Qv %*% xi_hat
    a4 <- sum(0.5 * (rk + nu) * vs)
    au <- 0.5 * (nu+q.rand) * vu-(0.5*nu + a)*log(b + 0.5*nu*exp(vu))
    a5 <- - sum(((0.5 * nu) + a) *  log(0.5*(nu * exp(vs)) + b))
    a6 <- (0.5 * nu) * vphi - ((0.5 * nu) + a) *  log(0.5*(nu * exp(vphi)) + b)
    value <- a1+a2+a3+a4+a5+au+a6
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

  
  Cvxi_mode = Cv_xi(xi_estim, Cv_mode)
  
  Sigma <- nearPD(Matrix::solve(-Matrix::Matrix(Hess.logpxi_stable(xi = xi_estim,
                                                                   Cv = Cv_mode,
                                                                   Qv = Qv_mode, Cvxi_mode))))$mat
  

  
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
        
        Cvxi <- Cv_xi(xi_hat, Cv_init)
        
        dlambda_beta <- Matrix::t(Cv_init)%*%(Matrix::Diagonal(x=exp(as.numeric(Matrix::Matrix(Cv_init %*% xi_hat)+log(offset)))))
        dlambda_beta_weighted = dlambda_beta%*%weight_matrix
        
        value1_hessian = dlambda_beta_weighted%*%Matrix::Diagonal(x = as.numeric(-y*(Cvxi)^(-2))) %*% Matrix::t(dlambda_beta_weighted)
        
        weight_hess = Matrix::Matrix(y*(Cvxi)^(-1)-1)
        nr_rep =rep(as.numeric(weight_hess), times = num_area$tot_points)
        first_part = Matrix::Diagonal(x=sqrt(exp(as.numeric(Matrix::Matrix(Cv_init %*% xi_hat)+log(offset)))))%*%Cv_init
        first_part_weight = Matrix::Diagonal(x=sqrt(abs(as.numeric(Weight_matrix_data$weight))))%*%first_part
        first_part_weight_scaled <- first_part_weight * as.numeric(nr_rep)
        value2_hessian <- Matrix::crossprod(first_part_weight_scaled, first_part_weight)
        
        value_hessian_init = value1_hessian+value2_hessian
        
        a2=0.5*as.numeric(Matrix::determinant(covs$Covw(c(v_init, v_init_3_rho)), logarithm = TRUE)$modulus)
        
        
        v_mode <- optim(par = v_init, 
                        fn = log_pv_partly, 
                        xi_hat = xi_hat,
                        Cvxi = Cvxi,
                        vphi = v_init_3_rho,
                        method = "Nelder-Mead", 
                        control = list(fnscale = -1, 
                                       reltol = 1e-8))$par
        
        
        Cv_mode <- C(c(v_mode, v_init_3_rho))
        Qv_mode <- Q(c(v_mode, v_init_3_rho))

        
        # Mode a posteriori estimate for xi
        xi_estim <- Laplace_stable(xi_hat, 
                                   Cv = Cv_mode, 
                                   Qv = Qv_mode)
        
        Cvxi_mode = Cv_xi(xi_hat, Cv_mode)
        
        dlambda_beta <- Matrix::t(Cv_mode)%*%(Matrix::Diagonal(x=exp(as.numeric(Matrix::Matrix(Cv_mode %*% xi_hat)+log(offset)))))
        dlambda_beta_weighted = dlambda_beta%*%weight_matrix
        
        value1_hessian = dlambda_beta_weighted%*%Matrix::Diagonal(x = as.numeric(-y*(Cvxi_mode)^(-2))) %*% Matrix::t(dlambda_beta_weighted)
        
        weight_hess = Matrix::Matrix(y*(Cvxi_mode)^(-1)-1)
        nr_rep =rep(as.numeric(weight_hess), times = num_area$tot_points)
        first_part = Matrix::Diagonal(x=sqrt(exp(as.numeric(Matrix::Matrix(Cv_mode %*% xi_hat)+log(offset)))))%*%Cv_mode
        first_part_weight = Matrix::Diagonal(x=sqrt(abs(as.numeric(Weight_matrix_data$weight))))%*%first_part
        first_part_weight_scaled <- first_part_weight * as.numeric(nr_rep)
        value2_hessian <- Matrix::crossprod(first_part_weight_scaled, first_part_weight)
        
        value_hessian_init = value1_hessian+value2_hessian
      
        log_xi_v[v_init_3] = log_pv_partly(v_mode, xi_hat, Cvxi_mode, v_init_3_rho)
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
    Cvxi_mode = Cv_xi(xi_estim, Cv_mode)
    
    Sigma <- nearPD(Matrix::solve(-Matrix::Matrix(Hess.logpxi_stable(xi = xi_estim,
                                                                     Cv = Cv_mode,
                                                                     Qv = Qv_mode, Cvxi_mode))))$mat
    
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
      
      source("functions/Krig_exact.R")
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
          
          Cvxi <- Cv_xi(xi_hat, Cv_init)
          
          dlambda_beta <- Matrix::t(Cv_init)%*%(Matrix::Diagonal(x=exp(as.numeric(Matrix::Matrix(Cv_init %*% xi_hat)+log(offset)))))
          dlambda_beta_weighted = dlambda_beta%*%weight_matrix
          
          value1_hessian = dlambda_beta_weighted%*%Matrix::Diagonal(x = as.numeric(-y*(Cvxi)^(-2))) %*% Matrix::t(dlambda_beta_weighted)
          
          weight_hess = Matrix::Matrix(y*(Cvxi)^(-1)-1)
          nr_rep =rep(as.numeric(weight_hess), times = num_area$tot_points)
          first_part = Matrix::Diagonal(x=sqrt(exp(as.numeric(Matrix::Matrix(Cv_init %*% xi_hat)+log(offset)))))%*%Cv_init
          first_part_weight = Matrix::Diagonal(x=sqrt(abs(as.numeric(Weight_matrix_data$weight))))%*%first_part
          first_part_weight_scaled <- first_part_weight * as.numeric(nr_rep)
          value2_hessian <- Matrix::crossprod(first_part_weight_scaled, first_part_weight)
          
          value_hessian_init = value1_hessian+value2_hessian
          
          a2=0.5*as.numeric(Matrix::determinant(covs$Covw(c(v_init, v_init_3_rho)), logarithm = TRUE)$modulus)
          
          
          v_mode <- optim(par = v_init, 
                          fn = log_pv_partly, 
                          xi_hat = xi_hat,
                          Cvxi = Cvxi,
                          vphi = v_init_3_rho,
                          method = "Nelder-Mead", 
                          control = list(fnscale = -1, 
                                         reltol = 1e-8))$par
          
          
          Cv_mode <- C(c(v_mode, v_init_3_rho))
          Qv_mode <- Q(c(v_mode, v_init_3_rho))
          
          # Mode a posteriori estimate for xi
          xi_estim <- Laplace_stable(xi_hat, 
                              Cv = Cv_mode, 
                              Qv = Qv_mode)
          
          Cvxi_mode = Cv_xi(xi_hat, Cv_mode)
          
          dlambda_beta <- Matrix::t(Cv_mode)%*%(Matrix::Diagonal(x=exp(as.numeric(Matrix::Matrix(Cv_mode %*% xi_hat)+log(offset)))))
          dlambda_beta_weighted = dlambda_beta%*%weight_matrix
          
          value1_hessian = dlambda_beta_weighted%*%Matrix::Diagonal(x = as.numeric(-y*(Cvxi_mode)^(-2))) %*% Matrix::t(dlambda_beta_weighted)
          
          weight_hess = Matrix::Matrix(y*(Cvxi_mode)^(-1)-1)
          nr_rep =rep(as.numeric(weight_hess), times = num_area$tot_points)
          first_part = Matrix::Diagonal(x=sqrt(exp(as.numeric(Matrix::Matrix(Cv_mode %*% xi_hat)+log(offset)))))%*%Cv_mode
          first_part_weight = Matrix::Diagonal(x=sqrt(abs(as.numeric(Weight_matrix_data$weight))))%*%first_part
          first_part_weight_scaled <- first_part_weight * as.numeric(nr_rep)
          value2_hessian <- Matrix::crossprod(first_part_weight_scaled, first_part_weight)
          
          value_hessian_init = value1_hessian+value2_hessian
          
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
      Cvxi_mode = Cv_xi(xi_estim, Cv_mode)
      
      Sigma <- nearPD(Matrix::solve(-Matrix::Matrix(Hess.logpxi_stable(xi = xi_estim,
                                                         Cv = Cv_mode,
                                                         Qv = Qv_mode, Cvxi_mode))))$mat

      

      
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
                 "exact" = T)
                 
}


