#######################################################################
#                   Fit model function                                #
#######################################################################
# This function should be called to fit the model
# Input of the function
#   - data: result of prepare_data function
#   - cov.method: "exponential", "matern", "circular" or "spherical" 
#   - knots: Spatial knots (result of find_knots function or a matrix with x and y coordinates as columns)
#   - offset: Offset per area for approximate method and offset per grid cell for exact method
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

Krig.fit <- function(data, cov.method, knots = knots, offset = offset,
                      multiple.init = T, nr_cores = 9, parallel = F, exact = F, spat_init = NULL,
                      n_init = 25, sigma_init = NULL, pen_init = NULL, B = 30){
  
  if (exact){
    source("functions/Krig_exact.R")
    fit <- Krig.exact(data, cov.method, knots, offset, multiple.init, nr_cores, parallel, spat_init,
                      n_init, sigma_init, pen_init, B)
  }else{
    source("functions/Krig_approx.R")
    fit <- Krig.approx(data, cov.method, knots, offset, multiple.init, nr_cores, parallel, spat_init,
                       n_init, sigma_init, pen_init, B)
  }

  return(fit)

}
