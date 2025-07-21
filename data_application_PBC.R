rm(list = ls())

#Source the functions
file_list = c("functions/Krig_fit.R","functions/cov.spatial.R", "functions/help_functions.R")
for (file in file_list) {
  source(file)
}


library(classInt)
library(RColorBrewer)
library(sf)
library(sp)
library(dplyr)
library(tmap)
library(raster)
library(ggplot2)
library(exactextractr)
library(mgcv)

wd <- getwd()

# Additional functions disaggregation package, SDALGCP and ATP
source("functions/Other functions/Disaggregation/disaggregation_help_functions.R")


# If SDALGCP not available for your version of R, reading dm_functions.R is sufficient
library(SDALGCP)
#source("functions/Other functions/SDALGCP/dm_functions.R")

setwd('functions/Other functions/ATA-Poisson-CoKriging-main/R')
source('ATA-ATP-Functions.R')

setwd(wd)



################################################################################
###################        Data specifications     #############################
################################################################################

# Read data
load("data/PBCshp.RData") # Liver cases
load('data/pop_den.RData') # Population density in raster of 300 x 300 meters
PBCshp_map = st_as_sf(PBCshp)

# Plot data
data_plot<-as(PBCshp, 'Spatial')
nclr <- 4        #### number of colors to be used
plotvar <-  PBCshp$X
plotclr <- rev(brewer.pal(nclr,"RdYlGn"))         #### define colors
class <- classIntervals(plotvar, nclr, style="fixed",fixedBreaks=c(0,1,3,5,Inf))    ### define categories
colcode <- findColours(class, plotclr)
plot(PBCshp, col=colcode, border=colcode)
legend('bottomleft',legend=names(attr(colcode, "table")),
       fill=attr(colcode, "palette"), cex=2, bty="n",y.intersp=1)



################################################################################
###################    Analysis specifications     #############################
################################################################################

cov.method <- "exponential"


# Find knots with space filling algorithm
knots = find_knots(PBCshp_map)

######### Grid for predictions
# Get the bounding box of the map
bbox <- st_bbox(PBCshp_map)

# Create a grid of points over the bounding box
grid_points <- expand.grid(
  x = seq(bbox["xmin"], bbox["xmax"], length.out = 200),
 y = seq(bbox["ymin"], bbox["ymax"], length.out = 200)
)


unioned_map = st_union(PBCshp_map)
outer_boundary = st_cast(unioned_map, "POLYGON")

grid_intersect = st_intersection(st_as_sf(grid_points[,1:2], coords = c("x","y"), crs = st_crs(outer_boundary)),
                                 outer_boundary)
grid<- as.data.frame(do.call(rbind, lapply(grid_intersect, st_coordinates)))
w1.0 <- grid[,1]
w2.0 <- grid[,2]



################################################################################
###################      Data analysis weighted    #############################
################################################################################

# Rasterize covariate IMD (same raster as population pop_den)
covariate_raster = rasterize(as(PBCshp_map,"Spatial"), pop_den, field = "IMD", 
                                         fun = mean, background = median(PBCshp_map$IMD))

# Find intersecting cells for each area
intersecting_cells_coordinates <- find_intersecting_cells(map = PBCshp_map, rstr_pop_density = pop_den,
                                                          rstr_covariates = list(linear = covariate_raster))

# Prepare data for method weighted SSDAM
dat_approx_weighted <- prepare_data(y = PBCshp_map$X, map = PBCshp_map,
                                   intersecting_cells = intersecting_cells_coordinates, # output from find_intersecting_cells
                                    exact = F, weighted = T)

# Area size population
true_pop <- exact_extract(pop_den,PBCshp_map, 'sum', progress = F)



# Initial values sigma
model <- glm(PBCshp_map$X ~PBCshp_map$IMD + offset(log(true_pop)), family = "poisson")
start_sigma = c(1/exp(1),min(mean(model$residuals^2),1))    



# Fit the model
mtime_area = proc.time()
model_approx_weighted <- Krig.fit(data = dat_approx_weighted, # Output from prepare_data
                             cov.method = cov.method,
                             knots = knots, # Output from find_knots
                             offset = true_pop, exact = F, n_init = 30, parallel = T,
                             sigma_init = start_sigma)
time_area = (proc.time()-mtime_area)[3]


################################################################################
################    Interpretation result weighted    ##########################
################################################################################

#### Plot spatial effect on new grid
spat_approx_weighted = predict_spatial(newdata = grid[,1:2], model = model_approx_weighted, plot = F,
                                  map = PBCshp_map)



dat1 = grid %>% mutate(Z = exp(spat_approx_weighted$fit))
bound <- raster::aggregate(PBCshp)
r1 <- raster::mask(raster::rasterFromXYZ(dat1), bound)
sp::spplot(r1, sp.layout=bound, at = seq(0, 3.8, length.out = 17),
           colorkey = list(labels = list(at = c(0,1,1.5,2,2.5,3,3.5)), space = "right"),
           as.table = TRUE)

#### Plot continuous incidence
pred_continuous_approx_weighted = predict_incidence_continuous(grid[,1:2], 
                                               model_approx_weighted, plot = F,
                                               rstr_covariates = list(linear = covariate_raster))




dat1_cont = grid %>% mutate(Z = exp(pred_continuous_approx_weighted$fit))
bound <- raster::aggregate(PBCshp)
r1_cont <- raster::mask(raster::rasterFromXYZ(dat1_cont), bound)
sp::spplot(r1_cont, sp.layout=bound, at = seq(0, 0.004, length.out = 17),
           colorkey = list(labels = list(at = c(0,0.001,0.002,0.003,0.004)), space = "right"),
           as.table = TRUE)


#### Estimated total count on area level
pred_inc_approx_weighted = predict_discrete_mean(PBCshp_map, model_approx_weighted, plot = F)


PBCshp_map$est_inc <- exp(pred_inc_approx_weighted$fit)/true_pop

#map the incidence
spplot(as(PBCshp_map,"Spatial"), c('est_inc'), 
       strip=lattice::strip.custom(factor.levels=c('est_inc')),
       at = seq(0, 0.0019, length.out = 17),
       colorkey = list(labels = list(at = c(0,0.0005,0.0010,0.0015)), space = "right"),
       as.table = TRUE)




################################################################################
##############      Data analysis weighted smooth     ##########################
################################################################################
covariate_raster = rasterize(as(PBCshp_map,"Spatial"), pop_den, field = "IMD", 
                             fun = mean, background = median(PBCshp_map$IMD))

# Smooth covariate
intersecting_cells_coordinates_smooth <- find_intersecting_cells(map = PBCshp_map, rstr_pop_density = pop_den,
                                                          rstr_covariates = list(smooth = covariate_raster))

dat_approx_weighted_smooth <- prepare_data(y = PBCshp_map$X, map = PBCshp_map,
                                    intersecting_cells = intersecting_cells_coordinates_smooth,
                                    exact = F, weighted = T)

IMD_smooth <- ps(PBCshp_map$IMD, df = 10)
model_smooth <- gam(PBCshp_map$X ~IMD_smooth + offset(log(true_pop)), family = "poisson",
                    paraPen = list(IMD_smooth = list(as.matrix(attributes(IMD_smooth)$S))))
start_sigma_smooth = c(1/exp(1),min(mean(model_smooth$residuals^2),1))    
start_pen_smooth = log(model_smooth$sp)


# Fit the model
mtime_area_smooth = proc.time()
model_approx_weighted_smooth <- Krig.fit(data = dat_approx_weighted_smooth,
                                  cov.method = cov.method,
                                  knots = knots,
                                  offset = true_pop, exact = F, n_init = 30, parallel = T,
                                  sigma_init = start_sigma_smooth, pen_init = start_pen_smooth, B = 10)
time_area_smooth = (proc.time()-mtime_area_smooth)[3]

#### Plot spatial effect
spat_approx_weighted_smooth = predict_spatial(newdata = grid[,1:2], model = model_approx_weighted_smooth, plot = F,
                                       limits = c(-1,0.7), map = PBCshp_map)



dat1_smooth = grid %>% mutate(Z = exp(spat_approx_weighted_smooth$fit))
bound <- raster::aggregate(PBCshp)
r1_smooth <- raster::mask(raster::rasterFromXYZ(dat1_smooth), bound)
sp::spplot(r1_smooth, sp.layout=bound, at = seq(0, 3.8, length.out = 17),
           colorkey = list(labels = list(at = c(0,1,1.5,2,2.5,3,3.5)), space = "right"),
           as.table = TRUE)


#### Plot continuous incidence
pred_continuous_approx_weighted_smooth = predict_incidence_continuous(grid[,1:2], 
                                                               model_approx_weighted_smooth, plot = F,
                                                               rstr_covariates = list(smooth = covariate_raster))


dat1_smooth_cont = grid %>% mutate(Z = exp(pred_continuous_approx_weighted_smooth$fit))
bound <- raster::aggregate(PBCshp)
r1_smooth_cont <- raster::mask(raster::rasterFromXYZ(dat1_smooth_cont), bound)
sp::spplot(r1_smooth_cont, sp.layout=bound, at = seq(0, 0.004, length.out = 17),
           colorkey = list(labels = list(at = c(0,0.001,0.002,0.003,0.004)), space = "right"),
           as.table = TRUE)



#### Plot estimated incidence
pred_inc_approx_weighted_smooth = predict_discrete_mean(PBCshp_map, model_approx_weighted_smooth, plot = F)


PBCshp_map$est_inc_smooth <- exp(pred_inc_approx_weighted_smooth$fit)/true_pop



#map the incidence
spplot(as(PBCshp_map,"Spatial"), c('est_inc_smooth'), 
       strip=lattice::strip.custom(factor.levels=c('est_inc_smooth')),
       at = seq(0, 0.0019, length.out = 17),
       colorkey = list(labels = list(at = c(0,0.0005,0.0010,0.0015)), space = "right"),
       as.table = TRUE)


#################  smooth component ###########
################################################
ngrid <- 100
xgrid <- seq(min(covariate_raster@data@values, na.rm=T), max(covariate_raster@data@values, na.rm=T), length.out = ngrid)

pred_smooth_weighted_RR <- predict_smooth_RR(xgrid = xgrid, model = model_approx_weighted_smooth, var = "smooth")

pdf("Figures/effect_IMD.pdf",height=6,width=8)
ggplot(pred_smooth_weighted_RR, aes(x = xgrid, y = fit)) +
  geom_line(color = "black") +  # Main effect line
  geom_line(aes(y = lower), color = "black", linetype = "dotted") +  # Lower CI
  geom_line(aes(y = upper), color = "black", linetype = "dotted") +  # Upper CI
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +  # Reference line at y=1
  labs(x ="IMD", y = "RR", title = "Effect of IMD on Outcome") +
  theme_minimal()
dev.off()


################################################################################
###################      Data analysis unweighted    #############################
################################################################################
covariate_raster = terra::rast(rasterize(as(PBCshp_map,"Spatial"), pop_den, field = "IMD", 
                                         fun = mean, background = median(PBCshp_map$IMD)))

intersecting_cells_coordinates <- find_intersecting_cells(map = PBCshp_map, rstr_pop_density = pop_den,
                                                          rstr_covariates = list(linear = covariate_raster))

dat_approx_unweighted <- prepare_data(y = PBCshp_map$X, map = PBCshp_map,
                                    intersecting_cells = intersecting_cells_coordinates,
                                    exact = F, weighted = F)

true_pop <- exact_extract(pop_den,PBCshp_map, 'sum', progress = F)


# Fit the model
mtime_area_unweighted = proc.time()
model_approx_unweighted <- Krig.fit(data = dat_approx_unweighted,
                                  cov.method = cov.method,
                                  knots = knots,
                                  offset = true_pop, exact = F, n_init = 30, parallel = T,
                                  sigma_init = start_sigma)
time_area_unweighted = (proc.time()-mtime_area_unweighted)[3]


################################################################################
################    Interpretation result unweighted    ########################
################################################################################

#### Plot spatial effect
spat_approx_unweighted = predict_spatial(newdata = grid[,1:2], model = model_approx_unweighted, plot = F,
                                       limits = c(-1,0.7), map = PBCshp_map)



dat2 = grid %>% mutate(Z = exp(spat_approx_unweighted$fit))
bound <- raster::aggregate(PBCshp)
r2 <- raster::mask(raster::rasterFromXYZ(dat2), bound)
sp::spplot(r2, sp.layout=bound, at = seq(0, 3.8, length.out = 17),
           colorkey = list(labels = list(at = c(0,1,1.5,2,2.5,3,3.5)), space = "right"),
           as.table = TRUE)

#### Plot continuous incidence
pred_continuous_approx_unweighted = predict_incidence_continuous(grid[,1:2], 
                                                                      model_approx_unweighted, plot = F,
                                                                      rstr_covariates = list(linear = covariate_raster))


dat2_cont = grid %>% mutate(Z = exp(pred_continuous_approx_unweighted$fit))
bound <- raster::aggregate(PBCshp)
r2_cont <- raster::mask(raster::rasterFromXYZ(dat2_cont), bound)
sp::spplot(r2_cont, sp.layout=bound, at = seq(0, 0.004, length.out = 17),
           colorkey = list(labels = list(at = c(0,0.001,0.002,0.003,0.004)), space = "right"),
           as.table = TRUE)


#### Plot estimated incidence
pred_inc_approx_unweighted = predict_discrete_mean(PBCshp_map, model_approx_unweighted, plot = F)
PBCshp_map$est_inc_unweighted <- exp(pred_inc_approx_unweighted$fit)/true_pop



#map the incidence
spplot(as(PBCshp_map,"Spatial"), c('est_inc_unweighted'), 
       strip=lattice::strip.custom(factor.levels=c('est_inc_unweighted')),
       at = seq(0, 0.0019, length.out = 17),
       colorkey = list(labels = list(at = c(0,0.0005,0.0010,0.0015)), space = "right"),
       as.table = TRUE)



################################################################################
###################        Data analysis exact     #############################
################################################################################
covariate_raster = terra::rast(rasterize(as(PBCshp_map,"Spatial"), pop_den, field = "IMD", 
                                         fun = mean, background = median(PBCshp_map$IMD)))

intersecting_cells_coordinates <- find_intersecting_cells(map = PBCshp_map, rstr_pop_density = pop_den,
                                                          rstr_covariates = list(linear = covariate_raster))

dat_exact <- prepare_data(y = PBCshp_map$X, map = PBCshp_map,
                                      intersecting_cells = intersecting_cells_coordinates,
                                      exact = T)

true_pop <- exact_extract(pop_den,PBCshp_map, 'sum', progress = F)


# Fit the model
mtime_area_exact = proc.time()
model_exact <- Krig.fit(data = dat_exact,
                                    cov.method = cov.method,
                                    knots = knots,
                                    offset = dat_exact$rstr_centers$pop, exact = T, n_init = 30, 
                                    parallel = T, sigma_init = start_sigma)
time_area_exact = (proc.time()-mtime_area_exact)[3]


################################################################################
################       Interpretation result exact      ########################
################################################################################

#### Plot spatial effect
spat_exact = predict_spatial(newdata = grid[,1:2], model = model_exact, plot = F,
                                         limits = c(-1,0.7), map = PBCshp_map)



dat3 = grid %>% mutate(Z = exp(spat_exact$fit))
bound <- raster::aggregate(PBCshp)
r3 <- raster::mask(raster::rasterFromXYZ(dat3), bound)
sp::spplot(r3, sp.layout=bound, at = seq(0, 3.8, length.out = 17),
           colorkey = list(labels = list(at = c(0,1,1.5,2,2.5,3,3.5)), space = "right"),
           as.table = TRUE)

#### Plot continuous incidence
pred_continuous_exact = predict_incidence_continuous(grid[,1:2], 
                                                      model_exact, plot = F,
                                                      rstr_covariates = list(linear = covariate_raster))


dat3_cont = grid %>% mutate(Z = exp(pred_continuous_exact$fit))
bound <- raster::aggregate(PBCshp)
r3_cont <- raster::mask(raster::rasterFromXYZ(dat3_cont), bound)
sp::spplot(r3_cont, sp.layout=bound, at = seq(0, 0.004, length.out = 17),
           colorkey = list(labels = list(at = c(0,0.001,0.002,0.003,0.004)), space = "right"),
           as.table = TRUE)

#### Plot estimated incidence
pred_inc_exact = predict_discrete_mean(PBCshp_map, model_exact, plot = F)
PBCshp_map$est_inc_exact <- exp(pred_inc_exact$fit)/true_pop



#map the incidence
spplot(as(PBCshp_map,"Spatial"), c('est_inc_exact'), 
       strip=lattice::strip.custom(factor.levels=c('est_inc_exact')),
       at = seq(0, 0.0019, length.out = 17),
       colorkey = list(labels = list(at = c(0,0.0005,0.0010,0.0015)), space = "right"),
       as.table = TRUE)





################################################################################
################### Comparing results with Johnson #############################
################################################################################


# BYM model
library(INLA)
library(spdep)

nb <- poly2nb(PBCshp_map, queen=T)
adj <- unlist(nb)
nb2INLA("adj",nb)
g <- inla.read.graph(filename="adj")

region <- 1:dim(PBCshp_map)[1]
formula <- y ~ xlin + f(region, model = "bym2", graph = g, param = c(2,0.5))
my.data.inla = data.frame(y = PBCshp_map$X)
my.data.inla$xlin = PBCshp_map$IMD
my.data.inla$area_pop = PBCshp_map$pop


mtime_BYM <- proc.time()
bym.model <- inla(formula, family = "Poisson", data = my.data.inla, offset = log(area_pop))
time_BYM <- (proc.time()-mtime_BYM)[3]



######## Predict
PBCshp_map$est_BYM <- as.numeric(as.numeric(bym.model$summary.fitted.values$mean)/PBCshp_map$pop)
spplot(as(PBCshp_map,"Spatial"), c('est_BYM'), 
       strip=lattice::strip.custom(factor.levels=c('est_BYM')),
       at = seq(0, 0.0021, length.out = 17),
       colorkey = list(labels = list(at = c(0,0.0005,0.0010,0.0015)), space = "right"),
       as.table = TRUE)



# SDALGCP model

#create the formula
data_lgcp <- data.frame(y = PBCshp_map$X, xlin = PBCshp_map$IMD, area_pop = PBCshp_map$pop)
FORM <- y ~ xlin +  offset(log(area_pop))

# create the set of phi's to use
phi <- seq(50, 2500, length.out = 100)

#set the mcmc parameters
control.mcmc <- controlmcmcSDA(n.sim = 1100000, burnin = 100000, thin= 100, h=1.65/(545^(1/6)),
                               c1.h = 0.01, c2.h = 1e-04)


start.time <- Sys.time()

s = pop_den
s[is.na(s[])] <- 0 
my_est <- SDALGCPMCML(data=data_lgcp, formula=FORM, my_shp=PBCshp, delta=300, phi=phi, method=1, pop_shp=s, 
                      weighted=TRUE,  plot=FALSE, par0=NULL, control.mcmc=control.mcmc)

# save(my_est, file="result_SDALGCP.RData")

######## Predict
Con_pred <- SDALGCPPred(para_est=my_est,  cellsize=300,  continuous=TRUE, pred.loc = grid, control.mcmc = control.mcmc)
# save(Con_pred, file="application/Johnson/models/result_pred_SDALGCP.RData")
end.time <- Sys.time()
time.taken1 <- end.time - start.time
time.taken1

summary(my_est)
confint(my_est)


PBCshp_map$est_SDALGCP <- Con_pred$my_shp$pMean_RR
PBCshp_map$sd_SDALGCP <- Con_pred$my_shp$pSD_RR

spplot(as(PBCshp_map,"Spatial"), c('est_SDALGCP'), 
       strip=lattice::strip.custom(factor.levels=c('est_SDALGCP')),
       at = seq(0, 0.0019, length.out = 17),
       colorkey = list(labels = list(at = c(0,0.0005,0.0010,0.0015)), space = "right"),
       as.table = TRUE)

dat4 <- data.frame(X=Con_pred$pred.loc[,1], Y=Con_pred$pred.loc[,2], Z = as.numeric(Con_pred$pred))
bound <- raster::aggregate(PBCshp)
r4 <- raster::mask(raster::rasterFromXYZ(dat4), bound)
sp::spplot(r4, sp.layout=bound, at = seq(0, 3.8, length.out = 17),
           colorkey = list(labels = list(at = c(0,1,1.5,2,2.5,3,3.5)), space = "right"),
           as.table = TRUE)



# Disaggregation package
library(disaggregation)
polygon_shapefile = PBCshp_map[,c("LSOA04CD", "X")]
polygon_shapefile <- polygon_shapefile%>% rename(ID = LSOA04CD)
covariate_raster = terra::rast(rasterize(as(PBCshp_map,"Spatial"), pop_den, field = "IMD", 
                                        fun = mean, background = median(PBCshp_map$IMD)))

dis_data = disaggregation::prepare_data(polygon_shapefile = polygon_shapefile,
                        covariate_rasters = covariate_raster,
                        aggregation_raster = terra::rast(pop_den),
                        id_var = "ID", response_var = "X", na_action = T)

fitted_disaggregation = disag_model(data = dis_data, iterations = 1000, family = "poisson", link = "log",
                                    hess_control_ndeps = 1e-10)

######## Predict

spatial_prediction = disaggr_predict_spatial(fitted_disaggregation, grid)

dat5 <- data.frame(X=grid[,1], Y=grid[,2], Z = exp(as.numeric(spatial_prediction)))
bound <- raster::aggregate(PBCshp)
r5 <- raster::mask(raster::rasterFromXYZ(dat5), bound)
sp::spplot(r5, sp.layout=bound, at = seq(0, 3.8, length.out = 17),
           colorkey = list(labels = list(at = c(0,1,1.5,2,2.5,3,3.5)), space = "right"),
           as.table = TRUE)


true_cov <- as.data.frame(cbind(grid, raster::extract(rstr_linear, grid)))
names(true_cov) <- c("x","y","xlin")

pred_continuous_disaggr = disaggr_pred_continuous_incidence(fitted_disaggregation,  
                  true_cov, predict_iid = T)

dat5_cont = grid %>% mutate(Z = exp(pred_continuous_disaggr))
bound <- raster::aggregate(PBCshp)
r5_cont <- raster::mask(raster::rasterFromXYZ(dat5_cont), bound)
sp::spplot(r5_cont, sp.layout=bound, at = seq(0, 0.004, length.out = 17),
           colorkey = list(labels = list(at = c(0,0.001,0.002,0.003,0.004)), space = "right"),
           as.table = TRUE)


#ATP

discretize.points = intersecting_cells_coordinates %>%
  rename(areaId = ID, ptx = x, pty = y, weight = tot_weight) %>%
  dplyr::select(areaId, ptx, pty, weight)

xy = intersecting_cells_coordinates %>% group_by(ID) %>%
  summarize(weightedpop(as.vector(x),
                        as.vector(y),
                        as.vector(tot_weight)))

data.test = data.frame(areaId = xy$ID, centx = xy$x, centy = xy$y,
                       counts = PBCshp_map$X, size = apply(matrix(PBCshp_map$pop), 1,int.pop)) 

createDataset <- function(data.test, discretePoints){
  rslt <- list(areaValues = data.test, discretePoints = discretePoints)
  class(rslt) <- c("list", "discreteArea")
  return(rslt)
}


data_atp = createDataset(data.test, discretize.points)

mtime_atp <- proc.time()
point_variogram <- deconvPointVgm(data_atp, model="Exp", rd = 1, fig = T)

predictions_atp <- atpKriging(data_atp, unknown = grid, point_variogram,
                              nmax = 5, showProgress = T)

time_atp <- (proc.time()-mtime_atp)[3]


dat6_cont = grid %>% mutate(Z = predictions_atp$pred) %>%
  mutate(Z = ifelse(Z<0,0,Z))
bound <- raster::aggregate(PBCshp)
r6_cont <- raster::mask(raster::rasterFromXYZ(dat6_cont), bound)
sp::spplot(r6_cont, sp.layout=bound, at = seq(0, 0.004, length.out = 17),
           colorkey = list(labels = list(at = c(0,0.001,0.002,0.003,0.004)), space = "right"),
           as.table = TRUE)




################### Plotting all models
# Map spatial effect
pdf("Figures/spatial_effect.pdf",width=5, height = 15)
rrr <- raster::stack(r3, r1, r2, r1_smooth, r4, r5)
sp::spplot(rrr, colorkey=list(space="right"), 
           names.attr = c('SSDEM', 'SSDAM I', 'SSDAM II', 'SSDAM I smooth', 'SDALGCP', 'disaggr.'),
           as.table = TRUE, layout = c(1,6), sp.layout=bound)
dev.off()

# Map continuous incidence
pdf("Figures/cont_incidence.pdf",width=5, height = 15)
rrr <- raster::stack(r3_cont, r1_cont, r2_cont, r1_smooth_cont, r5_cont, r6_cont)
sp::spplot(rrr, at = seq(0, 0.0019, length.out = 17),
           colorkey = list(labels = list(at = c(0,0.0005,0.0010,0.0015,0.0020)),space="right"),
           names.attr = c('SSDEM', 'SSDAM I', 'SSDAM II', 'SSDAM I smooth', 'disaggr.', 'atp'),
           as.table = TRUE, layout = c(1,6), sp.layout=bound)
dev.off()

# Map incidence
pdf("Figures/incidence_effect.pdf",width=5, height = 15)
spplot(as(PBCshp_map,"Spatial"), at = seq(0, 0.0019, length.out = 17),
       colorkey = list(labels = list(at = c(0,0.0005,0.0010,0.0015,0.0020)),space = "right"),
       c('est_inc_exact','est_inc','est_inc_unweighted','est_inc_smooth', 'est_SDALGCP', 'est_BYM'), 
       strip=lattice::strip.custom(factor.levels=c('SSDEM', 'SSDAM I', 'SSDAM II', 'SSDAM I smooth', 'SDALGCP', 'BYM')),
       layout = c(1,6), as.table = TRUE)
dev.off()

#################### Correlation measures
# Discrete
library(GGally)
cor_data_area <- data.frame("SSDEM" = PBCshp_map$est_inc_exact,
                            "SSDAM I" = PBCshp_map$est_inc,
                            "SSDAM II" = PBCshp_map$est_inc_unweighted,
                            "SSDAM I smooth" = PBCshp_map$est_inc_smooth,
                            "SDALGCP" = PBCshp_map$est_SDALGCP, 
                            "BYM" = PBCshp_map$est_BYM)

limitRange <- function(data, mapping, ...) { 
  ggplot(data = data, mapping = mapping, ...) + 
    geom_point(...) + 
    geom_abline(slope = 1, intercept = 0, linetype = 1, col = "red")+
    scale_y_continuous(limits = c(0, 0.0019)) +
    scale_x_continuous(limits = c(0, 0.0019)) 
}

pdf("Figures/incidence_cor.pdf", width = 12, height = 12)
ggpairs(cor_data_area ,columns = 1:6, 
        title = "Scatter Plot Matrix for area-level incidence", 
        axisLabels = "show", lower = list(continuous = limitRange)) 
dev.off()


# spatial
cor_data_continuous <- data.frame("SSDEM" = dat3$Z, "SSDAM I" = dat1$Z, "SSDAM II" = dat2$Z,
                                  "SSDAM I smooth" = dat1_smooth$Z,
                                  "SDALGCP" = dat4$Z, "disaggr." = dat5$Z)

limitRange <- function(data, mapping, ...) { 
  ggplot(data = data, mapping = mapping, ...) + 
    geom_point(...) + 
    geom_abline(slope = 1, intercept = 0, linetype = 1, col = "red")+
    scale_y_continuous(limits = c(0, 4)) +
    scale_x_continuous(limits = c(0, 4)) 
}

pdf("Figures/spatial_cor.pdf", width = 12, height = 12)
ggpairs(cor_data_continuous ,columns = 1:6, 
        title = "Scatter Plot Matrix for grid-level spatial effect", 
        axisLabels = "show", lower = list(continuous = limitRange)) 
dev.off()


# Continuous
cor_data_continuous_inc <- data.frame("SSDEM" = dat3_cont$Z, "SSDAM I" = dat1_cont$Z, "SSDAM II" = dat2_cont$Z,
                                  "SSDAM I smooth" = dat1_smooth_cont$Z, "disaggr." = dat5_cont$Z,
                                  "atp" = dat6_cont$Z)

limitRange <- function(data, mapping, ...) { 
  ggplot(data = data, mapping = mapping, ...) + 
    geom_point(...) + 
    geom_abline(slope = 1, intercept = 0, linetype = 1, col = "red")+
    scale_y_continuous(limits = c(0, 0.003)) +
    scale_x_continuous(limits = c(0, 0.003)) 
}

pdf("Figures/cont_inc_cor.pdf", width = 12, height = 12)
ggpairs(cor_data_continuous_inc ,columns = 1:6, 
        title = "Scatter Plot Matrix for grid-level incidence", 
        axisLabels = "show", lower = list(continuous = limitRange)) 
dev.off()

################# Compare covariance function ##################################
ind_cov = seq(0,40000, length = 500)

rho_approx_weighted = exp(model_approx_weighted$v_mode[length(model_approx_weighted$v_mode)])
sigma_approx_weighted  = 1/exp(model_approx_weighted$v_mode[length(model_approx_weighted$v_mode)-1])
cov_approx_weighted = sigma_approx_weighted * exp(-rho_approx_weighted* ind_cov)



rho_approx_unweighted = exp(model_approx_unweighted$v_mode[length(model_approx_unweighted$v_mode)])
sigma_approx_unweighted  = 1/exp(model_approx_unweighted$v_mode[length(model_approx_unweighted$v_mode)-1])
cov_approx_unweighted = sigma_approx_unweighted * exp(-rho_approx_unweighted* ind_cov)



rho_approx_weighted_smooth = exp(model_approx_weighted_smooth$v_mode[length(model_approx_weighted_smooth$v_mode)])
sigma_approx_weighted_smooth  = 1/exp(model_approx_weighted_smooth$v_mode[length(model_approx_weighted_smooth$v_mode)-1])
cov_approx_weighted_smooth = sigma_approx_weighted_smooth * exp(-rho_approx_weighted_smooth* ind_cov)




rho_exact = exp(model_exact$v_mode[length(model_exact$v_mode)])
sigma_exact  = 1/exp(model_exact$v_mode[length(model_exact$v_mode)-1])
cov_exact = sigma_exact * exp(-rho_exact* ind_cov)



range_disaggr = exp(fitted_disaggregation$opt$par["log_rho"])
sigma_disaggr = exp(fitted_disaggregation$opt$par["log_sigma"])^2
nu_disaggr = 1
rho_disaggr = sqrt(8)/range_disaggr #similar to rho_approx
cov_disaggr = sigma_disaggr * rho_disaggr * ind_cov * besselK(rho_disaggr * ind_cov, nu_disaggr)



range_sdalgcp = my_est$phi_opt
sigma_sdalgcp = my_est$sigma2_opt
rho_sdalgcp = 1/range_sdalgcp 
cov_sdalgcp = sigma_sdalgcp * exp(-rho_sdalgcp*ind_cov)



mean_inc = mean(PBCshp_map$X/true_pop)
sigma_atp <- log(point_variogram$pointVariogram$psill/mean_inc^2+1)
range_atp <- point_variogram$pointVariogram$range
rho_atp <- 1/range_atp
cov_atp <- sigma_atp * exp(-rho_atp * ind_cov)


data_cor <- data.frame(x = rep(ind_cov, 7),
                       y = c(cov_exact/sigma_exact, cov_approx_weighted/sigma_approx_weighted,
                             cov_approx_unweighted/sigma_approx_unweighted,
                             cov_approx_weighted_smooth/sigma_approx_weighted_smooth,
                             cov_disaggr/sigma_disaggr, cov_sdalgcp/sigma_sdalgcp,
                             cov_atp/sigma_atp),
                       Method = rep(c("SSDEM", "SSDAM I", "SSDAM II", "SSDAM I smooth",
                                  "disaggr.", "SDALGCP", "atp"), each = 500))

pdf("Figures/cov_uk.pdf", width = 8, height = 4)
ggplot(data = data_cor, aes(x = x, y = y, colour = Method, linetype = Method)) + 
  geom_line() +
  theme_bw()+
  labs(x = "distance", y = "cor")
dev.off()

################# Compare fixed effects ########################################


# SSDEM

est_ssdem <- model_exact$xi_estim[2]
sd_ssdem <- sqrt(model_exact$Sigma[2,2])


# Samples from posterior distribution
set.seed(123)
posterior_samples_ssdem <- rnorm(n = 1000000, est_ssdem, sd_ssdem)



# SSDAM I

est_ssdamI <- model_approx_weighted$xi_estim[2]
sd_ssdamI <- sqrt(model_approx_weighted$Sigma[2,2])

# Samples from posterior distribution
set.seed(123)
posterior_samples_ssdamI <- rnorm(n = 1000000, est_ssdamI, sd_ssdamI)



# SSDAM II

est_ssdamII <- model_approx_unweighted$xi_estim[2]
sd_ssdamII <- sqrt(model_approx_unweighted$Sigma[2,2])

# Samples from posterior distribution
set.seed(123)
posterior_samples_ssdamII <- rnorm(n = 1000000, est_ssdamII, sd_ssdamII)



# Disaggregation model
est_disaggr <- fitted_disaggregation$opt$par["layer"]
sd_disaggr <- sqrt(fitted_disaggregation$sd_out$cov.fixed[2,2])

# Samples from posterior distribution
set.seed(123)
posterior_samples_disaggr <- rnorm(n = 1000000, est_disaggr, sd_disaggr)


# SDALGCP model
est_lgcp <- my_est$beta_opt[2]
sd_lgcp <- sqrt(my_est$cov[2,2])

set.seed(123)
posterior_samples_sdalgcp <- rnorm(n = 1000000, est_lgcp, sd_lgcp)



# BYM model
est_bym <- bym.model$summary.results["xlin","Mean"]
posterior_samples_bym <- as.numeric(bym.model$samples$beta[,2])



pdf("Figures/IMD_effect.pdf", width = 8, height = 4)
ggplot() +
  geom_density(aes(x = posterior_samples_disaggr, col = "disaggr.")) +
  geom_density(aes(x = posterior_samples_ssdem, col = "SSDEM")) +
  geom_density(aes(x = posterior_samples_ssdamI, col = "SSDAM I")) +
  geom_density(aes(x = posterior_samples_ssdamII, col = "SSDAM II")) +
  geom_density(aes (x = posterior_samples_bym, col = "BYM")) +
  geom_density(aes(x = posterior_samples_sdalgcp, col = "SSDALGCP")) +
  labs(x = "Estimate IMD effect") +
  theme_bw()
dev.off()





# Compare to other correlation functions
model_approx_matern <- Krig.fit(data = dat_approx_weighted,
                                  cov.method = "matern",
                                  knots = knots,
                                  offset = true_pop, exact = F, n_init = 30, parallel = T,
                                  sigma_init = start_sigma)


rho_approx_matern = exp(model_approx_matern$v_mode[length(model_approx_matern$v_mode)])
sigma_approx_matern  = 1/exp(model_approx_matern$v_mode[length(model_approx_matern$v_mode)-1])
cov_approx_matern = exp(-rho_approx_matern * ind_cov)*(1+rho_approx_matern * ind_cov) * sigma_approx_matern



spat_approx_matern = predict_spatial(newdata = grid[,1:2], model = model_approx_matern, plot = F,
                                       limits = c(-1,0.7), map = PBCshp_map)



dat_matern = grid %>% mutate(Z = exp(spat_approx_matern$fit))
r_matern <- raster::mask(raster::rasterFromXYZ(dat_matern), bound)



model_approx_circular <- Krig.fit(data = dat_approx_weighted,
                                cov.method = "circular",
                                knots = knots,
                                offset = true_pop, exact = F, parallel = T,
                                sigma_init = start_sigma)

rho_approx_circular = exp(model_approx_circular$v_mode[length(model_approx_circular$v_mode)])
sigma_approx_circular  = 1/exp(model_approx_circular$v_mode[length(model_approx_circular$v_mode)-1])
cov_approx_circular = NULL

for (iter in 1:length(ind_cov)){
  theta = rho_approx_circular*ind_cov[iter]
  cov_approx_circular[iter] = sigma_approx_circular*ifelse(ind_cov[iter] <= 1/rho_approx_circular,
                                  1-2/pi*(theta*sqrt(1-theta^2)+asin(theta)),0)
  
}

spat_approx_circular = predict_spatial(newdata = grid[,1:2], model = model_approx_circular, plot = F,
                                     limits = c(-1,0.7), map = PBCshp_map)

dat_circular = grid %>% mutate(Z = exp(spat_approx_circular$fit))
r_circular <- raster::mask(raster::rasterFromXYZ(dat_circular), bound)


model_approx_spherical <- Krig.fit(data = dat_approx_weighted,
                                cov.method = "spherical",
                                knots = knots,
                                offset = true_pop, exact = F, parallel = T,
                                sigma_init = start_sigma)

rho_approx_spherical = exp(model_approx_spherical$v_mode[length(model_approx_spherical$v_mode)])
sigma_approx_spherical  = 1/exp(model_approx_spherical$v_mode[length(model_approx_spherical$v_mode)-1])
cov_approx_spherical = NULL

for (iter in 1:length(ind_cov)){
  cov_approx_spherical[iter] = sigma_approx_spherical*ifelse(ind_cov[iter]<= 1/rho_approx_spherical,
                                  1-1.5*rho_approx_spherical*ind_cov[iter]+
                                    0.5*rho_approx_spherical^3 * ind_cov[iter]^3,0)
}

spat_approx_spherical = predict_spatial(newdata = grid[,1:2], model = model_approx_spherical, plot = F,
                                     limits = c(-1,0.7), map = PBCshp_map)

dat_spherical = grid %>% mutate(Z = exp(spat_approx_spherical$fit))
r_spherical <- raster::mask(raster::rasterFromXYZ(dat_spherical), bound)


data_cor_dif <- data.frame(x = rep(ind_cov, 4),
                       y = c(cov_approx_weighted/sigma_approx_weighted,
                             cov_approx_matern/sigma_approx_matern,
                             cov_approx_circular/sigma_approx_circular,
                             cov_approx_spherical/sigma_approx_spherical),
                       Method = rep(c("exponential", "Matérn", "circular", "spherical"), each = 500))

pdf("Figures/cov_different.pdf", width = 8, height = 4)
ggplot(data = data_cor_dif, aes(x = x, y = y, colour = Method, linetype = Method)) + 
  geom_line() +
  theme_bw()+
  labs(x = "distance", y = "cor") +
  xlim(c(0,5000))
dev.off()


# Map spatial effect
pdf("Figures/spatial_effect_comparison.pdf",width=5, height = 11)
rrr <- raster::stack(r1, r_matern, r_circular, r_spherical)
sp::spplot(rrr, colorkey=list(space="right"), 
           names.attr = c('exponential', 'Matérn', 'circular', 'spherical'),
           as.table = TRUE, layout = c(1,4), sp.layout=bound)
dev.off()






# Compare to other correlation functions for exact method
model_exact_matern <- Krig.fit(data = dat_exact,
                                cov.method = "matern",
                                knots = knots,
                                offset = dat_exact$rstr_centers$pop, exact = T, n_init = 30, parallel = T,
                                sigma_init = start_sigma)


rho_exact_matern = exp(model_exact_matern$v_mode[length(model_exact_matern$v_mode)])
sigma_exact_matern  = 1/exp(model_exact_matern$v_mode[length(model_exact_matern$v_mode)-1])
cov_exact_matern = exp(-rho_exact_matern * ind_cov)*(1+rho_exact_matern * ind_cov) * sigma_exact_matern



spat_exact_matern = predict_spatial(newdata = grid[,1:2], model = model_exact_matern, plot = F,
                                     limits = c(-1,0.7), map = PBCshp_map)



dat3_matern = grid %>% mutate(Z = exp(spat_exact_matern$fit))
r3_matern <- raster::mask(raster::rasterFromXYZ(dat3_matern), bound)



model_exact_circular <- Krig.fit(data = dat_exact,
                                  cov.method = "circular",
                                  knots = knots,
                                  offset = dat_exact$rstr_centers$pop, exact = T, parallel = T,
                                  sigma_init = start_sigma)

rho_exact_circular = exp(model_exact_circular$v_mode[length(model_exact_circular$v_mode)])
sigma_exact_circular  = 1/exp(model_exact_circular$v_mode[length(model_exact_circular$v_mode)-1])
cov_exact_circular = NULL

for (iter in 1:length(ind_cov)){
  theta = rho_exact_circular*ind_cov[iter]
  cov_exact_circular[iter] = sigma_exact_circular*ifelse(ind_cov[iter] <= 1/rho_exact_circular,
                                                           1-2/pi*(theta*sqrt(1-theta^2)+asin(theta)),0)
  
}

spat_exact_circular = predict_spatial(newdata = grid[,1:2], model = model_exact_circular, plot = F,
                                       limits = c(-1,0.7), map = PBCshp_map)

dat3_circular = grid %>% mutate(Z = exp(spat_exact_circular$fit))
r3_circular <- raster::mask(raster::rasterFromXYZ(dat3_circular), bound)


model_exact_spherical <- Krig.fit(data = dat_exact,
                                   cov.method = "spherical",
                                   knots = knots,
                                   offset = dat_exact$rstr_centers$pop, exact = T, parallel = T,
                                   sigma_init = start_sigma)
rho_exact_spherical = exp(model_exact_spherical$v_mode[length(model_exact_spherical$v_mode)])
sigma_exact_spherical  = 1/exp(model_exact_spherical$v_mode[length(model_exact_spherical$v_mode)-1])
cov_exact_spherical = NULL

for (iter in 1:length(ind_cov)){
  cov_exact_spherical[iter] = sigma_exact_spherical*ifelse(ind_cov[iter]<= 1/rho_exact_spherical,
                                                             1-1.5*rho_exact_spherical*ind_cov[iter]+
                                                               0.5*rho_exact_spherical^3 * ind_cov[iter]^3,0)
}

spat_exact_spherical = predict_spatial(newdata = grid[,1:2], model = model_exact_spherical, plot = F,
                                        limits = c(-1,0.7), map = PBCshp_map)

dat3_spherical = grid %>% mutate(Z = exp(spat_exact_spherical$fit))
r3_spherical <- raster::mask(raster::rasterFromXYZ(dat3_spherical), bound)


data_cor_dif_exact <- data.frame(x = rep(ind_cov, 4),
                           y = c(cov_exact/sigma_exact,
                                 cov_exact_matern/sigma_exact_matern,
                                 cov_exact_circular/sigma_exact_circular,
                                 cov_exact_spherical/sigma_exact_spherical),
                           Method = rep(c("exponential", "Matérn", "circular", "spherical"), each = 500))

pdf("Figures/cov_different_exact.pdf", width = 8, height = 4)
ggplot(data = data_cor_dif_exact, aes(x = x, y = y, colour = Method, linetype = Method)) + 
  geom_line() +
  theme_bw()+
  labs(x = "distance", y = "cor") +
  xlim(c(0,5000))
dev.off()


# Map spatial effect
pdf("Figures/spatial_effect_comparison_exact.pdf",width=5, height = 11)
rrr_exact <- raster::stack(r3, r3_matern, r3_circular, r3_spherical)
sp::spplot(rrr_exact, colorkey=list(space="right"), 
           names.attr = c('exponential', 'Matérn', 'circular', 'spherical'),
           as.table = TRUE, layout = c(1,4), sp.layout=bound)
dev.off()








