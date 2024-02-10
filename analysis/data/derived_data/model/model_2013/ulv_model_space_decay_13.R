## ----libraries---------------------------------------------------------------

library(pacman)
pacman::p_load(INLA, tidyverse, sf, lubridate, parallel, here)


## ----data--------------------------------------------------------------------------------------------------

data <- readRDS(here("analysis", "data", "derived_data", "variables", "variables_2013", "vbles.space.decay.13.rds"))

## ----coordinates------------------------------------------------------------------------------------------------

# extract coordinates from data that will be used as vertices for the mesh

data.coord <- cbind(
  x = unlist(map(data$geometry, 1)),
  y = unlist(map(data$geometry, 2))
)

## ----create INLA mesh--------------------------------------------------------------------------------------------

# creates a polygon around the coordinates to allow smaller triangles inside polygon and larger triangles outside
# this saves computational power in the area with larger triangles
bnd <- inla.nonconvex.hull(data.coord, concave = 200)

#  mesh that allows discretization of the random field
mesh <- inla.mesh.2d(data.coord, boundary = bnd, max.edge = c(50, 100), cutoff = 20, min.angle = 30)
# plot of the mesh for this sample data
plot(mesh, asp = 1, main = "")
points(data.coord, pch = 16, col = "red", cex = 0.3, lwd = 1)


## ----projection A matrix---------------------------------------------------------------------------

# projection matrix A that projects the continuous Gaussian random field from the data points to the mesh nodes.
a1 <- inla.spde.make.A(mesh = mesh, loc = data.coord)

## ----SPDE with PC priors---------------------------------------------------------------------------

# builds the SPDE model
spde <- inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(10, 0.01), # the probability that the range is less than 10 (m) is 0.01
  prior.sigma = c(1, 0.01)
) # the probability that variance (on the log scale) is more that 1 is 0.01)

# make the index set for the SPDE model
sp.index <- inla.spde.make.index(
  name = "spatial.field", # name of the effect
  n.spde = spde$n.spde
) # number of vertices of the SPDE model



## ----models including a b variable to measure the effect of sprays in surrounding households-------------------------------------------------------------------------

#' Run INLA models including a variable that measures the effect of sprays in households surrounding household i
#'
#' This function takes a vector of variables and a data set as input and returns a list of INLA models where each model includes 
#' one variable that measures the effect of sprays in households surrounding household i
#' and one variable that was selected as best representing the effect of sprays within household i (g20_day.max)
#'
#' @param df A data set.
#' @param vble.list A vector of variables
#' @param ncores A numeric value of the number of computing cores to be used for this operation
#' @return A list of INLA models
#'

fx.model.2fe<- function(df, vble.list, ncores){
  fx.model<- function(df, vble.list){
  stack <- inla.stack(
  data= list(y= list(df[["AA_TOTAL"]])), 
  A = list(a1,1), 
  effects= list(c(sp.index, list(Intercept= 1)),
                # vble = one variable from a list of candidate variables that measure the effect of sprays in households surrounding household i
                list(vble= df[[vble.list]], 
                    # g20_day.max= the variable from the best-fitting model from the within-household model comparison
                    # this accounts for sprays within household i
                  g20_day.max= df[["g20_day.max"]],
                  location_code = df[["location_code"]], 
                     month= df[["month"]])),
  tag= "est")

  
formula <- y ~ -1 + Intercept + g20_day.max + vble + f(location_code, model ="iid")+ f(month, model ="rw1")+ f(spatial.field, model = spde)

m <- inla(formula, 
           data = inla.stack.data(stack, 
                                  spde=spde), 
           family= "nbinomial",
           control.family=list(link='log'),
           control.compute = list(cpo= TRUE, dic= TRUE, waic=TRUE), 
           control.predictor = list(link=1, A= inla.stack.A(stack), 
                                    compute=T)
           )
      c(m, fixed.effect.s= vble.list)
}

mclapply(X= vble.list  , FUN= fx.model, df= df , mc.cores = ncores)

}

## ----run surrounding house models--------------------------------------------------------------------------------------

cores= 2

# vector of variables that measures the effect of sprays in households surrounding household i
l.variables.decay <- c("i_total",
                       "g5_total", "g25_total", "g50_total", "g75_total", "g100_total",
                       "g125_total", "g150_total", "g200_total", "g250_total", "g300_total",
                       "e0.0025_total", "e0.0035_total", "e0.005_total", "e0.0075_total",
                       "e0.01_total", "e0.0125_total", "e0.02_total", "e0.045_total", "e0.2_total")
  

m.sur.list.2013 <- fx.model.2fe(df= data, vble.list= l.variables.decay, ncores= cores)

saveRDS(m.sur.list.2013 , here("analysis", "data", "derived_data", "model", "model_2013", "model_outputs_2013", "m.space.decay.13.rds"))

