## ----libraries---------------------------------------------------------------

library(pacman)
pacman::p_load(INLA, tidyverse, sf, lubridate, parallel, here)


## ----data--------------------------------------------------------------------------------------------------

data <- readRDS(here("analysis", "data", "derived_data", "variables", "variables_2013", "vbles.time.13.rds"))

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
#plot(mesh, asp = 1, main = "")
#points(data.coord, pch = 16, col = "red", cex = 0.3, lwd = 1)


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



## ----m1.1 month--------------------------------------------------------------------------------------------

# estimation INLA data stack
stack <- inla.stack(
  # response variable number of Ae. aegypti
  data = list(y = list(data[["AA_TOTAL"]])),

  # the projector matrix (for the spatial effect) and a linear vector (1) for the other effects
  # the effects are organised in a list of lists.
  A = list(a1, 1),

  # spatial effect that accounts for spatial autocorrelation and the intercept
  effects = list(
    c(sp.index, list(Intercept = 1)),
    # and then all the other effects, in this case month
    list(month = data[["month"]])
  ),
  tag = "est"
)

# f(month, model ="rw1") is a random walk of order one (RW1) effect for month
formula <- y ~ -1 + Intercept + f(month, model = "rw1") +
  # f(spatial.field, model = spde) accounts for spatial autocorrelation using SPDE
  f(spatial.field, model = spde)

m1.1 <- inla(formula,
  data = inla.stack.data(stack,
    spde = spde
  ),
  # negative binomial distribution
  family = "nbinomial",
  control.family = list(link = "log"),
  # add WAIC and other INLA result outputs
  control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE),
  control.predictor = list(
    link = 1, A = inla.stack.A(stack),
    compute = T
  )
)


## ----m1.2 household------------------------------------------------------------------------------------

stack <- inla.stack(
  # response variable number of Ae. aegypti
  data = list(y = list(data[["AA_TOTAL"]])),

  # the projector matrix (for the spatial effect) and a linear vector (1) for the other effects
  # the effects are organised in a list of lists.
  A = list(a1, 1),

  # spatial effect that accounts for spatial autocorrelation and the intercept
  effects = list(
    c(sp.index, list(Intercept = 1)),
    # and then all the other effects, in this case household
    list(location_code = data[["location_code"]])
  ),
  tag = "est"
)

# f(location_code, model ="iid") is an independent, identically distributed (iid) random effect that accounts for baseline abundance each i
formula <- y ~ -1 + Intercept + f(location_code, model = "iid") +
  # f(spatial.field, model = spde) accounts for spatial autocorrelation using SPDE
  f(spatial.field, model = spde)
# the spatial effect is specified using the spde tag (which is why we don't use the "" for it)

m1.2 <- inla(formula,
  data = inla.stack.data(stack,
    spde = spde
  ),
  # negative binomial distribution
  family = "nbinomial",
  control.family = list(link = "log"),
  # add WAIC and other INLA result outputs
  control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE),
  control.predictor = list(
    link = 1, A = inla.stack.A(stack),
    compute = T
  )
)


## ----m1.3 month + household----------------------------------------------------------------------------

stack <- inla.stack(
  # response variable number of Ae. aegypti
  data = list(y = list(data[["AA_TOTAL"]])),

  # the projector matrix (for the spatial effect) and a linear vector (1) for the other effects
  # the effects are organised in a list of lists.
  A = list(a1, 1),

  # spatial effect that accounts for spatial autocorrelation and the intercept
  effects = list(
    c(sp.index, list(Intercept = 1)),
    # and then all the other effects, in this case household and month
    list(
      location_code = data[["location_code"]],
      month = data[["month"]]
    )
  ),
  tag = "est"
)

# this model includes all three random effects: spatial autocorrelation, random walk 1 by month, and repeated measures by household
formula <- y ~ -1 + Intercept + f(location_code, model = "iid") + f(month, model = "rw1") +
  f(spatial.field, model = spde)
# the spatial effect is specified using the spde tag (which is why we don't use the "" for it)

m1.3 <- inla(formula,
  data = inla.stack.data(stack,
    spde = spde
  ),
  family = "nbinomial",
  # negative binomial distribution
  control.family = list(link = "log"),
  # add WAIC and other INLA result outputs
  control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE),
  control.predictor = list(
    link = 1, A = inla.stack.A(stack),
    compute = T
  )
)


## ----m.random----------------------------------------------------------------------------------------------

m.random <- list(m1.1, m1.2, m1.3)

saveRDS(m.random, here("analysis", "data", "derived_data", "model", "model_2013", "model_outputs_2013", "m.random.13.rds"))
