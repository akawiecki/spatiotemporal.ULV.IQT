
library(pacman)
pacman::p_load(INLA, tidyverse, sf, lubridate, parallel, here)


## ----data--------------------------------------------------------------------------------------------------------------

data.bio <-  readRDS(here("analysis", "data", "derived_data", "variables", "variables_2014", "vbles.space.ring.bio.14.rds"))%>% 
  mutate(month= as.factor(month))

data.100m <-  readRDS(here("analysis", "data", "derived_data", "variables", "variables_2014", "vbles.space.ring.100m.14.rds"))%>% 
  mutate(month= as.factor(month))

## ----coordinates------------------------------------------------------------------------------------------------

# extract coordinates from data that will be used as vertices for the mesh
# here we use data.bio but we could have used data.100m interchangeably as both data sets have the same household location points
data.coord <- cbind(
  x = unlist(map(data.bio$geometry, 1)),
  y = unlist(map(data.bio$geometry, 2))
)

## ----create INLA mesh--------------------------------------------------------------------------------------------

# creates a polygon around the coordinates to allow smaller triangles inside polygon and larger triangles outside
# this saves computational power in the area with larger triangles
bnd <- inla.nonconvex.hull(data.coord, concave = 200)

#  mesh that allows discretization of the random field
mesh <- inla.mesh.2d(data.coord, boundary = bnd, max.edge = c(50, 100), cutoff = 20, min.angle = 30)


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

# run models with proportions of households within a ring of a given distance from household i that had been sprayed in the week before t 

## ----biological plausibility rings------------------------------------------------------------------------------------------------------------------------


stack <- inla.stack(
  data= list(y= list(data.bio[["AA_TOTAL"]])), 
  #response variable
  A = list(a1,1), 
  # the projector matrix (for the spatial effect) and a linear vector (1) for the other effects
  # the effects are organised in a list of lists. 
  effects= list(c(sp.index, list(Intercept= 1)),
                list(p30m= data.bio[["perd.s.ring_0-30m"]], 
                     p100m= data.bio[["perd.s.ring_31-100m"]],
                     p300m= data.bio[["perd.s.ring_101-300m"]],
                     ptotalm= data.bio[["perd.s.ring_>300m"]],
                     g20_day.max= data.bio[["g20_day.max"]],
                     location_code = data.bio[["location_code"]], 
                     month= data.bio[["month"]])),
  #and then all the other effects
  tag= "est")
# The tag specify the name of this stack
#and then all the other effects

formula <- y ~ -1 + Intercept + g20_day.max + p30m + p100m+ p300m + ptotalm+ 
  f(location_code, model ="iid")+ f(month, model ="rw1")+ f(spatial.field, model = spde) 
# the spatial effect is specified using the spde tag (which is why we don't use the "" for it)

m.all.perd.bio <- c(inla(formula, 
                         data = inla.stack.data(stack, 
                                                spde=spde), 
                         family= "nbinomial",
                         control.family=list(link='log'),
                         control.compute = list(cpo= TRUE, dic= TRUE, waic=TRUE), 
                         control.predictor = list(link=1, A= inla.stack.A(stack), 
                                                  compute=T)
),fixed.effect.s= "g20_day.max + p30m + p100m + p300m + ptotalm")

stack <- inla.stack(
  data= list(y= list(data.bio[["AA_TOTAL"]])), 
  #response variable
  A = list(a1,1), 
  # the projector matrix (for the spatial effect) and a linear vector (1) for the other effects
  # the effects are organised in a list of lists. 
  effects= list(c(sp.index, list(Intercept= 1)),
                list(p30m= data.bio[["perd.s.ring_0-30m"]], 
                     p100m= data.bio[["perd.s.ring_31-100m"]],
                     p300m= data.bio[["perd.s.ring_101-300m"]],
                     g20_day.max= data.bio[["g20_day.max"]],
                     location_code = data.bio[["location_code"]], 
                     month= data.bio[["month"]])),
  #and then all the other effects
  tag= "est")
# The tag specify the name of this stack
#and then all the other effects

formula <- y ~ -1 + Intercept + g20_day.max + p30m + p100m+  p300m+ 
  f(location_code, model ="iid")+ f(month, model ="rw1")+ f(spatial.field, model = spde) 
# the spatial effect is specified using the spde tag (which is why we don't use the "" for it)

m.300.perd.bio <- c(inla(formula, 
                         data = inla.stack.data(stack, 
                                                spde=spde), 
                         family= "nbinomial",
                         control.family=list(link='log'),
                         control.compute = list(cpo= TRUE, dic= TRUE, waic=TRUE), 
                         control.predictor = list(link=1, A= inla.stack.A(stack), 
                                                  compute=T)
),fixed.effect.s= "g20_day.max + p30m + p100m + p300m")



stack <- inla.stack(
  data= list(y= list(data.bio[["AA_TOTAL"]])), 
  #response variable
  A = list(a1,1), 
  # the projector matrix (for the spatial effect) and a linear vector (1) for the other effects
  # the effects are organised in a list of lists. 
  effects= list(c(sp.index, list(Intercept= 1)),
                list(p30m= data.bio[["perd.s.ring_0-30m"]], 
                     p100m= data.bio[["perd.s.ring_31-100m"]],
                     g20_day.max= data.bio[["g20_day.max"]],
                     location_code = data.bio[["location_code"]], 
                     month= data.bio[["month"]])),
  #and then all the other effects
  tag= "est")
# The tag specify the name of this stack
#and then all the other effects

formula <- y ~ -1 + Intercept + g20_day.max + p30m + p100m + 
  f(location_code, model ="iid")+ f(month, model ="rw1")+ f(spatial.field, model = spde) 
# the spatial effect is specified using the spde tag (which is why we don't use the "" for it)

m.100.perd.bio <- c(inla(formula, 
                         data = inla.stack.data(stack, 
                                                spde=spde), 
                         family= "nbinomial",
                         control.family=list(link='log'),
                         control.compute = list(cpo= TRUE, dic= TRUE, waic=TRUE), 
                         control.predictor = list(link=1, A= inla.stack.A(stack), 
                                                  compute=T)
),fixed.effect.s= "g20_day.max + p30m + p100m")

stack <- inla.stack(
  data= list(y= list(data.bio[["AA_TOTAL"]])), 
  #response variable
  A = list(a1,1), 
  # the projector matrix (for the spatial effect) and a linear vector (1) for the other effects
  # the effects are organised in a list of lists. 
  effects= list(c(sp.index, list(Intercept= 1)),
                list(p30m= data.bio[["perd.s.ring_0-30m"]], 
                     g20_day.max= data.bio[["g20_day.max"]],
                     location_code = data.bio[["location_code"]], 
                     month= data.bio[["month"]])),
  #and then all the other effects
  tag= "est")
# The tag specify the name of this stack
#and then all the other effects

formula <- y ~ -1 + Intercept + g20_day.max + p30m +
  f(location_code, model ="iid")+ f(month, model ="rw1")+ f(spatial.field, model = spde) 
# the spatial effect is specified using the spde tag (which is why we don't use the "" for it)

m.30.perd.bio <- c(inla(formula, 
                        data = inla.stack.data(stack, 
                                               spde=spde), 
                        family= "nbinomial",
                        control.family=list(link='log'),
                        control.compute = list(cpo= TRUE, dic= TRUE, waic=TRUE), 
                        control.predictor = list(link=1, A= inla.stack.A(stack), 
                                                 compute=T)
),fixed.effect.s= "g20_day.max + p30m")



m.rings.bio <- list(m.all.perd.bio, m.300.perd.bio, m.100.perd.bio,m.30.perd.bio)

## ----model all perd.100m------------------------------------------------------------------------------------------------------------------------

stack <- inla.stack(
  data= list(y= list(data.100m[["AA_TOTAL"]])), 
  #response variable
  A = list(a1,1), 
  # the projector matrix (for the spatial effect) and a linear vector (1) for the other effects
  # the effects are organised in a list of lists. 
  effects= list(c(sp.index, list(Intercept= 1)),
                list(
                  p100m= data.100m[["perd.s.ring_0-100m"]],
                  p200m= data.100m[["perd.s.ring_101-200m"]],
                  p300m= data.100m[["perd.s.ring_201-300m"]],
                  p400m= data.100m[["perd.s.ring_301-400m"]],
                  p500m= data.100m[["perd.s.ring_401-500m"]],
                  ptotalm= data.100m[["perd.s.ring_>500m"]],
                  g20_day.max= data.100m[["g20_day.max"]],
                  location_code = data.100m[["location_code"]], 
                  month= data.100m[["month"]])),
  #and then all the other effects
  tag= "est")
# The tag specify the name of this stack
#and then all the other effects

formula <- y ~ -1 + Intercept + g20_day.max + p100m + p200m + p300m + p400m + p500m +ptotalm+ 
  f(location_code, model ="iid")+ f(month, model ="rw1")+ f(spatial.field, model = spde) 
# the spatial effect is specified using the spde tag (which is why we don't use the "" for it)

m.all.perd.100m <- c(inla(formula, 
                     data = inla.stack.data(stack, 
                                            spde=spde), 
                     family= "nbinomial",
                     control.family=list(link='log'),
                     control.compute = list(cpo= TRUE, dic= TRUE, waic=TRUE), 
                     control.predictor = list(link=1, A= inla.stack.A(stack), 
                                              compute=T)
),fixed.effect.s= "g20_day.max + p100m + p200m + p300m + p400m + p500m +ptotalm")


stack <- inla.stack(
  data= list(y= list(data.100m[["AA_TOTAL"]])), 
  #response variable
  A = list(a1,1), 
  # the projector matrix (for the spatial effect) and a linear vector (1) for the other effects
  # the effects are organised in a list of lists. 
  effects= list(c(sp.index, list(Intercept= 1)),
                list(
                  p100m= data.100m[["perd.s.ring_0-100m"]],
                  p200m= data.100m[["perd.s.ring_101-200m"]],
                  p300m= data.100m[["perd.s.ring_201-300m"]],
                  p400m= data.100m[["perd.s.ring_301-400m"]],
                  p500m= data.100m[["perd.s.ring_401-500m"]],
                  g20_day.max= data.100m[["g20_day.max"]],
                  location_code = data.100m[["location_code"]], 
                  month= data.100m[["month"]])),
  #and then all the other effects
  tag= "est")
# The tag specify the name of this stack
#and then all the other effects

formula <- y ~ -1 + Intercept + g20_day.max + p100m + p200m + p300m + p400m + p500m +
  f(location_code, model ="iid")+ f(month, model ="rw1")+ f(spatial.field, model = spde) 
# the spatial effect is specified using the spde tag (which is why we don't use the "" for it)

m.500m.perd.100m <- c(inla(formula, 
                 data = inla.stack.data(stack, 
                                        spde=spde), 
                 family= "nbinomial",
                 control.family=list(link='log'),
                 control.compute = list(cpo= TRUE, dic= TRUE, waic=TRUE), 
                 control.predictor = list(link=1, A= inla.stack.A(stack), 
                                          compute=T)
),fixed.effect.s= "g20_day.max + p100m + p200m + p300m + p400m + p500m")


stack <- inla.stack(
  data= list(y= list(data.100m[["AA_TOTAL"]])), 
  #response variable
  A = list(a1,1), 
  # the projector matrix (for the spatial effect) and a linear vector (1) for the other effects
  # the effects are organised in a list of lists. 
  effects= list(c(sp.index, list(Intercept= 1)),
                list(
                  p100m= data.100m[["perd.s.ring_0-100m"]],
                  p200m= data.100m[["perd.s.ring_101-200m"]],
                  p300m= data.100m[["perd.s.ring_201-300m"]],
                  p400m= data.100m[["perd.s.ring_301-400m"]],
                  g20_day.max= data.100m[["g20_day.max"]],
                  location_code = data.100m[["location_code"]], 
                  month= data.100m[["month"]])),
  #and then all the other effects
  tag= "est")
# The tag specify the name of this stack
#and then all the other effects

formula <- y ~ -1 + Intercept + g20_day.max + p100m + p200m + p300m + p400m +
  f(location_code, model ="iid")+ f(month, model ="rw1")+ f(spatial.field, model = spde) 
# the spatial effect is specified using the spde tag (which is why we don't use the "" for it)

m.400m.perd.100m <- c(inla(formula, 
                 data = inla.stack.data(stack, 
                                        spde=spde), 
                 family= "nbinomial",
                 control.family=list(link='log'),
                 control.compute = list(cpo= TRUE, dic= TRUE, waic=TRUE), 
                 control.predictor = list(link=1, A= inla.stack.A(stack), 
                                          compute=T)
),fixed.effect.s= "g20_day.max + p100m + p200m + p300m + p400m")


stack <- inla.stack(
  data= list(y= list(data.100m[["AA_TOTAL"]])), 
  #response variable
  A = list(a1,1), 
  # the projector matrix (for the spatial effect) and a linear vector (1) for the other effects
  # the effects are organised in a list of lists. 
  effects= list(c(sp.index, list(Intercept= 1)),
                list(
                  p100m= data.100m[["perd.s.ring_0-100m"]],
                  p200m= data.100m[["perd.s.ring_101-200m"]],
                  p300m= data.100m[["perd.s.ring_201-300m"]],
                  g20_day.max= data.100m[["g20_day.max"]],
                  location_code = data.100m[["location_code"]], 
                  month= data.100m[["month"]])),
  #and then all the other effects
  tag= "est")
# The tag specify the name of this stack
#and then all the other effects

formula <- y ~ -1 + Intercept + g20_day.max + p100m + p200m + p300m +
  f(location_code, model ="iid")+ f(month, model ="rw1")+ f(spatial.field, model = spde) 
# the spatial effect is specified using the spde tag (which is why we don't use the "" for it)

m.300m.perd.100m <- c(inla(formula, 
                 data = inla.stack.data(stack, 
                                        spde=spde), 
                 family= "nbinomial",
                 control.family=list(link='log'),
                 control.compute = list(cpo= TRUE, dic= TRUE, waic=TRUE), 
                 control.predictor = list(link=1, A= inla.stack.A(stack), 
                                          compute=T)
),fixed.effect.s= "g20_day.max + p100m + p200m + p300m")


stack <- inla.stack(
  data= list(y= list(data.100m[["AA_TOTAL"]])), 
  #response variable
  A = list(a1,1), 
  # the projector matrix (for the spatial effect) and a linear vector (1) for the other effects
  # the effects are organised in a list of lists. 
  effects= list(c(sp.index, list(Intercept= 1)),
                list(
                  p100m= data.100m[["perd.s.ring_0-100m"]],
                  p200m= data.100m[["perd.s.ring_101-200m"]],
                  g20_day.max= data.100m[["g20_day.max"]],
                  location_code = data.100m[["location_code"]], 
                  month= data.100m[["month"]])),
  #and then all the other effects
  tag= "est")
# The tag specify the name of this stack
#and then all the other effects

formula <- y ~ -1 + Intercept + g20_day.max + p100m + p200m +
  f(location_code, model ="iid")+ f(month, model ="rw1")+ f(spatial.field, model = spde) 
# the spatial effect is specified using the spde tag (which is why we don't use the "" for it)

m.200m.perd.100m <- c(inla(formula, 
                 data = inla.stack.data(stack, 
                                        spde=spde), 
                 family= "nbinomial",
                 control.family=list(link='log'),
                 control.compute = list(cpo= TRUE, dic= TRUE, waic=TRUE), 
                 control.predictor = list(link=1, A= inla.stack.A(stack), 
                                          compute=T)
),fixed.effect.s= "g20_day.max + p100m + p200m")


stack <- inla.stack(
  data= list(y= list(data.100m[["AA_TOTAL"]])), 
  #response variable
  A = list(a1,1), 
  # the projector matrix (for the spatial effect) and a linear vector (1) for the other effects
  # the effects are organised in a list of lists. 
  effects= list(c(sp.index, list(Intercept= 1)),
                list(
                  p100m= data.100m[["perd.s.ring_0-100m"]],
                  g20_day.max= data.100m[["g20_day.max"]],
                  location_code = data.100m[["location_code"]], 
                  month= data.100m[["month"]])),
  #and then all the other effects
  tag= "est")
# The tag specify the name of this stack
#and then all the other effects

formula <- y ~ -1 + Intercept + g20_day.max + p100m +
  f(location_code, model ="iid")+ f(month, model ="rw1")+ f(spatial.field, model = spde) 
# the spatial effect is specified using the spde tag (which is why we don't use the "" for it)

m.100m.perd.100m <- c(inla(formula, 
                 data = inla.stack.data(stack, 
                                        spde=spde), 
                 family= "nbinomial",
                 control.family=list(link='log'),
                 control.compute = list(cpo= TRUE, dic= TRUE, waic=TRUE), 
                 control.predictor = list(link=1, A= inla.stack.A(stack), 
                                          compute=T)
),fixed.effect.s= "g20_day.max + p100m")


m.rings.100m <- list(m.all.perd.100m, m.500m.perd.100m , m.400m.perd.100m, m.300m.perd.100m, m.200m.perd.100m, m.100m.perd.100m)


m.sur.rings.list.2014 <- flatten(list(m.rings.bio, m.rings.100m))

saveRDS(m.sur.rings.list.2014 , here("analysis", "data", "derived_data", "model", "model_2014", "model_outputs_2014", "m.space.rings.14.rds"))

