## ----libraries---------------------------------------------------------------

library(pacman)
pacman::p_load(INLA, tidyverse, sf, lubridate, parallel, here)


## ----data--------------------------------------------------------------------------------------------------

data <- readRDS(here("analysis", "data", "derived_data", "variables", "variables_2014", "vbles.time.14.rds"))

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


## ---- models including multiple yes/no variables that indicated if household i had been sprayed in the previous 1-6 weeks------------------------------------------------------------------------------------------------------------


stack <- inla.stack(
  # response variable number of Ae. aegypti
  data = list(y = list(data[["AA_TOTAL"]])),
  
  
  # the effects are organised in a list of lists.
  A = list(a1, 1),
  
  # spatial effect that accounts for spatial autocorrelation and the intercept
  effects = list(
    c(sp.index, list(Intercept = 1)),
    # and then all the other effects, in this a combination of spray.1week.lag and spray.2week.lag
    list(
      spray.1week.lag = data[["spray.1week.lag"]],
      spray.2week.lag = data[["spray.2week.lag"]],
      location_code = data[["location_code"]],
      month = data[["month"]]
    )
  ),
  tag = "est"
)

formula <- y ~ -1 + Intercept + spray.1week.lag + spray.2week.lag + f(location_code, model = "iid") + f(month, model = "rw1") + f(spatial.field, model = spde)

m.2.1 <- c(
  inla(formula,
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
  ),
  # add an element to the INLA results that names the fixed effects included in the model
  fixed.effect.s = "spray.1week.lag+ spray.2week.lag"
)


stack <- inla.stack(
  data = list(y = list(data[["AA_TOTAL"]])),
  A = list(a1, 1),
  effects = list(
    c(sp.index, list(Intercept = 1)),
    list(
      spray.1week.lag = data[["spray.1week.lag"]],
      spray.2week.lag = data[["spray.2week.lag"]],
      spray.3week.lag = data[["spray.3week.lag"]],
      location_code = data[["location_code"]],
      month = data[["month"]]
    )
  ),
  tag = "est"
)

formula <- y ~ -1 + Intercept + spray.1week.lag + spray.2week.lag + spray.3week.lag + f(location_code, model = "iid") + f(month, model = "rw1") + f(spatial.field, model = spde)

m.2.2 <- c(inla(formula,
                data = inla.stack.data(stack,
                                       spde = spde
                ),
                family = "nbinomial",
                control.family = list(link = "log"),
                control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE),
                control.predictor = list(
                  link = 1, A = inla.stack.A(stack),
                  compute = T
                )
), fixed.effect.s = "spray.1week.lag+ spray.2week.lag + spray.3week.lag")


stack <- inla.stack(
  data = list(y = list(data[["AA_TOTAL"]])),
  A = list(a1, 1),
  effects = list(
    c(sp.index, list(Intercept = 1)),
    list(
      spray.1week.lag = data[["spray.1week.lag"]],
      spray.2week.lag = data[["spray.2week.lag"]],
      spray.3week.lag = data[["spray.3week.lag"]],
      spray.4week.lag = data[["spray.4week.lag"]],
      location_code = data[["location_code"]],
      month = data[["month"]]
    )
  ),
  tag = "est"
)

formula <- y ~ -1 + Intercept + spray.1week.lag + spray.2week.lag + spray.3week.lag + spray.4week.lag + f(location_code, model = "iid") + f(month, model = "rw1") + f(spatial.field, model = spde)

m.2.3 <- c(inla(formula,
                data = inla.stack.data(stack,
                                       spde = spde
                ),
                family = "nbinomial",
                control.family = list(link = "log"),
                control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE),
                control.predictor = list(
                  link = 1, A = inla.stack.A(stack),
                  compute = T
                )
), fixed.effect.s = "spray.1week.lag+ spray.2week.lag + spray.3week.lag +spray.4week.lag")


stack <- inla.stack(
  data = list(y = list(data[["AA_TOTAL"]])),
  A = list(a1, 1),
  effects = list(
    c(sp.index, list(Intercept = 1)),
    list(
      spray.1week.lag = data[["spray.1week.lag"]],
      spray.2week.lag = data[["spray.2week.lag"]],
      spray.3week.lag = data[["spray.3week.lag"]], 
      spray.4week.lag = data[["spray.4week.lag"]],
      spray.5week.lag = data[["spray.5week.lag"]],
      location_code = data[["location_code"]],
      month = data[["month"]]
    )
  ),
  tag = "est"
)

formula <- y ~ -1 + Intercept + spray.1week.lag + spray.2week.lag + spray.3week.lag + spray.4week.lag + spray.5week.lag + f(location_code, model = "iid") + f(month, model = "rw1") + f(spatial.field, model = spde)

m.2.4 <- c(inla(formula,
                data = inla.stack.data(stack,
                                       spde = spde
                ),
                family = "nbinomial",
                control.family = list(link = "log"),
                control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE),
                control.predictor = list(
                  link = 1, A = inla.stack.A(stack),
                  compute = T
                )
), fixed.effect.s = "spray.1week.lag+ spray.2week.lag + spray.3week.lag +spray.4week.lag + spray.5week.lag")

stack <- inla.stack(
  data = list(y = list(data[["AA_TOTAL"]])),
  A = list(a1, 1),
  effects = list(
    c(sp.index, list(Intercept = 1)),
    list(
      spray.1week.lag = data[["spray.1week.lag"]],
      spray.2week.lag = data[["spray.2week.lag"]],
      spray.3week.lag = data[["spray.3week.lag"]],
      spray.4week.lag = data[["spray.4week.lag"]],
      spray.5week.lag = data[["spray.5week.lag"]],
      spray.6week.lag = data[["spray.6week.lag"]],
      location_code = data[["location_code"]],
      month = data[["month"]]
    )
  ),
  tag = "est"
)

formula <- y ~ -1 + Intercept + spray.1week.lag + spray.2week.lag + spray.3week.lag + spray.4week.lag + spray.5week.lag + spray.6week.lag + f(location_code, model = "iid") + f(month, model = "rw1") + f(spatial.field, model = spde)

m.2.5 <- c(inla(formula,
                data = inla.stack.data(stack,
                                       spde = spde
                ),
                family = "nbinomial",
                control.family = list(link = "log"),
                control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE),
                control.predictor = list(
                  link = 1, A = inla.stack.A(stack),
                  compute = T
                )
), fixed.effect.s = "spray.1week.lag+ spray.2week.lag + spray.3week.lag +spray.4week.lag + spray.5week.lag+ spray.6week.lag")


# produce a list of INLA models that include multiple yes/no variables that indicated if household i had been sprayed in the previous 1-6 weeks
m.week.lag <- list(m.2.1, m.2.2, m.2.3, m.2.4, m.2.5)


## ----models including a single variable that measured the effect of sprays within household i------------------------------------------------------------------------------------------

cores <- 2

#' Run INLA models including one variable from a vector of variables
#'
#' This function takes a vector of variables and a data set as input and returns a list of INLA models where each model includes one variable as a fixed effect.
#'
#' @param df A data set.
#' @param vble.list A vector of variables
#' @param ncores A numeric value of the number of computing cores to be used for this operation
#' @return A list of INLA models
#'
fx.model.1fe <- function(df, vble.list, ncores) {
  fx.model <- function(df, vble.list) {
    stack <- inla.stack(
      data = list(y = list(df[["AA_TOTAL"]])),
      A = list(a1, 1),
      effects = list(
        c(sp.index, list(Intercept = 1)),
        list(
          vble = df[[vble.list]],
          location_code = df[["location_code"]],
          month = df[["month"]]
        )
      ),
      tag = "est"
    )
    
    formula <- y ~ -1 + Intercept + vble + f(location_code, model = "iid") + f(month, model = "rw1") + f(spatial.field, model = spde)
    
    m <- inla(formula,
              data = inla.stack.data(stack,
                                     spde = spde
              ),
              family = "nbinomial",
              control.family = list(link = "log"),
              control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE),
              control.predictor = list(
                link = 1, A = inla.stack.A(stack),
                compute = T
              )
    )
    c(m, fixed.effect.s = vble.list)
  }
  
  mclapply(X = vble.list, FUN = fx.model, df = df, mc.cores = ncores)
}



## ----run models with 1 fixed effect for spray effect---------------------------------------------------------------------------


l.vbles.time <- list(
  ## discrete and continuous variables
  
  # house in spray zone/house in buffer zone
  "zone",
  # how many times was house i sprayed before the date of the adult survey t? 1,2...6
  "n.spray.by.date",
  # did any sprays occur by the date of the  adult survey t yes/no?
  "spray.by.date",
  # days since the most recent spray
  "day.min",
  # weeks since the most recent spray
  "w.min",
  
  # individual yes/no variables that indicated if household i had been sprayed in the previous 1-6 weeks
  "spray.1week.lag", "spray.2week.lag", "spray.3week.lag", "spray.4week.lag",
  "spray.5week.lag", "spray.6week.lag",
  
  ## variables weighting the spray effect using a decay function
  
  # cumulative inverse weight of all previous sprays
  "i_day",
  # inverse weight of the most recent spray
  "i_day.max",
  
  # cumulative exponential weight of all previous sprays (with varying k)
  "e0.005_day", "e0.01_day", "e0.015_day", "e0.02_day",
  "e0.03_day", "e0.04_day", "e0.06_day", "e0.1_day",
  "e0.2_day", "e0.4_day", "e1_day",
  # exponential weight of the most recent spray (with varying k)
  "e0.005_day.max", "e0.01_day.max", "e0.015_day.max", "e0.02_day.max",
  "e0.03_day.max", "e0.04_day.max", "e0.06_day.max", "e0.1_day.max",
  "e0.2_day.max", "e0.4_day.max", "e1_day.max",
  
  # cumulative gaussian weight of all previous sprays (with varying sigma)
  "g1_day", "g3_day", "g5_day", "g7_day",
  "g10_day", "g15_day", "g20_day", "g25_day",
  "g30_day", "g35_day", "g40_day", "g50_day",
  "g60_day", "g80_day",
  # gaussian weight of the most recent spray (with varying sigma)
  "g1_day.max", "g3_day.max", "g5_day.max",
  "g7_day.max", "g10_day.max", "g15_day.max",
  "g20_day.max", "g25_day.max", "g30_day.max",
  "g35_day.max", "g40_day.max", "g50_day.max",
  "g60_day.max", "g80_day.max"
)

# apply function to the data and produce a list of INLA models that include the individual variables
m.time.list <- fx.model.1fe(data, l.vbles.time, ncores = cores)

# create a list that includes all variables that measure the effect of sprays within the household
m.time.list <- flatten(list(m.time.list, m.week.lag))


saveRDS(m.time.list, here("analysis", "data", "derived_data", "model", "model_2014", "model_outputs_2014", "m.time.14.rds"))

