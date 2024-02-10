
library(pacman)
pacman::p_load(tidyverse, sf, lubridate, parallel, here)


## ----data---------------------------------------------------------------------------------------------------------------------------------------------

vbles.time.13 <- readRDS(here("analysis", "data", "derived_data", "variables", "variables_2013", "vbles.time.13.rds")) %>% 
  select(c("location_code","month","date", "AA_TOTAL","g20_day.max", "cycle_1", "cycle_2", "cycle_3", "cycle_4", 
           "cycle_5", "cycle_6", "geometry"))

## ----new surrounding houses function 2022-------------------------------------------------------------------------------------------------------------

cores <- 2

#' Calculate weights of effects of sprays in neighboring households (j) for every household i
#'
#' This function takes a data set of households with the date adult surveys and spray events occurred for each household as input
#' and returns the added weight of sprays in the surrounding houses (j) given their distance from household i (in m)
#' and also how many days before the adult survey in household i each spray in household j occurred.
#'
#' @param df A data frame with household points and dates of spray and adult survey events.
#' @param y A numeric value of the number of computing cores to be used for this operation.
#' @return A data frame with the cumulative spray weight of all surrounding houses for each household i and time t

fx.w <- function(df, ncores) {
  fx.dist.w <- function(x) {
    # center.geometry = geolocated point of household i
    center.geometry <- df[["geometry"]][x]
    # center.location = location_code for household i
    center.location <- df[["location_code"]][x]
    # center.date = date of adult survey for household i (t)
    center.date <- df[["date"]][[x]]

    # for every household i, we calculate the distance of sprays in time and space for every other household in the data frame (j)
    df.w <- distinct(df %>%
      select(c(
        "location_code", "cycle_1", "cycle_2", "cycle_3", "cycle_4",
        "cycle_5", "cycle_6", "geometry"
      ))) %>%
      # dist = the distance (in m) between the household i and every other household j in the study, or dij.
      mutate(dist = round(as.numeric(st_distance(center.geometry, geometry)))) %>%
      # once we calculate the distance we don't need it to be an sf object anymore
      # for the following analyses it's easier to work with data frame rather than spatial objects
      st_drop_geometry() %>%
      # diff_c[1-6] = what is the difference in days between the entomological collection date and the date the household was sprayed?
      # date entomological survey - each of the spray dates in the 6 spray cycles

      # add 1 to all time differences because if the spray occurred on the same date as the entomological collection, the time difference is 0

      # in S-2013, adult surveys were typically carried out during the spray period on Monday afternoons just prior to the initiation of each spray cycle (Gunning2018)
      # thus, the entomological collection occurred BEFORE the spray events in the 2013 experiments,
      # because we want to look at the effect of the spray on the number of Ae. aegypti,
      # we want to remove the effect of sprays occurring before the adult collection (as these sprays could not affect the number of mosquitoes)
      # diff_c[1-6]= 1 means the spray occurred on same day, diff_c[1-6] >1 means spray was in the past, diff_c[1-6] <1 means spray is in the future
      # to make sure sprays occurred before or on the same day as the adult collections don't have a spray effect weight in the model,
      # we turn values that are <1 or equal to 1 into NA
      mutate(
        diff_c1 = ifelse((center.date - cycle_1 + 1) <= 1, NA, (center.date - cycle_1 + 1)),
        diff_c2 = ifelse((center.date - cycle_2 + 1) <= 1, NA, (center.date - cycle_2 + 1)),
        diff_c3 = ifelse((center.date - cycle_3 + 1) <= 1, NA, (center.date - cycle_3 + 1)),
        diff_c4 = ifelse((center.date - cycle_4 + 1) <= 1, NA, (center.date - cycle_4 + 1)),
        diff_c5 = ifelse((center.date - cycle_5 + 1) <= 1, NA, (center.date - cycle_5 + 1)),
        diff_c6 = ifelse((center.date - cycle_6 + 1) <= 1, NA, (center.date - cycle_6 + 1))
      ) %>%
      # remove household i from the calculations (as we want to add the spray weights of the surrounding houses not counting sprays in household i)
      filter(location_code != center.location)

    df.w <- df.w %>%
      # day.min= days since the most recent spray for each surrounding household j 
      mutate(day.min = do.call(pmin, c(df.w[, c("diff_c1", "diff_c2", "diff_c3", "diff_c4", "diff_c5", "diff_c6")], list(na.rm = TRUE)))) %>%
      # create a variable that calculates the weighted measure of “time since spray event”
      # using the weighting scheme selected from the best model that measured the within-household spray effect,
      # in this case the weight of time since the most recent spray event according to a Gaussian function with 𝜎 = 20
      mutate(g20_day.min = exp(-(day.min^2) / (2 * 20^2))) %>%
      # assign weights to the distance (in m) between the household i and every surrounding household j, or dij,
      # using the inverse, gaussian and exponential decay functions,
      # accounting for spray events in household j prior to t by
      # multiplying the weighted measure of distance by the weighted measure of “time since spray event” using a Gaussian function with 𝜎 = 20

      # the assigned weights for every household j are then added
      # to create a single weight value of surrounding sprays per household i and time t
      mutate(

        # decay function used to assign a weight to dij:

        # INVERSE : dist^(-1)
        i_total = sum((1 / dist) * g20_day.min, na.rm = TRUE),

        # GAUSSIAN : exp(-(dist^2)/(2*sigma^2))
        # sigma= 5, 25, 50, 75, 100, 125, 150, 200, 250, 300

        g5_total = sum(exp(-(dist^2) / (2 * 5^2)) * g20_day.min, na.rm = TRUE),
        g25_total = sum(exp(-(dist^2) / (2 * 25^2)) * g20_day.min, na.rm = TRUE),
        g50_total = sum(exp(-(dist^2) / (2 * 50^2)) * g20_day.min, na.rm = TRUE),
        g75_total = sum(exp(-(dist^2) / (2 * 75^2)) * g20_day.min, na.rm = TRUE),
        g100_total = sum(exp(-(dist^2) / (2 * 100^2)) * g20_day.min, na.rm = TRUE),
        g125_total = sum(exp(-(dist^2) / (2 * 125^2)) * g20_day.min, na.rm = TRUE),
        g150_total = sum(exp(-(dist^2) / (2 * 150^2)) * g20_day.min, na.rm = TRUE),
        g200_total = sum(exp(-(dist^2) / (2 * 200^2)) * g20_day.min, na.rm = TRUE),
        g250_total = sum(exp(-(dist^2) / (2 * 250^2)) * g20_day.min, na.rm = TRUE),
        g300_total = sum(exp(-(dist^2) / (2 * 300^2)) * g20_day.min, na.rm = TRUE),

        # EXPONENTIAL : exp(-(k)*(dist))
        # k=0.0025, 0.0035, 0.005, 0.0075, 0.01, 0.0125, 0.02, 0.045, 0.2

        e0.0025_total = sum(exp(-(0.0025) * (dist)) * g20_day.min, na.rm = TRUE),
        e0.0035_total = sum(exp(-(0.0035) * (dist)) * g20_day.min, na.rm = TRUE),
        e0.005_total = sum(exp(-(0.005) * (dist)) * g20_day.min, na.rm = TRUE),
        e0.0075_total = sum(exp(-(0.0075) * (dist)) * g20_day.min, na.rm = TRUE),
        e0.01_total = sum(exp(-(0.01) * (dist)) * g20_day.min, na.rm = TRUE),
        e0.0125_total = sum(exp(-(0.0125) * (dist)) * g20_day.min, na.rm = TRUE),
        e0.02_total = sum(exp(-(0.02) * (dist)) * g20_day.min, na.rm = TRUE),
        e0.045_total = sum(exp(-(0.045) * (dist)) * g20_day.min, na.rm = TRUE),
        e0.2_total = sum(exp(-(0.2) * (dist)) * g20_day.min, na.rm = TRUE)
      ) %>%
      select(c(
        "i_total", "g5_total", "g25_total",
        "g50_total", "g75_total", "g100_total", "g125_total", "g150_total",
        "g200_total", "g250_total", "g300_total", "e0.0025_total", "e0.0035_total",
        "e0.005_total", "e0.0075_total", "e0.01_total", "e0.0125_total",
        "e0.02_total", "e0.045_total", "e0.2_total"
      )) %>%
      unique() %>%
      # add identifiers that indicates which household i and adult survey time t the weights correspond to
      # these weights are the sum of all spray weights for houses surrounding household i at time t
      mutate(
        center.location = center.location,
        center.date = center.date
      )
  }


  dist.w <- mclapply(1:nrow(df), fx.dist.w, mc.cores = ncores) %>% bind_rows()

  df.dist.w <- df %>%
    left_join(dist.w, by = c("location_code" = "center.location", "date" = "center.date"))
}


w.2013 <- fx.w(vbles.time.13, cores)

saveRDS(w.2013, here("analysis", "data", "derived_data", "variables", "variables_2013", "vbles.space.decay.13.rds"))
