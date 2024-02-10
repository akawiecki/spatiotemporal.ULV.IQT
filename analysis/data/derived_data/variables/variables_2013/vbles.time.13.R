## ----libraries-------------------------------------------------------------------------------------------------------------------------------

library(pacman)
pacman::p_load(tidyverse, sf, lubridate, parallel, here)


## ----read in data-------------------------------------------------------------------------------------------------------------------------------------------
data.raw <- readRDS(here("analysis", "data", "raw_data", "data.sample.rds"))

data <- data.raw %>%
  # takes away the spatial attributes
  # necessary for the rowSums function further down, to add inverse days across cols
  st_drop_geometry()

## ----create rds data files----------------------------------------------------------------------------------------------------------------------------------

data.13 <- data %>%
  # n_trt = how many spray events occurred in each household?
  mutate(n_trt = rowSums(!is.na(data[, c("cycle_1", "cycle_2", "cycle_3", "cycle_4", "cycle_5", "cycle_6")]))) %>%
  mutate(
    month = month(date),
    month.f = as.factor(month(date)),
    week = isoweek(date),
    week.f = as.factor(isoweek(date)),
    zone = as.factor(zone),
    lat = unlist(map(data.raw$geometry, 1)),
    long = unlist(map(data.raw$geometry, 2))
  ) %>%
  mutate(across(c("cycle_1", "cycle_2", "cycle_3", "cycle_4", "cycle_5", "cycle_6"), ymd)) %>%
  mutate(across(c("exp", "e_circuit", "moment", "org", "zone", ), as.factor)) %>%
  # selects only data from the 2013 experiment
  filter(exp == "2013") %>%
  # sprayed =  was house i sprayed previously? yes/no
  mutate(sprayed = if_else(n_trt == 0, 0, 1)) %>%
  # diff_c[1-6] = what is the difference in days between the entomological collection date and the date the household was sprayed?
  # date entomological survey (t) - each of the spray dates in the 6 spray cycles

  # add 1 to all time differences because if the spray occurred on the same date as the entomological collection, the time difference is 0

  # in S-2013, adult surveys were typically carried out during the spray period on Monday afternoons just prior to the initiation of each spray cycle (Gunning2018)
  # thus, the entomological collection occurred BEFORE the spray events in the 2013 experiments,
  # because we want to look at the effect of the spray on the number of Ae. aegypti,
  # we want to remove the effect of sprays occurring before the adult collection (as these sprays could not affect the number of mosquitoes)
  # diff_c[1-6]= 1 means the spray occurred on same day, diff_c[1-6] >1 means spray was in the past, diff_c[1-6] <1 means spray is in the future
  # to make sure sprays occurred before or on the same day as the adult collections don't have a spray effect weight in the model,
  # we turn values that are <1 or equal to 1 into NA
  mutate(
    diff_c1 = ifelse((date - cycle_1 + 1) <= 1, NA, (date - cycle_1 + 1)),
    diff_c2 = ifelse((date - cycle_2 + 1) <= 1, NA, (date - cycle_2 + 1)),
    diff_c3 = ifelse((date - cycle_3 + 1) <= 1, NA, (date - cycle_3 + 1)),
    diff_c4 = ifelse((date - cycle_4 + 1) <= 1, NA, (date - cycle_4 + 1)),
    diff_c5 = ifelse((date - cycle_5 + 1) <= 1, NA, (date - cycle_5 + 1)),
    diff_c6 = ifelse((date - cycle_6 + 1) <= 1, NA, (date - cycle_6 + 1))
  ) %>%
  # w_c[1-6] = how many weeks ago was the spray?
  # because we're adding +1 to the date difference we add +1 to the 7 so the actual difference in days is divided by 7
  mutate(
    w_c1 = diff_c1 / (7 + 1),
    w_c2 = diff_c2 / (7 + 1),
    w_c3 = diff_c3 / (7 + 1),
    w_c4 = diff_c4 / (7 + 1),
    w_c5 = diff_c5 / (7 + 1),
    w_c6 = diff_c6 / (7 + 1)
  ) %>% 
  
  # assign a weight to the difference in time between the spray event and the adult survey
  # based on the inverse decay function
  ### INVERSE : days^(-1) ###
  
  dplyr::mutate(
    i_day1 = 1 / (as.numeric(diff_c1)),
    i_day2 = 1 / (as.numeric(diff_c2)),
    i_day3 = 1 / (as.numeric(diff_c3)),
    i_day4 = 1 / (as.numeric(diff_c4)),
    i_day5 = 1 / (as.numeric(diff_c5)),
    i_day6 = 1 / (as.numeric(diff_c6)),
  )


data.13 <- data.13 %>%
  
  mutate(
    # day.min= days since the most recent spray
    day.min = do.call(pmin, c(data.13[, c("diff_c1", "diff_c2", "diff_c3", "diff_c4", "diff_c5", "diff_c6")], list(na.rm = TRUE))),
    # w.min= weeks since the most recent spray
    w.min = do.call(pmin, c(data.13[, c("w_c1", "w_c2", "w_c3", "w_c4", "w_c5", "w_c6")], list(na.rm = TRUE))),
    # i_day.max= inverse weight of the most recent spray
    i_day.max = do.call(pmax, c(data.13[, c("i_day1", "i_day2", "i_day3", "i_day4", "i_day5", "i_day6")], list(na.rm = TRUE)))
  ) %>%
  mutate(
    # assign a weight to the difference in time between the spray event and the adult survey
    # based on a decay function: 
    
    ### EXPONENTIAL: exp(-(k)*(days)) ###

    # k <- 0.005 0.010 0.015 0.020 0.030 0.040 0.060 0.100 0.200 0.400 1.000
    e0.005_day1 = exp(-(0.005) * diff_c1),
    e0.005_day2 = exp(-(0.005) * diff_c2),
    e0.005_day3 = exp(-(0.005) * diff_c3),
    e0.005_day4 = exp(-(0.005) * diff_c4),
    e0.005_day5 = exp(-(0.005) * diff_c5),
    e0.005_day6 = exp(-(0.005) * diff_c6),
    
    e0.01_day1 = exp(-(0.01) * diff_c1),
    e0.01_day2 = exp(-(0.01) * diff_c2),
    e0.01_day3 = exp(-(0.01) * diff_c3),
    e0.01_day4 = exp(-(0.01) * diff_c4),
    e0.01_day5 = exp(-(0.01) * diff_c5),
    e0.01_day6 = exp(-(0.01) * diff_c6),
    
    e0.015_day1 = exp(-(0.015) * diff_c1),
    e0.015_day2 = exp(-(0.015) * diff_c2),
    e0.015_day3 = exp(-(0.015) * diff_c3),
    e0.015_day4 = exp(-(0.015) * diff_c4),
    e0.015_day5 = exp(-(0.015) * diff_c5),
    e0.015_day6 = exp(-(0.015) * diff_c6),
    
    e0.02_day1 = exp(-(0.02) * diff_c1),
    e0.02_day2 = exp(-(0.02) * diff_c2),
    e0.02_day3 = exp(-(0.02) * diff_c3),
    e0.02_day4 = exp(-(0.02) * diff_c4),
    e0.02_day5 = exp(-(0.02) * diff_c5),
    e0.02_day6 = exp(-(0.02) * diff_c6),
    
    e0.03_day1 = exp(-(0.03) * diff_c1),
    e0.03_day2 = exp(-(0.03) * diff_c2),
    e0.03_day3 = exp(-(0.03) * diff_c3),
    e0.03_day4 = exp(-(0.03) * diff_c4),
    e0.03_day5 = exp(-(0.03) * diff_c5),
    e0.03_day6 = exp(-(0.03) * diff_c6),
    
    e0.04_day1 = exp(-(0.04) * diff_c1),
    e0.04_day2 = exp(-(0.04) * diff_c2),
    e0.04_day3 = exp(-(0.04) * diff_c3),
    e0.04_day4 = exp(-(0.04) * diff_c4),
    e0.04_day5 = exp(-(0.04) * diff_c5),
    e0.04_day6 = exp(-(0.04) * diff_c6),
    
    e0.06_day1 = exp(-(0.06) * diff_c1),
    e0.06_day2 = exp(-(0.06) * diff_c2),
    e0.06_day3 = exp(-(0.06) * diff_c3),
    e0.06_day4 = exp(-(0.06) * diff_c4),
    e0.06_day5 = exp(-(0.06) * diff_c5),
    e0.06_day6 = exp(-(0.06) * diff_c6),
    
    e0.1_day1 = exp(-(0.1) * diff_c1),
    e0.1_day2 = exp(-(0.1) * diff_c2),
    e0.1_day3 = exp(-(0.1) * diff_c3),
    e0.1_day4 = exp(-(0.1) * diff_c4),
    e0.1_day5 = exp(-(0.1) * diff_c5),
    e0.1_day6 = exp(-(0.1) * diff_c6),
    
    e0.2_day1 = exp(-(0.2) * diff_c1),
    e0.2_day2 = exp(-(0.2) * diff_c2),
    e0.2_day3 = exp(-(0.2) * diff_c3),
    e0.2_day4 = exp(-(0.2) * diff_c4),
    e0.2_day5 = exp(-(0.2) * diff_c5),
    e0.2_day6 = exp(-(0.2) * diff_c6),
    
    e0.4_day1 = exp(-(0.4) * diff_c1),
    e0.4_day2 = exp(-(0.4) * diff_c2),
    e0.4_day3 = exp(-(0.4) * diff_c3),
    e0.4_day4 = exp(-(0.4) * diff_c4),
    e0.4_day5 = exp(-(0.4) * diff_c5),
    e0.4_day6 = exp(-(0.4) * diff_c6),
    
    e1_day1 = exp(-(1) * diff_c1),
    e1_day2 = exp(-(1) * diff_c2),
    e1_day3 = exp(-(1) * diff_c3),
    e1_day4 = exp(-(1) * diff_c4),
    e1_day5 = exp(-(1) * diff_c5),
    e1_day6 = exp(-(1) * diff_c6),

    ### GAUSSIAN: exp(-(days^2)/(2*sigma^2) ###

    ### sigma <- 1  3  5  7 10 15 20 25 30 35 40 50 60 80

    g1_day1 = exp(-(diff_c1^2) / (2 * 1^2)),
    g1_day2 = exp(-(diff_c2^2) / (2 * 1^2)),
    g1_day3 = exp(-(diff_c3^2) / (2 * 1^2)),
    g1_day4 = exp(-(diff_c4^2) / (2 * 1^2)),
    g1_day5 = exp(-(diff_c5^2) / (2 * 1^2)),
    g1_day6 = exp(-(diff_c6^2) / (2 * 1^2)),
    
    g3_day1 = exp(-(diff_c1^2) / (2 * 3^2)),
    g3_day2 = exp(-(diff_c2^2) / (2 * 3^2)),
    g3_day3 = exp(-(diff_c3^2) / (2 * 3^2)),
    g3_day4 = exp(-(diff_c4^2) / (2 * 3^2)),
    g3_day5 = exp(-(diff_c5^2) / (2 * 3^2)),
    g3_day6 = exp(-(diff_c6^2) / (2 * 3^2)),
    
    g5_day1 = exp(-(diff_c1^2) / (2 * 5^2)),
    g5_day2 = exp(-(diff_c2^2) / (2 * 5^2)),
    g5_day3 = exp(-(diff_c3^2) / (2 * 5^2)),
    g5_day4 = exp(-(diff_c4^2) / (2 * 5^2)),
    g5_day5 = exp(-(diff_c5^2) / (2 * 5^2)),
    g5_day6 = exp(-(diff_c6^2) / (2 * 5^2)),
    
    g7_day1 = exp(-(diff_c1^2) / (2 * 7^2)),
    g7_day2 = exp(-(diff_c2^2) / (2 * 7^2)),
    g7_day3 = exp(-(diff_c3^2) / (2 * 7^2)),
    g7_day4 = exp(-(diff_c4^2) / (2 * 7^2)),
    g7_day5 = exp(-(diff_c5^2) / (2 * 7^2)),
    g7_day6 = exp(-(diff_c6^2) / (2 * 7^2)),
    
    g10_day1 = exp(-(diff_c1^2) / (2 * 10^2)),
    g10_day2 = exp(-(diff_c2^2) / (2 * 10^2)),
    g10_day3 = exp(-(diff_c3^2) / (2 * 10^2)),
    g10_day4 = exp(-(diff_c4^2) / (2 * 10^2)),
    g10_day5 = exp(-(diff_c5^2) / (2 * 10^2)),
    g10_day6 = exp(-(diff_c6^2) / (2 * 10^2)),
    
    g15_day1 = exp(-(diff_c1^2) / (2 * 15^2)),
    g15_day2 = exp(-(diff_c2^2) / (2 * 15^2)),
    g15_day3 = exp(-(diff_c3^2) / (2 * 15^2)),
    g15_day4 = exp(-(diff_c4^2) / (2 * 15^2)),
    g15_day5 = exp(-(diff_c5^2) / (2 * 15^2)),
    g15_day6 = exp(-(diff_c6^2) / (2 * 15^2)),
    
    g20_day1 = exp(-(diff_c1^2) / (2 * 20^2)),
    g20_day2 = exp(-(diff_c2^2) / (2 * 20^2)),
    g20_day3 = exp(-(diff_c3^2) / (2 * 20^2)),
    g20_day4 = exp(-(diff_c4^2) / (2 * 20^2)),
    g20_day5 = exp(-(diff_c5^2) / (2 * 20^2)),
    g20_day6 = exp(-(diff_c6^2) / (2 * 20^2)),
    
    g25_day1 = exp(-(diff_c1^2) / (2 * 25^2)),
    g25_day2 = exp(-(diff_c2^2) / (2 * 25^2)),
    g25_day3 = exp(-(diff_c3^2) / (2 * 25^2)),
    g25_day4 = exp(-(diff_c4^2) / (2 * 25^2)),
    g25_day5 = exp(-(diff_c5^2) / (2 * 25^2)),
    g25_day6 = exp(-(diff_c6^2) / (2 * 25^2)),
    
    g30_day1 = exp(-(diff_c1^2) / (2 * 30^2)),
    g30_day2 = exp(-(diff_c2^2) / (2 * 30^2)),
    g30_day3 = exp(-(diff_c3^2) / (2 * 30^2)),
    g30_day4 = exp(-(diff_c4^2) / (2 * 30^2)),
    g30_day5 = exp(-(diff_c5^2) / (2 * 30^2)),
    g30_day6 = exp(-(diff_c6^2) / (2 * 30^2)),
    
    g35_day1 = exp(-(diff_c1^2) / (2 * 35^2)),
    g35_day2 = exp(-(diff_c2^2) / (2 * 35^2)),
    g35_day3 = exp(-(diff_c3^2) / (2 * 35^2)),
    g35_day4 = exp(-(diff_c4^2) / (2 * 35^2)),
    g35_day5 = exp(-(diff_c5^2) / (2 * 35^2)),
    g35_day6 = exp(-(diff_c6^2) / (2 * 35^2)),
    
    g40_day1 = exp(-(diff_c1^2) / (2 * 40^2)),
    g40_day2 = exp(-(diff_c2^2) / (2 * 40^2)),
    g40_day3 = exp(-(diff_c3^2) / (2 * 40^2)),
    g40_day4 = exp(-(diff_c4^2) / (2 * 40^2)),
    g40_day5 = exp(-(diff_c5^2) / (2 * 40^2)),
    g40_day6 = exp(-(diff_c6^2) / (2 * 40^2)),
    
    g50_day1 = exp(-(diff_c1^2) / (2 * 50^2)),
    g50_day2 = exp(-(diff_c2^2) / (2 * 50^2)),
    g50_day3 = exp(-(diff_c3^2) / (2 * 50^2)),
    g50_day4 = exp(-(diff_c4^2) / (2 * 50^2)),
    g50_day5 = exp(-(diff_c5^2) / (2 * 50^2)),
    g50_day6 = exp(-(diff_c6^2) / (2 * 50^2)),
    
    g60_day1 = exp(-(diff_c1^2) / (2 * 60^2)),
    g60_day2 = exp(-(diff_c2^2) / (2 * 60^2)),
    g60_day3 = exp(-(diff_c3^2) / (2 * 60^2)),
    g60_day4 = exp(-(diff_c4^2) / (2 * 60^2)),
    g60_day5 = exp(-(diff_c5^2) / (2 * 60^2)),
    g60_day6 = exp(-(diff_c6^2) / (2 * 60^2)),
    
    g80_day1 = exp(-(diff_c1^2) / (2 * 80^2)),
    g80_day2 = exp(-(diff_c2^2) / (2 * 80^2)),
    g80_day3 = exp(-(diff_c3^2) / (2 * 80^2)),
    g80_day4 = exp(-(diff_c4^2) / (2 * 80^2)),
    g80_day5 = exp(-(diff_c5^2) / (2 * 80^2)),
    g80_day6 = exp(-(diff_c6^2) / (2 * 80^2))
  )

data.13 <- data.13 %>%
  # n.spray.by.date = how many times was house i sprayed before the date of the adult survey t? 1,2...6
  mutate(n.spray.by.date = 6 - rowSums(is.na(data.13[, c("diff_c1", "diff_c2", "diff_c3", "diff_c4", "diff_c5", "diff_c6")]))) %>%
  # spray.by.date = did any sprays occur by the date of the  adult survey t yes/no?
  mutate(spray.by.date = case_when(
    n.spray.by.date > 0 ~ 1,
    n.spray.by.date == 0 ~ 0
  )) %>%
  # f(x)sigma_day = cumulative weight of all past sprays for a given household and date (sum of all past spray weights)
  mutate(
    i_day = rowSums(data.13[, c("i_day1", "i_day2", "i_day3", "i_day4", "i_day5", "i_day6")], na.rm = TRUE),
    e0.005_day = rowSums(data.13[, c("e0.005_day1", "e0.005_day2", "e0.005_day3", "e0.005_day4", "e0.005_day5", "e0.005_day6")], na.rm = TRUE),
    e0.01_day = rowSums(data.13[, c("e0.01_day1", "e0.01_day2", "e0.01_day3", "e0.01_day4", "e0.01_day5", "e0.01_day6")], na.rm = TRUE),
    e0.015_day = rowSums(data.13[, c("e0.015_day1", "e0.015_day2", "e0.015_day3", "e0.015_day4", "e0.015_day5", "e0.015_day6")], na.rm = TRUE),
    e0.02_day = rowSums(data.13[, c("e0.02_day1", "e0.02_day2", "e0.02_day3", "e0.02_day4", "e0.02_day5", "e0.02_day6")], na.rm = TRUE),
    e0.03_day = rowSums(data.13[, c("e0.03_day1", "e0.03_day2", "e0.03_day3", "e0.03_day4", "e0.03_day5", "e0.03_day6")], na.rm = TRUE),
    e0.04_day = rowSums(data.13[, c("e0.04_day1", "e0.04_day2", "e0.04_day3", "e0.04_day4", "e0.04_day5", "e0.04_day6")], na.rm = TRUE),
    e0.06_day = rowSums(data.13[, c("e0.06_day1", "e0.06_day2", "e0.06_day3", "e0.06_day4", "e0.06_day5", "e0.06_day6")], na.rm = TRUE),
    e0.1_day = rowSums(data.13[, c("e0.1_day1", "e0.1_day2", "e0.1_day3", "e0.1_day4", "e0.1_day5", "e0.1_day6")], na.rm = TRUE),
    e0.2_day = rowSums(data.13[, c("e0.2_day1", "e0.2_day2", "e0.2_day3", "e0.2_day4", "e0.2_day5", "e0.2_day6")], na.rm = TRUE),
    e0.4_day = rowSums(data.13[, c("e0.4_day1", "e0.4_day2", "e0.4_day3", "e0.4_day4", "e0.4_day5", "e0.4_day6")], na.rm = TRUE),
    e1_day = rowSums(data.13[, c("e1_day1", "e1_day2", "e1_day3", "e1_day4", "e1_day5", "e1_day6")], na.rm = TRUE),
    g1_day = rowSums(data.13[, c("g1_day1", "g1_day2", "g1_day3", "g1_day4", "g1_day5", "g1_day6")], na.rm = TRUE),
    g3_day = rowSums(data.13[, c("g3_day1", "g3_day2", "g3_day3", "g3_day4", "g3_day5", "g3_day6")], na.rm = TRUE),
    g5_day = rowSums(data.13[, c("g5_day1", "g5_day2", "g5_day3", "g5_day4", "g5_day5", "g5_day6")], na.rm = TRUE),
    g7_day = rowSums(data.13[, c("g7_day1", "g7_day2", "g7_day3", "g7_day4", "g7_day5", "g7_day6")], na.rm = TRUE),
    g10_day = rowSums(data.13[, c("g10_day1", "g10_day2", "g10_day3", "g10_day4", "g10_day5", "g10_day6")], na.rm = TRUE),
    g15_day = rowSums(data.13[, c("g15_day1", "g15_day2", "g15_day3", "g15_day4", "g15_day5", "g15_day6")], na.rm = TRUE),
    g20_day = rowSums(data.13[, c("g20_day1", "g20_day2", "g20_day3", "g20_day4", "g20_day5", "g20_day6")], na.rm = TRUE),
    g25_day = rowSums(data.13[, c("g25_day1", "g25_day2", "g25_day3", "g25_day4", "g25_day5", "g25_day6")], na.rm = TRUE),
    g30_day = rowSums(data.13[, c("g30_day1", "g30_day2", "g30_day3", "g30_day4", "g30_day5", "g30_day6")], na.rm = TRUE),
    g35_day = rowSums(data.13[, c("g35_day1", "g35_day2", "g35_day3", "g35_day4", "g35_day5", "g35_day6")], na.rm = TRUE),
    g40_day = rowSums(data.13[, c("g40_day1", "g40_day2", "g40_day3", "g40_day4", "g40_day5", "g40_day6")], na.rm = TRUE),
    g50_day = rowSums(data.13[, c("g50_day1", "g50_day2", "g50_day3", "g50_day4", "g50_day5", "g50_day6")], na.rm = TRUE),
    g60_day = rowSums(data.13[, c("g60_day1", "g60_day2", "g60_day3", "g60_day4", "g60_day5", "g60_day6")], na.rm = TRUE),
    g80_day = rowSums(data.13[, c("g80_day1", "g80_day2", "g80_day3", "g80_day4", "g80_day5", "g80_day6")], na.rm = TRUE)
  ) %>%
  # f(x)sigma_day.max = weight of the most recent past spray
  mutate(
    e0.005_day.max = do.call(pmax, c(data.13[, c("e0.005_day1", "e0.005_day2", "e0.005_day3", "e0.005_day4", "e0.005_day5", "e0.005_day6")], list(na.rm = TRUE))),
    e0.01_day.max = do.call(pmax, c(data.13[, c("e0.01_day1", "e0.01_day2", "e0.01_day3", "e0.01_day4", "e0.01_day5", "e0.01_day6")], list(na.rm = TRUE))),
    e0.015_day.max = do.call(pmax, c(data.13[, c("e0.015_day1", "e0.015_day2", "e0.015_day3", "e0.015_day4", "e0.015_day5", "e0.015_day6")], list(na.rm = TRUE))),
    e0.02_day.max = do.call(pmax, c(data.13[, c("e0.02_day1", "e0.02_day2", "e0.02_day3", "e0.02_day4", "e0.02_day5", "e0.02_day6")], list(na.rm = TRUE))),
    e0.03_day.max = do.call(pmax, c(data.13[, c("e0.03_day1", "e0.03_day2", "e0.03_day3", "e0.03_day4", "e0.03_day5", "e0.03_day6")], list(na.rm = TRUE))),
    e0.04_day.max = do.call(pmax, c(data.13[, c("e0.04_day1", "e0.04_day2", "e0.04_day3", "e0.04_day4", "e0.04_day5", "e0.04_day6")], list(na.rm = TRUE))),
    e0.06_day.max = do.call(pmax, c(data.13[, c("e0.06_day1", "e0.06_day2", "e0.06_day3", "e0.06_day4", "e0.06_day5", "e0.06_day6")], list(na.rm = TRUE))),
    e0.1_day.max = do.call(pmax, c(data.13[, c("e0.1_day1", "e0.1_day2", "e0.1_day3", "e0.1_day4", "e0.1_day5", "e0.1_day6")], list(na.rm = TRUE))),
    e0.2_day.max = do.call(pmax, c(data.13[, c("e0.2_day1", "e0.2_day2", "e0.2_day3", "e0.2_day4", "e0.2_day5", "e0.2_day6")], list(na.rm = TRUE))),
    e0.4_day.max = do.call(pmax, c(data.13[, c("e0.4_day1", "e0.4_day2", "e0.4_day3", "e0.4_day4", "e0.4_day5", "e0.4_day6")], list(na.rm = TRUE))),
    e1_day.max = do.call(pmax, c(data.13[, c("e1_day1", "e1_day2", "e1_day3", "e1_day4", "e1_day5", "e1_day6")], list(na.rm = TRUE))),
    g1_day.max = do.call(pmax, c(data.13[, c("g1_day1", "g1_day2", "g1_day3", "g1_day4", "g1_day5", "g1_day6")], list(na.rm = TRUE))),
    g3_day.max = do.call(pmax, c(data.13[, c("g3_day1", "g3_day2", "g3_day3", "g3_day4", "g3_day5", "g3_day6")], list(na.rm = TRUE))),
    g5_day.max = do.call(pmax, c(data.13[, c("g5_day1", "g5_day2", "g5_day3", "g5_day4", "g5_day5", "g5_day6")], list(na.rm = TRUE))),
    g7_day.max = do.call(pmax, c(data.13[, c("g7_day1", "g7_day2", "g7_day3", "g7_day4", "g7_day5", "g7_day6")], list(na.rm = TRUE))),
    g10_day.max = do.call(pmax, c(data.13[, c("g10_day1", "g10_day2", "g10_day3", "g10_day4", "g10_day5", "g10_day6")], list(na.rm = TRUE))),
    g15_day.max = do.call(pmax, c(data.13[, c("g15_day1", "g15_day2", "g15_day3", "g15_day4", "g15_day5", "g15_day6")], list(na.rm = TRUE))),
    g20_day.max = do.call(pmax, c(data.13[, c("g20_day1", "g20_day2", "g20_day3", "g20_day4", "g20_day5", "g20_day6")], list(na.rm = TRUE))),
    g25_day.max = do.call(pmax, c(data.13[, c("g25_day1", "g25_day2", "g25_day3", "g25_day4", "g25_day5", "g25_day6")], list(na.rm = TRUE))),
    g30_day.max = do.call(pmax, c(data.13[, c("g30_day1", "g30_day2", "g30_day3", "g30_day4", "g30_day5", "g30_day6")], list(na.rm = TRUE))),
    g35_day.max = do.call(pmax, c(data.13[, c("g35_day1", "g35_day2", "g35_day3", "g35_day4", "g35_day5", "g35_day6")], list(na.rm = TRUE))),
    g40_day.max = do.call(pmax, c(data.13[, c("g40_day1", "g40_day2", "g40_day3", "g40_day4", "g40_day5", "g40_day6")], list(na.rm = TRUE))),
    g50_day.max = do.call(pmax, c(data.13[, c("g50_day1", "g50_day2", "g50_day3", "g50_day4", "g50_day5", "g50_day6")], list(na.rm = TRUE))),
    g60_day.max = do.call(pmax, c(data.13[, c("g60_day1", "g60_day2", "g60_day3", "g60_day4", "g60_day5", "g60_day6")], list(na.rm = TRUE))),
    g80_day.max = do.call(pmax, c(data.13[, c("g80_day1", "g80_day2", "g80_day3", "g80_day4", "g80_day5", "g80_day6")], list(na.rm = TRUE)))
  )

data.13 <- data.13 %>%
  # spray.[1-6]week.lag = was househhold i sprayed in the indicated previous week?
  # we don't count sprays on the same day because the entomological collection occurred in the morning before the spray
  mutate(
    # make an interval for the days 1-6 weeks before
    interval1 = interval(ymd(date - 7), ymd(date - 1)),
    interval2 = interval(ymd(date - 14), ymd(date - 8)),
    interval3 = interval(ymd(date - 21), ymd(date - 15)),
    interval4 = interval(ymd(date - 28), ymd(date - 22)),
    interval5 = interval(ymd(date - 35), ymd(date - 29)),
    interval6 = interval(ymd(date - 42), ymd(date - 36)),
    
    # if the date of a spray event falls within one of the intervals of a previous week, it assigns a value of 1
    spray.1week.lag = if_else(
      (cycle_1 %within% interval1 |
        cycle_2 %within% interval1 |
        cycle_3 %within% interval1 |
        cycle_4 %within% interval1 |
        cycle_5 %within% interval1 |
        cycle_6 %within% interval1) == TRUE, 1, 0
    ),
    spray.2week.lag = if_else(
      (cycle_1 %within% interval2 |
        cycle_2 %within% interval2 |
        cycle_3 %within% interval2 |
        cycle_4 %within% interval2 |
        cycle_5 %within% interval2 |
        cycle_6 %within% interval2) == TRUE, 1, 0
    ),
    spray.3week.lag = if_else(
      (cycle_1 %within% interval3 |
        cycle_2 %within% interval3 |
        cycle_3 %within% interval3 |
        cycle_4 %within% interval3 |
        cycle_5 %within% interval3 |
        cycle_6 %within% interval3) == TRUE, 1, 0
    ),
    spray.3week.lag = if_else(
      (cycle_1 %within% interval3 |
        cycle_2 %within% interval3 |
        cycle_3 %within% interval3 |
        cycle_4 %within% interval3 |
        cycle_5 %within% interval3 |
        cycle_6 %within% interval3) == TRUE, 1, 0
    ),
    spray.4week.lag = if_else(
      (cycle_1 %within% interval4 |
        cycle_2 %within% interval4 |
        cycle_3 %within% interval4 |
        cycle_4 %within% interval4 |
        cycle_5 %within% interval4 |
        cycle_6 %within% interval4) == TRUE, 1, 0
    ),
    spray.5week.lag = if_else(
      (cycle_1 %within% interval5 |
        cycle_2 %within% interval5 |
        cycle_3 %within% interval5 |
        cycle_4 %within% interval5 |
        cycle_5 %within% interval5 |
        cycle_6 %within% interval5) == TRUE, 1, 0
    ),
    spray.6week.lag = if_else(
      (cycle_1 %within% interval6 |
        cycle_2 %within% interval6 |
        cycle_3 %within% interval6 |
        cycle_4 %within% interval6 |
        cycle_5 %within% interval6 |
        cycle_6 %within% interval6) == TRUE, 1, 0
    )
  )


data.13<- data.13 %>%
  # remove redundant variables
  select(!ends_with("w1") & !ends_with("w2") & !ends_with("w3") & !ends_with("w4") & !ends_with("w5") & !ends_with("w6") &
    !ends_with("_day1") & !ends_with("_day2") & !ends_with("_day3") & !ends_with("_day4") & !ends_with("_day5") & !ends_with("_day6") &
    !starts_with("w_")) %>%
  # turn Nas into 0s
  mutate(across(ends_with("_day") | ends_with(".max") | ends_with(".lag") | ends_with(".min") | ends_with(".date"), ~ replace(., is.na(.), 0)))

# turn data.13 into an sf object again
data.13 <- st_as_sf(data.13, coords = c("lat", "long"), crs = st_crs(32718))

saveRDS(data.13, here("analysis", "data", "derived_data", "variables","variables_2013", "vbles.time.13.rds"))

