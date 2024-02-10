library(pacman)
pacman::p_load(tidyverse, sf, lubridate, parallel, here)


## ----data---------------------------------------------------------------------------------------------------------------------------------------------

vbles.time.14 <- readRDS(here("analysis", "data", "derived_data", "variables", "variables_2014", "vbles.time.14.rds"))

gis.buff.14 <- readRDS(here("analysis", "data", "raw_data", "gis.buff.sample.14.rds"))

# extract the data frame of study households from the gis.buff.14 data frame that includes all households within 1000m buffer of every study household
data.14 <- gis.buff.14 %>% 
  filter(!is.na(exp))

## ----biological plausibility rings ------------------------------------------------------------------------------------------------------------

cores= 2
  
  fx.perd.bio <- function (df, ncores) {
    
    fx.dist.w <- function (x) {
      # assign household i's geolocated point, location code and date of adult survey to objects 
      # to be used below
      center.geometry <- df[["geometry"]][x]
      center.location <- df[["location_code"]][x]
      center.date <- df[["date"]][[x]]
      center.moment <- df[["moment"]][[x]]
      
      df.w <- distinct(gis.buff.14 %>% 
                         select(c("location_code","cycle_1", "cycle_2", "cycle_3", "cycle_4", 
                                  "cycle_5", "cycle_6", "geometry"))) %>% 
        # dist = the distance (in m) between the household i and every other household j in the study, or dij.
        mutate(dist = round(as.numeric(st_distance(center.geometry, geometry)))) %>%
        # once we calculate the distance we don't need it to be an sf object anymore
        # for the following analyses it's easier to work with data frame rather than spatial objects
        st_drop_geometry() %>%
        # diff_c[1-6] = what is the difference in days between the entomological collection date and the date the household was sprayed?
        # date entomological survey - each of the spray dates in the 6 spray cycles
        
        # add 1 to all time differences because if the spray occurred on the same date as the entomological collection, the time difference is 0
        # and if we leave it as 0, it will not count as a spray weight
        
        # in L-2014, adult surveys were typically carried out 1 to 4 days after a house was sprayed (Gunning 2018).
        # diff_c[1-6]= 1 means the spray occurred on same day, diff_c[1-6] >1 means spray was in the past, diff_c[1-6] <1 means spray is in the future
        # to make sure sprays occurred before the adult collections don't have a spray effect weight in the model
        # we turn the values that are < 1 into NA
      mutate(
        diff_c1=ifelse((center.date - cycle_1+1) <1, NA, (center.date - cycle_1+1)),
        diff_c2=ifelse((center.date - cycle_2+1) <1, NA, (center.date - cycle_2+1)),
        diff_c3=ifelse((center.date - cycle_3+1) <1, NA, (center.date - cycle_3+1)),
        diff_c4=ifelse((center.date - cycle_4+1) <1, NA, (center.date - cycle_4+1)),
        diff_c5=ifelse((center.date - cycle_5+1) <1, NA, (center.date - cycle_5+1)),
        diff_c6=ifelse((center.date - cycle_6+1) <1, NA, (center.date - cycle_6+1))) %>% 
        # remove household i from the calculations (as we want to add the spray weights of the surrounding houses not counting sprays in household i)
        filter(location_code != center.location) %>% 
      # create distance rings around household i
        mutate(ring = case_when(
          between(dist, 0, 30) ~ "0-30m",
          between(dist, 31, 100) ~ "31-100m", 
          between(dist, 101, 300) ~ "101-300m", 
          between(dist, 301,1000) ~ ">300m"
          
        )) %>% 
        mutate(ring= factor(ring, levels=c("0-30m", "31-100m", "101-300m",">300m")))%>% 
        filter(!is.na(ring))
      
      df.w <- df.w %>% 
        # spray.past.w = was household j sprayed in the 7 days (1 week) previous to household i's adult survey time t?
        # we added +1 to all differences in time so 0 to 7 days before becomes 1 to 8 days before
        mutate(spray.past.w = if_any(c(diff_c1, diff_c2, diff_c3, diff_c4, diff_c5, diff_c6), ~ (between(.x, 1, 8)))) %>% 
        group_by(ring) %>% 
        
        mutate(
          #n.ring= total number of households in ring (including all households within 1000m buffer)
          n.ring= n(), 
          #n.s.ring= number of households in ring that were sprayed in the past 7 days
          n.s.ring= sum(spray.past.w ==TRUE, na.rm = TRUE)) %>% 
        
        group_by(ring) %>% 
        reframe(
          # perd.s.ring = proportion of households within a ring of a given distance from household i that had been sprayed in the week before t 
          # in other words, out of 10 households in each ring, how many were sprayed in the week before t?
          perd.s.ring= round(n.s.ring*10/n.ring, 2),
        )%>% 
        unique() %>% 
        pivot_wider(names_from = ring, names_prefix ="perd.s.ring_" , values_from = c( "perd.s.ring")) %>% 
        mutate(center.location= center.location,
               center.date=center.date
        )
    }
    
    dist.w <- lapply(1:nrow(df), fx.dist.w) %>% bind_rows()
    
    df.dist.w <- df %>% 
      left_join(dist.w , by=c("location_code" = "center.location", "date" = "center.date"))
    
  }


w.2014.bio <- fx.perd.bio(data.14,cores)
  

w.2014.bio <- w.2014.bio %>% 
  left_join(vbles.time.14[, c("location_code", "month","date","g20_day.max")] %>% st_drop_geometry(), by= c("location_code", "date"))

saveRDS(w.2014.bio , here("analysis", "data", "derived_data", "variables", "variables_2014", "vbles.space.ring.bio.14.rds"))



## ----evenly spaced 100m interval rings ------------------------------------------------------------------------------------------------------------


fx.perd.100m <- function (df, ncores) {
  
  fx.dist.w <- function (x) {
    
    center.geometry <- df[["geometry"]][x]
    center.location <- df[["location_code"]][x]
    center.date <- df[["date"]][[x]]
    center.moment <- df[["moment"]][[x]]
    
    df.w <- distinct(gis.buff.14 %>% 
                       select(c("location_code","cycle_1", "cycle_2", "cycle_3", "cycle_4", 
                                "cycle_5", "cycle_6", "geometry"))) %>% 
      # dist = the distance (in m) between the household i and every other household j in the study, or dij.
      mutate(dist = round(as.numeric(st_distance(center.geometry, geometry)))) %>%
      # once we calculate the distance we don't need it to be an sf object anymore
      # for the following analyses it's easier to work with data frame rather than spatial objects
      st_drop_geometry() %>%
      # diff_c[1-6] = what is the difference in days between the entomological collection date and the date the household was sprayed?
      # date entomological survey - each of the spray dates in the 6 spray cycles
      
      # add 1 to all time differences because if the spray occurred on the same date as the entomological collection, the time difference is 0
      # and if we leave it as 0, it will not count as a spray weight
      
      # in L-2014, adult surveys were typically carried out 1 to 4 days after a house was sprayed (Gunning 2018).
      # diff_c[1-6]= 1 means the spray occurred on same day, diff_c[1-6] >1 means spray was in the past, diff_c[1-6] <1 means spray is in the future
      # to make sure sprays occurred before the adult collections don't have a spray effect weight in the model
      # we turn the values that are < 1 into NA
      mutate(
        diff_c1=ifelse((center.date - cycle_1+1) <1, NA, (center.date - cycle_1+1)),
        diff_c2=ifelse((center.date - cycle_2+1) <1, NA, (center.date - cycle_2+1)),
        diff_c3=ifelse((center.date - cycle_3+1) <1, NA, (center.date - cycle_3+1)),
        diff_c4=ifelse((center.date - cycle_4+1) <1, NA, (center.date - cycle_4+1)),
        diff_c5=ifelse((center.date - cycle_5+1) <1, NA, (center.date - cycle_5+1)),
        diff_c6=ifelse((center.date - cycle_6+1) <1, NA, (center.date - cycle_6+1))) %>% 
      # remove household i from the calculations (as we want to add the spray weights of the surrounding houses not counting sprays in household i)
      filter(location_code != center.location) %>% 
      #filter out the center house
      mutate(ring = case_when(
        #evenly spaced rings by 1000m, from 1-500m 
        between(dist, 0, 100) ~ "0-100m", 
        between(dist, 101, 200) ~ "101-200m", 
        between(dist, 201, 300) ~ "201-300m", 
        between(dist, 301, 400) ~ "301-400m", 
        between(dist, 401, 500) ~ "401-500m",
        between(dist, 501,1000) ~ ">500m"
      )) %>% 
      mutate(ring= factor(ring, levels=c("0-100m", "101-200m", "201-300m", "301-400m","401-500m", ">500m"))) %>% 
      filter(!is.na(ring))
    
    df.w <- df.w %>% 
      # spray.past.w = was household j sprayed in the 7 days (1 week) previous to household i's adult survey time t?
      # we added +1 to all differences in time so 0 to 7 days before becomes 1 to 8 days before
      mutate(spray.past.w = if_any(c(diff_c1, diff_c2, diff_c3, diff_c4, diff_c5, diff_c6), ~ (between(.x, 1, 8)))) %>% 
      group_by(ring) %>% 
      
      mutate(
        #n.ring= total number of households in ring (including all households within 1000m buffer)
        n.ring= n(), 
        #n.s.ring= number of households in ring that were sprayed in the past 7 days
        n.s.ring= sum(spray.past.w ==TRUE, na.rm = TRUE)) %>% 
      
      group_by(ring) %>% 
      reframe(
        # perd.s.ring = proportion of households within a ring of a given distance from household i that had been sprayed in the week before t 
        # in other words, out of 10 households in each ring, how many were sprayed in the week before t?
        perd.s.ring= round(n.s.ring*10/n.ring, 2),
      )%>% 
      unique() %>% 
      pivot_wider(names_from = ring, names_prefix ="perd.s.ring_" , values_from = c( "perd.s.ring")) %>% 
      mutate(center.location= center.location,
             center.date=center.date
      )
  }
  
  dist.w <- lapply(1:nrow(df), fx.dist.w) %>% bind_rows()
  
  df.dist.w <- df %>% 
    left_join(dist.w , by=c("location_code" = "center.location", "date" = "center.date"))
  
}


w.2014.100m <- fx.perd.100m(data.14,cores)


w.2014.100m<- w.2014.100m %>% 
  left_join(vbles.time.14[, c("location_code","month", "date","g20_day.max")] %>% st_drop_geometry(), by= c("location_code", "date")) 

saveRDS(w.2014.100m , here("analysis", "data", "derived_data", "variables", "variables_2014", "vbles.space.ring.100m.14.rds"))


