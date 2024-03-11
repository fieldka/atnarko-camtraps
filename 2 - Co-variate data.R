
############################################
# Script 2 - Atnarko Camera Traps
# Join co-variate data:
# Salmon, berries, land-humans, water level
############################################

# setwd() 

# Load Packages
list.of.packages <- c("lubridate","dplyr", "ggplot2")
# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

## Salmon 
salm <- read.csv("salmon drift counts 2020-2021.csv")

# Join salmon to cam dataset and pull in coordinates. 
salm <- left_join(cam, salm, by = "riv.section.weekid")
sta_cords <- dplyr::select(sta, Latitude, Longitude, Deployment.Location.ID)
salm <- left_join(salm, sta_cords, by = "Deployment.Location.ID")
salm$date_of_video <- as.Date(salm$date_of_video)

# Create a dataframe that is only yearweeks, siteweek_IDs, riv.section.weekid, coordinates, and salmon.
salm <- salm[, c('year', 'week.x', 'yearweek_id.x','site_name', 'siteweek_id', "Latitude", "Longitude", "Sock_AL_EstTotal", "Coho_AL_EstTotal", "Pink_AL_EstTotal", "Chum_AL_EstTotal", "Chinook_AL_EstTotal", "Est_all_species", "River_section_ID.x", "riv.section.weekid")]
salm$yearweek_id <- salm$yearweek_id.x

# Collapse the dataset to have a row for every unique siteweek_id.
collapse <- salm[!duplicated(salm$siteweek_id),]
collapse <- collapse[order(collapse$yearweek_id),]  

# Data for distance from camera site to river mouth
dist <- read.csv("dist from cam site to river mouth.csv")

# Divide each river segment length by 2 for center point. Then, extract the camera trap location that falls closest to the center point of each river segment. 

# Belarko drift: Belarko site to talchacko confluence
# River segment length: 4.4km

# Mid: Cynthia's to Belarko site: 
23.3-4.4 #subtract the distance to Belarko site from the distance to Atnarko campground
# River segment length: 18.9km

# Talus Drift: Talus site to Cynthia's
28.7-23.3
# River segment length: 5.4km

# Above Talus Drift: Stillwater Forest to Talus
35.5-28.7
# River segment length: 6.8km

# Result:
# Belarko: Big Tweeds
# Mid: Fish Bone
# Talus: Sockeye island
# Above Talus: Line Cabin

collapse <- collapse[!duplicated(collapse$riv.section.weekid),]
collapse <- collapse[, c('year', 'week.x', 'yearweek_id',"Sock_AL_EstTotal", "Coho_AL_EstTotal", "Pink_AL_EstTotal", "Chum_AL_EstTotal", "Chinook_AL_EstTotal", "Est_all_species", "riv.section.weekid", "River_section_ID.x")]

# Assign distances from camera trap mid-points to each River_section_ID. These are the distance from Big Tweeds, Fish Bone, Sockeye Island, and Line Cabin to the confluence. 

collapse <- collapse %>%
  mutate(dist.along.riv.to.confluence =
           case_when(River_section_ID.x == "Belarko Drift" ~ 2.34,
                     River_section_ID.x == "Mid" ~ 13.07,
                     River_section_ID.x == "Talus Drift" ~ 26.77,
                     River_section_ID.x == "Above Talus Drift" ~ 31.53))

collapse <- collapse %>%
  mutate(order_riv_sections =
           case_when(River_section_ID.x == "Belarko Drift" ~ "A",
                     River_section_ID.x == "Mid" ~ "B",
                     River_section_ID.x == "Talus Drift" ~ "C",
                     River_section_ID.x == "Above Talus Drift" ~ "D"))

collapse <- collapse %>%  
  group_by(yearweek_id) %>% 
  arrange((order_riv_sections), .by_group = TRUE)

# Impute (2020 and 2021)
# Three steps:
# 1) Temporally impute for weeks when Belarko and/or Talus counts are missing (as described in methods)
# 2) Spatially impute for the mid section of the river (as described in methods)
# 3) Spatially impute for the above section of the river (as described in methods)

# 1) Temporally impute first for weeks when there is only one count done in Belarko and/or Talus. If data are missing for 1 week, add the two, within-river segment values from the bracketing weeks, and divide them by 2. E.g., 202037 (salmon value from 202036 + salmon value from 202038)/2. This was done manually in excel. Note: 202036 mid-upper river segment (Talus drift) data were imputed, and raw data were not used as per poor counting conditions (pers. comm. with Nuxalk Fisheries and Wildlife staff).  

collapse <- read.csv("temporally imputed_salmon drift counts 2020-2021.csv")

# 2) Spatially impute values using weighted averaging. 

# All species
impute <- collapse %>%
  arrange(yearweek_id, order_riv_sections) %>%
  group_by(yearweek_id) %>%
  mutate(
    Est_all_species = ifelse(
      is.na(Est_all_species) & River_section_ID.x == "Mid", 
      ((lag(Est_all_species) * (1/(dist.along.riv.to.confluence - lag(dist.along.riv.to.confluence)))) + 
         (lead(Est_all_species) * (1/(lead(dist.along.riv.to.confluence) - dist.along.riv.to.confluence)))) / ((1/(dist.along.riv.to.confluence - lag(dist.along.riv.to.confluence))) + ((1/(lead(dist.along.riv.to.confluence) - dist.along.riv.to.confluence)))), Est_all_species)
  )

# Sockeye
impute <- impute %>%
  arrange(yearweek_id, order_riv_sections) %>%
  group_by(yearweek_id) %>%
  mutate(
    Sock_AL_EstTotal = ifelse(
      is.na(Sock_AL_EstTotal) & River_section_ID.x == "Mid", 
      ((lag(Sock_AL_EstTotal) * (1/(dist.along.riv.to.confluence - lag(dist.along.riv.to.confluence)))) + 
         (lead(Sock_AL_EstTotal) * (1/(lead(dist.along.riv.to.confluence) - dist.along.riv.to.confluence)))) / ((1/(dist.along.riv.to.confluence - lag(dist.along.riv.to.confluence))) + ((1/(lead(dist.along.riv.to.confluence) - dist.along.riv.to.confluence)))), Sock_AL_EstTotal)
  )

# Coho
impute <- impute %>%
  arrange(yearweek_id, order_riv_sections) %>%
  group_by(yearweek_id) %>%
  mutate(
    Coho_AL_EstTotal = ifelse(
      is.na(Coho_AL_EstTotal) & River_section_ID.x == "Mid", 
      ((lag(Coho_AL_EstTotal) * (1/(dist.along.riv.to.confluence - lag(dist.along.riv.to.confluence)))) + 
         (lead(Coho_AL_EstTotal) * (1/(lead(dist.along.riv.to.confluence) - dist.along.riv.to.confluence)))) / ((1/(dist.along.riv.to.confluence - lag(dist.along.riv.to.confluence))) + ((1/(lead(dist.along.riv.to.confluence) - dist.along.riv.to.confluence)))), Coho_AL_EstTotal)
  )

# Pink
impute <- impute %>%
  arrange(yearweek_id, order_riv_sections) %>%
  group_by(yearweek_id) %>%
  mutate(
    Pink_AL_EstTotal = ifelse(
      is.na(Pink_AL_EstTotal) & River_section_ID.x == "Mid", 
      ((lag(Pink_AL_EstTotal) * (1/(dist.along.riv.to.confluence - lag(dist.along.riv.to.confluence)))) + 
         (lead(Pink_AL_EstTotal) * (1/(lead(dist.along.riv.to.confluence) - dist.along.riv.to.confluence)))) / ((1/(dist.along.riv.to.confluence - lag(dist.along.riv.to.confluence))) + ((1/(lead(dist.along.riv.to.confluence) - dist.along.riv.to.confluence)))), Pink_AL_EstTotal)
  )

# Pink in scientific notation.
impute$Pink_AL_EstTotal <- sprintf("%.3f", impute$Pink_AL_EstTotal)
impute$Pink_AL_EstTotal <- as.integer(impute$Pink_AL_EstTotal)
impute$Pink_AL_EstTotal <- round(impute$Pink_AL_EstTotal, digits = 0)

# Chum
impute <- impute %>%
  arrange(yearweek_id, order_riv_sections) %>%
  group_by(yearweek_id) %>%
  mutate(
    Chum_AL_EstTotal = ifelse(
      is.na(Chum_AL_EstTotal) & River_section_ID.x == "Mid", 
      ((lag(Chum_AL_EstTotal) * (1/(dist.along.riv.to.confluence - lag(dist.along.riv.to.confluence)))) + 
         (lead(Chum_AL_EstTotal) * (1/(lead(dist.along.riv.to.confluence) - dist.along.riv.to.confluence)))) / ((1/(dist.along.riv.to.confluence - lag(dist.along.riv.to.confluence))) + ((1/(lead(dist.along.riv.to.confluence) - dist.along.riv.to.confluence)))), Chum_AL_EstTotal)
  )

# Chinook
impute <- impute %>%
  arrange(yearweek_id, order_riv_sections) %>%
  group_by(yearweek_id) %>%
  mutate(
    Chinook_AL_EstTotal = ifelse(
      is.na(Chinook_AL_EstTotal) & River_section_ID.x == "Mid", 
      ((lag(Chinook_AL_EstTotal) * (1/(dist.along.riv.to.confluence - lag(dist.along.riv.to.confluence)))) + 
         (lead(Chinook_AL_EstTotal) * (1/(lead(dist.along.riv.to.confluence) - dist.along.riv.to.confluence)))) / ((1/(dist.along.riv.to.confluence - lag(dist.along.riv.to.confluence))) + ((1/(lead(dist.along.riv.to.confluence) - dist.along.riv.to.confluence)))), Chinook_AL_EstTotal)
  )

# 3) Spatially impute for the above Talus section as described in methods. 

# All species
impute_weekID <- function(River_section_ID.x, Est_all_species, dist.along.riv.to.confluence) {
  diff_belmid_salm <- Est_all_species[River_section_ID.x == 'Belarko Drift'] - Est_all_species[River_section_ID.x == 'Mid']
  diff_midtal_salm <- Est_all_species[River_section_ID.x == 'Mid'] - Est_all_species[River_section_ID.x == 'Talus Drift']
  
  diff_belmid_dist <- dist.along.riv.to.confluence[River_section_ID.x == 'Belarko Drift'] - dist.along.riv.to.confluence[River_section_ID.x == 'Mid']
  diff_midtal_dist <- dist.along.riv.to.confluence[River_section_ID.x == 'Mid'] - dist.along.riv.to.confluence[River_section_ID.x == 'Talus Drift']
  diff_talabove_dist <- dist.along.riv.to.confluence[River_section_ID.x == 'Talus Drift'] - dist.along.riv.to.confluence[River_section_ID.x == 'Above Talus Drift']
  
  normalized_diff_belmid <- abs(diff_belmid_salm / diff_belmid_dist)
  normalized_diff_midtal <- abs(diff_midtal_salm / diff_midtal_dist)
  
  avg_normalized_diff <- (normalized_diff_belmid + normalized_diff_midtal) / 2
  
  Est_all_species_shift <- abs(avg_normalized_diff * diff_talabove_dist)
  
  if (Est_all_species[River_section_ID.x == 'Belarko Drift'] < Est_all_species[River_section_ID.x == 'Talus Drift']) {
    return(Est_all_species[River_section_ID.x == 'Talus Drift'] + Est_all_species_shift)
  }
  return(Est_all_species[River_section_ID.x == 'Talus Drift'] - Est_all_species_shift)
}

impute_subset <- impute %>%
  group_by(yearweek_id) %>%
  filter(sum(is.na(Est_all_species))<2)

all_imputed_values <-
  impute_subset %>%
  group_by(yearweek_id) %>%
  summarize(imputed = impute_weekID(River_section_ID.x, Est_all_species, dist.along.riv.to.confluence)) 

# Substitute the entries of all_imputed_values for the site D entries in the original data frame:
impute_subset$Est_all_species[impute_subset$River_section_ID.x == 'Above Talus Drift'] <- all_imputed_values$imputed
column <- c("Est_all_species")
replace_w_zero <- function(x) {
  ifelse(x < 0, yes = 0, no = x)
}

impute_subset <- mutate_at(
  impute_subset,
  column,
  replace_w_zero
)

# Sockeye
impute_weekID <- function(River_section_ID.x, Sock_AL_EstTotal, dist.along.riv.to.confluence) {
  diff_belmid_salm <- Sock_AL_EstTotal[River_section_ID.x == 'Belarko Drift'] - Sock_AL_EstTotal[River_section_ID.x == 'Mid']
  diff_midtal_salm <- Sock_AL_EstTotal[River_section_ID.x == 'Mid'] - Sock_AL_EstTotal[River_section_ID.x == 'Talus Drift']
  
  diff_belmid_dist <- dist.along.riv.to.confluence[River_section_ID.x == 'Belarko Drift'] - dist.along.riv.to.confluence[River_section_ID.x == 'Mid']
  diff_midtal_dist <- dist.along.riv.to.confluence[River_section_ID.x == 'Mid'] - dist.along.riv.to.confluence[River_section_ID.x == 'Talus Drift']
  diff_talabove_dist <- dist.along.riv.to.confluence[River_section_ID.x == 'Talus Drift'] - dist.along.riv.to.confluence[River_section_ID.x == 'Above Talus Drift']
  
  normalized_diff_belmid <- abs(diff_belmid_salm / diff_belmid_dist)
  normalized_diff_midtal <- abs(diff_midtal_salm / diff_midtal_dist)
  
  avg_normalized_diff <- (normalized_diff_belmid + normalized_diff_midtal) / 2
  
  Sock_AL_EstTotal_shift <- abs(avg_normalized_diff * diff_talabove_dist)
  
  if (Sock_AL_EstTotal[River_section_ID.x == 'Belarko Drift'] < Sock_AL_EstTotal[River_section_ID.x == 'Talus Drift']) {
    return(Sock_AL_EstTotal[River_section_ID.x == 'Talus Drift'] + Sock_AL_EstTotal_shift)
  }
  return(Sock_AL_EstTotal[River_section_ID.x == 'Talus Drift'] - Sock_AL_EstTotal_shift)
}

impute_subset <- impute_subset %>%
  group_by(yearweek_id) %>%
  filter(sum(is.na(Sock_AL_EstTotal))<2)

all_imputed_values <-
  impute_subset %>%
  group_by(yearweek_id) %>%
  summarize(imputed = impute_weekID(River_section_ID.x, Sock_AL_EstTotal, dist.along.riv.to.confluence)) 

# Substitute the entries of all_imputed_values for the site D entries in the original data frame:
impute_subset$Sock_AL_EstTotal[impute_subset$River_section_ID.x == 'Above Talus Drift'] <- all_imputed_values$imputed
column <- c("Sock_AL_EstTotal")
replace_w_zero <- function(x) {
  ifelse(x < 0, yes = 0, no = x)
}

impute_subset <- mutate_at(
  impute_subset,
  column,
  replace_w_zero
)

# Coho
impute_weekID <- function(River_section_ID.x, Coho_AL_EstTotal, dist.along.riv.to.confluence) {
  diff_belmid_salm <- Coho_AL_EstTotal[River_section_ID.x == 'Belarko Drift'] - Coho_AL_EstTotal[River_section_ID.x == 'Mid']
  diff_midtal_salm <- Coho_AL_EstTotal[River_section_ID.x == 'Mid'] - Coho_AL_EstTotal[River_section_ID.x == 'Talus Drift']
  
  diff_belmid_dist <- dist.along.riv.to.confluence[River_section_ID.x == 'Belarko Drift'] - dist.along.riv.to.confluence[River_section_ID.x == 'Mid']
  diff_midtal_dist <- dist.along.riv.to.confluence[River_section_ID.x == 'Mid'] - dist.along.riv.to.confluence[River_section_ID.x == 'Talus Drift']
  diff_talabove_dist <- dist.along.riv.to.confluence[River_section_ID.x == 'Talus Drift'] - dist.along.riv.to.confluence[River_section_ID.x == 'Above Talus Drift']
  
  normalized_diff_belmid <- abs(diff_belmid_salm / diff_belmid_dist)
  normalized_diff_midtal <- abs(diff_midtal_salm / diff_midtal_dist)
  
  avg_normalized_diff <- (normalized_diff_belmid + normalized_diff_midtal) / 2
  
  Coho_AL_EstTotal_shift <- abs(avg_normalized_diff * diff_talabove_dist)
  
  if (Coho_AL_EstTotal[River_section_ID.x == 'Belarko Drift'] < Coho_AL_EstTotal[River_section_ID.x == 'Talus Drift']) {
    return(Coho_AL_EstTotal[River_section_ID.x == 'Talus Drift'] + Coho_AL_EstTotal_shift)
  }
  return(Coho_AL_EstTotal[River_section_ID.x == 'Talus Drift'] - Coho_AL_EstTotal_shift)
}

impute_subset <- impute_subset %>%
  group_by(yearweek_id) %>%
  filter(sum(is.na(Coho_AL_EstTotal))<2)

all_imputed_values <-
  impute_subset %>%
  group_by(yearweek_id) %>%
  summarize(imputed = impute_weekID(River_section_ID.x, Coho_AL_EstTotal, dist.along.riv.to.confluence)) 

# Substitute the entries of all_imputed_values for the site D entries in the original data frame:
impute_subset$Coho_AL_EstTotal[impute_subset$River_section_ID.x == 'Above Talus Drift'] <- all_imputed_values$imputed
column <- c("Coho_AL_EstTotal")
replace_w_zero <- function(x) {
  ifelse(x < 0, yes = 0, no = x)
}

impute_subset <- mutate_at(
  impute_subset,
  column,
  replace_w_zero
)

# Pink
impute_weekID <- function(River_section_ID.x, Pink_AL_EstTotal, dist.along.riv.to.confluence) {
  diff_belmid_salm <- Pink_AL_EstTotal[River_section_ID.x == 'Belarko Drift'] - Pink_AL_EstTotal[River_section_ID.x == 'Mid']
  diff_midtal_salm <- Pink_AL_EstTotal[River_section_ID.x == 'Mid'] - Pink_AL_EstTotal[River_section_ID.x == 'Talus Drift']
  
  diff_belmid_dist <- dist.along.riv.to.confluence[River_section_ID.x == 'Belarko Drift'] - dist.along.riv.to.confluence[River_section_ID.x == 'Mid']
  diff_midtal_dist <- dist.along.riv.to.confluence[River_section_ID.x == 'Mid'] - dist.along.riv.to.confluence[River_section_ID.x == 'Talus Drift']
  diff_talabove_dist <- dist.along.riv.to.confluence[River_section_ID.x == 'Talus Drift'] - dist.along.riv.to.confluence[River_section_ID.x == 'Above Talus Drift']
  
  normalized_diff_belmid <- abs(diff_belmid_salm / diff_belmid_dist)
  normalized_diff_midtal <- abs(diff_midtal_salm / diff_midtal_dist)
  
  avg_normalized_diff <- (normalized_diff_belmid + normalized_diff_midtal) / 2
  
  Pink_AL_EstTotal_shift <- abs(avg_normalized_diff * diff_talabove_dist)
  
  if (Pink_AL_EstTotal[River_section_ID.x == 'Belarko Drift'] < Pink_AL_EstTotal[River_section_ID.x == 'Talus Drift']) {
    return(Pink_AL_EstTotal[River_section_ID.x == 'Talus Drift'] + Pink_AL_EstTotal_shift)
  }
  return(Pink_AL_EstTotal[River_section_ID.x == 'Talus Drift'] - Pink_AL_EstTotal_shift)
}

impute_subset <- impute_subset %>%
  group_by(yearweek_id) %>%
  filter(sum(is.na(Pink_AL_EstTotal))<2)

all_imputed_values <-
  impute_subset %>%
  group_by(yearweek_id) %>%
  summarize(imputed = impute_weekID(River_section_ID.x, Pink_AL_EstTotal, dist.along.riv.to.confluence)) 

# Substitute the entries of all_imputed_values for the site D entries in the original data frame:
impute_subset$Pink_AL_EstTotal[impute_subset$River_section_ID.x == 'Above Talus Drift'] <- all_imputed_values$imputed
column <- c("Pink_AL_EstTotal")
replace_w_zero <- function(x) {
  ifelse(x < 0, yes = 0, no = x)
}

impute_subset <- mutate_at(
  impute_subset,
  column,
  replace_w_zero
)

# Pink in scientific notation. Change this then check the values are the same
impute_subset$Pink_AL_EstTotal <- sprintf("%.3f", impute_subset$Pink_AL_EstTotal)
impute_subset$Pink_AL_EstTotal <- as.integer(impute_subset$Pink_AL_EstTotal)
impute_subset$Pink_AL_EstTotal <- round(impute_subset$Pink_AL_EstTotal, digits = 0)

# Chum
impute_weekID <- function(River_section_ID.x, Chum_AL_EstTotal, dist.along.riv.to.confluence) {
  diff_belmid_salm <- Chum_AL_EstTotal[River_section_ID.x == 'Belarko Drift'] - Chum_AL_EstTotal[River_section_ID.x == 'Mid']
  diff_midtal_salm <- Chum_AL_EstTotal[River_section_ID.x == 'Mid'] - Chum_AL_EstTotal[River_section_ID.x == 'Talus Drift']
  
  diff_belmid_dist <- dist.along.riv.to.confluence[River_section_ID.x == 'Belarko Drift'] - dist.along.riv.to.confluence[River_section_ID.x == 'Mid']
  diff_midtal_dist <- dist.along.riv.to.confluence[River_section_ID.x == 'Mid'] - dist.along.riv.to.confluence[River_section_ID.x == 'Talus Drift']
  diff_talabove_dist <- dist.along.riv.to.confluence[River_section_ID.x == 'Talus Drift'] - dist.along.riv.to.confluence[River_section_ID.x == 'Above Talus Drift']
  
  normalized_diff_belmid <- abs(diff_belmid_salm / diff_belmid_dist)
  normalized_diff_midtal <- abs(diff_midtal_salm / diff_midtal_dist)
  
  avg_normalized_diff <- (normalized_diff_belmid + normalized_diff_midtal) / 2
  
  Chum_AL_EstTotal_shift <- abs(avg_normalized_diff * diff_talabove_dist)
  
  if (Chum_AL_EstTotal[River_section_ID.x == 'Belarko Drift'] < Chum_AL_EstTotal[River_section_ID.x == 'Talus Drift']) {
    return(Chum_AL_EstTotal[River_section_ID.x == 'Talus Drift'] + Chum_AL_EstTotal_shift)
  }
  return(Chum_AL_EstTotal[River_section_ID.x == 'Talus Drift'] - Chum_AL_EstTotal_shift)
}

impute_subset <- impute_subset %>%
  group_by(yearweek_id) %>%
  filter(sum(is.na(Chum_AL_EstTotal))<2)

all_imputed_values <-
  impute_subset %>%
  group_by(yearweek_id) %>%
  summarize(imputed = impute_weekID(River_section_ID.x, Chum_AL_EstTotal, dist.along.riv.to.confluence)) 

# Substitute the entries of all_imputed_values for the site D entries in the original data frame:
impute_subset$Chum_AL_EstTotal[impute_subset$River_section_ID.x == 'Above Talus Drift'] <- all_imputed_values$imputed
column <- c("Chum_AL_EstTotal")
replace_w_zero <- function(x) {
  ifelse(x < 0, yes = 0, no = x)
}

impute_subset <- mutate_at(
  impute_subset,
  column,
  replace_w_zero
)

# Chinook
impute_weekID <- function(River_section_ID.x, Chinook_AL_EstTotal, dist.along.riv.to.confluence) {
  diff_belmid_salm <- Chinook_AL_EstTotal[River_section_ID.x == 'Belarko Drift'] - Chinook_AL_EstTotal[River_section_ID.x == 'Mid']
  diff_midtal_salm <- Chinook_AL_EstTotal[River_section_ID.x == 'Mid'] - Chinook_AL_EstTotal[River_section_ID.x == 'Talus Drift']
  
  diff_belmid_dist <- dist.along.riv.to.confluence[River_section_ID.x == 'Belarko Drift'] - dist.along.riv.to.confluence[River_section_ID.x == 'Mid']
  diff_midtal_dist <- dist.along.riv.to.confluence[River_section_ID.x == 'Mid'] - dist.along.riv.to.confluence[River_section_ID.x == 'Talus Drift']
  diff_talabove_dist <- dist.along.riv.to.confluence[River_section_ID.x == 'Talus Drift'] - dist.along.riv.to.confluence[River_section_ID.x == 'Above Talus Drift']
  
  normalized_diff_belmid <- abs(diff_belmid_salm / diff_belmid_dist)
  normalized_diff_midtal <- abs(diff_midtal_salm / diff_midtal_dist)
  
  avg_normalized_diff <- (normalized_diff_belmid + normalized_diff_midtal) / 2
  
  Chinook_AL_EstTotal_shift <- abs(avg_normalized_diff * diff_talabove_dist)
  
  if (Chinook_AL_EstTotal[River_section_ID.x == 'Belarko Drift'] < Chinook_AL_EstTotal[River_section_ID.x == 'Talus Drift']) {
    return(Chinook_AL_EstTotal[River_section_ID.x == 'Talus Drift'] + Chinook_AL_EstTotal_shift)
  }
  return(Chinook_AL_EstTotal[River_section_ID.x == 'Talus Drift'] - Chinook_AL_EstTotal_shift)
}

impute_subset <- impute_subset %>%
  group_by(yearweek_id) %>%
  filter(sum(is.na(Chinook_AL_EstTotal))<2)

all_imputed_values <-
  impute_subset %>%
  group_by(yearweek_id) %>%
  summarize(imputed = impute_weekID(River_section_ID.x, Chinook_AL_EstTotal, dist.along.riv.to.confluence)) 

# Substitute the entries of all_imputed_values for the site D entries in the original data frame:
impute_subset$Chinook_AL_EstTotal[impute_subset$River_section_ID.x == 'Above Talus Drift'] <- all_imputed_values$imputed
column <- c("Chinook_AL_EstTotal")
replace_w_zero <- function(x) {
  ifelse(x < 0, yes = 0, no = x)
}

impute_subset <- mutate_at(
  impute_subset,
  column,
  replace_w_zero
)

# Now we have a dataframe with imputed values for weeks when both Bel and Tal counts were done. 
# Merge imputed DF to impute_subset DF now that we have imputed Above Talus Drift Values

impute <- left_join(impute, impute_subset, by = "riv.section.weekid")

impute <- impute[, c('year.x', "yearweek_id.x","Sock_AL_EstTotal.y", "Coho_AL_EstTotal.y", "Pink_AL_EstTotal.y", "Chum_AL_EstTotal.y", "Chinook_AL_EstTotal.y", "River_section_ID.x.x", "dist.along.riv.to.confluence.x","riv.section.weekid", "order_riv_sections.x")]

impute$River_section_ID.x.x <- factor(impute$River_section_ID.x.x, levels = c("Belarko Drift", "Mid", "Talus Drift", "Above Talus Drift"))

# Plots for salmon counts and biomass
# Add a week column
n_last <- 2
impute$week <- substr(impute$yearweek_id.x, nchar(impute$yearweek_id.x) - n_last + 1, nchar(impute$yearweek_id.x))

# Biomass calculations. For each species, we used average mass reported in Groot and Margolis (1991), and assumed a 1:1 sex ratio.

# Subset pink, sockeye, and Chinook
impute <- impute %>%
  dplyr::select(-Coho_AL_EstTotal.y, -Chum_AL_EstTotal.y)
impute <- impute %>%
  mutate(sock.biomass = impute$Sock_AL_EstTotal.y * 2.7,
         spring.biomass = impute$Chinook_AL_EstTotal.y * 13.6)

# Different biomass values for pinks depending on the year
impute <- impute %>%
  mutate(pink.biomass =
           case_when(year.x == "2020" ~ impute$Pink_AL_EstTotal.y * 1.7,
                     year.x == "2021" ~ impute$Pink_AL_EstTotal.y * 2.5))

final <- impute %>%
  mutate(total.biomass = rowSums(impute[,c("sock.biomass",
                                           "spring.biomass",
                                           "pink.biomass")]))

# Impute for 2019
# Calculate the NuSEDS to count ratio as described in methods. 
salm2020 <- filter(final, final$year.x == "2020")

# Sum the counts
salm2020.sum <- salm2020 %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

# 2020 drift count data by species (total within the year). 
# Sockeye: 164.4321
# Pink: 358382
# Chinook: 5782.568
# These include counts that occurred between Belarko to Fisheries Pool, between Talus and Cynthia's pool, and counts that were imputed using the average ratio approach above. There are decimal points because imputed values from above script are included. 

# 2020 NuSEDS estimates (Method: Expert Opinion)
# Dataset: open data available from https://open.canada.ca/data/en/dataset/c48669a3-045b-400d-b730-48aafe8c5ee6
# Sockeye: 8427
# Pink: 458612
# Chinook: 19644

# 2020 NuSEDS to 2020 Drift Count ratios

# Sockeye:
164.4321/8427
# 0.01951253

# Pink: 
358382/458612
# 0.7814492

# Chinook:
5782.568/19644
# 0.2943682

# Estimate 2019 count data for each species using these ratios

# 2019 NuSEDS data
# Sockeye: 4598
# Pink: 20550
# Chinook: 12017

# Estimate 2019 count data for each species using these ratios

# Sockeye: 
4598 * 0.1951253
# 897.1861

# Pink
20550 * 0.7814492
# 16058.78

# Chinook
12017 * 0.2943682
# 3537.423

# Calculate the proportion of counts that occurred in a given site and in a given week for each species. 
library(tidyverse)

# Same site, different weeks
# E.g., in week 1 at Belarko, sockeye counts were x% of their total count across all weeks at Belarko. 
d <- final
d <- na.omit(d)
d <- d %>%
  dplyr::select(year.x, yearweek_id.x, Sock_AL_EstTotal.y, Pink_AL_EstTotal.y, Chinook_AL_EstTotal.y, River_section_ID.x.x)

d <- d %>% 
  group_by(year.x, River_section_ID.x.x) %>% 
  mutate(across(Sock_AL_EstTotal.y:Chinook_AL_EstTotal.y, 
                ~ prop.table (.x),
                .names = "{.col}_prop"))

n_last <- 2
d$week <- substr(d$yearweek_id.x, nchar(d$yearweek_id.x) - n_last + 1, nchar(d$yearweek_id.x))

d <- d %>%
  group_by(week,River_section_ID.x.x) %>%
  mutate(mean.sock = mean(Sock_AL_EstTotal.y_prop)) %>%
  mutate(mean.pink = mean(Pink_AL_EstTotal.y_prop)) %>%
  mutate(mean.chin = mean(Chinook_AL_EstTotal.y_prop))

# Take the estimated count data for 2019, and use the proportions at each siteweek to impute for 2019. 
# Create a site-week ID, excluding year. 
# For every site-week, there will be a mean proportion for each species. 
# Multiply the estimated 2019 counts by this proportion to impute for each week. 
d$riversection_weekID <- paste(d$River_section_ID.x.x, d$week, sep = "")

# Proportions should add up to one across sites, within years. 
test <- subset(d, River_section_ID.x.x == "Talus Drift")
test <- subset(test, year.x == "2020")
sum(test$Sock_AL_EstTotal.y_prop)

d <- d %>%
  dplyr::select(mean.chin, mean.sock, mean.pink, riversection_weekID)

# create columns for estimated 2019 counts.
d$sock2019 <- 897.1861
d$pink2019 <- 16058.78
d$chin2019 <- 3537.423
d$year <- 2019

# Remove duplicates
d <- d[!duplicated(d$riversection_weekID),]

# Impute by multiplying the total counts for 2019 by the mean proportion
d$imputed.sock2019 <- d$mean.sock * d$sock2019
d$imputed.pink2019 <- d$mean.pink * d$pink2019
d$imputed.chin2019 <- d$mean.chin * d$chin2019
d$yearweek_id.x <- paste(d$year, d$week, sep = "")
d$riversection_weekid <- paste(d$River_section_ID.x.x, d$yearweek_id.x, sep = "")
d <- d %>%
  dplyr::select(year, yearweek_id.x, imputed.sock2019, imputed.chin2019, imputed.pink2019, riversection_weekid)

# Calculate biomass
d$sock.biomass <- d$imputed.sock2019 *2.7
d$spring.biomass <- d$imputed.chin2019 *13.6
d$pink.biomass <- d$imputed.pink2019 * 2.5
d$total.biomass = rowSums(d[,c("sock.biomass","spring.biomass","pink.biomass")])

# Join 2019 to the other 2 years
salmon2020.2021 <- final %>%
  dplyr::select(year.x, yearweek_id.x, riv.section.weekid, total.biomass, River_section_ID.x.x)

salmon2020.2021 <- salmon2020.2021 %>%
  rename(year = year.x,
         riversection_weekid = riv.section.weekid)

salmon2020.2021$yearweek_id.x <- as.character(salmon2020.2021$yearweek_id.x)

salmon2019 <- d %>%
  dplyr::select(year, yearweek_id.x, riversection_weekid, total.biomass)

salmon <- bind_rows(salmon2019, salmon2020.2021)

n_last <- 2
salmon$week <- substr(salmon$yearweek_id.x, nchar(salmon$yearweek_id.x) - n_last + 1, nchar(salmon$yearweek_id.x))

# Figure S6
salmon_impute <- ggplot(salmon, aes(x = week, y = total.biomass, group = River_section_ID.x.x)) +
  geom_line(aes(colour = River_section_ID.x.x), size = 1) + 
  scale_colour_brewer(name = "River Segment", 
                      palette = "Set1",
                      labels = c("Lower", "Mid-Lower", "Mid-Upper", "Upper")) +
  facet_grid(rows = vars(year), scales = "free") +
  labs(x = "Week", y = "Total salmon biomass (kg)") +
  geom_vline(data=filter(salmon, year=="2019"), aes(xintercept=2)) + 
  geom_vline(data=filter(salmon, year=="2019"), aes(xintercept=11), linetype = "dotted") + 
  geom_vline(data=filter(salmon, year=="2020"), aes(xintercept=2)) + 
  geom_vline(data=filter(salmon, year=="2020"), aes(xintercept=11), linetype = "dotted") + 
  geom_vline(data=filter(salmon, year=="2021"), aes(xintercept=1)) + 
  geom_vline(data=filter(salmon, year=="2021"), aes(xintercept=10), linetype = "dotted") + 
  theme_classic()

salmon_impute
ggsave(file="salmon_impute.png", width=6, height=4, dpi=300)

salmon_impute_bars <- 
  ggplot(salmon, aes(x = week, y = total.biomass, group = River_section_ID.x.x)) +
  geom_col(aes(fill = River_section_ID.x.x), size = 1) + 
  scale_fill_brewer(name = "River Segment", 
                    palette = "Set1",
                    labels = c("Lower", "Mid-Lower", "Mid-Upper", "Upper")) +
  facet_grid(rows = vars(year)) +
  labs(x = "Week", y = "Salmon biomass (kg)") +
  theme_classic()

# Call in berry data 
# The number of species ripe is the sum of the number of unique ripe species within a revisit from berry_surveys.csv.
berries <- read.csv("berries_summarized.csv")

# Call in human data 
humans <- read.csv("humans.csv")
humans_hist <- na.omit(humans) #one empty row
humans_hist <- humans_hist %>%
  mutate(Treatment =
           case_when(site_name == "Stillwater Trail" ~ "No tour",
                     site_name == "Mosher Trail" ~ "Land-based",
                     site_name == "Fisheries Pool" ~ "Land- and boat-based"))
# SI Figure
hum.hist <- humans_hist %>%
  ggplot(aes(x = total_estimated_humans,fill = Treatment)) +
  geom_histogram(alpha=0.6, position = 'identity') +
  facet_wrap(Treatment ~., scales = "free_x") +
  scale_fill_brewer(palette = "Dark2") +
  labs(fill="Spatial treatment") +
  ylab(label = "Frequency") +
  xlab(label = "Visitors") +
  theme_classic() +
  theme(legend.position = "none") 
hum.hist
ggsave(file="hum.hist.png", width=6, height=4, dpi=300)

# Humans over time
hum2 <- humans
hum2 <- na.omit(hum2) #one empty row
hum2 <- hum2%>%
  mutate(Treatment =
           case_when(site_name == "Stillwater Trail" ~ "No tour",
                     site_name == "Mosher Trail" ~ "Land-based",
                     site_name == "Fisheries Pool" ~ "Land- and boat-based"))
hum2$year <- as.factor(hum2$year)
hum2$week <- as.factor(hum2$week)
hum2$Treatment <- as.factor(hum2$Treatment)

# SI Figure
ggplot(hum2, aes(x = week, y = total_estimated_humans, fill = Treatment)) +
  geom_bar(alpha=0.6,stat = "identity") +
  facet_grid(rows=vars(year), scales = "free_y") +
  labs(x = "Week",
       y = "Visitors") +
  scale_fill_brewer(palette = "Dark2") +
  theme_classic() +
  theme(legend.title = element_blank())
ggsave(file="visitors.over.time.png", width=6, height=4, dpi=300)

# Call in hydrometric data
# Real-Time Hydrometric Data Graph for ATNARKO RIVER NEAR THE MOUTH (08FB006) [BC]
# Source: Government of Canada
# https://wateroffice.ec.gc.ca/report/real_time_e.html?stn=08FB006 
hydro <- read.csv("hydro.csv")
hydro$date <- hydro$Timestamp..UTC.08.00.
hydro <- hydro %>%
  group_by(year = year(date), week = isoweek(date)) %>% 
  mutate(yearweek_id = paste(year,week,sep = ""))
hydro <- hydro %>%
  group_by(year, week, yearweek_id) %>%
  summarise(weekly_mean_waterlevel.m = mean(Value))
hydro$yearweek_id <- as.character(hydro$yearweek_id)
number_ticks <- function(n) {function(limits) pretty(limits, n)}

# SI Figure
ggplot(hydro) +
  geom_point(aes(x = week, y = weekly_mean_waterlevel.m)) +
  facet_grid(row = vars(year)) +
  labs(y = "Mean water level (m)", x = "Week") +
  scale_x_continuous(breaks=number_ticks(10)) +
  geom_vline(data=filter(salmon, year=="2019"), aes(xintercept=33)) + 
  geom_vline(data=filter(salmon, year=="2019"), aes(xintercept=42), linetype = "dotted") + 
  geom_vline(data=filter(salmon, year=="2020"), aes(xintercept=33)) + 
  geom_vline(data=filter(salmon, year=="2020"), aes(xintercept=42), linetype = "dotted") + 
  geom_vline(data=filter(salmon, year=="2021"), aes(xintercept=32)) + 
  geom_vline(data=filter(salmon, year=="2021"), aes(xintercept=41), linetype = "dotted") + 
  theme_classic()
ggsave(file="hydro.png", width=6, height=4, dpi=300)
