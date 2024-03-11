
#############################################
# Script 1
# Data preparation and exploration; includes script to subset study period of interest, 
# identify camera activity gaps, and calculate diel period metric. 
#############################################

# setwd() 

# This script produces the following two files:
# 1) cam.sum (sum denotes summarise): This is the dataset to be used for the detection rate analysis
# 2) cam.multi (multi denotes multinomial): This is the dataset to be used for the multinomial analysis 

# Load Packages
list.of.packages <- c("leaflet", "colortools", "kriging", "corrplot", "lubridate", "ggplot2", "tidyr", "hms", "dplyr", "circular", "suncalc", "ISOweek", "igraph", "tidyverse", "data.table")
# Check you have them and load them
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

cam <- read.csv("atnarko_detection_data.csv")
sta <- read.csv("atnarko_station_data.csv")
eff <- read.csv("atnarko_deployment_data.csv")

cam <- cam %>%
  mutate(date_of_video = ymd(date_of_video)) %>%
  filter(between(yday(date_of_video), 227, 288))

# For each site and each year, subset a date range with the following conditions:
# The range will begin from the latest date prior to August 15th, and will end at the earliest date after October 15th. These are our effort data.

eff <- eff %>%
  group_by(Deployment.Location.ID, year) %>%
  filter(between(Camera.Deployment.Begin.Date,
                 max(Camera.Deployment.Begin.Date[
                   Camera.Deployment.Begin.Date < as.Date(ISOdate(year, 8, 15))]),
                 min(Camera.Deployment.Begin.Date[
                   Camera.Deployment.Begin.Date > as.Date(ISOdate(year, 10, 15))]))) %>%
  ungroup

# Data checks: code revised from Christopher Beirne, Cole Burton's Wildlife Coexistence Lab at UBC, and the WildCAM Network.
# Timezone [Use UTC if your cameras do not correct for daylight saving time, if they do use the timezone where the data were collected]
tz <- "UTC"

# 1) dat$Blank must be logical
is.logical(cam$Blank)
# If this is FALSE convert this column to TRUE/FALSE
# cam$Blank <- as.logical(cam$Blank)

# 2) All dates must be in YYYY-MM-DD in 'eff' and YYYY-MM-DD HH:MM:SS in 'cam'
# If the following return NA, change formatting
as.Date(ymd_hms(eff$Camera.Deployment.Begin.Date[1], truncated=3))

# Format date
eff$Camera.Deployment.Begin.Date <- format(as.Date(eff$Camera.Deployment.Begin.Date, format = "%Y-%m-%d"), "%Y-%m-%d")

# If the following return NA, change formatting
ymd_hms(cam$Date_Time.Captured[1],truncated=2)

# 3) the dates in 'eff$Camera.Deployment.End.Date' must be the when the camera fails, not when you check the camera. If the camera fails (due to damage or full sd card), use the last day it functions.

# 4) Ensure your species names are consistent - check in the list below
as.data.frame(table(cam$Species))

# 5) Ensure Number.of.Animals and doesn't have any non-numeric data. The following should return TRUE
is.numeric(cam$Number.of.Bears)

# 6) Ensure all deployment dates are before retrieval dates for each deployment
# Logic = are the stations active for 0 or more days -> all should read TRUE
table((strptime(eff$Camera.Deployment.End.Date, "%Y-%m-%d", tz="UTC")-strptime(eff$Camera.Deployment.Begin.Date, "%Y-%m-%d", tz="UTC"))>=0)

# 7) Do you have lat/long data for all of your sites you have effort data for? If yes, the value should be 0
length(setdiff(eff$Deployment.Location.ID, sta$Deployment.Location.ID))

# Prepare dates
ymd(eff$Camera.Deployment.Begin.Date[1])
ymd_hms(cam$Date_Time.Captured[1],truncated=2)
eff$Camera.Deployment.Begin.Date <- strptime(as.Date(ymd_hms(eff$Camera.Deployment.Begin.Date, truncated=3, tz=tz)), "%Y-%m-%d", tz=tz)
eff$Camera.Deployment.End.Date   <- strptime(as.Date(ymd_hms(eff$Camera.Deployment.End.Date, truncated=3, tz=tz)), "%Y-%m-%d", tz=tz)
eff$Days <- as.numeric(round(difftime(eff$Camera.Deployment.End.Date, eff$Camera.Deployment.Begin.Date, units="days"),1))
cam$Date_Time.Captured <- ymd_hms(cam$Date_Time.Captured, truncated=2, tz=tz)

# Report camera inactivity
# Assign the number of days within a week the cameras were active.
# First apply independence rule

# Subset grizzly bears only
cam %>%
  group_by(species.1) %>%
  tally

cam <- subset(cam, species.1 != "B")
cam <- subset(cam, species.1 != "U")

cam <- cam %>%
  arrange(Deployment.Location.ID, Date_Time.Captured) %>%
  group_by(Deployment.Location.ID) %>%
  mutate(diff = c(Inf, diff(as.numeric(Date_Time.Captured)/60))) %>%
  ungroup

nr <- nrow(cam)
g <- make_empty_graph(n = nr)
wx <- which(cam$diff < 30)
g <- add_edges(g, c(rbind(wx - 1, wx)))
cam$membership <- components(g)$membership

cam <- cam %>%
  group_by(membership) %>%
  slice(1) %>%
  ungroup

# Report camera inactivity 
eff$Camera.Deployment.Begin.Date <- as.Date(as.POSIXct((eff$Camera.Deployment.Begin.Date), format = "%Y-%m-%d"))
eff$Camera.Deployment.End.Date <- as.Date(as.POSIXct((eff$Camera.Deployment.End.Date), format = "%Y-%m-%d"))

camact <- eff %>%  
  dplyr::arrange(Camera.Deployment.Begin.Date) %>%
  # unnest day sequence from start to end into df 
  dplyr::group_by(rn = dplyr::row_number()) %>% 
  dplyr::mutate(dates = list(seq.Date(from = Camera.Deployment.Begin.Date, to = Camera.Deployment.End.Date, by = "days"))) %>% 
  unnest(dates) %>%
  dplyr::ungroup() %>%
  # right join list of all dates with iso week and year
  dplyr::right_join(dplyr::tibble(dates = seq.Date(from = min(eff$Camera.Deployment.Begin.Date), max(eff$Camera.Deployment.End.Date), by = "days")) %>%
                      dplyr::mutate(year.2 = lubridate::year(dates),
                                    iso_week = lubridate::isoweek(dates)),
                    by = "dates") %>%
  # fill up the rn in case it is zero with a number that is larger all rns
  dplyr::mutate(rn = ifelse(is.na(rn), nrow(eff) + 1, rn)) %>%
  # summarize data
  dplyr::group_by(year.2, Deployment.Location.ID, iso_week, rn) %>%
  dplyr::summarize(bdate =  min(Camera.Deployment.Begin.Date, na.rm = TRUE), 
                   edate = min(Camera.Deployment.End.Date, na.rm = TRUE), 
                   days = sum(ifelse(is.na(Camera.Deployment.Begin.Date), 0, 1))) %>%
  dplyr::ungroup() %>%
  # get lowest sequential numbering per week 
  dplyr::group_by(year.2,Deployment.Location.ID,iso_week) %>%
  dplyr::slice_min(order_by = rn, n = 1) %>%
  dplyr::ungroup()

# Omit NAs (which are rows in which the study was not active)
camact <- camact[!is.na(camact$Deployment.Location.ID),]

# Convert iso_week column to include the year
camact <- camact %>%
  mutate(isoweek = paste(year.2, iso_week, sep = ""))

# Create a list of active days per site
active <- camact %>%
  
  # get days from begin.date through end.date per row
  # (these will be stored as numerical instead of dates)
  rowwise() %>%
  mutate(days.in.range = list(ymd(bdate):ymd(edate))) %>%
  
  # concatenate per group
  group_by(Deployment.Location.ID) %>%
  summarise(days.in.range = list(unlist(days.in.range)))

# change from data frame to list
active.list <- active$days.in.range
names(active.list) <- active$Deployment.Location.ID

# for each row in camact, get ISO week days and tally how many match with active.list

result <- camact %>%
  
  # create (temporary) columns for properly formatted ISO week, and the iso week first and last day
  mutate(.isoweek = gsub('(\\d{4})(\\d+)', '\\1-W\\2', isoweek),
         .isoweekfirst = ISOweek2date(paste0(.isoweek,'-',1)),
         .isoweeklast = ISOweek2date(paste0(.isoweek,'-',7))) %>%
  
  # create lists of days in ISO week (`.isodays`)
  rowwise() %>%
  mutate(.isodays = list(.isoweekfirst:.isoweeklast) ) %>%
  ungroup() %>%
  
  # for each row, check how many values in .isodays match with those for the respective group in `active.list`
  mutate(active.days.in.iso.week = unlist(map2(.isodays,
                                               Deployment.Location.ID,
                                               function(.isodays, Deployment.Location.ID) sum(.isodays %in% active.list[[Deployment.Location.ID]])))) 

# Add a column that is just the site number
n_last <- 2
result$site_number <- substr(result$Deployment.Location.ID, nchar(result$Deployment.Location.ID) - n_last + 1, nchar(result$Deployment.Location.ID))

# Create a histogram that looks at the frequency of weeks by the number of days within a week that a camera was active
result <- result %>%
  filter(iso_week == "33" | 
           iso_week == "34" |
           iso_week == "35" |
           iso_week =="36" |
           iso_week == "37"|
           iso_week == "38" |
           iso_week == "39" |
           iso_week == "40" |
           iso_week == "41" |
           iso_week == "42")

result <- result %>%
  mutate(Treatment =
           case_when(Deployment.Location.ID == "2019AT01" ~ "Reference",
                     Deployment.Location.ID == "2019AT03" ~ "Reference",
                     Deployment.Location.ID == "2019AT05" ~ "Reference",
                     Deployment.Location.ID == "2019AT06" ~ "Reference",
                     Deployment.Location.ID == "2019AT07" ~ "Reference",
                     Deployment.Location.ID == "2019AT08" ~ "Reference",
                     Deployment.Location.ID == "2019AT10" ~ "Landbased",
                     Deployment.Location.ID == "2019AT11" ~ "Landbased",
                     Deployment.Location.ID == "2019AT12" ~ "Landbased",
                     Deployment.Location.ID == "2019AT13" ~ "Landbased",
                     Deployment.Location.ID == "2019AT15" ~ "Landbased",
                     Deployment.Location.ID == "2019AT17" ~ "Landbased",
                     Deployment.Location.ID == "2019AT18" ~ "Landbased",
                     Deployment.Location.ID == "2019AT20" ~ "Landbased",
                     Deployment.Location.ID == "2019AT21" ~ "Landbased",
                     Deployment.Location.ID == "2019AT21B" ~ "Landbased",
                     Deployment.Location.ID == "2019AT23" ~ "Landbased",
                     Deployment.Location.ID == "2019AT24" ~ "Platform/Drift",
                     Deployment.Location.ID == "2019AT24B" ~ "Platform/Drift",
                     Deployment.Location.ID == "2019AT25" ~ "Platform/Drift",
                     Deployment.Location.ID == "2019AT26" ~ "Platform/Drift",
                     Deployment.Location.ID == "2019AT27" ~ "Platform/Drift",
                     Deployment.Location.ID == "2019AT28" ~ "Landbased",
                     Deployment.Location.ID == "2019AT30" ~ "Platform/Drift",
                     Deployment.Location.ID == "2019AT31" ~ "Reference",
                     Deployment.Location.ID == "2019AT32" ~ "Reference",
                     Deployment.Location.ID == "2019AT33" ~ "Reference",
                     Deployment.Location.ID == "2019AT36" ~ "Platform/Drift",
                     Deployment.Location.ID == "2019AT37" ~ "Platform/Drift",
                     Deployment.Location.ID == "2020AT01" ~ "Reference",
                     Deployment.Location.ID == "2020AT03" ~ "Reference",
                     Deployment.Location.ID == "2020AT05" ~ "Reference",
                     Deployment.Location.ID == "2020AT06" ~ "Reference",
                     Deployment.Location.ID == "2020AT07" ~ "Reference",
                     Deployment.Location.ID == "2020AT08" ~ "Reference",
                     Deployment.Location.ID == "2020AT10" ~ "Landbased",
                     Deployment.Location.ID == "2020AT11" ~ "Landbased",
                     Deployment.Location.ID == "2020AT12" ~ "Landbased",
                     Deployment.Location.ID == "2020AT13" ~ "Landbased",
                     Deployment.Location.ID == "2020AT15" ~ "Landbased",
                     Deployment.Location.ID == "2020AT17" ~ "Landbased",
                     Deployment.Location.ID == "2020AT18" ~ "Landbased",
                     Deployment.Location.ID == "2020AT20" ~ "Landbased",
                     Deployment.Location.ID == "2020AT21" ~ "Landbased",
                     Deployment.Location.ID == "2020AT21B" ~ "Landbased",
                     Deployment.Location.ID == "2020AT23" ~ "Landbased",
                     Deployment.Location.ID == "2020AT24" ~ "Platform/Drift",
                     Deployment.Location.ID == "2020AT24B" ~ "Platform/Drift",
                     Deployment.Location.ID == "2020AT25" ~ "Platform/Drift",
                     Deployment.Location.ID == "2020AT26" ~ "Platform/Drift",
                     Deployment.Location.ID == "2020AT27" ~ "Platform/Drift",
                     Deployment.Location.ID == "2020AT28" ~ "Landbased",
                     Deployment.Location.ID == "2020AT30" ~ "Platform/Drift",
                     Deployment.Location.ID == "2020AT31" ~ "Reference",
                     Deployment.Location.ID == "2020AT32" ~ "Reference",
                     Deployment.Location.ID == "2020AT33" ~ "Reference",
                     Deployment.Location.ID == "2020AT36" ~ "Platform/Drift",
                     Deployment.Location.ID == "2020AT37" ~ "Platform/Drift",
                     Deployment.Location.ID == "2021AT01" ~ "Reference",
                     Deployment.Location.ID == "2021AT03" ~ "Reference",
                     Deployment.Location.ID == "2021AT05" ~ "Reference",
                     Deployment.Location.ID == "2021AT06" ~ "Reference",
                     Deployment.Location.ID == "2021AT07" ~ "Reference",
                     Deployment.Location.ID == "2021AT08" ~ "Reference",
                     Deployment.Location.ID == "2021AT10" ~ "Landbased",
                     Deployment.Location.ID == "2021AT11" ~ "Landbased",
                     Deployment.Location.ID == "2021AT12" ~ "Landbased",
                     Deployment.Location.ID == "2021AT13" ~ "Landbased",
                     Deployment.Location.ID == "2021AT15" ~ "Landbased",
                     Deployment.Location.ID == "2021AT17" ~ "Landbased",
                     Deployment.Location.ID == "2021AT18" ~ "Landbased",
                     Deployment.Location.ID == "2021AT20" ~ "Landbased",
                     Deployment.Location.ID == "2021AT21" ~ "Landbased",
                     Deployment.Location.ID == "2021AT21B" ~ "Landbased",
                     Deployment.Location.ID == "2021AT23" ~ "Landbased",
                     Deployment.Location.ID == "2021AT24" ~ "Platform/Drift",
                     Deployment.Location.ID == "2021AT24B" ~ "Platform/Drift",
                     Deployment.Location.ID == "2021AT25" ~ "Platform/Drift",
                     Deployment.Location.ID == "2021AT26" ~ "Platform/Drift",
                     Deployment.Location.ID == "2021AT27" ~ "Platform/Drift",
                     Deployment.Location.ID == "2021AT28" ~ "Landbased",
                     Deployment.Location.ID == "2021AT30" ~ "Platform/Drift",
                     Deployment.Location.ID == "2021AT31" ~ "Reference",
                     Deployment.Location.ID == "2021AT32" ~ "Reference",
                     Deployment.Location.ID == "2021AT33" ~ "Reference",
                     Deployment.Location.ID == "2021AT36" ~ "Platform/Drift",
                     Deployment.Location.ID == "2021AT37" ~ "Platform/Drift"))

# Frequency of site-weeks that had 1-7 active camera days

# SI Figure
ggplot(result, aes(x = active.days.in.iso.week)) +
  geom_histogram(aes(y=(..count..)/sum(..count..)),
                 position = "identity",binwidth=1, fill="grey", color="black", ) +
 geom_text(aes(y = ((..count../sum(..count..))),
               label = scales::percent((round(..count../sum(..count..)*100))/100, scale_factor = 1)),
           stat = "count", vjust = -.3)+
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Days camera active in a site-week", y = "Percent") +
  scale_x_continuous(breaks=1:7) +
  theme_classic()
ggsave(file="Day.cam.active.all.zone.png", width=6, height=4, dpi=300)

# Subset each treatment 
reference <- result %>%
  filter(Treatment == "Reference")
land <- result %>%
  filter(Treatment == "Landbased")
plat <- result %>%
  filter(Treatment == "Platform/Drift")

# SI Figure
ggplot(reference, aes(x = active.days.in.iso.week)) +
  geom_histogram(aes(y=(..count..)/sum(..count..)),
                 position = "identity",binwidth=1, fill="grey", color="black", ) +
  geom_text(aes(y = ((..count../sum(..count..))),
                label = scales::percent((round(..count../sum(..count..)*100))/100, scale_factor = 1)),
            stat = "count", vjust = -0.5)+ 
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Days camera active in a site-week", y = "Percent") +
  scale_x_continuous(breaks=1:7) +
  theme_classic()
ggsave(file="Days.cam.active.no.tour.zone.png", width=6, height=4, dpi=300)

# SI Figure
ggplot(land, aes(x = active.days.in.iso.week)) +
  geom_histogram(aes(y=(..count..)/sum(..count..)),
                 position = "identity",binwidth=1, fill="grey", color="black", ) +
  geom_text(aes(y = ((..count../sum(..count..))),
                label = scales::percent((round(..count../sum(..count..)*100))/100, scale_factor = 1)),
            stat = "count", vjust = -0.5)+ 
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Days camera active in a site-week", y = "Percent") +
  scale_x_continuous(breaks=1:7) +
  theme_classic()
ggsave(file="Days.cam.active.landbased.zone.png", width=6, height=4, dpi=300)

# SI Figure
ggplot(plat, aes(x = active.days.in.iso.week)) +
  geom_histogram(aes(y=(..count..)/sum(..count..)),
                 position = "identity",binwidth=1, fill="grey", color="black", ) +
  geom_text(aes(y = ((..count../sum(..count..))),
                label = scales::percent((round(..count../sum(..count..)*100))/100, scale_factor = 1)),
            stat = "count", vjust = -0.5)+ 
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Days camera active in a site-week", y = "Percent") +
  scale_x_continuous(breaks=1:7) +
  theme_classic()
ggsave(file="Days.cam.active.land.boat.zone.png", width=6, height=4, dpi=300)

# Summarise bear detections per site-week
# This is the dataset we use for the detection rate analysis 
cam.sum <- cam %>%
  group_by(year, week, site_name, Treatment, riv.section.weekid, siteweek_id) %>%
  summarise(detections.per.siteweek = sum(bear_present))
cam.sum$Treatment <- factor(cam.sum$Treatment, levels = c("Platform/Drift", "Landbased", "Reference"))

ggplot(cam.sum) +
  geom_point(aes(x = week, y = detections.per.siteweek)) +
  facet_grid(row = vars(year), cols = vars(Treatment), scales = "free_x") +
  labs(y = "Detections per site-week", x = "Week") +
  theme_classic()  

# Prepare a separate dataset for the multinomial analysis.
cam.multi <- cam

# Check correct date format
cam.multi$Date_Time.Captured <- as.POSIXct(cam.multi$Date_Time.Captured, format = "%Y-%m-%d %H:%M:%S")

cam.multi$time <- as_hms(cam.multi$Date_Time.Captured)
class(cam.multi$time)
cam.multi$date_of_video <- as.Date(cam.multi$date_of_video, "%Y-%m-%d") 
class(cam.multi$date_of_video)

# Define days 7am-7am
cam.multi <- cam.multi %>% 
  mutate(time_of_day = case_when(
    time >= lubridate::hms("00:00:00") & time < lubridate::hms("07:00:00") ~ 1,
    time >= lubridate::hms("07:00:00")  & time <= lubridate::hms("23:59:59") ~ 2,
  )
  )

cam.multi <- cam.multi %>%
  mutate(Date_7am = case_when(
    time_of_day == 2 ~ date_of_video,
    time_of_day == 1 ~ date_of_video - 1,
  ))

# Get sunset and sunrise times
sun <- getSunlightTimes(date = cam.multi$date_of_video, lat = 52, lon = -126, tz = "US/Pacific")
sun <- sun[c(1, 8, 9)]
sun$daylightHours <- sun$sunsetStart - sun$sunriseEnd
names(sun)[names(sun) == "date"] <- "date_of_video"
sun <- sun[!duplicated(sun$date_of_video),]
cam.multi <- left_join(cam.multi, sun, by = "date_of_video")

# Assign the dark point. Dark point is defined by the mid-point between sunset and sunrise. 
cam.multi$nightlightHours <- 24 - cam.multi$daylightHours
cam.multi$DarkPoint <- (cam.multi$nightlightHours/2) + cam.multi$sunsetStart
cam.multi$DarkPoint <- as_hms(cam.multi$DarkPoint)

# Get mean dark point across the study period 
seconds_to_period(mean(period_to_seconds(lubridate::hms(cam.multi$DarkPoint))))
# [1] "1H 20M 50
cam.multi$DarkPoint <- as_hms("01:20:50")

# Calculate diurnality metric. Subtract the time of detection from the dark point. The lower the output, the more nocturnal. This output is the number of hours from the darkest point in the night. 
# Code for how close the detection is to the dark point time, no matter the date. 
cam.multi <- cam.multi %>%
  mutate(diff   = lubridate::hms(cam.multi$time)-lubridate::hms(cam.multi$DarkPoint),
         TimeDiff = abs(as.numeric(if_else(diff > hours(12), diff-hours(24), diff))/3600)
  )

#Re-categorize age sex categories
cam.multi <- cam.multi %>%
  mutate(AgeSex = case_when(
    age1 == "SUBADULT" ~ "Sub-adult",
    sex1 == "FWY" ~ "Female with young", 
    age1 == "ADULT" & sex1 == "M" ~ "Adult male",
    age1 == "ADULT" & sex1 == "F" ~ "Adult female",
    age1 == "ADULT" & sex1 == "U" ~ "Unknown adult",
    sex1 == "U" & sex2 == "U" ~ "2 Adults",
    sex1 == "F" & sex2 == "U" ~ "2 Adults",
    sex1 == "M" & sex2 == "U" ~ "2 Adults",
    sex1 == "F" & sex2 == "M" ~ "2 Adults",
    sex1 == "M" & sex2 == "F" ~ "2 Adults",
  ))

# Prop male vs. weekly detections 
plot <- cam.multi %>%
  mutate(male = ifelse(AgeSex == "Adult male", 1, 0))

plot <- plot %>%
  mutate(male = ifelse(is.na(male), 0, male))

plot <- plot %>%
  group_by(yearweek_id) %>%
  summarize(
    weekly_detection_rate = sum(!is.na(AgeSex)),
    prop.male = sum(male) / sum(!is.na(AgeSex)))

# SI Figure 
ggplot(plot, aes(x = weekly_detection_rate, y = prop.male)) +
  geom_point() +
  labs(
    x = "Weekly detections",
    y = "Proportion male") + 
  geom_smooth(method = "lm", se = TRUE, color = "blue")+
  theme_classic()
ggsave(file="prop.male.by.det.rates.png", width=6, height=4, dpi=300)