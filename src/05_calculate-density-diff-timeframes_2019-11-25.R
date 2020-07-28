#-------------------------------------------------------------------------------

# Title: White-tail deer density for Mel et Maud
# Created: November 25, 2019
# Author: Marcus Becker

# Objectives: Calculate density of white-tailed deer for each 2018 ABMI and CMU
#             camera deployment at five specific time periods.

#-------------------------------------------------------------------------------

# Load packages

library(tidyverse)

#-------------------------------------------------------------------------------

# Import data

# Set path to Mel's public drive folder
public_tomel <- "Z:/toMelanie/WTD-Density/"

# Lookup vector for CMU deployment abbreviations
cmu_abb <- read_csv(paste0(public_tomel, "cmu-dep-abb.csv")) %>%
  pull(abb)

# Lookup table for distance groups
df_dist_groups <- read_csv(paste0(public_tomel, "species-distance-groups.csv"))

# Camera operating days (All)
time.by.day <- read_csv(paste0(public_tomel, "camera-operating-days-all.csv"))

# Series data
df_series_all <- read_csv(paste0(public_tomel, "series-data-all.csv"))

# Detection distance modeling
# load(paste0(public_tomel,"Detection distances by site species and season.rData"))
load("G:/Shared drives/ABMI Camera Mammals/data/processed/detection-distance/predictions/Detection distances by site species and season_2020-05-27.rData")

#-------------------------------------------------------------------------------

# Set parameters

cam_fov_ang <- 42

#-------------------------------------------------------------------------------

# Generate time-by-day summaries for each requested time period.
df_tbd_summary <- time.by.day %>%
  as_tibble(rownames = "DeploymentYear") %>%
  # Numbers of columns correspond to Julian day
  mutate(total = rowSums(select(., -DeploymentYear)),
         # Spring is April 1 to June 17
         total.spring = rowSums(select(., V91:V168)),
         # Winter months
         total.nov = rowSums(select(., V305:V334)),
         total.dec = rowSums(select(., V335:V365)),
         total.jan = rowSums(select(., V1:V31)),
         total.feb = rowSums(select(., V32:V59)),
         total.mar = rowSums(select(., V60:V90)))%>%
  select(DeploymentYear, total:total.mar) %>%
  # Filter only for ABMI and CMU deployments from 2018
  filter(str_detect(DeploymentYear, "2018|2017|2016|2015|2019")) %>%
  filter(str_detect(DeploymentYear, "^ABMI"),
         !str_detect(DeploymentYear, "OG"))

# Calculate total duration spent on camera by deployment, year, and species:
df_series_summary <- df_tt %>%
  rename(DeploymentYear = name_year) %>%
  # Correct number of individuals issue
  # mutate(series_total_time = mean_animals * series_total_time) %>%
  # Filter only for ABMI and CMU deployments from 2018
  filter(str_detect(DeploymentYear, "2018|2017|2016|2015|2019")) %>%
  filter(str_detect(DeploymentYear, "^ABMI"),
         !str_detect(DeploymentYear, "OG")) %>%
  # Calculate Julian day and time period of each series
  mutate(julian = as.numeric(format(date_detected, "%j")),
         season = case_when(
           julian >= 1 & julian <= 31 ~ "January",
           julian >= 32 & julian <= 59 ~ "February",
           julian >= 60 & julian <= 90 ~ "March",
           julian >= 305 & julian <= 334 ~ "November",
           julian >= 335 & julian <= 365 ~ "December",
           julian >= 91 & julian <= 168 ~ "Spring",
           TRUE ~ "Other")) %>%
  mutate_at(c("DeploymentYear", "common_name", "season"), factor) %>%
  # Important to add .drop = F argument in group_by so all combinations are preserved.
  group_by(DeploymentYear, common_name, season, .drop = FALSE) %>%
  # Calculate total duration in front of camera
  summarise(total_duration = sum(series_total_time)) %>%
  ungroup() %>%
  mutate_if(is.factor, as.character) %>%
  # Join operating days
  left_join(df_tbd_summary, by = "DeploymentYear") %>%
  filter(!season == "Other")
  # ^Note: Missing info from deployments w/o any images of native mammals. Added below.

# Character vector of all native species (41 in total)
chr_species <- as.character(sort(unique(df_series_summary$common_name)))

# Add deployments w/o any images of native mammals
df_series_nonative <- df_tbd_summary %>%
  # Preserve only those Deployment-Years w/o images of native mammals.
  anti_join(df_series_summary, by = "DeploymentYear") %>% # 65 Deployment-Years
  # Expand to include all combinations of DeploymentYear, season, and common_name.
  expand(DeploymentYear,
         season = c("January", "February", "March", "November", "December", "Spring"),
         common_name = chr_species) %>%
  # Join in operating days.
  left_join(df_tbd_summary, by = "DeploymentYear") %>%
  # Total_Duration is 0 because these species weren't seen.
  mutate(total_duration = 0)

# Bind together the two df's of series (w/ and w/o native mammals)
df_series_full <- df_series_summary %>%
  bind_rows(df_series_nonative) %>%
  arrange(DeploymentYear, common_name, season) %>%
  # Combine total.x's into `days` (time period specified as `season`)
  mutate(days = case_when(
    season == "January" ~ total.jan,
    season == "February" ~ total.feb,
    season == "March" ~ total.mar,
    season == "November" ~ total.nov,
    season == "December" ~ total.dec,
    season == "Spring" ~ total.spring)) %>%
  select(-c(total:total.mar))

# Add detection distance information
# Tidy first:
df_detdist <- dd %>%
  as.data.frame() %>%
  rownames_to_column(var = "DeploymentYear") %>%
  gather(key = "SpGroupSeason", value = "detdist", Bear.Summer:`Bighorn sheep.Winter`) %>%
  mutate(SpGroupSeason = str_replace_all(SpGroupSeason, "[//(//)// ]", "")) %>%
  # Create two new columns: Detection Distance Group and Season, sep by "."
  separate(SpGroupSeason, into = c("detdistgroup", "season")) %>%
  mutate(detdistgroup = str_replace(detdistgroup, "wapiti", "")) %>%
  # Select only 2018 ABMI & CMU deployments
  filter(str_detect(DeploymentYear, "2018|2017|2016|2015|2019")) %>%
  filter(str_detect(DeploymentYear, "^ABMI"),
         !str_detect(DeploymentYear, "OG")) %>%
  # Create new time periods - kind of a hacky way to do this.
  crossing(new.season = c("Spring", "January", "February",  "March", "November", "December")) %>%
  # Set criteria. Note that I'm using the Winter models for Jan,Feb,Mar,Nov,Dec, Summer for Spring.
  filter((season == "Winter" & new.season == "January") |
         (season == "Winter" & new.season == "February") |
         (season == "Winter" & new.season == "March") |
         (season == "Winter" & new.season == "November") |
         (season == "Winter" & new.season == "December") |
         (season == "Summer" & new.season == "Spring")) %>%
  select(DeploymentYear, detdistgroup, season = new.season, detdist) %>%
  mutate(DeploymentYear = toupper(DeploymentYear))

# Combine, and then we'll have all the ingredients we need for density calc.
df_dens_ing <- df_series_full %>%
  mutate(DeploymentYear = toupper(DeploymentYear)) %>%
  left_join(df_dist_groups, by = "common_name") %>%
  rename(detdistgroup = dist_group) %>%
  mutate(detdistgroup = ifelse(detdistgroup == "Bighorn sheep", "Bighornsheep", detdistgroup),
         detdistgroup = ifelse(detdistgroup == "Mule deer", "Muledeer", detdistgroup)) %>%
  left_join(df_detdist, by = c("DeploymentYear", "season", "detdistgroup")) %>%
  select(DeploymentYear, common_name, detdistgroup, season, detdist, everything())

# Calculate density
df_density <- df_dens_ing %>%
  # Calculate effort per site, in 100-m^2 * days
  mutate(effort = days * (detdist^2 * pi * (cam_fov_ang/360)) / 100,
         # Calculate seconds (s) of animal presence per effort
         cpue = total_duration / effort,
         # Convert to per km^2
         cpue_km2 = cpue / 60 / 60 / 24 * 10000)

# Finally, filter for species of interest
df_density_wtd <- df_density %>%
  filter(common_name == "White-tailed Deer") %>%
  mutate(season = factor(season,
                         levels = c("November", "December", "January",
                                    "February", "March", "Spring"))) %>%
  select(DeploymentYear, common_name, season, days, density_km2 = cpue_km2) %>%
  arrange(DeploymentYear, season) %>%
  na_if("NaN") %>%
  na_if("Inf")

# Lure adjustment
df_lure <- read_csv("G:/Shared drives/ABMI Camera Mammals/data/lookup/lure/all-cam-lure_2020-06-01.csv") %>%
  mutate(DeploymentYear = paste0(name, "_", year)) %>%
  select(DeploymentYear, lure) %>%
  distinct()

all <- read_csv("G:/Shared drives/ABMI Camera Mammals/results/density/Marcus/abmi-all-years_density_2020-06-09.csv") %>%
  filter(str_detect(name_year, "ABMI")) %>%
  filter(common_name == "White-tailed Deer") %>%
  group_by(name_year) %>%
  summarise(density = mean(density_km2, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(df_lure, by = "name_year") %>%
  na_if("Inf") %>%
  na_if("NaN") %>%
  group_by(lure) %>%
  summarise(density = mean(density, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(!is.na(lure))

lure.effect <- 0.96

df_density_wtd_lure <- df_density_wtd %>%
  left_join(df_lure, by = "DeploymentYear") %>%
  mutate(density_km2_adj = ifelse(lure == "Yes",
         density_km2 / lure.effect, density_km2)) %>%
  select(DeploymentYear, common_name, season, days, lure, density_km2, density_km2_adj)

# Write csv file
write_csv(df_density_wtd, paste0(public_tomel, "density-wtd-abmicmu-2015-18_2020-06-15.csv"))

#-----------------------------------------------------------------------------------------------------------------------

# Further exploration

df_density_wtd <- read_csv(paste0(public_tomel, "density-wtd-abmi-2015-19_2020-06-17.csv"))

tbd <- read_csv("G:/Shared drives/ABMI Camera Mammals/data/processed/time-by-day/abmi-all-years_tbd-summary_2020-06-08.csv")

deployments <- df_density_wtd %>%
  select(name_year = DeploymentYear) %>%
  distinct() %>%
  pull()

check1 <- df_abmi_native %>%
  filter(name_year %in% deployments,
         common_name == "White-tailed Deer")

check2 <- tbd %>%
  filter(name_year %in% deployments)

mean(check2$total_days, na.rm = TRUE)

check3 <- check1 %>%
  mutate(julian = as.numeric(format(date_detected, "%j"))) %>%
  filter(julian > 305 | julian < 168)



































