#-----------------------------------------------------------------------------------------------------------------------

# Title: Calculate density for each camera deployment
# Author: Dave Huggard, Marcus Becker

# Previous scripts: 01_clean-raw-data, 02_gap-classes, 04_camera-timeframes, 05_calc-time-in-camera-fov

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)

# Set path to Google Drive
root <- "G:/Shared drives/ABMI Camera Mammals/"

#-----------------------------------------------------------------------------------------------------------------------

# Load time in front of camera data:
df_tt_full <- read_csv(paste0(root, "data/processed/time-in-cam-fov/abmi-all-years_cam-time_2020-06-25.csv"))

# Load effective detection distance (EDD) modeling:
load(paste0(root,"data/processed/detection-distance/predictions/Detection distances by site species and season.rdata"))

# Lookup:

df_dist_groups <- read_csv(paste0(root, "data/lookup/species-distance-groups.csv"))

# Set parameters:

# Camera field of view angle
cam_fov_angle <- 42

#-----------------------------------------------------------------------------------------------------------------------

# Step 1. Append EDD information

df_detdist <- dd %>%
  as.data.frame() %>%
  rownames_to_column(var = "name_year") %>%
  gather(key = "SpGroupSeason", value = "detdist", Bear.Summer:`Bighorn sheep.Winter`) %>%
  mutate(SpGroupSeason = str_replace_all(SpGroupSeason, "[//(//)// ]", "")) %>%
  # Create two new columns: Detection Distance Group and Season, sep by "."
  separate(SpGroupSeason, into = c("dist_group", "season")) %>%
  mutate(dist_group = str_replace(dist_group, "wapiti", ""),
         season = tolower(season))

df_dens_ing <- df_tt_full %>%
  left_join(df_dist_groups, by = "common_name") %>%
  left_join(df_detdist, by = c("name_year", "dist_group", "season"))

#-----------------------------------------------------------------------------------------------------------------------

# Step 2. Calculate density

df_density <- df_dens_ing %>%
  mutate(effort = total_season_days * (detdist ^ 2 * pi * (cam_fov_angle / 360)) / 100,
         # Catch per unit effort
         cpue = total_duration / effort,
         # captch per unit effort in km2
         cpue_km2 = cpue / 60 / 60 / 25 * 10000) %>%
  select(name_year:total_season_days, density_km2 = cpue_km2)

#-----------------------------------------------------------------------------------------------------------------------

# Write results

write_csv(df_density, paste0(root, "results/density/Marcus/abmi-all-years_density_", Sys.Date(), ".csv"))

#-----------------------------------------------------------------------------------------------------------------------
