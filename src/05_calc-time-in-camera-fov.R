#-----------------------------------------------------------------------------------------------------------------------

# Title: Identify series, add probabilistic gap assignment, and calculate time in front of camera
# Author: Dave Huggard, Marcus Becker

# Previous scripts: 01_clean-raw-data, 02_gap-classes, 04_camera-timeframes

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(readr)
library(dplyr)
library(tidyr)

# Set path to Google Drive
root <- "G:/Shared drives/ABMI Camera Mammals/"

#-----------------------------------------------------------------------------------------------------------------------

# Previously processed data:

# 1. Probabilistic gaps
df_leave_prob_pred <- read_csv(paste0(root, "data/processed/probabilistic-gaps/gap-leave-prob_predictions_2020-06-25.csv"))
# 2. Time between photos
df_tbp <- read_csv(paste0(root, "data/processed/time-btwn-images/Table of time between photos within series by species incl 2019 May 2020.csv")) %>%
  rename(common_name = Species)
# 3. Time by day summary
df_tbd <- read_csv(paste0(root, "data/processed/time-by-day/abmi-all-years_tbd-summary_2020-06-08.csv"))
# 4. Gap classes
df_gap <- read_csv(paste0(root, "data/processed/probabilistic-gaps/abmi-all-years_gap-class_2020-05-30.csv"))

# Lookup:

# Native species
load(paste0(root, "data/lookup/wt_native_sp.RData"))
# Gap groups
df_gap_groups <- read_csv(paste0(root, "data/lookup/species-gap-groups.csv"))

# Parameters:

# Seasonal start/end dates (julian day):
summer.start.j <- 106 # April 16
summer.end.j <- 288 # October 15

# Load tag data
df_abmi_all <- read_csv(paste0(root, "data/base/clean/abmi-all-years_all-data_clean_2020-06-08.csv"))

#-----------------------------------------------------------------------------------------------------------------------

# Step 1. Identify Series

# Filter the whole dataset for native species only, and within the camera field of view
df_abmi_native <- df_abmi_all %>%
  filter(common_name %in% native_sp,
         field_of_view == "WITHIN") %>%
  unite(col = "name_year", name, year, sep = "_", remove = FALSE)

# Save native only data
write_csv(df_abmi_native, paste0(root, "data/base/clean/abmi-all-years_native-sp_clean_", Sys.Date(), ".csv"))

df_series <- df_abmi_native %>%
  # Join gap class
  left_join(df_gap, by = c("name", "date_detected", "common_name")) %>%
  # Order observations
  arrange(name, year, date_detected, common_name) %>%
  # Identify series and gaps requiring probabilistic time assignment
  mutate(series_num = 0,
         # Lagged time stamp
         date_detected_previous = lag(date_detected),
         # Calculate difference in time between ordered images
         diff_time = as.numeric(date_detected - date_detected_previous),
         # Lagged species
         common_name_previous = lag(common_name),
         # Was is a different species?
         diff_sp = ifelse(common_name != common_name_previous, TRUE, FALSE),
         # Lagged deployment
         name_previous = lag(name),
         # Was is a different deployment?
         diff_name = ifelse(name != name_previous, TRUE, FALSE),
         # Flag gaps that will need checking
         gap_check = ifelse(diff_name == FALSE & diff_sp == FALSE & (diff_time <= 120 & diff_time >= 20), 1, 0),
         # Lagged gap class
         gap_class_previous = replace_na(lag(gap_class), ""),
         # Identify new series, based on being different deployment, species, greater than 120 seconds, and approp gaps
         diff_series = ifelse(diff_name == TRUE | diff_sp == TRUE | diff_time > 120 | (gap_class_previous == "L" | gap_class_previous == "N"), 1, 0),
         # Number series
         series_num = c(0, cumsum(diff_series[-1])),
         # Flag gaps that require probabilistic time assignment
         gap_prob = replace_na(ifelse(gap_check == 1 & (gap_class_previous == "" | gap_class_previous == "U"), 1, 0), 0)) %>%
  # Join gap group lookup table
  left_join(df_gap_groups, by = "common_name") %>%
  # Join gap leaving predictions
  left_join(df_leave_prob_pred, by = c("gap_group", "diff_time")) %>%
  # Adjust time difference between ordered images that require probabilistic time assignment
  mutate(pred = replace_na(pred, 1),
         diff_time_adj = round(ifelse(gap_prob == 1, diff_time * (1 - pred), diff_time), digits = 2))

#-----------------------------------------------------------------------------------------------------------------------

# Step 2. Calculate time between photos (tbp), by species. This is ABMI-only, but we use a greater pool of data below.

df_tbp_abmi <- df_series %>%
  mutate(series_num_previous = lag(series_num)) %>%
  # Remove first image from each series
  filter(series_num == series_num_previous) %>%
  group_by(common_name) %>%
  # Calculate average tbp and number of images from each species
  summarise(tbp = mean(diff_time),
            sample_size = n())

# Write results
write_csv(df_tbp_abmi, paste0(root, "data/processed/time-btwn-images/abmi-all-years_tbp_", Sys.Date(), ".csv"))

#-----------------------------------------------------------------------------------------------------------------------

# Step 3. Calculate total time in front of the camera, by series (tts = Total Time by Series)

# Start with series that have >1 image
df_tts_multiple <- df_series %>%
  # Remove first image from dataframe (as to not include diff_time in time totals for the series)
  mutate(series_num_previous = lag(series_num)) %>%
  filter(series_num_previous == series_num) %>%
  group_by(series_num) %>%
  summarise(n_images = n(),
            total_time = sum(diff_time_adj))

# Next, single-image series'
df_tts_single <- df_series %>%
  group_by(series_num) %>%
  summarise(n_images = n()) %>%
  # Keep only series with 1 image
  filter(n_images == 1) %>%
  mutate(total_time = 0)

# Bind together
df_tts_all <- df_tts_multiple %>%
  bind_rows(df_tts_single) %>%
  arrange(series_num)

# Add time between photos, accounting for the average number of animals in each photo of a series
df_tts_final <- df_series %>%
  # Join in average tbp, by species; note that this is slightly different from results calculated above (incl. non-ABMI)
  left_join(df_tbp, by = "common_name") %>%
  select(series_num, common_name, number_individuals, time_btwn_images = TBP) %>%
  group_by(series_num) %>%
  # Number of individuals in a series
  mutate(avg_individuals = mean(number_individuals)) %>%
  ungroup() %>%
  select(-number_individuals) %>%
  distinct() %>%
  left_join(df_tts_all, by = "series_num") %>%
  mutate(series_total_time = (total_time + time_btwn_images) * avg_individuals) %>%
  select(series_num, common_name, n_images, series_total_time)

#-----------------------------------------------------------------------------------------------------------------------

# Step 4. Calculate total time in front of camera, by deployment, year, and species (tt = total time)

df_tt <- df_series %>%
  group_by(series_num) %>%
  arrange(date_detected) %>%
  filter(row_number() == 1) %>%
  left_join(df_tts_final, by = c("series_num", "common_name")) %>%
  select(series_num, name_year, name, year, date_detected, common_name, series_total_time) %>%
  ungroup() %>%
  mutate(julian = as.numeric(format(date_detected, "%j")),
         season = ifelse(julian >= summer.start.j & julian <= summer.end.j, "summer", "winter")) %>%
  mutate_at(c("name_year", "common_name", "season"), factor) %>%
  group_by(name_year, common_name, season, .drop = FALSE) %>%
  summarise(total_duration = sum(series_total_time)) %>%
  ungroup() %>%
  mutate_if(is.factor, as.character) %>%
  left_join(df_tbd, by = "name_year")

sp <- as.character(sort(unique(df_tt$common_name)))

# For deployments with no images of native animals (nn = no natives)

# Vector of all deployments in the ABMI projects:
dep <- df_abmi_all %>%
  select(name, year) %>%
  distinct() %>%
  unite(col = "name_year", name, year, sep = "_", remove = TRUE) %>%
  pull()

df_tt_nn <- df_tbd %>%
  # Filter out non-ABMI deployments in df_tbd
  filter(name_year %in% dep) %>%
  # Retrieve only those that had no images of native species
  anti_join(df_tt, by = "name_year") %>%
  expand(name_year, season = c("summer", "winter"), common_name = sp) %>%
  # Re-join time-by-day information
  left_join(df_tbd, by = "name_year") %>%
  # Add total_duration column, which is zero in these cases
  mutate(total_duration = 0)

df_tt_full <- df_tt %>%
  bind_rows(df_tt_nn) %>%
  arrange(name_year, common_name, season) %>%
  mutate(total_season_days = ifelse(season == "summer", total_summer_days, total_winter_days)) %>%
  select(name_year:season, total_season_days, total_duration)

#-----------------------------------------------------------------------------------------------------------------------

# Write results

write_csv(df_tt_full, paste0(root, "data/processed/time-in-cam-fov/abmi-all-years_cam-time_", Sys.Date(), ".csv"))

#-----------------------------------------------------------------------------------------------------------------------
















