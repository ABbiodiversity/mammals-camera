#-----------------------------------------------------------------------------------------------------------------------

# Title: Modeling leave probability based on gap tags
# Author: Dave Huggard, Marcus Becker

# Previous scripts: 01_clean-raw-data, 02_gap-classes

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(readr)
library(tidyr)
library(dplyr)
library(purrr)

# Set path to Google Drive
root <- "G:/Shared drives/ABMI Camera Mammals/"

#-----------------------------------------------------------------------------------------------------------------------

# Previously processed data
df_gap <- read_csv(paste0(root, "data/processed/probabilistic-gaps/abmi-all-years_gap-class_2020-05-30.csv"))

# Lookup table
df_gap_groups <- read_csv(paste0(root, "data/lookup/species-gap-groups.csv"))

# Load tag data
df_abmi_all <- read_csv(paste0(root, "data/base/clean/abmi-all-years_all-data_clean_2020-06-08.csv"))

#-----------------------------------------------------------------------------------------------------------------------

# Filter the whole dataset for native species only, and within the camera field of view
df_abmi_native <- df_abmi_all %>%
  filter(common_name %in% native_sp,
         field_of_view == "WITHIN") %>%
  unite(col = "name_year", name, year, sep = "_", remove = FALSE)

#-----------------------------------------------------------------------------------------------------------------------

# Probabilistic Gaps

# Set up data for modeling leaving probability using observations from gap analysis:
df_leave_prob <- df_abmi_native %>%
  # Join gap class
  left_join(df_gap, by = c("name", "date_detected", "common_name")) %>%
  # Order observations
  arrange(name, year, date_detected, common_name) %>%
  # Identify series and gaps requiring probabilistic time assignment
  # Lagged time stamp
  mutate(date_detected_previous = lag(date_detected),
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
         gap_class_previous = lag(gap_class)) %>%
  # Join gap group lookup table
  left_join(df_gap_groups, by = "common_name") %>%
  # Isolate observations requiring a gap check
  filter(!is.na(gap_class_previous) & gap_check == "1") %>%
  select(common_name, gap_group, gap_class_previous, diff_time) %>%
  # Determine whether the animal left the field of view or not (i.e. tagged as 'P' (present) by gap analysis)
  mutate(left = ifelse(gap_class_previous == "P", 0, 1)) %>%
  select(-gap_class_previous)

# Modeling by gap group, followed by predictions using those models:
df_leave_prob_pred <- df_leave_prob %>%
  filter(!is.na(gap_group)) %>%
  group_by(gap_group) %>%
  nest() %>%
  # Define model for each gap group, and then make predictions
  mutate(model = map(.x = data, ~ smooth.spline(x = .$diff_time, y = .$left, df = 3)),
         pred = map(.x = model, ~ predict(., x = 20:120))) %>%
  select(gap_group, pred) %>%
  unnest_wider(pred) %>%
  unnest(cols = c(x, y)) %>%
  rename(diff_time = x, pred = y) %>%
  ungroup()

#-----------------------------------------------------------------------------------------------------------------------

# Write results

write_csv(df_leave_prob, paste0(root, "data/processed/probabilistic-gaps/gap-leave-prob_raw-data_", Sys.Date(), ".csv"))

write_csv(df_leave_prob_pred, paste0(root, "data/processed/probabilistic-gaps/gap-leave-prob_predictions_", Sys.Date(), ".csv"))

#-----------------------------------------------------------------------------------------------------------------------




