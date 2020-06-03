#-----------------------------------------------------------------------------------------------------------------------

# Title: Calculate camera deployment operating times by day
# Authors: Dave Huggard, Marcus Becker

# Previous scripts: 01_clean-raw-data

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(readr)
library(dplyr)
library(lubridate)

# Set path to Google Drive
root <- "G:/Shared drives/ABMI Camera Mammals/"

#-----------------------------------------------------------------------------------------------------------------------

# Load tag data
df_all <- read_csv(paste0(root, "data/base/clean/abmi-all-years_all-data_clean_2020-05-30.csv")) %>%
  mutate(year = ifelse(is.na(year), 2019, year))

# Camera start and end times
df_cam_range <- read_csv(paste0(root, "data/lookup/start-end/all-cam_startend_2020-06-01.csv"))

#-----------------------------------------------------------------------------------------------------------------------

# Truncate time ranges if there is a field of view issues (i.e. END without a subsequent START)
df_cam_range_trunc <- df_all %>%
  select(name, year, date_detected, common_name, field_of_view) %>%
  filter(field_of_view == "END - Last Good Image in FOV" | field_of_view == "START - First Good Image in FOV") %>%
  arrange(name, date_detected) %>%
  mutate(starts_again = ifelse(lead(field_of_view) == "START - First Good Image in FOV", 1, 0)) %>%
  filter(field_of_view == "END - Last Good Image in FOV", starts_again == "0") %>%
  select(name, year, updated_latest = date_detected)

# Update deployment operating time ranges
df_cam_range_upd <- df_cam_range %>%
  left_join(df_cam_range_trunc, by = c("name", "year")) %>%
  mutate(end_date_time = if_else(is.na(updated_latest), end_date_time, updated_latest)) %>%
  select(-updated_latest)

# Create intermediate End-Start pairs dataframe, formated as END / START in subsequent rows
inter <- df_all %>%
  filter(field_of_view == "START - First Good Image in FOV" | field_of_view == "END - Last Good Image in FOV") %>%
  group_by(name) %>%
  tally() %>%
  filter(n > 1) %>%
  select(name) %>%
  pull()

df_inter_pairs <- df_all %>%
  filter(field_of_view == "START - First Good Image in FOV" | field_of_view == "END - Last Good Image in FOV") %>%
  filter(name %in% inter) %>%
  arrange(name, year, date_detected) %>%
  select(name, year, date_detected, field_of_view) %>%
  mutate(starts_again = ifelse(lead(field_of_view) == "START - First Good Image in FOV", 1, 0)) %>%
  filter(!(field_of_view == "END - Last Good Image in FOV" & starts_again == "0")) %>%
  select(-starts_again) %>%
  group_split(field_of_view) %>%
  bind_cols() %>%
  mutate(time_diff = date_detected1 - date_detected) %>%
  filter(time_diff > 12) %>%
  select(-c(time_diff, name1))

ends <- inter_pairs %>%
  select(name, date_detected, field_of_view)

df_inter_pairs <- df_inter_pairs %>%
  select(name, date_detected = date_detected1, field_of_view = field_of_view1) %>%
  bind_rows(ends) %>%
  arrange(name, date_detected)






