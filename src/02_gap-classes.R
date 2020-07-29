#-----------------------------------------------------------------------------------------------------------------------

# Title: Retrieve gap class designations from old ABMI camera tagging system, add N gap classes following NONE images
# Authors: Dave Huggard, Marcus Becker

# Previous scripts: 01_clean-raw-data

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(readr)
library(dplyr)
library(lubridate)

# Load native species list
load(paste0(root, "data/lookup/wt_native_sp.RData"))

# Set path to Google Drive
root <- "G:/Shared drives/ABMI Camera Mammals/"

#-----------------------------------------------------------------------------------------------------------------------

# Retrieve gap class designations from previous data (years 2013 through 2018)

# Note: These tags were added with previous ABMI tagging system.

df_gap <-
  # Read in previous data (years 2013 through 2018)
  read_csv(paste0(root, "data/base/raw/previous/ALL_native-mammals_2019-09-01.csv"),
                         col_types = cols(distance = col_character(),
                                          number_during_gap = col_number(),
                                          number_individuals = col_character()),
                         na = "") %>%
  filter(!is.na(gap_class)) %>%
  mutate(date_detected = ymd_hms(date_time_taken)) %>%
  select(name = deployment, date_detected, common_name, gap_class)

#-----------------------------------------------------------------------------------------------------------------------

# Add 'N' gap classes for images following a 'NONE' image for ABMI 2019

# Note: a 'NONE' image is used to demarcate when a series should be truncated because an animal left the field of view.

df_abmi_2019 <- read_csv(paste0(root, "data/base/clean/abmi-2019_all-data_clean_2020-06-02.csv"))

df_abmi_2019_ngap <- df_abmi_2019 %>%
  select(name, date_detected, common_name) %>%
  arrange(name, date_detected) %>%
  # Create gap_class column
  mutate(common_name_next = lead(common_name),
         gap_class = ifelse(common_name != "NONE" & common_name_next == "NONE", "N", NA)) %>%
  # Include only N gap class for native mammals
  filter(gap_class == "N",
         common_name %in% native_sp) %>%
  select(-common_name_next)

# Combine with df_gap, write to csv
df_gap <- bind_rows(df_gap, df_abmi_2019_ngap) %>%
  write_csv(paste0(root, "data/processed/probabilistic-gaps/abmi-all-years_gap-class_", Sys.Date(), ".csv"))

#-----------------------------------------------------------------------------------------------------------------------
