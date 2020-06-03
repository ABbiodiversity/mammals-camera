#-----------------------------------------------------------------------------------------------------------------------

# Title: Download and clean raw ABMI data from WildTrax in preparation for downstream density estimation.
# Author: Marcus Becker

# Previous scripts: None

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
library(lubridate)
library(keyring)

# Load functions
source("./src/functions/abmi-camera-mammals-fns.R")

# Set path to Google Drive
root <- "G:/Shared drives/ABMI Camera Mammals/"

#-----------------------------------------------------------------------------------------------------------------------

# Download raw data

# 1. Establish connection to WildTrax database

connection <- wt_conn(username = key_get("username", keyring = "wildtrax"),
                      password = key_get("password", keyring = "wildtrax"))

keyring_lock("wildtrax")

# 2. Query the database for ABMI camera data

project_ids <- c(205, 194, 197, 193, 195, 214)

# Query raw data:
df_abmi_all_raw <- map_df(.x = project_ids,
                      .f = ~ wt_cam_report_tag_summary(conn = connection, proj_id = .x, native_only = FALSE))

# Disconnect from the database:
dbDisconnect(connection)

#-----------------------------------------------------------------------------------------------------------------------

# ABMI 2019

# 1. Standardize deployment names
# 2. Remove duplicated deployments
# 3. Combine tags of same species species in same image

duplicates <- read_csv(paste0(root, "data/2019-temporary/abmi-2019-removed-duplicates.csv"))

df_abmi_2019 <- df_abmi_all_raw %>%
  filter(project == "ABMI 2019") %>%
  # Fix year variable
  mutate(year = 2019) %>%
  # Clean deployment names
  mutate(name = ifelse(str_detect(name, "OG"), name, str_remove(name, "ABMI-")),
         alt_cam_mod = case_when(
           str_detect(name, "HF") ~ "HF2",
           str_detect(name, "HP") ~ "HP2X",
           str_detect(name, "CUDDE") ~ "CUDDE",
           TRUE ~ "Original"),
         name = ifelse(str_detect(name, "HF|CUDDE|HP"), str_remove(name, "-HF2|-CUDDE|-HP2X"), name),
         alt_cam_mod = ifelse(str_detect(name, "-1$"), paste0(alt_cam_mod, "-1"), alt_cam_mod),
         name = ifelse(str_detect(name, "-1$") & !str_detect(name, "OG"), str_sub(name, end = -3), name),
         name = str_pad(name, width = 7, side = "left", pad = "0"),
         name = ifelse(str_detect(name, "OG"), name, paste0("ABMI-", name)),
         name = ifelse(str_detect(alt_cam_mod, "CUDDE|HF|HP"), paste0(name, "_", alt_cam_mod), name),
         name = ifelse(name == "ABMI-638-61-2", "OG-ABMI-638-61-2", name)) %>%
  select(-alt_cam_mod) %>%
  # Remove duplicated deployments
  mutate(id_by = str_trim(id_by, side = "both")) %>%
  anti_join(duplicates, by = c("name", "id_by")) %>%
  # Remove duplicated observations
  distinct() %>%
  # Amalgamate tags of same species in same image
  mutate(number_individuals = as.numeric(ifelse(number_individuals == "VNA", 1, number_individuals))) %>%
  group_by(name, date_detected, common_name) %>%
  mutate(number_individuals = sum(number_individuals),
         age_class = paste0(age_class, collapse = ", "),
         sex = paste0(sex, collapse = ", ")) %>%
  distinct(name, date_detected, common_name, number_individuals, .keep_all = TRUE) %>%
  ungroup() %>%
  # Fix wrong year in 701-61-89
  mutate(date_detected = ymd_hms(date_detected),
         date_detected = if_else(name == "OG-ABMI-701-61-89", date_detected %m-% years(2), date_detected))

# Save cleaned data
write_csv(df_abmi_2019, paste0(root, "data/base/clean/abmi-2019_all-data_clean_", Sys.Date(), ".csv"))

#-----------------------------------------------------------------------------------------------------------------------

# ABMI All Years
df_abmi_allyears <- df_abmi_all_raw %>%
  filter(!project == "ABMI 2019") %>%
  mutate(date_detected = ymd_hms(date_detected),
         number_individuals = as.numeric(ifelse(number_individuals == "VNA", 1, number_individuals))) %>%
  # Amalgamate tags of same species in same image
  group_by(name, date_detected, common_name) %>%
  mutate(number_individuals = sum(number_individuals),
         age_class = paste0(age_class, collapse = ", "),
         sex = paste0(sex, collapse = ", ")) %>%
  distinct(name, date_detected, common_name, number_individuals, .keep_all = TRUE) %>%
  ungroup()

# Join the years together
df_abmi_allyears <- df_abmi_allyears %>%
  # Rejoin the cleaned up version of ABMI 2019
  bind_rows(df_abmi_2019) %>%
  # Keep only necessary columns for downstream density estimation:
  select(project, name, year, date_detected, common_name, number_individuals, age_class, sex, field_of_view)

# Save cleaned data
write_csv(df_abmi_allyears, paste0(root, "data/base/clean/abmi-all-years_all-data_clean_", Sys.Date(), ".csv"))

#-----------------------------------------------------------------------------------------------------------------------
















