#-----------------------------------------------------------------------------------------------------------------------

# Title: Adjust density estimates by lure status
# Author: Dave Huggard, Marcus Becker

# Previous scripts: 01_clean-raw-data, 02_gap-classes, 04_camera-timeframes, 05_calc-time-in-camera-fov

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(readr)
library(dplyr)
library(stringr)

# Set path to Google Drive
root <- "G:/Shared drives/ABMI Camera Mammals/"

# Lure lookup
df_lure <- read_csv(paste0(root, "data/lookup/lure/all-cam-lure_2020-06-01.csv"))

# Read in density data:
df_density <- read_csv(paste0(root, "results/density/Marcus/abmi-all-years_density_2020-06-09.csv"))

#-----------------------------------------------------------------------------------------------------------------------

# Filter for only ABMI core site estimates:
df_density_abmi_core <- df_density %>%
  filter(str_detect(name_year, "ABMI"),
         !str_detect(name_year, "OG"))

# Calculate lure:unlured density ratios:
df_lure_ratios <- df_lure %>%
  unite(col = "name_year", name, year, sep = "_") %>%
  right_join(df_density_abmi_core, by = "name_year") %>%
  group_by(common_name, lure) %>%
  filter(!is.na(lure),
         !density_km2 == "Inf") %>%
  summarise(avg_density = mean(density_km2, na.rm = TRUE)) %>%
  summarise(ratio = avg_density[lure == "Yes"] / avg_density[lure == "No"])

# Adjust density estimates:
df_density_adj <- df_lure %>%
  unite(col = "name_year", name, year, sep = "_") %>%
  right_join(df_density, by = "name_year") %>%
  left_join(df_lure_ratios, by = "common_name") %>%
  mutate(density_km2_lure_adj = ifelse(lure == "Yes", density_km2 / ratio, density_km2)) %>%
  select(-ratio)

#-----------------------------------------------------------------------------------------------------------------------

# Write results

write_csv(df_lure_ratios, paste0(root, "data/processed/lure/species-lure-ratios_", Sys.Date(), ".csv"))

write_csv(df_density_adj, paste0(root, "results/density/Marcus/abmi-all-years_density-lure-adj_", Sys.Date(), ".csv"))

#-----------------------------------------------------------------------------------------------------------------------




