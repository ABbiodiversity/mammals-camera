# ------------------------------------------------------------------------------

# Title: Build dataframe of all deployments with images
# Created: July 12, 2019
# Author: Marcus Becker

# Objectives: Read in camera data from all years and all sources, combine in one
#             dataframe for further downstream processing. Compare to Dave H's
#             original df of this.

# ------------------------------------------------------------------------------

# Load packages
library(tidyverse)
library(fs)
library(vroom)

# Build relative path if accessing data through S: drive
# s_drive <- "S:/github-repos-data/SC-Mammals/"
# setwd(s_drive)

# Import data
paths <- dir_ls("data/base/raw/", glob = "*.csv")
names <- c("original", "abmi1314", "abmi15", "abmi16", "abmi17", "abmi18",
           "cmu17", "cmu18", "other1", "other2")

ls_cam_data <- map(.x = paths, vroom::vroom) %>% # "vroom vroom" - Mazda, probably.
  # "original" is Dave H's df
  set_names(names)

# Lookup table for CMU deployment abbreviations
cmu_abb <- read_csv("./data/lookup/cmu-dep-abb.csv") %>%
  pull(abb)

#-------------------------------------------------------------------------------

# Process and compare data, build df of all deployments with images

# 2013-14 ABMI

df_orig_1314 <- ls_cam_data[["original"]] %>%
  filter(Year == "2013" | Year == "2014") %>%
  filter(str_detect(Deployment, "ABMI")) # 40 deployments

df_1314 <- ls_cam_data[["abmi1314"]] %>%
  # filter(tagged == TRUE) %>% I guess we're okay with time lapse images after all?
  group_by(name, deployment_year) %>%
  summarise(Images = as.numeric(n())) %>%
  # Some wonky deployment names
  mutate(Deployment = gsub("_", "-", name)) %>%
  separate(Deployment, into = c("Dep_Who", "Dep_Code", "Dep_Corner"), sep = "-", remove = FALSE) %>%
  mutate(Dep_Corner = case_when(
    Dep_Corner == 3 ~ "SE",
    Dep_Corner == 5 ~ "NE",
    Dep_Corner == 7 ~ "NW",
    Dep_Corner == 9 ~ "SW",
    TRUE ~ Dep_Corner)) %>%
  mutate(Deployment = paste(Dep_Who, Dep_Code, Dep_Corner, sep = "-")) %>%
  ungroup() %>%
  select(Deployment, Year = deployment_year, Images)

all_equal(df_orig_1314, df_1314) # Same.

rm(df_orig_1314)

# 2015 ABMI

df_orig_15 <- ls_cam_data[["original"]] %>%
  filter(Year == "2015") %>%
  filter(str_detect(Deployment, "ABMI") | str_detect(Deployment, "OG-OGC")) # 648 deployments

df_15 <- ls_cam_data[["abmi15"]] %>%
  group_by(Name, deployment_year) %>%
  summarise(Images = as.numeric(n())) %>%
  rename(Deployment = Name, Year = deployment_year) %>%
  ungroup()

all_equal(df_orig_15, df_15) # Same.

rm(df_orig_15)

# 2016 ABMI

df_orig_16 <- ls_cam_data[["original"]] %>%
  filter(Year == "2016") %>%
  filter(str_detect(Deployment, "ABMI|OG-CITSCI")) # 647 deployments

df_16 <- ls_cam_data[["abmi16"]] %>%
  filter(tagged == TRUE) %>% # Suddenly we aren't okay with lapse images ...
  group_by(name, deployment_year) %>%
  summarise(Images = as.numeric(n())) %>%
  rename(Deployment = name, Year = deployment_year) %>%
  ungroup()

all_equal(df_16, df_orig_16) # Same.

rm(df_orig_16)

# 2017 ABMI

df_orig_17 <- ls_cam_data[["original"]] %>%
  filter(Year == "2017") %>%
  filter(str_detect(Deployment, "ABMI|Reconyx|OG-AAC|OG-EI")) # 917 deployments

df_17 <- ls_cam_data[["abmi17"]] %>%
  filter(tagged == TRUE) %>% # Weird.
  group_by(name, deployment_year) %>%
  summarise(Images = as.numeric(n())) %>%
  rename(Deployment = name, Year = deployment_year) %>%
  ungroup()

all_equal(df_orig_17, df_17) # Same.

rm(df_orig_17)

# 2018 ABMI

df_orig_18 <- ls_cam_data[["original"]] %>%
  filter(Year == "2018") %>%
  filter(str_detect(Deployment, "ABMI|RIVR")) # 819 deployments

df_18 <- ls_cam_data[["abmi18"]] %>%
  filter(tagged == TRUE) %>%
  group_by(name, deployment_year) %>%
  summarise(Images = as.numeric(n())) %>%
  rename(Deployment = name, Year = deployment_year) %>%
  ungroup()

all_equal(df_orig_18, df_18) # Same.

rm(df_orig_18)

# CMU 2017

df_orig_cmu17 <- ls_cam_data[["original"]] %>%
  filter(Year == "2017") %>%
  filter(str_detect(Deployment, paste(cmu_abb, collapse = "|"))) # 145 deployments

df_cmu17 <- ls_cam_data[["cmu17"]] %>%
  # filter(tagged == TRUE) %>% What? Not okay with them again?
  mutate(Deployment = gsub("-", "", name)) %>%
  group_by(Deployment, deployment_year) %>%
  summarise(Images = as.numeric(n())) %>%
  rename(Year = deployment_year) %>%
  ungroup()

all_equal(df_orig_cmu17, df_cmu17) # Same.

rm(df_orig_cmu17)

# CMU 2018

df_orig_cmu18 <- ls_cam_data[["original"]] %>%
  filter(Year == "2018") %>%
  filter(str_detect(Deployment, paste(cmu_abb, collapse = "|"))) # 270 deployments

df_cmu18 <- ls_cam_data[["cmu18"]] %>%
  # filter(tagged == TRUE) %>%
  group_by(name, deployment_year) %>%
  summarise(Images = as.numeric(n())) %>%
  ungroup() %>%
  rename(Year = deployment_year, Deployment = name) %>%
  mutate(Deployment = gsub("-", "", Deployment))

all_equal(df_cmu18, df_orig_cmu18) # Same.

# ^ This was a weird one for naming issues. For example, a deployment may be named FMM-1 and FMM1.
# When we get rid of the hyphen earlier, this aggregates these 'two' deployments (are they really two?)
# I think they're just one deployment.

rm(df_orig_cmu18)

# Other data 1 - Various ABMI OGW and RIVR

df_orig_other1 <- ls_cam_data[["original"]] %>%
  filter(Year == "2018" | Year == "2016") %>%
  filter(str_detect(Deployment, "OGW-CITSCI|OG-RIVR")) # 41 deployments

df_other1 <- ls_cam_data[["other1"]] %>%
  filter(tagged == TRUE) %>%
  group_by(name, deployment_year) %>%
  summarise(Images = as.numeric(n())) %>%
  ungroup() %>%
  rename(Year = deployment_year, Deployment = name)

# Let's get rid of the RIVR deployments - duplicated in 2016 ABMI data.

df_other1 <- df_other1 %>%
  filter(!str_detect(Deployment, "OG-RIVR"))

all_equal(df_other1, df_orig_other1) # Not same anymore.

# ^ This is weird too. OG-RIVR were already included in the 2016 ABMI data, but OGW-CITSCI were not included in the 2018 ABMI data.

rm(df_orig_other1)

# Other data 2

df_orig_other2 <- ls_cam_data[["original"]] %>%
  filter(str_detect(Deployment, "BG|NEXEN|WHEC|TSP"))

df_other2 <- ls_cam_data[["other2"]] %>%
  filter(tagged == TRUE) %>%
  group_by(name, deployment_year) %>%
  summarise(Images = as.numeric(n())) %>%
  ungroup() %>%
  rename(Year = deployment_year, Deployment = name) %>%
  mutate(Deployment = ifelse(str_detect(Deployment, "^WHEC"), tolower(Deployment), Deployment),
         Deployment = gsub("whec", "WHEC", Deployment))

all_equal(df_orig_other2, df_other2) # Same

rm(df_orig_other2)

#-------------------------------------------------------------------------------

# Unite them all under a common banner

df_dep_img <- bind_rows(
  df_1314,
  df_15,
  df_16,
  df_17,
  df_18,
  df_cmu17,
  df_cmu18,
  df_other1,
  df_other2
)

df_dep_img_orig <- ls_cam_data[["original"]]

all_equal(df_dep_img, df_dep_img_orig) # TRUE

write_csv(df_dep_img, "./data/processed/All deployments with images July 2019 ALL.csv")




