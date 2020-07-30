#-------------------------------------------------------------------------------

# Title: Summarise Camera Operating Time
# Created: October 1, 2019
# Author: Dave Huggard, Marcus Becker

# Objectives: Calculate the amount of time in days that each deployment has been
#             operating, broken down by season.

#-------------------------------------------------------------------------------

# Load packages

library(tidyverse)
library(lubridate)
library(mgcv)

#-------------------------------------------------------------------------------

# Import data

# Path to data folder on ABMI Science Centre S: drive
abmisc <- "S:/github-repos-data/SC-Camera-Mammals/data/"

# Deployments with images of native mammals
df_dep_native <- read_csv(paste0(abmisc,"processed/All deployments with images July 2019 ALL.csv")) %>%
  mutate(DeploymentYear = paste(Deployment, Year, sep = "_"))

# All deployments start and end times (regardless of whether images were captured)
df_cam_time <- read_csv(paste0(abmisc,"base/Start end times Lure Feb 2019 ALL.csv")) %>%
  # Filter out NAs and null values - flag as something to come back to.
  filter(!is.na(EndTime) & EndTime != "#VALUE!" & EndTime != "#N/A")

# Intermediate end/start pairs (cameras w/ two operating periods)
df_cam_time_inter <- read_csv(paste0(abmisc,"base/Intermediate END START pairs Feb 2019 ALL.csv"))

#-------------------------------------------------------------------------------

# Parameters

# Calculate total days from Jan 1, 2009 to the latest deployment end time.
start <- as.Date("2009-01-01")
end <- max(as.Date(df_cam_time$EndTime))

interval <- start %--% end
days <- as.duration(interval) / ddays(1)

#-------------------------------------------------------------------------------

# Clean data

df_cam_time <- df_cam_time %>%
  mutate(StartTime = ymd_hms(StartTime),
         EndTime = ymd_hms(EndTime),
         DeploymentYear = paste(Deployment, Year, sep = "_")) %>%
  # Join df_dep_native - use semi_join to just keep matching values.
  semi_join(df_dep_native, by = "DeploymentYear")

df_cam_time_inter <- df_cam_time_inter %>%
  mutate(StartEndTime1 = ymd_hms(date_time_original),
         DeploymentYear = paste(Deployment, Year, sep = "_")) %>%
  semi_join(df_dep_native, by = "DeploymentYear")

# Create vector of DeploymentYear's
sitelist <- unique(as.character(df_cam_time$DeploymentYear))

# Create array to track each Julian day that cameras were operating
time.by.day<-array(0,c(length(sitelist),366))

# Create array to track each day since 2009 that cameras were operating
time.since.2009 <- array(0, c(length(sitelist), days))

# Populate array
for (i in 1:length(sitelist)) {
  df_cam_time_1 <- df_cam_time[df_cam_time$DeploymentYear == sitelist[i],]
  for (k in 1:nrow(df_cam_time_1)) {
    # For days since 2009. floor() is used so that first & last day of operation is included.
    j1 <- floor(julian(df_cam_time_1$StartTime[k],
                       origin = strptime("2009-01-01 00:00:00",
                                         format="%Y-%m-%d %H:%M:%S")))
    j2 <- floor(julian(df_cam_time_1$EndTime[k],
                       origin = strptime("2009-01-01 00:00:00",
                                         format="%Y-%m-%d %H:%M:%S")))
    if (!is.na(j1) & !is.na(j2)) {
      time.since.2009[i, (j1:j2)] <- 1
    }
  }
  # To take off time(s) when camera wasn't working.
  if (sitelist[i] %in% df_cam_time_inter$DeploymentYear) {
    df_inter_1 <- df_cam_time_inter[as.character(df_cam_time_inter$DeploymentYear) == as.character(sitelist[i]),]
    for (j in seq(1, (nrow(df_inter_1) - 1), 2)) { # Assumes all extra times are formatted as END/START pairs
      # Use ceiling() so that day of failure is excluded from operating days
      j1 <- ceiling(julian(df_inter_1$StartEndTime1[j],
                           origin = strptime("2009-01-01 00:00:00",
                                             format="%Y-%m-%d %H:%M:%S")))
      j2 <- ceiling(julian(df_inter_1$StartEndTime1[j],
                           origin = strptime("2009-01-01 00:00:00",
                                             format="%Y-%m-%d %H:%M:%S")))
      if (j2 > j1) time.since.2009[i, j1:(j2-1)] <- 0
    }
  }
}

# Then calculate time.by.day from time.since.2009

days.per.year <- c(365,365,365,366,365,365,365,366,365,365,365)

Jan1 <- cumsum(c(1, days.per.year))

yrnum <- julday <- NULL

for (i in 1:ncol(time.since.2009)) {
  yrnum[i] <- sum(i >= Jan1)
  julday[i] <- i - Jan1[yrnum[i]] + 1
}
for (i in 1:ncol(time.by.day)) {
  time.by.day[,i] <- rowSums(time.since.2009[,which(julday == i)])
}

rownames(time.by.day) <- rownames(time.since.2009) <- sitelist

# Export data

time.by.day %>% as_tibble(rownames = "DeploymentYear") %>%
  write_csv(paste0(abmisc,"processed/Camera operating days Feb 2019 ALL.csv"))

time.since.2009 %>% as_tibble(rownames = "DeploymentYear") %>%
  write_csv(paste0(abmisc,"processed/Camera operating days since 2009 Feb 2019 ALL.csv"))

# A bit more processing.

df_tbd_summary <- time.by.day %>%
  as_tibble(rownames = "DeploymentYear") %>%
  mutate(total = rowSums(select(., -DeploymentYear)),
         total.summer = rowSums(select(., -c(DeploymentYear:V106, V289:V366))),
         total.winter = rowSums(select(., -c(DeploymentYear, V107:V288)))) %>%
  select(DeploymentYear, total, total.summer, total.winter)

write_csv(df_tbd_summary, paste0(abmisc,"processed/Summary of camera operating days Feb 2019 ALL.csv")





