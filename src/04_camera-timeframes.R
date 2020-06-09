#-----------------------------------------------------------------------------------------------------------------------

# Title: Calculate camera deployment operating times by day
# Authors: Dave Huggard, Marcus Becker

# Previous scripts: 01_clean-raw-data

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(readr)
library(dplyr)
library(lubridate)
library(tibble)

# Set path to Google Drive
root <- "G:/Shared drives/ABMI Camera Mammals/"

#-----------------------------------------------------------------------------------------------------------------------

# Load tag data
df_all <- read_csv(paste0(root, "data/base/clean/abmi-all-years_all-data_clean_2020-06-02.csv"))

# Camera start and end times
df_cam_range <- read_csv(paste0(root, "data/lookup/start-end/all-cam_startend_2020-06-01.csv"))

#-----------------------------------------------------------------------------------------------------------------------

# Temporary: fix number of individuals issue. (to be moved to script 01_clean-raw-data when WT query works again)
df_all <- df_all %>%
  group_by(name, date_detected, common_name) %>%
  mutate(number_individuals = sum(number_individuals)) %>%
  distinct(name, date_detected, common_name, number_individuals, .keep_all = TRUE) %>%
  ungroup()

write_csv(df_all, paste0(root, "data/base/clean/abmi-all-years_all-data_clean_", Sys.Date(), ".csv"))

#-----------------------------------------------------------------------------------------------------------------------

# Truncate time ranges if there is a field of view issues (i.e. END without a subsequent START):
df_cam_range_trunc <- df_all %>%
  select(name, year, date_detected, common_name, field_of_view) %>%
  filter(field_of_view == "END - Last Good Image in FOV" | field_of_view == "START - First Good Image in FOV") %>%
  group_by(name) %>%
  arrange(date_detected) %>%
  mutate(starts_again = ifelse(lead(field_of_view) == "START - First Good Image in FOV", 1, 0)) %>%
  filter(field_of_view == "END - Last Good Image in FOV", is.na(starts_again)) %>%
  select(name, year, updated_latest = date_detected)

# Update deployment operating time ranges with new end_date_time if applicable:
df_cam_range_upd <- df_cam_range %>%
  left_join(df_cam_range_trunc, by = c("name", "year")) %>%
  mutate(end_date_time = if_else(is.na(updated_latest), end_date_time, updated_latest),
         name_year = paste0(name, "_", year)) %>%
  select(name_year, start_date_time, end_date_time)

# Create intermediate End-Start pairs dataframe, formated as END / START in subsequent rows:
inter <- df_all %>%
  filter(field_of_view == "START - First Good Image in FOV" | field_of_view == "END - Last Good Image in FOV") %>%
  mutate(name_year = paste0(name, "_", year)) %>%
  group_by(name_year) %>%
  tally() %>%
  filter(n > 1) %>%
  select(name_year) %>%
  pull()

df_inter_pairs <- df_all %>%
  filter(field_of_view == "START - First Good Image in FOV" | field_of_view == "END - Last Good Image in FOV") %>%
  mutate(name_year = paste0(name, "_", year)) %>%
  filter(name_year %in% inter) %>%
  arrange(name_year, date_detected) %>%
  select(name_year, date_detected, field_of_view) %>%
  group_by(name_year) %>%
  # This code is gross.
  mutate(starts_again = ifelse(lead(field_of_view) == "START - First Good Image in FOV" & field_of_view == "END - Last Good Image in FOV", 1, NA),
         restart = ifelse(lag(starts_again) == "1" & lag(field_of_view) == "END - Last Good Image in FOV", 1, NA)) %>%
  filter(starts_again == "1" | restart == "1") %>%
  select(-c(starts_again, restart)) %>%
  ungroup() %>%
  group_split(field_of_view) %>%
  bind_cols() %>%
  mutate(time_diff = difftime(date_detected1, date_detected, units = "hours")) %>%
  filter(time_diff > 12) %>%
  select(-c(time_diff, name_year1))

ends <- df_inter_pairs %>%
  select(name_year, date_detected, field_of_view)

df_inter_pairs <- df_inter_pairs %>%
  select(name_year, date_detected = date_detected1, field_of_view = field_of_view1) %>%
  bind_rows(ends) %>%
  arrange(name_year, date_detected)

#-----------------------------------------------------------------------------------------------------------------------

start <- as.Date("2009-01-01")
end <- df_cam_range_upd %>% filter(!is.na(end_date_time)) %>% pull(end_date_time) %>% max()
interval <- start %--% end
days <- ceiling(as.duration(interval) / ddays(1))

# Create vector of deployments
dep <- df_cam_range_upd %>% filter(!is.na(end_date_time)) %>% pull(name_year)

# Date ranges, no NAs
ranges <- df_cam_range_upd %>% filter(!is.na(end_date_time))

# Build arrays
time.by.day <- array(0, c(length(dep), 366))
time.since.2009 <- array(0, c(length(dep), days))

# Populate arrays
for (i in 1:length(dep)) {
  df <- ranges[ranges$name_year == dep[i],]
  for (k in 1:nrow(df)) {
    # For days since 2009. floor() is used so that first & last day of operation is included.
    j1 <- floor(julian(df$start_date_time[k], origin = strptime("2009-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")))
    j2 <- floor(julian(df$end_date_time[k], origin = strptime("2009-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")))
    if (!is.na(j1) & !is.na(j2)) {
      time.since.2009[i, (j1:j2)] <- 1
    }
  }
  # To take off time(s) when camera wasn't working.
  if (dep[i] %in% df_inter_pairs$name_year) {
    df1<- df_inter_pairs[as.character(df_inter_pairs$name_year) == as.character(dep[i]),]
    for (j in seq(1, (nrow(df1) - 1), 2)) { # Assumes all extra times are formatted as END/START pairs
      # Use ceiling() so that day of failure is excluded from operating days
      j1 <- ceiling(julian(df1$date_detected[j], origin = strptime("2009-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")))
      j2 <- floor(julian(df1$date_detected[j + 1], origin = strptime("2009-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%S")))
      if (j2 > j1) time.since.2009[i, j1:(j2-1)] <- 0
    }
  }
}

days.per.year <- c(365,365,365,366,365,365,365,366,365,365,365,366)
Jan1 <- cumsum(c(1, days.per.year))

yrnum <- julday <- NULL

for (i in 1:ncol(time.since.2009)) {
  yrnum[i] <- sum(i >= Jan1)
  julday[i] <- i - Jan1[yrnum[i]] + 1
}
for (i in 1:ncol(time.by.day)) {
  time.by.day[,i] <- rowSums(time.since.2009[,which(julday == i)])
}

rownames(time.by.day) <- rownames(time.since.2009) <- dep
columns <- as.character(1:366)

# Summarise time-by-day for each camera deployment
df_abmi_tbd_summary <- time.by.day %>%
  as_tibble(rownames = "name_year", .name_repair = ~ columns) %>%
  mutate(total_days = rowSums(select(., -name_year)),
         total_summer_days = rowSums(select(., -c(name_year:106, 289:366))),
         total_winter_days = rowSums(select(., -c(name_year, 107:288)))) %>%
  select(name_year, total_days, total_summer_days, total_winter_days)

# Write results
write_csv(df_abmi_tbd_summary, paste0(root, "data/processed/time-by-day/abmi-all-years_tbd-summary_", Sys.Date(), ".csv"))




























