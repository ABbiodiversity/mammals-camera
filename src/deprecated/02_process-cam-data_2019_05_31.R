# ------------------------------------------------------------------------------

# Title: Processing Camera Data
# Created: May 30, 2019
# Author: Dave Huggard, Marcus Becker

# Objectives: Process and tidy ABMI camera data for use in downstream analyses

# ------------------------------------------------------------------------------

# Load packages
library(tidyverse)
library(lubridate)

# Import camera data
df_cam_all_18 <- read_csv("S:/github-repos-data/SC-Mammals/data/base/raw/cam-data-abmi_2018.csv")

# Explore data

df_deploy <- df_cam_all_18 %>%
  group_by(deployment_id) %>%
  count() # 820 deployments.

summary(df_cam_all)

#-------------------------------------------------------------------------------

# Step 1 - Reduce initial data file and correct wrong dates

# Scan for wrong dates
df_deptime <- df_cam_all %>%
  mutate(date_time_original = ymd_hms(date_time_original)) %>%
  group_by(deployment_id) %>%
  slice(which.min(date_time_original)) %>%
  select(deployment_id, name, date_time_original) %>%
  arrange(date_time_original)
  # Seems like there are 3 cameras that are "too early" and 1 that is "too late"
  # Doesn't quite jibe with Dave's table in Access - doesn't add 365 days to one
  # of the three camera (ABMI-Cudde-504-NE)
  # To-do: Need to add or subtract days (365) where appropriate.

# Omit time lapse mode images

df_cam_all %>% group_by(trigger_mode) %>% count # n=5 CodeLoc Not Entered

df_cam_md <- df_cam_all %>%
  filter(trigger_mode != "Time Lapse")
  # Haven't added/subtracted days yet.

# Extra Start and End times

df_cam_all %>% group_by(out_of_range) %>% count

df_range <- df_cam_all %>%
  filter(out_of_range == "End" | out_of_range == "Start")

write_csv(df_range, "S:/github-repos-data/SC-Mammals/data/processed/extra-start-end/extra-times_2018.csv")

# All deployments with images

df_cam_img <- df_cam_all %>%
  filter(tagged == TRUE) %>%
  group_by(deployment_id, name, deployment_year) %>%
  summarise(Images = n())

write_csv(df_cam_img, "./data/processed/All deployments with images 2018.csv")

#-------------------------------------------------------------------------------

# Step 2 - Add N gap classes for following NONE

# Import reduced dataset (just motion detected) directly from Access database

df_reduced <- read_csv("C:/Users/mabec/Documents/R/SC-Mammals/data/base/Camtraps_abmi_data_2018_2.csv")
# Note: same obs as df_cam_md from above

df_reduced_gap <- df_cam_md %>%
  mutate(name = as.character(name)) %>%
  mutate(date_time_original1 = ymd_hms(date_time_original)) %>%
  arrange(name, date_time_original1) %>%
  mutate(common_name1 = lead(common_name)) %>%
  select(id:common_name, common_name1, date_time_original, date_time_original1,
         gap_class, tsn_id:distance) %>%
  mutate(gap_class = ifelse(is.na(gap_class) & common_name != "NONE" &
                            common_name1 == "NONE", "N",
                            as.character(gap_class))) %>%
  select(-common_name1)

test$gap_class<-ifelse(test$gap_class=="" & d1$common_name!="NONE" & d1$common_name1=="NONE","N",as.character(d1$gap_class))

#-------------------------------------------------------------------------------

# Step 3 - Produce native-mammals only file and other outputs from Access

# 3.1 - Summarise all species found in current dataset

# Import species groups look-up table
df_species_groups <- read_csv("./data/lookup/camera-species-groups.csv")

# Summarize species categories
df_species <- df_reduced_gap %>%
  group_by(scientific_name, common_name) %>%
  summarise(count = n()) %>%
  filter(common_name %in% df_species_groups$common_name) # All species included

# 3.2 - Filter for just native mammals
df_native <- df_reduced_gap %>%
  left_join(df_species_groups, by = "common_name") %>%
  filter(Type == "Mammal")# Have an extra six observations.

# 3.3 - Multiple animal records with comments to check

# 3.4 - Distance output

# Import camera info look-up table
df_cam_info <- read_csv("./data/lookup/camera-info-2018.csv")

df_cam_info_l <- df_cam_info %>%
  select(deployment, Lure) %>%
  rename(name = deployment)

df_distance <- df_native %>%
  left_join(df_cam_info_l, by = "name") %>%
  filter(Lure == "No") %>%
  select(name, deployment_id, deployment_year, date_time_original1, common_name,
         sex, age_class, number_individuals, distance, Lure)

write_csv(df_distance, "./data/processed/Distance output for native mammals unlured 2018.csv")
# Supposed to append this to results from other years as well.

# 3.5 - Run cow summaries

# 3.6 - Add to camera start/end times data file
# Import start/end times and lure information for each deployment
# Includes both ABMI and non-ABMI cameras
# Straight from Access database? Hmmm.
df_cam_time <- read_csv("./data/base/Start end times Lure Feb 2019 ALL.csv")

#-------------------------------------------------------------------------------

# Step 4 - Reconcile deployment names in different input files

#-------------------------------------------------------------------------------

# Step 5 - Process data with second-pass information

# Import native mammal data from all years
df_native_all <- read_csv("./data/base/native-mammals-2019.csv",
                          col_types = cols(distance = col_character(),
                                          number_during_gap = col_number(),
                                          number_individuals = col_character()),
                          # Important for 'distance' variable
                          na = "")

# Checks for consistency between df_native and df_native_all (2018 data)
test1 <- df_native %>%
  group_by(name) %>%
  summarise(count = n()) %>%
  rename(deployment = name)

test2 <- df_native_all %>%
  filter(Year == "2018") %>%
  group_by(deployment) %>%
  summarise(count = n()) %>%
  left_join(test1, by = "deployment") # Looks consistent

# Summarise species seen
df_species_all <- df_native_all %>%
  group_by(common_name) %>%
  summarise(count = n()) %>%
  arrange(desc(count))
  # Lots of white-tailed deer.
  # Potential issue: 'Deer' category (vs white-tailed or mule specified)

df_native_all_1 <- df_native_all %>%
  mutate(common_name = ifelse(common_name == "Mule Deer", "Mule deer",
                              common_name)) %>%
  mutate(deployment = ifelse(str_detect(deployment, "^WH"),
                             tolower(deployment),
                             deployment),
         deployment = ifelse(str_detect(deployment, "^w"),
                             gsub("^whec", "WHEC", deployment), deployment)) %>%
  # Lose 'VNA' for some elk obs
  mutate(number_individuals = as.numeric(ifelse(is.na(number_individuals) | 0,
                                                1,
                                                number_individuals))) %>%
  # No 'OGW-CITSCI-1538-31' deployment, oddly.
  mutate(distance = str_replace_all(distance, "NA", "X"),
         distance = str_replace_all(distance, "NA, NA", "X, X"),
         distance = str_replace_all(distance, "NA, NA, NA", "X, X, X")) %>%
  # Make columns of number of individuals at each pole position
  mutate(PoleAt = str_count(distance, "A"),
         PoleBehind = str_count(distance, "B"),
         PoleFront = str_count(distance, "F"),
         PoleIC = str_count(distance, "IC"),
         PoleIP = str_count(distance, "IP"),
         PoleNA = str_count(distance, "X"))

# 5.2 - Summarize individuals per image (just for reporting)

# 5.3 - Identify series

df_series <- df_native_all_1 %>%
  mutate(Time1 = ymd_hms(date_time_taken)) %>%
  arrange(deployment, Year, common_name, Time1) %>%
  mutate(SeriesNum = 0,
         GapCheck = "") %>%
  mutate(Time2 = lag(Time1),
         # Calculate difference in time (seconds)
         DiffTime = as.numeric(Time1 - Time2), # Don't need to worry about NAs
         common_name1 = lag(common_name),
         DiffSp = ifelse(common_name != common_name1, TRUE, FALSE),
         deployment1 = lag(deployment),
         DiffSite = ifelse(deployment != deployment1, TRUE, FALSE),
         DiffSeries = ifelse(DiffSite == TRUE | DiffSp == TRUE | DiffTime > 120,
                             1, 0),
         SeriesNum = c(0, cumsum(DiffSeries[-1])), # Off by 1 from Dave's script
         GapCheck = ifelse(DiffSite == FALSE & DiffSp == FALSE &
                             (DiffTime <= 120 & DiffTime >= 20), 1, 0),
         gap_class = lag(gap_class))

# 5.4 - Plot gap length versus left/not



# 5.5 - Gap-leaving probability models for each gap group

# Import gap group data
df_gap <- read_csv("./data/base/Gap groups.csv")

df_gapcheck <- df_series %>%
  filter(!is.na(gap_class) & GapCheck == "1")

p.gap<-array(NA,c(max(df_gap$GapGroup),120))

for (i in 1:max(df_gap$GapGroup)) {
  sp.list<-df_gap$Sp[df_gap$GapGroup==i]  # Species that are in that gap group
  d5<-df_gapcheck[df_gapcheck$common_name %in% sp.list,]  # All records for those species
  p.gap[i,20:120]<- predict(smooth.spline(d5$DiffTime[d5$GapCheck==1],
                          ifelse(d5$gap_class[d5$GapCheck==1]=="P",0,1),df=3),
                          x=20:120)$y
  # Smooth spline used to model probability of leaving.
  # Assume only code "P" used for animals that stayed
}

row.names(p.gap)<-c("Most ungulates","Moose","Small carnivores","Canids cougar",
                    "Bears","Small mammals")

# 5.6 - Re-run series numbering to include checked gap information

# This is for when 20-120 second gaps have been checked from images, so don't
# need probabilistic gap-leaving.

df_series_2 <- df_series %>%
  mutate(gap_class = replace_na(gap_class, "")) %>%
  mutate(DiffSeries1 = ifelse(
    DiffSp == TRUE | DiffSite == TRUE | DiffTime > 120 |
      (gap_class == "L" | gap_class == "N"), 1, 0)) %>%
  mutate(SeriesNum = c(0, cumsum(DiffSeries1[-1]))) %>%
  # Gaps needing probabilistic time assignment
  mutate(ProbGap = ifelse(GapCheck == 1 & (gap_class == "" | gap_class == "U"),
                          1, 0)) %>%
  mutate(ProbGap = replace_na(ProbGap, 0))

# Remove hyphens in CMU deployments

test <- df_series_2 %>%
  filter(str_detect(deployment, "CAL|CHR|DEW|FMM|LLB|LRN|MAC|MCC|WAB")) %>%
  group_by(deployment) %>%
  tally
  # Hyphens don't seem to be present in CMU deployments.

# Save df_series_2
write_csv(df_series_2, "./data/processed/Mammal images with series and checkable gaps Feb 2019 ALL.csv")

# 5.7 - Total time for each series

# Time between photos, by species
df_tbp <- df_series_2 %>%
  mutate(SeriesNum1 = lag(SeriesNum)) %>%
  filter(SeriesNum == SeriesNum1) %>%
  group_by(common_name) %>%
  summarise(time_btwn_photos = mean(DiffTime),
            sample_size = n())

write_csv(df_tbp, "./data/processed/Table of time between photos within series by species.csv")

pgap <- p.gap %>%
  set_rownames(c("Most ungulates","Moose","Small carnivores","Canids cougar",
                 "Bears","Small mammals")) %>%
  as.data.frame() %>%
  rownames_to_column(var = "GapGroup") %>%
  gather(DiffTime, Prob, V1:V120) %>%
  mutate(DiffTime = as.numeric(str_remove(DiffTime, "V")))

df_gap1 <- df_gap %>%
  rename(common_name = Sp) %>%
  mutate(GapGroup = case_when(
    GapGroup == 1 ~ "Most ungulates",
    GapGroup == 2 ~ "Moose",
    GapGroup == 3 ~ "Small carnivores",
    GapGroup == 4 ~ "Canids cougar",
    GapGroup == 5 ~ "Bears",
    GapGroup == 6 ~ "Small mammals"
  )) %>%
  select(-nGap)

df_series_3 <- df_series_2 %>%
  left_join(df_gap1, by = "common_name") %>%
  left_join(pgap, by = c("GapGroup", "DiffTime")) %>%
  mutate(Prob = replace_na(Prob, 1),
         # Adjusted for probability of leaving
         DiffTime_adj = ifelse(ProbGap == 1, DiffTime * (1 - Prob), DiffTime))

# Time for each series
df_tfs_1 <- df_series_3 %>%
  mutate(SeriesNum1 = lag(SeriesNum)) %>%
  filter(SeriesNum1 == SeriesNum) %>%
  group_by(SeriesNum) %>%
  summarise(count = n(),
            total_time = sum(DiffTime_adj))

df_tfs_2 <- df_series_3 %>%
  mutate(SeriesNum1 = lag(SeriesNum)) %>%
  group_by(SeriesNum) %>%
  summarise(count = n()) %>%
  filter(count == 1) %>%
  mutate(total_time = 0)

df_tfs <- df_tfs_1 %>%
  bind_rows(df_tfs_2) %>%
  arrange(SeriesNum)
  # Not a perfect solution, suspect I'll find a better one eventually.

# Add time-between-photos to each series
df_time_final <- df_series_3 %>%
  left_join(df_tbp, by = "common_name") %>%
  select(SeriesNum, common_name, time_btwn_photos) %>%
  distinct() %>%
  left_join(df_tfs, by = "SeriesNum") %>%
  mutate(total_time_tbp = total_time + time_btwn_photos) %>%
  select(SeriesNum, common_name, count, total_time_tbp)

# 5.8 - Mean number of animals per series

df_n_ani <- df_series_3 %>%
  mutate(number_individuals = replace_na(number_individuals, 1)) %>%
  group_by(SeriesNum) %>%
  # Calculate number of photos in each series; mean number of animals
  summarise(n_photos = n(),
            mean_animals = mean(number_individuals))

df_series_4 <- df_series_3 %>%
  group_by(SeriesNum) %>%
  arrange(Time1) %>%
  filter(row_number() == 1) %>%
  left_join(df_n_ani, by = "SeriesNum") %>%
  left_join(df_time_final, by = c("SeriesNum", "common_name")) %>%
  rename(Deployment = deployment) %>%
  # Lure information
  left_join(df_cam_time, by = c("Deployment", "Year")) %>%
  dplyr::select(Deployment, Year, Lure, Time1, SeriesNum, n_photos, common_name,
                sex, age_class, multiple_animals, mean_animals, total_time_tbp)

# 5.9 - Summarise time that the cameras are operating

# Deployments with images - processed from other script (for now)

df_dwti <- read_csv("./data/processed/All deployments with images July 2019 ALL.csv") %>%
  mutate(DeploymentYear = paste(Deployment, Year, sep = "_"))

# All deployments start and end times, regardless of whether there are images
df_cam_time <- read_csv("./data/base/Start end times Lure Feb 2019 ALL.csv") # n = 5324

df_cam_time <- df_cam_time %>%
  mutate(StartTime = ymd_hms(StartTime),
         EndTime = ymd_hms(EndTime),
         DeploymentYear = paste(Deployment, Year, sep = "_")) %>%
  # Semi join to keep all rows from cam_time with matching values in dwti
  semi_join(df_dwti, by = "DeploymentYear")

sitelist <- unique(as.character(df_cam_time$DeploymentYear))

# Read in intermediate end start pairs csv file (from Dave)
df_inter <- read_csv("./data/base/Intermediate END START pairs Feb 2019 ALL.csv")

df_inter <- df_inter %>%
  mutate(DeploymentYear = paste(Deployment, Year, sep = "_")) %>%
  mutate(StartEndTime1 = ymd_hms(date_time_original)) %>%
  semi_join(df_dwti, by = "DeploymentYear")
  # Looks like we're missing WHEC cameras from 2010. Weird, because how do we know about the END START pairs then

# To track each julian day that cameras were operating
time.by.day<-array(0,c(length(sitelist),366))

# Instead of adding 365 each year, we can define a variable at the start of the script.
time.since.2009<-array(0,c(length(sitelist),365+365+365+366+365+365+365+366+365+365+365))

# Now, for some for loop magic.

for (i in 1:length(sitelist)) {
  df_cam_time_1 <- df_cam_time[df_cam_time$DeploymentYear == sitelist[i],]
  for (k in 1:nrow(df_cam_time_1)) {
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
  if (sitelist[i] %in% df_inter$DeploymentYear) {
    df_inter_1 <- df_inter[as.character(df_inter$DeploymentYear) == as.character(sitelist[i]),]
    for (j in seq(1, (nrow(df_inter_1) - 1), 2)) {
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
  write_csv("./data/processed/Camera operating days Feb 2019 ALL.csv")

time.since.2009 %>% as_tibble(rownames = "DeploymentYear") %>%
  write_csv("./data/processed/Camera operating days since 2009 Feb 2019 ALL.csv")

# Remove mammal records after final end or in intermediate end-start periods

df_series_4 <- df_series_4 %>%
  mutate(DeploymentYear = paste(Deployment, Year, sep = "_"))

x <- rep(0, nrow(df_series_4))

for (i in 1:nrow(df_series_4)) {
  j <- floor(julian(df_series_4$Time1[i], origin = strptime("2009-01-01 00:00:00",format="%Y-%m-%d %H:%M:%S")))
  k <- match(df_series_4$DeploymentYear[i], rownames(time.since.2009))
  if (is.na(k) | is.na(j) | time.since.2009[k,j] == 0) x[i] <- 1
}

df_series_5 <- df_series_4[x==0,] # Starting to deviate from Dave a bit here.

# Export data

write_csv(df_series_5, "./data/processed/Mammals by series Feb 2019 ALL.csv")

# Summary table of series (including off-grid, other projects) - for reporting.

df_series_2 <- read_csv("./data/processed/Mammal images with series and checkable gaps Feb 2019 ALL.csv")

df_series_2_img <- df_series_2 %>%
  group_by(common_name, Year) %>%
  count(name = "Images")

df_summary <- df_series_5 %>%
  filter(Year >= "2015") %>%
  group_by(common_name, Year) %>%
  summarise(Series = n(),
            Sites = n_distinct(Deployment),
            Total_Indiv = sum(mean_animals),
            Total_Dur = sum(total_time_tbp),
            Mean_Dur = Total_Dur / Series) %>%
  left_join(df_series_2_img, by = c("common_name", "Year"))

write_csv(df_summary, "./data/processed/Series summary table.csv")

# Summarize total time, winter, and summer time

time.by.day <- read_csv("./data/processed/Camera operating days Feb 2019 ALL.csv")

df_tbd_summary <- time.by.day %>%
  mutate(total = rowSums(select(., -DeploymentYear)),
         total.summer = rowSums(select(., -c(DeploymentYear:V106, V289:V366))),
         total.winter = rowSums(select(., -c(DeploymentYear, V107:V288)))) %>%
  select(DeploymentYear, total, total.summer, total.winter)

write_csv(df_tbd_summary, "./data/processed/Summary of camera operating days Feb 2019 ALL.csv")

# Plot deployment time figures - for reporting

x<-as.character(time.by.day$DeploymentYear)
year<-as.numeric(substr(x,nchar(x)-3,nchar(x)))
ordernum<-ifelse(substr(x,1,4)=="ABMI",ifelse(substr(x,1,6)=="ABMI-W",2,1),3)  # Put ABMI first, then ABMI-W, then OG and other studies
time.by.day1<-time.by.day[order(year,ordernum,x),]  # Sort by year, grouping, alphabetical deployment name
x<-as.character(time.by.day1$DeploymentYear)  # And re-do these variables in the new order
year<-as.numeric(substr(x,nchar(x)-3,nchar(x)))

for (yr in 1:4) {
  fname<-paste("C:/Users/mabec/Documents/R/SC-Mammals/beta/Deployment Operating Days ",yr+2014,".png",sep="")
  png(fname,width=500,height=800)
  t<-time.by.day1[year==yr+2014,]
  x<-as.character(t$DeploymentYear)  # And re-do these variables in the new order
  col1<-ifelse(substr(x,1,4)=="ABMI",ifelse(substr(x,1,6)=="ABMI-W","blue2","green2"),"orange2")
  y<-1:nrow(t)
  plot(0,0,xlim=c(1,366),ylim=c(1,max(y)),xlab="",xaxt="n",typ="n",yaxt="n")
  axis(side=1,at=cumsum(c(0,31+28,31+30,31+30,31+31,30+31)),lab=c("Jan","Mar","May","Jul","Sep","Nov"),tcl=0.015,cex.axis=1.2)
  for (i in 1:length(y)) points(2:367,rep(max(y)-y[i]+1,366),pch=15,cex=sign(as.numeric(t[i,2:367]))*0.3,col=col1[i])  # ABMI at top, to others at bottom
  abline(v=106,lty=2)  # April 16
  abline(v=288,lty=2)  # Oct 15
  text(53,-20,"Winter",adj=0.5)
  text(328,-20,"Winter",adj=0.5)
  text(198,-20,"Summer",adj=0.5)
  mtext(side=3,at=1,yr+2014,cex=1.3,adj=0)
  graphics.off()
}

# Histograms of operating hours -> maybe this would be a good place for walk().
par(mfrow=c(2,2))
x<-hist(c(q$time.total),plot=FALSE)
hist(q$time.total[q$Year==2015],xlab="Operating time (days)",main="2015",col="grey79",cex.lab=1.3,breaks=x$breaks)
hist(q$time.total[q$Year==2016],xlab="Operating time (days)",main="2016",col="grey79",cex.lab=1.3,breaks=x$breaks)
hist(q$time.total[q$Year==2017],xlab="Operating time (days)",main="2017",col="grey79",cex.lab=1.3,breaks=x$breaks)
hist(q$time.total[q$Year==2018],xlab="Operating time (days)",main="2018",col="grey79",cex.lab=1.3,breaks=x$breaks)
par(mfrow=c(2,2))
x<-hist(c(q$time.winter),plot=FALSE)
hist(q$time.winter[q$Year==2015],xlab="Operating time - Winter",main="2015",col="grey79",cex.lab=1.3,breaks=x$breaks)
hist(q$time.winter[q$Year==2016],xlab="Operating time - Winter",main="2016",col="grey79",cex.lab=1.3,breaks=x$breaks)
hist(q$time.winter[q$Year==2017],xlab="Operating time - Winter",main="2017",col="grey79",cex.lab=1.3,breaks=x$breaks)
hist(q$time.winter[q$Year==2018],xlab="Operating time - Winter",main="2018",col="grey79",cex.lab=1.3,breaks=x$breaks)
par(mfrow=c(2,2))
x<-hist(c(q$time.summer),plot=FALSE)
hist(q$time.summer[q$Year==2015],xlab="Operating time - Summer",main="2015",col="grey79",cex.lab=1.3,breaks=x$breaks)
hist(q$time.summer[q$Year==2016],xlab="Operating time - Summer",main="2016",col="grey79",cex.lab=1.3,breaks=x$breaks)
hist(q$time.summer[q$Year==2017],xlab="Operating time - Summer",main="2017",col="grey79",cex.lab=1.3,breaks=x$breaks)
hist(q$time.summer[q$Year==2018],xlab="Operating time - Summer",main="2018",col="grey79",cex.lab=1.3,breaks=x$breaks)

# 5.10 - Summarise distance wrt pole into detection distance models

# No NAs, blanks.

# Lure information for all deployments.
df_cam_lure <- read_csv("./data/lookup/camera-lure-info-all.csv")

df_distance_all <- df_native_all_1 %>%
  select(deployment:common_name, sex, age_class, number_individuals, distance:PoleNA) %>%
  filter(!is.na(distance)) %>%
  left_join(df_cam_lure, by = c("deployment", "Year")) %>%
  mutate(Julian = as.numeric(format(ymd_hms(date_time_taken), "%j"))) %>%
  filter(Lure == "n")
  # Note that this is different from Dave's output. Not sure why. Seems like he's missing some observations.

# Species distance groups
df_dist_groups <- read_csv("./data/lookup/species-distance-groups.csv")

# Add veg+HF information for camera point

df_veghf_dist <- read_csv("./data/lookup/Combined vegHF soilHF and detection distance veg for cameras Feb 2019.csv") %>%
  select(deployment, Year, DeploymentYear, VegForDetectionDistance) %>%
  mutate(VegHF = ifelse(VegForDetectionDistance == "WetShrub", "Shrub", VegForDetectionDistance)) %>%
  select(-VegForDetectionDistance) %>%
  right_join(df_distance_all, by = c("deployment", "Year")) %>%
  mutate(Time1 = ymd_hms(date_time_taken)) %>%
  left_join(df_dist_groups, by = "common_name")

# Make combined veg types -> seems like different categories/classifications to try in the models.

df_veghf_dist$VegHF1<-"Wet"
df_veghf_dist$VegHF1<-ifelse(df_veghf_dist$VegHF=="Conif" | df_veghf_dist$VegHF=="Decid", "ConifDecid", df_veghf_dist$VegHF1)
df_veghf_dist$VegHF1<-ifelse(df_veghf_dist$VegHF=="Grass" | df_veghf_dist$VegHF=="Shrub", "GrassShrub", df_veghf_dist$VegHF1)
df_veghf_dist$VegHF1<-ifelse(df_veghf_dist$VegHF=="HF","HF",df_veghf_dist$VegHF1)

df_veghf_dist$VegHF2<-"GrassWater"
df_veghf_dist$VegHF2<-ifelse(df_veghf_dist$VegHF=="Conif" | df_veghf_dist$VegHF=="Decid" | df_veghf_dist$VegHF=="WetTreed","Treed",df_veghf_dist$VegHF2)
df_veghf_dist$VegHF2<-ifelse(df_veghf_dist$VegHF=="Shrub","Shrub",df_veghf_dist$VegHF2)
df_veghf_dist$VegHF2<-ifelse(df_veghf_dist$VegHF=="HF","HF",df_veghf_dist$VegHF2)

df_veghf_dist$VegHF3<-"Wet"
df_veghf_dist$VegHF3<-ifelse(df_veghf_dist$VegHF=="Conif" | df_veghf_dist$VegHF=="Decid","ConifDecid",df_veghf_dist$VegHF3)
df_veghf_dist$VegHF3<-ifelse(df_veghf_dist$VegHF=="Grass" | df_veghf_dist$VegHF=="Shrub" | df_veghf_dist$VegHF=="HF","GrassShrubHF",df_veghf_dist$VegHF3)

df_veghf_dist$VegHF4<-"GrassWaterHF"
df_veghf_dist$VegHF4<-ifelse(df_veghf_dist$VegHF=="Conif" | df_veghf_dist$VegHF=="Decid" | df_veghf_dist$VegHF=="WetTreed","Treed",df_veghf_dist$VegHF4)
df_veghf_dist$VegHF4<-ifelse(df_veghf_dist$VegHF=="Shrub","Shrub",df_veghf_dist$VegHF4)

df_veghf_dist$VegHF5<-"GrassWaterShrubHF"
df_veghf_dist$VegHF5<-ifelse(df_veghf_dist$VegHF=="Conif" | df_veghf_dist$VegHF=="Decid" | df_veghf_dist$VegHF=="WetTreed","Treed",df_veghf_dist$VegHF5)

# Create lookup table for pole distance model predictions.
df_pole_veg_lookup <- df_veghf_dist %>%
  select(VegHF, VegHF1:VegHF5) %>%
  filter_all(all_vars(!is.na(.))) %>%
  distinct()

write_csv(df_pole_veg_lookup, "./data/lookup/Lookup for pole distance model predictions.csv")

# Estimate time between images behind and in front of the pole, by species

# Something to come back to one day - just auxiliary info for now.

d1<-df_veghf_dist[order(df_veghf_dist$DeploymentYear,df_veghf_dist$common_name,df_veghf_dist$Time1),]
d1.1<-d1[2:nrow(d1),]  # Offset one row
d1<-d1[-nrow(d1),]
time.diff<-difftime(d1.1$Time1,d1$Time1,units="secs")
time.diff.front<-time.diff.behind<-time.diff.front.n<-time.diff.behind.n<-NULL  # One value for each SpGroup
for (sp in 1:length(SpTable)) {

  # Time differences for consecutive records of a species in that group where both records were in front
  time.diff.front1<-time.diff[d1.1$SpGroup==SpTable[sp] & d1.1$DeploymentYear==d1$DeploymentYear & d1.1$common_name==d1$common_name & time.diff<120 & time.diff>0 & d1.1$PoleFront==d1.1$number_individuals & d1$PoleFront==d1$number_individuals]

  # Time differences for consecutive records of a species in that group where both records were behind
  time.diff.behind1<-time.diff[d1.1$SpGroup==SpTable[sp] & d1.1$DeploymentYear==d1$DeploymentYear & d1.1$common_name==d1$common_name & time.diff<120 & time.diff>0 & d1.1$PoleBehind==d1.1$number_individuals & d1$PoleBehind==d1$number_individuals]

  time.diff.front[sp]<-mean(time.diff.front1,na.rm=T)
  time.diff.behind[sp]<-mean(time.diff.behind1,na.rm=T)

  time.diff.front.n[sp]<-length(time.diff.front1)
  time.diff.behind.n[sp]<-length(time.diff.behind1)
}

q<-data.frame(SpGroup=SpTable,TimeDiffFront=time.diff.front,TimeDiffBehind=time.diff.behind,TimeDiffFront.n<-time.diff.front.n,TimeDiffBehind.n=time.diff.behind.n)

# 5.11 - Detection Distance Models

library(mgcv)

df_pole_veg_lookup <- read_csv("./data/lookup/Lookup for pole distance model predictions.csv")

# Only use records where number_individuals = PoleBehind or PoleFront

df_detection <- df_veghf_dist %>%
  mutate(n = rowSums(select(., PoleAt:PoleIP))) %>%
  filter(!is.na(VegHF)) %>%
  filter(n == PoleBehind | n == PoleFront) %>%
  mutate(Season = as.factor(ifelse(Julian >= summer.start.j & Julian <= summer.end.j,
                         "Summer",
                         "Winter")))

SpTable <- unique(df_detection$SpGroup)

for (sp in 1:length(SpTable)) {

  print(paste(sp,length(SpTable),SpTable[sp],date()))

  d.sp<-df_detection[df_detection$SpGroup==SpTable[sp],]

  # Number of informative trials for that record CHANGE THIS TO SIGN - Doesn't matter how many animals - only one trips the camera.
  d.sp$n<-d.sp$PoleFront+d.sp$PoleBehind

  d.sp$pBehind<-d.sp$PoleBehind/d.sp$n

  m<-list(NULL)
  m[[1]]<-try(gam(pBehind~1,weights=d.sp$n,data=d.sp,family="binomial"))
  if (SpTable[sp]!="Pronghorn" & SpTable[sp]!="Bighorn sheep") {
    m[[2]]<-try(gam(pBehind~as.factor(VegHF),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[3]]<-try(gam(pBehind~as.factor(VegHF1),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[4]]<-try(gam(pBehind~as.factor(VegHF2),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[5]]<-try(gam(pBehind~as.factor(VegHF3),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[6]]<-try(gam(pBehind~as.factor(VegHF4),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[7]]<-try(gam(pBehind~as.factor(VegHF5),weights=d.sp$n,data=d.sp,family="binomial"))

    m[[8]]<-try(gam(pBehind~Season,weights=d.sp$n,data=d.sp,family="binomial"))

    m[[9]]<-try(gam(pBehind~as.factor(VegHF)+as.factor(Season),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[10]]<-try(gam(pBehind~as.factor(VegHF1)+as.factor(Season),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[11]]<-try(gam(pBehind~as.factor(VegHF2)+as.factor(Season),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[12]]<-try(gam(pBehind~as.factor(VegHF3)+as.factor(Season),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[13]]<-try(gam(pBehind~as.factor(VegHF4)+as.factor(Season),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[14]]<-try(gam(pBehind~as.factor(VegHF5)+as.factor(Season),weights=d.sp$n,data=d.sp,family="binomial"))
  }

  nModels<-length(m)

  # Additional modification to not count models where the minimum number of pole positions in a veg-HF type is <20
  min.n<-c(min(by(d.sp$n,d.sp$VegHF,sum)),min(by(d.sp$n,d.sp$VegHF1,sum)),min(by(d.sp$n,d.sp$VegHF2,sum)),min(by(d.sp$n,d.sp$VegHF3,sum)),min(by(d.sp$n,d.sp$VegHF4,sum)),min(by(d.sp$n,d.sp$VegHF5,sum)))
  min.n<-ifelse(is.na(min.n),0,min.n)

  # BIC calculation
  bic.ta<-rep(999999999,(nModels))
  for (i in 1:(nModels)) {
    if (!is.null(m[[i]]) & class(m[[i]])[1]!="try-error") {  # last part is to not used non-converged models, unless none converged
      bic.ta[i]<-BIC(m[[i]])
    }
  }

  if (SpTable[sp]!="Pronghorn" & SpTable[sp]!="Bighorn sheep") {
    bic.ta[c(2:7)]<-bic.ta[c(2:7)]+(min.n<20)*999999999  # Inflate BIC if <20 records in any one vegHF type
    bic.ta[c(9:14)]<-bic.ta[c(9:14)]+(min.n<20)*999999999  # Inflate BIC if <20 records in any one vegHF type
  }

  if (SpTable[sp]=="Bear" | SpTable[sp]=="Bighorn sheep" | SpTable[sp]=="Elk (wapiti)") bic.ta[8:14]<-999999999  # Too few winter data to fit properly
  bic.delta<-bic.ta-min(bic.ta)
  bic.exp<-exp(-1/2*bic.delta)
  bic.wt<-bic.exp/sum(bic.exp)
  best.model<-which.max(bic.wt)

  # Figures predicting for each veg type and season
  hab.list<-levels(as.factor(d.sp$VegHF))
  v.l1<-df_pole_veg_lookup[df_pole_veg_lookup$VegHF %in% hab.list,]
  if (sum(d.sp$VegHF=="Water")==0) v.l1[7,1]<-"WetGrass"  # For species that have no data in water
  p.veg<-array(0,c(length(hab.list),nModels))  # veg+HF types, for each model
  p.season<-array(0,c(2,nModels))  # Two seasons
  p.veg[,1]<-p.season[,1]<-mean(predict(m[[1]]))  # Null model

  if (SpTable[sp]!="Pronghorn" & SpTable[sp]!="Bighorn sheep") {
    for (i in 2:nModels) {
      if(class(m[[i]])[1]!="try-error" & bic.wt[i]>0.001 ) {
        p.veg[,i]<-predict(m[[i]],newdata=data.frame(v.l1,Season="Summer"))
        p.season[,i]<-predict(m[[i]],newdata=data.frame(v.l1[1,],Season=c("Summer","Winter")))
      }
    }
    p1.veg<-colSums(bic.wt*t(p.veg))
    p1.season<-colSums(bic.wt*t(p.season))
  } else {
    p1.veg<-t(p.veg)
    p1.season<-p.season
  }

  dist.veg<-5/sqrt(1-plogis(p1.veg)) # Time-between-images effect would be added in here, but no meaningful differences found
  dist.season<-5/sqrt(1-plogis(p1.season))  # Time-between-images effect would be added in here, but no meaningful differences found
  dist.veg<-ifelse(dist.veg>20,20,dist.veg)  # Can be inf when all locations are behind pole
  dist.season<-ifelse(dist.season>20,20,dist.season)  # Can be inf when all locations are behind pole
  ymax<-ifelse(max(c(dist.veg,dist.season))>10,20,10)
  hab.n<-by(d.sp$n,d.sp$VegHF,sum)
  fname<-paste("C:/Users/mabec/Documents/R/SC-Mammals/beta/Detection distance ",SpTable[sp],".jpg",sep="")
  jpeg(width=500,height=900,file=fname)
  par(mfrow=c(2,1),mai=c(1.8,0.8,0.3,0.3))
  # Figures - veg types
  x1<-barplot(dist.veg,ylim=c(0,ymax),xlab="",ylab="Effective distance (m)",cex.axis=1.3,cex.lab=1.4)
  mtext(side=1,at=x1,hab.list,las=2,line=1,cex=1.3)
  box(bty="l")
  text(x1[1],ymax*0.96,SpTable[sp],cex=1.5,adj=0)
  text(x1,rep(ymax*0.89,length(x1)),hab.n[match(v.l1$VegHF,names(hab.n))],cex=1.1)
  # Figures - season
  par(mai=c(1,1.8,0.1,1.3))
  x1<-barplot(dist.season,ylim=c(0,ymax),xlab="",ylab="Effective distance (m)",cex.axis=1.3,cex.lab=1.4)
  mtext(side=1,at=x1,c("Summer","Winter"),las=1,line=1,cex=1.3)
  box(bty="l")
  graphics.off()

  # Save models and bic wt
  fname<-paste("C:/Users/mabec/Documents/R/SC-Mammals/beta/Detection distance models ",SpTable[sp],".rdata",sep="")
  save(file=fname,m,bic.wt,hab.list)

}  # Next SpGroup

# 5.12 Lure Effects

df_series_lure <- df_series_5 %>%
  ungroup() %>%
  # Don't use off-grids, wetlands, or non-ABMI projects
  filter(str_detect(Deployment, "^ABMI"),
         !str_detect(Deployment, "^ABMI-W"),
         !is.na(Lure),
         # Only from 2015 onwards
         Year >= 2015)

df_cam_time_1 <- df_cam_time %>%
  ungroup() %>%
  filter(str_detect(Deployment, "^ABMI"),
         !str_detect(Deployment, "^ABMI-W"),
         Year >= 2015) %>%
  # to add qualifying deployments without native mammals
  anti_join(df_series_lure, by = "DeploymentYear") %>%
  select(Deployment, Year, Lure) %>%
  distinct() %>%
  group_by(Year, Lure) %>%
  tally() %>%
  ungroup()

df_lure_summary <- df_series_lure %>%
  ungroup() %>%
  select(Deployment, Year, Lure) %>%
  distinct() %>%
  group_by(Year, Lure) %>%
  tally() %>%
  ungroup() %>%
  left_join(df_cam_time_1, by = c("Year", "Lure")) %>%
  mutate(total = n.x + n.y) %>%
  select(Year, Lure, total) %>%
  arrange(Lure, Year)

df_series_lure <- df_series_lure %>%
  left_join(df_lure_summary, by = c("Year", "Lure"))

df_lure_effects <- df_series_lure %>%
  group_by(common_name, Year, Lure) %>%
  summarise(sum_dur = sum(total_time_tbp),
            length_dur = length(total_time_tbp)) %>%
  left_join(df_lure_summary, by = c("Year", "Lure")) %>%
  mutate(mean_dur = sum_dur / total,
         mean_rec = length_dur / total,
         mean_dur_per_rec = mean_dur / mean_rec)

write_csv(df_lure_effects, "./data/processed/Series summary unlured versus lured Feb 2019 ALL.csv")

# Bootstrapped version to check for changes in lure effect in different years

q<-by(df_series_lure$n_photos,df_series_lure$common_name,sum)  # Figure out top species
sp.top<-names(sort(-q))[1:23]
sp.top<-sp.top[-which(sp.top=="Bighorn sheep")]  # Only one year
sp.top<-sp.top[-which(sp.top=="Richardson's Ground Squirrel")]
sp.top<-sp.top[-which(sp.top=="Bison")]  # Only one year
d2<-df_series_lure[df_series_lure$common_name %in% sp.top,]  # Limit to top species only
d2$Deployment<-as.character(d2$Deployment)
d2$common_name <-as.factor(as.character(d2$common_name))  # And make sure that excluded species are excluded
d2<-d2[substr(d2$Deployment,1,4)=="ABMI" & d2$Year>=2015,]  # Only use ABMI on-grids with paired design
d2<-d2[substr(d2$Deployment,1,6)!="ABMI-W",]  # And don't use wetlands for same reason
d2<-d2[!is.na(d2$Lure),]
d2$Deployment<-as.character(d2$Deployment)  # This ensures that excluded deployments are not included
d2$Site<-substr(d2$Deployment,6,nchar(d2$Deployment)-3)  # This assumes all sites are ABMI-[Site]-{NW,NE,SW,SE}
# Add qualifying deployments with no native mammals (in sp.top)
t2<-df_cam_time_1[is.na(match(df_cam_time_1$DeploymentYear,d2$DeploymentYear)),]   # t1 prepared above
t3<-data.frame(Deployment=t2$Deployment,Year=t2$Year,Lure=t2$Lure,Time1=t2$StartTime,SeriesNum=0,n_photos=0,common_name=NA,sex=NA,age_class=NA,multiple_animals=NA,mean_animals=0,total_time_tbp=0,DeploymentYear=t2$DeploymentYear,total=NA,Site=NA)
d2<-rbind(d2,t3)
niter<-1000
lure.bs<-array(NA,c(length(sp.top),4,2,niter))  # For each species, year, percent of duration or nPhotos in lured sites, for each iteration
for (yr in 1:4) {
  d1<-d2[d2$Year==2014+yr,]
  site.list<-unique(d1$Site)
  for (iter in 1:niter) {
    if((iter-1)/100==floor((iter-1)/100)) print(paste(yr,iter,date()))
    s<-sample(1:length(site.list),length(site.list),replace=TRUE)
    i<-unlist(lapply(site.list[s],function(a) which(d1$Site %in% a)))
    d1.bs<-d1[i,]
    q.u<-by(d1.bs$total_time_tbp[d1.bs$Lure=="n"],d1.bs$common_name[d1.bs$Lure=="n"],sum)/length(unique(d1.bs$Deployment[d1.bs$Lure=="n"]))  # Average duration per site
    q.l<-by(d1.bs$total_time_tbp[d1.bs$Lure=="y"],d1.bs$common_name[d1.bs$Lure=="y"],sum)/length(unique(d1.bs$Deployment[d1.bs$Lure=="y"]))
    n.u<-by(d1.bs$total_time_tbp[d1.bs$Lure=="n"],d1.bs$common_name[d1.bs$Lure=="n"],length)/length(unique(d1.bs$Deployment[d1.bs$Lure=="n"]))  # Average records per site
    n.l<-by(d1.bs$total_time_tbp[d1.bs$Lure=="y"],d1.bs$common_name[d1.bs$Lure=="y"],length)/length(unique(d1.bs$Deployment[d1.bs$Lure=="y"]))
    lure.bs[,yr,1,iter]<-q.l/(q.u+q.l)*100
    lure.bs[,yr,2,iter]<-n.l/(n.u+n.l)*100
  }
}
lure.bs.sum<-array(NA,c(length(sp.top),4,2,7))
for (sp in 1:length(sp.top)) {
  for (yr in 1:4) {
    for (x in 1:2) {
      lure.bs.sum[sp,yr,x,]<-quantile(lure.bs[sp,yr,x,],c(0.025,0.05,0.1,0.5,0.9,0.95,0.975),na.rm=TRUE)
    }
  }
}
dimnames(lure.bs)[[1]]<-dimnames(lure.bs.sum)[[1]]<-sort(sp.top)
dimnames(lure.bs)[[2]]<-dimnames(lure.bs.sum)[[2]]<-2015:2018
dimnames(lure.bs)[[3]]<-dimnames(lure.bs.sum)[[3]]<-c("duration","photos")

# Figure of results

sp.order<-c(1,8,7,3,14,19,6,9,2,10,5,11,17,20,13,16,18,12,15)  # Bears, dogs, cats+weasels, ungulates, others
x<-c(1,2,3.5,4.5,5.5,7,8,9,10.5,12,13,14,15,16,17,18.5,19.5,20.5)
col1<-rep(rainbow(10,v=0.8),2)
png(file="C:/Users/mabec/Documents/R/SC-Mammals/beta/Lure effect over years.png",height=600,width=900)
par(mai=c(3,0.9,0.3,0.3))
plot(0,0,xlim=range(x),ylim=log(c(0.2,21)),xlab="",ylab="Lured:Unlured",xaxt="n",yaxt="n",typ="n")
axis(side=2,at=log(c(0.1,0.2,0.25,1/3,1/2,1,2,3,4,5,7,10,20)),lab=rep("",13),tck=1,col="grey80")
axis(side=2,at=log(c(0.1,0.2,0.25,1/3,1/2,1,2,3,4,5,7,10,20)),lab=c(0.1,0.2,0.25,0.33,0.5,1,2,3,4,5,7,10,20),tck=0.015,cex.lab=1.3)
for (i in 1:length(x)) {
  y<-(lure.bs.sum[sp.order[i],,1,]/100)/(1-lure.bs.sum[sp.order[i],,1,]/100)  # Convert percent of time in lured to lured/unlured
  points(x[i]+c(-0.3,-0.1,0.1,0.3),log(y[,4]),pch=18,cex=1.7,typ="p",lwd=2,col=col1[i])
  for (j in 1:4) lines(rep(x[i]+c(-0.3,-0.1,0.1,0.3)[j],2),log(y[j,c(2,6)]),col=col1[i])
  mtext(side=1,at=x[i],dimnames(lure.bs.sum)[[1]][sp.order[i]],las=2,cex=1.3,adj=1,line=0.5,col=col1[i])
}
graphics.off()

# 5.13 Camera year calibration

# Wasn't done for 2018.




























