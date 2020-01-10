#-------------------------------------------------------------------------------

# Title: Calculating Mammal Density
# Created: August 23, 2019
# Author: Dave Huggard, Marcus Becker

# Objectives: Calculate density of mammal species for each camera deployment.
#             Density is later used in habitat modeling.

#-------------------------------------------------------------------------------

# Load packages

library(tidyverse)
library(lubridate)
library(mgcv)

# Path to data folder on ABMI Science Centre S: drive
abmisc_data <- "S:/github-repos-data/SC-Camera-Mammals/data/"
# Path to results folder on ABMI Science Centre S: drive
abmisc_results <- "S:/github-repos-data/SC-Camera-Mammals/results/"

#-------------------------------------------------------------------------------

# Import data

# Veg+HF information for each deployment:
df_dep_veghf <- read_csv("./data/lookup/Combined vegHF soilHF and detection distance veg for cameras Feb 2019.csv")

# Deployment-year operating days (total, summer, winter):
df_depy_od <- read_csv("./data/processed/Summary of camera operating days Feb 2019 ALL.csv")

# Deployment lure status:
df_depy_lure <- read_csv("./data/lookup/camera-lure-info-all.csv") %>%
  distinct()

# Native mammals by series:
df_series <- read_csv("./data/processed/Mammals by series Feb 2019 ALL.csv") %>%
  filter(!common_name == "Falcons and allies")

# Species detection distance groups:
df_det_groups <- read_csv("./data/lookup/species-det-distance-groups.csv")

# Look up table for VegHF detection distance categories:
df_pole_veg_lookup <- read_csv("./data/lookup/Lookup for pole distance model predictions.csv")

# Camera study designation:
df_study <- read_csv("./data/lookup/camera-study.csv")

#-------------------------------------------------------------------------------

# Set parameters:

# Summer start/end dates (Julian)
summer.start.j <- 106 # April 16
summer.end.j <- 288 # October 15

# Reconyx camera field-of-view angle
cam_fov_ang <- 42


#-------------------------------------------------------------------------------

# 1. Calculate detection distances - Species groups x Deployment x Season x Veghf type

df_dep_veghf_dist <- df_dep_veghf %>%
  # WetShrub not used in modeling.
  mutate(VegHF = ifelse(VegForDetectionDistance == "WetShrub", "Shrub", VegForDetectionDistance)) %>%
  left_join(df_pole_veg_lookup, by = "VegHF")

test <- df_dep_veghf %>%
  # WetShrub not used in modeling.
  mutate(VegHF = ifelse(VegForDetectionDistance == "WetShrub", "Shrub", VegForDetectionDistance)) %>%
  left_join(df_pole_veg_lookup, by = "VegHF")

SpTable<-c("Bear","BigForestCarnivores","BigForestUngulates","CoyoteFox",
           "Elk (wapiti)","Lynx","Mule deer","SmallForest","SmallOpen","WTDeer",
           "Pronghorn","Bighorn sheep")  # Update this if more modeling groups added.
           # ^ The last (2) are species where the full set of models did not run - treated separately below.

pred.season<-c("Summer","Winter")

site.list<-unique(df_dep_veghf_dist$DeploymentYear)

# Path to detection distance models.
fname.models <- paste0(abmisc_results, "modeling/Detection distance models ")

dd<-dd.lci<-dd.uci<-array(NA,c(length(site.list),length(SpTable),2))

for (sp in c(1:(length(SpTable) - 2))) {  # The last (2) species do not have the full set of models - treated separately.
  fname<-paste(fname.models, SpTable[sp], ".rdata", sep="")
  load(fname)  # Model m, bic.wt
  hab.list<-as.character(m[[2]]$xlevels[[1]])  # The shared subset of habitat types in both VegHF models
  s1<- df_dep_veghf_dist  # Missing data for some veg types for some species, so need to collapse vegHF types in s
  s1$VegHF<-as.character(s1$VegHF)
  s1$VegHF1<-as.character(s1$VegHF1)
  s1$VegHF2<-as.character(s1$VegHF2)
  s1$VegHF3<-as.character(s1$VegHF3)
  s1$VegHF4<-as.character(s1$VegHF4)
  # Substitutions for when a habitat type does not show up in a species' model.
  # Need to revise these by trial and error each year.
  if (("Water" %in% hab.list)==FALSE) s1$VegHF<-ifelse(s1$VegHF=="Water","WetGrass",s1$VegHF)
  if (("WetTreed" %in% hab.list)==FALSE) s1$VegHF<-ifelse(s1$VegHF=="WetTreed","WetGrass",s1$VegHF)
  if (("WetGrass" %in% hab.list)==FALSE) s1$VegHF<-ifelse(s1$VegHF=="WetGrass","Grass",s1$VegHF)
  if (("Grass" %in% hab.list)==FALSE) s1$VegHF<-ifelse(s1$VegHF=="Grass","WetGrass",s1$VegHF)
  if (("Shrub" %in% hab.list)==FALSE) {
    s1$VegHF<-ifelse(s1$VegHF=="Shrub","Grass",s1$VegHF)
    s1$VegHF2<-ifelse(s1$VegHF2=="Shrub","GrassWater",s1$VegHF2)
    s1$VegHF4<-ifelse(s1$VegHF4=="Shrub","GrassWaterHF",s1$VegHF4)
  }
  if ("wet" %in% m[[10]]$xlevels==FALSE) {
    s1$VegHF1<-ifelse(s1$VegHF1=="Wet","GrassShrub",s1$VegHF1)
    s1$VegHF3<-ifelse(s1$VegHF3=="Wet","GrassShrubHF",s1$VegHF3)
  }
  for (j in 1:2) {  # Predictions for the two seasons
    p<-p.se<-array(0,c(length(m),length(site.list)))
    p1<-predict(m[[1]],se.fit=TRUE)
    p[1,]<-rep(p1$fit[1],length(site.list))
    p.se[1,]<-rep(p1$se.fit[1],length(site.list))
    for (k in 2:length(m)) {
      p1<-predict(m[[k]],newdata=data.frame(s1[,c("VegHF","VegHF1","VegHF2","VegHF3","VegHF4","VegHF5")],Season=pred.season[j]),se.fit=TRUE)
      # The SE produces CI's on the detection distance, but not currently used (would be used if we did full CI's on estimates)
      p[k,]<-p[k,]+p1$fit
      p.se[k,]<-p.se[k,]+p1$se.fit
    }
    p.all<-colSums(p*bic.wt)  # BIC weighted average
    p.se.all<-colSums(bic.wt*sqrt(p.se^2+(p-rep(p.all,each=length(m)))^2))
    # Convert probability of being behind 5m pole into effective detection distance
    dd[,sp,j]<-5/sqrt(1-plogis(p.all))
    dd.lci[,sp,j]<-5/sqrt(1-plogis(p.all-2*p.se.all))
    dd.uci[,sp,j]<-5/sqrt(1-plogis(p.all+2*p.se.all))
  }  # Next time period j
}

# Separate predictions for pronghorn and bighorn, because only have null model
for (sp in (length(SpTable)-1):length(SpTable)) {
  fname<-paste(fname.models,SpTable[sp],".rdata",sep="")
  load(fname)  # Model m, bic.wt, hab.list
  s1<-s  # Missing data for some veg types for some species, so need to collapse vegHF types in s
  for (j in 1:2) {
    p1<-predict(m[[1]],se.fit=TRUE)  # Only null model
    dd[,sp,j]<-rep(5/sqrt(1-plogis(p1$fit[1])),length(site.list))  # Only null model
    dd.lci[,sp,j]<-rep(5/sqrt(1-plogis(p1$fit[1]-2*p1$se.fit[1])),length(site.list))  # Only null model
    dd.uci[,sp,j]<-rep(5/sqrt(1-plogis(p1$fit[1]+2*p1$se.fit[1])),length(site.list))  # Only null model
  }  # Next time period j
}  # Next of these two species

dimnames(dd)[[1]]<-dimnames(dd.lci)[[1]]<-dimnames(dd.uci)[[1]]<-site.list
dimnames(dd)[[2]]<-dimnames(dd.lci)[[2]]<-dimnames(dd.uci)[[2]]<-SpTable
dimnames(dd)[[3]]<-dimnames(dd.lci)[[3]]<-dimnames(dd.uci)[[3]]<-pred.season

#-------------------------------------------------------------------------------

# 2. Tidy data

# Calculate total duration spent on camera by deployment, year, and species:
df_series_2 <- df_series %>%
  # Calculte Julian days and season
  mutate(Julian = as.numeric(format(ymd_hms(Time1), "%j")),
         Season = ifelse(Julian >= summer.start.j & Julian <= summer.end.j,
                         "Summer",
                         "Winter")) %>%
  mutate_at(c("DeploymentYear", "common_name", "Season"), factor) %>%
  # Important to add .drop = F argument in group_by so all combinations are preserved.
  group_by(DeploymentYear, common_name, Season, .drop = FALSE) %>%
  summarise(Total_Duration = sum(total_time_tbp)) %>%
  ungroup() %>%
  mutate_if(is.factor, as.character) %>%
  # Join operating days and lure status
  left_join(df_depy_od, by = "DeploymentYear") %>%
  left_join(df_depy_lure, by = "DeploymentYear")
  # ^Note: Missing info from deployments w/o any images of native mammals.

# Character vector of all native species (43 in total)
chr_species <- as.character(sort(unique(df_series_2$common_name)))

# Add deployments w/o any images of native mammals
df_series_3 <- df_depy_od %>%
  # Preserve only those Deployment-Years w/o images of native mammals.
  anti_join(df_series_2, by = "DeploymentYear") %>% # 519 Deployment-Years
  # Expand to include all combinations of DeploymentYear, Season, and Species.
  expand(DeploymentYear, Season = c("Winter", "Summer"), common_name = chr_species) %>%
  # Join in operating days.
  left_join(df_depy_od, by = "DeploymentYear") %>%
  # Total_Duration is 0 because these species weren't seen.
  mutate(Total_Duration = 0) %>%
  # Join in lure status.
  left_join(df_depy_lure, by = "DeploymentYear")

# Bind togther the two df's of series (w/ and w/o native mammals)
df_series_full <- df_series_2 %>%
  bind_rows(df_series_3) %>%
  arrange(DeploymentYear, common_name, Season) %>%
  # Combine total.summer and total.winter into `seasonal` (Season already specified)
  mutate(seasonal = ifelse(Season == "Summer", total.summer, total.winter)) %>%
  select(DeploymentYear:Season, total, seasonal, Lure, Total_Duration)

# Load detection distance modeling from above
load("./beta/det-dist/Detection distances by site species and season.rData")

# Detection distance for each DeploymentYear, SpeciesGroup, and Season combination, tidied.
df_detdist <- dd %>%
  as.data.frame() %>%
  rownames_to_column(var = "DeploymentYear") %>%
  gather(key = "SpGroupSeason", value = "detdist", Bear.Summer:`Bighorn sheep.Winter`) %>%
  mutate(SpGroupSeason = str_replace_all(SpGroupSeason, "[//(//)// ]", "")) %>%
  # Create two new columns: Detection Distance Group and Season, sep by "."
  separate(SpGroupSeason, into = c("DetDistGroup", "Season")) %>%
  mutate(DetDistGroup = str_replace(DetDistGroup, "wapiti", ""))

# Join detection distance information to df_series_full.
# Now contains all the ingredients to calculate density at each deployment.
df_dens_ing <- df_series_full %>%
  left_join(df_det_groups, by = "common_name") %>%
  left_join(df_detdist, by = c("DeploymentYear", "Season", "DetDistGroup")) %>%
  # ^ Missing some detdist info ... which I traced back to missing veghf info.
  # Something to look into.
  select(DeploymentYear, common_name, DetDistGroup, Season, detdist, everything())

#-------------------------------------------------------------------------------

# 3. Calculate Density

df_dens <- df_dens_ing %>%
  # Calculate effort per site, in 100-m^2 * days
  mutate(effort = seasonal * (detdist^2 * pi * (cam_fov_ang/360)) / 100,
  # Calculate seconds (s) of animal presence per effort
         cpue = Total_Duration / effort,
  # Convert to per km^2
         cpue_km2 = cpue / 60 / 60 / 24 * 10000)

write_csv(df_dens, "./data/processed/Camera mammal data processed Feb 2019.csv")

#-------------------------------------------------------------------------------

# 4. Summaries by study - This is done as part of deciding which studies to use in the habitat modeling.

df_dens <- read_csv("./data/processed/Camera mammal data processed Feb 2019.csv")

# Study list - update if needed.
Study.list<-c("ABMI","ABMI-W","OG-ABMI","CMU","BigGrid","OG-CITSCI","OGW-CITSCI",
              "OG-RIVR","OG-EI","OG-AAC","OG-OGC","Nexen","WHEC")


# General summary of effort in each study.
df_summary_effort <- df_dens %>%
  left_join(df_study, by = "DeploymentYear") %>%
  select(DeploymentYear, Season, total, seasonal, Study) %>%
  distinct() %>%
  group_by(Study, Season) %>%
  summarise(total.days = sum(total),
            total.seasonal = sum(seasonal),
            deployments = n_distinct(DeploymentYear)) %>%
  mutate(DaysPerDep = round(total.days / deployments, digits = 2),
         SeasonalDaysPerDep = round(total.seasonal / deployments, digits = 2))
# ^ Note: Not quite clear how it's presented. To fix later.




s1 <- site_metadata

























