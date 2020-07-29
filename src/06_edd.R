#-----------------------------------------------------------------------------------------------------------------------

# Title: Effective detection distance (EDD) modeling
# Author: Dave Huggard, Marcus Becker

# Previous scripts: 01_clean-raw-data

#-----------------------------------------------------------------------------------------------------------------------

# Attach packages
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(mgcv)

# Set path to Google Drive
root <- "G:/Shared drives/ABMI Camera Mammals/"

# Lure lookup
df_lure <- read_csv(paste0(root, "data/lookup/lure/all-cam-lure_2020-06-01.csv"))

# Veg/HF lookup
df_veghf <- read_csv(paste0(root,"data/lookup/Combined vegHF soilHF and detection distance veg for cameras May 2020.csv"))

# EDD modeling species groups
df_edd_groups <- read_csv(paste0(root, "data/lookup/species-distance-groups.csv"))

# Veg/HF groupings to try in modeling
df_veghf_groups <- read_csv(paste0(root,"data/lookup/veg-pole-distance.csv"))

#-----------------------------------------------------------------------------------------------------------------------

# Retrieve pole position information from previous data (years 2013 through 2018)

# Note: these tags were added with the previous ABMI tagging system

df_pole <- read_csv(
  paste0(root,"data/base/raw/previous/ALL_native-mammals_2019-09-01.csv"),
         col_types = cols(distance = col_character(),
                          number_during_gap = col_number(),
                          number_individuals = col_character()),
         na = "") %>%
  # Only observations with pole information
  filter(!is.na(distance)) %>%
  # Make columns of number of individuals at each pole position
  mutate(at_pole = str_count(distance, "A"),
         behind_pole = str_count(distance, "B"),
         front_pole = str_count(distance, "F"),
         ic_pole = str_count(distance, "IC"),
         ip_pole = str_count(distance, "IP"),
         na_pole = str_count(distance, "NA")) %>%
  # Standardize variable names
  select(name = deployment, year = Year, date_detected = date_time_taken, common_name,
         number_individuals, distance:na_pole) %>%
  # Join lure information
  left_join(df_lure, by = c("name", "year")) %>%
  # Filter out lured deployments; only unlured are used in the EDD modeling
  filter(lure == "No") %>%
  # Create variable for julian date
  mutate(julian = as.numeric(format(ymd_hms(date_detected), "%j")))

# Retrieve Veg and HF information for each deployment
df_pole_veghf <- df_veghf %>%
  select(name = deployment, year = Year, VegForDetectionDistance) %>%
  # Get rid of WetShrub - not used for modeling.
  mutate(VegHF = ifelse(VegForDetectionDistance == "WetShrub",
                        "Shrub", VegForDetectionDistance)) %>%
  select(-VegForDetectionDistance) %>%
  # Join back to pole information
  right_join(df_pole, by = c("name", "year")) %>%
  mutate(date_detected = ymd_hms(date_detected)) %>%
  # Join in EDD species grouping lookup
  left_join(df_edd_groups, by = "common_name") %>%
  # Join alternative VegHF categories to try in the modeling
  left_join(df_veghf_groups, by = "VegHF") %>%
  filter(!is.na(VegHF)) %>%
  # Sum total individuals in each of at_pole through ip_pole
  mutate(n = rowSums(select(., at_pole:ip_pole))) %>%
  # Only use records where number_individuals = behind_pole or front_pole
  filter(n == behind_pole | n == front_pole) %>%
  # Create new season variable based on julian day, and calculate the proportion of individuals behind pole
  mutate(season = as.factor(ifelse(julian  >= summer.start.j & julian <= summer.end.j, "summer", "winter")),
         prop_behind = behind_pole / n) %>%
  select(name, year, common_name, dist_group, n, prop_behind, VegHF, VegHF1:VegHF5, season)

# Save data

#-----------------------------------------------------------------------------------------------------------------------

# Step 2. EDD modeling.

# Detach a couple packages - they wreak some havoc by using tibbles.
unloadNamespace(c("tidyr", "dplyr"))

# Species groups
SpTable <- unique(df_pole_veghf$dist_group)

# Loop through each species group, try each model variation, make predictions, output figures
for (sp in 1:length(SpTable)) {

  print(paste(sp, length(SpTable), SpTable[sp], date()))

  d.sp <- df_pole_veghf[df_pole_veghf$dist_group == SpTable[sp], ]

  m <- list(NULL)
  m[[1]]<-try(gam(prop_behind~1,weights=d.sp$n,data=d.sp,family="binomial"))
  if (SpTable[sp]!="Pronghorn" & SpTable[sp]!="Bighorn sheep") {
    m[[2]]<-try(gam(prop_behind~as.factor(VegHF),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[3]]<-try(gam(prop_behind~as.factor(VegHF1),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[4]]<-try(gam(prop_behind~as.factor(VegHF2),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[5]]<-try(gam(prop_behind~as.factor(VegHF3),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[6]]<-try(gam(prop_behind~as.factor(VegHF4),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[7]]<-try(gam(prop_behind~as.factor(VegHF5),weights=d.sp$n,data=d.sp,family="binomial"))

    m[[8]]<-try(gam(prop_behind~season,weights=d.sp$n,data=d.sp,family="binomial"))

    m[[9]]<-try(gam(prop_behind~as.factor(VegHF)+as.factor(season),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[10]]<-try(gam(prop_behind~as.factor(VegHF1)+as.factor(season),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[11]]<-try(gam(prop_behind~as.factor(VegHF2)+as.factor(season),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[12]]<-try(gam(prop_behind~as.factor(VegHF3)+as.factor(season),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[13]]<-try(gam(prop_behind~as.factor(VegHF4)+as.factor(season),weights=d.sp$n,data=d.sp,family="binomial"))
    m[[14]]<-try(gam(prop_behind~as.factor(VegHF5)+as.factor(season),weights=d.sp$n,data=d.sp,family="binomial"))
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
  hab.list <- levels(as.factor(d.sp$VegHF))
  v.l1<-df_veghf_groups[df_veghf_groups$VegHF %in% hab.list,]
  if (sum(d.sp$VegHF=="Water")==0) v.l1[7,1]<-"WetGrass"  # For species that have no data in water
  p.veg<-array(0,c(length(hab.list),nModels))  # veg+HF types, for each model
  p.season<-array(0,c(2,nModels))  # Two seasons
  p.veg[,1]<-p.season[,1]<-mean(predict(m[[1]]))  # Null model

  if (SpTable[sp]!="Pronghorn" & SpTable[sp]!="Bighorn sheep") {
    for (i in 2:nModels) {
      if(class(m[[i]])[1]!="try-error" & bic.wt[i]>0.001 ) {
        p.veg[,i]<-predict(m[[i]],newdata=data.frame(v.l1, season="summer"))
        p.season[,i]<-predict(m[[i]],newdata=data.frame(v.l1[1,], season=c("summer","winter")))
      }
    }
    p1.veg<-colSums(bic.wt*t(p.veg))
    p1.season<-colSums(bic.wt*t(p.season))
  } else {
    p1.veg<-t(p.veg)
    p1.season<-p.season
  }
  dist.veg<-5/sqrt(1-plogis(p1.veg))  # Time-between-images effect would be added in here, but no meaningful differences found
  dist.season<-5/sqrt(1-plogis(p1.season))  # Time-between-images effect would be added in here, but no meaningful differences found
  dist.veg<-ifelse(dist.veg>20,20,dist.veg)  # Can be inf when all locations are behind pole
  dist.season<-ifelse(dist.season>20,20,dist.season)  # Can be inf when all locations are behind pole
  ymax<-ifelse(max(c(dist.veg,dist.season))>10,20,10)
  hab.n<-by(d.sp$n,d.sp$VegHF,sum)
  fname<-paste0(root, "data/processed/detection-distance/figures/Detection distance ",SpTable[sp],".jpg",sep="")
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
  fname<-paste0(root, "data/processed/detection-distance/group-models/Detection distance models ",SpTable[sp],".rdata",sep="")
  save(file=fname,m,bic.wt,hab.list)

}  # Next SpGroup

