# Calculate detection distance for each species group and days operating for each camera site,
# and use these with the image data summarized to series to calculate density

# Detection distance models are created in a separate script,
# as part of general second-pass camera data processing

library(mgcv)  # For GAM predictions for detection distances

# Calculate detection distance for each species group in each season

# Need veg+HF information for camera point
# s<-read.csv("C:/Dave/ABMI/Data/Site info/2018/Combined vegHF soilHF and detection distance veg for cameras Feb 2019.csv")  # Processed to broad veg + HF categories in Excel
s <- read.csv("./data/lookup/Combined vegHF soilHF and detection distance veg for cameras Feb 2019.csv")

s$VegHF<-NULL  # Replace full veg+HF with broad classes for detection distance
names(s)[which(names(s)=="VegForDetectionDistance")]<-"VegHF"
s$VegHF<-as.character(s$VegHF)
s$VegHF<-ifelse(s$VegHF=="WetShrub","Shrub",s$VegHF)  # WetShrub not used in modeling
# Make combined veg types
s$VegHF1<-"Wet"
s$VegHF1<-ifelse(s$VegHF=="Conif" | s$VegHF=="Decid","ConifDecid",s$VegHF1)
s$VegHF1<-ifelse(s$VegHF=="Grass" | s$VegHF=="Shrub","GrassShrub",s$VegHF1)
s$VegHF1<-ifelse(s$VegHF=="HF","HF",s$VegHF1)
s$VegHF2<-"GrassWater"
s$VegHF2<-ifelse(s$VegHF=="Conif" | s$VegHF=="Decid" | s$VegHF=="WetTreed","Treed",s$VegHF2)
s$VegHF2<-ifelse(s$VegHF=="Shrub","Shrub",s$VegHF2)
s$VegHF2<-ifelse(s$VegHF=="HF","HF",s$VegHF2)
s$VegHF3<-"Wet"
s$VegHF3<-ifelse(s$VegHF=="Conif" | s$VegHF=="Decid","ConifDecid",s$VegHF3)
s$VegHF3<-ifelse(s$VegHF=="Grass" | s$VegHF=="Shrub" | s$VegHF=="HF","GrassShrubHF",s$VegHF3)
s$VegHF4<-"GrassWaterHF"
s$VegHF4<-ifelse(s$VegHF=="Conif" | s$VegHF=="Decid" | s$VegHF=="WetTreed","Treed",s$VegHF4)
s$VegHF4<-ifelse(s$VegHF=="Shrub","Shrub",s$VegHF4)
s$VegHF5<-"GrassWaterShrubHF"
s$VegHF5<-ifelse(s$VegHF=="Conif" | s$VegHF=="Decid" | s$VegHF=="WetTreed","Treed",s$VegHF5)
# Use that vegHF info to predict detection distances for each camera, species group and season
SpTable<-c("Bear","BigForestCarnivores","BigForestUngulates","CoyoteFox","Elk (wapiti)","Lynx","Mule deer","SmallForest","SmallOpen","WTDeer","Pronghorn","Bighorn sheep")  # Update this if more modeling groups added.  The last (2) are species where the full set of models did not run - treated separately below
pred.season<-c("Summer","Winter")
site.list<-unique(s$DeploymentYear)

#fname.models<-"C:/Dave/ABMI/Cameras/2018 analysis/Distance models/Detection distance models "
fname.models <- "./beta/Detection distance models "

# Detection distance by site, species group and season
dd<-dd.lci<-dd.uci<-array(NA,c(length(site.list),length(SpTable),2))

for (sp in c(1:(length(SpTable)-2))) {  # The last (2) species do not have the full set of models - treated separately.
  fname<-paste(fname.models,SpTable[sp],".rdata",sep="")
  load(fname)

  # Model m, bic.wt
  hab.list<-as.character(m[[2]]$xlevels[[1]])  # The shared subset of habitat types in both VegHF models
  s1<-s  # Missing data for some veg types for some species, so need to collapse vegHF types in s
  s1$VegHF<-as.character(s1$VegHF)
  s1$VegHF1<-as.character(s1$VegHF1)
  s1$VegHF2<-as.character(s1$VegHF2)
  s1$VegHF3<-as.character(s1$VegHF3)
  s1$VegHF4<-as.character(s1$VegHF4)
  if (("Water" %in% hab.list)==FALSE) s1$VegHF<-ifelse(s1$VegHF=="Water","WetGrass",s1$VegHF)  # Substitutions for when a habitat type does not show up in a species' model.  Need to revise these by trial and error each year
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

      p1<-predict(m[[k]],
                  newdata=data.frame(s1[,c("VegHF","VegHF1","VegHF2","VegHF3","VegHF4","VegHF5")],
                                     Season=pred.season[j]),se.fit=TRUE)

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

#save(dd,dd.lci,dd.uci,SpTable,site.list,file="C:/Dave/ABMI/Cameras/2018 analysis/Detection distances by site species and season.rData")
save(dd, dd.lci, dd.uci, SpTable, site.list, file = "./beta/det-dist/Detection distances by site species and season.rData")


# Plot percentage of cameras operating at each day -
# originally used to find a good common period, but that approach is not being used any more (seasons instead)

tt<-read.csv("C:/Dave/ABMI/Cameras/2018 analysis/Start end times/Camera operating days Feb 2019 ALL.csv")
names(tt)[1]<-"DeploymentYear"  # The other 365 columns are whether or not that site was surveyed on that day
cpd<-colSums(tt[,2:366])/nrow(tt)*100
dplot(1:365,cpd,lwd=3,xlab="",ylab="Cameras operating (%)",xaxt="n",typ="l",ylim=c(0,100))
ml<-c(31,28,31,30,31,30,31,31,30,31,30,31)
x<-c(0,cumsum(ml[1:11]))+1
axis(side=1,at=x,lab=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),tck=0.015)
start1<-min(which(cpd>70))  # Period when >70% operating
end1<-max(which(cpd>70))
lines(rep(start1,2),c(-10,cpd[start1]),lty=2)
lines(rep(end1,2),c(-10,cpd[start1]),lty=2)
abline(70,0,lty=2)
text(start1,cpd[start1]+5,format(as.Date(start1,origin=as.Date("2016-12-31")),"%b %d"),cex=1.2)
text(end1,cpd[end1]+5,format(as.Date(end1,origin=as.Date("2016-12-31")),"%b %d"),cex=1.2)
print(paste("% of camera-days in operating period",round(sum(tt[,(start1+1):(end1+1)])/sum(tt[,2:366])*100,2)))

#-------------------------------------------------------------------------------

# Then process data files together with effort - this is done for all the mammal data

# First, summarize mammal series data file to one row per deployment, with columns for total duration of each species*season

# d<-read.csv("C:/Dave/ABMI/Cameras/2018 analysis/Mammals by series Feb 2019 ALL.csv")
d <- read_csv("./data/processed/Mammals by series Feb 2019 ALL.csv")

d<-d[d$common_name!="Falcons and allies",]  # Don't know how that one record got in here!?

d$Time1<-strptime(d$Time1,format="%Y-%m-%d %H:%M:%S")

SpAll<-as.character(sort(unique(d$common_name)))

d$DeploymentYear<-paste(d$Deployment,d$Year,sep="_")

DepYearList<-as.character(sort(unique(d$DeploymentYear)))

d.d<-array(0,c(length(DepYearList),length(SpAll)*2))  # To store total durations for each deployment_year * species (for two seasons)

rownames(d.d)<-DepYearList
colnames(d.d)<-paste(rep(SpAll,2),rep(c("Summer","Winter"),each=length(SpAll)),sep=".")

d$common_name<-as.character(d$common_name)  # Otherwise, puts results below in the wrong places if a species or deployment is missing from d1

d$DeploymentYear<-as.character(d$DeploymentYear)

j.summer.start<-106  # April 16
j.summer.end<-288  # Oct 15

for (i in 1:nrow(d)) {  # Run through each series and add to appropriate Deployment_Year and Species.Season
  jul<-as.numeric(strftime(d$Time1[i],format="%j"))
  seas<-ifelse(jul>=j.summer.start & jul<=j.summer.end,"Summer","Winter")
  label<-paste(d$common_name[i],seas,sep=".")
  d.d[d$DeploymentYear[i],label]<-d.d[d$DeploymentYear[i],label]+d$total_time_tbp[i]
}

# And compile deployment, year and lure information for each row
i<-match(DepYearList,d$DeploymentYear)  # The first occurrence of each DeploymentYear in d
Dep1<-d$Deployment[i]
Year1<-d$Year[i]
Lured1<-d$Lure[i]
qq<-data.frame(DeploymentYear=DepYearList,Deployment=Dep1,Year=Year1,Lured=Lured1,d.d)
names(qq)<-gsub("\\.","",names(qq))
names(qq)<-gsub("wapiti","",names(qq))

write.table(qq,"./beta/Dave-Objects/Camera data by site duration Feb 2019.csv",sep=",",row.names=FALSE)

#-------------------------------------------------------------------------------

# Then calculate effort for each deployment and season
d.d<-read.csv("./beta/Dave-Objects/Camera data by site duration Feb 2019.csv")  # This just converts into proper dataframe, gets rid of rownames

#t<-read.csv("C:/Dave/ABMI/Cameras/2018 analysis/Start end times/Camera operating days Feb 2019 ALL.csv")
t <- read.csv("./data/processed/Camera operating days Feb 2019 ALL.csv")

names(t)[which(names(t)=="X")]<-"DeploymentYear"

t$TotalDays<-rowSums(t[,2:367])  # Originally substracted 1 day for the half days at start and end, but took that out because it causes infinities for cameras that run <1 day but get images

t$SummerDays<-rowSums(t[,(j.summer.start+1):(j.summer.end+1)])  # Days operating during summer (assumes one initial column only, for DeploymentYear)

t$SummerDays<-ifelse(t$SummerDays<0,0,t$SummerDays)

t$WinterDays<-t$TotalDays-t$SummerDays

d<-merge(d.d,t[,c("DeploymentYear","TotalDays","SummerDays","WinterDays")],all.y=T)
# Note: includes cameras with no photos  Should check for loss of data

# Add deployment, year, lured for deployments that don't have any mammals, and convert NA's to 0's
d$DeploymentYear<-as.character(d$DeploymentYear)
d$Deployment<-substr(d$DeploymentYear,1,nchar(d$DeploymentYear)-5)
d$Year<-substr(d$DeploymentYear,nchar(d$DeploymentYear)-3,nchar(d$DeploymentYear))

# qqq<-read.csv("C:/Dave/ABMI/Cameras/2018 analysis/Start end times/Start end times Lure Feb 2019 ALL.csv")
qqq <- read.csv("./data/lookup/camera-startend-lure-all.csv")

qqq$DeploymentYear<-paste(qqq$Deployment,qqq$Year,sep="_")
lured1<-as.character(qqq$Lure[match(d$DeploymentYear,qqq$DeploymentYear)])
d$Lured<-as.character(d$Lured)
d$Lured<-ifelse(is.na(d$Lured),lured1,d$Lured)
d$Lured<-ifelse(is.na(d$Lured),"n",d$Lured)  # A few remaining unaccounted for deployments, mostly off-grids, so assume no lure
FirstSpCol<-which(names(d)=="BadgerSummer")  # Check if species list changes
LastSpCol<-which(names(d)=="WoodlandCaribouWinter")  # Check if species list changes
for (i in FirstSpCol:LastSpCol) 	d[,i]<-ifelse(is.na(d[,i]),0,d[,i])

# Add detection distances by spgroup and time period
# load("C:/Dave/ABMI/Cameras/2018 analysis/Detection distances by site species and season.rData")  # dd,dd.lci,dd.uci,SpTable,site.list.  Includes sites with metadata but that weren't collected or haven't been tagged (or not properly)
load("./beta/det-dist/Detection distances by site species and season.rData")

q<-data.frame(DeploymentYear=dimnames(dd)[[1]],dd[,,1],dd.lci[,,1],dd.uci[,,1],dd[,,2],dd.lci[,,2],dd.uci[,,2])

names(q)<-gsub("\\.","",names(q))
names(q)<-gsub("wapiti","",names(q))
names(q)[2:(dim(dd)[2]+1)]<-paste(names(q)[2:(dim(dd)[2]+1)],"Summer",sep=".")
names(q)<-gsub("1",".Summer.LCI",names(q))
names(q)<-gsub("2",".Summer.UCI",names(q))
names(q)<-gsub("3",".Winter",names(q))
names(q)<-gsub("4",".Winter.LCI",names(q))
names(q)<-gsub("5",".Winter.UCI",names(q))
names(q)[2:ncol(q)]<-paste("DetDist",names(q)[2:ncol(q)],sep=".")
d<-merge(d,q,by="DeploymentYear")  # Should check for loss of data

#taxa.file<-read.csv("C:/Dave/ABMI/Data/Mammals/2015/Mammal taxa.csv")  # To get proper species names
taxa.file <- read.csv("./data/lookup/Mammal taxa.csv")

#-------------------------------------------------------------------------------

# Process species duration and effort variables into CPUE (density) value
FirstSpCol<-which(names(d)=="BadgerSummer")  # Check if species list changes
LastSpCol<-which(names(d)=="WoodlandCaribouSummer")  # Check if species list changes ("Summer" here, because just want one copy of each species' name for list of all species)
SpAll<-names(d)[FirstSpCol:LastSpCol]
SpAll<-gsub("Summer","",SpAll)

for (sp in 1:length(SpAll)) {
  SpGroup<-taxa.file$DetDistGroup[match(SpAll[sp],taxa.file$Label)]
  # Summer
  i<-which(regexpr("DetDist",names(d))>-1 & regexpr(SpGroup,names(d))>-1 & (regexpr(".LCI",names(d))==-1 & regexpr(".UCI",names(d))==-1 & regexpr("Summer",names(d))>-1))
  detdist<-d[,i]
  effort<-d$SummerDays*detdist^2*3.14159/8/100  # Effort per site, in 100m2*days.
  label<-paste(SpAll[sp],"Summer",sep="")
  CPUE<-d[,label]/effort  # Units here are seconds (of animal presence) per 100m2*day.
  CPUE<-CPUE/60/60/24*10000  # Convert to /km2
  names(d)[which(names(d)==label)]<-paste(label,"Duration",sep=".")  # Rename these, to use the simple species label for the CPUE value
  d[,label]<-CPUE  # And create this variable as the main value for the species



  # Winter
  i<-which(regexpr("DetDist",names(d))>-1 & regexpr(SpGroup,names(d))>-1 & (regexpr(".LCI",names(d))==-1 & regexpr(".UCI",names(d))==-1 & regexpr("Winter",names(d))>-1))
  detdist<-d[,i]
  effort<-d$WinterDays*detdist^2*3.14159/8/100  # Effort per site, in 100m2*days.
  label<-paste(SpAll[sp],"Winter",sep="")
  CPUE<-d[,label]/effort  # Units here are seconds (of animal presence) per 100m2*day.
  CPUE<-CPUE/60/60/24*10000  # Convert to /km2
  names(d)[which(names(d)==label)]<-paste(label,"Duration",sep=".")  # Rename these, to use the simple species label for the CPUE value
  d[,label]<-CPUE  # And create this variable as the main value for the species
}


write.table(d,file="./beta/Dave-Objects/Camera mammal data processed Feb 2019.csv",sep=",",row.names=F)

#-------------------------------------------------------------------------------

# Summaries by study. This is done as part of deciding which studies to use in the habitat modeling
d<-read.csv("./beta/Dave-Objects/Camera mammal data processed Feb 2019.csv")

Study.list<-c("ABMI","ABMI-W","OG-ABMI","CMU","BigGrid","OG-CITSCI","OGW-CITSCI","OG-RIVR","OG-EI","OG-AAC","OG-OGC","Nexen","WHEC")  # Update this and lines below if more studies are added, or site names change
d$Study<-"ABMI"
d$Study<-ifelse(substr(d$Deployment,1,6)=="ABMI-W","ABMI-W",d$Study)
d$Study<-ifelse(substr(d$Deployment,1,4)=="WHEC","WHEC",d$Study)
d$Study<-ifelse(substr(d$Deployment,1,3) %in% c("CAL","CHR","DEW","FMM","LLB","LRN","MAC","MCC","WAB","LIB") == TRUE,"CMU",d$Study)
d$Study<-ifelse(substr(d$Deployment,1,2)=="BG","BigGrid",d$Study)
d$Study<-ifelse(substr(d$Deployment,1,3)=="NEX","Nexen",d$Study)
d$Study<-ifelse(substr(d$Deployment,1,6)=="OG-AAC","OG-AAC",d$Study)
d$Study<-ifelse(substr(d$Deployment,1,6)=="OG-ABM","OG-ABMI",d$Study)
d$Study<-ifelse(substr(d$Deployment,1,9)=="OG-CITSCI","OG-CITSCI",d$Study)
d$Study<-ifelse(substr(d$Deployment,1,9)=="OGW-CITSC","OGW-CITSCI",d$Study)
d$Study<-ifelse(substr(d$Deployment,1,5)=="OG-EI","OG-EI",d$Study)
d$Study<-ifelse(substr(d$Deployment,1,6)=="OG-OGC","OG-OGC",d$Study)
d$Study<-ifelse(substr(d$Deployment,1,7)=="OG-RIVR","OG-RIVR",d$Study)
d$Study<-as.factor(d$Study)  # To ensure 0's show up

# General summary of effort in each study
q<-data.frame(Study=sort(unique(d$Study)),Deployments=as.numeric(by(d$DeploymentYear,d$Study,length)), DeploymentsSummer=as.numeric(by(d$DeploymentYear[d$SummerDays>0],d$Study[d$SummerDays>0],length)), DeploymentsWinter=as.numeric(by(d$DeploymentYear[d$WinterDays>0],d$Study[d$WinterDays>0],length)),
              TotalDays=as.numeric(by(d$TotalDays,d$Study,function(x) sum(x,na.rm=TRUE))), SummerDays=as.numeric(by(d$SummerDays,d$Study,function(x) sum(x,na.rm=TRUE))), WinterDays=as.numeric(by(d$WinterDays,d$Study,function(x) sum(x,na.rm=TRUE))))
q$TotalDaysPerDep<-q$TotalDays/q$Deployments
q$SummerDaysPerDep<-q$SummerDays/q$Deployments
q$WinterDaysPerDep<-q$WinterDays/q$Deployments
q<-q[match(Study.list,q$Study),]
write.table(q,file="./beta/Dave-Objects/Study summary - effort.csv",sep=",",row.names=FALSE)

# Summary of main species in each study - summer
SpTable<-c("WhitetailedDeer","Muledeer","Moose","Elk","WoodlandCaribou","Pronghorn","BlackBear","Coyote","GrayWolf","CanadaLynx","Fisher","Marten","SnowshoeHare","WhitetailedJackRabbit")
study.list<-sort(unique(d$Study))
d1<-d[,c("Study",paste(SpTable,"Summer",sep=""))]
q.summer<-data.frame(Study=sort(unique(d$Study)))
dens.summer.sum<-array(NA,c(length(study.list),length(SpTable),3))  # Mean, 5%CI, 95%CI for each study and species
for (sp in 1:length(SpTable)) {
  for (s in 1:length(study.list)) {
    x<-d1[d1$Study==study.list[s],paste(SpTable[sp],"Summer",sep="")]
    x<-x[!is.na(x)]  # NA's where that season not sampled at that deployment
    p<-sum(sign(x))/length(x)  # Proportion occurrences
    x1<-x[x>0]  # Abundance given presence
    if (length(x1)>1) {
      x.bs<-rbinom(10000,length(x),p)/length(x) * exp(rnorm(10000,mean(log(x1)),sd(log(x1))/sqrt(length(x1))))
      x.bs<-x.bs*mean(x)/mean(x.bs)  # Adjust for any log bias
      dens.summer.sum[s,sp,]<-c(mean(x),quantile(x.bs,c(0.05,0.95)))
    }
    else {
      dens.summer.sum[s,sp,]<-c(mean(x),NA,NA)
    }
  }
}

dimnames(dens.summer.sum)<-list(study.list,SpTable,c("Mean","LCI","UCI"))

dens.summer.sum<-dens.summer.sum[match(Study.list,dimnames(dens.summer.sum)[[1]]),,]

write.table(dens.summer.sum,file="./beta/Dave-Objects/Study summary - summer densities.csv",sep=",",col.names=NA)

# Summary of main species in each study - winter
SpTable<-c("WhitetailedDeer","Muledeer","Moose","Elk","WoodlandCaribou","Pronghorn","BlackBear","Coyote","GrayWolf","CanadaLynx","Fisher","Marten","SnowshoeHare","WhitetailedJackRabbit")
study.list<-sort(unique(d$Study))
d1<-d[,c("Study",paste(SpTable,"Winter",sep=""))]
q.winter<-data.frame(Study=sort(unique(d$Study)))
dens.winter.sum<-array(NA,c(length(study.list),length(SpTable),3))  # Mean, 5%CI, 95%CI for each study and species
for (sp in 1:length(SpTable)) {
  for (s in 1:length(study.list)) {
    x<-d1[d1$Study==study.list[s],paste(SpTable[sp],"Winter",sep="")]
    x<-x[!is.na(x)]  # NA's where that season not sampled at that deployment
    p<-sum(sign(x))/length(x)  # Proportion occurrences
    x1<-x[x>0]  # Abundance given presence
    if (length(x1)>1) {
      x.bs<-rbinom(10000,length(x),p)/length(x) * exp(rnorm(10000,mean(log(x1)),sd(log(x1))/sqrt(length(x1))))
      x.bs<-x.bs*mean(x)/mean(x.bs)  # Adjust for any log bias
      dens.winter.sum[s,sp,]<-c(mean(x),quantile(x.bs,c(0.05,0.95)))
    }
    else {
      dens.winter.sum[s,sp,]<-c(mean(x),NA,NA)
    }
  }
}

dimnames(dens.winter.sum)<-list(study.list,SpTable,c("Mean","LCI","UCI"))

dens.winter.sum<-dens.winter.sum[match(Study.list,dimnames(dens.winter.sum)[[1]]),,]

write.table(dens.winter.sum,file="C:/Dave/ABMI/Cameras/2018 analysis/Summaries/Study summary - winter densities.csv",sep=",",col.names=NA)

# Summarize locations, incl. maps
s<-read.csv("./data/lookup/Deployment locations all Feb 2019.csv")  # Summarized from larger meta-data file

s$Public.Latitude<-as.numeric(as.character(s$Public.Latitude))
s$Public.Longitude<-as.numeric(as.character(s$Public.Longitude))
s$NearestSite<-as.numeric(as.character(s$NearestSite))

# To get natural regions via nearest ABMI sites
s1<-read.csv("C:/Dave/ABMI/Data/Site info/Site summary with climate.csv")

# To get natural regions via nearest ABMI sites
for (i in 1:nrow(s)) {  # Fill in missing nearest ABMI sites based on lat longs
  if (is.na(s$NearestSite[i]) & !is.na(s$Lat[i])) s$NearestSite[i]<-s1$SITE_ID[which.min(sqrt((s$Lat[i]-s1$PUBLIC_LATTITUDE)^2 + (s$Long[i]-s1$PUBLIC_LONGITUDE)^2))]
}

s$NR<-s1$NATURAL_REGIONS[match(s$NearestSite,s1$SITE_ID)]

s$DeploymentYear<-paste(s$Site.Name,s$Year,sep="_")

names(s)[which(names(s)=="Public.Latitude")]<-"Lat"
names(s)[which(names(s)=="Public.Longitude")]<-"Long"

d$DeploymentYear[is.na(match(d$DeploymentYear,s$DeploymentYear))]# Check that there aren't many of these

d<-merge(d,s[,c("DeploymentYear","Lat","Long","NR")])

# Plot each study location
load("./data/lookup/kgrid_table_km.Rdata")  # For map background using NR's

c1<-c(rgb(0.7,0.9,0.9),rgb(0.8,1,0.9),rgb(0.8,0.9,0.7),rgb(1,0.8,1),rgb(1,0.9,0.7),rgb(1,1,0.9))[match(kgrid$NRNAME,c("Canadian Shield","Boreal","Foothills","Rocky Mountain","Parkland","Grassland"))]
for (i in 1:length(Study.list)) {
  d1<-d[d$Study==Study.list[i],]
  fname<-paste("./beta/study-maps/Study maps ",Study.list[i],".png",sep="")
  png(file=fname,width=600,height=800)
  plot(kgrid$POINT_X,kgrid$POINT_Y,pch=15,cex=1,col=c1,xlim=range(d$Long,na.rm=T),ylim=range(d$Lat,na.rm=T))
  points(d1$Long,d1$Lat,cex=1,lwd=2)
  mtext(side=3,at=-120,adj=0,Study.list[i],cex=1.3)
  graphics.off()
}

# Extra cluster diagrams of studies based on densities.  Not saved, just cut and paste images if needed
plot(hclust(dist(dens.summer.sum[,,1],method="canberra")))
xxx<-dens.winter.sum[,,1]  # Winter takes a few more steps, because not all studies had winter samples
xxx<-xxx[!is.na(xxx[,1]),]
xxx<-xxx[rowSums(xxx)>0,]
plot(hclust(dist(xxx,method="canberra")))


# Summaries of species' densities by season
d<-read.csv("./beta/Dave-Objects/Camera mammal data processed Feb 2019.csv")

d<-d[substr(d$Deployment,1,4)=="ABMI" & substr(d$Deployment,1,6)!="ABMI-W",] # ABMI on-grid terrestrial sites only

SpTable<-c("WhitetailedDeer","Muledeer","Elk","Moose","WoodlandCaribou","GrayWolf","Coyote","CanadaLynx","Fisher","Marten","SnowshoeHare","WhitetailedJackRabbit")

dens.season.sum<-array(NA,c(length(SpTable),2,3))  # For each species, summer/winter, {mean, lci,uci}
dimnames(dens.season.sum)[[1]]<-SpTable
dimnames(dens.season.sum)[[2]]<-c("Summer","Winter")
dimnames(dens.season.sum)[[3]]<-c("Mean","LCI","UCI")
for (sp in 1:length(SpTable)) {
  y.s<-d[,paste(SpTable[sp],"Summer",sep="")]
  y.w<-d[,paste(SpTable[sp],"Winter",sep="")]
  y.s<-na.omit(y.s)
  y.w<-na.omit(y.w)
  # Put confidence intervals on each season
  p.s<-sum(sign(y.s))/length(y.s)  # Proportion occurrences
  agp.s<-y.s[y.s>0]  # Abundance given presence
  # Resampling CI's for two-part binomial presence/absence and log-normal abundance-given-presence.
  # Doesn't work if <2 occurrences.
  y.s.bs<-rbinom(10000,length(y.s),p.s)/length(y.s) * exp(rnorm(10000,mean(log(agp.s)),sd(log(agp.s))/sqrt(length(agp.s))))
  y.s.bs<-y.s.bs*mean(y.s)/mean(y.s.bs)  # Adjust for any log bias
  dens.season.sum[sp,1,]<-c(mean(y.s),quantile(y.s.bs,c(0.05,0.95)))
  p.w<-sum(sign(y.w))/length(y.w)  # Proportion occurrences
  agp.w<-y.w[y.w>0]  # Abundance given presence
  y.w.bs<-rbinom(10000,length(y.w),p.w)/length(y.w) * exp(rnorm(10000,mean(log(agp.w)),sd(log(agp.w))/sqrt(length(agp.w))))
  y.w.bs<-y.w.bs*mean(y.w)/mean(y.w.bs)  # Adjust for any log bias
  dens.season.sum[sp,2,]<-c(mean(y.w),quantile(y.w.bs,c(0.05,0.95)))
}

x<-c(4,3,2,1,5,7,8,9,10,11,13,14)
fname<-"./beta/seasonal-densities/Seasonal density.png"
png(file=fname,width=500,height=600)
par(mai=c(2.1,0.9,0.3,0.3))
plot(x,dens.season.sum[,1,1],xlab="",ylab="Density (/km2)",pch=18,cex=2,col="red3",xaxt="n")
axis(side=1,at=x,lab=rep("",length(x)),tcl=0.015)
points(x+0.3,dens.season.sum[,2,1],pch=18,cex=2,col="blue2")
for (i in 1:length(x)) {
  lines(rep(x[i],2),dens.season.sum[i,1,2:3],col="red4")
  lines(rep(x[i],2)+0.3,dens.season.sum[i,2,2:3],col="blue2")
}
mtext(side=1,at=x,SpTable,line=0.5,las=2,cex=1.2)
graphics.off()


# Summaries by year
d<-read.csv("./beta/Dave-Objects/Camera mammal data processed Feb 2019.csv")

d<-d[substr(d$Deployment,1,4)=="ABMI" & substr(d$Deployment,1,6)!="ABMI-W",] # ABMI on-grid terrestrial sites only

SpTable<-c("WhitetailedDeer","Muledeer","Elk","Moose","WoodlandCaribou","GrayWolf","Coyote","CanadaLynx","Fisher","Marten","SnowshoeHare","WhitetailedJackRabbit")

dens.year.sum<-array(NA,c(length(SpTable),4,3))  # For each species, 2015:2018, {mean, lci,uci}
dimnames(dens.year.sum)[[1]]<-SpTable
dimnames(dens.year.sum)[[2]]<-2015:2018
dimnames(dens.year.sum)[[3]]<-c("Mean","LCI","UCI")
for (sp in 1:length(SpTable)) {
  for (yr in 1:4) {
    d1<-d[d$Year==yr+2014,]
    d1<-d1[d1$SummerDays>=10 & d1$WinterDays>=10,] # Only use deployments with enough summer and winter densities
    y.s<-d1[,paste(SpTable[sp],"Summer",sep="")]
    y.w<-d1[,paste(SpTable[sp],"Winter",sep="")]
    y<-(y.s+y.w)/2  # Equally weighted seasons
    # Put confidence intervals on
    p<-sum(sign(y))/length(y)  # Proportion occurrences
    agp<-y[y>0]  # Abundance given presence
    y.bs<-rbinom(10000,length(y),p)/length(y) * exp(rnorm(10000,mean(log(agp)),sd(log(agp))/sqrt(length(agp))))
    y.bs<-y.bs*mean(y)/mean(y.bs)  # Adjust for any log bias
    dens.year.sum[sp,yr,]<-c(mean(y),quantile(y.bs,c(0.05,0.95)))
  }  # Next yr
}  # Next species

fname<-"./beta/seasonal-densities/Yearly density.png"
png(file=fname,width=700,height=800)
par(mfrow=c(4,3),mai=c(0.8,0.9,0.3,0.3))
for (sp in 1:length(SpTable)) {
  plot(2015:2018,dens.year.sum[sp,,1],xlab="",ylab="Density (/km2)",pch=18,cex=2,col="blue2",xaxt="n",ylim=c(0,max(dens.year.sum[sp,,])),xlim=c(2014.8,2018) )
  axis(side=1,at=2015:2018,cex.axis=1.2,tcl=0.015)
  for (i in 1:4) 	lines(rep(i+2014,2),dens.year.sum[sp,i,2:3],col="blue2")
  mtext(side=3,at=2015,SpTable[sp],cex=1.1,adj=0)
}
graphics.off()




