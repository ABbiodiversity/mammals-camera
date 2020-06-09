# R Processing of camera data for native mammals with second pass information
# 2018 version
# This step converts the individual images into series, runs detection distance models,
# figures out times of series (incl. probablistic gap-leaving), summarizes times cameras operating,
# calculates lure effects, compares camera models or years, and does summary output for those steps for reporting

# Convert distance codes to columns
d<-read.csv("./beta/Dave-Objects/Native mammals only Feb 2019 ALL.csv")

# Combined files from each year.  Make sure that seconds have not been lost in date_time field in Excel.
d$common_name<-as.character(d$common_name)
d$common_name<-ifelse(d$common_name=="Mule Deer","Mule deer",d$common_name)  # Check for any other multiple-format names
d$deployment<-as.character(d$deployment)
# Change all WHEC letter additions to lower case (inconsistent in data)
i<-which(substr(d$deployment,1,4)=="WHEC")
d$deployment[i]<-gsub("C","c",d$deployment[i])
d$deployment[i]<-gsub("D","d",d$deployment[i])
d$deployment[i]<-gsub("N","n",d$deployment[i])
d$deployment[i]<-gsub("U","u",d$deployment[i])
d$deployment[i]<-gsub("WHEc","WHEC",d$deployment[i])
# Other changes
d$number_individuals<-as.numeric(as.character(d$number_individuals))
d$number_individuals<-ifelse(is.na(d$number_individuals),1,d$number_individuals)  # Assume 1 individual if not recorded
d$number_individuals<-ifelse(d$number_individuals==0,1,d$number_individuals)  # I checked a few of these, and they all had 1 individual
# Delete and modify incorrect 2015 records [OG-CITSCI-1537-31 also has 2015 dates, but no native mammals.  Others have been fixed in the database]
d$date_time_taken<-as.character(d$date_time_taken)
i<-which(d$deployment=="OGW-CITSCI-1538-31")
d$date_time_taken[i]<-paste("2016",substr(d$date_time_taken[i],5,nchar(d$date_time_taken[i])),sep="")
# Convert NA to "X" to avoid confusion with R's NA
d$distance<-as.character(d$distance)
d$distance[which(is.na(d$distance))]<-"X"
d$distance[which(d$distance=="NA, NA")]<-"X, X"
d$distance[which(d$distance=="NA, NA, NA")]<-"X, X, X"
table(d$distance)  # Check that no other odd codes here

# Make columns of number of individuals at each pole position
q<-gregexpr("A",d$distance)
w<-do.call(rbind, lapply(q, function(x) length(x[x>-1])))
d$PoleAt<-w
q<-gregexpr("B",d$distance)
w<-do.call(rbind, lapply(q, function(x) length(x[x>-1])))
d$PoleBehind<-w
q<-gregexpr("F",d$distance)
w<-do.call(rbind, lapply(q, function(x) length(x[x>-1])))
d$PoleFront<-w
q<-gregexpr("IC",d$distance)
w<-do.call(rbind, lapply(q, function(x) length(x[x>-1])))
d$PoleIC<-w
q<-gregexpr("IP",d$distance)
w<-do.call(rbind, lapply(q, function(x) length(x[x>-1])))
d$PoleIP<-w
q<-gregexpr("X",d$distance)
w<-do.call(rbind, lapply(q, function(x) length(x[x>-1])))
d$PoleNA<-w

write.table(d,file="C:/Dave/ABMI/Data/Mammals/2018/Native mammals only with distances Feb 2019 ALL.csv",sep=",",row.names=F)

# Summarize individuals per image.  This is just for reporting
q<-table(d$common_name,ifelse(d$number_individuals>10,10,d$number_individuals))
q1<-data.frame(Species=row.names(q),array(q,dim(q)),Total=as.numeric(by(d$number_individuals,d$common_name,sum)),Mean=as.numeric(by(d$number_individuals,d$common_name,mean)))
write.table(q1,file="C:/Dave/ABMI/Cameras/2018 analysis/Individuals per photo/Number of individuals per image by species Feb 2019 ALL.csv",sep=",",col.names=NA)

# R Identify series.
# Use gaps-in-series results to modify series where available.
# Also, plot gap-in-series versus whether the animal left or not.
# NOTE: (Almost) all gaps were done in 2015, but only select ones were done in 2016-2018
d$gap_class<-as.character(d$gap_class)
d$Time1<-strptime(d$date_time_taken,format="%Y-%m-%d %H:%M:%S")
d1<-d[order(d$deployment,d$Year,d$common_name,d$Time1),]  # Sort by site, year, species, time
d1$SeriesNum<-0  # Create this field for later use
d1$GapCheck<-""  # This variable was to flag 20-120 seconds for checking images - not done any more.
                # But still used below to indicate which gaps need to have the probablistic gap-leaving models applied
d2<-d1[-1,]
d3<-d1[-nrow(d1),]  # Offset by one record, to allow comparison of time between subsequent images
d1$DiffTime<-c(99999,as.vector(difftime(d2$Time1,d3$Time1,units="secs")))
# Look for NA DiffTimes and fill these in by hand if necessary with time from the previous image.  (I don't know why R can't calculate time differences for a few pairs of records...)
i<-which(is.na(d1$DiffTime))
data.frame(i,d1$Time1[i-1],d1$Time1[i])
d1$DiffTime[which(is.na(d1$DiffTime))]<-c(9999,20,9999,9999,1,2,rep(9999,4),30,5,2,rep(9999,13))  # These numbers are generated manually by checking the pairs of records with NA DiffTimes, and are in the same order as i
# Then use changes in species or site, or long times between images to identify the start of each new series
d1$DiffSp<-c(TRUE,d2$common_name!=d3$common_name)
d1$DiffSite<-c(TRUE,d2$deployment!=d3$deployment)
x<-sign(d1$DiffSp+d1$DiffSite+(d1$DiffTime>120))  # 0 if no change in species, site and time gap<120 seconds; 1 otherwise
d1$SeriesNum<-cumsum(x)
d1$GapCheck<- 1-sign(d1$DiffSp+d1$DiffSite+((d1$DiffTime>120) | (d1$DiffTime<20)) )
# 1 if all false - i.e., not different species, not different site, and time difference is between 20 and 120 seconds
d1$gap_class<-c("",d1$gap_class[-nrow(d1)])
# Shift these down one, because GapCheck, DiffTime are assigned to image after gap, while gap_class in database is assigned to image before gap

# Plot gap length versus left/not
jpeg(file="C:/Users/mabec/Documents/R/SC-Mammals/beta/Gap check versus gap length 2.jpg",width=600,height=500)

par(mai=c(1,1,0.3,2))

d4<-d1[d1$gap_class!="" & d1$GapCheck==1,]

plot(jitter(d4$DiffTime[d4$GapCheck==1]),jitter(ifelse(d4$gap_class[d4$GapCheck==1]=="P",0,1),0.2),cex=0.5,col="#0000FF20",xlab="Gap length (s)",ylab="Probability of Leaving",yaxt="n")

axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1),lab=c("Stayed","0.2","0.4","0.6","0.8","Left"),cex.axis=1.3,las=2)

lines(smooth.spline(d4$DiffTime[d4$GapCheck==1],ifelse(d4$gap_class[d4$GapCheck==1]=="P",0,1),df=3),lwd=2)

sp<-"White-tailed Deer"
c1<-"chocolate"
m<-smooth.spline(d4$DiffTime[d4$GapCheck==1 & d4$common_name==sp],ifelse(d4$gap_class[d4$GapCheck==1 & d4$common_name==sp]=="P",0,1),df=3)
lines(predict(m),lwd=1,col=c1)
mtext(side=4,at=predict(m)$y[length(predict(m)$y)]-0.01,sp,col=c1,las=1,cex=1.2)
sp<-"Mule deer"
c1<-"chocolate4"
m<-smooth.spline(d4$DiffTime[d4$GapCheck==1 & d4$common_name==sp],ifelse(d4$gap_class[d4$GapCheck==1 & d4$common_name==sp]=="P",0,1),df=3)
lines(predict(m),lwd=1,col=c1)
mtext(side=4,at=predict(m)$y[length(predict(m)$y)],sp,col=c1,las=1,cex=1.2)
sp<-"Moose"
c1<-"black"
m<-smooth.spline(d4$DiffTime[d4$GapCheck==1 & d4$common_name==sp],ifelse(d4$gap_class[d4$GapCheck==1 & d4$common_name==sp]=="P",0,1),df=3)
lines(predict(m),lwd=1,col=c1)
mtext(side=4,at=predict(m)$y[length(predict(m)$y)],sp,col=c1,las=1,cex=1.2)
sp<-"Elk (wapiti)"
c1<-"brown2"
m<-smooth.spline(d4$DiffTime[d4$GapCheck==1 & d4$common_name==sp],ifelse(d4$gap_class[d4$GapCheck==1 & d4$common_name==sp]=="P",0,1),df=3)
lines(predict(m),lwd=1,col=c1)
mtext(side=4,at=predict(m)$y[length(predict(m)$y)]-0.02,sp,col=c1,las=1,cex=1.2)
sp<-"Pronghorn"
c1<-"orange2"
m<-smooth.spline(d4$DiffTime[d4$GapCheck==1 & d4$common_name==sp],ifelse(d4$gap_class[d4$GapCheck==1 & d4$common_name==sp]=="P",0,1),df=3)
lines(predict(m),lwd=1,col=c1)
mtext(side=4,at=predict(m)$y[length(predict(m)$y)],sp,col=c1,las=1,cex=1.2)
sp<-"Black Bear"
c1<-"darkred"
m<-smooth.spline(d4$DiffTime[d4$GapCheck==1 & d4$common_name==sp],ifelse(d4$gap_class[d4$GapCheck==1 & d4$common_name==sp]=="P",0,1),df=3)
lines(predict(m),lwd=1,col=c1)
mtext(side=4,at=predict(m)$y[length(predict(m)$y)],sp,col=c1,las=1,cex=1.2)
sp<-"Gray Wolf"
c1<-"gray50"
m<-smooth.spline(d4$DiffTime[d4$GapCheck==1 & d4$common_name==sp],ifelse(d4$gap_class[d4$GapCheck==1 & d4$common_name==sp]=="P",0,1),df=3)
lines(predict(m),lwd=1,col=c1)
mtext(side=4,at=predict(m)$y[length(predict(m)$y)]+0.03,sp,col=c1,las=1,cex=1.2)
sp<-"Coyote"
c1<-"slateblue3"
m<-smooth.spline(d4$DiffTime[d4$GapCheck==1 & d4$common_name==sp],ifelse(d4$gap_class[d4$GapCheck==1 & d4$common_name==sp]=="P",0,1),df=3)
lines(predict(m),lwd=1,col=c1)
mtext(side=4,at=predict(m)$y[length(predict(m)$y)],sp,col=c1,las=1,cex=1.2)
sp<-"Canada Lynx"
c1<-"cyan1"
m<-smooth.spline(d4$DiffTime[d4$GapCheck==1 & d4$common_name==sp],ifelse(d4$gap_class[d4$GapCheck==1 & d4$common_name==sp]=="P",0,1),df=3)
lines(predict(m),lwd=1,col=c1)
mtext(side=4,at=predict(m)$y[length(predict(m)$y)]+0.02,sp,col=c1,las=1,cex=1.2)
sp<-"Fisher"
c1<-"green2"
m<-smooth.spline(d4$DiffTime[d4$GapCheck==1 & d4$common_name==sp],ifelse(d4$gap_class[d4$GapCheck==1 & d4$common_name==sp]=="P",0,1),df=3)
lines(predict(m),lwd=1,col=c1)
mtext(side=4,at=predict(m)$y[length(predict(m)$y)]+0.03,sp,col=c1,las=1,cex=1.2)
sp<-"Marten"
c1<-"green4"
m<-smooth.spline(d4$DiffTime[d4$GapCheck==1 & d4$common_name==sp],ifelse(d4$gap_class[d4$GapCheck==1 & d4$common_name==sp]=="P",0,1),df=3)
lines(predict(m),lwd=1,col=c1)
mtext(side=4,at=predict(m)$y[length(predict(m)$y)],sp,col=c1,las=1,cex=1.2)
sp<-"Snowshoe Hare"
c1<-"blue"
m<-smooth.spline(d4$DiffTime[d4$GapCheck==1 & d4$common_name==sp],ifelse(d4$gap_class[d4$GapCheck==1 & d4$common_name==sp]=="P",0,1),df=3)
lines(predict(m),lwd=1,col=c1)
mtext(side=4,at=predict(m)$y[length(predict(m)$y)],sp,col=c1,las=1,cex=1.2)
sp<-"Porcupine"
c1<-"turquoise4"
m<-smooth.spline(d4$DiffTime[d4$GapCheck==1 & d4$common_name==sp],ifelse(d4$gap_class[d4$GapCheck==1 & d4$common_name==sp]=="P",0,1),df=3)
lines(predict(m),lwd=1,col=c1)
mtext(side=4,at=predict(m)$y[length(predict(m)$y)]+0.02,sp,col=c1,las=1,cex=1.2)
graphics.off()

# Gap-leaving probability models for each gap group
g<-read.csv("C:/Dave/ABMI/Cameras/2018 analysis/Series and gaps/Gap groups.csv")
# Prepared by listing species from d4 above and adding gap groups manually
p.gap<-array(NA,c(max(g$GapGroup),120))
# Probability gap left for each gap group and 1-120 seconds (only filled in for 20-120 seconds)
for (i in 1:max(g$GapGroup)) {
  sp.list<-g$Sp[g$GapGroup==i]  # Species that are in that gap group
  d5<-d4[d4$common_name %in% sp.list,]  # All records for those species
  p.gap[i,20:120]<-predict(smooth.spline(d5$DiffTime[d5$GapCheck==1],ifelse(d5$gap_class[d5$GapCheck==1]=="P",0,1),df=3),x=20:120)$y
  # Smooth spline used to model probability of leaving.  Assume only code "P" used for animals that stayed
}

row.names(p.gap)<-c("Most ungulates","Moose","Small carnivores","Canids cougar","Bears","Small mammals")

save(p.gap,file="C:/Dave/ABMI/Cameras/2018 analysis/Series and gaps/Gap models 2018.RData")  # For use by other summaries of camera data

# Re-run series numbering to include checked gap information.
# This is for when 20-120 gaps have been checked from images, so don't need to use probablistic gap-leaving
d1$SeriesNum<-0
x<-sign(d1$DiffSp+d1$DiffSite+(d1$gap_class=="L" | d1$gap_class=="N")+(d1$DiffTime>120))
# 0 if no change in species, site, gap class is not L (=left) or N (=followed by NONE record) and time gap<120 seconds; 1 otherwise
# Note: The series stays the same if it is a checkable gap but there was no check or check was "U" (unclear).
# In that case, the total gap time is subtracted and the probablistic gap time is added below
d1$SeriesNum<-cumsum(x)
d1$DiffSeries1 <- x
d1$ProbGap<-ifelse(d1$GapCheck==1 & (d1$gap_class=="" | d1$gap_class=="U"),1,0)  # Gaps needing probablistic time assignment
# Get rid of hyphen in CMU deployments - not used consistently
d1$deployment<-as.character(d1$deployment)
cmu<-which(substr(d1$deployment,1,3) %in% c("CAL","CHR","DEW","FMM","LLB","LRN","MAC","MCC","WAB")==TRUE)
d1$deployment[cmu]<-gsub("-","",d1$deployment[cmu])
write.table(d1,file="C:/Dave/ABMI/Cameras/2018 analysis/Series and gaps/Mammals images with series and checkable gaps Feb 2019 ALL.csv",sep=",",col.names=NA)

# Total time for each series
d1<-read.csv("C:/Dave/ABMI/Cameras/2018 analysis/Series and gaps/Mammals images with series and checkable gaps Feb 2019 ALL.csv")
d1$Time1<-strptime(d1$Time1,format="%Y-%m-%d %H:%M:%S")  # Will read in as factor otherwise
g<-read.csv("C:/Dave/ABMI/Cameras/2018 analysis/Series and gaps/Gap groups.csv")
# Needed to assign probability of leaving to probablistic gaps for each species based on its gap group
# Average time between photos within series (TBP), by species
# - add this much to each series (half at start, half at end) to account for time before first photo and after last photo
#   (incl. for singleton series)
i1<-c(0,d1$SeriesNum[-nrow(d1)])  # Series numbers shifted down one
j<-which(d1$SeriesNum!=i1 & d1$ProbGap!=1 )
# These are the rows in d1 where new series start AND rows with a probablistic gap.
# Don't use these DiffTimes because they are not within the series (they are the times since the previous series)
# OR they are uncertain gaps
TBP<-by(d1$DiffTime[-j],d1$common_name[-j],mean)  # Time Between Photos = Average time difference within series, by species
TBP.n<-by(d1$DiffTime[-j],d1$common_name[-j],length)  # Sample size for TBP
q<-data.frame(Species=names(TBP),TBP=as.numeric(TBP),TBP.n=as.numeric(TBP.n))
write.table(q,file="C:/Dave/ABMI/Cameras/2018 analysis/Series and gaps/Table of time between photos within series by species.csv",sep=",",row.names=F)
# Time for each series (TFS).
# These start out as 0 for singletons, or total time between photos for series >1, then extra time is added for the start and end (including for singletons).
# Then time is taken off for probability that animal left a probablistic gap.
all.series<-sort(unique(d1$SeriesNum))
TFS<-rep(0,length(all.series))  # Default is 0, for singletons
TFS.temp<-by(d1$DiffTime[-j],d1$SeriesNum[-j],sum)  # This is only for series with >1 photo
TFS[match(names(TFS.temp),all.series)]<-as.numeric(TFS.temp)
TFS.extra<-as.numeric(TBP)[match(d1$common_name,names(TBP))]  # For extra time at start and end, including for singletons
TFS.extra<-ifelse(is.na(TFS.extra),mean(as.numeric(TBP),na.rm=T),TFS.extra)
# NA when there aren't any series >1 for a species.  Use overall average in that case
TFS.extra<-by(TFS.extra,d1$SeriesNum,mean)  # Summarize by series
TFI<-as.numeric(TFS+TFS.extra)
# At this point, TFI includes the full gap time for probablistic gaps, so subtract those.
i<-which(d1$ProbGap==1)
group<-g$GapGroup[match(d1$common_name[i],g$Sp)]
PGT<-NULL
TFI1<-TFI  # To check the next step
for (j in 1:length(i)) {  # Cycle through each series with a probalistic gap
  PGT[j]<-p.gap[group[j],d1$DiffTime[i][j]] * d1$DiffTime[i][j]  # The probability-of-leaving-adjusted  gap time for probablistic gaps - this is the time to subtract
  TFI[d1$SeriesNum[i[j]]]<-TFI[d1$SeriesNum[i[j]]]-PGT[j]  # Subtract that time from the appropriate series
}
dplot(TFI1,TFI,xlab="Total series time (s)",ylab="Series time adjusted for probablistic gaps (s)",cex=0.6)
# Points should be on or below 1:1 line (on for no probablistic gaps in that series, below if there are probablistic gaps)

# Compile data.frame with one record for each mammal series
# Mean number of animals and total number of photos per series.
# Note: Mean number of animals is not time-weighted - i.e., each image in a series is weighted equally
mean.animals<-as.numeric(by(d1$number_individuals,d1$SeriesNum,function(x) mean(ifelse(is.na(x),1,x))))
n.photos<-as.numeric(by(d1$SeriesNum,d1$SeriesNum,length))
# Get camera location, species, date for each series, by finding first record for each
d1<-d1[order(d1$SeriesNum,d1$Time1),]
same.series<-c("F",ifelse(d1$SeriesNum[2:nrow(d1)]==d1$SeriesNum[1:(nrow(d1)-1)],"T","F"))
# Marks each record as being the same series as the previous one (T), or not (F)
d2<-d1[same.series=="F",]  # Use only the first record of each series
d2<-d2[order(d2$SeriesNum),]
q<-(by(d1$SeriesNum,d1$SeriesNum,length))
dplot(1:nrow(d2),match(names(q),d2$SeriesNum),cex=0.3)
# Make sure this is a straight line - i.e., that n.photos and d2 have same series in same order
d2$Photos<-n.photos
d2$Duration<-TFI
d2$Photos<-n.photos
d2$Individuals<-mean.animals
# Get Lured information from separate camera data file
d3<-read.csv("C:/Dave/ABMI/Cameras/2018 analysis/Start end times/Start end times Lure Feb 2019 ALL.csv")
# Exported from Access - separately for each year, then combined
names(d3)[which(names(d3)=="Deployment")]<-"deployment"
d2<-merge(d2,d3,all.x=T)
d<-data.frame(Deployment=d2$deployment,Year=d2$Year,Lured=d2$Lure,DateTime=d2$Time1,Series=d2$SeriesNum,Photos=d2$Photos,Species=d2$common_name,Sex=d2$sex,Age=d2$age_class,Multiple=d2$multiple_animals,Individuals=d2$Individuals,Duration=d2$Duration)
# d is data file of series

# Summarize time cameras operating. Used below to get rid of mammal series from times when the camera wasn't operating properly
t<-read.csv("C:/Dave/ABMI/Cameras/2018 analysis/Start end times/Start end times Lure Feb 2019 ALL.csv")
t <- read.csv(paste0(abmisc,"processed/dave-objects/Start end times Lure Feb 2019 ALL.csv"))
t$DeploymentYear<-paste(t$Deployment,t$Year,sep="_")
# Truncate to just cameras that have images
# dwti <- read.csv("C:/Dave/ABMI/Data/Mammals/2018/All deployments with images Feb 2019 ALL.csv")
dwti<-read.csv(paste0(abmisc,"processed/dave-objects/All deployments with images July 2019 ALL.csv"))
# Deployments with images, compiled from images database.
# Needed, because there is metadata from many non-ABMI dpeloyments that have not been tagged
dwti$DeploymentYear<-paste(dwti$Deployment,dwti$Year,sep="_")
t<-t[!is.na(match(t$DeploymentYear,dwti$DeploymentYear)),]  # Get rid of any rows in t that do not have tagged images
site.list<-unique(as.character(t$DeploymentYear))

t$StartTime1<-strptime(t$StartTime,format="%Y-%m-%d %H:%M:%S")  # Will read in as factor otherwise
t$EndTime1<-strptime(t$EndTime,format="%Y-%m-%d %H:%M:%S")  # Will read in as factor otherwise
t.extra<-read.csv(paste0(abmisc,"processed/dave-objects/Intermediate END START pairs Feb 2019 ALL.csv"))
t.extra$DeploymentYear<-paste(t.extra$Deployment,t.extra$Year,sep="_")
t.extra<-t.extra[!is.na(match(t.extra$DeploymentYear,dwti$DeploymentYear)),]
# Get rid of any rows in t.extra that do not have tagged images
# (shouldn't be any of these, because wouldn't know about intermediate END/STOP pairs if the images haven't been tagged, but you never know)

t.extra$StartEndTime1<-strptime(t.extra$date_time_original,format="%Y-%m-%d %H:%M:%S")  # Will read in as factor otherwise

# To track each julian day that cameras were operating
time.by.day<-array(0,c(length(site.list),366))

# To track days since Jan 1 2009 for deployments in all years.  Add extra days each year
# Do time since 2009 first, then calculate time-by-day from that, to deal with deployments that cross Jan 1
time.since.2009<-array(0,c(length(site.list),365+365+365+366+365+365+365+366+365+365+365))

for (i in 1:length(site.list)) {
  t1<-t[t$DeploymentYear==site.list[i],]
  for (k in 1:nrow(t1)) {  # There can be more than one line here, when the same camera is dpeloyed more than once in a year
    j1<-floor(julian(t1$StartTime1[k],origin=strptime("2009-01-01 00:00:00",format="%Y-%m-%d %H:%M:%S")))  # For days since 2009.  "floor" so that first day is included in operating days
    j2<-floor(julian(t1$EndTime1[k],origin=strptime("2009-01-01 00:00:00",format="%Y-%m-%d %H:%M:%S")))  # "floor" so that last day (but not following day) is included in operating days
    if (!is.na(j1) & !is.na(j2)) {   # Some blank lines in start/end date file for undeployed or uncollected cameras
      time.since.2009[i,(j1:j2)]<-1
    }
  } # Next line of start and end times for that DeploymentYear
  if (site.list[i] %in% t.extra$DeploymentYear) {  # Take off time(s) when camera wasn't working
    t.extra1<-t.extra[as.character(t.extra$DeploymentYear)==as.character(site.list[i]),]
    for (j in seq(1,(nrow(t.extra1)-1),2)) {  # Assumes all extra times are formatted as END/START row pairs
      j1<-ceiling(julian(t.extra1$StartEndTime1[j],origin=strptime("2009-01-01 00:00:00",format="%Y-%m-%d %H:%M:%S")))  # For days since 2009.  "ceiling" so that day of failure is not excluded from operating days
      j2<-floor(julian(t.extra1$StartEndTime1[j+1],origin=strptime("2009-01-01 00:00:00",format="%Y-%m-%d %H:%M:%S")))  # "floor" so that day of recovery is included in operating days
      if (j2>j1) time.since.2009[i,j1:(j2-1)]<-0  # If these aren't the same or consecutive days, exclude days between.  This keeps the day of failure and day of recovery as operating days
    }  # Next intermediate END/START pair for that deployment
  }  # End if for END-START pairs in that deployment
}  # Next DeploymentYear

# Then calculate time.by.day from time.since.2009
days.per.year<-c(365,365,365,366,365,365,365,366,365,365,365)  # 2009 to 2019. Add extra days each year
Jan1<-cumsum(c(1,days.per.year))  # Day of Jan 1 from 2009 to 2020
yrnum<-julday<-NULL  # Years from 2009 and Julian day for each column
for (i in 1:ncol(time.since.2009)) {
  yrnum[i]<-sum(i>=Jan1)
  julday[i]<- i-Jan1[yrnum[i]]+1
}
for (i in 1:ncol(time.by.day)) {
  time.by.day[,i]<-rowSums(time.since.2009[,which(julday==i)])
}
rownames(time.by.day)<-rownames(time.since.2009)<-site.list
# Export
write.table(time.by.day,file="C:/Dave/ABMI/Cameras/2018 analysis/Start end times/Camera operating days Feb 2019 ALL.csv",sep=",",col.names=NA)
write.table(time.since.2009,file="C:/Dave/ABMI/Cameras/2018 analysis/Start end times/Camera operating days since 2009 Feb 2019 ALL.csv",sep=",",col.names=NA)

# Then remove mammal records after final end or in intermediate end-start periods - using operating days in time.since.2009
#d$DateTime<-strptime(d$DateTime,format="%Y-%m-%d %H:%M:%S")  # May be factor otherwise
d$DeploymentYear<-paste(d$Deployment,d$Year,sep="_")
x<-rep(0,nrow(d))  # 0 for keep, change to 1 for delete
for (i in 1:nrow(d)) {
  j<-floor(julian(d$DateTime[i],origin=strptime("2009-01-01 00:00:00",format="%Y-%m-%d %H:%M:%S")))  # For days since 2009
  k<-match(d$DeploymentYear[i],rownames(time.since.2009))
  if (is.na(k) | is.na(j) | time.since.2009[k,j]==0) x[i]<-1
  # Mammal record when camera not operating.
  # NA's for deployments that are not in meta-data (and hence not in time.since.2009) <-- delete these.
  # I think is.na(j) happens with illegitimate times - when day-light savings is happening
}
# Worth checking records with x==1 here - make sure they are all legitimate records that should be deleted.
d<-d[x==0,]
write.table(d,file="C:/Dave/ABMI/Cameras/2018 analysis/Mammals by series Feb 2019 ALL.csv",sep=",",row.names=F)

# Summary table of series (includes off-grid, other projects) - for reporting
for (y in 2015:2018) {  # Ignoring haphazard cameras of earlier years
  d.y<-d[d$Year==y,]
  q<-table(d.y$Species)
  td<-as.numeric(by(d.y$Duration,d.y$Species,sum))
  dur.mean<-td/as.numeric(q)
  sites<-as.numeric(by(d.y$Deployment,d.y$Species,function(x) length(unique(x))))
  indiv<-as.numeric(by(d.y$Individuals,d.y$Species,sum))
  images<-as.numeric(table(d1$common_name[d1$Year==y]))  # d1 is still the original native mammal image dataset
  if (y==2015) q1<-data.frame(Species=names(q),Year=y,Sites=sites,Series=as.numeric(q),Individuals=indiv,Images=images,Duration.total=td,Duration.mean=dur.mean)
  if (y>2015) q1<-rbind(q1,data.frame(Species=names(q),Year=y,Sites=sites,Series=as.numeric(q),Individuals=indiv,Images=images,Duration.total=td,Duration.mean=dur.mean))
}
write.table(q1,file="C:/Dave/ABMI/Cameras/2018 analysis/Summaries/Series summary table.csv",sep=",",row.names=F)

# Summarize total time, and winter, summer time
time.by.day<-read.csv("C:/Dave/ABMI/Cameras/2018 analysis/Start end times/Camera operating days Feb 2019 ALL.csv")
time.by.day$X<-as.character(time.by.day$X)
q<-gregexpr("_",time.by.day$X)  # Row names get read in as column X
deployment<-substr(time.by.day$X,1,as.numeric(lapply(q,function(x) x[length(x)]))-1)
year<-substr(time.by.day$X,as.numeric(lapply(q,function(x) x[length(x)]))+1,nchar(time.by.day$X))
time.by.day$X<-NULL
j.summer.start<-106  # April 16
j.summer.end<-288  # Oct 15
time.total<-rowSums(time.by.day)
time.summer<-rowSums(time.by.day[,j.summer.start:j.summer.end])
time.winter<-rowSums(time.by.day[,c(1:(j.summer.start-1),(j.summer.end+1):ncol(time.by.day))])
q<-data.frame(Deployment=deployment,Year=year,time.total,time.summer,time.winter)
write.table(q,file="C:/Dave/ABMI/Cameras/2018 analysis/Start end times/Summary of camera operating times Feb 2019 ALL.csv",sep=",",row.names=FALSE)

# Plot deployment time figures - for reporting
time.by.day<-read.csv(file="C:/Dave/ABMI/Cameras/2018 analysis/Start end times/Camera operating days Feb 2019 ALL.csv")
names(time.by.day)[1]<-"DeploymentYear"
x<-as.character(time.by.day$DeploymentYear)
year<-as.numeric(substr(x,nchar(x)-3,nchar(x)))
ordernum<-ifelse(substr(x,1,4)=="ABMI",ifelse(substr(x,1,6)=="ABMI-W",2,1),3)  # Put ABMI first, then ABMI-W, then OG and other studies
time.by.day1<-time.by.day[order(year,ordernum,x),]  # Sort by year, grouping, alphabetical deployment name
x<-as.character(time.by.day1$DeploymentYear)  # And re-do these variables in the new order
year<-as.numeric(substr(x,nchar(x)-3,nchar(x)))
for (yr in 1:4) {
  fname<-paste("C:/Dave/ABMI/Cameras/2018 analysis/Start end times/Deployment operating days ",yr+2014,".png",sep="")
  png(fname,width=500,height=800)
  t<-time.by.day1[year==yr+2014,]
  x<-as.character(t$DeploymentYear)  # And re-do these variables in the new order
  col1<-ifelse(substr(x,1,4)=="ABMI",ifelse(substr(x,1,6)=="ABMI-W","blue2","green2"),"orange2")
  y<-1:nrow(t)
  dplot(0,0,xlim=c(1,366),ylim=c(1,max(y)),xlab="",xaxt="n",typ="n",yaxt="n")
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

# Histograms of operating hours - for reporting, but not saved, so just cut and paste images if needed
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


# Summarize distance wrt pole into detection distance models
d<-read.csv("C:/Dave/ABMI/Data/Mammals/2018/Distance output for native mammals unlured 2015 to 2018.csv")  # Exported from Access.  Only records with non-NA non-blank Distance field.  Lured sites removed.
d$number_individuals<-ifelse(is.na(d$number_individuals),1,d$number_individuals)
# Convert distance coding to number of animals in each distance class
x<-gregexpr("B",d$distance)
d$Behind<-unlist(lapply(x,function(x) sum(sign(x))))
d$Behind<-ifelse(d$Behind==-1,0,d$Behind)
x<-gregexpr("F",d$distance)
d$Front<-unlist(lapply(x,function(x) sum(sign(x))))
d$Front<-ifelse(d$Front==-1,0,d$Front)
x<-gregexpr("A",d$distance)
d$At<-unlist(lapply(x,function(x) sum(sign(x))))
d$At<-ifelse(d$At==-1,0,d$At)
x<-gregexpr("IC",d$distance)
d$IC<-unlist(lapply(x,function(x) sum(sign(x))))
d$IC<-ifelse(d$IC==-1,0,d$IC)
x<-gregexpr("IP",d$distance)
d$IP<-unlist(lapply(x,function(x) sum(sign(x))))
d$IP<-ifelse(d$IP==-1,0,d$IP)
d$n<-d$Front+d$Behind+d$At+d$IP+d$IC

d$Julian<-as.numeric(strftime(d$date_time_taken,format="%j"))
names(d)[which(names(d)=="deployment")]<-"Deployment"
d$Time1<-strptime(d$date_time_taken,format="%Y-%m-%d %H:%M:%S")  # Will read in as factor otherwise

# Add veg+HF information for camera point
s<-read.csv("C:/Dave/ABMI/Data/Site info/2018/Combined vegHF soilHF and detection distance veg for cameras Feb 2019.csv")
# Detection distance veg is in this file, processed in Excel to standard categories
s<-s[,c("deployment","Year","DeploymentYear","VegForDetectionDistance")]
# Use only the reduced veg here
s$VegForDetectionDistance<-ifelse(as.character(s$VegForDetectionDistance)=="WetShrub","Shrub",as.character(s$VegForDetectionDistance))  # WetShrub category not used consistently
names(s)[which(names(s)=="VegForDetectionDistance")]<-"VegHF"
names(s)[which(names(s)=="deployment")]<-"Deployment"
d$DeploymentYear<-paste(d$Deployment,d$Year,sep="_")
d<-merge(d,s[,c("DeploymentYear","VegHF")],by="DeploymentYear",all.x=TRUE)
d$VegHF<-as.character(d$VegHF)  # Get rid of factor levels that aren't used
# Make combined veg types
d$VegHF1<-"Wet"
d$VegHF1<-ifelse(d$VegHF=="Conif" | d$VegHF=="Decid", "ConifDecid", d$VegHF1)
d$VegHF1<-ifelse(d$VegHF=="Grass" | d$VegHF=="Shrub", "GrassShrub", d$VegHF1)
d$VegHF1<-ifelse(d$VegHF=="HF","HF",d$VegHF1)
d$VegHF2<-"GrassWater"
d$VegHF2<-ifelse(d$VegHF=="Conif" | d$VegHF=="Decid" | d$VegHF=="WetTreed","Treed",d$VegHF2)
d$VegHF2<-ifelse(d$VegHF=="Shrub","Shrub",d$VegHF2)
d$VegHF2<-ifelse(d$VegHF=="HF","HF",d$VegHF2)
d$VegHF3<-"Wet"
d$VegHF3<-ifelse(d$VegHF=="Conif" | d$VegHF=="Decid","ConifDecid",d$VegHF3)
d$VegHF3<-ifelse(d$VegHF=="Grass" | d$VegHF=="Shrub" | d$VegHF=="HF","GrassShrubHF",d$VegHF3)
d$VegHF4<-"GrassWaterHF"
d$VegHF4<-ifelse(d$VegHF=="Conif" | d$VegHF=="Decid" | d$VegHF=="WetTreed","Treed",d$VegHF4)
d$VegHF4<-ifelse(d$VegHF=="Shrub","Shrub",d$VegHF4)
d$VegHF5<-"GrassWaterShrubHF"
d$VegHF5<-ifelse(d$VegHF=="Conif" | d$VegHF=="Decid" | d$VegHF=="WetTreed","Treed",d$VegHF5)

# Combine rare species with more common ones
# Edit this if new species appear
d$SpGroup<-"SmallOpen"  # Badger, Raccoon, Groundhog, Richardson's Ground Squirrel, Striped Skunk, Voles, White-tailed Jack Rabbit, Beaver
d$SpGroup<-ifelse(d$common_name=="Fisher" | d$common_name=="Marten" | d$common_name=="Porcupine" | d$common_name=="Red Squirrel" | d$common_name=="Snowshoe Hare" | d$common_name=="Mink" | d$common_name=="Weasels and Ermine" | d$common_name=="Northern Flying Squirrel","SmallForest",d$SpGroup)
d$SpGroup<-ifelse(d$common_name=="Canada Lynx" | d$common_name=="Bobcat","Lynx",d$SpGroup)
d$SpGroup<-ifelse(d$common_name=="Bison" | d$common_name=="Moose" | d$common_name=="Woodland Caribou", "BigForestUngulates",d$SpGroup)
d$SpGroup<-ifelse(d$common_name=="Black Bear" | d$common_name=="Grizzly bear","Bear",d$SpGroup)
d$SpGroup<-ifelse(d$common_name=="White-tailed Deer" | d$common_name=="Deer","WTDeer",d$SpGroup)
d$SpGroup<-ifelse(d$common_name=="Gray Wolf" | d$common_name=="Cougar" | d$common_name=="Wolverine" | d$common_name=="Wloves, Coyotes and Allies","BigForestCarnivores",d$SpGroup)
d$SpGroup<-ifelse(d$common_name=="Coyote" | d$common_name=="Red fox" | d$common_name=="Foxes" | d$common_name=="Swift fox","CoyoteFox",d$SpGroup)
d$SpGroup<-ifelse(d$common_name=="Mule deer" | d$common_name=="Elk (wapiti)" | d$common_name=="Pronghorn" | d$common_name=="Bighorn sheep",as.character(d$common_name),d$SpGroup)
d$SpGroup<-ifelse(d$common_name=="Mountain goat","Bighorn sheep",d$SpGroup)
table(d$common_name,d$SpGroup)  # Check that these make sense - especially species in default SmallOpen group
SpTable<-unique(d$SpGroup)

# Estimate time between images behind and in front of pole, by species
# This is just extra information, to think about whether detection distance modeling should include time between images in each pole position, rather than just the number of images
d1<-d[order(d$DeploymentYear,d$common_name,d$Time1),]
d1.1<-d1[2:nrow(d1),]  # Offset one row
d1<-d1[-nrow(d1),]
time.diff<-difftime(d1.1$Time1,d1$Time1,units="secs")
time.diff.front<-time.diff.behind<-time.diff.front.n<-time.diff.behind.n<-NULL  # One value for each SpGroup
for (sp in 1:length(SpTable)) {
  time.diff.front1<-time.diff[d1.1$SpGroup==SpTable[sp] & d1.1$DeploymentYear==d1$DeploymentYear & d1.1$common_name==d1$common_name & time.diff<120 & time.diff>0 & d1.1$Front==d1.1$n & d1$Front==d1$n]  # Time differences for consecutive records of a species in that group where both records were in front
  time.diff.behind1<-time.diff[d1.1$SpGroup==SpTable[sp] & d1.1$DeploymentYear==d1$DeploymentYear & d1.1$common_name==d1$common_name & time.diff<120 & time.diff>0 & d1.1$Behind==d1.1$n & d1$Behind==d1$n]  # Time differences for consecutive records of a species in that group where both records were behind
  time.diff.front[sp]<-mean(time.diff.front1,na.rm=T)
  time.diff.behind[sp]<-mean(time.diff.behind1,na.rm=T)
  time.diff.front.n[sp]<-length(time.diff.front1)
  time.diff.behind.n[sp]<-length(time.diff.behind1)
}
q<-data.frame(SpGroup=SpTable,TimeDiffFront=time.diff.front,TimeDiffBehind=time.diff.behind,TimeDiffFront.n<-time.diff.front.n,TimeDiffBehind.n=time.diff.behind.n)
write.table(q,file="C:/Dave/ABMI/Cameras/2018 analysis/Distance models/Time between images in front of and behind pole.csv",sep=",",row.names=FALSE)

# Detection distance models

library(mgcv)  # Note: Not using random effects for site, because initial test showed little effect on coefficients

v.l<-read.csv("C:/Dave/ABMI/Cameras/2015 analysis/Distance models/Lookup for pole distance model predictions.csv")  # Lookup table for broader vegHF1, vegHF2, etc types for each vegHF type

# Only use records where d$n==d$Front | d$n==d$Behind, because if there are mixed distances for a group, we can't know which animal(s) triggered or would have triggered the camera
d<-d[d$n==d$Front | d$n==d$Behind,]
d<-d[!is.na(d$Deployment),]  # Bunch of completely empty records generated, don't know why
d$Season<-ifelse(d$Julian>=j.summer.start & d$Julian<=j.summer.end,"Summer","Winter")
d$Season<-as.factor(d$Season)

for (sp in 1:length(SpTable)) {

  print(paste(sp,length(SpTable),SpTable[sp],date()))

  d.sp<-d[d$SpGroup==SpTable[sp],]

  # Number of informative trials for that record CHANGE THIS TO SIGN - Doesn't matter how many animals - only one trips the camera.
  d.sp$n<-d.sp$Front+d.sp$Behind

  d.sp$pBehind<-d.sp$Behind/d.sp$n

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
  dist.veg<-5/sqrt(1-plogis(p1.veg))  # Time-between-images effect would be added in here, but no meaningful differences found
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

# Lure effect
# Do this using ABMI sites where lure and unlured are matched.  There are close to equal numbers of lured and unlured deployments
d<-read.csv(paste0(abmisc,"processed/dave-objects/Mammals by series Feb 2019 ALL.csv"))
d<-d[substr(d$Deployment,1,4)=="ABMI" & d$Year>=2015,]  # Don't use off-grids, or non-ABMI projects - unlured only mostly, or all lured for WHEC
d<-d[substr(d$Deployment,1,6)!="ABMI-W",]  # And don't use wetlands for same reason
d<-d[!is.na(d$Lure),]
z<-by(d$Lure,paste(d$Deployment,d$Year,sep="_"),function(x) x[1])  # Check number of lured and unlured per year
Year<-substr(names(z),nchar(names(z))-3,nchar(names(z)))
x<-table(Year,as.character(z))
# Add deployments with no native mammals - t is metadata file processed above
t$Deployment<-as.character(t$Deployment)
t1<-t[substr(t$Deployment,1,4)=="ABMI" & t$Year>=2015,]
t1<-t1[substr(t1$Deployment,1,6)!="ABMI-W",]
t2<-t1[is.na(match(t1$DeploymentYear,d$DeploymentYear)),]  # Qualifying deployments with no native mammals
t2$Lure<-as.character(t2$Lure)
x1<-table(t2$Year,t2$Lure)
x<-x+x1  # Add to totals by year and lure status
# Then summarize for each year
d1<-d[d$Year==2015,]
q.u<-by(d1$total_time_tbp[d1$Lure=="n"],d1$common_name[d1$Lure=="n"],sum)/x[1,1]  # Mean total_time_tbp per site
q.l<-by(d1$total_time_tbp[d1$Lure=="y"],d1$common_name[d1$Lure=="y"],sum)/x[1,2]
n.u<-by(d1$total_time_tbp[d1$Lure=="n"],d1$common_name[d1$Lure=="n"],length)/x[1,1]  # Mean records per site
n.l<-by(d1$total_time_tbp[d1$Lure=="y"],d1$common_name[d1$Lure=="y"],length)/x[1,2]
q<-data.frame(Year=2015,common_name=names(q.l),total_time_tbp.UnLure=as.numeric(q.u),total_time_tbp.Lure=as.numeric(q.l),n.UnLure=as.numeric(n.u),n.Lure=as.numeric(n.l),Mean.total_time_tbp.UnLure=as.numeric(q.u)/as.numeric(n.u),Mean.total_time_tbp.Lure=as.numeric(q.l)/as.numeric(n.l))
d1<-d[d$Year==2016,]
q.u<-by(d1$total_time_tbp[d1$Lure=="n"],d1$common_name[d1$Lure=="n"],sum)/x[2,1]  # Mean total_time_tbp per site
q.l<-by(d1$total_time_tbp[d1$Lure=="y"],d1$common_name[d1$Lure=="y"],sum)/x[2,2]
n.u<-by(d1$total_time_tbp[d1$Lure=="n"],d1$common_name[d1$Lure=="n"],length)/x[2,1]  # Mean records per site
n.l<-by(d1$total_time_tbp[d1$Lure=="y"],d1$common_name[d1$Lure=="y"],length)/x[2,2]
q<-rbind(q,data.frame(Year=2016,common_name=names(q.l),total_time_tbp.UnLure=as.numeric(q.u),total_time_tbp.Lure=as.numeric(q.l),n.UnLure=as.numeric(n.u),n.Lure=as.numeric(n.l),Mean.total_time_tbp.UnLure=as.numeric(q.u)/as.numeric(n.u),Mean.total_time_tbp.Lure=as.numeric(q.l)/as.numeric(n.l)))
d1<-d[d$Year==2017,]
q.u<-by(d1$total_time_tbp[d1$Lure=="n"],d1$common_name[d1$Lure=="n"],sum)/x[3,1]  # Mean total_time_tbp per site
q.l<-by(d1$total_time_tbp[d1$Lure=="y"],d1$common_name[d1$Lure=="y"],sum)/x[3,2]
n.u<-by(d1$total_time_tbp[d1$Lure=="n"],d1$common_name[d1$Lure=="n"],length)/x[3,1]  # Mean records per site
n.l<-by(d1$total_time_tbp[d1$Lure=="y"],d1$common_name[d1$Lure=="y"],length)/x[3,2]
q<-rbind(q,data.frame(Year=2017,common_name=names(q.l),total_time_tbp.UnLure=as.numeric(q.u),total_time_tbp.Lure=as.numeric(q.l),n.UnLure=as.numeric(n.u),n.Lure=as.numeric(n.l),Mean.total_time_tbp.UnLure=as.numeric(q.u)/as.numeric(n.u),Mean.total_time_tbp.Lure=as.numeric(q.l)/as.numeric(n.l)))
d1<-d[d$Year==2018,]
q.u<-by(d1$total_time_tbp[d1$Lure=="n"],d1$common_name[d1$Lure=="n"],sum)/x[4,1]  # Mean total_time_tbp per site
q.l<-by(d1$total_time_tbp[d1$Lure=="y"],d1$common_name[d1$Lure=="y"],sum)/x[4,2]
n.u<-by(d1$total_time_tbp[d1$Lure=="n"],d1$common_name[d1$Lure=="n"],length)/x[4,1]  # Mean records per site
n.l<-by(d1$total_time_tbp[d1$Lure=="y"],d1$common_name[d1$Lure=="y"],length)/x[4,2]
q<-rbind(q,data.frame(Year=2018,common_name=names(q.l),total_time_tbp.UnLure=as.numeric(q.u),total_time_tbp.Lure=as.numeric(q.l),n.UnLure=as.numeric(n.u),n.Lure=as.numeric(n.l),Mean.total_time_tbp.UnLure=as.numeric(q.u)/as.numeric(n.u),Mean.total_time_tbp.Lure=as.numeric(q.l)/as.numeric(n.l)))
write.table(q, file=paste0(abmisc,"processed/dave-objects/Series summary unlured versus lured Feb 2019 ALL.csv"), sep = ",", row.names = FALSE)

# And bootstrapped version to check for changes in lure effect in different years
q<-by(d$n_photos,d$common_name,sum)  # Figure out top species
sp.top<-names(sort(-q))[1:23]
sp.top<-sp.top[-which(sp.top=="Bighorn sheep")]  # Only one year
sp.top<-sp.top[-which(sp.top=="Richardson's Ground Squirrel")]
sp.top<-sp.top[-which(sp.top=="Bison")]  # Only one year
d2<-d[d$common_name%in% sp.top,]  # Limit to top species only
d2$Deployment<-as.character(d2$Deployment)
d2$common_name<-as.factor(as.character(d2$common_name))  # And make sure that excluded species are excluded
d2<-d2[substr(d2$Deployment,1,4)=="ABMI" & d2$Year>=2015,]  # Only use ABMI on-grids with paired design
d2<-d2[substr(d2$Deployment,1,6)!="ABMI-W",]  # And don't use wetlands for same reason
d2<-d2[!is.na(d2$Lure),]
d2$Deployment<-as.character(d2$Deployment)  # This ensures that excluded deployments are not included
d2$Site<-substr(d2$Deployment,6,nchar(d2$Deployment)-3)  # This assumes all sites are ABMI-[Site]-{NW,NE,SW,SE}
# Add qualifying deployments with no native mammals (in sp.top)
t2<-t1[is.na(match(t1$DeploymentYear,d2$DeploymentYear)),]   # t1 prepared above
t3<-data.frame(Deployment=t2$Deployment,Year=t2$Year,Lure=t2$Lure,Time1=t2$StartTime,SeriesNum=0,n_photos=0,common_name=NA,sex=NA,age_class=NA,multiple_animals=NA,mean_animals=0,total_time_tbp=0,DeploymentYear=t2$DeploymentYear,Site=NA)
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
    q.u<-by(d1.bs$total_time_tbp[d1.bs$Lure=="n"],d1.bs$common_name[d1.bs$Lure=="n"],sum)/length(unique(d1.bs$Deployment[d1.bs$Lure=="n"]))  # Average total_time_tbp per site
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
png(file="S:/github-repos-data/SC-Camera-Mammals/results/figures/Lure effect over years.png",height=600,width=900)
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

# Section below has not been updated in 2018
# Camera year calibration (2017 and 2018 only)  NOTE: Because mostly 2 cameras at 2018 sites - one lured, one not - there are few additional comparisons provided by 2018 data
d<-read.csv("C:/Dave/ABMI/Cameras/2018 analysis/Mammals by series 2015 to 2018 ALL.csv")
d$Deployment<-as.character(d$Deployment)  # This ensures that excluded deployments are not included
d$DeploymentYear<-paste(d$Deployment,d$Year,sep="_")
d<-d[substr(d$Deployment,1,4)=="ABMI",]  # Don't use off-grids or non-ABMI - no paired design
d<-d[substr(d$Deployment,1,6)!="ABMI-W",]  # Don't use wetlands - no paired design
d<-d[substr(d$Deployment,nchar(d$Deployment),nchar(d$Deployment))!="b",]  # Don't use b sites - not paired
d<-d[d$Year==2017 | d$Year==2018,]  # Paired comparison only in 2017 and 2018
c1<-read.csv("C:/Dave/ABMI/Data/Mammals/2018/Camera deployments with year camera bought ALL.csv")  # Assembled in Excel
c1$DeploymentYear<-paste(c1$Deployment,c1$Year,sep="_")
c1<-c1[duplicated(c1$DeploymentYear)==FALSE,]  # Sometimes 2 or more record per deployment (when there was a break in the middle)
d2<-merge(d,c1[,c("DeploymentYear","YearBought")],by="DeploymentYear")
d2<-d2[d2$Lured=="n",]  # Pairing was only done somewhat rigorously for unlured sites
d2$YearBought<-as.factor(as.character(d2$YearBought))  # This ensures that all years show up for each species below, even if the species was not detected in qualifying deployments one year
d2$Site<-substr(d2$Deployment,6,nchar(d2$Deployment)-3)  # This assumes all sites are ABMI-[Site]-{NW,NE,SW,SE}
# Find top ten species (number of photos)
q<-by(d2$Photos,d2$Species,sum)
sp.top<-names(sort(-q))[1:19]
sp.top<-sp.top[-which(sp.top=="Woodland Caribou")]  # Exclude, because none in 2013 cameras (due to small sample, presumably)
# Figure out which sites have pairs to be compared (this doesn't account for sites with no native mammals at all)
q<-table(d2$Site,d2$YearBought)
q<-sign(q[rowSums(sign(q))==2,])
sites.13.14<-rownames(q)[which(q[,"2013"]==1 & q[,"2014"]==1)]  # Sites at which there is a paired comparison of 2013 and 2014 bought cameras
sites.14.16<-rownames(q)[which(q[,"2016"]==1 & q[,"2014"]==1)]  # Sites at which there is a paired comparison of 2014 and 2016 bought cameras
niter<-1000
# 2013 versus 2014
d3<-d2[d2$Site %in% sites.13.14,]
d3$YearBought<-as.factor(as.character(d3$YearBought))  # This ensures that all years show up for each species below, even if the species was not detected in qualifying deployments one year
pc.13.14<-array(NA,c(length(sp.top),niter))
for (iter in 1:niter) {
  if ((iter-1)/100==floor((iter-1)/100)) print(paste(iter,niter,date()))
  s<-sample(1:length(sites.13.14),length(sites.13.14),replace=TRUE)
  i<-unlist(lapply(sites.13.14[s],function(a) which(d3$Site %in% a)))
  d3.bs<-d3[i,]  # Bootstrap resample of sites with that comparison
  for (sp in 1:length(sp.top)) {
    d4<-d3.bs[d3.bs$Sp==sp.top[sp],]
    if (nrow(d4)>0) x<-by(d4$Photos,d4$YearBought,sum)/length(sites.13.14) # Average number of images per deployment by YearBought (including 0's)
    x<-ifelse(is.na(x),0,x)
    pc.13.14[sp,iter]<-x[1]/sum(x)*100  # Percentage of photos that are from the first YearBought
  }  # Next sp
}  # Next iter
for (sp in 1:length(sp.top)) {
  if (sp==1) {
    qq.13.14<-data.frame(Sp=sp.top[sp],t(quantile(pc.13.14[sp,],c(0.025,0.05,0.1,0.5,0.9,0.95,0.975))))
  } else {
    qq.13.14<-rbind(qq.13.14,data.frame(Sp=sp.top[sp],t(quantile(pc.13.14[sp,],c(0.025,0.05,0.1,0.5,0.9,0.95,0.975)))))
  }
}
names(qq.13.14)<-c("Sp","pc.13.vs.14.q2.5","pc.13.vs.14.q5","pc.13.vs.14.q10","pc.13.vs.14.median","pc.13.vs.14.q90","pc.13.vs.14.q95","pc.13.vs.14.q97.5")
# 2014 versus 2016
d3<-d2[d2$Site %in% sites.14.16,]
d3$YearBought<-as.factor(as.character(d3$YearBought))  # This ensures that all years show up for each species below, even if the species was not detected in qualifying deployments one year
pc.14.16<-array(NA,c(length(sp.top),niter))
for (iter in 1:niter) {
  if ((iter-1)/100==floor((iter-1)/100)) print(paste(iter,niter,date()))
  s<-sample(1:length(sites.14.16),length(sites.14.16),replace=TRUE)
  i<-unlist(lapply(sites.14.16[s],function(a) which(d3$Site %in% a)))
  d3.bs<-d3[i,]  # Bootstrap resample of sites with that comparison
  for (sp in 1:length(sp.top)) {
    d4<-d3.bs[d3.bs$Sp==sp.top[sp],]
    if (nrow(d4)>0) x<-by(d4$Photos,d4$YearBought,sum)/length(sites.14.16) # Average number of images per deployment by YearBought (including 0's)
    x<-ifelse(is.na(x),0,x)
    pc.14.16[sp,iter]<-x[1]/sum(x)*100  # Percentage of photos that are from the first YearBought
  }  # Next sp
}  # Next iter
for (sp in 1:length(sp.top)) {
  if (sp==1) {
    qq.14.16<-data.frame(Sp=sp.top[sp],t(quantile(pc.14.16[sp,],c(0.025,0.05,0.1,0.5,0.9,0.95,0.975))))
  } else {
    qq.14.16<-rbind(qq.14.16,data.frame(Sp=sp.top[sp],t(quantile(pc.14.16[sp,],c(0.025,0.05,0.1,0.5,0.9,0.95,0.975)))))
  }
}
names(qq.14.16)<-c("Sp","pc.14.vs.16.q2.5","pc.14.vs.16.q5","pc.14.vs.16.q10","pc.14.vs.16.median","pc.14.vs.16.q90","pc.14.vs.16.q95","pc.14.vs.16.q97.5")
save(file="C:/Dave/ABMI/Cameras/2018 analysis/Calibration/R objects Calibration results for camera bought year 2018.rData",qq.14.16,qq.13.14,sites.14.16,sites.13.14)
# Figures
qq.13.14<-qq.13.14[qq.13.14$Sp!="Deer",]
qq.14.16<-qq.14.16[qq.14.16$Sp!="Deer",]
sp.order<-c(7,12,6,13,14,11,3,4,2,1,5,8,9,10)  # Bears, dogs, cats, ungulates, others
x<-c(1,2.5,3.5,4.5,6,7.5,9,10,11,12,13,14.5,15.5,16.5)
col1<-rainbow(15,v=0.8)
png(file="C:/Dave/ABMI/Cameras/2018 analysis/Calibration/Camera bought-year effect 2018.png",height=750,width=750)
par(mai=c(3,0.9,0.3,0.3))
dplot(0,0,xlim=range(x),ylim=log(c(0.01,100)),xlab="",ylab="Model Year1 : Model Year2",xaxt="n",yaxt="n",typ="n")
axis(side=2,at=log(c(0.01,0.02,0.1,0.2,0.25,1/3,1/2,1,2,3,4,5,7,10,20,50,100)),lab=rep("",17),tck=1,col="grey80")
axis(side=2,at=log(c(0.01,0.02,0.1,0.2,0.25,1/3,1/2,1,2,3,4,5,7,10,20,50,100)),lab=c(0.01,0.02,0.1,0.2,0.25,1/3,1/2,1,2,3,4,5,7,10,20,50,100),tck=0.015,cex.lab=1.3)
for (i in 1:length(x)) {
  y<-as.numeric((qq.13.14[sp.order[i],2:7]/100)/(1-qq.13.14[sp.order[i],2:7]/100))  # Convert percent of time in first year to first:second year ratio
  y<-ifelse(y==Inf,100,y)
  y<-ifelse(y==0,0.01,y)
  points(x[i]-0.2,log(y[4]),pch=18,cex=1.7,typ="p",lwd=2,col=col1[i])
  lines(rep(x[i]-0.2,2),log(y[c(2,6)]),col=col1[i])
  y<-as.numeric((qq.14.16[sp.order[i],2:7]/100)/(1-qq.14.16[sp.order[i],2:7]/100))  # Convert percent of time in first year to first:second year ratio
  y<-ifelse(y==Inf,100,y)
  y<-ifelse(y==0,0.01,y)
  points(x[i]+0.2,log(y[4]),pch=18,cex=1.7,typ="p",lwd=2,col=col1[i])
  lines(rep(x[i]+0.2,2),log(y[c(2,6)]),col=col1[i])
  mtext(side=1,at=x[i],qq.14.16$Sp[sp.order[i]],las=2,cex=1.3,adj=1,line=0.5,col=col1[i])
}
graphics.off()












