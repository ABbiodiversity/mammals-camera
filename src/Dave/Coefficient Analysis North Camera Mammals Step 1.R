# R ABMI coefficients analysis for mammals in the north, using camera results
# This is step 1, process data files to be used by step 2 (simple run with diagnostic figures).
# Step 3, a bootstrapped version, has not been updated for the camera mammals, due to lack of any interest in those results.
# The data file has one row per camera, with summaries of the species' photos, series and duration, as well as seasonal times camera was operating and detection distances per species group.  Pre-processed in Access and R.
# Also, processing of km2 grid for the province is done separately, mainly because of memory limits

datafile<-"C:/Dave/ABMI/Cameras/2018 analysis/Camera mammal data processed Feb 2019.csv"  # Processed file from Access, preliminary R scrips, etc.
vegfile<-"C:/Dave/ABMI/Data/Site info/2018/Combined vegHF soilHF and detection distance veg for cameras Feb 2019.csv"  #  File simplified and corrected from GIS raw point summary, and information added from images for missed deployments
HFgroupfile<-"C:/Dave/ABMI/Data/Site info/2018/lookup-hf-class.csv"  # Lookup table for HF to HF group from Peter - modified for what I think are this year's HF categories
pmfile<-"C:/Dave/ABMI/Cameras/Coefficients/2018/Analysis north/Prediction matrix for mammal models NORTH 2018.csv"  # Prediction matrix - shows which (grouped) veg types and HF to use to predict which coefficients.  Including for P/A A|P models. Needs to be edited in Excel when veg or HF groupings change in models.
siteinfofile<-"C:/Dave/ABMI/Data/Site info/Site summary with climate.csv"
deploymentlocationsfile<-"C:/Dave/ABMI/Data/Mammals/2018/Deployment locations all Feb 2019.csv"  # Needed to fill in missing nearest sites below
dataset.out<-"C:/Dave/ABMI/Cameras/Coefficients/2018/Analysis north/R Dataset SpTable for ABMI North mammal coefficients 2018.RData"  # To save combined species-veg-HF file, plus SpTable.

# Read data file
d1<-read.csv(datafile)
d1$Deployment<-as.character(d1$Deployment)
# Get rid of unnecessary columns
j<-which(regexpr(".Duration",names(d1))>0)
d1<-d1[,-j]
j<-which(regexpr("DetDist",names(d1))>0)
d1<-d1[,-j]
# And eliminate studies that don't use ABMI protocols or are otherwise weird
study.list<-c("ABMI","ABMI-W","OG-ABMI","CMU","BigGrid","OG-OGC","OG-EI")  # Qualifying studies.  EI study further edited below, to get rid of roadside deployments.  Need to change this and lines below if additional studies are added, or site naming changes
d1$Study<-"ABMI"
d1$Study<-ifelse(substr(d1$Deployment,1,6)=="ABMI-W","ABMI-W",d1$Study)
d1$Study<-ifelse(substr(d1$Deployment,1,4)=="WHEC","WHEC",d1$Study)
d1$Study<-ifelse(substr(d1$Deployment,1,3) %in% c("CAL","CHR","DEW","FMM","LLB","LRN","MAC","MCC","WAB") == TRUE,"CMU",d1$Study)
d1$Study<-ifelse(substr(d1$Deployment,1,2)=="BG","BigGrid",d1$Study)
d1$Study<-ifelse(substr(d1$Deployment,1,3)=="NEX","Nexen",d1$Study)
d1$Study<-ifelse(substr(d1$Deployment,1,6)=="OG-AAC","OG-AAC",d1$Study)
d1$Study<-ifelse(substr(d1$Deployment,1,6)=="OG-ABM","OG-ABMI",d1$Study)
d1$Study<-ifelse(substr(d1$Deployment,1,9)=="OG-CITSCI","OG-CITSCI",d1$Study)
d1$Study<-ifelse(substr(d1$Deployment,1,9)=="OGW-CITSC","OGW-CITSCI",d1$Study)
d1$Study<-ifelse(substr(d1$Deployment,1,5)=="OG-EI","OG-EI",d1$Study)
d1$Study<-ifelse(substr(d1$Deployment,1,6)=="OG-OGC","OG-OGC",d1$Study)
d1$Study<-ifelse(substr(d1$Deployment,1,7)=="OG-RIVR","OG-RIVR",d1$Study)
d1<-d1[d1$Study %in% study.list,]
i<-which(d1$Study=="OG-EI" & substr(d1$Deployment,nchar(d1$Deployment),nchar(d1$Deployment))=="1")  # Also get rid of EI edge locations
d1<-d1[-i,]

# Read in Veg + HF, process
v<-read.csv(vegfile)
v$Deployment<-as.character(v$deployment)
v$deployment<-NULL
v$VegHF<-as.character(v$VegHF)  # This column includes corrections from comparing the GIS info to the images
v$VegHF<-ifelse(substr(v$Deployment,1,6)=="ABMI-W","WetlandMargin",v$VegHF)  # Set all ABMI-W sites to WetlandMargin - specific stratum targeted by that study's design
# Set up columns, 1 per veg or HF type
VegHF.list<-sort(unique(v$VegHF))
q<-array(0,c(nrow(v),length(VegHF.list)))
colnames(q)<-VegHF.list
v<-cbind(v,q)
for (i in VegHF.list) v[,i]<-ifelse(v$VegHF==i,1,0)  # Assign proportion of area = 1 to appropriate veg or HF column for each site
# Group HF types
HFgl<-read.csv(HFgroupfile)  # HF group lookup
HFgroups<-unique(HFgl$UseInAnalysis)
for (i in HFgroups) {
  HFtypesingroup<-colnames(v)[na.omit(match(HFgl$HF_GROUP[HFgl$UseInAnalysis==i],colnames(v)))]
  if (i=="HFor") HFtypesingroup<-c(HFtypesingroup,names(v)[substr(names(v),1,2)=="CC"])  # Add the individual cutblock types to HFor group
  if (length(HFtypesingroup)>1) {
    v<-cbind(v,rowSums(v[,HFtypesingroup])) # Add a column for each HF group, summing the component HF types
  } else {
    v<-cbind(v,v[,HFtypesingroup])  # If there is only one HF type in that HF group, just add it
  }
  if (length(HFtypesingroup)==0) HFgroups<-HFgroups[-which(HFgroups==i)]  # Remove that group from the list if there are one of that HF type in the camera sample
}
names(v)<-c(names(v)[1:(ncol(v)-length(HFgroups))],as.character(HFgroups))
# And broader groups - SuccAlien
HFgroups<-unique(HFgl$SuccAlien)
for (i in HFgroups) {
  HFtypesingroup<-colnames(v)[na.omit(match(HFgl$HF_GROUP[HFgl$SuccAlien==i],colnames(v)))]
  if (i=="Succ") HFtypesingroup<-c(HFtypesingroup,names(v)[substr(names(v),1,2)=="CC"])  # Add the individual cutblock types to Succ group
  if (length(HFtypesingroup)>1) {
    v<-cbind(v,rowSums(v[,HFtypesingroup])) # Add a column for each HF group, summing the component HF types
  } else {
    v<-cbind(v,v[,HFtypesingroup])  # If there is only one HF type in that HF group, just add it
  }
  if (length(HFtypesingroup)==0) HFgroups<-HFgroups[-which(HFgroups==i)]  # Remove that group from the list if there are one of that HF type in the camera sample
}
names(v)<-c(names(v)[1:(ncol(v)-length(HFgroups))],as.character(HFgroups))
# And broader groups - NonlinLin
HFgroups<-unique(HFgl$NonlinLin)
for (i in HFgroups) {
  HFtypesingroup<-colnames(v)[na.omit(match(HFgl$HF_GROUP[HFgl$NonlinLin==i],colnames(v)))]
  if (i=="Nonlin") HFtypesingroup<-c(HFtypesingroup,names(v)[substr(names(v),1,2)=="CC"])  # Add the individual cutblock types to Nonlin group
  if (length(HFtypesingroup)>1) {
    v<-cbind(v,rowSums(v[,HFtypesingroup])) # Add a column for each HF group, summing the component HF types
  } else {
    v<-cbind(v,v[,HFtypesingroup])  # If there is only one HF type in that HF group, just add it
  }
}
names(v)<-c(names(v)[1:(ncol(v)-length(HFgroups))],as.character(HFgroups))
# CHECK for na in (any) veg class - make sure these deployments make sense
# Also check that the max of broad HF classes is never more than 1 - if so, there is some problem in the lookup (often a duplicated name)

# Add lat long, natural regions, etc. to sites
# First add in nearest.sites for non-ABMI studies that do not have that info in the deployment name
s1<-read.csv(deploymentlocationsfile)  # Summarized from larger meta-data file
s1$Lat<-as.numeric(as.character(s1$Public.Latitude))
s1$Long<-as.numeric(as.character(s1$Public.Longitude))
s1$DeploymentYear<-paste(s1$Site.Name,s1$Year,sep="_")
d1<-merge(d1,s1[,c("DeploymentYear","NearestSite","Long","Lat")])  # Check for loss of info
d1$NearestSite<-as.numeric(as.character(d1$NearestSite))  # Warning message about NA's is okay - trying to change some "" records for NA's here
s<-read.csv(siteinfofile)  # To get ABMI site lat longs, and then natural regions via nearest ABMI sites
names(s)[which(names(s)=="PUBLIC_LONGITUDE")]<-"Long"
names(s)[which(names(s)=="PUBLIC_LATTITUDE")]<-"Lat"
names(s)[which(names(s)=="NATURAL_REGIONS")]<-"NR"
names(s)[which(names(s)=="NATURAL_SUBREGIONS")]<-"NSR"
names(s)[which(names(s)=="LANDUSE_FRAMEWORK")]<-"LUF"
for (i in 1:nrow(d1)) {  # Fill in missing nearest ABMI sites
  if (is.na(d1$NearestSite[i]) & !is.na(d1$Lat[i])) d1$NearestSite[i]<-s$SITE_ID[which.min(sqrt((d1$Long[i]-s$Long)^2 + (d1$Lat[i]-s$Lat)^2))]
}
q<-merge(d1,s[,c("SITE_ID","NR","NSR","LUF","AHM","PET","FFP","MAP","MAT","MCMT","MWMT")],by.x="NearestSite",by.y="SITE_ID")  # Check for loss of sites
q$TrueLat<-q$Lat
q$Lat<-ifelse(q$Lat<51.5,51.5,q$Lat)  # Truncated latitude for spatial modeling
# Then merge in veg
v$DeploymentYear<-as.character(v$DeploymentYear)
q$DeploymentYear<-as.character(q$DeploymentYear)
q0<-merge(q,v[,-(1:4)],by="DeploymentYear")  # Check for loss of sites.
d<-q0
names(d)[which(names(d)=="Deployment.x")]<-"Deployment"
d$Deployment.y<-NULL

# Combine indistinguishable wetland classes and group age classes
# Note: All age class 0 have been allocated an age class by lookng at the images, unless there were no images at that deployment (in which case, they should not be in d any more anyway)
i<-which(substr(names(d),1,8)=="TreedBog")
d$TreedBog<-rowSums(d[,i])
i<-which(substr(names(d),1,8)=="TreedFen")
d$TreedFen<-rowSums(d[,i])
d$ShrubbyBogFen<-d$ShrubbyBog+d$ShrubbyFen
i<-which(substr(names(d),1,10)=="TreedSwamp")
d$TreedSwamp<-rowSums(d[,i])
for (i1 in i[order(-i)]) d[,i1]<-NULL # Get rid of these so they are not used by mistake
i<-which(substr(names(d),1,7)=="CCDecid" | substr(names(d),1,7)=="CCMixed")
d$CCDecidMixed<-rowSums(d[,i])
i<-which(substr(names(d),1,8)=="CCSpruce")
d$CCSpruce<-rowSums(d[,i])
i<-which(substr(names(d),1,6)=="CCPine")
d$CCPine<-rowSums(d[,i])
d$CCConif<-d$CCSpruce+d$CCPine
d$Marsh<-d$Marsh+d$GraminoidFen  # Change these all to marsh  Note: Could also be GrassyBog, if this isn't converted to another type in the original data file
d$GraminoidFen<-rep(0,nrow(d))  # Set to 0 so not two types for those sites  (and also GrassyBog if that is included in previous line)
# Truncate age distributions at class 8 (so this class becomes >140yr), because no ABMI data in older pine, deciduous, or mixed, and marginal in other upland and wet conifers.  Also CC's
# Check for other stand types with age 9 data in future iterations
d$Spruce8<-d$Spruce8+d$Spruce9
d$Spruce9<-NULL  # Get rid of those old classes to avoid confusion
d$Decid8<-d$Decid8+d$Decid9
d$Decid9<-NULL  # Get rid of those old classes to avoid confusion
d$Pine8<-d$Pine8+d$Pine9
d$Pine9<-NULL  # Get rid of those old classes to avoid confusion
d$Mixedwood8<-d$Mixedwood8+d$Mixedwood9
d$Mixedwood9<-NULL  # Get rid of those old classes to avoid confusion
d$TreedFen8<-d$TreedFen8+d$TreedFen9
d$TreedFen9<-NULL  # Get rid of those old classes to avoid confusion
d$TreedBog8<-d$TreedBog8+d$TreedBog9
d$TreedBog9<-NULL  # Get rid of those old classes to avoid confusion
d$CCPine1<-d$CCPine1+d$CCPine2+d$CCPine3  # Too few of class 2 and 3.  Check for other classes in updates
d$CCPine2<-NULL  # Get rid of those old classes to avoid confusion
d$CCPine3<-NULL  # Get rid of those old classes to avoid confusion
d$CCMixedwood2<-d$CCMixedwood2+d$CCMixedwood3+d$CCMixedwood4  # Too few of class 3 and 4. Check for other classes in updates
d$CCMixedwood3<-NULL  # Get rid of those old classes to avoid confusion
d$CCMixedwood4<-NULL  # Get rid of those old classes to avoid confusion
d$CCSpruce2<-d$CCSpruce2+d$CCSpruce3+d$CCSpruce4  # Too few of class 3 and 4. Check for other classes in updates
d$CCSpruce3<-NULL  # Get rid of those old classes to avoid confusion
d$CCSpruce4<-NULL  # Get rid of those old classes to avoid confusion
# Make combined classes for entire stand types
d$Spruce<-d$SpruceR+d$Spruce1+d$Spruce2+d$Spruce3+d$Spruce4+d$Spruce5+d$Spruce6+d$Spruce7+d$Spruce8
d$Pine<-d$PineR+d$Pine1+d$Pine2+d$Pine3+d$Pine4+d$Pine5+d$Pine6+d$Pine7+d$Pine8  # Check if all these are in input data - causes error if some type is not found in input file - and check that there are no new ones in future iterations
d$Decid<-d$DecidR+d$Decid1+d$Decid2+d$Decid3+d$Decid4+d$Decid5+d$Decid6+d$Decid7+d$Decid8  # Check if all these are in input data - causes error if some type is not found in input file
d$Mixedwood<-d$MixedwoodR+d$Mixedwood2+d$Mixedwood3+d$Mixedwood4+d$Mixedwood5+d$Mixedwood6+d$Mixedwood7+d$Mixedwood8
d$TreedBogFen<-d$TreedFen+d$TreedBogR+d$TreedBog1+d$TreedBog2+d$TreedBog3+d$TreedBog4+d$TreedBog5+d$TreedBog6+d$TreedBog7+d$TreedBog8

# Add combined veg or HF variables for more general models
d$DecidMixed<-d$Decid+d$Mixedwood
d$UpCon<-d$Spruce+d$Pine
d$CCAll<-d$CCDecidMixed+d$CCSpruce
d$Cult<-d$Crop+d$RoughPasture+d$TamePasture
d$ShrubbyWet<-d$ShrubbyBogFen+d$ShrubbySwamp
d$OpenWet<-d$ShrubbyWet+d$Marsh
d$TreedWet<-d$TreedBogFen+d$TreedSwamp
d$GrassShrub<-d$GrassHerb+d$Shrub
d$TreedAll<-d$DecidMixed+d$UpCon+d$TreedWet
d$OpenAll<-d$GrassShrub+d$OpenWet
d$Boreal<-d$TreedAll+d$OpenWet  # Everything except GrassShrub
d$THF<-d$Alien+d$Succ
d$Lowland<-d$TreedWet+d$OpenWet
d$UplandForest<-d$UpCon+d$DecidMixed  # Upland, except GrassShrub
d$Upland<-1-d$Lowland-d$THF-d$WetlandMargin
# CHECK for max of some of these groups being >1 - indicates duplication somewhere, which needs to be corrected

# Remove sites not in North (boreal or foothills)
d<-d[!is.na(d$Water),]  # Saskatchewan + a few oddballs
d$UseAsNorth<-ifelse(d$NR=="Boreal" | d$NR=="Canadian Shield" | d$NR=="Foothills" | d$NR=="Rocky Mountain" | d$NR=="Parkland","Y","N")
d$UseAsNorth<-ifelse(d$Water==1,"N",d$UseAsNorth)
d<-d[d$UseAsNorth=="Y",]

# Add weights for each record - because some sites sampled >1 time
# Eliminate summer records for any deployment with <10 days summer sampling, and ditto for winter.  And eliminate entirely any deployments that don't qualify in either season
d<-d[d$SummerDays>=10 | d$WinterDays>=10,]
i<-which(d$SummerDays<10)
j<-which(regexpr("Summer",names(d))>0)
j<-j[-1]  # First column is for Days - leave as is
d[i,j]<-NA  # Convert SummerDays to NA if SummerDays is <10 (and not already NA)
i<-which(d$WinterDays<10)
j<-which("Winter" %in% names(d)==TRUE)
j<-j[-1]  # First column is for Days - leave as is
d[i,j]<-NA
# Then weights based on number of qualifying visits (for each season)
q<-by(ifelse(d$SummerDays>=10,1,0),d$Deployment,sum)
d$wt.s<-1/as.numeric(q[match(d$Deployment,names(q))])
d$wt.s<-ifelse(d$wt.s==Inf,0,d$wt.s)  # no weight to deployments that have never been sampled (10+ days) in summer
q<-by(ifelse(d$WinterDays>=10,1,0),d$Deployment,sum)
d$wt.w<-1/as.numeric(q[match(d$Deployment,names(q))])
d$wt.w<-ifelse(d$wt.w==Inf,0,d$wt.w)  # no weight to deployments that have never been sampled (10+ days) in summer

# Set up list of species to analyse, separately for summer and winter
# Includes main analysis list as SpTable, and larger group to do use/availability for as SpTable.ua
FirstSpCol.s<-which(names(d)=="BadgerSummer")  # Find species names.  Check if species list changes
LastSpCol.s<-which(names(d)=="WoodlandCaribouSummer")
FirstSpCol.w<-which(names(d)=="BadgerWinter")  # Find species names
LastSpCol.w<-which(names(d)=="WoodlandCaribouWinter")
# Summer species table
SpTable.s<-SpTable.s.ua<-names(d)[FirstSpCol.s:LastSpCol.s] # All speciesXseasons
SpTable.s<-SpTable.s[regexpr("Summer",SpTable.s)>0]  # and Summer only
occ1.s<-NULL  # Figure out total number of occurrences of each species
for (i in 1:length(SpTable.s)) occ1.s[i]<-sum(sign(d[,SpTable.s[i]])*d$wt.s,na.rm=TRUE)
SpTable.s<-SpTable.s[-which(occ1.s<20)]  # Omit species with <20 occurrences
occ1.s<-NULL  # Figure out total number of occurrences of each species
for (i in 1:length(SpTable.s.ua)) occ1.s[i]<-sum(sign(d[,SpTable.s.ua[i]])*d$wt.s,na.rm=TRUE)  # Do not exclude parkland here,  - now used in use/availability in north (otherwise, no figures for veg ua of parkland species)
SpTable.s.ua<-SpTable.s.ua[-which(occ1.s<3)]  # Omit species with <3 occurrences for use/availability summaries
# and winter species table
SpTable.w<-SpTable.w.ua<-names(d)[FirstSpCol.w:LastSpCol.w] # All speciesXseasons
SpTable.w<-SpTable.w[regexpr("Winter",SpTable.w)>0]  # and Winter only
occ1.w<-NULL  # Figure out total number of occurrences of each species
for (i in 1:length(SpTable.w)) occ1.w[i]<-sum(sign(d[,SpTable.w[i]])*d$wt.w,na.rm=TRUE)
SpTable.w<-SpTable.w[-which(occ1.w<20)]  # Omit species with <20 occurrences
occ1.w<-NULL  # Figure out total number of occurrences of each species
for (i in 1:length(SpTable.w.ua)) occ1.w[i]<-sum(sign(d[,SpTable.w.ua[i]])*d$wt.w,na.rm=TRUE)  # Do not exclude parkland here,  - now used in use/availability in north (otherwise, no figures for veg ua of parkland species)
SpTable.w.ua<-SpTable.w.ua[-which(occ1.w<3)]  # Omit species with <3 occurrences for use/availability summaries

# Read prediction matrix into pm, so that it will be available in next step
pm<-read.csv(pmfile)

save(file=dataset.out,d,FirstSpCol.s,LastSpCol.s,FirstSpCol.w,LastSpCol.w,SpTable.s,SpTable.s.ua,SpTable.w,SpTable.w.ua,pm)


