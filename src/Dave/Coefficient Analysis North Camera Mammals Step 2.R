# R ABMI coefficients analysis for mammals in the north
# This version includes three options for spatial/climate residual mapping
# It also uses only the best veg+HF and age model for map predictions and figures
# Data has already been processed to density (/km2)
# Run step 1 first to set up the data files, and also km2 processing

library(mgcv)  # For binomial GAM
library(mapproj)  # For projected maps
library(binom)  # For exact binomial confidence intervals (not currently being used)
library(pROC)  # For AUC of ROC

# Set up file names for various outputs
fname.fig<-"C:/Dave/ABMI/Cameras/Coefficients/2018/Analysis north/Figures/Best model/"  # Start of figure file, one per species for veg+HF.  Subdirectory, species name and extension will be added
fname.map<-"C:/Dave/ABMI/Cameras/Coefficients/2018/Analysis north/Maps/"  # Start of map file, one per species for, spatial+climate map, current, reference and difference maps.  Subdirectory, species name and extension will be added.  Note: the web uses only whole-province maps, so naming here can be non-official
fname.sumout<-"C:/Dave/ABMI/Cameras/Coefficients/2018/Analysis north/Coefficient tables/Mammal coefficients North Feb 2019 Best model"  # For exporting coefficient tables
fname.Robjects<-"C:/Dave/ABMI/Cameras/Coefficients/2018/Analysis north/R objects/R objects North mammal coefficients Feb 2019 Best model"  # Start of file name to save each species' models
fname.Robjects.sum<-"C:/Dave/ABMI/Cameras/Coefficients/2018/Analysis north/R objects/R objects North mammals Coefficient tables Feb 2019 Best model.Rdata"
fname.km2summaries<-"C:/Dave/ABMI/Cameras/Coefficients/2018/Analysis north/Km2 summaries/Best model/Km2 North reference and current Feb 2019" # Start of file name for km2 grid reference and current output for each species
fname.useavail<-"C:/Dave/ABMI/Cameras/Coefficients/2018/Analysis north/Figures/Best model/Use availability/"  # Subdirectory for use-availability figures and table
km2file.out<-"C:/Dave/ABMI/Data/Km2 grid/2018/R km2 grid current and backfilled processed 2016 NORTH.Rdata"    # File with processed km2 grid files

# File name for previously processed datafile+veg+HF, km2 grid info, SpTable, and look-up matrix for veg groups
dataset.out<-"C:/Dave/ABMI/Cameras/Coefficients/2018/Analysis north/R Dataset SpTable for ABMI North mammal coefficients 2018.RData"
load(dataset.out)
SpTable.all<-sort(unique(c(gsub("Summer","",SpTable.s),gsub("Winter","",SpTable.w))))
# Correct mistakes in data found in analysis
i<-which(d$DeploymentYear=="ABMI-6-SW_2017")
d$SnowshoeHareSummer[i]<-d$SnowshoeHareSummer[i]+d$WhitetailedJackRabbitSummer[i]
d$WhitetailedJackRabbitSummer[i]<-0

# Specify space/climate model options for each species in SpTable.all  1=based on total abundance, 2=based on presence/absence only, 3=none.  Based on results of cross-validation at site and regional level.
sc.option<-c(1,2,3,2,2,2,3,3,3,3,1,1,3,3,2,2,2,3,2,2,3,3)  # These are the decisions from the previous year
sc.option<-ifelse(sc.option==3,2,sc.option)  # For this year, first try running all species with at least space/climate modeling of occurrence
data.frame(SpTable.all,sc.option)  # Check that these are in right order

taxa.file<-read.csv("C:/Dave/ABMI/Data/Mammals/2015/Mammal taxa.csv")  # To get proper species names
sp.names.s<-paste( as.character(taxa.file$Species[match(SpTable.s,paste(taxa.file$Label,"Summer",sep=""))]), "Summer")
sp.names11.s<-paste( as.character(taxa.file$Species[match(SpTable.s.ua,paste(taxa.file$Label,"Summer",sep=""))]), "Summer")
sp.names.w<-paste( as.character(taxa.file$Species[match(SpTable.w,paste(taxa.file$Label,"Winter",sep=""))]), "Winter")
sp.names11.w<-paste( as.character(taxa.file$Species[match(SpTable.w.ua,paste(taxa.file$Label,"Winter",sep=""))]), "Winter")
sp.names.all<-as.character(taxa.file$Species[match(SpTable.all,taxa.file$Label)])

# Figure parameters set up based on any possible variable in the veg+HF models
fp<-read.csv("C:/Dave/ABMI/Cameras/Coefficients/2018/Analysis north/Veg HF plot lookup.csv")  # Figure parameters, set up for all possible variables
fp$col1<-as.character(fp$col1)

load(km2file.out)
km2$Intercept<-1
km2.b$Intercept<-1

load("C:/Dave/ABMI/Data/Km2 grid/R object water dominated km2.rdata")  # Water-dominated km2 rasters to add to maps km2.water

# Set up km2 prediction matrices
# For climate+space effects only - for diagnostic mapping of these effects.  NSR groupings are being used experimentally
km2.sc<-data.frame(km2[,c("Intercept","Lat","Long","TrueLat","AHM","PET","FFP","MAP","MAT","MCMT","MWMT","NSR1Parkland","NSR1DryMixedwood","NSR1CentralMixedwood","NSR1Foothills","NSR1North","NSR1Shield","NSR1Mountain")],
                   Lat2=km2$Lat^2,Lat3=km2$Lat^3,Long2=km2$Long^2,LatLong=km2$Lat*km2$Long,Lat2Long2=km2$Lat^2*km2$Long^2,LongMAT=km2$Long*km2$MAT,MAPPET=km2$MAP*km2$PET,MATAHM=km2$MAT*km2$AHM,MAPFFP=km2$MAP*km2$FFP,MAT2=km2$MAT*(km2$MAT+10),MWMT2=km2$MWMT^2)
vnames.b<-c("SpruceR","Spruce1","Spruce2","Spruce3","Spruce4","Spruce5","Spruce6","Spruce7","Spruce8","PineR","Pine1","Pine2","Pine3","Pine4","Pine5","Pine6","Pine7","Pine8","DecidR","Decid1","Decid2","Decid3","Decid4","Decid5","Decid6","Decid7","Decid8",
            "MixedwoodR","Mixedwood1","Mixedwood2","Mixedwood3","Mixedwood4","Mixedwood5","Mixedwood6","Mixedwood7","Mixedwood8","TreedBogR","TreedBog1","TreedBog2","TreedBog3","TreedBog4","TreedBog5","TreedBog6","TreedBog7","TreedBog8",
            "TreedFen","TreedSwamp","Shrub","GrassHerb","GraminoidFen","Marsh","ShrubbyBog","ShrubbyFen","ShrubbySwamp","WetlandMargin")  # Variables to use in reference predictions
DF.res.null<-data.frame(Intercept=1,Lat=mean(d$Lat),Long=mean(d$Long),AHM=mean(d$AHM),PET=mean(d$PET),FFP=mean(d$FFP),MAP=mean(d$MAP),MAT=mean(d$MAT),MCMT=mean(d$MCMT),MWMT=mean(d$MWMT),
                        Lat2=mean(d$Lat)^2,Long=mean(d$Long)^2,LatLong=mean(d$Lat)*mean(d$Long),MAPPET=mean(d$MAP)*mean(d$PET),MATAHM=mean(d$MAT)*mean(d$AHM),MAPFFP=mean(d$MAP)*mean(d$FFP),MAT2=mean(d$MAT*(d$MAT+10)),MWMT2=mean(d$MWMT)^2,
                        NSR1Parkland=1/7,NSR1DryMixedwood=1/7,NSR1CentralMixedwood=1/7,NSR1Foothills=1/7,NSR1North=1/7,NSR1Shield=1/7,NSR1Mountain=1/7)  # Means, used to plot residual relationships

# Set up information for mapping
city.y<-c(51,53,49,56,58,55)+c(3,33,42,44,31,10)/60
city.x<- -c(114,113,112,111,117,118)-c(5,30,49,23,8,48)/60
city<-c("Calgary","Edmonton","Lethbridge","Fort McMurray","High Level","Grande Prairie")
m<-mapproject(c(km2.sc$Long,city.x,km2.water$Long),c(km2.sc$TrueLat,city.y,km2.water$Lat),projection="albers",par=c(49,60))  # These have to be projected together, otherwise get slightly different results, don't know why
km2.sc$proj.x<-m$x[1:nrow(km2.sc)]  # Projected x-value for km2 raster longitudes
km2.sc$proj.y<-m$y[1:nrow(km2.sc)]  # Projected y-value for km2 raster latitudes
m.city.x<-m$x[(nrow(km2.sc)+1):(nrow(km2.sc)+length(city))]  # Projected x for cities locations
m.city.y<-m$y[(nrow(km2.sc)+1):(nrow(km2.sc)+length(city))]  # Projected y for cities locations
m.water<-data.frame(x=m$x[(nrow(km2.sc)+length(city)+1):length(m$x)],y=m$y[(nrow(km2.sc)+length(city)+1):length(m$y)])
m.title<-mapproject(-115,60.5,projection="albers",par=c(49,60))  # Just coordinates for the map title
c1<-rev(c("#D73027","#FC8D59","#FEE090","#E0F3F8","#91BFDB","#4575B4"))  # Colour gradient for reference and current
c2<-colorRampPalette(c1, space = "rgb") # Function to interpolate among these colours for reference and current
d1<-c("#C51B7D","#E9A3C9","#FDE0EF","#E6F5D0","#A1D76A","#4D9221")  # Colour gradient for difference map
d2<-colorRampPalette(d1, space = "rgb") # Function to interpolate among these colours for difference map

# Set up variables to store results
# .pa used for presence/absence components, .agp for abundance-given-presence, .sc for space/climate, .s for summer, .w for winter
vnames.pa<-c("SpruceR","Spruce1","Spruce2","Spruce3","Spruce4","Spruce5","Spruce6","Spruce7","Spruce8","PineR","Pine1","Pine2","Pine3","Pine4","Pine5","Pine6","Pine7","Pine8","DecidR","Decid1","Decid2","Decid3","Decid4","Decid5","Decid6","Decid7","Decid8",
             "MixedwoodR","Mixedwood1","Mixedwood2","Mixedwood3","Mixedwood4","Mixedwood5","Mixedwood6","Mixedwood7","Mixedwood8","TreedBogR","TreedBog1","TreedBog2","TreedBog3","TreedBog4","TreedBog5","TreedBog6","TreedBog7","TreedBog8",
             "TreedFen","TreedSwamp","Shrub","GrassHerb","GraminoidFen","Marsh","ShrubbyBog","ShrubbyFen","ShrubbySwamp",
             "CCSpruceR","CCSpruce1","CCSpruce2","CCSpruce3","CCSpruce4","CCPineR","CCPine1","CCPine2","CCPine3","CCPine4","CCDecidR","CCDecid1","CCDecid2","CCDecid3","CCDecid4","CCMixedwoodR","CCMixedwood1","CCMixedwood2","CCMixedwood3","CCMixedwood4",
             "Crop","TamePasture","RoughPasture","UrbInd","RuralResInd","Wells","SoftLin","HardLin","WetlandMargin")  # Names of coefficients
Coef.pa.s<-Coef.pa.s.se<-array(0,c(length(SpTable.s),length(vnames.pa)))  # Set up tables to store coefficients and their SE's - Summer
Coef.pa.w<-Coef.pa.w.se<-array(0,c(length(SpTable.w),length(vnames.pa)))  # Set up tables to store coefficients and their SE's - Winter
PooledStandCoef.pa.s<-PooledStandCoef.pa.s.se<-array(0,c(length(SpTable.s),5))  # Extra matrix for pooled coefficients for each species for 5 stand types that have ages estimated - Summer
PooledStandCoef.pa.w<-PooledStandCoef.pa.w.se<-array(0,c(length(SpTable.w),5))  # Extra matrix for pooled coefficients for each species for 5 stand types that have ages estimated - Winter
colnames(Coef.pa.s)<-colnames(Coef.pa.s.se)<-colnames(Coef.pa.w)<-colnames(Coef.pa.w.se)<-vnames.pa
colnames(PooledStandCoef.pa.s)<-colnames(PooledStandCoef.pa.w)<-c("Decid","Mixedwood","Pine","Spruce","TreedBog")  # To store pooled coefficients for these stand types that are later divided into age classes
rownames(Coef.pa.s)<-rownames(Coef.pa.s.se)<-rownames(PooledStandCoef.pa.s)<-SpTable.s
rownames(Coef.pa.w)<-rownames(Coef.pa.w.se)<-rownames(PooledStandCoef.pa.w)<-SpTable.w
vnames.agp<-c("Spruce","Pine","Decid","Mixedwood","TreedBog",
              "TreedFen","TreedSwamp","Shrub","GrassHerb","GraminoidFen","Marsh","ShrubbyBog","ShrubbyFen","ShrubbySwamp",
              "CCSpruceR","CCSpruce1","CCSpruce2","CCSpruce3","CCSpruce4","CCPineR","CCPine1","CCPine2","CCPine3","CCPine4","CCDecidR","CCDecid1","CCDecid2","CCDecid3","CCDecid4","CCMixedwoodR","CCMixedwood1","CCMixedwood2","CCMixedwood3","CCMixedwood4",
              "Crop","TamePasture","RoughPasture","UrbInd","RuralResInd","Wells","SoftLin","HardLin","WetlandMargin")  # Names of coefficients
Coef.agp.s<-Coef.agp.s.se<-array(0,c(length(SpTable.s),length(vnames.agp)))  # Set up tables to store coefficients and their SE's - Summer
Coef.agp.w<-Coef.agp.w.se<-array(0,c(length(SpTable.w),length(vnames.agp)))  # Set up tables to store coefficients and their SE's - Winter
colnames(Coef.agp.s)<-colnames(Coef.agp.s.se)<-colnames(Coef.agp.w)<-colnames(Coef.agp.w.se)<-vnames.agp
rownames(Coef.agp.s)<-rownames(Coef.agp.s.se)<-SpTable.s
rownames(Coef.agp.w)<-rownames(Coef.agp.w.se)<-SpTable.w
Coef.mean.s<-Coef.lci.s<-Coef.uci.s<-array(0,c(length(SpTable.s),length(vnames.pa)))  # To store combined presence * abundance|presence estimates and their CI's (CI's because not symmetrical SE's) - Summer
Coef.mean.w<-Coef.lci.w<-Coef.uci.w<-array(0,c(length(SpTable.w),length(vnames.pa)))  # To store combined presence * abundance|presence estimates and their CI's (CI's because not symmetrical SE's) - Winter
colnames(Coef.mean.s)<-colnames(Coef.lci.s)<-colnames(Coef.uci.s)<-colnames(Coef.mean.w)<-colnames(Coef.lci.w)<-colnames(Coef.uci.w)<-vnames.pa
rownames(Coef.mean.s)<-rownames(Coef.lci.s)<-rownames(Coef.uci.s)<-SpTable.s
rownames(Coef.mean.w)<-rownames(Coef.lci.w)<-rownames(Coef.uci.w)<-SpTable.w
Coef.pa.all<-Coef.pa.lci.all<-Coef.pa.uci.all<-Coef.mean.all<-Coef.lci.all<-Coef.uci.all<-array(0,c(length(SpTable.all),length(vnames.pa)))  # To store combined presence * abundance|presence estimates and their CI's (CI's because not symmetrical SE's) - average of Summer and Winter
colnames(Coef.pa.all)<-colnames(Coef.pa.lci.all)<-colnames(Coef.pa.uci.all)<-colnames(Coef.mean.all)<-colnames(Coef.lci.all)<-colnames(Coef.uci.all)<-vnames.pa
rownames(Coef.pa.all)<-rownames(Coef.pa.lci.all)<-rownames(Coef.pa.uci.all)<-rownames(Coef.mean.all)<-rownames(Coef.lci.all)<-rownames(Coef.uci.all)<-SpTable.all
vnames.sc<-c("Intercept","Lat","Long","LatLong","Lat2","Lat3","Long2","Lat2Long2","PET","AHM","MAT","FFP","MAP","MAPFFP","MAPPET","MATAHM","MWMT","MCMT","MWMT2","MAT2","LongMAT","NSR1Parkland","NSR1DryMixedwood","NSR1CentralMixedwood","NSR1Foothills","NSR1North","NSR1Shield","NSR1Mountain")
Res.coef<-array(0,c(length(SpTable.all),length(vnames.sc))) # Done for both seasons together
colnames(Res.coef)<-vnames.sc
rownames(Res.coef)<-SpTable.all
lure.pa.s<-lure.agp.s<-lure.pa.w<-lure.agp.w<-NULL  # To store lure coefficients for pres/abs, abundance | presence
aic.age.s<-array(NA,c(length(SpTable.s),3))  # AIC for age models using 5 stand separately, 2 pairs + black spruce, all stands together - Summer
aic.age.w<-array(NA,c(length(SpTable.w),3))  # AIC for age models using 5 stand separately, 2 pairs + black spruce, all stands together - Winter
cutoff.pc.for.age<-0.1  # Set the cut-off for the minimum proportion of a stand type for a site to be included in the age analysis for a stand type (currently the same for each stand type - could be flexible).  Note: weighting is proportional to proportion of site, so low values can be used - those sites just won't contribute much.  Not that any of this matters with point veg+HF...
# Arrays to save aic weights for each age model, and for which sets of models are best.
aic.wt.pa.save.s<-array(NA,c(length(SpTable.all),16))  # To save aic weights for each species, presence/absence - Summer.  Need to change the second dimension if number of models changes
aic.wt.pa.save.w<-array(NA,c(length(SpTable.all),16))  # To save aic weights for each species, presence/absence - Winter.  Need to change the second dimension if number of models changes
aic.wt.agp.save.s<-array(NA,c(length(SpTable.all),17))  # To save aic weights for each species, abundance|presence.  Include null. - Summer.  Need to change the second dimension if number of models changes
aic.wt.agp.save.w<-array(NA,c(length(SpTable.all),17))  # To save aic weights for each species, abundance|presence.  Include null. - Winter.  Need to change the second dimension if number of models changes
rownames(aic.wt.pa.save.s)<-rownames(aic.wt.agp.save.s)<-SpTable.all
rownames(aic.wt.pa.save.w)<-rownames(aic.wt.agp.save.w)<-SpTable.all
aic.wt.age.save.s<-array(NA,c(length(SpTable.all),3))  # This saves the AIC weights for the 3 combining age models for each species. - Summer
aic.wt.age.save.w<-array(NA,c(length(SpTable.all),3))  # This saves the AIC weights for the 3 combining age models for each species. - Winter
aic.wt.age.models.save.s<-array(NA,c(length(SpTable.all),8))  # This saves the AIC weights for the spline model (relative to the null) for each of the eight stand types or combinations - Summer
aic.wt.age.models.save.w<-array(NA,c(length(SpTable.all),8))  # This saves the AIC weights for the spline model (relative to the null) for each of the eight stand types or combinations - Winter
dimnames(aic.wt.age.save.s)[[1]]<-dimnames(aic.wt.age.models.save.s)[[1]]<-SpTable.all
dimnames(aic.wt.age.save.w)[[1]]<-dimnames(aic.wt.age.models.save.w)[[1]]<-SpTable.all
dimnames(aic.wt.age.save.s)[[2]]<-dimnames(aic.wt.age.save.w)[[2]]<-c("Separate","Intermediate","Combined")
dimnames(aic.wt.age.models.save.s)[[2]]<-dimnames(aic.wt.age.models.save.w)[[2]]<-c("Spruce","Pine","Decid","Mixedwood","TreedBog","UpCon","DecidMixed","All")
auc.fit<-NULL  # AUC for the presence/absence fit, one value for each species

# Extra list of variables for direct summary of densities by veg+HF class.  Need to eliminate unsampled types.  This summary is just meant as a test of the resonableness of the fitted coefficients.
vnames.sum<-vnames.pa[vnames.pa %in% names(d)]
n.sum<-colSums(sign(d[,vnames.sum]))
vnames.sum<-vnames.sum[which(colSums(d[,vnames.sum])>0)]  # A couple veg+HF types have columns, but no samples in used dataset
n.sum<-colSums(sign(d[,vnames.sum]))
which(colSums(sign(d[d$SummerDays>=10,vnames.sum]))==0)  # Check if any of these veg+HF types were not sampled in Summer
which(colSums(sign(d[d$WinterDays>=10,vnames.sum]))==0)  # Check if any of these veg+HF types were not sampled in winter
simple.sum.mean.s<-simple.sum.n.s<-simple.sum.se.s<-array(NA,c(length(SpTable.s),length(vnames.sum)))  # Summer
simple.sum.mean.w<-simple.sum.n.w<-simple.sum.se.w<-array(NA,c(length(SpTable.w),length(vnames.sum)))  # Winter
colnames(simple.sum.mean.s)<-colnames(simple.sum.n.s)<-colnames(simple.sum.se.s)<-colnames(simple.sum.mean.w)<-colnames(simple.sum.n.w)<-colnames(simple.sum.se.w)<-vnames.sum
rownames(simple.sum.mean.s)<-rownames(simple.sum.n.s)<-rownames(simple.sum.se.s)<-SpTable.s
rownames(simple.sum.mean.w)<-rownames(simple.sum.n.w)<-rownames(simple.sum.se.w)<-SpTable.w

# Loop through species
# For each species, models are fit for summer (and coefficient figures plotted), then for winter (and coefficient figures plotted), then results are combined (and plotted) and space/climate residual models are fit (and maps are plotted)
for (sp2 in 1:length(SpTable.all)) {
  d1<-data.frame(DeploymentYear=d$DeploymentYear,Lured=d$Lured,count.summer=NA,count.winter=NA,p.pa.summer=NA,p.agp.summer=NA,p.summer=NA,p.pa.winter=NA,p.agp.winter=NA,p.winter=NA)  # Set up data.frame to keep counts and predictions for the species for each season - to use in SC modeling using residuals after year-round average prediction
  # Summer
  if (paste(SpTable.all[sp2],"Summer",sep="") %in% SpTable.s == TRUE) {
    sp<-which(SpTable.s==paste(SpTable.all[sp2],"Summer",sep=""))
    print(paste(sp,length(SpTable.s),SpTable.s[sp],date()))
    d.sp<-d[,c(1:(FirstSpCol.s-1),(LastSpCol.w+1):ncol(d),which(colnames(d)==SpTable.s[sp]))]  # Extract site descriptors and just the target species.  Assumes species columns are ASummer, AWinter...ZSummer, ZWinter
    colnames(d.sp)[ncol(d.sp)]="Count"  # Change the species count column name to "Count"
    d.sp<-d.sp[!is.na(d.sp$Count) & d.sp$SummerDays>10,]  # Omit deployments with no sampling in that season or too few days
    d.sp$Lured<-as.character(d.sp$Lured)
    q<-by(sign(d.sp$Count[d.sp$Study=="ABMI"]),d.sp$Lured[d.sp$Study=="ABMI"],mean)  # Only using (paired) ABMI sites
    lure.pa.s[sp2]<-q["y"]/q["n"]

    # S0. Simple summary of density in each veg+HF type that has been sampled
    lure1<-mean(d.sp$Count[d.sp$Study=="ABMI" & d.sp$Lured=="y"])/mean(d.sp$Count[d.sp$Study=="ABMI" & d.sp$Lured=="n"])  # Lure effect on total abundance
    for (i in 1:length(vnames.sum)) {
      x<-d.sp$Count[d.sp[,vnames.sum[i]]==1]
      l<-d.sp$Lured[d.sp[,vnames.sum[i]]==1]
      x<-ifelse(l=="y",x/lure1,x)
      simple.sum.mean.s[sp,i]<-mean(x)
      simple.sum.n.s[sp,i]<-sum(sign(x))
      if (n.sum[i]>1) simple.sum.se.s[sp,i]<-sd(x)/sqrt(n.sum[sp])
    }

    # S1. Fit models to broad veg types and HF types - presence/absence
    # S1.1 Fit broad veg models
    m.pa.s<-list(NULL)
    pCount<-sign(d.sp$Count)/ifelse(d.sp$Lured=="y",lure.pa.s[sp2],1)  # Standardize all to no-lure
    pCount<-pCount/max(pCount)  # This is in case the lure effect is <1
    d.sp$pCount.pa<-pCount  # Used in the age models below
    m.pa.s[[1]]<-try(glm(pCount~Decid+Mixedwood+Pine+Spruce+TreedBog+TreedFen+TreedSwamp+GrassHerb+Shrub+Marsh+ShrubbySwamp+ShrubbyBogFen+CCDecidR+CCDecid1+CCDecid2+CCMixedwoodR+CCMixedwood1+CCMixedwood2+CCPineR+CCPine1+CCSpruceR+CCSpruce1+CCSpruce2+SoftLin+TamePasture+RoughPasture+Wells+RuralResInd+WetlandMargin+SummerDays,family="binomial",data=d.sp,weights=d.sp$wt.s))  # Intercept is crop. No UrbInd sampled
    m.pa.s[[2]]<-try(glm(pCount~Decid+Mixedwood+Pine+Spruce+TreedBog+TreedFen+TreedSwamp+GrassHerb+Shrub+Marsh+ShrubbySwamp+ShrubbyBogFen+CCDecidMixed+CCPine+CCSpruce+SoftLin+TamePasture+RoughPasture+Wells+RuralResInd+WetlandMargin+SummerDays,family="binomial",data=d.sp,weights=d.sp$wt.s))  # Intercept is crop. No UrbInd sampled
    m.pa.s[[3]]<-try(glm(pCount~Decid+Mixedwood+Pine+Spruce+TreedBog+TreedFen+TreedSwamp+GrassHerb+Shrub+Marsh+ShrubbySwamp+ShrubbyBogFen+CCDecidMixed+CCPine+CCSpruce+SoftLin+WetlandMargin+SummerDays,family="binomial",data=d.sp,weights=d.sp$wt.s))  # Intercept is alien
    m.pa.s[[4]]<-try(glm(pCount~Decid+Mixedwood+Pine+Spruce+TreedBog+TreedFen+TreedSwamp+GrassHerb+Shrub+Marsh+ShrubbySwamp+ShrubbyBogFen+CCAll+SoftLin+WetlandMargin+SummerDays,family="binomial",data=d.sp,weights=d.sp$wt.s))  # Intercept is alien.  Weak data for cutblocks
    m.pa.s[[5]]<-try(glm(pCount~Decid+Mixedwood+Pine+Spruce+TreedBogFen+TreedSwamp+GrassHerb+Shrub+OpenWet+CCAll+SoftLin+WetlandMargin+SummerDays,family="binomial",data=d.sp,weights=d.sp$wt.s))  # Intercept is alien.
    m.pa.s[[6]]<-try(glm(pCount~Decid+Mixedwood+Pine+Spruce+TreedWet+GrassShrub+OpenWet++CCDecidR+CCDecid1+CCDecid2+CCMixedwoodR+CCMixedwood1+CCMixedwood2+CCPineR+CCPine1+CCSpruceR+CCSpruce1+CCSpruce2+SoftLin+WetlandMargin+SummerDays,family="binomial",data=d.sp,weights=d.sp$wt.s))  # Intercept is alien
    m.pa.s[[7]]<-try(glm(pCount~Decid+Mixedwood+Pine+Spruce+TreedWet+GrassShrub+OpenWet+CCDecidMixed+CCConif+SoftLin+WetlandMargin+SummerDays,family="binomial",data=d.sp,weights=d.sp$wt.s))  # Intercept is alien
    m.pa.s[[8]]<-try(glm(pCount~DecidMixed+UpCon+TreedWet+GrassShrub+OpenWet+CCDecidMixed+CCConif+SoftLin+WetlandMargin+SummerDays,family="binomial",data=d.sp,weights=d.sp$wt.s))  # Intercept is alien
    m.pa.s[[9]]<-try(glm(pCount~DecidMixed+UpCon+TreedWet+GrassShrub+OpenWet+CCAll+SoftLin+TamePasture+RoughPasture+Wells+RuralResInd+WetlandMargin+SummerDays,family="binomial",data=d.sp,weights=d.sp$wt.s))  # Intercept is crop
    m.pa.s[[10]]<-try(glm(pCount~DecidMixed+UpCon+TreedWet+GrassShrub+OpenWet+CCAll+SoftLin+WetlandMargin+SummerDays,family="binomial",data=d.sp,weights=d.sp$wt.s))  # Intercept is alien
    m.pa.s[[11]]<-try(glm(pCount~DecidMixed+UpCon+TreedWet+GrassShrub+OpenWet+Succ+WetlandMargin+SummerDays,family="binomial",data=d.sp,weights=d.sp$wt.s))  # Intercept is alien
    m.pa.s[[12]]<-try(glm(pCount~Upland+Lowland+CCAll+SoftLin+WetlandMargin+SummerDays,family="binomial",data=d.sp,weights=d.sp$wt.s))  # Intercept is alien
    m.pa.s[[13]]<-try(glm(pCount~TreedAll+OpenAll+Succ+WetlandMargin+SummerDays,family="binomial",data=d.sp,weights=d.sp$wt.s))  # Intercept is alien
    m.pa.s[[14]]<-try(glm(pCount~UplandForest+Lowland+GrassShrub+Succ+WetlandMargin+SummerDays,family="binomial",data=d.sp,weights=d.sp$wt.s))  # Intercept is alien
    m.pa.s[[15]]<-try(glm(pCount~Upland+Lowland+Succ+WetlandMargin+SummerDays,family="binomial",data=d.sp,weights=d.sp$wt.s))  # Intercept is alien
    m.pa.s[[16]]<-try(glm(pCount~Boreal+GrassShrub+Succ+WetlandMargin+SummerDays,family="binomial",data=d.sp,weights=d.sp$wt.s))  # "Boreal" is everything native except GrassShrub. Intercept is alien
    # AIC calculation  (I'm using AIC here, because this is primarily for prediction, rather than finding a minimal best model)
    nModels.pa<-length(m.pa.s)
    aic.ta<-rep(999999999,(nModels.pa))
    for (i in 1:(nModels.pa)) {
      if (!is.null(m.pa.s[[i]]) & class(m.pa.s[[i]])[1]!="try-error") {  # last part is to not used non-converged models, unless none converged
        aic.ta[i]<-AICc(m.pa.s[[i]])
      }
    }
    aic.delta<-aic.ta-min(aic.ta)
    aic.exp<-exp(-1/2*aic.delta)
    aic.wt.pa<-aic.exp/sum(aic.exp)
    best.model.pa<-which.max(aic.wt.pa)
    aic.wt.pa.save.s[sp2,]<-aic.wt.pa

    # S1.2 Then do abundance | presence model
    d.p<-d.sp[d.sp$Count>0,]  # Just presence records
    q<-by(d.p$Count[d.sp$Study=="ABMI"],d.p$Lured[d.sp$Study=="ABMI"],mean)
    lure.agp.s[sp2]<-q["y"]/q["n"]
    pCount<-d.p$Count/ifelse(d.p$Lured=="y",lure.agp.s[sp2],1)
    j<-0
    m.agp.s<-list(NULL)
    mnums<-NULL  # To keep track of which models were actually fit
    for (i in 1:nModels.pa) {
      x<-colSums(d.p[,attr(m.pa.s[[i]]$terms,"term.labels")])  # Minimum number of presence records in each veg+HF type included in the model
      if (min(x)>3 & nrow(d.p)-sum(x[-length(x)])) {  # Only use the model if each type is represented by >3 presence records (incl. intercept)
        j<-j+1
        m.agp.s[[j]]<-glm(m.pa.s[[i]]$formula,data=d.p,family=gaussian(link="log"),weights=d.p$wt)
        mnums<-c(mnums,i)
      }
    }
    m.agp.s[[j+1]]<-glm(pCount~SummerDays,data=d.p,family=gaussian(link="log"),weights=d.p$wt)  # And add a null
    mnums<-c(mnums,nModels.pa+1)  # Add the null model to list - 1 beyond the last pa model
    # AIC calculation  (I'm using AIC here, because this is primarily for prediction, rather than finding a minimal best model)
    nModels.agp<-length(m.agp.s)
    aic.ta<-rep(999999999,(nModels.agp))
    for (i in 1:(nModels.agp)) {
      if (!is.null(m.agp.s[[i]]) & class(m.agp.s[[i]])[1]!="try-error") {  # last part is to not used non-converged models, unless none converged
        aic.ta[i]<-AICc(m.agp.s[[i]])
      }
    }
    aic.delta<-aic.ta-min(aic.ta)
    aic.exp<-exp(-1/2*aic.delta)
    aic.wt.agp<-aic.exp/sum(aic.exp)
    best.model.agp<-which.max(aic.wt.agp)
    aic.wt.agp.save.s[sp2,mnums]<-aic.wt.agp

    # S1.3 Predict from models for 100% of each veg type and for each site - pres/abs - using only the best model
    Intercept<-c("Crop","Crop",rep("Alien",6),"Crop",rep("Alien",7))[best.model.pa]  # Change this if models change
    terms.pa<-c(attr(m.pa.s[[best.model.pa]]$terms,"term.labels"),Intercept)
    Coef1.pa<-Coef1.pa.se<-rep(NA,length(terms.pa))
    names(Coef1.pa)<-names(Coef1.pa.se)<-terms.pa
    for (i in 1:length(terms.pa)) {
      pm1<-rep(0,length(terms.pa))
      pm1[i]<-1
      names(pm1)<-terms.pa
      pm1<-data.frame(t(pm1))
      p<-predict(m.pa.s[[best.model.pa]],newdata=data.frame(SummerDays=100,pm1),se.fit=T)
      Coef1.pa[i]<-plogis(p$fit)  # Ordinal scale
      Coef1.pa.se[i]<-p$se.fit  # logit scale
    }
    # Adjust so that mean prediction at standardized SummerDays across all qualifying sites = mean observed count.  This is to compensate for inaccuracies due to fitting SummerDays coefficients
    pCount<-sign(d.sp$Count)/ifelse(d.sp$Lured=="y",lure.pa.s[sp2],1)  # Standardize all to no-lure
    pCount<-pCount/max(pCount)  # This is in case the lure effect is <1
    p.adj<-predict(m.pa.s[[best.model.pa]],newdata=data.frame(SummerDays=100,d.sp))  # Prediction for each data point, but at standardized days
    Coef1.pa<-plogis(qlogis(Coef1.pa)+qlogis(mean(pCount))-qlogis(mean(plogis(p.adj))))  # Adjustment made on logit scale.  Var doesn't change on logit scale(?)
    # And predict for each site - used for age model below
    d.sp$p<-predict(m.pa.s[[best.model.pa]])  # Logit scale

    # S1.4 Predict from models for 100% of each veg type and for each site - abund|pres
    terms.pa1<-terms.pa[terms.pa!="SummerDays"]
    Coef1.agp<-Coef1.agp.se<-rep(NA,length(terms.pa1))  # For agp coefficients averaged to terms in best pa model
    names(Coef1.agp)<-names(Coef1.agp.se)<-terms.pa1
    i<-best.model.agp
    p<-predict(m.agp.s[[i]],newdata=data.frame(pm,SummerDays=100),se.fit=T)  # For each type
    tTypeMean.agp<-p$fit
    tTypeVar.agp<-p$se.fit^2  # Using variance, for consistency with other scripts at this point
    names(tTypeMean.agp)<-names(tTypeVar.agp)<-pm$VegType
    # Adjust so that mean prediction at standardized SummerDays across all qualifying sites = mean observed count.  This is to compensate for inaccuracies due to fitting SummerDays coefficients
    pCount<-d.p$Count/ifelse(d.p$Lured=="y",lure.agp.s[sp2],1)
    p.adj<-predict(m.agp.s[[i]],newdata=data.frame(SummerDays=100,d.sp))  # Prediction for each data point, but at standardized days
    tTypeMean.agp<-log(exp(tTypeMean.agp)+mean(pCount)-mean(exp(p.adj)))  # Var doesn't change on log scale
    # Then average those for each broader group included in the best pa model
    for (i in 1:length(terms.pa1)) {
      j<-pm$VegType[which(pm[,terms.pa1[i]]==1)]  # Names of fine hab+HF types included in that broader group
      x<-tTypeMean.agp[as.character(j)]
      x.var<-tTypeVar.agp[as.character(j)]
      Coef1.agp[i]<-mean(x)  # Simple mean, log scale
      Coef1.agp.se[i]<-sqrt(mean(x.var))  # This is for straight-up mean, log scale
    }
    # Put non-stand-age coefficients in Coef matrix (and SE's)
    Coef.agp.s[sp,]<-exp(tTypeMean.agp[na.omit(match(colnames(Coef.agp.s),names(tTypeMean.agp)))])  # On ordinal scale
    Coef.agp.s.se[sp,]<-sqrt(tTypeVar.agp[na.omit(match(colnames(Coef.agp.s),names(tTypeMean.agp)))])  # On log scale

    # S2. Run models for age within each stand type - pres/abs only currently
    # S2.1. Set up separate dataframes for sites containing a minimum amount of each broad stand type
    d.spruce<-d.pine<-d.decid<-d.mixed<-d.treedbog<-NULL
    for (i in 0:8) {  # Add sites with each age class of the stand type to a separate data frames for each stand type
      i1<-ifelse(i==0,"R",i)  # For variable name
      i2<-ifelse(i==0,0.5,i)  # For twenty-year age
      cn<-paste("Spruce",i1,sep="")  # Col name
      if(cn %in% names(d.sp)) if (sum(d.sp[,cn]>cutoff.pc.for.age)>0) d.spruce<-rbind(d.spruce,data.frame(pCount=d.sp[d.sp[,cn]>cutoff.pc.for.age,"pCount.pa"], age=i2, wt1=d.sp[d.sp[,cn]>cutoff.pc.for.age,cn]*d.sp$wt.s[d.sp[,cn]>cutoff.pc.for.age], p=d.sp$p[d.sp[,cn]>cutoff.pc.for.age]))  # Weight is the ppoportion of the site of that age class and stand type, multiplied by the original weight (which accounts for revisited sites)
      cn<-paste("Pine",i1,sep="")  # Col name
      if(cn %in% names(d.sp)) if (sum(d.sp[,cn]>cutoff.pc.for.age)>0) d.pine<-rbind(d.pine,data.frame(pCount=d.sp[d.sp[,cn]>cutoff.pc.for.age,"pCount.pa"], age=i2, wt1=d.sp[d.sp[,cn]>cutoff.pc.for.age,cn]*d.sp$wt.s[d.sp[,cn]>cutoff.pc.for.age], p=d.sp$p[d.sp[,cn]>cutoff.pc.for.age]))
      cn<-paste("Decid",i1,sep="")  # Col name
      if(cn %in% names(d.sp)) if (sum(d.sp[,cn]>cutoff.pc.for.age)>0) d.decid<-rbind(d.decid,data.frame(pCount=d.sp[d.sp[,cn]>cutoff.pc.for.age,"pCount.pa"], age=i2, wt1=d.sp[d.sp[,cn]>cutoff.pc.for.age,cn]*d.sp$wt.s[d.sp[,cn]>cutoff.pc.for.age], p=d.sp$p[d.sp[,cn]>cutoff.pc.for.age]))
      cn<-paste("Mixedwood",i1,sep="")  # Col name
      if(cn %in% names(d.sp)) if (sum(d.sp[,cn]>cutoff.pc.for.age)>0) d.mixed<-rbind(d.mixed,data.frame(pCount=d.sp[d.sp[,cn]>cutoff.pc.for.age,"pCount.pa"], age=i2, wt1=d.sp[d.sp[,cn]>cutoff.pc.for.age,cn]*d.sp$wt.s[d.sp[,cn]>cutoff.pc.for.age], p=d.sp$p[d.sp[,cn]>cutoff.pc.for.age]))
      cn<-paste("TreedBog",i1,sep="")  # Col name
      if(cn %in% names(d.sp)) if (sum(d.sp[,cn]>cutoff.pc.for.age)>0) d.treedbog<-rbind(d.treedbog,data.frame(pCount=d.sp[d.sp[,cn]>cutoff.pc.for.age,"pCount.pa"], age=i2, wt1=d.sp[d.sp[,cn]>cutoff.pc.for.age,cn]*d.sp$wt.s[d.sp[,cn]>cutoff.pc.for.age], p=d.sp$p[d.sp[,cn]>cutoff.pc.for.age]))
    }
    # Note: Not using grass and shrub here, because including parklands, so most of these are prairies, not recent fires
    # Combined data frames for models using more than one stand type
    d.upcon<-rbind(cbind(d.spruce,Sp="Spruce"),cbind(d.pine,Sp="Pine"))
    d.decidmixed<-rbind(cbind(d.decid,Sp="Decid"),cbind(d.mixed,Sp="Mixed"))
    d.all<-rbind(cbind(d.spruce,Sp="Spruce"),cbind(d.pine,Sp="Pine"),cbind(d.decid,Sp="Decid"),cbind(d.mixed,Sp="Mixed"),cbind(d.treedbog,Sp="TreedBog"))
    # S2.2. Then fit models of age functions
    # Note: Days not included in these models, because its effect on presence/absence is already included in the prediction p used in the offset
    m.age.s<-m.null.age.s<-list(NULL)  # Age models for each of the 5 stand types plus 3 combinations.  Null model is for AIC comparison to smoothing spline model (which can be worse, even if it has >1 equiv DF)
    m.null.age.s[[1]]<-gam(pCount~1+offset(p),data=d.spruce,family="binomial",weights=d.spruce$wt1)
    if (sum(sign(d.spruce$pCount))>4) {
      m.age.s[[1]]<-gam(pCount~s(sqrt(age),k=4,m=2)+offset(p),data=d.spruce,family="binomial",weights=d.spruce$wt1)  # Fit spline through age data for that stand type if enough records
    } else {
      m.age.s[[1]]<-m.null.age.s[[1]]  # Use constant-only if too few records in that stand type
    }
    m.null.age.s[[2]]<-gam(pCount~1+offset(p),data=d.pine,family="binomial",weights=d.pine$wt1)
    if (sum(sign(d.pine$pCount))>4) {
      m.age.s[[2]]<-gam(pCount~s(sqrt(age),k=4,m=2)+offset(p),data=d.pine,family="binomial",weights=d.pine$wt1)
    } else {
      m.age.s[[2]]<-m.null.age.s[[2]]  # Use constant-only if too few records in that stand type
    }
    m.null.age.s[[3]]<-gam(pCount~1+offset(p),data=d.decid,family="binomial",weights=d.decid$wt1)
    if (sum(sign(d.decid$pCount))>4) {
      m.age.s[[3]]<-gam(pCount~s(sqrt(age),k=4,m=2)+offset(p),data=d.decid,family="binomial",weights=d.decid$wt1)
    } else {
      m.age.s[[3]]<-m.null.age.s[[3]]  # Use constant-only if too few records in that stand type
    }
    m.null.age.s[[4]]<-gam(pCount~1+offset(p),data=d.mixed,family="binomial",weights=d.mixed$wt1)
    if (sum(sign(d.mixed$pCount))>4) {
      m.age.s[[4]]<-gam(pCount~s(sqrt(age),k=4,m=2)+offset(p),data=d.mixed,family="binomial",weights=d.mixed$wt1)
    } else {
      m.age.s[[4]]<-m.null.age.s[[4]]  # Use constant-only if too few records in that stand type
    }
    m.null.age.s[[5]]<-gam(pCount~1+offset(p),data=d.treedbog,family="binomial",weights=d.treedbog$wt1)
    if (sum(sign(d.treedbog$pCount))>4) {
      m.age.s[[5]]<-gam(pCount~s(sqrt(age),k=4,m=2)+offset(p),data=d.treedbog,family="binomial",weights=d.treedbog$wt1)
    } else {
      m.age.s[[5]]<-m.null.age.s[[5]]  # Use constant-only if too few records in that stand type
    }
    m.null.age.s[[6]]<-gam(pCount~1+offset(p),data=d.upcon,family="binomial",weights=d.upcon$wt1)
    if (sum(sign(d.upcon$pCount))>4) {
      m.age.s[[6]]<-gam(pCount~s(sqrt(age),k=4,m=2)+offset(p),data=d.upcon,family="binomial",weights=d.upcon$wt1)
    } else {
      m.age.s[[6]]<-m.null.age.s[[6]]  # Use constant-only if too few records in that stand type
    }
    m.null.age.s[[7]]<-gam(pCount~1+offset(p),data=d.decidmixed,family="binomial",weights=d.decidmixed$wt1)
    if (sum(sign(d.decidmixed$pCount))>4) {
      m.age.s[[7]]<-gam(pCount~s(sqrt(age),k=4,m=2)+offset(p),data=d.decidmixed,family="binomial",weights=d.decidmixed$wt1)
    } else {
      m.age.s[[7]]<-m.null.age.s[[7]]  # Use constant-only if too few records in that stand type
    }
    m.null.age.s[[8]]<-gam(pCount~1+offset(p),data=d.all,family="binomial",weights=d.all$wt1)
    if (sum(sign(d.all$pCount))>4) {
      m.age.s[[8]]<-gam(pCount~s(sqrt(age),k=4,m=2)+offset(p),data=d.all,family="binomial",weights=d.all$wt1)
    } else {
      m.age.s[[8]]<-m.null.age.s[[8]]  # Use constant-only if too few records in that stand type
    }
    # AIC for options - each separate, pairwise combinations, combine all
    aic.age.s[sp,1]<-AIC(m.age.s[[1]])+AIC(m.age.s[[2]])+AIC(m.age.s[[3]])+AIC(m.age.s[[4]])+AIC(m.age.s[[5]])-8  # Minus 8 for the 4 extra variances
    aic.age.s[sp,2]<-AIC(m.age.s[[5]])+AIC(m.age.s[[6]])+AIC(m.age.s[[7]])-4
    aic.age.s[sp,3]<-AIC(m.age.s[[8]])
    aic.delta<-aic.age.s[sp,]-min(aic.age.s[sp,])
    aic.exp<-exp(-1/2*aic.delta)
    aic.wt.age<-aic.exp/sum(aic.exp)
    aic.wt.age.models<-NULL
    for (i in 1:8) {  # AIC weight for spline (versus null) for each age grouping
      q<-c(AICc(m.age.s[[i]]),AICc(m.null.age.s[[i]]))
      q<-q-min(q)
      aic.wt.age.models[i]<-exp(-1/2*q[1])/sum(exp(-1/2*q))
    }
    names(aic.wt.age.models)<-c("Spruce","Pine","Decid","Mixedwood","TreedBog","DecidMixed","UpCon","TreedAll")
    # S2.3 Then make age predictions only if the spline model is better than the null, and only for (groupings of) stand types in the pa model
    age.flag<-rep(0,10)  # Flag for separate age predictions for each stand type {Spruce, Pine, Deciduous, Mixedwood, TreedBog, TreedWet (includes TreedSwamp), DecidMixed, UpCon, UplandForest, TreedAll) (0=no, 1=yes)
    p.age<-p.age.se<-array(NA,c(10,9))  # Predictions for 10 stand types, 9 age classes
    rownames(p.age)<-rownames(p.age.se)<-names(age.flag)<-c("Spruce","Pine","Decid","Mixedwood","TreedBog","TreedWet","DecidMixed","UpCon","UplandForest","TreedAll")
    if ("Spruce" %in% terms.pa) {  # Do spruce separately
      if (aic.wt.age.models["Spruce"]>0.5) {
        age.flag[1]<-1
        p<-predict(m.age.s[[1]],newdata=data.frame(age=c(0.5,1:8),p=qlogis(Coef1.pa["Spruce"])),se.fit=T)
        p.age["Spruce",]<-plogis(p$fit)
        p.age.se["Spruce",]<-p$se.fit  # Still on link scale
      }
    }
    if ("Pine" %in% terms.pa) {  # Do pine separately
      if (aic.wt.age.models["Pine"]>0.5) {
        age.flag["Pine"]<-1
        p<-predict(m.age.s[[2]],newdata=data.frame(age=c(0.5,1:8),p=qlogis(Coef1.pa["Pine"])),se.fit=T)
        p.age["Pine",]<-plogis(p$fit)
        p.age.se["Pine",]<-p$se.fit  # Still on link scale
      }
    }
    if ("Decid" %in% terms.pa) {  # Do deciduous separately
      if (aic.wt.age.models["Decid"]>0.5) {
        age.flag["Decid"]<-1
        p<-predict(m.age.s[[3]],newdata=data.frame(age=c(0.5,1:8),p=qlogis(Coef1.pa["Decid"])),se.fit=T)
        p.age["Decid",]<-plogis(p$fit)
        p.age.se["Decid",]<-p$se.fit  # Still on link scale
      }
    }
    if ("Mixedwood" %in% terms.pa) {  # Do mixedwood separately
      if (aic.wt.age.models["Mixedwood"]>0.5) {
        age.flag["Mixedwood"]<-1
        p<-predict(m.age.s[[4]],newdata=data.frame(age=c(0.5,1:8),p=qlogis(Coef1.pa["Mixedwood"])),se.fit=T)
        p.age["Mixedwood",]<-plogis(p$fit)
        p.age.se["Mixedwood",]<-p$se.fit  # Still on link scale
      }
    }
    if ("TreedBog" %in% terms.pa) {  # Do treedbog separately
      if (aic.wt.age.models["TreedBog"]>0.5) {
        age.flag["TreedBog"]<-1
        p<-predict(m.age.s[[5]],newdata=data.frame(age=c(0.5,1:8),p=qlogis(Coef1.pa["TreedBog"])),se.fit=T)
        p.age["TreedBog",]<-plogis(p$fit)
        p.age.se["TreedBog",]<-p$se.fit  # Still on link scale
      }
    }
    if ("TreedWet" %in% terms.pa) {  # Do TreedWet separately, using TreedBog model
      if (aic.wt.age.models["TreedBog"]>0.5) {
        age.flag["TreedWet"]<-1
        p<-predict(m.age.s[[5]],newdata=data.frame(age=c(0.5,1:8),p=qlogis(Coef1.pa["TreedWet"])),se.fit=T)
        p.age["TreedWet",]<-plogis(p$fit)
        p.age.se["TreedWet",]<-p$se.fit  # Still on link scale
      }
    }
    if ("UpCon" %in% terms.pa) {  # Do spruce+pine together
      if (aic.wt.age.models["UpCon"]>0.5) {
        age.flag["UpCon"]<-1
        p<-predict(m.age.s[[6]],newdata=data.frame(age=c(0.5,1:8),p=qlogis(Coef1.pa["UpCon"])),se.fit=T)
        p.age["UpCon",]<-plogis(p$fit)
        p.age.se["UpCon",]<-p$se.fit  # Still on link scale
      }
    }
    if ("DecidMixed" %in% terms.pa) {  # Do deciduous+mixed together
      if (aic.wt.age.models["DecidMixed"]>0.5) {
        age.flag["DecidMixed"]<-1
        p<-predict(m.age.s[[7]],newdata=data.frame(age=c(0.5,1:8),p=qlogis(Coef1.pa["DecidMixed"])),se.fit=T)
        p.age["DecidMixed",]<-plogis(p$fit)
        p.age.se["DecidMixed",]<-p$se.fit  # Still on link scale
      }
    }
    if ("UplandForest" %in% terms.pa)  { # Do those four stand types together, using "all" age model
      if (aic.wt.age.models["TreedAll"]>0.5) {
        age.flag["UplandForest"]<-1
        p<-predict(m.age.s[[8]],newdata=data.frame(age=c(0.5,1:8),p=qlogis(Coef1.pa["UplandForest"])),se.fit=T)
        p.age["UplandForest",]<-plogis(p$fit)
        p.age.se["UplandForest",]<-p$se.fit  # Still on link scale
      }
    }
    if ("TreedAll" %in% terms.pa)  {  # Do all five together
      if (aic.wt.age.models["TreedAll"]>0.5) {
        age.flag["TreedAll"]<-1
        p<-predict(m.age.s[[8]],newdata=data.frame(age=c(0.5,1:8),p=qlogis(Coef1.pa["TreedAll"])),se.fit=T)
        p.age["TreedAll",]<-plogis(p$fit)
        p.age.se["TreedAll",]<-p$se.fit  # Still on link scale
      }
    }
    # Save aic weights
    aic.wt.age.save.s[sp2,]<-aic.wt.age
    aic.wt.age.models.save.s[sp2,]<-aic.wt.age.models

    # S3.Assemble presence/absence coefficients for each veg+HF type, including expanding out stand types that have ages, to populate entire Coef.pa matrix, then multiply by appropriate Coef.agp term to generate Coef.mean matrix
    for (i in 1:length(terms.pa1)) {
      j<-as.character(pm$VegType[pm[,terms.pa1[i]]==1])  # The fine veg+HF types that are covered by the (potentially broader) variable included in the best model
      for (j1 in 1:length(j)) {
        if (j[j1] %in% c("Spruce","Pine","Decid","Mixedwood","TreedBog")) {
          # Stand types that have age classes in full coefficients
          x<-paste(j[j1],c("R","1","2","3","4","5","6","7","8"),sep="")  # Names in Coef.pa and Coef.mean with ages
          if ((age.flag[terms.pa1[i]]==0)==TRUE | is.na(age.flag[terms.pa1[i]])) {   # But these broad types do not have age classes in the best model
            Coef.pa.s[sp,x]<-Coef1.pa[terms.pa1[i]]  # Fill in all ages with value for that (broad) stand type
            Coef.pa.s.se[sp,x]<-Coef1.pa.se[terms.pa1[i]]  # Fill in all ages with value for that (broad) stand type
            Coef.mean.s[sp,x]<-Coef.pa.s[sp,x]*Coef.agp.s[sp,j[j1]]  # Only the one value for agp
            Coef.lci.s[sp,x]<-plogis(qlogis(Coef.pa.s[sp,x])-Coef.pa.s.se[sp,x]*1.28) * exp(log(Coef.agp.s[sp,j[j1]])-Coef.agp.s.se[sp,j[j1]]*1.28)  # Using 10% intervals for each to multiply to 5% intervals (assuming independence - checked empirically)
            Coef.uci.s[sp,x]<-plogis(qlogis(Coef.pa.s[sp,x])+Coef.pa.s.se[sp,x]*1.28) * exp(log(Coef.agp.s[sp,j[j1]])+Coef.agp.s.se[sp,j[j1]]*1.28)  # Using 10% intervals for each to multiply to 5% intervals (assuming independence - checked empirically)
          } else {  # These broad types do have age classes in the best model
            Coef.pa.s[sp,x]<-p.age[terms.pa1[i],]  # Use separate age values for that (broad) stand type
            Coef.pa.s.se[sp,x]<-p.age.se[terms.pa1[i],]  # Use separate age values for that (broad) stand type
            Coef.mean.s[sp,x]<-Coef.pa.s[sp,x]*Coef.agp.s[sp,j[j1]]  # Only the one value for agp
            Coef.lci.s[sp,x]<-plogis(qlogis(Coef.pa.s[sp,x])-Coef.pa.s.se[sp,x]*1.28) * exp(log(Coef.agp.s[sp,j[j1]])-Coef.agp.s.se[sp,j[j1]]*1.28)  # Using 10% intervals for each to multiply to 5% intervals (assuming independence - checked empirically)
            Coef.uci.s[sp,x]<-plogis(qlogis(Coef.pa.s[sp,x])+Coef.pa.s.se[sp,x]*1.28) * exp(log(Coef.agp.s[sp,j[j1]])+Coef.agp.s.se[sp,j[j1]]*1.28)  # Using 10% intervals for each to multiply to 5% intervals (assuming independence - checked empirically)
          }
        } else {
          # Stand types that do not have age classes in full coefficients
          Coef.pa.s[sp,j[j1]]<-Coef1.pa[terms.pa1[i]]  # Fill in the non-age coefficient
          Coef.pa.s.se[sp,j[j1]]<-Coef1.pa.se[terms.pa1[i]]  # Fill in the non-age coefficient
          Coef.mean.s[sp,j[j1]]<-Coef.pa.s[sp,j[j1]]*Coef.agp.s[sp,j[j1]]
          Coef.lci.s[sp,j[j1]]<-plogis(qlogis(Coef.pa.s[sp,j[j1]])-Coef.pa.s.se[sp,j[j1]]*1.28) * exp(log(Coef.agp.s[sp,j[j1]])-Coef.agp.s.se[sp,j[j1]]*1.28)  # Using 10% intervals for each to multiply to 5% intervals (assuming independence - checked empirically)
          Coef.uci.s[sp,j[j1]]<-plogis(qlogis(Coef.pa.s[sp,j[j1]])+Coef.pa.s.se[sp,j[j1]]*1.28) * exp(log(Coef.agp.s[sp,j[j1]])+Coef.agp.s.se[sp,j[j1]]*1.28)  # Using 10% intervals for each to multiply to 5% intervals (assuming independence - checked empirically)
        }  # End if for stand types with age classes
      }  # Next fine veg+HF type within broader class
    }  # Next broader class in best model
    # Further adjustments so that mean(Coef.mean) = mean(observed density).  Adjusts for geometric mean issue, and also for any remaining issues with Days adjustments
    pCount<-d.sp$Count/(lure.pa.s[sp2]*lure.agp.s[sp2])
    z<-colnames(Coef.mean.s)[!is.na(match(colnames(Coef.mean.s),colnames(d.sp)))]  # Coef names that are in d.sp - excludes added categories not in data, like old clearcuts
    z<-z[!is.na(Coef.mean.s[sp,z])]
    q<-colSums(t(d.sp[,z]) * Coef.mean.s[sp,z])  # Prediction for each deployment, based on the coefficients and its veg composition
    Coef.mean.s[sp,]<-Coef.mean.s[sp,]*mean(pCount)/mean(q)
    Coef.lci.s[sp,]<-Coef.lci.s[sp,]*mean(pCount)/mean(q)
    Coef.uci.s[sp,]<-Coef.uci.s[sp,]*mean(pCount)/mean(q)

    # S3.2 Not using cutblock convergence for best-model analysis (i.e., for "honest" coefficient figures), but it is used below for the final all-coefficients exported
    # Predictions for each site, as offsets below
    d.sp$t.p.pa<-colSums(Coef1.pa[terms.pa1]*t(d.sp[,terms.pa1]))  # Prediction of presence/absence at each site - ordinal scale
    d.sp$t.p.agp<-colSums(Coef1.agp[terms.pa1]*t(d.sp[,terms.pa1]))  # Prediction of abundance|presence at each site - log scale
    d.sp$p.ta<-d.sp$t.p.pa*exp(d.sp$t.p.agp)  # Prediction (ordinal scale) of total abundance at each site

    # S3.3 Store values for Summer
    j<-match(d.sp$DeploymentYear,d1$DeploymentYear)
    d1$count.summer[j]<-d.sp$Count
    d1$p.pa.summer[j]<-d.sp$t.p.pa
    d1$p.agp.summer[j]<-exp(d.sp$t.p.agp)
    d1$p.summer[j]<-d.sp$p.ta

    # S4. Coefficient figures
    # S4.1. Combined total abundance
    # Custom figure for each best model and age options
    x<-col1<-y<-y.lci<-y.uci<-w<-class<-space1<-NULL
    for (i in 1:length(terms.pa1)) {
      if (age.flag[terms.pa1[i]]==0 | is.na(age.flag[terms.pa1[i]])) {  # Variable without age info
        j<-match(terms.pa1[i],fp$Class)
        x<-c(x,fp$x[j])
        veg.to.use<-as.character(pm$VegType[pm[,terms.pa1[i]]==1])[1]  # To use Coef.mean, find the first fine veg+HF types that is included in the (potentially broader) variable being plotted
        if (veg.to.use %in% c("Spruce","Pine","Decid","Mixedwood","TreedBog")) veg.to.use<-paste(veg.to.use,"R",sep="")  # If no age relationship, use the first value for that fine stand type (don't add the R for types that don't have age classes at all)
        y<-c(y,Coef.mean.s[sp,veg.to.use])
        y.lci<-c(y.lci,Coef.lci.s[sp,veg.to.use])
        y.uci<-c(y.uci,Coef.uci.s[sp,veg.to.use])
        col1<-c(col1,fp$col1[j])
        w<-c(w,fp$width[j])
        class<-c(class,terms.pa1[i])
        space1<-c(space1,fp$spaceafter[j])
      } else {  # Variable with age info
        j<-match(paste(terms.pa1[i],"R",sep=""),fp$Class):match(paste(terms.pa1[i],"8",sep=""),fp$Class)
        x<-c(x,fp$x[j])
        veg.to.use<-as.character(pm$VegType[pm[,terms.pa1[i]]==1])[1]  # To use Coef.mean, find the first fine veg+HF types that is included in the (potentially broader) variable being plotted
        veg.to.use<-paste(veg.to.use,c("R","1","2","3","4","5","6","7","8"),sep="")  # If age relationship, use each value for that fine stand type
        y<-c(y,Coef.mean.s[sp,veg.to.use])
        y.lci<-c(y.lci,Coef.lci.s[sp,veg.to.use])
        y.uci<-c(y.uci,Coef.uci.s[sp,veg.to.use])
        col1<-c(col1,fp$col1[j])
        w<-c(w,fp$width[j])
        class<-c(class,as.character(fp$Class[j]))
        space1<-c(space1,fp$spaceafter[j])
      }
    }
    ord<-order(x)  # Sort all by x
    y<-y[ord]
    y.lci<-y.lci[ord]
    y.uci<-y.uci[ord]
    col1<-col1[ord]
    w<-w[ord]
    class<-class[ord]
    space1<-space1[ord]
    x<-x[ord]
    # Rectify spaces between x's
    for (i in 1:(length(x)-1)) {
      for (j in (i+1):length(x)) x[j]<-x[j]+space1[i]-(x[i+1]-x[i])  # Alter all subsequent positions accordingly
    }
    # Make bar plot
    ymax<-min(max(y.uci,na.rm=TRUE),2*max(y,na.rm=TRUE))  # This keeps the figures readable when there are extreme UCI's
    space<-c(1,x[-1]-x[-length(x)])-0.99  # The spacing between bars
    density<-ifelse(substr(class,1,2)=="CC",50,NA)
    fname<-paste(fname.fig,"Veg+HF figure best model ",SpTable.s[sp],".png",sep="")
    png(file=fname,width=ifelse(length(y)>5,1500,1000),height=700)
    par(mai=c(1.9,1,0.2,0.3))
    x1<-barplot(y,space=space,width=w,border="white",col="grey30",ylim=c(0,ymax),yaxt="n",ylab="Relative abundance",col.lab="grey50",cex.lab=1.2,axisnames=F)[,1]  # To get strips on CC bars
    abline(h=pretty(c(0,ymax)),col="grey80")
    x1<-barplot(y,space=space,width=w,border="white",col="grey30",ylim=c(0,ymax),yaxt="n",ylab="Relative abundance",col.lab="grey50",cex.lab=1.2,axisnames=F,add=TRUE)[,1]  # To get strips on CC bars, and to put bars in front of horizontal axis lines
    x1<-barplot(y,space=space,width=w,border="white",density=density,col=col1,ylim=c(0,ymax),yaxt="n",ylab="Relative abundance",col.lab="grey50",cex.lab=1.2,axisnames=F,add=TRUE)[,1]
    axis(side=2,tck=0.02,cex.axis=0.9,col.axis="grey50",col.ticks="grey50",las=2,at=pretty(c(0,ymax)))
    box(bty="l",col="grey50")
    for (i in 1:length(x1)) {
      lines(rep(x1[i],2),c(y[i],y.uci[i]),col=col1[i])
      lines(rep(x1[i],2),c(y[i],y.lci[i]),col="grey90")
    }
    for (i in 1:length(class)) {  # Label x axis
      if (substr(class[i],1,2)=="CC") {
        mtext(side=1,line=0,at=x1[i],"cut",col=col1[i],cex=0.8)
        if (class[i]=="CCAll") mtext(side=1,at=x1[i],line=1,"All",col=col1[i])
        if (class[i]=="CCConif") mtext(side=1,at=x1[i],line=1,"Conif",col=col1[i])
      } else {
        if (substr(class[i],nchar(class[i]),nchar(class[i])) %in% c("R","1","2","3","4","5","6","7","8")==TRUE) {
          if (substr(class[i],nchar(class[i]),nchar(class[i]))=="R") mtext(side=1,line=0,at=x1[i]-0.5,"0",col=col1[i],cex=0.8)
          if (substr(class[i],nchar(class[i]),nchar(class[i]))=="2") mtext(side=1,line=0,at=x1[i]-0.5,"20",col=col1[i],cex=0.8)
          if (substr(class[i],nchar(class[i]),nchar(class[i]))=="4") {
            mtext(side=1,line=0,at=x1[i]-0.5,"60",col=col1[i],cex=0.8)
            mtext(side=1,line=1,at=x1[i],substr(class[i],1,nchar(class[i])-1),col=col1[i],cex=1.2)
          }
          if (substr(class[i],nchar(class[i]),nchar(class[i]))=="6") mtext(side=1,line=0,at=x1[i]-0.5,"100",col=col1[i],cex=0.8)
          if (substr(class[i],nchar(class[i]),nchar(class[i]))=="8") mtext(side=1,line=0,at=x1[i]-0.5,"140",col=col1[i],cex=0.8)
        } else {
          if (class[i]=="Crop") class[i]<-"Crop+"  # Because this also includes other alienating
          mtext(side=1,at=x1[i],line=1,adj=ifelse(w[i]>2 | length(y)<9,0.5,1),las=ifelse(w[i]>2 | length(y)<9,1,2),class[i],col=col1[i],cex=1.2)
        }  # End if for aged treed type
      } # End if for CC
    }
    mtext(side=3,at=x1[1],adj=0,paste(sp.names.s[sp],"- North"),col="grey30",cex=1.2)
    text(max(x1),ymax*0.98,paste("Detected at",sum(sign(d.sp$Count)),"of",nrow(d.sp),"Summer camera locations"),cex=1.1,adj=1,col="grey40") # Add sample size
    graphics.off()
    sp.s<-sp
  }  # End if for that species being in the Summer species table

  # WINTER
  if (paste(SpTable.all[sp2],"Winter",sep="") %in% SpTable.w == TRUE) {
    sp<-which(SpTable.w==paste(SpTable.all[sp2],"Winter",sep=""))
    print(paste(sp,length(SpTable.w),SpTable.w[sp],date()))
    d.sp<-d[,c(1:(FirstSpCol.s-1),(LastSpCol.w+1):ncol(d),which(colnames(d)==SpTable.w[sp]))]  # Extract site descriptors and just the target species.  Assumes species columns are ASummer, AWinter...ZSummer, ZWinter
    colnames(d.sp)[ncol(d.sp)]="Count"  # Change the species count column name to "Count"
    d.sp<-d.sp[!is.na(d.sp$Count) & d.sp$WinterDays>=10,]  # Omit deployments with no sampling in that season or too few days
    d.sp$Lured<-as.character(d.sp$Lured)
    q<-by(sign(d.sp$Count[d.sp$Study=="ABMI"]),d.sp$Lured[d.sp$Study=="ABMI"],mean)
    lure.pa.w[sp2]<-q["y"]/q["n"]

    # W0. Simple summary of density in each veg+HF type that has been sampled
    lure1<-mean(d.sp$Count[d.sp$Study=="ABMI" & d.sp$Lured=="y"])/mean(d.sp$Count[d.sp$Study=="ABMI" & d.sp$Lured=="n"])  # Lure effect on total abundance
    for (i in 1:length(vnames.sum)) {
      x<-d.sp$Count[d.sp[,vnames.sum[i]]==1]
      l<-d.sp$Lured[d.sp[,vnames.sum[i]]==1]
      x<-ifelse(l=="y",x/lure1,x)
      simple.sum.mean.w[sp,i]<-mean(x)
      simple.sum.n.w[sp,i]<-sum(sign(x))
      if (n.sum[i]>1) simple.sum.se.w[sp,i]<-sd(x)/sqrt(n.sum[sp])
    }

    # W1. Fit models to broad veg types and HF types - presence/absence
    # W1.1 Fit broad veg models
    m.pa.w<-list(NULL)
    pCount<-sign(d.sp$Count)/ifelse(d.sp$Lured=="y",lure.pa.w[sp2],1)  # Standardize all to no-lure
    pCount<-pCount/max(pCount)  # This is in case the lure effect is <1
    d.sp$pCount.pa<-pCount  # Used in the age models below
    m.pa.w[[1]]<-try(glm(pCount~Decid+Mixedwood+Pine+Spruce+TreedBog+TreedFen+TreedSwamp+GrassHerb+Shrub+Marsh+ShrubbySwamp+ShrubbyBogFen+CCDecidR+CCDecid1+CCDecid2+CCMixedwoodR+CCMixedwood1+CCMixedwood2+CCPineR+CCPine1+CCSpruceR+CCSpruce1+CCSpruce2+SoftLin+TamePasture+RoughPasture+Wells+RuralResInd+WetlandMargin+WinterDays,family="binomial",data=d.sp,weights=d.sp$wt.w))  # Intercept is crop. No UrbInd sampled. HardLin assumed 0 because only 2 sampled.
    m.pa.w[[2]]<-try(glm(pCount~Decid+Mixedwood+Pine+Spruce+TreedBog+TreedFen+TreedSwamp+GrassHerb+Shrub+Marsh+ShrubbySwamp+ShrubbyBogFen+CCDecidMixed+CCPine+CCSpruce+SoftLin+TamePasture+RoughPasture+Wells+RuralResInd+WetlandMargin+WinterDays,family="binomial",data=d.sp,weights=d.sp$wt.w))  # Intercept is crop. No UrbInd sampled
    m.pa.w[[3]]<-try(glm(pCount~Decid+Mixedwood+Pine+Spruce+TreedBog+TreedFen+TreedSwamp+GrassHerb+Shrub+Marsh+ShrubbySwamp+ShrubbyBogFen+CCDecidMixed+CCPine+CCSpruce+SoftLin+WetlandMargin+WinterDays,family="binomial",data=d.sp,weights=d.sp$wt.w))  # Intercept is alien
    m.pa.w[[4]]<-try(glm(pCount~Decid+Mixedwood+Pine+Spruce+TreedBog+TreedFen+TreedSwamp+GrassHerb+Shrub+Marsh+ShrubbySwamp+ShrubbyBogFen+CCAll+SoftLin+WetlandMargin+WinterDays,family="binomial",data=d.sp,weights=d.sp$wt.w))  # Intercept is alien.  Weak data for cutblocks
    m.pa.w[[5]]<-try(glm(pCount~Decid+Mixedwood+Pine+Spruce+TreedBogFen+TreedSwamp+GrassHerb+Shrub+OpenWet+CCAll+SoftLin+WetlandMargin+WinterDays,family="binomial",data=d.sp,weights=d.sp$wt.w))  # Intercept is alien.
    m.pa.w[[6]]<-try(glm(pCount~Decid+Mixedwood+Pine+Spruce+TreedWet+GrassShrub+OpenWet++CCDecidR+CCDecid1+CCDecid2+CCMixedwoodR+CCMixedwood1+CCMixedwood2+CCPineR+CCPine1+CCSpruceR+CCSpruce1+CCSpruce2+SoftLin+WetlandMargin+WinterDays,family="binomial",data=d.sp,weights=d.sp$wt.w))  # Intercept is alien
    m.pa.w[[7]]<-try(glm(pCount~Decid+Mixedwood+Pine+Spruce+TreedWet+GrassShrub+OpenWet+CCDecidMixed+CCConif+SoftLin+WetlandMargin+WinterDays,family="binomial",data=d.sp,weights=d.sp$wt.w))  # Intercept is alien
    m.pa.w[[8]]<-try(glm(pCount~DecidMixed+UpCon+TreedWet+GrassShrub+OpenWet+CCDecidMixed+CCConif+SoftLin+WetlandMargin+WinterDays,family="binomial",data=d.sp,weights=d.sp$wt.w))  # Intercept is alien
    m.pa.w[[9]]<-try(glm(pCount~DecidMixed+UpCon+TreedWet+GrassShrub+OpenWet+CCAll+SoftLin+TamePasture+RoughPasture+Wells+RuralResInd+WetlandMargin+WinterDays,family="binomial",data=d.sp,weights=d.sp$wt.w))  # Intercept is crop
    m.pa.w[[10]]<-try(glm(pCount~DecidMixed+UpCon+TreedWet+GrassShrub+OpenWet+CCAll+SoftLin+WetlandMargin+WinterDays,family="binomial",data=d.sp,weights=d.sp$wt.w))  # Intercept is alien
    m.pa.w[[11]]<-try(glm(pCount~DecidMixed+UpCon+TreedWet+GrassShrub+OpenWet+Succ+WetlandMargin+WinterDays,family="binomial",data=d.sp,weights=d.sp$wt.w))  # Intercept is alien
    m.pa.w[[12]]<-try(glm(pCount~Upland+Lowland+CCAll+SoftLin+WetlandMargin+WinterDays,family="binomial",data=d.sp,weights=d.sp$wt.w))  # Intercept is alien
    m.pa.w[[13]]<-try(glm(pCount~TreedAll+OpenAll+Succ+WetlandMargin+WinterDays,family="binomial",data=d.sp,weights=d.sp$wt.w))  # Intercept is alien
    m.pa.w[[14]]<-try(glm(pCount~UplandForest+Lowland+GrassShrub+Succ+WetlandMargin+WinterDays,family="binomial",data=d.sp,weights=d.sp$wt.w))  # Intercept is alien
    m.pa.w[[15]]<-try(glm(pCount~Upland+Lowland+Succ+WetlandMargin+WinterDays,family="binomial",data=d.sp,weights=d.sp$wt.w))  # Intercept is alien
    m.pa.w[[16]]<-try(glm(pCount~Boreal+GrassShrub+Succ+WetlandMargin+WinterDays,family="binomial",data=d.sp,weights=d.sp$wt.w))  # "Boreal" is everything native except GrassShrub. Intercept is alien
    # AIC calculation  (I'm using AIC here, because this is primarily for prediction, rather than finding a minimal best model)
    nModels.pa<-length(m.pa.w)
    aic.ta<-rep(999999999,(nModels.pa))
    for (i in 1:(nModels.pa)) {
      if (!is.null(m.pa.w[[i]]) & class(m.pa.w[[i]])[1]!="try-error") {  # last part is to not used non-converged models, unless none converged
        aic.ta[i]<-AICc(m.pa.w[[i]])
      }
    }
    aic.delta<-aic.ta-min(aic.ta)
    aic.exp<-exp(-1/2*aic.delta)
    aic.wt.pa<-aic.exp/sum(aic.exp)
    best.model.pa<-which.max(aic.wt.pa)
    aic.wt.pa.save.w[sp2,]<-aic.wt.pa

    # W1.2 Then do abundance | presence model
    d.p<-d.sp[d.sp$Count>0,]  # Just presence records
    q<-by(d.p$Count[d.sp$Study=="ABMI"],d.p$Lured[d.sp$Study=="ABMI"],mean)
    lure.agp.w[sp2]<-q["y"]/q["n"]
    pCount<-d.p$Count/ifelse(d.p$Lured=="y",lure.agp.w[sp2],1)
    j<-0
    m.agp.w<-list(NULL)
    mnums<-NULL  # To keep track of which models were actually fit
    for (i in 1:nModels.pa) {
      x<-colSums(d.p[,attr(m.pa.w[[i]]$terms,"term.labels")])  # Minimum number of presence records in each veg+HF type included in the model
      if (min(x)>3 & nrow(d.p)-sum(x[-length(x)])) {  # Only use the model if each type is represented by >3 presence records (incl. intercept)
        j<-j+1
        m.agp.w[[j]]<-glm(m.pa.w[[i]]$formula,data=d.p,family=gaussian(link="log"),weights=d.p$wt)
        mnums<-c(mnums,i)
      }
    }
    m.agp.w[[j+1]]<-glm(pCount~WinterDays,data=d.p,family=gaussian(link="log"),weights=d.p$wt)  # And add a null
    mnums<-c(mnums,nModels.pa+1)  # Add the null model to list - 1 beyond the last pa model
    # AIC calculation  (I'm using AIC here, because this is primarily for prediction, rather than finding a minimal best model)
    nModels.agp<-length(m.agp.w)
    aic.ta<-rep(999999999,(nModels.agp))
    for (i in 1:(nModels.agp)) {
      if (!is.null(m.agp.w[[i]]) & class(m.agp.w[[i]])[1]!="try-error") {  # last part is to not used non-converged models, unless none converged
        aic.ta[i]<-AICc(m.agp.w[[i]])
      }
    }
    aic.delta<-aic.ta-min(aic.ta)
    aic.exp<-exp(-1/2*aic.delta)
    aic.wt.agp<-aic.exp/sum(aic.exp)
    best.model.agp<-which.max(aic.wt.agp)
    aic.wt.agp.save.w[sp2,mnums]<-aic.wt.agp

    # W1.3 Predict from models for 100% of each veg type and for each site - pres/abs - using only the best model
    Intercept<-c("Crop","Crop",rep("Alien",6),"Crop",rep("Alien",7))[best.model.pa]  # Change this if models change
    terms.pa<-c(attr(m.pa.w[[best.model.pa]]$terms,"term.labels"),Intercept)
    Coef1.pa<-Coef1.pa.se<-rep(NA,length(terms.pa))
    names(Coef1.pa)<-names(Coef1.pa.se)<-terms.pa
    for (i in 1:length(terms.pa)) {
      pm1<-rep(0,length(terms.pa))
      pm1[i]<-1
      names(pm1)<-terms.pa
      pm1<-data.frame(t(pm1))
      p<-predict(m.pa.w[[best.model.pa]],newdata=data.frame(WinterDays=100,pm1),se.fit=T)
      Coef1.pa[i]<-plogis(p$fit)  # Ordinal scale
      Coef1.pa.se[i]<-p$se.fit  # logit scale
    }
    # Adjust so that mean prediction at standardized WinterDays across all qualifying sites = mean observed count.  This is to compensate for inaccuracies due to fitting WinterDays coefficients
    pCount<-sign(d.sp$Count)/ifelse(d.sp$Lured=="y",lure.pa.w[sp2],1)  # Standardize all to no-lure
    pCount<-pCount/max(pCount)  # This is in case the lure effect is <1
    p.adj<-predict(m.pa.w[[best.model.pa]],newdata=data.frame(WinterDays=100,d.sp))  # Prediction for each data point, but at standardized days
    Coef1.pa<-plogis(qlogis(Coef1.pa)+qlogis(mean(pCount))-qlogis(mean(plogis(p.adj))))  # Adjustment made on logit scale.  Var doesn't change on logit scale(?)
    # And predict for each site - used for age model below
    d.sp$p<-predict(m.pa.w[[best.model.pa]])  # Logit scale

    # W1.4 Predict from models for 100% of each veg type and for each site - abund|pres
    terms.pa1<-terms.pa[terms.pa!="WinterDays"]
    Coef1.agp<-Coef1.agp.se<-rep(NA,length(terms.pa1))  # For agp coefficients averaged to terms in best pa model
    names(Coef1.agp)<-names(Coef1.agp.se)<-terms.pa1
    i<-best.model.agp
    p<-predict(m.agp.w[[i]],newdata=data.frame(WinterDays=100,pm),se.fit=T)  # For each type
    tTypeMean.agp<-p$fit
    tTypeVar.agp<-p$se.fit^2  # Using variance, for consistency with other scripts at this point
    names(tTypeMean.agp)<-names(tTypeVar.agp)<-pm$VegType
    # Adjust so that mean prediction at standardized WinterDays across all qualifying sites = mean observed count.  This is to compensate for inaccuracies due to fitting WinterDays coefficients
    pCount<-d.p$Count/ifelse(d.p$Lured=="y",lure.agp.w[sp2],1)
    p.adj<-predict(m.agp.w[[i]],newdata=data.frame(WinterDays=100,d.sp))  # Prediction for each data point, but at standardized days
    tTypeMean.agp<-log(exp(tTypeMean.agp)+mean(pCount)-mean(exp(p.adj)))  # Var doesn't change on log scale
    # Then average those for each broader group included in the best pa model
    for (i in 1:length(terms.pa1)) {
      j<-pm$VegType[which(pm[,terms.pa1[i]]==1)]  # Names of fine hab+HF types included in that broader group
      x<-tTypeMean.agp[as.character(j)]
      x.var<-tTypeVar.agp[as.character(j)]
      Coef1.agp[i]<-mean(x)  # Simple mean, log scale
      Coef1.agp.se[i]<-sqrt(mean(x.var))  # This is for straight-up mean, log scale
    }
    # Put non-stand-age coefficients in Coef matrix (and SE's)
    Coef.agp.w[sp,]<-exp(tTypeMean.agp[na.omit(match(colnames(Coef.agp.w),names(tTypeMean.agp)))])  # On ordinal scale
    Coef.agp.w.se[sp,]<-sqrt(tTypeVar.agp[na.omit(match(colnames(Coef.agp.w),names(tTypeMean.agp)))])  # On log scale

    # W2. Run models for age within each stand type - pres/abs only currently
    # W2.1. Set up separate dataframes for sites containing a minimum amount of each broad stand type
    d.spruce<-d.pine<-d.decid<-d.mixed<-d.treedbogfen<-NULL
    for (i in 0:8) {  # Add sites with each age class of the stand type to a separate data frames for each stand type
      i1<-ifelse(i==0,"R",i)  # For variable name
      i2<-ifelse(i==0,0.5,i)  # For twenty-year age
      cn<-paste("Spruce",i1,sep="")  # Col name
      if(cn %in% names(d.sp)) if (sum(d.sp[,cn]>cutoff.pc.for.age)>0) d.spruce<-rbind(d.spruce,data.frame(pCount=d.sp[d.sp[,cn]>cutoff.pc.for.age,"pCount.pa"], age=i2, wt1=d.sp[d.sp[,cn]>cutoff.pc.for.age,cn]*d.sp$wt.w[d.sp[,cn]>cutoff.pc.for.age], p=d.sp$p[d.sp[,cn]>cutoff.pc.for.age]))  # Weight is the ppoportion of the site of that age class and stand type, multiplied by the original weight (which accounts for revisited sites)
      cn<-paste("Pine",i1,sep="")  # Col name
      if(cn %in% names(d.sp)) if (sum(d.sp[,cn]>cutoff.pc.for.age)>0) d.pine<-rbind(d.pine,data.frame(pCount=d.sp[d.sp[,cn]>cutoff.pc.for.age,"pCount.pa"], age=i2, wt1=d.sp[d.sp[,cn]>cutoff.pc.for.age,cn]*d.sp$wt.w[d.sp[,cn]>cutoff.pc.for.age], p=d.sp$p[d.sp[,cn]>cutoff.pc.for.age]))
      cn<-paste("Decid",i1,sep="")  # Col name
      if(cn %in% names(d.sp)) if (sum(d.sp[,cn]>cutoff.pc.for.age)>0) d.decid<-rbind(d.decid,data.frame(pCount=d.sp[d.sp[,cn]>cutoff.pc.for.age,"pCount.pa"], age=i2, wt1=d.sp[d.sp[,cn]>cutoff.pc.for.age,cn]*d.sp$wt.w[d.sp[,cn]>cutoff.pc.for.age], p=d.sp$p[d.sp[,cn]>cutoff.pc.for.age]))
      cn<-paste("Mixedwood",i1,sep="")  # Col name
      if(cn %in% names(d.sp)) if (sum(d.sp[,cn]>cutoff.pc.for.age)>0) d.mixed<-rbind(d.mixed,data.frame(pCount=d.sp[d.sp[,cn]>cutoff.pc.for.age,"pCount.pa"], age=i2, wt1=d.sp[d.sp[,cn]>cutoff.pc.for.age,cn]*d.sp$wt.w[d.sp[,cn]>cutoff.pc.for.age], p=d.sp$p[d.sp[,cn]>cutoff.pc.for.age]))
      cn<-paste("TreedBog",i1,sep="")  # Col name
      if(cn %in% names(d.sp)) if (sum(d.sp[,cn]>cutoff.pc.for.age)>0) d.treedbog<-rbind(d.treedbog,data.frame(pCount=d.sp[d.sp[,cn]>cutoff.pc.for.age,"pCount.pa"], age=i2, wt1=d.sp[d.sp[,cn]>cutoff.pc.for.age,cn]*d.sp$wt.w[d.sp[,cn]>cutoff.pc.for.age], p=d.sp$p[d.sp[,cn]>cutoff.pc.for.age]))
    }
    # Note: Not using grass and shrub here, because including parklands, so most of these are prairies, not recent fires
    # Combined data frames for models using more than one stand type
    d.upcon<-rbind(cbind(d.spruce,Sp="Spruce"),cbind(d.pine,Sp="Pine"))
    d.decidmixed<-rbind(cbind(d.decid,Sp="Decid"),cbind(d.mixed,Sp="Mixed"))
    d.all<-rbind(cbind(d.spruce,Sp="Spruce"),cbind(d.pine,Sp="Pine"),cbind(d.decid,Sp="Decid"),cbind(d.mixed,Sp="Mixed"),cbind(d.treedbog,Sp="TreedBog"))
    # 2.2. Then fit models of age functions
    # Note: Days not included in these models, because its effect on presence/absence is already included in the prediction p used in the offset
    m.age.w<-m.null.age.w<-list(NULL)  # Age models for each of the 5 stand types plus 3 combinations.  Null model is for AIC comparison to smoothing spline model (which can be worse, even if it has >1 equiv DF)
    m.null.age.w[[1]]<-gam(pCount~1+offset(p),data=d.spruce,family="binomial",weights=d.spruce$wt1)
    if (sum(sign(d.spruce$pCount))>4) {
      m.age.w[[1]]<-gam(pCount~s(sqrt(age),k=4,m=2)+offset(p),data=d.spruce,family="binomial",weights=d.spruce$wt1)  # Fit spline through age data for that stand type if enough records
    } else {
      m.age.w[[1]]<-m.null.age.w[[1]]  # Use constant-only if too few records in that stand type
    }
    m.null.age.w[[2]]<-gam(pCount~1+offset(p),data=d.pine,family="binomial",weights=d.pine$wt1)
    if (sum(sign(d.pine$pCount))>4) {
      m.age.w[[2]]<-gam(pCount~s(sqrt(age),k=4,m=2)+offset(p),data=d.pine,family="binomial",weights=d.pine$wt1)
    } else {
      m.age.w[[2]]<-m.null.age.w[[2]]  # Use constant-only if too few records in that stand type
    }
    m.null.age.w[[3]]<-gam(pCount~1+offset(p),data=d.decid,family="binomial",weights=d.decid$wt1)
    if (sum(sign(d.decid$pCount))>4) {
      m.age.w[[3]]<-gam(pCount~s(sqrt(age),k=4,m=2)+offset(p),data=d.decid,family="binomial",weights=d.decid$wt1)
    } else {
      m.age.w[[3]]<-m.null.age.w[[3]]  # Use constant-only if too few records in that stand type
    }
    m.null.age.w[[4]]<-gam(pCount~1+offset(p),data=d.mixed,family="binomial",weights=d.mixed$wt1)
    if (sum(sign(d.mixed$pCount))>4) {
      m.age.w[[4]]<-gam(pCount~s(sqrt(age),k=4,m=2)+offset(p),data=d.mixed,family="binomial",weights=d.mixed$wt1)
    } else {
      m.age.w[[4]]<-m.null.age.w[[4]]  # Use constant-only if too few records in that stand type
    }
    m.null.age.w[[5]]<-gam(pCount~1+offset(p),data=d.treedbog,family="binomial",weights=d.treedbog$wt1)
    if (sum(sign(d.treedbog$pCount))>4) {
      m.age.w[[5]]<-gam(pCount~s(sqrt(age),k=4,m=2)+offset(p),data=d.treedbog,family="binomial",weights=d.treedbog$wt1)
    } else {
      m.age.w[[5]]<-m.null.age.w[[5]]  # Use constant-only if too few records in that stand type
    }
    m.null.age.w[[6]]<-gam(pCount~1+offset(p),data=d.upcon,family="binomial",weights=d.upcon$wt1)
    if (sum(sign(d.upcon$pCount))>4) {
      m.age.w[[6]]<-gam(pCount~s(sqrt(age),k=4,m=2)+offset(p),data=d.upcon,family="binomial",weights=d.upcon$wt1)
    } else {
      m.age.w[[6]]<-m.null.age.w[[6]]  # Use constant-only if too few records in that stand type
    }
    m.null.age.w[[7]]<-gam(pCount~1+offset(p),data=d.decidmixed,family="binomial",weights=d.decidmixed$wt1)
    if (sum(sign(d.decidmixed$pCount))>4) {
      m.age.w[[7]]<-gam(pCount~s(sqrt(age),k=4,m=2)+offset(p),data=d.decidmixed,family="binomial",weights=d.decidmixed$wt1)
    } else {
      m.age.w[[7]]<-m.null.age.w[[7]]  # Use constant-only if too few records in that stand type
    }
    m.null.age.w[[8]]<-gam(pCount~1+offset(p),data=d.all,family="binomial",weights=d.all$wt1)
    if (sum(sign(d.all$pCount))>4) {
      m.age.w[[8]]<-gam(pCount~s(sqrt(age),k=4,m=2)+offset(p),data=d.all,family="binomial",weights=d.all$wt1)
    } else {
      m.age.w[[8]]<-m.null.age.w[[8]]  # Use constant-only if too few records in that stand type
    }
    # AIC for options - each separate, pairwise combinations, combine all
    aic.age.w[sp,1]<-AIC(m.age.w[[1]])+AIC(m.age.w[[2]])+AIC(m.age.w[[3]])+AIC(m.age.w[[4]])+AIC(m.age.w[[5]])-8  # Minus 8 for the 4 extra variances
    aic.age.w[sp,2]<-AIC(m.age.w[[5]])+AIC(m.age.w[[6]])+AIC(m.age.w[[7]])-4
    aic.age.w[sp,3]<-AIC(m.age.w[[8]])
    aic.delta<-aic.age.w[sp,]-min(aic.age.w[sp,])
    aic.exp<-exp(-1/2*aic.delta)
    aic.wt.age<-aic.exp/sum(aic.exp)
    aic.wt.age.models<-NULL
    for (i in 1:8) {  # AIC weight for spline (versus null) for each age grouping
      q<-c(AICc(m.age.w[[i]]),AICc(m.null.age.w[[i]]))
      q<-q-min(q)
      aic.wt.age.models[i]<-exp(-1/2*q[1])/sum(exp(-1/2*q))
    }
    names(aic.wt.age.models)<-c("Spruce","Pine","Decid","Mixedwood","TreedBog","DecidMixed","UpCon","TreedAll")
    # W2.3 Then make age predictions only if the spline model is better than the null, and only for (groupings of) stand types in the pa model
    age.flag<-rep(0,10)  # Flag for separate age predictions for each stand type {Spruce, Pine, Deciduous, Mixedwood, TreedBog, TreedWet (includes TreedSwamp), DecidMixed, UpCon, UplandForest, TreedAll) (0=no, 1=yes)
    p.age<-p.age.se<-array(NA,c(10,9))  # Predictions for 10 stand types, 9 age classes
    rownames(p.age)<-rownames(p.age.se)<-names(age.flag)<-c("Spruce","Pine","Decid","Mixedwood","TreedBog","TreedWet","DecidMixed","UpCon","UplandForest","TreedAll")
    if ("Spruce" %in% terms.pa) {  # Do spruce separately
      if (aic.wt.age.models["Spruce"]>0.5) {
        age.flag[1]<-1
        p<-predict(m.age.w[[1]],newdata=data.frame(age=c(0.5,1:8),p=qlogis(Coef1.pa["Spruce"])),se.fit=T)
        p.age["Spruce",]<-plogis(p$fit)
        p.age.se["Spruce",]<-p$se.fit  # Still on link scale
      }
    }
    if ("Pine" %in% terms.pa) {  # Do pine separately
      if (aic.wt.age.models["Pine"]>0.5) {
        age.flag["Pine"]<-1
        p<-predict(m.age.w[[2]],newdata=data.frame(age=c(0.5,1:8),p=qlogis(Coef1.pa["Pine"])),se.fit=T)
        p.age["Pine",]<-plogis(p$fit)
        p.age.se["Pine",]<-p$se.fit  # Still on link scale
      }
    }
    if ("Decid" %in% terms.pa) {  # Do deciduous separately
      if (aic.wt.age.models["Decid"]>0.5) {
        age.flag["Decid"]<-1
        p<-predict(m.age.w[[3]],newdata=data.frame(age=c(0.5,1:8),p=qlogis(Coef1.pa["Decid"])),se.fit=T)
        p.age["Decid",]<-plogis(p$fit)
        p.age.se["Decid",]<-p$se.fit  # Still on link scale
      }
    }
    if ("Mixedwood" %in% terms.pa) {  # Do mixedwood separately
      if (aic.wt.age.models["Mixedwood"]>0.5) {
        age.flag["Mixedwood"]<-1
        p<-predict(m.age.w[[4]],newdata=data.frame(age=c(0.5,1:8),p=qlogis(Coef1.pa["Mixedwood"])),se.fit=T)
        p.age["Mixedwood",]<-plogis(p$fit)
        p.age.se["Mixedwood",]<-p$se.fit  # Still on link scale
      }
    }
    if ("TreedBog" %in% terms.pa) {  # Do TreedBog separately
      if (aic.wt.age.models["TreedBog"]>0.5) {
        age.flag["TreedBog"]<-1
        p<-predict(m.age.w[[5]],newdata=data.frame(age=c(0.5,1:8),p=qlogis(Coef1.pa["TreedBog"])),se.fit=T)
        p.age["TreedBog",]<-plogis(p$fit)
        p.age.se["TreedBog",]<-p$se.fit  # Still on link scale
      }
    }
    if ("TreedWet" %in% terms.pa) {  # Do TreedWet separately, using TreedBog model
      if (aic.wt.age.models["TreedBog"]>0.5) {
        age.flag["TreedWet"]<-1
        p<-predict(m.age.w[[5]],newdata=data.frame(age=c(0.5,1:8),p=qlogis(Coef1.pa["TreedWet"])),se.fit=T)
        p.age["TreedWet",]<-plogis(p$fit)
        p.age.se["TreedWet",]<-p$se.fit  # Still on link scale
      }
    }
    if ("UpCon" %in% terms.pa) {  # Do spruce+pine together
      if (aic.wt.age.models["UpCon"]>0.5) {
        age.flag["UpCon"]<-1
        p<-predict(m.age.w[[6]],newdata=data.frame(age=c(0.5,1:8),p=qlogis(Coef1.pa["UpCon"])),se.fit=T)
        p.age["UpCon",]<-plogis(p$fit)
        p.age.se["UpCon",]<-p$se.fit  # Still on link scale
      }
    }
    if ("DecidMixed" %in% terms.pa) {  # Do deciduous+mixed together
      if (aic.wt.age.models["DecidMixed"]>0.5) {
        age.flag["DecidMixed"]<-1
        p<-predict(m.age.w[[7]],newdata=data.frame(age=c(0.5,1:8),p=qlogis(Coef1.pa["DecidMixed"])),se.fit=T)
        p.age["DecidMixed",]<-plogis(p$fit)
        p.age.se["DecidMixed",]<-p$se.fit  # Still on link scale
      }
    }
    if ("UplandForest" %in% terms.pa)  { # Do those four stand types together, using "all" age model
      if (aic.wt.age.models["TreedAll"]>0.5) {
        age.flag["UplandForest"]<-1
        p<-predict(m.age.w[[8]],newdata=data.frame(age=c(0.5,1:8),p=qlogis(Coef1.pa["UplandForest"])),se.fit=T)
        p.age["UplandForest",]<-plogis(p$fit)
        p.age.se["UplandForest",]<-p$se.fit  # Still on link scale
      }
    }
    if ("TreedAll" %in% terms.pa)  {  # Do all five together
      if (aic.wt.age.models["TreedAll"]>0.5) {
        age.flag["TreedAll"]<-1
        p<-predict(m.age.w[[8]],newdata=data.frame(age=c(0.5,1:8),p=qlogis(Coef1.pa["TreedAll"])),se.fit=T)
        p.age["TreedAll",]<-plogis(p$fit)
        p.age.se["TreedAll",]<-p$se.fit  # Still on link scale
      }
    }
    # Save aic weights
    aic.wt.age.save.w[sp2,]<-aic.wt.age
    aic.wt.age.models.save.w[sp2,]<-aic.wt.age.models

    # W3. Assemble presence/absence coefficients for each veg+HF type, including expanding out stand types that have ages, to populate entire Coef.pa matrix, then multiply by appropriate Coef.agp term to generate Coef.mean matrix
    for (i in 1:length(terms.pa1)) {
      j<-as.character(pm$VegType[pm[,terms.pa1[i]]==1])  # The fine veg+HF types that are covered by the (potentially broader) variable included in the best model
      for (j1 in 1:length(j)) {
        if (j[j1] %in% c("Spruce","Pine","Decid","Mixedwood","TreedBog")) {
          # Stand types that have age classes in full coefficients
          x<-paste(j[j1],c("R","1","2","3","4","5","6","7","8"),sep="")  # Names in Coef.pa and Coef.mean with ages
          if ((age.flag[terms.pa1[i]]==0)==TRUE | is.na(age.flag[terms.pa1[i]])) {   # But these broad types do not have age classes in the best model
            Coef.pa.w[sp,x]<-Coef1.pa[terms.pa1[i]]  # Fill in all ages with value for that (broad) stand type
            Coef.pa.w.se[sp,x]<-Coef1.pa.se[terms.pa1[i]]  # Fill in all ages with value for that (broad) stand type
            Coef.mean.w[sp,x]<-Coef.pa.w[sp,x]*Coef.agp.w[sp,j[j1]]  # Only the one value for agp
            Coef.lci.w[sp,x]<-plogis(qlogis(Coef.pa.w[sp,x])-Coef.pa.w.se[sp,x]*1.28) * exp(log(Coef.agp.w[sp,j[j1]])-Coef.agp.w.se[sp,j[j1]]*1.28)  # Using 10% intervals for each to multiply to 5% intervals (assuming independence - checked empirically)
            Coef.uci.w[sp,x]<-plogis(qlogis(Coef.pa.w[sp,x])+Coef.pa.w.se[sp,x]*1.28) * exp(log(Coef.agp.w[sp,j[j1]])+Coef.agp.w.se[sp,j[j1]]*1.28)  # Using 10% intervals for each to multiply to 5% intervals (assuming independence - checked empirically)
          } else {  # These broad types do have age classes in the best model
            Coef.pa.w[sp,x]<-p.age[terms.pa1[i],]  # Use separate age values for that (broad) stand type
            Coef.pa.w.se[sp,x]<-p.age.se[terms.pa1[i],]  # Use separate age values for that (broad) stand type
            Coef.mean.w[sp,x]<-Coef.pa.w[sp,x]*Coef.agp.w[sp,j[j1]]  # Only the one value for agp
            Coef.lci.w[sp,x]<-plogis(qlogis(Coef.pa.w[sp,x])-Coef.pa.w.se[sp,x]*1.28) * exp(log(Coef.agp.w[sp,j[j1]])-Coef.agp.w.se[sp,j[j1]]*1.28)  # Using 10% intervals for each to multiply to 5% intervals (assuming independence - checked empirically)
            Coef.uci.w[sp,x]<-plogis(qlogis(Coef.pa.w[sp,x])+Coef.pa.w.se[sp,x]*1.28) * exp(log(Coef.agp.w[sp,j[j1]])+Coef.agp.w.se[sp,j[j1]]*1.28)  # Using 10% intervals for each to multiply to 5% intervals (assuming independence - checked empirically)
          }
        } else {
          # Stand types that do not have age classes in full coefficients
          Coef.pa.w[sp,j[j1]]<-Coef1.pa[terms.pa1[i]]  # Fill in the non-age coefficient
          Coef.pa.w.se[sp,j[j1]]<-Coef1.pa.se[terms.pa1[i]]  # Fill in the non-age coefficient
          Coef.mean.w[sp,j[j1]]<-Coef.pa.w[sp,j[j1]]*Coef.agp.w[sp,j[j1]]
          Coef.lci.w[sp,j[j1]]<-plogis(qlogis(Coef.pa.w[sp,j[j1]])-Coef.pa.w.se[sp,j[j1]]*1.28) * exp(log(Coef.agp.w[sp,j[j1]])-Coef.agp.w.se[sp,j[j1]]*1.28)  # Using 10% intervals for each to multiply to 5% intervals (assuming independence - checked empirically)
          Coef.uci.w[sp,j[j1]]<-plogis(qlogis(Coef.pa.w[sp,j[j1]])+Coef.pa.w.se[sp,j[j1]]*1.28) * exp(log(Coef.agp.w[sp,j[j1]])+Coef.agp.w.se[sp,j[j1]]*1.28)  # Using 10% intervals for each to multiply to 5% intervals (assuming independence - checked empirically)
        }  # End if for stand types with age classes
      }  # Next fine veg+HF type within broader class
    }  # Next broader class in best model
    # Further adjustments so that mean(Coef.mean) = mean(observed density).  Adjusts for geometric mean issue, and also for any remaining issues with Days adjustments
    pCount<-d.sp$Count/(lure.pa.w[sp2]*lure.agp.w[sp2])
    z<-colnames(Coef.mean.w)[!is.na(match(colnames(Coef.mean.w),colnames(d.sp)))]  # Coef names that are in d.sp - excludes added categories not in data, like old clearcuts
    z<-z[!is.na(Coef.mean.s[sp,z])]
    q<-colSums(t(d.sp[,z]) * Coef.mean.w[sp,z])  # Prediction for each deployment, based on the coefficients and its veg composition
    Coef.mean.w[sp,]<-Coef.mean.w[sp,]*mean(pCount)/mean(q)
    Coef.lci.w[sp,]<-Coef.lci.w[sp,]*mean(pCount)/mean(q)
    Coef.uci.w[sp,]<-Coef.uci.w[sp,]*mean(pCount)/mean(q)

    # W3.2 Not using cutblock convergence for best-model analysis (i.e., for "honest" coefficient figures), but it is used below for the final all-coefficients exported
    # Predictions for each site, as offsets below
    d.sp$t.p.pa<-colSums(Coef1.pa[terms.pa1]*t(d.sp[,terms.pa1]))  # Prediction of presence/absence at each site - ordinal scale
    d.sp$t.p.agp<-colSums(Coef1.agp[terms.pa1]*t(d.sp[,terms.pa1]))  # Prediction of abundance|presence at each site - log scale
    d.sp$p.ta<-d.sp$t.p.pa*exp(d.sp$t.p.agp)  # Prediction (ordinal scale) of total abundance at each site

    # W3.3 Store values for winter
    j<-match(d.sp$DeploymentYear,d1$DeploymentYear)
    d1$count.winter[j]<-d.sp$Count
    d1$p.pa.winter[j]<-d.sp$t.p.pa
    d1$p.agp.winter[j]<-exp(d.sp$t.p.agp)
    d1$p.winter[j]<-d.sp$p.ta

    # W4. Coefficient figures
    # W4.1. Combined total abundance
    # Custom figure for each best model and age options
    x<-col1<-y<-y.lci<-y.uci<-w<-class<-space1<-NULL
    for (i in 1:length(terms.pa1)) {
      if (age.flag[terms.pa1[i]]==0 | is.na(age.flag[terms.pa1[i]])) {  # Variable without age info
        j<-match(terms.pa1[i],fp$Class)
        x<-c(x,fp$x[j])
        veg.to.use<-as.character(pm$VegType[pm[,terms.pa1[i]]==1])[1]  # To use Coef.mean, find the first fine veg+HF types that is included in the (potentially broader) variable being plotted
        if (veg.to.use %in% c("Spruce","Pine","Decid","Mixedwood","TreedBog")) veg.to.use<-paste(veg.to.use,"R",sep="")  # If no age relationship, use the first value for that fine stand type (don't add the R for types that don't have age classes at all)
        y<-c(y,Coef.mean.w[sp,veg.to.use])
        y.lci<-c(y.lci,Coef.lci.w[sp,veg.to.use])
        y.uci<-c(y.uci,Coef.uci.w[sp,veg.to.use])
        col1<-c(col1,fp$col1[j])
        w<-c(w,fp$width[j])
        class<-c(class,terms.pa1[i])
        space1<-c(space1,fp$spaceafter[j])
      } else {  # Variable with age info
        j<-match(paste(terms.pa1[i],"R",sep=""),fp$Class):match(paste(terms.pa1[i],"8",sep=""),fp$Class)
        x<-c(x,fp$x[j])
        veg.to.use<-as.character(pm$VegType[pm[,terms.pa1[i]]==1])[1]  # To use Coef.mean, find the first fine veg+HF types that is included in the (potentially broader) variable being plotted
        veg.to.use<-paste(veg.to.use,c("R","1","2","3","4","5","6","7","8"),sep="")  # If age relationship, use each value for that fine stand type
        y<-c(y,Coef.mean.w[sp,veg.to.use])
        y.lci<-c(y.lci,Coef.lci.w[sp,veg.to.use])
        y.uci<-c(y.uci,Coef.uci.w[sp,veg.to.use])
        col1<-c(col1,fp$col1[j])
        w<-c(w,fp$width[j])
        class<-c(class,as.character(fp$Class[j]))
        space1<-c(space1,fp$spaceafter[j])
      }
    }
    ord<-order(x)  # Sort all by x
    y<-y[ord]
    y.lci<-y.lci[ord]
    y.uci<-y.uci[ord]
    col1<-col1[ord]
    w<-w[ord]
    class<-class[ord]
    space1<-space1[ord]
    x<-x[ord]
    # Rectify spaces between x's
    for (i in 1:(length(x)-1)) {
      for (j in (i+1):length(x)) x[j]<-x[j]+space1[i]-(x[i+1]-x[i])  # Alter all subsequent positions accordingly
    }
    # Make bar plot
    ymax<-min(max(y.uci,na.rm=TRUE),2*max(y,na.rm=TRUE))  # This keeps the figures readable when there are extreme UCI's
    space<-c(1,x[-1]-x[-length(x)])-0.99  # The spacing between bars
    density<-ifelse(substr(class,1,2)=="CC",50,NA)
    fname<-paste(fname.fig,"Veg+HF figure best model ",SpTable.w[sp],".png",sep="")
    png(file=fname,width=ifelse(length(y)>5,1500,1000),height=700)
    par(mai=c(1.9,1,0.2,0.3))
    x1<-barplot(y,space=space,width=w,border="white",col="grey30",ylim=c(0,ymax),yaxt="n",ylab="Relative abundance",col.lab="grey50",cex.lab=1.2,axisnames=F)[,1]  # To get strips on CC bars
    abline(h=pretty(c(0,ymax)),col="grey80")
    x1<-barplot(y,space=space,width=w,border="white",col="grey30",ylim=c(0,ymax),yaxt="n",ylab="Relative abundance",col.lab="grey50",cex.lab=1.2,axisnames=F,add=TRUE)[,1]  # To get strips on CC bars, and to put bars in front of horizontal axis lines
    x1<-barplot(y,space=space,width=w,border="white",density=density,col=col1,ylim=c(0,ymax),yaxt="n",ylab="Relative abundance",col.lab="grey50",cex.lab=1.2,axisnames=F,add=TRUE)[,1]
    axis(side=2,tck=0.02,cex.axis=0.9,col.axis="grey50",col.ticks="grey50",las=2,at=pretty(c(0,ymax)))
    box(bty="l",col="grey50")
    for (i in 1:length(x1)) {
      lines(rep(x1[i],2),c(y[i],y.uci[i]),col=col1[i])
      lines(rep(x1[i],2),c(y[i],y.lci[i]),col="grey90")
    }
    for (i in 1:length(class)) {  # Label x axis
      if (substr(class[i],1,2)=="CC") {
        mtext(side=1,line=0,at=x1[i],"cut",col=col1[i],cex=0.8)
        if (class[i]=="CCAll") mtext(side=1,at=x1[i],line=1,"All",col=col1[i])
        if (class[i]=="CCConif") mtext(side=1,at=x1[i],line=1,"Conif",col=col1[i])
      } else {
        if (substr(class[i],nchar(class[i]),nchar(class[i])) %in% c("R","1","2","3","4","5","6","7","8")==TRUE) {
          if (substr(class[i],nchar(class[i]),nchar(class[i]))=="R") mtext(side=1,line=0,at=x1[i]-0.5,"0",col=col1[i],cex=0.8)
          if (substr(class[i],nchar(class[i]),nchar(class[i]))=="2") mtext(side=1,line=0,at=x1[i]-0.5,"20",col=col1[i],cex=0.8)
          if (substr(class[i],nchar(class[i]),nchar(class[i]))=="4") {
            mtext(side=1,line=0,at=x1[i]-0.5,"60",col=col1[i],cex=0.8)
            mtext(side=1,line=1,at=x1[i],substr(class[i],1,nchar(class[i])-1),col=col1[i],cex=1.2)
          }
          if (substr(class[i],nchar(class[i]),nchar(class[i]))=="6") mtext(side=1,line=0,at=x1[i]-0.5,"100",col=col1[i],cex=0.8)
          if (substr(class[i],nchar(class[i]),nchar(class[i]))=="8") mtext(side=1,line=0,at=x1[i]-0.5,"140",col=col1[i],cex=0.8)
        } else {
          if (class[i]=="Crop") class[i]<-"Crop+"  # Because this also includes other alienating
          mtext(side=1,at=x1[i],line=1,adj=ifelse(w[i]>2 | length(y)<9,0.5,1),las=ifelse(w[i]>2 | length(y)<9,1,2),class[i],col=col1[i],cex=1.2)
        }  # End if for aged treed type
      } # End if for CC
    }
    mtext(side=3,at=x1[1],adj=0,paste(sp.names.w[sp],"- North"),col="grey30",cex=1.2)
    text(max(x1),ymax*0.98,paste("Detected at",sum(sign(d.sp$Count)),"of",nrow(d.sp),"Winter camera locations"),cex=1.1,adj=1,col="grey40") # Add sample size
    graphics.off()
    sp.w<-sp
  }  # End if for that species being in the winter species table

  # 5.1. Combine Summer and winter coefficients
  # Use these for mapping, below (pres/abs needed for cases where space/climate residual models are just for pres/abs)
  if (paste(SpTable.all[sp2],"Summer",sep="") %in% SpTable.s  & paste(SpTable.all[sp2],"Winter",sep="") %in% SpTable.w) {
    Coef.pa.all[sp2,]<-(Coef.pa.s[sp.s,]+Coef.pa.w[sp.w,])/2
    Coef.mean.all[sp2,]<-(Coef.mean.s[sp.s,]+Coef.mean.w[sp.w,])/2
    Coef.lci.all[sp2,]<-(Coef.lci.s[sp.s,]+Coef.lci.w[sp.w,])/2
    Coef.uci.all[sp2,]<-(Coef.uci.s[sp.s,]+Coef.uci.w[sp.w,])/2
  }
  if (paste(SpTable.all[sp2],"Summer",sep="") %in% SpTable.s  & (paste(SpTable.all[sp2],"Winter",sep="") %in% SpTable.w == FALSE)) {
    Coef.pa.all[sp2,]<-Coef.pa.s[sp.s,]
    Coef.mean.all[sp2,]<-Coef.mean.s[sp.s,]
    Coef.lci.all[sp2,]<-Coef.lci.s[sp.s,]
    Coef.uci.all[sp2,]<-Coef.uci.s[sp.s,]
  }
  if (paste(SpTable.all[sp2],"Winter",sep="") %in% SpTable.w  & (paste(SpTable.all[sp2],"Summer",sep="") %in% SpTable.s == FALSE)) {
    Coef.pa.all[sp2,]<-Coef.pa.w[sp.w,]
    Coef.mean.all[sp2,]<-Coef.mean.w[sp.w,]
    Coef.lci.all[sp2,]<-Coef.lci.w[sp.w,]
    Coef.uci.all[sp2,]<-Coef.uci.w[sp.w,]
  }
  # Do cutblock convergence for class 2 and unsampled class 3 and 4
  Coef.mean.all[sp2,"CCSpruce2"]<-Coef.mean.all[sp2,"CCSpruce2"]*(1-0.5)+Coef.mean.all[sp2,"Spruce2"]*0.5
  Coef.mean.all[sp2,"CCSpruce3"]<-Coef.mean.all[sp2,"CCSpruce3"]*(1-0.849)+Coef.mean.all[sp2,"Spruce3"]*0.849
  Coef.mean.all[sp2,"CCSpruce4"]<-Coef.mean.all[sp2,"CCSpruce4"]*(1-0.96)+Coef.mean.all[sp2,"Spruce4"]*0.96
  Coef.mean.all[sp2,"CCPine2"]<-Coef.mean.all[sp2,"CCPine2"]*(1-0.5)+Coef.mean.all[sp2,"Pine2"]*0.5
  Coef.mean.all[sp2,"CCPine3"]<-Coef.mean.all[sp2,"CCPine3"]*(1-0.849)+Coef.mean.all[sp2,"Pine3"]*0.849
  Coef.mean.all[sp2,"CCPine4"]<-Coef.mean.all[sp2,"CCPine4"]*(1-0.96)+Coef.mean.all[sp2,"Pine4"]*0.96
  Coef.mean.all[sp2,"CCMixedwood2"]<-Coef.mean.all[sp2,"CCMixedwood2"]*(1-0.705)+Coef.mean.all[sp2,"Mixedwood2"]*0.705
  Coef.mean.all[sp2,"CCMixedwood3"]<-Coef.mean.all[sp2,"CCMixedwood3"]*(1-0.912)+Coef.mean.all[sp2,"Mixedwood3"]*0.912
  Coef.mean.all[sp2,"CCMixedwood4"]<-Coef.mean.all[sp2,"CCMixedwood4"]*(1-0.97)+Coef.mean.all[sp2,"Mixedwood4"]*0.97
  Coef.mean.all[sp2,"CCDecid2"]<-Coef.mean.all[sp2,"CCDecid2"]*(1-0.705)+Coef.mean.all[sp2,"Decid2"]*0.705
  Coef.mean.all[sp2,"CCDecid3"]<-Coef.mean.all[sp2,"CCDecid3"]*(1-0.912)+Coef.mean.all[sp2,"Decid3"]*0.912
  Coef.mean.all[sp2,"CCDecid4"]<-Coef.mean.all[sp2,"CCDecid4"]*(1-0.97)+Coef.mean.all[sp2,"Decid4"]*0.97
  # And CI's
  Coef.lci.all[sp2,"CCSpruce2"]<-Coef.lci.all[sp2,"CCSpruce2"]*(1-0.5)+Coef.lci.all[sp2,"Spruce2"]*0.5
  Coef.lci.all[sp2,"CCSpruce3"]<-Coef.lci.all[sp2,"CCSpruce3"]*(1-0.849)+Coef.lci.all[sp2,"Spruce3"]*0.849
  Coef.lci.all[sp2,"CCSpruce4"]<-Coef.lci.all[sp2,"CCSpruce4"]*(1-0.96)+Coef.lci.all[sp2,"Spruce4"]*0.96
  Coef.lci.all[sp2,"CCPine2"]<-Coef.lci.all[sp2,"CCPine2"]*(1-0.5)+Coef.lci.all[sp2,"Pine2"]*0.5
  Coef.lci.all[sp2,"CCPine3"]<-Coef.lci.all[sp2,"CCPine3"]*(1-0.849)+Coef.lci.all[sp2,"Pine3"]*0.849
  Coef.lci.all[sp2,"CCPine4"]<-Coef.lci.all[sp2,"CCPine4"]*(1-0.96)+Coef.lci.all[sp2,"Pine4"]*0.96
  Coef.lci.all[sp2,"CCMixedwood2"]<-Coef.lci.all[sp2,"CCMixedwood2"]*(1-0.705)+Coef.lci.all[sp2,"Mixedwood2"]*0.705
  Coef.lci.all[sp2,"CCMixedwood3"]<-Coef.lci.all[sp2,"CCMixedwood3"]*(1-0.912)+Coef.lci.all[sp2,"Mixedwood3"]*0.912
  Coef.lci.all[sp2,"CCMixedwood4"]<-Coef.lci.all[sp2,"CCMixedwood4"]*(1-0.97)+Coef.lci.all[sp2,"Mixedwood4"]*0.97
  Coef.lci.all[sp2,"CCDecid2"]<-Coef.lci.all[sp2,"CCDecid2"]*(1-0.705)+Coef.lci.all[sp2,"Decid2"]*0.705
  Coef.lci.all[sp2,"CCDecid3"]<-Coef.lci.all[sp2,"CCDecid3"]*(1-0.912)+Coef.lci.all[sp2,"Decid3"]*0.912
  Coef.lci.all[sp2,"CCDecid4"]<-Coef.lci.all[sp2,"CCDecid4"]*(1-0.97)+Coef.lci.all[sp2,"Decid4"]*0.97
  Coef.uci.all[sp2,"CCSpruce2"]<-Coef.uci.all[sp2,"CCSpruce2"]*(1-0.5)+Coef.uci.all[sp2,"Spruce2"]*0.5
  Coef.uci.all[sp2,"CCSpruce3"]<-Coef.uci.all[sp2,"CCSpruce3"]*(1-0.849)+Coef.uci.all[sp2,"Spruce3"]*0.849
  Coef.uci.all[sp2,"CCSpruce4"]<-Coef.uci.all[sp2,"CCSpruce4"]*(1-0.96)+Coef.uci.all[sp2,"Spruce4"]*0.96
  Coef.uci.all[sp2,"CCPine2"]<-Coef.uci.all[sp2,"CCPine2"]*(1-0.5)+Coef.uci.all[sp2,"Pine2"]*0.5
  Coef.uci.all[sp2,"CCPine3"]<-Coef.uci.all[sp2,"CCPine3"]*(1-0.849)+Coef.uci.all[sp2,"Pine3"]*0.849
  Coef.uci.all[sp2,"CCPine4"]<-Coef.uci.all[sp2,"CCPine4"]*(1-0.96)+Coef.uci.all[sp2,"Pine4"]*0.96
  Coef.uci.all[sp2,"CCMixedwood2"]<-Coef.uci.all[sp2,"CCMixedwood2"]*(1-0.705)+Coef.uci.all[sp2,"Mixedwood2"]*0.705
  Coef.uci.all[sp2,"CCMixedwood3"]<-Coef.uci.all[sp2,"CCMixedwood3"]*(1-0.912)+Coef.uci.all[sp2,"Mixedwood3"]*0.912
  Coef.uci.all[sp2,"CCMixedwood4"]<-Coef.uci.all[sp2,"CCMixedwood4"]*(1-0.97)+Coef.uci.all[sp2,"Mixedwood4"]*0.97
  Coef.uci.all[sp2,"CCDecid2"]<-Coef.uci.all[sp2,"CCDecid2"]*(1-0.705)+Coef.uci.all[sp2,"Decid2"]*0.705
  Coef.uci.all[sp2,"CCDecid3"]<-Coef.uci.all[sp2,"CCDecid3"]*(1-0.912)+Coef.uci.all[sp2,"Decid3"]*0.912
  Coef.uci.all[sp2,"CCDecid4"]<-Coef.uci.all[sp2,"CCDecid4"]*(1-0.97)+Coef.uci.all[sp2,"Decid4"]*0.97

  # 5.2. Figure for overall coefficients
  # Do the next three lines here to use count.all in figure
  d1$count.all<-ifelse(!is.na(d1$count.summer) & !is.na(d1$count.winter),(d1$count.summer+d1$count.winter)/2, ifelse(is.na(d1$count.summer),0,d1$count.summer)+ifelse(is.na(d1$count.winter),0,d1$count.winter))  # Use unweighted average if both available, because predictions are for equal effort (100 days) per season.  Otherwise, use just whichever season is available
  d1$p.pa.all<-ifelse(!is.na(d1$count.summer) & !is.na(d1$count.winter),(d1$p.pa.summer+d1$p.pa.winter)/2, ifelse(is.na(d1$count.summer),0,d1$p.pa.summer)+ifelse(is.na(d1$count.winter),0,d1$p.pa.winter))  # Use (average) prediction(s) for whichever season(s) have/has available count
  d1$p.all<-ifelse(!is.na(d1$count.summer) & !is.na(d1$count.winter),(d1$p.summer+d1$p.winter)/2, ifelse(is.na(d1$count.summer),0,d1$p.summer)+ifelse(is.na(d1$count.winter),0,d1$p.winter))  # Use (average) prediction(s) for whichever season(s) have/has available count
  x<-col1<-y<-y.lci<-y.uci<-w<-class<-space1<-NULL
  for (i in 1:length(vnames.agp)) {  # vnames.agp doesn't have age classes for natural stands
    if (vnames.agp[i] %in% c("Spruce","Pine","Decid","Mixedwood","TreedBog") == FALSE) {  # Variable without age info
      j<-match(vnames.agp[i],fp$Class)
      x<-c(x,fp$x[j])
      veg.to.use<-as.character(pm$VegType[pm[,vnames.agp[i]]==1])[1]  # To use Coef.mean, find the first fine veg+HF types that is included in the (potentially broader) variable being plotted
      if (veg.to.use %in% c("Spruce","Pine","Decid","Mixedwood","TreedBog")) veg.to.use<-paste(veg.to.use,"R",sep="")  # If no age relationship, use the first value for that fine stand type (don't add the R for types that don't have age classes at all)
      y<-c(y,Coef.mean.all[sp2,veg.to.use])
      y.lci<-c(y.lci,Coef.lci.all[sp2,veg.to.use])
      y.uci<-c(y.uci,Coef.uci.all[sp2,veg.to.use])
      col1<-c(col1,fp$col1[j])
      w<-c(w,fp$width[j])
      class<-c(class,vnames.agp[i])
      space1<-c(space1,fp$spaceafter[j])
    } else {  # Variable with age info
      j<-match(paste(vnames.agp[i],"R",sep=""),fp$Class):match(paste(vnames.agp[i],"8",sep=""),fp$Class)
      x<-c(x,fp$x[j])
      veg.to.use<-as.character(pm$VegType[pm[,vnames.agp[i]]==1])[1]  # To use Coef.mean, find the first fine veg+HF types that is included in the (potentially broader) variable being plotted
      veg.to.use<-paste(veg.to.use,c("R","1","2","3","4","5","6","7","8"),sep="")  # If age relationship, use each value for that fine stand type
      y<-c(y,Coef.mean.all[sp2,veg.to.use])
      y.lci<-c(y.lci,Coef.lci.all[sp2,veg.to.use])
      y.uci<-c(y.uci,Coef.uci.all[sp2,veg.to.use])
      col1<-c(col1,fp$col1[j])
      w<-c(w,fp$width[j])
      class<-c(class,as.character(fp$Class[j]))
      space1<-c(space1,fp$spaceafter[j])
    }
  }
  ord<-order(x)  # Sort all by x
  y<-y[ord]
  y.lci<-y.lci[ord]
  y.uci<-y.uci[ord]
  col1<-col1[ord]
  w<-w[ord]
  class<-class[ord]
  space1<-space1[ord]
  x<-x[ord]
  # Rectify spaces between x's
  for (i in 1:(length(x)-1)) {
    for (j in (i+1):length(x)) x[j]<-x[j]+space1[i]-(x[i+1]-x[i])  # Alter all subsequent positions accordingly
  }
  # Make bar plot
  ymax<-min(max(y.uci,na.rm=TRUE),2*max(y,na.rm=TRUE))  # This keeps the figures readable when there are extreme UCI's
  space<-c(1,x[-1]-x[-length(x)])-0.99  # The spacing between bars
  density<-ifelse(substr(class,1,2)=="CC",50,NA)
  fname<-paste(fname.fig,"Veg+HF figure best model ",SpTable.all[sp2],".png",sep="")
  png(file=fname,width=ifelse(length(y)>5,1500,1000),height=700)
  par(mai=c(1.9,1,0.2,0.3))
  x1<-barplot(y,space=space,width=w,border="white",col="grey30",ylim=c(0,ymax),yaxt="n",ylab="Relative abundance",col.lab="grey50",cex.lab=1.2,axisnames=F)[,1]  # To get strips on CC bars
  abline(h=pretty(c(0,ymax)),col="grey80")
  x1<-barplot(y,space=space,width=w,border="white",col="grey30",ylim=c(0,ymax),yaxt="n",ylab="Relative abundance",col.lab="grey50",cex.lab=1.2,axisnames=F,add=TRUE)[,1]  # To get strips on CC bars, and to put bars in front of horizontal axis lines
  x1<-barplot(y,space=space,width=w,border="white",density=density,col=col1,ylim=c(0,ymax),yaxt="n",ylab="Relative abundance",col.lab="grey50",cex.lab=1.2,axisnames=F,add=TRUE)[,1]
  axis(side=2,tck=0.02,cex.axis=0.9,col.axis="grey50",col.ticks="grey50",las=2,at=pretty(c(0,ymax)))
  box(bty="l",col="grey50")
  for (i in 1:length(x1)) {
    lines(rep(x1[i],2),c(y[i],y.uci[i]),col=col1[i])
    lines(rep(x1[i],2),c(y[i],y.lci[i]),col="grey90")
  }
  for (i in 1:length(class)) {  # Label x axis
    if (substr(class[i],1,2)=="CC") {
      if(substr(class[i],nchar(class[i]),nchar(class[i]))=="2") mtext(side=1,line=0,at=x1[i],"cut",col=col1[i],cex=0.8)
      if (class[i]=="CCAll") mtext(side=1,at=x1[i],line=1,"All",col=col1[i])
      if (class[i]=="CCConif") mtext(side=1,at=x1[i],line=1,"Conif",col=col1[i])
    } else {
      if (substr(class[i],nchar(class[i]),nchar(class[i])) %in% c("R","1","2","3","4","5","6","7","8")==TRUE) {
        if (substr(class[i],nchar(class[i]),nchar(class[i]))=="R") mtext(side=1,line=0,at=x1[i]-0.5,"0",col=col1[i],cex=0.8)
        if (substr(class[i],nchar(class[i]),nchar(class[i]))=="2") mtext(side=1,line=0,at=x1[i]-0.5,"20",col=col1[i],cex=0.8)
        if (substr(class[i],nchar(class[i]),nchar(class[i]))=="4") {
          mtext(side=1,line=0,at=x1[i]-0.5,"60",col=col1[i],cex=0.8)
          mtext(side=1,line=1,at=x1[i],substr(class[i],1,nchar(class[i])-1),col=col1[i],cex=1.2)
        }
        if (substr(class[i],nchar(class[i]),nchar(class[i]))=="6") mtext(side=1,line=0,at=x1[i]-0.5,"100",col=col1[i],cex=0.8)
        if (substr(class[i],nchar(class[i]),nchar(class[i]))=="8") mtext(side=1,line=0,at=x1[i]-0.5,"140",col=col1[i],cex=0.8)
      } else {
        if (class[i]=="Crop") class[i]<-"Crop+"  # Because this also includes other alienating
        mtext(side=1,at=x1[i],line=1,adj=ifelse(w[i]>2 | length(y)<9,0.5,1),las=ifelse(w[i]>2 | length(y)<9,1,2),class[i],col=col1[i],cex=1.2)
      }  # End if for aged treed type
    } # End if for CC
  }
  mtext(side=3,at=x1[1],adj=0,paste(sp.names.all[sp2],"- North"),col="grey30",cex=1.2)
  text(max(x1),ymax*0.98,paste("Detected at",sum(sign(d1$count.all)),"of",nrow(d1),"camera locations"),cex=1.1,adj=1,col="grey40") # Add sample size
  graphics.off()

  # 6. Residual variation due to location and climate
  # These models are fit to the residual variation after the average Winter and winter predictions (or whichever one(s) are available).
  # Average Winter and winter counts and predictions, and add explanatory variables from d to d1.  Check first that DeploymentYears in d1 and d are still in same order
  dplot(1:nrow(d),match(d$DeploymentYear,d1$DeploymentYear),cex=0.3)  # Needs to be a straight 1:1 line
  d1<-cbind(d1,d[,c("Lat","TrueLat","Long","PET","AHM","MAT","FFP","MAP","MWMT","MCMT","NSR1")])
  if (sc.option[sp2]==2 | sc.option[sp2]==1) { # Presence/absence model only.  It is run for option 1 here solely to do fit AUC; the total abundance models replace these in the next section for option 1.
    # Climate and spatial variable sets are tried
    dtemp<-d1[substr(d1$DeploymentYear,1,4)=="ABMI" & substr(d1$DeploymentYear,1,6)!="ABMI-W",]  # To calculate lure effect only with ABMI (notABMI-W) sites
    lure1<-mean(sign(dtemp$count.all[dtemp$Lured=="y"]))/mean(sign(dtemp$count.all[dtemp$Lured=="n"]))  # Lure effect on presence/absence
    pCount1<-sign(d1$count.all)/ifelse(d1$Lured=="y",lure1,1)
    pCount1<-pCount1/max(pCount1)  # In case lure effect is <1
    m.sc.pa<-list(NULL)
    wt1<-d1$wt  # Need to do this to use model.matrix below for plotting
    m.sc.pa[[1]]<-try(glm(pCount1~offset(qlogis(0.998*d1$p.pa.all+0.001))-1,data=d1,family="binomial",weights=wt1))
    m.sc.pa[[2]]<-try(update(m.sc.pa[[1]],.~.+Lat+Long+Lat:Long))  # Note that "Lat" here is the truncated latitude, where southern points are treated as being further north.  ("TrueLat" is the actual latitude of the site)
    m.sc.pa[[3]]<-try(update(m.sc.pa[[1]],.~.+Lat+Long+Lat:Long+I(Lat^2)+I(Long^2)))
    m.sc.pa[[4]]<-try(update(m.sc.pa[[1]],.~.+Lat+Long+Lat:Long+I(Lat^2)+I(Long^2)+I(Lat^3)))
    m.sc.pa[[5]]<-try(update(m.sc.pa[[1]],.~.+PET))  # Was EREF, but not available in current climate variable summary, so all EREF changed to PET
    m.sc.pa[[6]]<-try(update(m.sc.pa[[1]],.~.+AHM))
    m.sc.pa[[7]]<-try(update(m.sc.pa[[1]],.~.+MAT))
    m.sc.pa[[8]]<-try(update(m.sc.pa[[1]],.~.+FFP))
    m.sc.pa[[9]]<-try(update(m.sc.pa[[1]],.~.+MAP+FFP))
    m.sc.pa[[10]]<-try(update(m.sc.pa[[1]],.~.+MAP+FFP+MAP:FFP))
    m.sc.pa[[11]]<-try(update(m.sc.pa[[1]],.~.+MAT+MAP+PET+AHM))
    m.sc.pa[[12]]<-try(update(m.sc.pa[[1]],.~.+MAT+MAP+PET+AHM+PET:MAP+MAT:AHM))
    m.sc.pa[[13]]<-try(update(m.sc.pa[[1]],.~.+MAT+MAP))
    m.sc.pa[[14]]<-try(update(m.sc.pa[[1]],.~.+MWMT+MCMT))
    m.sc.pa[[15]]<-try(update(m.sc.pa[[1]],.~.+AHM+PET))
    m.sc.pa[[16]]<-try(update(m.sc.pa[[1]],.~.+MAT+I(MAT*(MAT+10)) ))  #
    m.sc.pa[[17]]<-try(update(m.sc.pa[[1]],.~.+MWMT+MCMT+FFP+MAT))
    m.sc.pa[[18]]<-try(update(m.sc.pa[[1]],.~.+Lat+Long+Lat:Long+I(Lat^2)+I(Long^2)+I(Lat^2*Long^2)))
    m.sc.pa[[19]]<-try(update(m.sc.pa[[1]],.~.+MAT+I(MAT*(MAT+10))+I(Long*MAT) ))
    #		for (i in 1:19) m.sc.pa[[i+19]]<-update(m.sc.pa[[i]],.~.+NSR1)
    # Models with spatial and climate variables together not currently used, because highly correlated
    nModels.sc.pa<-length(m.sc.pa)
    # BIC calculation to select best covariate set Uses BIC for more conservative variable set
    bic.sc.pa<-rep(999999999,nModels.sc.pa)
    for (i in 1:(nModels.sc.pa)) {
      if (!is.null(m.sc.pa[[i]]) & class(m.sc.pa[[i]])[1]!="try-error") {
        bic.sc.pa[i]<-BIC(m.sc.pa[[i]])
      }
    }
    best.model.sc.pa<-which.min(bic.sc.pa)
  }

  # Do AUC of fit here, before the total abundance option is run for option 1
  if (sc.option[sp2]==1 | sc.option[sp2]==2) p<-plogis(predict(m.sc.pa[[best.model.sc.pa]]))
  if (sc.option[sp2]==3) p<-plogis(d1$p.pa.all)  # The original veg+HF only prediction if there is no sc model (this prediction is on the ordinal scale, despite the "t." in the name...)
  auc.fit[sp2]<-auc(roc(sign(d1$count.all),p))  # No correction for lure here.

  if (sc.option[sp2]==1) { # Model full total abundance
    # Climate and spatial variables sets are tried
    dtemp<-d1[substr(d1$DeploymentYear,1,4)=="ABMI" & substr(d1$DeploymentYear,1,6)!="ABMI-W",]  # To calculate lure effect only with ABMI (notABMI-W) sites
    lure1<-mean(dtemp$count.all[dtemp$Lured=="y"])/mean(dtemp$count.all[dtemp$Lured=="n"])  # Lure effect on total abundance
    pCount1<-d1$count.all/ifelse(d1$Lured=="y",lure1,1)
    log.offset<-min(pCount1[pCount1>0])/2
    pCount1<-log(pCount1+log.offset)  # Using log(x+offset)
    log.adj<-mean(log(d1$p.all+log.offset)-pCount1)  # Compensate for geometric mean issue
    pCount1<-pCount1+log.adj
    m.sc.pa<-list(NULL)
    wt1<-d1$wt  # Need to do this to use model.matrix below for plotting
    m.sc.pa[[1]]<-try(glm(pCount1~offset(log(d1$p.all+log.offset))-1,data=d1,family=gaussian,weights=wt1))
    m.sc.pa[[2]]<-try(update(m.sc.pa[[1]],.~.+Lat+Long+Lat:Long))  # Note that "Lat" here is the truncated latitude, where southern points are treated as being further north.  ("TrueLat" is the actual latitude of the site)
    m.sc.pa[[3]]<-try(update(m.sc.pa[[1]],.~.+Lat+Long+Lat:Long+I(Lat^2)+I(Long^2)))
    m.sc.pa[[4]]<-try(update(m.sc.pa[[1]],.~.+Lat+Long+Lat:Long+I(Lat^2)+I(Long^2)+I(Lat^3)))
    m.sc.pa[[5]]<-try(update(m.sc.pa[[1]],.~.+PET))
    m.sc.pa[[6]]<-try(update(m.sc.pa[[1]],.~.+AHM))
    m.sc.pa[[7]]<-try(update(m.sc.pa[[1]],.~.+MAT))
    m.sc.pa[[8]]<-try(update(m.sc.pa[[1]],.~.+FFP))
    m.sc.pa[[9]]<-try(update(m.sc.pa[[1]],.~.+MAP+FFP))
    m.sc.pa[[10]]<-try(update(m.sc.pa[[1]],.~.+MAP+FFP+MAP:FFP))
    m.sc.pa[[11]]<-try(update(m.sc.pa[[1]],.~.+MAT+MAP+PET+AHM))
    m.sc.pa[[12]]<-try(update(m.sc.pa[[1]],.~.+MAT+MAP+PET+AHM+PET:MAP+MAT:AHM))
    m.sc.pa[[13]]<-try(update(m.sc.pa[[1]],.~.+MAT+MAP))
    m.sc.pa[[14]]<-try(update(m.sc.pa[[1]],.~.+MWMT+MCMT))
    m.sc.pa[[15]]<-try(update(m.sc.pa[[1]],.~.+AHM+PET))
    m.sc.pa[[16]]<-try(update(m.sc.pa[[1]],.~.+MAT+I(MAT*(MAT+10)) ))  #
    m.sc.pa[[17]]<-try(update(m.sc.pa[[1]],.~.+MWMT+MCMT+FFP+MAT))
    m.sc.pa[[18]]<-try(update(m.sc.pa[[1]],.~.+Lat+Long+Lat:Long+I(Lat^2)+I(Long^2)+I(Lat^2*Long^2)))
    m.sc.pa[[19]]<-try(update(m.sc.pa[[1]],.~.+MAT+I(MAT*(MAT+10))+I(Long*MAT) ))  # Version including Long did not work as well
    #		for (i in 1:19) m.sc.pa[[i+19]]<-update(m.sc.pa[[i]],.~.+NSR1)
    # Models with spatial and climate variables together not currently used, because highly correlated
    nModels.sc<-length(m.sc.pa)
    # BIC calculation to select best covariate set Uses BIC for more conservative variable set
    bic.sc<-rep(999999999,nModels.sc)
    for (i in 1:(nModels.sc)) {
      if (!is.null(m.sc.pa[[i]]) & class(m.sc.pa[[i]])[1]!="try-error") {
        bic.sc[i]<-BIC(m.sc.pa[[i]])
      }
    }
    best.model.sc.pa<-which.min(bic.sc)
  }
  # No models if option=3

  # Change models if necessary when a non-best model is more accurate than the best model
  if (SpTable.all[sp2]=="Muledeer") best.model.sc.pa<-11

  c1<-coef(m.sc.pa[[best.model.sc.pa]])  # Variable names in best sc model

  # And post-hoc modifications of coefficients when necessary
  if (SpTable.all[sp2]=="Grizzlybear") c1<-c(c1,"NSR1North"= -16, "NSR1Shield"= -16, "NSR1CentralMixedwood" = -2)  # Add to (partially) censor those NSR1's

  # 6.2 Save coefficients for the subset of climate and/or spatial variables
  if (sc.option[sp2]==1 | sc.option[sp2]==2) {
    vnames<-names(c1)
    vnames1<-ifelse(vnames=="(Intercept)","Intercept",vnames)  # This is for the names used in km2.res (where "(", "^", etc. can't be used)
    vnames1<-ifelse(vnames1=="I(Lat^2)","Lat2",vnames1)
    vnames1<-ifelse(vnames1=="I(Lat^3)","Lat3",vnames1)
    vnames1<-ifelse(vnames1=="I(Long^2)","Long2",vnames1)
    vnames1<-ifelse(vnames1=="I(MAT * (MAT + 10))","MAT2",vnames1)
    vnames1<-ifelse(vnames1=="I(MWMT^2)","MWMT2",vnames1)
    vnames1<-ifelse(vnames1=="I(Lat^2 * Long^2)","Lat2Long2",vnames1)
    vnames1<-ifelse(vnames1=="I(Long * MAT)","LongMAT",vnames1)
    vnames1<-gsub(":","",vnames1)
    Res.coef[sp2,match(vnames1,colnames(Res.coef))]<-c1
    if (length(c1)==0) Res.coef[sp2,]<-0  # For case when models were fit, but best one was null
  }
  # If option=3, the res coefs are already 0

  # 6.3 Mapping of residual climate/spatial effect
  # This is just the additional climate and spatial effect, not the basic veg type effects
  # Done in the original way, since this is just for us
  if (sc.option[sp2]==1 | sc.option[sp2]==2) {
    vnames<-names(c1)  # Figure out the names of the included variables in the km2 raster data frame
    vnames<-ifelse(vnames=="(Intercept)","Intercept",vnames)
    vnames<-ifelse(vnames=="I(Lat^2)","Lat2",vnames)
    vnames<-ifelse(vnames=="I(Lat^3)","Lat3",vnames)
    vnames<-ifelse(vnames=="I(Long^2)","Long2",vnames)
    vnames<-ifelse(vnames=="I(MAT * (MAT + 10))","MAT2",vnames)
    vnames<-ifelse(vnames=="I(MWMT^2)","MWMT2",vnames)
    vnames<-ifelse(vnames=="I(Lat^2 * Long^2)","Lat2Long2",vnames)
    vnames<-ifelse(vnames=="I(Long * MAT)","LongMAT",vnames)
    vnames<-gsub(":","",vnames)
    km2.sc1<-km2.sc[,vnames]
    if (sc.option[sp2]==1) km2.p1<-exp(rowSums( t(c1*t(km2.sc1))))   # Predictions from just the climate and spatial part of the residual model for each km2 raster - log-link total abundance
    if (sc.option[sp2]==2) km2.p1<-plogis(rowSums( t(c1*t(km2.sc1))))   # Predictions from just the climate and spatial part of the residual model for each km2 raster - logit-link presence/absence
    km2.p1<-ifelse(km2.p1>quantile(km2.p1,0.99),quantile(km2.p1,0.99),km2.p1)
    km2.p1<-km2.p1/max(km2.p1)
    if (max(km2.p1)==min(km2.p1)) {
      km2.p1<-rep(0.5,length(km2.p1))  # make all white map if null model is best
    } else {
      km2.p1<-(km2.p1-min(km2.p1))/(max(km2.p1)-min(km2.p1))
    }
    r<-ifelse(km2.p1<0.5,1,(1-km2.p1)*2)[!is.na(km2.p1)]  # RGB for map (white at no change, more green for more positive residual effect, more red for more negative
    g<-ifelse(km2.p1>0.5,1,km2.p1*2)[!is.na(km2.p1)]  # RGB for map
    b<-(1-abs(km2.p1-0.5)*2)[!is.na(km2.p1)]  # RGB for map
    fname<-paste(fname.map,"Climate and spatial/",SpTable.all[sp2],".jpg",sep="")
    jpeg(filename=fname,width=5,height=8.7,units="in",res=300)
    dplot(km2.sc$Long[!is.na(km2.p1)],km2.sc$TrueLat[!is.na(km2.p1)],pch=15,cex=0.4,col=rgb(r,g,b),xlab="",ylab="")
    points(km2.water$Long,km2.water$Lat,pch=15,cex=0.3,col=rgb(0.5,0.4,0.9))
    points(d1$Long,d1$TrueLat,cex=1.5*sqrt(d1$count.all)+0.2,col="yellow",lwd=2)  # Add data points, size proportional to count
    title(paste("Climate+spatial",sp.names.all[sp2]))
    graphics.off()
  }

  # 6.4 Full map for North
  # This includes veg types, surrounding HF and additional climate and spatial effects
  km2.pveg<-colSums(Coef.mean.all[sp2,]*t(km2[,colnames(Coef.mean.all)]))   # Prediction based on veg types only.  Note: water, barren and mines not included, so they are treated as 0.
  km2.pres<-colSums(Res.coef[sp2,]*t(km2.sc[,colnames(Res.coef)]))  # Prediction of residual effect (note: Uses truncated latitude "Lat", but point is plotted at true latitude)
  km2.pveg.ref<-colSums(Coef.mean.all[sp2,vnames.b]*t(km2.b[,vnames.b]))   # Prediction based on non-HF veg types only for reference
  km2.pres.ref<-colSums(Res.coef[sp2,]*t(km2.sc[,colnames(Res.coef)]))  # Prediction of residual effect (note: Uses truncated latitude "Lat", but point is plotted at true latitude)
  if (sc.option[sp2]==1 | sc.option[sp2]==3) {  # Use this also for no model - all 0 coefficients become 1 multipliers on exp scale
    km2.p<-km2.pveg*exp(km2.pres) # Using simple multiplication of residual effect, to avoid offset problems
    km2.p.ref<-km2.pveg.ref*exp(km2.pres.ref)  # Using simple multiplication of residual effect, to avoid offset problems - here, for log-linked residual total abundance predictions
  }
  if (sc.option[sp2]==2) {  # Logit scale, and need to extract veg presence/absence and AGP components
    km2.pveg.pa<-colSums(Coef.pa.all[sp2,]*t(km2[,colnames(Coef.pa.all)]))   # Presence/absence prediction based on veg types only.  Note: water, barren and mines not included, so they are treated as 0.
    km2.pveg.agp<-ifelse(km2.pveg.pa<0.0001,0,km2.pveg/km2.pveg.pa)  # This is the abundance-given-presence prediction for each km2 raster, to be multiplied by the sc-adjusted presence/absence
    km2.p<-plogis(qlogis(0.998*km2.pveg.pa+0.001)+km2.pres) * km2.pveg.agp  # Using simple multiplication of residual effect, to avoid offset problems - here, for logit-linked residual presence/absence predictions.  And multiply by agp to get total abundance.
    km2.p<-ifelse(is.na(km2.p),0,km2.p)  # This is if the prediction is 0
    # And repeat for reference
    km2.pveg.pa.ref<-colSums(Coef.pa.all[sp2,vnames.b]*t(km2.b[,vnames.b]))   # Presence/absence prediction based on veg types only.  Note: water, barren and mines not included, so they are treated as 0.
    km2.pveg.agp.ref<-ifelse(km2.pveg.pa.ref<0.0001,0,km2.pveg.ref/km2.pveg.pa.ref)  # This is the abundance-given-presence prediction for each km2 raster, to be multiplied by the sc-adjusted presence/absence
    km2.p.ref<-plogis(qlogis(0.998*km2.pveg.pa.ref+0.001)+km2.pres.ref) * km2.pveg.agp.ref  # Using simple multiplication of residual effect, to avoid offset problems - here, for logit-linked residual presence/absence predictions. And multiply by agp to get total abundance.
    km2.p.ref<-ifelse(is.na(km2.p.ref),0,km2.p.ref)  # This is if the prediction is 0
  }
  x.ref1<-km2.p.ref  # Colour gradient direct with predicted abundance
  x.curr1<-km2.p
  #	x.curr1<-ifelse(x.curr1<min(pCount),min(pCount),x.curr1)  # Some extreme low values otherwise
  #	x.ref1<-ifelse(x.ref1<min(pCount),min(pCount),x.ref1)
  x.ref.trunc<-ifelse(x.ref1>quantile(c(x.ref1,x.curr1),0.99,na.rm=T),quantile(c(x.ref1,x.curr1),0.99,na.rm=T),x.ref1)  # Clip to 99 percentile, to prevent colour scaling problems with a few high values
  x.curr.trunc<-ifelse(x.curr1>quantile(c(x.ref1,x.curr1),0.99,na.rm=T),quantile(c(x.ref1,x.curr1),0.99,na.rm=T),x.curr1)  # Clip to 99 percentile, to prevent colour scaling problems with a few high values
  x.curr<-x.curr.trunc/max(c(x.ref.trunc,x.curr.trunc),na.rm=T)
  x.ref<-x.ref.trunc/max(c(x.ref.trunc,x.curr.trunc),na.rm=T)
  # 6.4.1 Current
  c3<-c2(1000)[1+999*x.curr]
  fname<-paste(fname.map,"Current/",SpTable.all[sp2],".jpg",sep="")
  jpeg(file=fname,width=600,height=1000)
  plot(km2.sc$proj.x,km2.sc$proj.y,pch=15,cex=0.2,col=c3,xaxt="n",yaxt="n",xlab="",ylab="",bty="n",xlim=range(c(km2.sc$proj.x,m.water$x)),ylim=range(c(km2.sc$proj.y,m.water$y)))
  points(m.water$x,m.water$y,pch=15,cex=0.2,col=rgb(0.4,0.3,0.8))
  #points(km2.sc$proj.x[km2$NatRegion=="Rocky Mountain"],km2.sc$proj.y[km2$NatRegion=="Rocky Mountain"],pch=15,cex=0.2,col="lightcyan4")
  #text(-0.025,-0.745,"Insufficient \n   data",col="white",cex=0.9)
  mtext(side=3,at=m.title$x,paste(sp.names.all[sp2],"Current"),adj=0.5,cex=1.4,col="grey40")
  mtext(side=3,at=m.title$x,line=-1,paste("Detected at",sum(sign(d1$count.all)),"of",nrow(d),"camera locations"),adj=0.5,cex=1.2,col="grey40")
  points(m.city.x,m.city.y,pch=18,col="grey10")
  text(m.city.x,m.city.y,city,cex=0.8,adj=-0.1,col="grey10")
  graphics.off()
  # 6.4.2 Reference
  c3<-c2(1000)[1+999*x.ref]
  fname<-paste(fname.map,"Reference/",SpTable.all[sp2],".jpg",sep="")
  jpeg(file=fname,width=600,height=1000)
  plot(km2.sc$proj.x,km2.sc$proj.y,pch=15,cex=0.2,col=c3,xaxt="n",yaxt="n",xlab="",ylab="",bty="n",xlim=range(c(km2.sc$proj.x,m.water$x)),ylim=range(c(km2.sc$proj.y,m.water$y)))
  points(m.water$x,m.water$y,pch=15,cex=0.2,col=rgb(0.4,0.3,0.8))
  #points(km2$proj.x[km2$NatRegion=="Rocky Mountain"],km2$proj.y[km2$NatRegion=="Rocky Mountain"],pch=15,cex=0.2,col="lightcyan4")
  #text(-0.025,-0.745,"Insufficient \n   data",col="white",cex=0.9)
  mtext(side=3,at=m.title$x,paste(sp.names.all[sp2],"Reference"),adj=0.5,cex=1.4,col="grey40")
  mtext(side=3,at=m.title$x,line=-1,paste("Detected at",sum(sign(d1$count.all)),"of",nrow(d),"camera locations"),adj=0.5,cex=1.2,col="grey40")
  points(m.city.x,m.city.y,pch=18,col="grey10")
  text(m.city.x,m.city.y,city,cex=0.8,adj=-0.1,col="grey10")
  graphics.off()
  # 6.4.3 Difference
  trunc<-quantile(c(km2.p,km2.p.ref),0.99,na.rm=T)
  diff<-(ifelse(km2.p>trunc,trunc,km2.p)-ifelse(km2.p.ref>trunc,trunc,km2.p.ref))/trunc
  diff<-sign(diff)*abs(diff)^0.5
  d3<-d2(1000)[500+499*diff]
  fname<-paste(fname.map,"Difference/",SpTable.all[sp2],".jpg",sep="")
  jpeg(file=fname,width=600,height=1000)
  plot(km2.sc$proj.x,km2.sc$proj.y,pch=15,cex=0.2,col=d3,xaxt="n",yaxt="n",xlab="",ylab="",bty="n",xlim=range(c(km2.sc$proj.x,m.water$x)),ylim=range(c(km2.sc$proj.y,m.water$y)))
  points(m.water$x,m.water$y,pch=15,cex=0.2,col=rgb(0.4,0.3,0.8))
  #points(km2$proj.x[km2$NatRegion=="Rocky Mountain"],km2$proj.y[km2$NatRegion=="Rocky Mountain"],pch=15,cex=0.2,col="lightcyan4")
  #text(-0.025,-0.745,"Insufficient \n   data",col="white",cex=0.9)
  mtext(side=3,at=m.title$x,paste(sp.names.all[sp2],"Difference"),adj=0.5,cex=1.4,col="grey40")
  mtext(side=3,at=m.title$x,line=-1,paste("Detected at",sum(sign(d1$count.all)),"of",nrow(d),"camera locations"),adj=0.5,cex=1.2,col="grey40")
  points(m.city.x,m.city.y,pch=18,col="grey10")
  text(m.city.x,m.city.y,city,cex=0.8,adj=-0.1,col="grey10")
  graphics.off()

  # 7. Output by species
  # 7.1 km2 raster reference and current abundance
  q<-data.frame(LinkID=km2$LinkID,Ref=km2.p.ref,Curr=km2.p)
  fname<-paste(fname.km2summaries," ",SpTable.all[sp2],".csv",sep="")
  write.table(q,file=fname,sep=",",row.names=FALSE)

  # 7.2 Save models, AIC wts, etc. for each species
  fname<-paste(fname.Robjects," ",SpTable.all[sp2],".Rdata",sep="")
  if (sc.option[sp2]==3) {  # Just so it has something to save in these cases
    m.sc.pa<-list(NULL)
    bic.sc.pa<-NA
  }
  if (paste(SpTable.all[sp2],"Summer",sep="") %in% SpTable.s  & paste(SpTable.all[sp2],"Winter",sep="") %in% SpTable.w) {
    save(file=fname,m.pa.s,m.pa.w,m.agp.s,m.agp.w,m.age.s,m.age.w,m.null.age.s,m.null.age.w,m.sc.pa,aic.wt.pa.save.s,aic.wt.pa.save.w,aic.wt.agp.save.s,aic.wt.agp.save.w,aic.age.s,aic.age.w,aic.wt.age.models.save.s,aic.wt.age.models.save.w,aic.wt.age.save.s,aic.wt.age.save.w,bic.sc.pa)
  }
  if (paste(SpTable.all[sp2],"Summer",sep="") %in% SpTable.s  & (paste(SpTable.all[sp2],"Winter",sep="") %in% SpTable.w == FALSE)) {
    save(file=fname,m.pa.s,m.agp.s,m.age.s,m.null.age.s,m.sc.pa,aic.wt.pa.save.s,aic.wt.agp.save.s,aic.age.s,aic.wt.age.models.save.s,aic.wt.age.save.s,bic.sc.pa)
  }
  if (paste(SpTable.all[sp2],"Winter",sep="") %in% SpTable.w  & (paste(SpTable.all[sp2],"Summer",sep="") %in% SpTable.s == FALSE)) {
    save(file=fname,m.pa.w,m.agp.w,m.age.w,m.null.age.w,m.sc.pa,aic.wt.pa.save.w,aic.wt.agp.save.w,aic.age.w,aic.wt.age.models.save.w,aic.wt.age.save.w,bic.sc.pa)
  }

}  # Next species

# 8. Export .csv files for website or other uses
# Veg and HF coefficients
lu.names<-read.csv("C:/Dave/ABMI/Cameras/Coefficients/2018/Analysis north/Lookup for coefficient table names North Mar 2019.csv")  # To translate names used for coefficients here to official names
Coef.mean.all<-cbind(Coef.mean.all,rep(0,nrow(Coef.mean.all)),rep(0,nrow(Coef.mean.all)))  # Add columns for bare and water
Coef.lci.all<-cbind(Coef.lci.all,rep(0,nrow(Coef.lci.all)),rep(0,nrow(Coef.lci.all)))  # Add columns for bare and water
Coef.uci.all<-cbind(Coef.uci.all,rep(0,nrow(Coef.uci.all)),rep(0,nrow(Coef.uci.all)))  # Add columns for bare and water
colnames(Coef.mean.all)[(ncol(Coef.mean.all)-1):ncol(Coef.mean.all)]<-c("Bare","Water")  # Currently not modeled so assumed to be 0
i<-match(lu.names$VegHF,colnames(Coef.mean.all))  # Check for NA's
Coef.official<-Coef.mean.all[,i]
i<-match(lu.names$VegHF,colnames(Coef.lci.all))  # Check for NA's
Coef.official.lci<-Coef.lci.all[,i]
i<-match(lu.names$VegHF,colnames(Coef.uci.all))  # Check for NA's
Coef.official.uci<-Coef.uci.all[,i]
colnames(Coef.official)<-colnames(Coef.official.lci)<-colnames(Coef.official.uci)<-lu.names$CoefName
rownames(Coef.official)<-rownames(Coef.official.lci)<-rownames(Coef.official.uci)<-SpTable.all
# and climate/space coefficients
Res.coef.official<-Res.coef  # Seems to be in right format already
# Save
fname<-paste(fname.sumout,"OFFICIAL coefficients.Rdata")
save(file=fname,Coef.official,Coef.official.lci,Coef.official.uci,Res.coef.official)  # Save to compile later with South results

# AUC fit
q<-data.frame(Sp=SpTable.all,AUC.fit=auc.fit)
write.table(q,file="C:/Dave/ABMI/Cameras/Coefficients/2018/Analysis north/AUC of fit for camera mammals North Mar 2019.csv",sep=",",row.names=F)

# Lure effects
fname<-"C:/Dave/ABMI/Cameras/Coefficients/2018/Analysis north/Lure effects North.csv"
q<-data.frame(Season="Summer",Measure="PresAbs",t(lure.pa.s))
q<-rbind(q,data.frame(Season="Summer",Measure="AGP",t(lure.agp.s)))
q<-rbind(q,data.frame(Season="Winter",Measure="PresAbs",t(lure.pa.w)))
q<-rbind(q,data.frame(Season="Winter",Measure="AGP",t(lure.agp.w)))
names(q)<-c("Season","Measure",SpTable.all)
write.table(q,file=fname,sep=",",row.names=FALSE)

# Veg+HF Model weights
fname<-"C:/Dave/ABMI/Cameras/Coefficients/2018/Analysis north/Model weights Veg+HF North.csv"
q<-data.frame(Season="Summer",Measure="PresAbs",Model=1:ncol(aic.wt.pa.save.s),t(aic.wt.pa.save.s))
q<-rbind(q,data.frame(Season="Summer",Measure="AGP",Model=1:ncol(aic.wt.agp.save.s),t(aic.wt.agp.save.s)))
q<-rbind(q,data.frame(Season="Winter",Measure="PresAbs",Model=1:ncol(aic.wt.pa.save.w),t(aic.wt.pa.save.w)))
q<-rbind(q,data.frame(Season="Winter",Measure="AGP",Model=1:ncol(aic.wt.agp.save.w),t(aic.wt.agp.save.w)))
names(q)<-c("Season","Measure","Model",SpTable.all)
write.table(q,file=fname,sep=",",row.names=FALSE)

# Age Model weights
fname<-"C:/Dave/ABMI/Cameras/Coefficients/2018/Analysis north/Model weights Age stand groupings.csv"
q<-data.frame(Season="Summer",Measure="PresAbs",Grouping=c("Separate","Intermediate","Combined"),t(aic.wt.age.save.s))
q<-rbind(q,data.frame(Season="Winter",Measure="PresAbs",Grouping=c("Separate","Intermediate","Combined"),t(aic.wt.age.save.w)))
names(q)<-c("Season","Measure","Grouping",SpTable.all)
write.table(q,file=fname,sep=",",row.names=FALSE)
fname<-"C:/Dave/ABMI/Cameras/Coefficients/2018/Analysis north/Model weights Age spline versus null.csv"
q<-data.frame(Season="Summer",Measure="PresAbs",VegType=c("Spruce","Pine","Decid","Mixedwood","TreedBog","UpCon","DecidMixed","All"),t(aic.wt.age.models.save.s))
q<-rbind(q,data.frame(Season="Winter",Measure="PresAbs",VegType=c("Spruce","Pine","Decid","Mixedwood","TreedBog","UpCon","DecidMixed","All"),t(aic.wt.age.models.save.w)))
write.table(q,file=fname,sep=",",row.names=FALSE)

# 9. Use/availability figures for all species with >3 detections (separate section, because larger species list)
# Uses Peter's terminology and lines from his script
library(RColorBrewer)
#d1<-d[d$NR=="Boreal" | d$NR=="Foothills" | d$NR=="Canadian Shield" ,]   # Excluding parkland for this
d1<-d  # NOT excluding parkland for this
x<-data.frame(Deciduous=d1$Decid,Mixedwood=d1$Mixedwood,WhiteSpruce=d1$Spruce,Pine=d1$Pine,BlackSpruce=d1$TreedBog,TreedFen=d1$TreedFen,Open=d1$GrassShrub,Wetland=d1$OpenWet,
              HFor=d1$HFor,Crop=d1$Crop,TameP=d1$TamePasture,RoughP=d1$RoughPasture,RuralResInd=d1$RuralResInd,Wells=d1$WellSite,UrbInd=d1$RuralResInd,HardLin=d1$HardLin,SoftLin=d1$SoftLin)  # No UrbInd sampled
pAvail<-colMeans(x*d1$TotalDays)  # Correction here is only for sampling length.  Detection distance not included, for simplicity
pAvail<-pAvail/sum(pAvail)  # Proportional availability of types
col<-brewer.pal(8, "Accent")[c(1,1,1,1, 2,2,3,5, 7,rep(8,8))]
op<-par(mar=c(6,4,2,2)+0.1, las=2)
SpTable.ua<-sort(unique(gsub("Winter","",SpTable.w.ua))) # Original SpTable.w.ua and SpTable.s.ua mixed up.  Just check that this is all species
SpTable.ua<-sort(unique(gsub("Summer","",SpTable.ua)))
WRSI<-rWRSI<-array(NA,c(length(SpTable.ua),length(col)))  # Store results for each species
colnames(WRSI)<-colnames(rWRSI)<-names(x)
for (sp in 1:length(SpTable.ua)) {
  y.w<-d[,paste(SpTable.ua[sp],"Winter",sep="")]
  y.s<-d[,paste(SpTable.ua[sp],"Summer",sep="")]
  y.ave<-(ifelse(is.na(y.w),0,y.w*d$WinterDays)+ifelse(is.na(y.s),0,y.s*d$SummerDays)) / (ifelse(is.na(y.s),0,d$SummerDays) + ifelse(is.na(y.w),0,d$WinterDays))
  fname<-paste(fname.fig,"Use availability/",SpTable.ua[sp],".png",sep="")
  png(file=fname,width=480,height=480)
  par(mar=c(6,4,2,2)+0.1, las=2)
  pUse<-colMeans(x * sign(y.ave))
  pUse<-pUse/sum(pUse)
  WRSI[sp,]<-pUse/pAvail  # Simple use/availability
  rWRSI[sp,]=(exp(2 * log(WRSI[sp,])) - 1)/(1 + exp(2 * log(WRSI[sp,])))  # Transform to -1 to 1 (from Peter)
  x1<-barplot(rWRSI[sp,], horiz=FALSE, ylab="Affinity",space=NULL, col=col, border=col, ylim=c(-1,1), axes=FALSE)
  axis(side=2)
  abline(h=0, col="red4", lwd=2)
  mtext(side=3,at=x1[1],adj=0,taxa.file$Species[match(SpTable.ua[sp],taxa.file$Label)],cex=1.2,col="grey40",las=1)
  # Add sample size
  text(max(x1),0.97,paste("Detected at",sum(sign(y.ave)),"of",nrow(d),"camera locations"),cex=1.2,adj=1,col="grey40")
  graphics.off()
}
par(op)
# Save table for all species
#ua<-data.frame(SpLabel=SpTable.ua,Species=taxa.file$Species[match(SpTable.ua,taxa.file$Label)],WRSI,rWRSI)  # Original version
#names(ua)<-c("SpLabel","Species",paste(names(x),"WRSI"),paste(names(x),"rWRSI"))  # Original version
UseavailNorth<-rWRSI  # Official version
rownames(UseavailNorth)<-SpTable.ua  # Official version
fname<-paste(fname.useavail,"UseavailNorth.Rdata",sep="")
save(UseavailNorth,file=fname)

# 10. Information for header file saying what information is available for what species - to be compiled with south
q<-data.frame(SpeciesID=SpTable.ua,ScientificName=NA,TSNID=NA,CommonName=taxa.file$Species[match(SpTable.ua,taxa.file$Label)],
              ModelNorth=!is.na(match(SpTable.ua,SpTable.all)),ModelSouth=NA,Nonnative=FALSE,
              LinkHabitatNorth=ifelse(is.na(match(SpTable.ua,SpTable.all)),NA,"Logit/Log"),LinkHabitatSouth=NA,
              LinkSpclimNorth=ifelse(is.na(match(SpTable.ua,SpTable.all)),NA,c("Log","Logit")[sc.option[match(SpTable.ua,SpTable.all)]]),LinkSpclimSouth=NA,
              ModelNorthWinter=!is.na(match(SpTable.ua,gsub("Winter","",SpTable.w))),ModelNorthSummer=!is.na(match(SpTable.ua,gsub("Summer","",SpTable.s))),
              UseavailNorth=TRUE,UseavailSouth=NA,
              ModelSouthWinter=NA,ModelSouthSummer=NA)
fname<-paste(fname.sumout,"Header table north.csv",sep="")
write.table(file=fname,q,sep=",",row.names=FALSE)

# 11. Extra - summary of average density of on-grid sites per year
OnOff<-ifelse(substr(d$Deployment,1,1)=="A","On","Off")
FirstSpCol<-which(names(d)=="BadgerSummer")
LastSpCol<-which(names(d)=="WoodlandCaribouWinter")
dens.sum<-dens.sd<-array(NA,c(LastSpCol-FirstSpCol+1,4))  # Mean density for each species, year
d3<-d[OnOff=="On" & (d$NR=="Boreal" | d$NR=="Foothills" | d$NR=="Canadian Shield"),]
SpList<-names(d3)[FirstSpCol:LastSpCol]
lure<-NULL
for (sp in 1:length(SpList)) {
  x<-d3[,SpList[sp]]
  lure[sp]<-mean(x[d3$Study=="ABMI" & d3$Lured=="y"],na.rm=TRUE)/mean(x[d3$Study=="ABMI" & d3$Lured=="n"],na.rm=TRUE)  # Standardize to unlured (yearly lure effect)
}
for (i in 1:4) {
  d2<-d3[d3$Year==c("2015","2016","2017","2018")[i],]
  for (sp in 1:length(SpList)) {
    x<-d2[,SpList[sp]]
    x<-ifelse(d2$Lured=="y",x/lure[sp],x)
    dens.sum[sp,i]<-mean(x,na.rm=TRUE)
    dens.sd[sp,i]<-sd(x,na.rm=TRUE)
  }
}
n.per.yr<-as.numeric(table(d3$Year[d3$Year>=2015]))
dens.se<-t(t(dens.sd)/sqrt(n.per.yr))*sqrt(4)  # This is to (conservatively) account for dependence among deployments at a site
rownames(dens.sum)<-rownames(dens.se)<-names(d)[FirstSpCol:LastSpCol]
colnames(dens.sum)<-2015:2018
colnames(dens.se)<-paste("SE",2015:2018,sep="")
q<-cbind(dens.sum,dens.se)
write.table(q,file="C:/Dave/ABMI/Cameras/Coefficients/2018/Analysis north/Species density summary Boreal Foothills Shield.csv",sep=",",col.names=NA)

# 12. Extra section - produces initial grouping of sites for cross-validation.  May have to be revised by eye.
# Ermias version
D.coord <- d[ , c("Long", "TrueLat")]
names(D.coord) <- c("longitude", "latitude")
require(fossil)
# geodist <- earth.dist (D.coord)  # Too slow!
geodist<-dist(data.frame(long=111*D.coord$longitude,lat=65*D.coord$latitude))
geo.gps <- hclust (geodist , method = "ward") #I have checked other algorithms and ward seems to give a better cluster
Group.geo <- cutree (geo.gps, k=24) #k is number of clusters#
plot(d$Long,d$TrueLat,pch=18,cex=2.5,col=rainbow(24,s=0.5)[Group.geo], main="Geographical distance based groupings")
text(d$Long,d$TrueLat,d$Site,cex=0.6)
ds<-data.frame(DeploymentYear=d$DeploymentYear,Group=Group.geo)
ds<-ds[duplicated(ds$DeploymentYear)==FALSE,]
write.table(file="C:/Dave/ABMI/Data/Site info/Groups for BS and subareas/Site groupings for North mammals 24 groups Mar 2019.csv",ds,sep=",",row.names=F)





