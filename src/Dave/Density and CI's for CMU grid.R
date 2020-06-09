# R Density and CI's for CMU grid
# This version does seasons separately, then averages them

# Read data file
datafile<-"./beta/data/CMU Camera mammal data processed Mar 2019.csv"
d<-read.csv(datafile)
d<-d[-which(d$Deployment=="LRN64"),]  # Missing most images
TotalDays.all<-by(d$TotalDays,d$DeploymentYear,sum)  # Summer and winter days added for each deployment-year
failed<-names(TotalDays.all)[TotalDays.all<100]
i<-which(d$DeploymentYear %in% failed)
d<-d[-i,]  # Cameras that failed - because of winter deployments, these have no non-snow images
d$Deployment<-as.character(d$Deployment)
d$DeploymentYear<-as.character(d$DeploymentYear)
# Get rid of unnecessary columns
j<-which(regexpr(".Duration",names(d))>0)
d<-d[,-j]
j<-which(regexpr("DetDist",names(d))>0)
d<-d[,-j]

# Add combined WTD + (unknown) Deer
d$DeerAllSummer<-d$WhitetailedDeerSummer+d$DeerSummer
d$DeerAllWinter<-d$WhitetailedDeerWinter+d$DeerWinter
FirstSpCol<-which(names(d)=="BlackBearSummer")
LastSpCol<-which(names(d)=="DeerAllWinter")
SpTable<-names(d)[FirstSpCol:LastSpCol]

d$Grid<-substr(d$Deployment,1,3)
GridList<-sort(unique(d$Grid))
YearList<-c(2017,2018)

pa<-agp<-agp.se<-dens<-dens.lci<-dens.uci<-array(NA,c(length(SpTable),length(GridList),length(YearList)))
dimnames(pa)[[1]]<-dimnames(agp)[[1]]<-dimnames(dens)[[1]]<-dimnames(dens.lci)[[1]]<-dimnames(dens.uci)[[1]]<-SpTable
dimnames(pa)[[2]]<-dimnames(agp)[[2]]<-dimnames(dens)[[2]]<-dimnames(dens.lci)[[2]]<-dimnames(dens.uci)[[2]]<-GridList
dimnames(pa)[[3]]<-dimnames(agp)[[3]]<-dimnames(dens)[[3]]<-dimnames(dens.lci)[[3]]<-dimnames(dens.uci)[[3]]<-YearList
d$WinterDays<-ifelse(d$WinterDays<10,NA,d$WinterDays)  # Don't use deployments if <10 days in a season
d$SummerDays<-ifelse(d$SummerDays<10,NA,d$SummerDays)  # Don't use deployments if <10 days in a season

for (sp in 1:length(SpTable)) {
  for (grid in 1:length(GridList)) {
    for (year in 1:length(YearList)) {
      d1<-d[d$Grid==GridList[grid] & d$Year==YearList[year],]
      if (regexpr("Winter",SpTable[sp])>0) d1<-d1[!is.na(d1$WinterDays),]  # Exclude deployments with <10 days in that season
      if (regexpr("Summer",SpTable[sp])>0) d1<-d1[!is.na(d1$SummerDays),]
      if (nrow(d1)>0) {
        y<-na.omit(d1[,SpTable[sp]])  # na.omit for cameras that did not operate in that season
        pa[sp,grid,year]<-sum(sign(y))/length(y)
        agp[sp,grid,year]<-mean(y[y>0])
        agp.se[sp,grid,year]<-sd(y[y>0])/sqrt(sum(y>0))
        if (sum(y)>0) {
          dens[sp,grid,year]<-pa[sp,grid,year]*agp[sp,grid,year]
        } else {
          dens[sp,grid,year]<-0  # Get NA's in agp and therefore in dens otherwise
        }
        if (sum(sign(y))>1) {
          y1<-rbinom(10000,25,pa[sp,grid,year])/25*exp(rnorm(10000,log(agp[sp,grid,year]),sqrt(agp.se[sp,grid,year]^2/agp[sp,grid,year]^2)))  # Simulate distribution, using agp on log scale (with SE from delta method)
          y1<-y1*dens[sp,grid,year]/mean(y1)  # Compensate for any difference from mean (should be small)
          dens.lci[sp,grid,year]<-quantile(y1,0.05)
          dens.uci[sp,grid,year]<-quantile(y1,0.95)
        }
      }  # End if for survey for that grid in that year
    }  # Next year
  }  # Next grid
}  # Next sp

# And average density across seasons
SpTable1<-unique(substr(row.names(dens),1,nchar(row.names(dens))-6))
pa.avg<-agp.avg<-dens.avg<-dens.lci.avg<-dens.uci.avg<-array(NA,c(length(SpTable1),dim(dens)[2],length(YearList)))
dimnames(pa.avg)[[1]]<-dimnames(agp.avg)[[1]]<-dimnames(dens.avg)[[1]]<-dimnames(dens.lci.avg)[[1]]<-dimnames(dens.uci.avg)[[1]]<-SpTable1
dimnames(pa.avg)[[2]]<-dimnames(agp.avg)[[2]]<-dimnames(dens.avg)[[2]]<-dimnames(dens.lci.avg)[[2]]<-dimnames(dens.uci.avg)[[2]]<-dimnames(dens)[[2]]
dimnames(pa.avg)[[3]]<-dimnames(agp.avg)[[3]]<-dimnames(dens.avg)[[3]]<-dimnames(dens.lci.avg)[[3]]<-dimnames(dens.uci.avg)[[3]]<-YearList
for (sp in 1:length(SpTable1)) {
  for (year in 1:length(YearList)) {
    i<-which(regexpr(SpTable1[sp],row.names(pa))>0)
    pa.avg[sp,,year]<-colMeans(pa[i,,year])
    i<-which(regexpr(SpTable1[sp],row.names(agp))>0)
    agp.avg[sp,,year]<-colMeans(agp[i,,year],na.rm=TRUE)
    i<-which(regexpr(SpTable1[sp],row.names(dens))>0)
    dens.avg[sp,,year]<-colMeans(dens[i,,year])
    i<-which(regexpr(SpTable1[sp],row.names(dens.lci))>0)
    dens.lci.avg[sp,,year]<-colMeans(ifelse(is.na(dens.lci[i,,year]),dens[i,,year],dens.lci[i,,year]))  # This is using the mean of seasons that do not have CI's --> overly narrow CI's
    #		dens.lci.avg[sp,,year]<-colMeans(dens.lci[i,,year],na.rm=TRUE)  # This is using just the CI's on one season if the other season does not have an estimate --> CI's not centred on mean
    i<-which(regexpr(SpTable1[sp],row.names(dens.uci))>0)
    dens.uci.avg[sp,,year]<-colMeans(ifelse(is.na(dens.uci[i,,year]),dens[i,,year],dens.uci[i,,year]))  # All/most? NA's come from ~0 estimates# This is using the mean of seasons that do not have CI's --> overly narrow CI's
    #		dens.uci.avg[sp,,year]<-colMeans(dens.uci[i,,year],na.rm=TRUE)  # This is using just the CI's on one season if the other season does not have an estimate  --> CI's not centred on mean
  }
}

# Full summary
q<-data.frame(Year=2017,Sp=c(SpTable,SpTable1),Occ=rbind(pa[,,1],pa.avg[,,1]),AGP=rbind(agp[,,1]*100,agp.avg[,,1]*100),Density=rbind(dens[,,1]*100,dens.avg[,,1]*100),Density.LCI=rbind(dens.lci[,,1]*100,dens.lci.avg[,,1]*100),Density.UCI=rbind(dens.uci[,,1]*100,dens.uci.avg[,,1]*100))
q<-rbind(q,data.frame(Year=2018,Sp=c(SpTable,SpTable1),Occ=rbind(pa[,,2],pa.avg[,,2]),AGP=rbind(agp[,,2]*100,agp.avg[,,2]*100),Density=rbind(dens[,,2]*100,dens.avg[,,2]*100),Density.LCI=rbind(dens.lci[,,2]*100,dens.lci.avg[,,2]*100),Density.UCI=rbind(dens.uci[,,2]*100,dens.uci.avg[,,2]*100)))
q<-q[order(q$Year,rownames(q)),]
write.table(q,file="C:/Dave/ABMI/Cameras/2018 analysis/CMU/Density summary CMU 2017 and 2018 Seasons Mar 2019.csv",sep=",",row.names=F)
# Summary to edit for plotting
q<-data.frame(Year=2017,Sp=c(SpTable,SpTable1),Density=rbind(dens[,,1]*100,dens.avg[,,1]*100),Density.LCI=rbind(dens.lci[,,1]*100,dens.lci.avg[,,1]*100),Density.UCI=rbind(dens.uci[,,1]*100,dens.uci.avg[,,1]*100))
q<-rbind(q,data.frame(Year=2018,Sp=c(SpTable,SpTable1),Density=rbind(dens[,,2]*100,dens.avg[,,2]*100),Density.LCI=rbind(dens.lci[,,2]*100,dens.lci.avg[,,2]*100),Density.UCI=rbind(dens.uci[,,2]*100,dens.uci.avg[,,2]*100)))
q<-q[order(q$Year,rownames(q)),]
write.table(q,file="C:/Dave/ABMI/Cameras/2018 analysis/CMU/CMU densities to plot 2017 and 2018 Seasons Mar 2019.csv",sep=",",row.names=F)






