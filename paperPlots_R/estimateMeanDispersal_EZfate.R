#code to plot mean and standard deviation of particle dispersal as a function
#of dispersal time. 

library(ggplot2)
library(tidyverse)
library(furrr)

#source routines to manipulate connectivity. Which one to use depends on
#which machine is being used...
source('../../../oceanparcels/github_4_EZfate/EZfate/code/connectivityUtilities.R')
#source('/Users/pringle/Documents/GitHub/EZfate/code/connectivityUtilities.R')
#source('/Users/pringle/Documents/GitHub/EZfate/code/connectivityUtilities_parallelSubset.R')

#===============================================================================
#write function to trim to habitat
trim2habitat<-function(Ein,depthMin,depthMax){
  
  if (length(latPoly)>0){
    limitPoly=st_polygon(list(cbind(lonPoly,latPoly)))
    limitPoly<-st_sfc(limitPoly,crs=4326) #make into a surface, 4326 is the code for the WGS84 datum in proj
    
    #subset to polygon
    Eout<-subsetConnectivity_byPolygon(Ein,limitPoly,trimTo=trimTo)
  } else {
    #nothing in latPoly, so don't subset
    Eout<-Ein
  }
  
  #now subset to locations in shallow water
  Eout<-subsetConnectivity_byGridValue(Eout,gridDepth,depthMin,depthMax,trimTo=trimTo) #water depth < 25 m
  
  return(Eout)
}

#===============================================================================
library('geosphere')
library('testit')
#write function to calculate Ladv and Ldiff from E
#in answer include initial lat, lon, Ladv, Ldiff
meanSquareDist <-function(lonLatPnt,lonVec,latVec,numVec) {
  #this function returns the mean of distances squared between lonLatPnt 
  #(two element vector, c(longitude,latitude) in degrees)
  #degrees) and all the points in the vectors lonVec and latVec. Units are km^2
  
  lonPnt=lonLatPnt[1]
  latPnt=lonLatPnt[2]
  
  assert('length of lonVec and latVec should be equal',length(lonVec)==length(latVec))
  lonLatMat=matrix(0,nrow=length(lonVec),ncol=2)
  lonLatMat[,1]=lonVec
  lonLatMat[,2]=latVec
  
  dist=distHaversine(lonLatMat,c(lonPnt,latPnt),r=6378.137)
  
  #must upweight the dist^2 by the number of particles that make that trip
  meanSqDist=sum((dist^2)*numVec)/(sum(numVec)) #biased estamator, but won't die if N=1
  
  #for debugging
  #print(c('Point is',c(lonPnt,latPnt),meanSqDist))
  
  return(meanSqDist)
}

findMeanAndSTD<-function(lonVec,latVec,numVec){
  #this function takes two vectors of equal size, the first of longitudes,
  #the second of latitudes, and returns a vector of the point that minimizes
  #the sum of square distances to all the points (the "mean" location) and 
  #the sqrt(sum of square distances) to all the points from this minimum. The
  #later is equivalent to a standard deviation.
  #
  #It does this by minimizing the function meanSquareDist as a function of location of 
  #the mean point to find the mean point
  
  #sign, since latVecIn and lonVecIn are vector in data frame, I need to unpack
  #them just to get what I want
  #latVec<-latVecIn[[1]]; lonVec<-lonVecIn[[1]]
  
  assert('length of lonVec and latVec should be equal',length(lonVec)==length(latVec))
  
  #find the lon and lat to start the optimization at by taking mean 
  #of lat and lon vectors
  latStart=mean(latVec)
  lonStart=mean(lonVec)
  
  #for debugging
  #print(paste(lonStart,latStart,length(lonVec)))
  
  if (length(lonVec)>0) {
    #for syntax of optim, refer to 
    #https://stackoverflow.com/questions/24623488/how-do-i-use-a-function-with-parameters-in-optim-in-r
    result=optim(par=c(lonStart,latStart),fn=meanSquareDist,lonVec=lonVec,latVec=latVec,numVec=numVec)
    return(list(meanPoint=result$par,stdDist=sqrt(result$value)))
  } else {
    return(list(meanPoint=c(NaN,NaN),stdDist=NaN))
  }
  
}

findLadvLdiff<-function(lonStart,latStart,lonVec,latVec,numVec,numLaunched){
  #this code takes a starting position (lonStart,latStart) and a vector of ending
  #longitudes lonVec and latitudes latVec, and computes the standard deviation
  #in distance of lonVec,latVec, and the distance from (lonStart,latStart) to the
  #average location of lonVec and latVec

  #calculate mean and std of lonVec,latVec
  jnk<-findMeanAndSTD(lonVec,latVec,numVec)
  Ldiff<-jnk$stdDist
  meanPoint<-jnk$meanPoint; meanLon<-meanPoint[1]; meanLat<-meanPoint[2]
  Ladv<-distHaversine(c(meanLon,meanLat),c(lonStart,latStart),r=6378.137) #in km
  fracReturn<-sum(numVec)/numLaunched
  
  return(data.frame(lonFrom=lonStart,latFrom=latStart,Ladv=Ladv,Ldiff=Ldiff,fracReturn))
}

#For a given PLD, calculate Ladv and Ldiff for all points 
#in the polygon described by latPoly on lonPoly. First get all the data

#what PLD to get?
PLD=16 # in days

#get data to analyze
regionName<-'theAmericas'
depth<-1
year<-'climatology'
verticalBehavior<-'fixed'
months<-c(4,5,6)
minPLD<-PLD; maxPLD<-minPLD

tic(paste('working on connectivity data for PLD=',PLD,'in'))
print(paste('...getting month',months[1]))
E<-getConnectivityData(regionName,depth,year,verticalBehavior,months[1],minPLD,maxPLD)
for (month in months[2:length(months)]){
  print(paste('...getting and combining month',month))
  E<-combineConnectivityData(E,getConnectivityData(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD))
}
toc()
Estart<-E 

#now iterate over the regions we are analyzing
dispersalData=data.frame()
for (whichCase in c(1,2)) {
  
  E<-Estart #we will trim E below, so we need to start with a clean copy
  
  #set up cases 
  if (whichCase==1) {
    #define what region and times to analyze
    #define region to be subset to
    caseName<-'East Coast USA'
    trimTo<-TRUE
    lonMin=-81.5; lonMax=-59.0 ; latMin=25.0 ; latMax=46.0
    latPoly<-c(latMin,latMax,latMax,latMin,latMin)
    lonPoly<-c(lonMin,lonMin,lonMax,lonMax,lonMin)
  } else {
    #define what region and times to analyze
    #define region to be subset to
    caseName<-'West Coast USA'
    trimTo<-TRUE
    lonMin=-128.5; lonMax=-116.0 ; latMin=32.5 ; latMax=49.0
    latPoly<-c(latMin,latMax,latMax,latMin,latMin)
    lonPoly<-c(lonMin,lonMin,lonMax,lonMax,lonMin)
  }
  
  print(paste('working on',caseName))
  tic(paste('subset connectivity data for PLD=',PLD,'in'))
  #subset to region. WATCH OUT, ISOBATH RANGE OF RELEASE IS ENCODED HERE.
  #add lat lon
  E<-addLatLon(E)
  print('WARNING, SUBSETTING TO DEPTH, NOT DISTANCE FROM SHORE!')
  E<-trim2habitat(E,depthMin=0.0, depthMax=25.0) #for 1m release
  toc()
  
  #apply findLadvLdiff() to each row of E
  tic(paste('calculated Ladv and Ldiff for PLD=',PLD,'in'))
  if (TRUE){
    #if TRUE, run in parallel for speed
    dispersalDataNow<-future_pmap_dfr(E[c('lonFrom','latFrom','lonTo','latTo','numTo','numLaunched')],function(lonFrom,latFrom,lonTo,latTo,numTo,numLaunched){findLadvLdiff(lonFrom,latFrom,lonTo,latTo,numTo,numLaunched)})
  } else {
    #else if FALSE run in serial, to make debugging easier
    dispersalDataNow<-pmap_dfr(E[c('lonFrom','latFrom','lonTo','latTo','numTo','numLaunched')],function(lonFrom,latFrom,lonTo,latTo,numTo,numLaunched){findLadvLdiff(lonFrom,latFrom,lonTo,latTo,numTo,numLaunched)})
  }
  toc()
  
  #add PLD and caseName to dispersalData, 
  dispersalDataNow$PLD<-PLD
  dispersalDataNow$caseName<-caseName
  
  #add to dispersalData
  dispersalData<-rbind(dispersalData,dispersalDataNow)
  
}

#calculate mean velocity in cm/s by converting Ladv from km to meters and dividing by PLD in seconds
dispersalData$meanSpeed<-(dispersalData$Ladv*1000*100)/(dispersalData$PLD*86400)

#get case means
jnk<-dispersalData %>% group_by(caseName) 
caseMeans<-jnk %>% summarize(meanV=mean(meanSpeed,na.rm=TRUE))
caseSTD<-jnk %>% summarize(stdV=sd(meanSpeed,na.rm=TRUE))
print(caseMeans)
print(caseSTD)

#make histogram of velocity
jnk<-ggplot(dispersalData,aes(x=meanSpeed,fill=caseName))+geom_density(
                                                                       na.rm=TRUE,position='identity',alpha=0.5)+
  geom_vline(data=caseMeans,aes(xintercept=meanV,color=caseName),linetype='dashed',linewidth=1)+
  xlab('speed, cm/s')+ggtitle(paste('Distribution of mean speed for PLD=',PLD,'days'))
print(jnk)


