# example for calculating neibourhood (DBH/distance) from tree census data
# load packages
library(bigmemory)
library(bigtabulate)
library(doParallel)
library(foreach)
library(ggplot2)
library(ggforce)
library(here)

#define all functions
source(here("HOI_cal_allpair.R"))

#load data
load(here("cbstest.RData"))

#before running, check the spatial distribution, how many individuals and species are within focals' radius  
ggplot(cbstest)+geom_point(aes(x=gx,y=gy,color=Accepted_species_tnrs))+geom_circle(aes(x0=gx,y0=gy,r=10,color=Accepted_species_tnrs))+coord_fixed()#10m
ggplot(cbstest)+geom_point(aes(x=gx,y=gy,color=Accepted_species_tnrs))+geom_circle(aes(x0=gx,y0=gy,r=20,color=Accepted_species_tnrs))+coord_fixed()#20m
ggplot(cbstest)+geom_point(aes(x=gx,y=gy,color=Accepted_species_tnrs))+geom_circle(aes(x0=gx,y0=gy,r=50,color=Accepted_species_tnrs))+coord_fixed()#50m

#make a list when you put more sites
plotl<-list(cbs=cbstest)
for(i in c("cbs")){
  di.cal(plotname=i,plotlist=plotl)#here I run it in the radius of 10, 20, 50 meters as well;
}

#check the ouput csv files in the cbs/out 
