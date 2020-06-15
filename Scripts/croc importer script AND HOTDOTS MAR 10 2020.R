library(tidyverse)
library(Morpho)
library(geomorph)
library(Rvcg)
library(paleomorph)
library(EMMLi)
library(qgraph)
library(ape)
library(geiger)
library(abind)
#load the script for resampling functions
source("~/Google Drive/utlity_scripts/resamplingsV3.R")
#set the working directory
setwd("~/Dropbox/Crocs/pts3")


#curvedata must be a csv with columns 'lm1' 'lm2' and 'ptswanted'
curvedata<-read_csv("~/Google Drive/NHM/crocs/data/croc_curves_july.csv")
#list which landmarks are fixed and will not slide
fixed <- c(1:102)#load in the list of pts files from the current working directory
ptslist <- dir(pattern='.pts', recursive=F)  
pts_tibble <- as_tibble()#initialize an empty tibble to store the pts data
filenames <- ptslist
ntaxa <- length(filenames) #define the number of taxa based on how many pts files are in the working directory

#define the curves using the 'curvedata' object
subsampled.curve<-sapply(paste("SC",c(1:nrow(curvedata)),sep=""),function(x) NULL)
#subsampled.curve <- list()
subsampled.curve[[1]]<-as.integer(c(curvedata$lm1[1],((length(fixed)+1):(length(fixed)+curvedata$ptswanted[1])),curvedata$lm2[1]))
for (i in 2:length(subsampled.curve)){
  subsampled.curve[[i]]<-as.integer(c(curvedata$lm1[i],((max(unlist(subsampled.curve))+1):(max(unlist(subsampled.curve))+curvedata$ptswanted[i])),curvedata$lm2[i]))
}

subsampled.curve.in<-subsampled.curve
for(i in 1:length(subsampled.curve.in)){
  subsampled.curve.in[[i]]<-tail(head(subsampled.curve[[i]],-1),-1)
}


slidings.sub<-c((max(fixed)+1):max(unlist(subsampled.curve)))


## IMPORT DATA ##
{
  for(i in 1:length(ptslist))
  {
    specimen.tmp <- as_tibble(read.table(file=ptslist[i],skip=2,header=F,sep="")) #import a single specimen
    specimen.tmp <- specimen.tmp %>% mutate(., V1 = as.character(V1)) #convert the first row, the lm names, to characters instead of factors
    specimen.tmp <- specimen.tmp %>% mutate(.,spec.id=filenames[i]) #add a column with the specimen name
    pts_tibble<-bind_rows(pts_tibble,specimen.tmp) #paste it to the end of the tibble with the rest of the specimens 
  }
  #this will give a warning message but its nothing to worry about
  pts_tibble <- pts_tibble %>% separate(., V1, into = c("class", "id","sub_lm"), remove = FALSE)
  
  pts_tibble <- pts_tibble %>% mutate(., id = as.factor(id)) %>%
    rename(., index = V1, X = V2, Y = V3, Z = V4) #rename the coordinate data columns
  print("SERIOUSLY DONT WORRY ABOUT THE WARNINGS BELOW")
}

#make a list of how many sliding semilandmarks on curves there are
curvepoints <- sum(curvedata$ptswanted)
#convert the tibble to a 3D array compatable with geomorph
pts_tibble_tmp <- pts_tibble%>%filter(.,class=="S")%>%group_by(spec.id)%>%select(.,X,Y,Z)%>%nest()%>%transpose()
ptsarray_tmp <- array(dim=c(length(fixed),3,ntaxa))



#####some code to check whether you have any specimens that have the wrong number of landmarks
test<-lapply(1:length(pts_tibble_tmp), function(k) dim(pts_tibble_tmp[[k]]$data)[1]) %>% unlist(.)
print("The following specimens have the wrong number of fixed points in them:")
print(filenames[which(test!=length(fixed))])


test2<-pts_tibble%>%filter(., class=="C")%>%group_by(spec.id)%>%select(.,id,X,Y,Z)%>%mutate(., id = as.numeric(as.character(id)))%>%nest()%>%transpose()
test2a<- lapply(1:length(test2), function(k) length(unique(test2[[k]]$data$id))) %>% unlist(.)
print("The following specimens have the wrong number of curves in them:")
filenames[which(test2a!=length(curvedata$Curves))]

for(i in 1:length(pts_tibble_tmp))
{
  ptsarray_tmp[,,i] <- as.matrix(select(pts_tibble_tmp[[i]]$data, c(X,Y,Z)))
}



##############################################################
####making all the curves have the correct number of landmarks
##############################################################
#make an empty array with the correct number of landmarks and specimens 
newpts <- array(data = NA, dim = c(length(fixed),3,ntaxa))
#give it dimension names based on your specimen list
dimnames(newpts)[3] <- list(substr(filenames,1,(nchar(filenames)-4)))
#fill in the fixed landmarks 
newpts[fixed,,] <- ptsarray_tmp

#OPTIONAL:
#check for outliers in your anatomical landmarks before subsampling
#find.outliers(newpts)

#right now I have this set to do 1 through 69 becuase there are 69 curves we want to use
#but in theory it should be 1:nrow(curvedata) to do this for all curves
for (which.curve in 1:nrow(curvedata)){
  this.curve <- array(data=NA, dim=c(curvedata$ptswanted[which.curve],3,ntaxa))
  for (which.spec in 1:length(filenames)){
    orig.curve <- pts_tibble %>% filter(.,spec.id==filenames[which.spec])%>%filter(., class=="C")%>%filter(., id==which.curve) %>% select(., X,Y,Z)
    orig.curve.anchors <- pts_tibble %>% filter(.,spec.id==filenames[which.spec])%>%slice(c(subsampled.curve[[which.curve]][1],last(subsampled.curve[[which.curve]]))) %>% select(., X,Y,Z)
    orig.curve <- rbind(orig.curve.anchors[1,],orig.curve,orig.curve.anchors[2,])
    new.curve <- cursub.interpo(orig.curve, curvedata$ptswanted[which.curve])
    this.curve[,,which.spec] <- as.matrix(new.curve)[2:(dim(new.curve)[1]-1),]
  }
  newpts <- abind(newpts, this.curve, along=1)
}




#separate your template points from the actual specimens
subsampled.lm <- newpts[,,-which(dimnames(newpts)[[3]]=="hemispheremodelHIGH2")]
template.lm <- newpts[,,which(dimnames(newpts)[[3]]=="hemispheremodelHIGH2")]
ntaxa <- ntaxa-1 #remove the hemisphere model from the taxon list


#OPTIONAL
#define which curves belong to which modules for a pretty plot
#curvemodules<-c()
#for (i in 1:nrow(curvedata)) {
#  x<-rep(curvedata$Colors[i],curvedata$ptswanted[i])
#  curvemodules<-c(curvemodules,x)
#}
#curvemodules<-as.factor(curvemodules)
#levels(curvemodules)<-c(1:21)
#curvemodules<-c(rep("black",length(fixed)),curvemodules)




#if you have any missing points, Checkpoint will set them to x,y,z=9999
#this makes for bad plots. switch them to NA
subsampled.lm[which(subsampled.lm==9999)]<-NA
subsampled.lm[which(is.nan(subsampled.lm))]<-NA
subsampled.lm[which(subsampled.lm==-Inf)]<-NA



#anatosuchus is the wrong size, scale it correctly please:
subsampled.lm[,,which(dimnames(subsampled.lm)[[3]]=="Anatoschus_Merged_Mirrored")]<-subsampled.lm[,,which(dimnames(subsampled.lm)[[3]]=="Anatoschus_Merged_Mirrored")]/10




#this is the landmark numbers which are patch points on your template in checkpoint
checkpoint.patch.nums<-c(103:699)

#loading your template
{
  template.lm.and.patch <- read_csv("/Volumes/R Felice 4T/crocs/updated template files/hemisphere template 2 CSV.csv", col_names = TRUE, skip = 1)
  template.lm.and.patch <- template.lm.and.patch[,4:6]
  template.mesh <- vcgImport(file="~/Dropbox/Crocs/pts3/ply/hemispheremodelHIGH.ply")# load surface scan from the template
  #####IMPORTANT:
  #####replace 
  patch <- as.matrix(template.lm.and.patch[checkpoint.patch.nums,])
}
# creation of a vector containing rowindices of sliding landmarks of surface


## CREATION of the *subsampled* ATLAS

atlas.subsampl<-createAtlas(template.mesh,as.matrix(template.lm),as.matrix(patch), corrCurves=subsampled.curve.in, patchCurves = NULL, keep.fix = fixed)
total.lm.num<-dim(template.lm)[[1]]+dim(patch)[[1]]


#Plot the atlas
open3d()
shade3d(atlas.subsampl$mesh,col="white")
spheres3d(atlas.subsampl$landmarks[fixed,],col="red",radius=0.04)
spheres3d(atlas.subsampl$landmarks[slidings.sub,],col="gold",radius=0.04)
spheres3d(atlas.subsampl$patch,col=4,radius=0.04)

#check LM again
#checkLM(subsampled.lm,atlas=atlas.subsampl,begin=42,pt.size = 2,path="./ply/",suffix=".ply",render="s")


#load Aki's function for redefining curves
ReCorrCurves <- function(curve.in, curvein.nos, rec1, n.fixed) {
  no.lm <- length(unique(unlist(subsampled.curve.in[curvein.nos])))
  new.no <- as.integer(seq(1:no.lm)+n.fixed)
  new.curvein <- list()
  counter <- 1
  for(i in curvein.nos) {
    new.curvein[[counter]] <- new.no[((sum(rec1[curvein.nos[1:counter]])-(rec1[i]-1)):sum(rec1[curvein.nos[1:counter]]))]
    counter <- counter + 1
  }
  return(new.curvein)
}


{
    premax_vent <- c(1,2,3)
    
    max_vent <- c(4,5,6)
    
    max_dorsal <- c(7,8,9)
    
    premax_dorsal <- c(10,11,12)
    
    nasal <- c(13,14,15)
    
    frontal <- c(17,18,19,20)-1
    
    parietal <- c(21,22,23,24)-1
    
    exocc <- c(25,26,27,28)-1
    
    occ.condyle <- c(29,30,31)-1
    
    basiocc <- c(32,33)-1
    
    basisph <- c(34,35,36)-1
    
    pterygoid <- c(37,38,39,40)-1
    
    palatine <- c(41,42,43,44)-1
    
    jugal <- c(45,46,47,48,49)-1
    
    squamosal <- c(50,51,52,53)-1
    
    lacrimal <- c(54,55,56,57)-1
    
    ecto.pter <- c(58,59,60,61)-1
    
    ptery.joint <- c(62,63)-1
    
    jawjoint <- c(64,65,66)-1
    
    postorb <- c(67:71)-1
    
}

#define the patches:
module_defs <- read_csv("~/Google Drive/NHM/crocs/data/croc module defs3-aug7.csv")
completeness <- read_csv("~/Google Drive/NHM/crocs/data/croc_module_presence_July_2019.csv")
completeness<-completeness %>% arrange(., Filename)


color.palette <- c( "#6a3d9a","dimgrey","#fb9a99",  "gold", "#009E73",  "#D55E00", "#CC79A7", "cyan2",  "#e31a1c", "#0072B2", "#b2df8a", "#E69F00",  "whitesmoke" ,  "deeppink",   "#a6cee3",   "#F0E442","blue","red","brown", "black")
color_table <- read_csv("~/Google Drive/_UCL/Projects/Archosaur Modularity/Colortable.csv")



rec1 <- rep(10,length(subsampled.curve))
rec1[c(5,6,7,36,46)] <- 20
rec1[c(14,21,28,30,32,37,43,45,50,52,54,63,64,66,67,68,70)] <- 5
rec1 <- curvedata$ptswanted

modulecolors<-as.factor(color_table$croc_colors)
levels(modulecolors)<-color.palette



lm_premax_d <- subsampled.curve[c(nasal, premax_dorsal)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_premax_d <- module_defs %>% filter(., description=="patch_premax_d") %>% pull(., lm)-102

atlas_premax_d<- createAtlas(mesh = template.mesh, 
                             landmarks = as.matrix(template.lm[lm_premax_d,]),
                             patch = as.matrix(patch[patch_premax_d,]),
                             corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                     curvein.nos = c(nasal, premax_dorsal), 
                                                     rec1 = rec1,
                                                     n.fixed = length(lm_premax_d[lm_premax_d<103])), 
                             patchCurves = NULL,
                             keep.fix=1:length(lm_premax_d[lm_premax_d<103]))
patched_premax_d <- placePatch(atlas_premax_d,
                               subsampled.lm[lm_premax_d,,which(completeness$premax_d==1)],
                               path="./ply",
                               prefix=NULL,
                               fileext=".ply",
                               ray=TRUE,
                               inflate=2,
                               tol=.8,
                               #relax.patch=FALSE,
                               relax.patch=TRUE,
                               keep.fix=1:length(lm_premax_d[lm_premax_d<103]),
                               #rhotol=NULL,
                               silent=FALSE, mc.cores=3)

#checkLM(patched_premax_d,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_premax_d,alpha=1,begin=1)

lm_premax_d2 <- subsampled.curve[c(nasal, premax_dorsal, max_vent, premax_vent, max_dorsal)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_premax_d2 <- module_defs %>% filter(., description=="patch_premax_d") %>% pull(., lm)-102

atlas_premax_d2 <- createAtlas(mesh = template.mesh, 
                             landmarks = as.matrix(template.lm[lm_premax_d2,]),
                             patch = as.matrix(patch[patch_premax_d2,]),
                             corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                     curvein.nos = c(nasal, premax_dorsal, max_vent, premax_vent, max_dorsal), 
                                                     rec1 = rec1,
                                                     n.fixed = length(lm_premax_d2[lm_premax_d2<103])), 
                             patchCurves = NULL,
                             keep.fix=1:length(lm_premax_d2[lm_premax_d2<103]))

patched_premax_d2 <- placePatch(atlas_premax_d2,
                               subsampled.lm[lm_premax_d2,,which(completeness$premax_d==1)],
                               path="./ply",
                               prefix=NULL,
                               fileext=".ply",
                               ray=TRUE,
                               inflate=2,
                               tol=.8,
                               #relax.patch=FALSE,
                               relax.patch=TRUE,
                               keep.fix=1:length(lm_premax_d2[lm_premax_d2<103]),
                               #rhotol=NULL,
                               silent=FALSE, mc.cores=3)
#checkLM(patched_premax_d2,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_premax_d2,alpha=1,begin=1)



lm_premax_d3 <- subsampled.curve[c(premax_dorsal)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_premax_d3 <- module_defs %>% filter(., description=="patch_premax_d") %>% pull(., lm)-102

atlas_premax_d3 <- createAtlas(mesh = template.mesh, 
                               landmarks = as.matrix(template.lm[lm_premax_d3,]),
                               patch = as.matrix(patch[patch_premax_d3,]),
                               corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                       curvein.nos = c(premax_dorsal), 
                                                       rec1 = rec1,
                                                       n.fixed = length(lm_premax_d3[lm_premax_d3<103])), 
                               patchCurves = NULL,
                               keep.fix=1:length(lm_premax_d3[lm_premax_d3<103]))

patched_premax_d3 <- placePatch(atlas_premax_d3,
                                subsampled.lm[lm_premax_d3,,which(completeness$premax_d==1)],
                                path="./ply",
                                prefix=NULL,
                                fileext=".ply",
                                ray=TRUE,
                                inflate=2,
                                tol=.8,
                                #relax.patch=FALSE,
                                relax.patch=TRUE,
                                keep.fix=1:length(lm_premax_d3[lm_premax_d3<103]),
                                #rhotol=NULL,
                                silent=FALSE, mc.cores=3) 
checkLM(patched_premax_d3,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_premax_d3,alpha=1,begin=33)

#################################
lm_max_d <- subsampled.curve[c(max_dorsal,nasal,max_vent,jugal)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_max_d <- module_defs %>% filter(., description=="patch_max_d") %>% pull(., lm)-102

atlas_max_d<- createAtlas(mesh = template.mesh, 
                             landmarks = as.matrix(template.lm[lm_max_d,]),
                             patch = as.matrix(patch[patch_max_d,]),
                             corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                     curvein.nos = c(max_dorsal,nasal,max_vent,jugal), 
                                                     rec1 = rec1,
                                                     n.fixed = length(lm_max_d[lm_max_d<103])), 
                             patchCurves = NULL,
                             keep.fix=1:length(lm_max_d[lm_max_d<103]))
patched_max_d <- placePatch(atlas_max_d,
                               subsampled.lm[lm_max_d,,which(completeness$max_d ==1)],
                               path="./ply",
                               prefix=NULL,
                               fileext=".ply",
                               ray=TRUE,
                               inflate=4,
                               tol=2,
                               #relax.patch=FALSE,
                               relax.patch=TRUE,
                               keep.fix=1:length(lm_max_d[lm_max_d<103]),
                               #rhotol=NULL,
                               silent=FALSE, mc.cores=3)

checkLM(patched_max_d,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_max_d,alpha=1,begin=33)


lm_max_d2 <- subsampled.curve[c(max_dorsal,premax_dorsal)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_max_d2 <- module_defs %>% filter(., description=="patch_max_d") %>% pull(., lm)-102

atlas_max_d2<- createAtlas(mesh = template.mesh, 
                          landmarks = as.matrix(template.lm[lm_max_d2,]),
                          patch = as.matrix(patch[patch_max_d2,]),
                          corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                  curvein.nos = c(max_dorsal,premax_dorsal), 
                                                  rec1 = rec1,
                                                  n.fixed = length(lm_max_d2[lm_max_d2<103])), 
                          patchCurves = NULL,
                          keep.fix=1:length(lm_max_d2[lm_max_d2<103]))

#just for spec 24 whcih has no max vent
patched_max_d2 <- placePatch(atlas_max_d2,
                            subsampled.lm[lm_max_d2,,which(completeness$max_d==1)],
                            path="./ply",
                            prefix=NULL,
                            fileext=".ply",
                            ray=TRUE,
                            inflate=2,
                            tol=.8,
                            #relax.patch=FALSE,
                            relax.patch=TRUE,
                            keep.fix=1:length(lm_max_d2[lm_max_d2<103]),
                            #rhotol=NULL,
                            silent=FALSE, mc.cores=3)

checkLM(patched_max_d2,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_max_d2,alpha=1,begin=10)
#############
lm_max_d3 <- subsampled.curve[c(max_dorsal, max_vent)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_max_d3 <- module_defs %>% filter(., description=="patch_max_d") %>% pull(., lm)-102

atlas_max_d3<- createAtlas(mesh = template.mesh, 
                           landmarks = as.matrix(template.lm[lm_max_d3,]),
                           patch = as.matrix(patch[patch_max_d3,]),
                           corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                   curvein.nos = c(max_dorsal, max_vent), 
                                                   rec1 = rec1,
                                                   n.fixed = length(lm_max_d3[lm_max_d3<103])), 
                           patchCurves = NULL,
                           keep.fix=1:length(lm_max_d3[lm_max_d3<103]))

#just for spec 24 whcih has no max vent
patched_max_d3 <- placePatch(atlas_max_d3,
                             subsampled.lm[lm_max_d3,,which(completeness$max_d==1)],
                             path="./ply",
                             prefix=NULL,
                             fileext=".ply",
                             ray=TRUE,
                             inflate=2,
                             tol=.8,
                             #relax.patch=FALSE,
                             relax.patch=TRUE,
                             keep.fix=1:length(lm_max_d3[lm_max_d3<103]),
                             #rhotol=NULL,
                             silent=FALSE, mc.cores=3)

checkLM(patched_max_d3,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_max_d3,alpha=1,begin=33)

#################################
lm_max_vent <- subsampled.curve[c(nasal, max_vent)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_max_vent <- module_defs %>% filter(., description=="patch_max_vent") %>% pull(., lm)-102

atlas_max_vent<- createAtlas(mesh = template.mesh, 
                          landmarks = as.matrix(template.lm[lm_max_vent,]),
                          patch = as.matrix(patch[patch_max_vent,]),
                          corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                  curvein.nos = c(nasal, max_vent), 
                                                  rec1 = rec1,
                                                  n.fixed = length(lm_max_vent[lm_max_vent<103])), 
                          patchCurves = NULL,
                          keep.fix=1:length(lm_max_vent[lm_max_vent<103]))
patched_max_vent <- placePatch(atlas_max_vent,
                            subsampled.lm[lm_max_vent,,which(completeness$max_vent==1)],
                            path="./ply",
                            prefix=NULL,
                            fileext=".ply",
                            ray=TRUE,
                            inflate=2,
                            tol=5,
                            #relax.patch=FALSE,
                            relax.patch=TRUE,
                            keep.fix=1:length(lm_max_vent[lm_max_vent<103]),
                            #rhotol=NULL,
                            silent=FALSE, mc.cores=3)

#checkLM(patched_max_vent,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_premax_d,alpha=1,begin=1)

######################
lm_premax_vent <- subsampled.curve[c(nasal, premax_vent)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_premax_vent <- module_defs %>% filter(., description=="patch_premax_vent") %>% pull(., lm)-102

atlas_premax_vent<- createAtlas(mesh = template.mesh, 
                             landmarks = as.matrix(template.lm[lm_premax_vent,]),
                             patch = as.matrix(patch[patch_premax_vent,]),
                             corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                     curvein.nos = c(nasal, premax_vent), 
                                                     rec1 = rec1,
                                                     n.fixed = length(lm_premax_vent[lm_premax_vent<103])), 
                             patchCurves = NULL,
                             keep.fix=1:length(lm_premax_vent[lm_premax_vent<103]))
patched_premax_vent <- placePatch(atlas_premax_vent,
                               subsampled.lm[lm_premax_vent,,which(completeness$premax_vent==1)],
                               path="./ply",
                               prefix=NULL,
                               fileext=".ply",
                               ray=TRUE,
                               inflate=2,
                               tol=5,
                               #relax.patch=FALSE,
                               relax.patch=TRUE,
                               keep.fix=1:length(lm_premax_vent[lm_premax_vent<103]),
                               #rhotol=NULL,
                               silent=FALSE, mc.cores=3)


######################
lm_nasal <- subsampled.curve[c(nasal)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_nasal <- module_defs %>% filter(., description=="patch_nasal") %>% pull(., lm)-102

atlas_nasal<- createAtlas(mesh = template.mesh, 
                             landmarks = as.matrix(template.lm[lm_nasal,]),
                             patch = as.matrix(patch[patch_nasal,]),
                             corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                     curvein.nos = c(nasal), 
                                                     rec1 = rec1,
                                                     n.fixed = length(lm_nasal[lm_nasal<103])), 
                             patchCurves = NULL,
                             keep.fix=1:length(lm_nasal[lm_nasal<103]))
patched_nasal <- placePatch(atlas_nasal,
                               subsampled.lm[lm_nasal,,which(completeness$nasal==1)],
                               path="./ply",
                               prefix=NULL,
                               fileext=".ply",
                               ray=TRUE,
                               inflate=2,
                               tol=.8,
                               #relax.patch=FALSE,
                               relax.patch=TRUE,
                               keep.fix=1:length(lm_nasal[lm_nasal<103]),
                               #rhotol=NULL,
                               silent=FALSE, mc.cores=3)

checkLM(patched_nasal,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_nasal,alpha=1,begin=1)


lm_nasal2 <- subsampled.curve[c(premax_vent, max_vent, nasal)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_nasal2 <- module_defs %>% filter(., description=="patch_nasal") %>% pull(., lm)-102

atlas_nasal2<- createAtlas(mesh = template.mesh, 
                          landmarks = as.matrix(template.lm[lm_nasal2,]),
                          patch = as.matrix(patch[patch_nasal2,]),
                          corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                  curvein.nos = c(premax_vent, max_vent, nasal), 
                                                  rec1 = rec1,
                                                  n.fixed = length(lm_nasal2[lm_nasal2<103])), 
                          patchCurves = NULL,
                          keep.fix=1:length(lm_nasal2[lm_nasal2<103]))
patched_nasal2 <- placePatch(atlas_nasal2,
                            subsampled.lm[lm_nasal2,,which(completeness$nasal==1)],
                            path="./ply",
                            prefix=NULL,
                            fileext=".ply",
                            ray=TRUE,
                            inflate=2,
                            tol=.8,
                            relax.patch=FALSE,
                            #relax.patch=TRUE,
                            keep.fix=1:length(lm_nasal2[lm_nasal2<103]),
                            #rhotol=NULL,
                            silent=FALSE, mc.cores=3)

checkLM(patched_nasal2,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_nasal2,alpha=1,begin=33)

lm_nasal3 <- subsampled.curve[c(premax_dorsal, max_dorsal, nasal)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_nasal3 <- module_defs %>% filter(., description=="patch_nasal") %>% pull(., lm)-102

atlas_nasal3<- createAtlas(mesh = template.mesh, 
                           landmarks = as.matrix(template.lm[lm_nasal3,]),
                           patch = as.matrix(patch[patch_nasal3,]),
                           corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                   curvein.nos = c(premax_dorsal, max_dorsal, nasal), 
                                                   rec1 = rec1,
                                                   n.fixed = length(lm_nasal3[lm_nasal3<103])), 
                           patchCurves = NULL,
                           keep.fix=1:length(lm_nasal3[lm_nasal3<103]))
patched_nasal3 <- placePatch(atlas_nasal3,
                             subsampled.lm[lm_nasal3,,which(completeness$nasal==1)],
                             path="./ply",
                             prefix=NULL,
                             fileext=".ply",
                             ray=TRUE,
                             inflate=2,
                             tol=.8,
                             relax.patch=FALSE,
                             #relax.patch=TRUE,
                             keep.fix=1:length(lm_nasal3[lm_nasal3<103]),
                             #rhotol=NULL,
                             silent=FALSE, mc.cores=3)

checkLM(patched_nasal3,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_nasal3,alpha=1,begin=33)

######################
lm_frontal <- subsampled.curve[c(frontal)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_frontal <- module_defs %>% filter(., description=="patch_frontal") %>% pull(., lm)-102

atlas_frontal<- createAtlas(mesh = template.mesh, 
                          landmarks = as.matrix(template.lm[lm_frontal,]),
                          patch = as.matrix(patch[patch_frontal,]),
                          corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                  curvein.nos = c(frontal), 
                                                  rec1 = rec1,
                                                  n.fixed = length(lm_frontal[lm_frontal<103])), 
                          patchCurves = NULL,
                          keep.fix=1:length(lm_frontal[lm_frontal<103]))
patched_frontal <- placePatch(atlas_frontal,
                            subsampled.lm[lm_frontal,,which(completeness$frontal==1)],
                            path="./ply",
                            prefix=NULL,
                            fileext=".ply",
                            ray=TRUE,
                            inflate=4,
                            tol=3,
                            #relax.patch=FALSE,
                            relax.patch=TRUE,
                            keep.fix=1:length(lm_frontal[lm_frontal<103]),
                            #rhotol=NULL,
                            silent=FALSE, mc.cores=3)




######################
lm_parietal <- subsampled.curve[c(parietal,premax_dorsal)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_parietal <- module_defs %>% filter(., description=="patch_parietal") %>% pull(., lm)-102

atlas_parietal<- createAtlas(mesh = template.mesh, 
                          landmarks = as.matrix(template.lm[lm_parietal,]),
                          patch = as.matrix(patch[patch_parietal,]),
                          corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                  curvein.nos = c(parietal,premax_dorsal), 
                                                  rec1 = rec1,
                                                  n.fixed = length(lm_parietal[lm_parietal<103])), 
                          patchCurves = NULL,
                          keep.fix=1:length(lm_parietal[lm_parietal<103]))
patched_parietal <- placePatch(atlas_parietal,
                            subsampled.lm[lm_parietal,,which(completeness$parietal==1)],
                            path="./ply",
                            prefix=NULL,
                            fileext=".ply",
                            ray=TRUE,
                            inflate=4,
                            tol=3,
                            #relax.patch=FALSE,
                            relax.patch=TRUE,
                            keep.fix=1:length(lm_parietal[lm_parietal<103]),
                            #rhotol=NULL,
                            silent=FALSE, mc.cores=3)

#checkLM(patched_parietal,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_parietal,alpha=1,begin=1)


######################
lm_exocc <- subsampled.curve[c(premax_dorsal, nasal, exocc)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_exocc <- module_defs %>% filter(., description=="patch_exocc") %>% pull(., lm)-102

atlas_exocc<- createAtlas(mesh = template.mesh, 
                          landmarks = as.matrix(template.lm[lm_exocc,]),
                          patch = as.matrix(patch[patch_exocc,]),
                          corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                  curvein.nos = c(premax_dorsal, nasal, exocc), 
                                                  rec1 = rec1,
                                                  n.fixed = length(lm_exocc[lm_exocc<103])), 
                          patchCurves = NULL,
                          keep.fix=1:length(lm_exocc[lm_exocc<103]))
patched_exocc <- placePatch(atlas_exocc,
                            subsampled.lm[lm_exocc,,which(completeness$supraocc==1)],
                            path="./ply",
                            prefix=NULL,
                            fileext=".ply",
                            ray=TRUE,
                            inflate=4,
                            tol=3,
                            #relax.patch=FALSE,
                            relax.patch=TRUE,
                            keep.fix=1:length(lm_exocc[lm_exocc<103]),
                            #rhotol=NULL,
                            silent=FALSE, mc.cores=3)

#checkLM(patched_exocc,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_exocc,alpha=1,begin=26)


lm_exocc2 <- subsampled.curve[c(exocc, parietal)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_exocc2 <- module_defs %>% filter(., description=="patch_exocc") %>% pull(., lm)-102

atlas_exocc2<- createAtlas(mesh = template.mesh, 
                          landmarks = as.matrix(template.lm[lm_exocc2,]),
                          patch = as.matrix(patch[patch_exocc2,]),
                          corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                  curvein.nos = c(exocc, parietal), 
                                                  rec1 = rec1,
                                                  n.fixed = length(lm_exocc2[lm_exocc2<103])), 
                          patchCurves = NULL,
                          keep.fix=1:length(lm_exocc2[lm_exocc2<103]))
patched_exocc2 <- placePatch(atlas_exocc2,
                            subsampled.lm[lm_exocc2,,which(completeness$supraocc==1)],
                            path="./ply",
                            prefix=NULL,
                            fileext=".ply",
                            ray=TRUE,
                            inflate=4,
                            tol=3,
                            #relax.patch=FALSE,
                            relax.patch=TRUE,
                            keep.fix=1:length(lm_exocc2[lm_exocc2<103]),
                            #rhotol=NULL,
                            silent=FALSE, mc.cores=3)

#checkLM(patched_exocc2,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_exocc2,alpha=1,begin=31)

######################
lm_occ.condyle <- subsampled.curve[c(occ.condyle)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_occ.condyle <- module_defs %>% filter(., description=="patch_occ_condyle") %>% pull(., lm)-102

atlas_occ.condyle<- createAtlas(mesh = template.mesh, 
                             landmarks = as.matrix(template.lm[lm_occ.condyle,]),
                             patch = as.matrix(patch[patch_occ.condyle,]),
                             corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                     curvein.nos = c(occ.condyle), 
                                                     rec1 = rec1,
                                                     n.fixed = length(lm_occ.condyle[lm_occ.condyle<103])), 
                             patchCurves = NULL,
                             keep.fix=1:length(lm_occ.condyle[lm_occ.condyle<103]))
patched_occ.condyle <- placePatch(atlas_occ.condyle,
                               subsampled.lm[lm_occ.condyle,,which(completeness$occ_condyle==1)],
                               path="./ply",
                               prefix=NULL,
                               fileext=".ply",
                               ray=TRUE,
                               inflate=2,
                               tol=5,
                               #relax.patch=FALSE,
                               relax.patch=TRUE,
                               keep.fix=1:length(lm_occ.condyle[lm_occ.condyle<103]),
                               #rhotol=NULL,
                               silent=FALSE, mc.cores=3)

######################
lm_basiocc <- subsampled.curve[c(basiocc)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_basiocc <- module_defs %>% filter(., description=="patch_basiocc") %>% pull(., lm)-102

atlas_basiocc<- createAtlas(mesh = template.mesh, 
                             landmarks = as.matrix(template.lm[lm_basiocc,]),
                             patch = as.matrix(patch[patch_basiocc,]),
                             corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                     curvein.nos = c(basiocc), 
                                                     rec1 = rec1,
                                                     n.fixed = length(lm_basiocc[lm_basiocc<103])), 
                             patchCurves = NULL,
                             keep.fix=1:length(lm_basiocc[lm_basiocc<103]))
patched_basiocc <- placePatch(atlas_basiocc,
                               subsampled.lm[lm_basiocc,,which(completeness$basiocc==1)],
                               path="./ply",
                               prefix=NULL,
                               fileext=".ply",
                               ray=TRUE,
                               inflate=2,
                               tol=5,
                               #relax.patch=FALSE,
                               relax.patch=TRUE,
                               keep.fix=1:length(lm_basiocc[lm_basiocc<103]),
                               #rhotol=NULL,
                               silent=FALSE, mc.cores=3)

######################
lm_pterygoid <- subsampled.curve[c(nasal,pterygoid)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_pterygoid <- module_defs %>% filter(., description=="patch_pterygoid") %>% pull(., lm)-102

atlas_pterygoid<- createAtlas(mesh = template.mesh, 
                             landmarks = as.matrix(template.lm[lm_pterygoid,]),
                             patch = as.matrix(patch[patch_pterygoid,]),
                             corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                     curvein.nos = c(nasal,pterygoid), 
                                                     rec1 = rec1,
                                                     n.fixed = length(lm_pterygoid[lm_pterygoid<103])), 
                             patchCurves = NULL,
                             keep.fix=1:length(lm_pterygoid[lm_pterygoid<103]))
patched_pterygoid <- placePatch(atlas_pterygoid,
                               subsampled.lm[lm_pterygoid,,which(completeness$pterygoid==1)],
                               path="./ply",
                               prefix=NULL,
                               fileext=".ply",
                               ray=TRUE,
                               inflate=5,
                               tol=4,
                               #relax.patch=FALSE,
                               relax.patch=TRUE,
                               keep.fix=1:length(lm_pterygoid[lm_pterygoid<103]),
                               #rhotol=NULL,
                               silent=FALSE, mc.cores=3)

patched_pterygoid_small <- placePatch(atlas_pterygoid,
                                subsampled.lm[lm_pterygoid,,which(completeness$pterygoid==1)],
                                path="./ply",
                                prefix=NULL,
                                fileext=".ply",
                                ray=TRUE,
                                inflate=2,
                                tol=.8,
                                #relax.patch=FALSE,
                                relax.patch=TRUE,
                                keep.fix=1:length(lm_pterygoid[lm_pterygoid<103]),
                                #rhotol=NULL,
                                silent=FALSE, mc.cores=3)

#checkLM(patched_pterygoid,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_pterygoid,alpha=1,begin=1)


######################
lm_palatine <- subsampled.curve[c(max_vent, nasal, palatine)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_palatine <- module_defs %>% filter(., description=="patch_palatine") %>% pull(., lm)-102

atlas_palatine<- createAtlas(mesh = template.mesh, 
                             landmarks = as.matrix(template.lm[lm_palatine,]),
                             patch = as.matrix(patch[patch_palatine,]),
                             corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                     curvein.nos = c(max_vent, nasal, palatine), 
                                                     rec1 = rec1,
                                                     n.fixed = length(lm_palatine[lm_palatine<103])), 
                             patchCurves = NULL,
                             keep.fix=1:length(lm_palatine[lm_palatine<103]))
patched_palatine <- placePatch(atlas_palatine,
                               subsampled.lm[lm_palatine,,which(completeness$palatine==1)],
                               path="./ply",
                               prefix=NULL,
                               fileext=".ply",
                               ray=TRUE,
                               inflate=4,
                               tol=3,
                               #relax.patch=FALSE,
                               relax.patch=TRUE,
                               keep.fix=1:length(lm_palatine[lm_palatine<103]),
                               #rhotol=NULL,
                               silent=FALSE, mc.cores=3)

#checkLM(patched_palatine,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_palatine,alpha=1,begin=1)



######################
lm_jugal <- subsampled.curve[c(jugal)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_jugal <- module_defs %>% filter(., description=="patch_jugal") %>% pull(., lm)-102

atlas_jugal<- createAtlas(mesh = template.mesh, 
                             landmarks = as.matrix(template.lm[lm_jugal,]),
                             patch = as.matrix(patch[patch_jugal,]),
                             corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                     curvein.nos = c(jugal), 
                                                     rec1 = rec1,
                                                     n.fixed = length(lm_jugal[lm_jugal<103])), 
                             patchCurves = NULL,
                             keep.fix=1:length(lm_jugal[lm_jugal<103]))
patched_jugal <- placePatch(atlas_jugal,
                               subsampled.lm[lm_jugal,,which(completeness$jugal==1)],
                               path="./ply",
                               prefix=NULL,
                               fileext=".ply",
                               ray=TRUE,
                               inflate=2,
                               tol=5,
                               #relax.patch=FALSE,
                               relax.patch=TRUE,
                               keep.fix=1:length(lm_jugal[lm_jugal<103]),
                               #rhotol=NULL,
                               silent=FALSE, mc.cores=3)

#checkLM(patched_jugal,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_jugal,alpha=1,begin=39)


######################
lm_squamosal <- subsampled.curve[c(frontal, jugal, squamosal)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_squamosal <- module_defs %>% filter(., description=="patch_squamosal") %>% pull(., lm)-102

atlas_squamosal<- createAtlas(mesh = template.mesh, 
                          landmarks = as.matrix(template.lm[lm_squamosal,]),
                          patch = as.matrix(patch[patch_squamosal,]),
                          corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                  curvein.nos = c(frontal, jugal, squamosal), 
                                                  rec1 = rec1,
                                                  n.fixed = length(lm_squamosal[lm_squamosal<103])), 
                          patchCurves = NULL,
                          keep.fix=1:length(lm_squamosal[lm_squamosal<103]))
patched_squamosal <- placePatch(atlas_squamosal,
                            subsampled.lm[lm_squamosal,,which(completeness$squamosal==1)],
                            path="./ply",
                            prefix=NULL,
                            fileext=".ply",
                            ray=TRUE,
                            inflate=5,
                            tol=4,
                            #relax.patch=FALSE,
                            relax.patch=TRUE,
                            keep.fix=1:length(lm_squamosal[lm_squamosal<103]),
                            #rhotol=NULL,
                            silent=FALSE, mc.cores=3)

#checkLM(patched_squamosal,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_squamosal,alpha=1,begin=40)


lm_squamosal2 <- subsampled.curve[c(frontal, squamosal)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_squamosal2 <- module_defs %>% filter(., description=="patch_squamosal") %>% pull(., lm)-102

atlas_squamosal2 <- createAtlas(mesh = template.mesh, 
                              landmarks = as.matrix(template.lm[lm_squamosal2,]),
                              patch = as.matrix(patch[patch_squamosal2,]),
                              corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                      curvein.nos = c(frontal, squamosal), 
                                                      rec1 = rec1,
                                                      n.fixed = length(lm_squamosal2[lm_squamosal2<103])), 
                              patchCurves = NULL,
                              keep.fix=1:length(lm_squamosal2[lm_squamosal2<103]))
patched_squamosal2 <- placePatch(atlas_squamosal2,
                                subsampled.lm[lm_squamosal2,,which(completeness$squamosal==1)],
                                path="./ply",
                                prefix=NULL,
                                fileext=".ply",
                                ray=TRUE,
                                inflate=2,
                                tol=.8,
                                #relax.patch=FALSE,
                                relax.patch=TRUE,
                                keep.fix=1:length(lm_squamosal2[lm_squamosal2<103]),
                                #rhotol=NULL,
                                silent=FALSE, mc.cores=3)

#checkLM(patched_squamosal2,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_squamosal2,alpha=1,begin=27)

######################
lm_lacrimal <- subsampled.curve[c(lacrimal)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_lacrimal <- module_defs %>% filter(., description=="patch_lac") %>% pull(., lm)-102

atlas_lacrimal<- createAtlas(mesh = template.mesh, 
                          landmarks = as.matrix(template.lm[lm_lacrimal,]),
                          patch = as.matrix(patch[patch_lacrimal,]),
                          corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                  curvein.nos = c(lacrimal), 
                                                  rec1 = rec1,
                                                  n.fixed = length(lm_lacrimal[lm_lacrimal<103])), 
                          patchCurves = NULL,
                          keep.fix=1:length(lm_lacrimal[lm_lacrimal<103]))
patched_lacrimal <- placePatch(atlas_lacrimal,
                            subsampled.lm[lm_lacrimal,,which(completeness$lacrimal==1)],
                            path="./ply",
                            prefix=NULL,
                            fileext=".ply",
                            ray=TRUE,
                            inflate=2,
                            tol=5,
                            #relax.patch=FALSE,
                            relax.patch=TRUE,
                            keep.fix=1:length(lm_lacrimal[lm_lacrimal<103]),
                            #rhotol=NULL,
                            silent=FALSE, mc.cores=3)


checkLM(patched_lacrimal,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_lacrimal,alpha=1,begin=1)

######################\
lm_lacrimal2 <- subsampled.curve[c(lacrimal,palatine)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_lacrimal2 <- module_defs %>% filter(., description=="patch_lac") %>% pull(., lm)-102

atlas_lacrimal2<- createAtlas(mesh = template.mesh, 
                             landmarks = as.matrix(template.lm[lm_lacrimal2,]),
                             patch = as.matrix(patch[patch_lacrimal2,]),
                             corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                     curvein.nos = c(lacrimal,palatine), 
                                                     rec1 = rec1,
                                                     n.fixed = length(lm_lacrimal2[lm_lacrimal2<103])), 
                             patchCurves = NULL,
                             keep.fix=1:length(lm_lacrimal2[lm_lacrimal2<103]))
patched_lacrimal2 <- placePatch(atlas_lacrimal2,
                               subsampled.lm[lm_lacrimal2,,which(completeness$lacrimal==1)],
                               path="./ply",
                               prefix=NULL,
                               fileext=".ply",
                               ray=TRUE,
                               inflate=2,
                               tol=.8,
                               #relax.patch=FALSE,
                               relax.patch=TRUE,
                               keep.fix=1:length(lm_lacrimal2[lm_lacrimal2<103]),
                               #rhotol=NULL,
                               silent=FALSE, mc.cores=3)

#checkLM(patched_lacrimal2,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_lacrimal,alpha=1,begin=15)


######################
lm_ecto.pter <- subsampled.curve[c(nasal, jugal, ecto.pter)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_ecto.pter <- module_defs %>% filter(., description=="patch_ectoptery") %>% pull(., lm)-102

atlas_ecto.pter<- createAtlas(mesh = template.mesh, 
                          landmarks = as.matrix(template.lm[lm_ecto.pter,]),
                          patch = as.matrix(patch[patch_ecto.pter,]),
                          corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                  curvein.nos = c(nasal, jugal, ecto.pter), 
                                                  rec1 = rec1,
                                                  n.fixed = length(lm_ecto.pter[lm_ecto.pter<103])), 
                          patchCurves = NULL,
                          keep.fix=1:length(lm_ecto.pter[lm_ecto.pter<103]))
patched_ecto.pter <- placePatch(atlas_ecto.pter,
                            subsampled.lm[lm_ecto.pter,,which(completeness$ectopt==1)],
                            path="./ply",
                            prefix=NULL,
                            fileext=".ply",
                            ray=TRUE,
                            inflate=2,
                            tol=.8,
                            #relax.patch=FALSE,
                            relax.patch=TRUE,
                            keep.fix=1:length(lm_ecto.pter[lm_ecto.pter<103]),
                            #rhotol=NULL,
                            silent=FALSE, mc.cores=3)

#checkLM(patched_ecto.pter,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_ecto.pter,alpha=1,begin=1)


######################
lm_ptery.joint <- subsampled.curve[c(nasal, ptery.joint)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_ptery.joint <- module_defs %>% filter(., description=="patch_ptery_joint") %>% pull(., lm)-102

atlas_ptery.joint<- createAtlas(mesh = template.mesh, 
                             landmarks = as.matrix(template.lm[lm_ptery.joint,]),
                             patch = as.matrix(patch[patch_ptery.joint,]),
                             corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                     curvein.nos = c(nasal, ptery.joint), 
                                                     rec1 = rec1,
                                                     n.fixed = length(lm_ptery.joint[lm_ptery.joint<103])), 
                             patchCurves = NULL,
                             keep.fix=1:length(lm_ptery.joint[lm_ptery.joint<103]))
patched_ptery.joint <- placePatch(atlas_ptery.joint,
                               subsampled.lm[lm_ptery.joint,,which(completeness$ptery_flange==1)],
                               path="./ply",
                               prefix=NULL,
                               fileext=".ply",
                               ray=TRUE,
                               inflate=2,
                               tol=.8,
                               relax.patch=FALSE,
                               #relax.patch=TRUE,
                               keep.fix=1:length(lm_ptery.joint[lm_ptery.joint<103]),
                               #rhotol=NULL,
                               silent=FALSE, mc.cores=3)

#checkLM(patched_ptery.joint,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_ptery.joint,alpha=1,begin=1)


######################
lm_jawjoint <- subsampled.curve[c(nasal, jawjoint)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_jawjoint <- module_defs %>% filter(., description=="patch_jawjoint") %>% pull(., lm)-102

atlas_jawjoint<- createAtlas(mesh = template.mesh, 
                             landmarks = as.matrix(template.lm[lm_jawjoint,]),
                             patch = as.matrix(patch[patch_jawjoint,]),
                             corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                     curvein.nos = c(nasal, jawjoint), 
                                                     rec1 = rec1,
                                                     n.fixed = length(lm_jawjoint[lm_jawjoint<103])), 
                             patchCurves = NULL,
                             keep.fix=1:length(lm_jawjoint[lm_jawjoint<103]))
patched_jawjoint <- placePatch(atlas_jawjoint,
                               subsampled.lm[lm_jawjoint,,which(completeness$jawjoint==1)],
                               path="./ply",
                               prefix=NULL,
                               fileext=".ply",
                               ray=TRUE,
                               inflate=5,
                               tol=4,
                               #relax.patch=FALSE,
                               relax.patch=TRUE,
                               keep.fix=1:length(lm_jawjoint[lm_jawjoint<103]),
                               #rhotol=NULL,
                               silent=FALSE, mc.cores=3)

#checkLM(patched_jawjoint,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_jawjoint,alpha=1,begin=1)


######################
lm_postorb <- subsampled.curve[c(nasal, postorb)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_postorb <- module_defs %>% filter(., description=="patch_postorb") %>% pull(., lm)-102

atlas_postorb<- createAtlas(mesh = template.mesh, 
                             landmarks = as.matrix(template.lm[lm_postorb,]),
                             patch = as.matrix(patch[patch_postorb,]),
                             corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                     curvein.nos = c(nasal, postorb), 
                                                     rec1 = rec1,
                                                     n.fixed = length(lm_postorb[lm_postorb<103])), 
                             patchCurves = NULL,
                             keep.fix=1:length(lm_postorb[lm_postorb<103]))
patched_postorb <- placePatch(atlas_postorb,
                               subsampled.lm[lm_postorb,,which(completeness$postorb==1)],
                               path="./ply",
                               prefix=NULL,
                               fileext=".ply",
                               ray=TRUE,
                               inflate=3,
                               tol=2,
                               #relax.patch=FALSE,
                               relax.patch=TRUE,
                               keep.fix=1:length(lm_postorb[lm_postorb<103]),
                               #rhotol=NULL,
                               silent=FALSE, mc.cores=3)


lm_postorb2 <- subsampled.curve[c(nasal, squamosal, parietal, palatine, postorb)]%>%unlist(.)%>%unique(.)%>%sort(.)

patch_postorb2 <- module_defs %>% filter(., description=="patch_postorb") %>% pull(., lm)-102

atlas_postorb2<- createAtlas(mesh = template.mesh, 
                            landmarks = as.matrix(template.lm[lm_postorb2,]),
                            patch = as.matrix(patch[patch_postorb2,]),
                            corrCurves=ReCorrCurves(curve.in = subsampled.curve.in,
                                                    curvein.nos = c(nasal, squamosal, parietal, palatine, postorb), 
                                                    rec1 = rec1,
                                                    n.fixed = length(lm_postorb2[lm_postorb2<103])), 
                            patchCurves = NULL,
                            keep.fix=1:length(lm_postorb2[lm_postorb2<103]))
patched_postorb2 <- placePatch(atlas_postorb2,
                              subsampled.lm[lm_postorb2,,which(completeness$postorb==1)],
                              path="./ply",
                              prefix=NULL,
                              fileext=".ply",
                              ray=TRUE,
                              inflate=3,
                              tol=2,
                              #relax.patch=FALSE,
                              relax.patch=TRUE,
                              keep.fix=1:length(lm_postorb2[lm_postorb2<103]),
                              #rhotol=NULL,
                              silent=FALSE, mc.cores=3)
#checkLM(patched_postorb2,path="./ply/",suffix=".ply",pt.size=2,render="s",atlas=atlas_postorb2,alpha=1,begin=27)


 
 
 ################################
 ##combine patches
 ################################
 patched.combined <- array(data = NA, dim = c(total.lm.num,3,ntaxa)) #empty array with the correct number of specimens and landmarks
 dimnames(patched.combined)[3]<-dimnames(subsampled.lm)[3] #assign taxon names to empty array
 patched.combined[c(1:dim(subsampled.lm)[1]),,]<-subsampled.lm #put the landmarks and semilandmarks in the array
 
 num.regular<-last(slidings.sub) #number of non patch points
 
 patched.combined[patch_frontal+num.regular,,which(completeness$frontal==1)]<-patched_frontal[-c(1:length(lm_frontal)),,]
 patched.combined[patch_parietal+num.regular,,which(completeness$parietal==1)]<-patched_parietal[-c(1:length(lm_parietal)),,]
 patched.combined[patch_jawjoint+num.regular,,which(completeness$jawjoint==1)]<-patched_jawjoint[-c(1:length(lm_jawjoint)),,]
 patched.combined[patch_lacrimal+num.regular,,which(completeness$lacrimal==1)]<-patched_lacrimal[-c(1:length(lm_lacrimal)),,]
 patched.combined[patch_max_d2+num.regular,,which(completeness$max_d==1)]<-patched_max_d2[-c(1:length(lm_max_d2)),,]
 patched.combined[patch_premax_d+num.regular,,which(completeness$premax_d==1)]<-patched_premax_d[-c(1:length(lm_premax_d)),,]
 patched.combined[patch_nasal3+num.regular,,which(completeness$nasal==1)]<-patched_nasal3[-c(1:length(lm_nasal3)),,]
 patched.combined[patch_basiocc+num.regular,,which(completeness$basiocc==1)]<-patched_basiocc[-c(1:length(lm_basiocc)),,]
 
 patched.combined[patch_exocc+num.regular,,which(completeness$supraocc==1)]<-patched_exocc[-c(1:length(lm_exocc)),,]
 patched.combined[patch_occ.condyle+num.regular,,which(completeness$occ_condyle==1)]<-patched_occ.condyle[-c(1:length(lm_occ.condyle)),,]

  patched.combined[patch_palatine+num.regular,,which(completeness$palatine==1)]<-patched_palatine[-c(1:length(lm_palatine)),,]
 patched.combined[patch_pterygoid+num.regular,,which(completeness$pterygoid==1)]<-patched_pterygoid[-c(1:length(lm_pterygoid)),,]
 patched.combined[patch_premax_vent+num.regular,,which(completeness$premax_vent==1)]<-patched_premax_vent[-c(1:length(lm_premax_vent)),,]
 patched.combined[patch_max_vent+num.regular,,which(completeness$max_vent==1)]<-patched_max_vent[-c(1:length(lm_max_vent)),,]
 patched.combined[patch_squamosal+num.regular,,which(completeness$squamosal==1)]<-patched_squamosal[-c(1:length(lm_squamosal)),,]
 patched.combined[patch_postorb2+num.regular,,which(completeness$postorb==1)]<-patched_postorb2[-c(1:length(lm_postorb2)),,]
 patched.combined[patch_jugal+num.regular,,which(completeness$jugal==1)]<-patched_jugal[-c(1:length(lm_jugal)),,]
 patched.combined[patch_ptery.joint+num.regular,,which(completeness$ptery_flange==1)]<-patched_ptery.joint[-c(1:length(lm_ptery.joint)),,]
 patched.combined[patch_ecto.pter+num.regular,,which(completeness$ectopt==1)]<-patched_ecto.pter[-c(1:length(lm_ecto.pter)),,]
 
 patched.combined[patch_max_d+num.regular,,which(dimnames(patched.combined)[[3]]=="Marilasuchus_cleaned retrodeformed and mirrored")]<-patched_max_d2[-c(1:length(lm_max_d2)),,which(dimnames(patched_max_d2)[[3]]=="Marilasuchus_cleaned retrodeformed and mirrored")]
 patched.combined[patch_squamosal2+num.regular,,which(dimnames(patched.combined)[[3]]=="Geosaurus(Cricosaurus)_araucanensis_MLP_72_IV_7_1" )]<-patched_squamosal2[-c(1:length(lm_squamosal2)),,which(dimnames(patched_squamosal2)[[3]]=="Geosaurus(Cricosaurus)_araucanensis_MLP_72_IV_7_1")]
 patched.combined[patch_pterygoid+num.regular,,which(dimnames(patched.combined)[[3]]=="Araripesuchus gomesii AMNH24450 skull 11022016")]<-patched_pterygoid_small[-c(1:length(lm_pterygoid)),,which(dimnames(patched_pterygoid_small)[[3]]=="Araripesuchus gomesii AMNH24450 skull 11022016"  )]
 patched.combined[patch_premax_d2+num.regular,,which(dimnames(patched.combined)[[3]]=="Sarcosuchus_imperator_TMP_2009_003_0005(cast)-lowres")]<-patched_premax_d2[-c(1:length(lm_premax_d2)),,which(dimnames(patched_premax_d2)[[3]]=="Sarcosuchus_imperator_TMP_2009_003_0005(cast)-lowres")]
 patched.combined[patch_nasal2+num.regular,,which(dimnames(patched.combined)[[3]]=="Araripesuchus_Skull_repaired")]<-patched_nasal2[-c(1:length(lm_nasal2)),,which(dimnames(patched_nasal2)[[3]]=="Araripesuchus_Skull_repaired")]
 patched.combined[patch_nasal2+num.regular,,which(dimnames(patched.combined)[[3]]=="Marilasuchus_cleaned retrodeformed and mirrored")]<-patched_nasal2[-c(1:length(lm_nasal2)),,which(dimnames(patched_nasal2)[[3]]=="Marilasuchus_cleaned retrodeformed and mirrored")]
 patched.combined[patch_nasal2+num.regular,,which(dimnames(patched.combined)[[3]]=="Araripesuchus gomesii AMNH24450 skull 11022016" )]<-patched_nasal2[-c(1:length(lm_nasal2)),,which(dimnames(patched_nasal2)[[3]]=="Araripesuchus gomesii AMNH24450 skull 11022016" )]

 patched.combined[patch_nasal+num.regular,,which(dimnames(patched.combined)[[3]]=="Gavialis_gangenticus_uf-herp-118998_M20938-39872" )]<-patched_nasal[-c(1:length(lm_nasal)),,which(dimnames(patched_nasal)[[3]]=="Gavialis_gangenticus_uf-herp-118998_M20938-39872" )]
 patched.combined[patch_nasal2+num.regular,,which(dimnames(patched.combined)[[3]]=="Geosaurus(Cricosaurus)_araucanensis_MLP_72_IV_7_1")]<-patched_nasal2[-c(1:length(lm_nasal2)),,which(dimnames(patched_nasal2)[[3]]=="Geosaurus(Cricosaurus)_araucanensis_MLP_72_IV_7_1")]
 patched.combined[patch_nasal2+num.regular,,which(dimnames(patched.combined)[[3]]=="Palagosaurus_cranium")]<-patched_nasal2[-c(1:length(lm_nasal2)),,which(dimnames(patched_nasal2)[[3]]=="Palagosaurus_cranium" )]
 
 patched.combined[patch_postorb+num.regular,,which(dimnames(patched.combined)[[3]]=="Geosaurus(Cricosaurus)_araucanensis_MLP_72_IV_7_1"  )]<-patched_postorb[-c(1:length(lm_postorb)),,which(dimnames(patched_postorb)[[3]]=="Geosaurus(Cricosaurus)_araucanensis_MLP_72_IV_7_1"  )]
 patched.combined[patch_premax_d2+num.regular,,which(dimnames(patched.combined)[[3]]=="Alligator_mississippiensis_AMNH_R_40582_craniumMIRRORED")]<-patched_premax_d2[-c(1:length(lm_premax_d2)),,which(dimnames(patched_premax_d2)[[3]]=="Alligator_mississippiensis_AMNH_R_40582_craniumMIRRORED")]
 patched.combined[patch_premax_d3+num.regular,,which(dimnames(patched.combined)[[3]]=="Crocodylus_acutus_AMNH_7857_cranium")]<-patched_premax_d3[-c(1:length(lm_premax_d3)),,which(dimnames(patched_premax_d3)[[3]]=="Crocodylus_acutus_AMNH_7857_cranium")]
 patched.combined[patch_premax_d3+num.regular,,which(dimnames(patched.combined)[[3]]=="Crocodylus_johnstoni_TMM_M-6807")]<-patched_premax_d3[-c(1:length(lm_premax_d3)),,which(dimnames(patched_premax_d3)[[3]]=="Crocodylus_johnstoni_TMM_M-6807")]
 patched.combined[patch_premax_d3+num.regular,,which(dimnames(patched.combined)[[3]]=="Crocodylus_mindorensis_FMNH-11136-7_Cranium2")]<-patched_premax_d3[-c(1:length(lm_premax_d3)),,which(dimnames(patched_premax_d3)[[3]]=="Crocodylus_mindorensis_FMNH-11136-7_Cranium2")]
 patched.combined[patch_premax_d3+num.regular,,which(dimnames(patched.combined)[[3]]=="Crocodylus_novaeguineae_AMNH_64425_cranium")]<-patched_premax_d3[-c(1:length(lm_premax_d3)),,which(dimnames(patched_premax_d3)[[3]]=="Crocodylus_novaeguineae_AMNH_64425_cranium")]
 patched.combined[patch_premax_d3+num.regular,,which(dimnames(patched.combined)[[3]]=="Crocodylus_raninus_AMNH_29294_craniumMIRRORED")]<-patched_premax_d3[-c(1:length(lm_premax_d3)),,which(dimnames(patched_premax_d3)[[3]]=="Crocodylus_raninus_AMNH_29294_craniumMIRRORED")]
 patched.combined[patch_premax_d3+num.regular,,which(dimnames(patched.combined)[[3]]=="Crocodylus_siamensis cranium BMNH 1920.1.20.1626")]<-patched_premax_d3[-c(1:length(lm_premax_d3)),,which(dimnames(patched_premax_d3)[[3]]=="Crocodylus_siamensis cranium BMNH 1920.1.20.1626")]
 patched.combined[patch_premax_d3+num.regular,,which(dimnames(patched.combined)[[3]]=="Gavialis_gangenticus_uf-herp-118998_M20938-39872")]<-patched_premax_d3[-c(1:length(lm_premax_d3)),,which(dimnames(patched_premax_d3)[[3]]=="Gavialis_gangenticus_uf-herp-118998_M20938-39872")]
 patched.combined[patch_premax_d3+num.regular,,which(dimnames(patched.combined)[[3]]=="Mecistops_cataphractus_AMNH_10075_cranium")]<-patched_premax_d3[-c(1:length(lm_premax_d3)),,which(dimnames(patched_premax_d3)[[3]]=="Mecistops_cataphractus_AMNH_10075_cranium")]
 patched.combined[patch_premax_d3+num.regular,,which(dimnames(patched.combined)[[3]]=="Gavialis_gangenticus_uf-herp-118998_M20938-39872")]<-patched_premax_d3[-c(1:length(lm_premax_d3)),,which(dimnames(patched_premax_d3)[[3]]=="Gavialis_gangenticus_uf-herp-118998_M20938-39872")]
 patched.combined[patch_premax_d3+num.regular,,which(dimnames(patched.combined)[[3]]=="Melanosuchus_niger_AMNH_58132_cranium")]<-patched_premax_d3[-c(1:length(lm_premax_d3)),,which(dimnames(patched_premax_d3)[[3]]=="Melanosuchus_niger_AMNH_58132_cranium")]
 patched.combined[patch_premax_d3+num.regular,,which(dimnames(patched.combined)[[3]]=="Palagosaurus_cranium")]<-patched_premax_d3[-c(1:length(lm_premax_d3)),,which(dimnames(patched_premax_d3)[[3]]=="Palagosaurus_cranium")]
 patched.combined[patch_premax_d3+num.regular,,which(dimnames(patched.combined)[[3]]=="Paleosuchus_trigonatus_AMNH_6691_cranium")]<-patched_premax_d3[-c(1:length(lm_premax_d3)),,which(dimnames(patched_premax_d3)[[3]]=="Paleosuchus_trigonatus_AMNH_6691_cranium")]
 patched.combined[patch_premax_d3+num.regular,,which(dimnames(patched.combined)[[3]]=="Tomistoma_schlegelii_MCZ_12459_cranium")]<-patched_premax_d3[-c(1:length(lm_premax_d3)),,which(dimnames(patched_premax_d3)[[3]]=="Tomistoma_schlegelii_MCZ_12459_cranium")]
 patched.combined[patch_lacrimal2+num.regular,,which(dimnames(patched.combined)[[3]]=="Crocodylus_affinus_UWBM_88113")]<-patched_lacrimal2[-c(1:length(lm_lacrimal2)),,which(dimnames(patched_lacrimal2)[[3]]=="Crocodylus_affinus_UWBM_88113")]
 patched.combined[patch_premax_d3+num.regular,,which(dimnames(patched.combined)[[3]]=="Crocodilus_acutus_cranium_MNHN_A_5310")]<-patched_premax_d3[-c(1:length(lm_premax_d3)),,which(dimnames(patched_premax_d3)[[3]]=="Crocodilus_acutus_cranium_MNHN_A_5310")]
 #patched.combined[patch_max_d3+num.regular,,which(dimnames(patched.combined)[[3]]=="Palagosaurus_cranium")]<-patched_max_d3[-c(1:length(lm_max_d3)),,which(dimnames(patched_max_d3)[[3]]=="Palagosaurus_cranium")]
 #patched.combined[patch_max_d3+num.regular,,which(dimnames(patched.combined)[[3]]=="Pholidosaurus_purbeckensis_cranium" )]<-patched_max_d3[-c(1:length(lm_max_d3)),,which(dimnames(patched_max_d3)[[3]]=="Pholidosaurus_purbeckensis_cranium" )]
 patched.combined[patch_exocc2+num.regular,,which(dimnames(patched.combined)[[3]]=="Palagosaurus_cranium" )]<-patched_exocc2[-c(1:length(lm_exocc2)),,which(dimnames(patched_exocc2)[[3]]=="Palagosaurus_cranium" )]
 patched.combined[patch_max_d+num.regular,,which(dimnames(patched.combined)[[3]]=="Tomistoma_schlegelii_MCZ_12459_cranium"   )]<-patched_max_d[-c(1:length(lm_max_d)),,which(dimnames(patched_max_d)[[3]]=="Tomistoma_schlegelii_MCZ_12459_cranium"   )]
 
 
 color.palette <- c( "#6a3d9a","dimgrey","#fb9a99",  "gold", "#009E73",  "#D55E00", "#CC79A7", "cyan2",  "#e31a1c", "#0072B2", "#b2df8a", "#E69F00",  "whitesmoke" ,  "deeppink",   "#a6cee3",   "#F0E442","blue","red","brown", "black")
 modulecolors1<-as.factor(module_defs$full_model)
 levels(modulecolors1)<-color.palette
 
 
 source('~/Google Drive/utlity_scripts/check.modules.R')
 check.modules(patched.combined,path="./ply/",suffix=".ply",pt.size=3,render="s",module.colors = modulecolors1,begin=30)
   
#save(patched.combined, file = "~/Google Drive/NHM/crocs/data/patched.crocs_sept_23_2019.R")
load("~/Google Drive/NHM/crocs/data/patched.crocs_sept_23_2019.R")
 
 
 
 module_defs_nosph<-module_defs[-which(module_defs$full_model==19),]
 
 
 #import meshes
 m1<-vcgImport("~/Google Drive/NHM/crocs/meshes/Alligator_mississippiensis_AMNH_R_40582_craniumforfigures.ply")
 m2<-vcgImport("~/Dropbox/Crocs/pts3/ply/Alligator_sinensis_AMNH_R_23900_cranium.ply")

 

 
 
  patched.combined[which(patched.combined==9999)]
 
  patched.combined2<-patched.combined
  #checkLM(patched.combined,atlas=atlas.subsampl,begin=1,pt.size = 2,path="./ply/",suffix=".ply",render="s")
  
  
  fixed.nums<-fixed
  curve.nums<-subsampled.curve.in %>% unlist(.) %>% unique(.) %>% sort(.)
  patch.nums<-c((last(curve.nums)+1):dim(patched.combined2)[[1]])
  
  
  dat.array<-patched.combined
  
  n <- dim(dat.array)[3]
  k <- dim(dat.array)[1]
  m <- dim(dat.array)[2]
  
  ignore = module_defs$row.id[which(module_defs$full_model ==19)]
  li <- length(ignore)
  lm.old <- c(1:k)[-ignore]
  mat.ptr <- matrix(c(1:(k-li),lm.old),k-li,2)
  ptr <- function(xo)	### define pointer function for indexing
  {
    if (length(which(ignore %in% xo))!= 0)
      xo <- xo[-which(xo %in% ignore)]
    for (i in 1:(k-li))
      xo[which(xo==mat.ptr[i,2])] <- mat.ptr[i,1]
    return(xo)
  }
  ### update outline indices
  outlines1 <- lapply(subsampled.curve,ptr)
  outlines1 <- outlines1[-c(33:35)]
  outlines.in1 <- lapply(subsampled.curve.in,ptr)
  
  
  surp1 <- ptr(patch.nums)
  
  SMvector1 <- ptr(c(curve.nums,patch.nums))
  
  
  #remove bsphenoid!
  lm_sphenoid <-  subsampled.curve[basisph]%>%unlist(.)%>%unique(.)%>%sort(.)
  lm_no_sphenoid <- c(1:max(slidings.sub))[-lm_sphenoid]
  curves_no_sphenoid <- ReCorrCurves(curve.in = subsampled.curve.in,
                                  
                                      curvein.nos = c(1:length(subsampled.curve.in))[-basisph], 
                                      rec1 = rec1,
                                      n.fixed = length(lm_no_sphenoid[lm_no_sphenoid<103]))
  

  
 
  #get complete specimens:
 good_spec_names<-completeness %>% na_if(., 0) %>% na.omit(.) %>% pull(., Filename)
 good_specs<- patched.combined[,,good_spec_names]
 lm_sphenoid <-  subsampled.curve[basisph]%>%unlist(.)%>%unique(.)%>%sort(.)
 
 ###########
 #temporarilly remove spec"Stangerochampsa_mccabei_TMP_1986_061_0001"  cuz it has messed up premax vent
 #temporarilly remove pelagosaurus typus  cuz it has messed up ectopterygoid
 
#good_specs<-good_specs[,,-c(30,38)]
 
#now remove sphenoid points
 
good_specs2 <- good_specs[-which(module_defs$full_model==19),,]

 good_specs2<-estimate.missing(good_specs2)
 
 modulecolors2<-as.factor(module_defs_nosph$full_model)
 levels(modulecolors2)<-color.palette
 
 check.modules(good_specs2,path="./ply/",suffix=".ply",pt.size=2,render="s",module.colors = modulecolors2,begin=1)
 


 {
   
   slided4.new <- slider3d(good_specs2, SMvector= SMvector1,
                                     outlines = outlines1, surp =surp1,
                                     sur.path = "./ply", sur.name = NULL, 
                                     meshlist = paste("./ply/",dimnames(good_specs2)[[3]],".ply",sep=""), ignore = NULL,
                                     sur.type = "ply", tol = 1e-10, deselect = FALSE, inc.check = FALSE,
                                     recursive = TRUE, iterations = 3, initproc = TRUE,
                                     pairedLM = 0, mc.cores = 4, bending=TRUE,
                                     fixRepro = FALSE,stepsize=0.2)
   dimnames(slided4.new[["dataslide"]])[3]<-dimnames(good_specs2)[3]
   
   #save(slided4.new,file="~/Google Drive/NHM/crocs/data/slid.crocs.sept23.R")
   #load("~/Google Drive/NHM/crocs/data/slid.crocs.sept23.R")
 }
 

x<-gpagen(slided4.new$dataslide)
plotTangentSpace(x$coords)
 
 

#mirror
#bilat.landmarks <- cbind(module_defs$row.id[which(module_defs$l_or_r=="r")],module_defs$row.id[which(module_defs$l_or_r=="l")])
bilat.landmarks <- cbind(which(module_defs_nosph$l_or_r=="r"),which(module_defs_nosph$l_or_r=="l"))

bilat.landmarks <- bilat.landmarks[c(1:41),]
specimens<-slided4.new$dataslide
midline<- which(module_defs_nosph$l_or_r=="m")#a list of which landmarks fall on the midline
#a list of the CORRESPONDING left side landmarks. (the order of right.lm and left.lm must represent bilateral pairs of landmarks)
surp.no.sphenoid <- c((max(unlist(curves_no_sphenoid))+1):dim(good_specs2)[1])

all.curves.present<- c(1:length(curves_no_sphenoid))
right.curves <- all.curves.present[-c(4,13,16,20,24,30,32)]#a list of the curves which are going to be mirrored to the other side (excludes midline curves)
right.curve.list<-unlist(curves_no_sphenoid[right.curves])
rightside<-c(bilat.landmarks[,1],right.curve.list,surp.no.sphenoid)
num.missing<-(length(rightside)-length(bilat.landmarks[,2]))
blanks<-c((dim(specimens)[1]+1):(dim(specimens)[1]+num.missing))         
leftside<-c(bilat.landmarks[,2],blanks)

module_defs_right <- module_defs_nosph %>% filter(., l_or_r != "l")

add_col_or_row = function(x, n = 1, add_col = T, fill = 0)
{
  m1 = matrix(x, ncol = if(add_col) nrow(x) * ncol(x) else nrow(x), byrow = T)
  m2 = matrix(fill, nrow = if(add_col) dim(x)[3] else prod(dim(x)[-1]), 
              ncol = if(add_col) nrow(x) * n else n)
  array(t(cbind(m1, m2)), 
        c(nrow(x) + ((!add_col) * n), ncol(x) + (add_col * n), dim(x)[3]))
}

specimens2<-add_col_or_row(specimens,n=num.missing,add_col=FALSE,fill=NA)
dimnames(specimens2)[3]<-dimnames(specimens)[3]
bilats<-cbind(rightside, leftside)
library(paleomorph)
newarray<-mirrorfill(specimens2,l1=as.integer(midline),l2=bilats)
dimnames(newarray)[3]<-dimnames(specimens)[3]
checkLM(newarray, begin=37,pt.size = 2,path="./ply/",suffix=".ply",render="s")


color_table <- read_csv("~/Google Drive/_UCL/Projects/Archosaur Modularity/Colortable.csv")
modulecolors<-as.factor(module_defs_right$full_model)
levels(modulecolors)<-color_table$croc_colors

open3d();shade3d(m1,color='#E4D1C0')
spheres3d(newarray[-leftside,,1],col=modulecolors,radius=1.5)

#croc.vent<-list(zoom=par3d()$zoom,
 #               userMatrix=par3d()$userMatrix,
  #              windowRect=par3d()$windowRect)

#croc.lateral<-list(zoom=par3d()$zoom,
 #               userMatrix=par3d()$userMatrix,
  #              windowRect=par3d()$windowRect)

#croc.dorsal<-list(zoom=par3d()$zoom,
 #               userMatrix=par3d()$userMatrix,
  #              windowRect=par3d()$windowRect)

#croc.views<-list(ventral=croc.vent, lateral=croc.lateral, dorsal=croc.dorsal)
#save(croc.views,file="~/Google Drive/_UCL/Projects/Archosaur Modularity/Skull Figures/croc.3d.params.R")
load("~/Google Drive/_UCL/Projects/Archosaur Modularity/Skull Figures/croc.3d.params.R")

open3d(zoom=croc.views$lateral$zoom,userMatrix = croc.views$lateral$userMatrix,windowRect = croc.views$lateral$windowRect)
shade3d(m1,color='#E4D1C0')
spheres3d(newarray[-leftside,,1],col=modulecolors,radius=2.5)
rgl.snapshot("~/Google Drive/_UCL/Projects/Archosaur Modularity/Skull Figures/croc_modules_lateral.png")

specsnosph<-good_specs2[-which(module_defs$full_model==19),,]
open3d();spheres3d(specsnosph[module_defs_nosph$full_model==20,,1],col="black",radius=2.5)
writePLY("/Users/felice/Data/making cover image/croc/module19.ply")



open3d(zoom=croc.views$lateral$zoom,userMatrix = croc.views$lateral$userMatrix,windowRect = croc.views$lateral$windowRect)
shade3d(m1,color='#E4D1C0')
spheres3d(slided4$dataslide[,,1],col="blue",radius=1.6)
spheres3d(slided4$dataslide[c(1:100),,1],col="red",radius=1.8)
spheres3d(slided4$dataslide[c(101:735),,1],col="gold",radius=1.8)

dim(specsnosph)
#procrustes

check.modules(newarray,path="./ply/",suffix=".ply",pt.size=2,render="s",module.colors = modulecolors,begin=1)


GPA_40_crocs <- gpagen(newarray)
GPA_38_crocs <- gpagen(newarray[,,-which(dimnames(newarray)[[3]]%in%c("Crocodylus_acutus_AMNH_7857_cranium","Crocodilus_biporcatus_cranium_MNHN_A_5316"))])

#remove lefts 
proc_aligned_rightside<-GPA_38_crocs$coords[-bilats[,2],,]*-1
plotTangentSpace(proc_aligned_rightside)

modulecolors<-as.factor(module_defs_right$full_model)
levels(modulecolors)<-color_table$croc_colors
spheres3d(proc_aligned_rightside[,,1],color=modulecolors,radius=.0007)
#get correlation matrix
corrmat1<-dotcorr(proc_aligned_rightside)
#run emmli:
em1<-EMMLi(corr = corrmat1, N_sample = 38, mod = as.data.frame(module_defs_right[,-c(1,2,4)]))
source("~/Google Drive/_UCL/Projects/Improving_EMMLi/V3/EMMLi_CVv3.R")
source("~/Google Drive/_UCL/Projects/Improving_EMMLi/V3/EIC_CVv3.R")
emv3<-EMMLi2(corr = corrmat1, N_sample = 38, mod = as.data.frame(module_defs_right[,-c(1,2,4)]))
modelsA<-data.frame(rep(module_defs_right%>%pull(.,description),each=3) , apply(as.data.table(module_defs_right[,-c(1:4)]), 2, rep, each=3))

emv32<-EIC_EMMLi(data = two.d.array(proc_aligned_rightside), nboot = 1000L, models =  modelsA)

nmodules=20
#rho list directs the code to the results. this will be the "rhos_best" thing for subsampled emmli
rholist<-t( em1$rho$`full_model.sep.Mod + sep.between`)
words<-strsplit(rownames(rholist), split = " ")
corrmat_new<-matrix(data = NA, nrow = nmodules, ncol = nmodules)
for (i in 1:(length(words)-1)){
  if (length(words[[i]]) == 2){
    corrmat_new[as.numeric(words[[i]][2]),as.numeric(words[[i]][2])] <- rholist[i,2]
  }
  if (length(words[[i]]) == 3){
    corrmat_new[as.numeric(words[[i]][1]),as.numeric(words[[i]][3])] <- rholist[i,2]
    corrmat_new[as.numeric(words[[i]][3]),as.numeric(words[[i]][1])] <- rholist[i,2]
  }
}
corrmat_new <- corrmat_new[-19,-19]
modnames<-c("premax_vent","max_vent","max_dorsal","premax_dorsal","nasal","frontal","parietal","supra_occ","occ_condyle","basiocc","pterygoid","palatine","jjoint","jugal","squamosal","lac/prefont","ectoptery","ptery_flange","postorbital") 
rownames(corrmat_new)<-colnames(corrmat_new)<-modnames
#write.csv(corrmat_new,file = "/Users/rnf/Google Drive/_UCL/Projects/Crocs 2019/Manuscript/Results Tables/Emmli_crocs_nonphylo_oct4.csv")
within<-diag(corrmat_new)
between<-corrmat_new
#load layout
layout<-as.matrix(read.csv("~/Google Drive/NHM/crocs/data/croc_grid_layout2.csv",header = FALSE))
qgraph(between, shape="circle", posCol="#056571",  labels = rownames(corrmat_new), vsize=(within)*10, diag = FALSE, title="rho values",layout=layout)


###
#phylo version

#load trees
tree1<-read.nexus("/Users/felice/Dropbox/Crocs/Croc_FBD_Project/Data/FBD_Bayestraits/Tree_1_FBD_MCC.nex")
tree2<-read.nexus("/Users/felice/Dropbox/Crocs/Croc_FBD_Project/Data/FBD_Bayestraits/Tree_2_FBD_MCC.nex")
tree3<-read.nexus("/Users/felice/Dropbox/Crocs/Croc_FBD_Project/Data/FBD_Bayestraits/Tree_3_FBD_MCC.nex")
tree4<-read.nexus("/Users/felice/Dropbox/Crocs/Croc_FBD_Project/Data/FBD_Bayestraits/Tree_4_FBD_MCC.nex")

treelist<-list(tree1,tree2,tree3,tree4)

#change the names of the data to match the tree:

for (i in 1:nrow(completeness)){
    dimnames(proc_aligned_rightside)[[3]][which(dimnames(proc_aligned_rightside)[[3]]==completeness$Filename[i])]<-completeness$Tip_label[i]
}

nc1 <- name.check(tree1,two.d.array(proc_aligned_rightside))

treelist2 <- lapply(treelist,drop.tip,nc1$tree_not_data)

for (i in 1:length(treelist2)){
  treelist2[[i]]$root.time <- max(nodeHeights(treelist2[[i]]))
}

geoscalePhylo(tree=treelist2[[4]],cex.tip=0.7)

coords.2d <- two.d.array(proc_aligned_rightside)

contrasts_tree1 <- apply(coords.2d, 2, FUN=pic, phy=treelist2[[1]]) %>% arrayspecs(., p=dim(proc_aligned_rightside)[1], k = 3)
contrasts_tree2 <- apply(coords.2d, 2, FUN=pic, phy=treelist2[[2]]) %>% arrayspecs(., p=dim(proc_aligned_rightside)[1], k = 3)
contrasts_tree3 <- apply(coords.2d, 2, FUN=pic, phy=treelist2[[3]]) %>% arrayspecs(., p=dim(proc_aligned_rightside)[1], k = 3)
contrasts_tree4 <- apply(coords.2d, 2, FUN=pic, phy=treelist2[[4]]) %>% arrayspecs(., p=dim(proc_aligned_rightside)[1], k = 3)
 
corrmat_contrasts_t1<-dotcorr(contrasts_tree1)
corrmat_contrasts_t2<-dotcorr(contrasts_tree2)
corrmat_contrasts_t3<-dotcorr(contrasts_tree3)
corrmat_contrasts_t4<-dotcorr(contrasts_tree4)

EMMLI_PIC_tree1 <- EMMLi(corr = corrmat_contrasts_t1, N_sample = 38, mod = as.data.frame(module_defs_right[,-c(1,2,4)])) 
EMMLI_PIC_tree2 <- EMMLi(corr = corrmat_contrasts_t2, N_sample = 38, mod = as.data.frame(module_defs_right[,-c(1,2,4)])) 
EMMLI_PIC_tree3 <- EMMLi(corr = corrmat_contrasts_t3, N_sample = 38, mod = as.data.frame(module_defs_right[,-c(1,2,4)])) 
EMMLI_PIC_tree4 <- EMMLi(corr = corrmat_contrasts_t4, N_sample = 38, mod = as.data.frame(module_defs_right[,-c(1,2,4)])) 


EMMLI_PIC_tree1A <- EMMLi(corr = corrmat_contrasts_t1, N_sample = 38, mod = as.data.frame(module_defs_right[,c(3,6)])) 
EMMLI_PIC_tree2A <- EMMLi(corr = corrmat_contrasts_t2, N_sample = 38, mod = as.data.frame(module_defs_right[,c(3,6)])) 
EMMLI_PIC_tree3A <- EMMLi(corr = corrmat_contrasts_t3, N_sample = 38, mod = as.data.frame(module_defs_right[,c(3,6)])) 
EMMLI_PIC_tree4A <- EMMLi(corr = corrmat_contrasts_t4, N_sample = 38, mod = as.data.frame(module_defs_right[,c(3,6)])) 

#this will force emmli to give me results for the hypothesis with merged occipitals and mergeed pterygoids
EMMLI_PIC_tree1.for.icvm <- EMMLi(corr = corrmat_contrasts_t1, N_sample = 38, mod = as.data.frame(module_defs_right[,c(3,14)])) 
EMMLI_PIC_tree2.for.icvm <- EMMLi(corr = corrmat_contrasts_t2, N_sample = 38, mod = as.data.frame(module_defs_right[,c(3,14)])) 
EMMLI_PIC_tree3.for.icvm <- EMMLi(corr = corrmat_contrasts_t3, N_sample = 38, mod = as.data.frame(module_defs_right[,c(3,14)])) 
EMMLI_PIC_tree4.for.icvm <- EMMLi(corr = corrmat_contrasts_t4, N_sample = 38, mod = as.data.frame(module_defs_right[,c(3,14)])) 

EMMLI_PIC_tree1sub <- subsampleEMMLi(landmarks = contrasts_tree1, models = as.data.frame(module_defs_right[,c(3,15,6)]), min_landmark = 10, fractions = 0.5863, nrep = 100) 
EMMLI_PIC_tree2sub <- subsampleEMMLi(landmarks = contrasts_tree2, models = as.data.frame(module_defs_right[,c(3,15,6)]), min_landmark = 10, fractions = 0.5863, nrep = 100) 
EMMLI_PIC_tree3sub <- subsampleEMMLi(landmarks = contrasts_tree3, models = as.data.frame(module_defs_right[,c(3,15,6)]), min_landmark = 10, fractions = 0.5863, nrep = 100) 
EMMLI_PIC_tree4sub <- subsampleEMMLi(landmarks = contrasts_tree4, models = as.data.frame(module_defs_right[,c(3,15,6)]), min_landmark = 10, fractions = 0.5863, nrep = 100) 


CR.PIC_tree1<-modularity.test(contrasts3d,partition.gp = new.hypothesis)
 
#rho list directs the code to the results. this will be the "rhos_best" thing for subsampled emmli
nmodules=20
rholist<-t(EMMLI_PIC_tree4.for.icvm$rho$`occipitals_merged_and_pterygoids_merged.sep.Mod + sep.between`)
words<-strsplit(rownames(rholist), split = " ")
corrmat_new<-matrix(data = NA, nrow = nmodules, ncol = nmodules)
for (i in 1:(length(words)-1)){
  if (length(words[[i]]) == 2){
    corrmat_new[as.numeric(words[[i]][2]),as.numeric(words[[i]][2])] <- rholist[i,2]
  }
  if (length(words[[i]]) == 3){
    corrmat_new[as.numeric(words[[i]][1]),as.numeric(words[[i]][3])] <- rholist[i,2]
    corrmat_new[as.numeric(words[[i]][3]),as.numeric(words[[i]][1])] <- rholist[i,2]
  }
}
corrmat_new <- corrmat_new[-19,-19]
modnames<-c("premax_vent","max_vent","max_dorsal","premax_dorsal","nasal","frontal","parietal","supra_occ","occ_condyle","basiocc","pterygoid","palatine","jjoint","jugal","squamosal","lac/prefont","ectoptery","ptery_flange","postorbital") 
rownames(corrmat_new)<-colnames(corrmat_new)<-modnames
#write.csv(corrmat_new,file = "~/Google Drive/NHM/crocs/results/Emmli_crocs_phylo_tree1.csv")
within<-diag(corrmat_new)
between<-corrmat_new
#load layout
qgraph(between, shape="circle", posCol="#056571",  labels = rownames(corrmat_new), vsize=(within)*10, diag = FALSE, 
       title="rho values croc tree4",layout=layout)

#export results for the occ and ptery merged 
nmodules=20
rholist<-t(EMMLI_PIC_tree3.for.icvm$rho$`occipitals_merged_and_pterygoids_merged.sep.Mod + sep.between`)
words<-strsplit(rownames(rholist), split = " ")
corrmat_new<-matrix(data = NA, nrow = nmodules, ncol = nmodules)
for (i in 1:(length(words)-1)){
  if (length(words[[i]]) == 2){
    corrmat_new[as.numeric(words[[i]][2]),as.numeric(words[[i]][2])] <- rholist[i,2]
  }
  if (length(words[[i]]) == 3){
    corrmat_new[as.numeric(words[[i]][1]),as.numeric(words[[i]][3])] <- rholist[i,2]
    corrmat_new[as.numeric(words[[i]][3]),as.numeric(words[[i]][1])] <- rholist[i,2]
  }
}
corrmat_new <- corrmat_new[-19,-19]
modnames<-c("premax_vent","max_vent","max_dorsal","premax_dorsal","nasal","frontal","parietal","supra_occ","occ_condyle","basiocc","pterygoid","palatine","jjoint","jugal","squamosal","lac/prefont","ectoptery","ptery_flange","postorbital") 
rownames(corrmat_new)<-colnames(corrmat_new)<-modnames
write.csv(corrmat_new,file = "~/Google Drive/NHM/crocs/results/Emmli_crocs_phylo_tree3.occ.ptery.merged.csv")
within<-diag(corrmat_new)
between<-corrmat_new
#load layout
qgraph(between, shape="circle", posCol="#056571",  labels = rownames(corrmat_new), vsize=(within)*10, diag = FALSE, 
       title="rho values croc tree4",layout=layout)






#FOR BIRD MODULARITY:
nmodules=8

rholist<-t(EMMLI_PIC_tree4sub$`0.5863`$rhos_best$`full_model_Birds.sep.Mod + sep.between`)
words<-strsplit(rownames(rholist), split = " ")
corrmat_new<-matrix(data = NA, nrow = nmodules, ncol = nmodules)
for (i in 1:(length(words)-1)){
  if (length(words[[i]]) == 2){
    corrmat_new[as.numeric(words[[i]][2]),as.numeric(words[[i]][2])] <- rholist[i,2]
  }
  if (length(words[[i]]) == 3){
    corrmat_new[as.numeric(words[[i]][1]),as.numeric(words[[i]][3])] <- rholist[i,2]
    corrmat_new[as.numeric(words[[i]][3]),as.numeric(words[[i]][1])] <- rholist[i,2]
  }
}

modnames<-c("Rostrum","Vault","Basiphenoid","Palate","Pterygoid","Naris","Occipital","Quadrate") 
rownames(corrmat_new)<-colnames(corrmat_new)<-modnames
corrmat_new <- corrmat_new[-3,-3]
#write.csv(corrmat_new,file = "~/Google Drive/NHM/crocs/results/Emmli_crocs_phylo_subampled_tree4.csv")
within<-diag(corrmat_new)
between<-corrmat_new
#load layout
#layout<-as.matrix(read.csv("~/Google Drive/NHM/crocs/data/bird layout for crocs.csv",header=FALSE))
qgraph(between, shape="circle", posCol="#056571",  labels = rownames(corrmat_new), vsize=(within)*10, diag = FALSE, 
       title="croc data, subsampled, bird modules tree4",layout=layout)

cs<-GPA_38_crocs$Csize
allom<-procD.lm(proc_aligned_rightside~log(cs))
summary(allom)
shape.resid <- arrayspecs(allom$residuals, p=dim(proc_aligned_rightside)[1], k=dim(proc_aligned_rightside)[2])
adj.shape <- shape.resid + array(mshape(proc_aligned_rightside), dim(proc_aligned_rightside))

coords.2d.adj<-two.d.array(adj.shape)

contrasts_tree1.adj <- apply(coords.2d.adj, 2, FUN=pic, phy=treelist2[[1]]) %>% arrayspecs(., p=dim(adj.shape)[1], k = 3)
contrasts_tree2.adj <- apply(coords.2d.adj, 2, FUN=pic, phy=treelist2[[2]]) %>% arrayspecs(., p=dim(adj.shape)[1], k = 3)
contrasts_tree3.adj <- apply(coords.2d.adj, 2, FUN=pic, phy=treelist2[[3]]) %>% arrayspecs(., p=dim(adj.shape)[1], k = 3)
contrasts_tree4.adj <- apply(coords.2d.adj, 2, FUN=pic, phy=treelist2[[4]]) %>% arrayspecs(., p=dim(adj.shape)[1], k = 3)

corrmat_contrasts_t1.adj<-dotcorr(contrasts_tree1.adj)
corrmat_contrasts_t2.adj<-dotcorr(contrasts_tree2.adj)
corrmat_contrasts_t3.adj<-dotcorr(contrasts_tree3.adj)
corrmat_contrasts_t4.adj<-dotcorr(contrasts_tree4.adj)

EMMLI_PIC_tree1.adj <- EMMLi(corr = corrmat_contrasts_t1.adj, N_sample = 38, mod = as.data.frame(module_defs_right[,-c(1,2,4)])) 


#rho list directs the code to the results. this will be the "rhos_best" thing for subsampled emmli
nmodules=20
rholist<-t(EMMLI_PIC_tree1.adj$rho$`full_model.sep.Mod + sep.between`)
words<-strsplit(rownames(rholist), split = " ")
corrmat_new<-matrix(data = NA, nrow = nmodules, ncol = nmodules)
for (i in 1:(length(words)-1)){
  if (length(words[[i]]) == 2){
    corrmat_new[as.numeric(words[[i]][2]),as.numeric(words[[i]][2])] <- rholist[i,2]
  }
  if (length(words[[i]]) == 3){
    corrmat_new[as.numeric(words[[i]][1]),as.numeric(words[[i]][3])] <- rholist[i,2]
    corrmat_new[as.numeric(words[[i]][3]),as.numeric(words[[i]][1])] <- rholist[i,2]
  }
}
corrmat_new <- corrmat_new[-19,-19]
modnames<-c("premax_vent","max_vent","max_dorsal","premax_dorsal","nasal","frontal","parietal","supra_occ","occ_condyle","basiocc","pterygoid","palatine","jjoint","jugal","squamosal","lac/prefont","ectoptery","ptery_flange","postorbital") 
rownames(corrmat_new)<-colnames(corrmat_new)<-modnames
#write.csv(corrmat_new,file = "~/Google Drive/NHM/crocs/results/Emmli_crocs_phylo_tree1SIZEADJUSTED.csv")
within<-diag(corrmat_new)
between<-corrmat_new
#load layout
qgraph(between, shape="circle", posCol="#056571",  labels = rownames(corrmat_new), vsize=(within)*10, diag = FALSE, 
       title="rho values croc tree1 size adjusted",layout=layout)





#HOTDOTS WITH THIS DATA:
 ##################
 source("~/Google Drive/utlity_scripts/RNF_toolbox.R")
 
 x.croc<-two.d.array(proc_aligned_rightside)

 phy.parts <- phylo.mat(x.croc, treelist2[[3]])
 invC <- phy.parts$invC
 D.mat <- phy.parts$D.mat
 C = phy.parts$C
 
 global<-sig.calc(x.croc,invC,D.mat,Subset=TRUE)
 
 global.array1<-arrayspecs(global$R,p=1291,k=3)
 rates.vector<-colSums(matrix(diag(global$R), nrow=3))
 module.id<-module_defs_right$full_model
 module.id2<-module_defs_right$occipitals_merged
 
 variances<-rowSums(apply(proc_aligned_rightside ,c(1,2),var))
 
 per.lm.rates<-cbind(rates.vector,module.id2,variances)
 rates.and.vars<-as.data.frame(per.lm.rates)
 
 
 cols1<-colorRampPalette(c("#6e016b","#0c2c84","#225ea8","#005a32","#ffff00","#fe9929","#fc4e2a","red"))
 
 cols<-cols1(100)
 
 #calculate log rates:
 x=(log10(rates.vector))
 xlims<-NULL
 tol <- 1e-06
 xlims <- range(x) + c(-tol, tol)
 nbin=100
 breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
 whichColor <- function(p, cols, breaks) {
   i <- 1
   while (p >= breaks[i] && p > breaks[i + 1]) i <- i +
       1
   cols[i]
 }
 ratecolors <- sapply(x, whichColor, cols = cols, breaks = breaks)
 
 #calculate log variance:
 
 x=(log10(variances))
 xlims<-NULL
 tol <- 1e-06
 xlims <- range(x) + c(-tol, tol)
 nbin=100
 breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
 whichColor <- function(p, cols, breaks) {
   i <- 1
   while (p >= breaks[i] && p > breaks[i + 1]) i <- i +
       1
   cols[i]
 }
 dispcolors <- sapply(x, whichColor, cols = cols, breaks = breaks)
 
 
 open3d()
 mfrow3d(nr=1,nc=2)
 spheres3d(slided4.new$dataslide[-which(module_defs_nosph$l_or_r=="l"),,1],radius = 2.5,col=dispcolors)
 shade3d(m1,col="white")
 title3d(main="variance")
 next3d()
 spheres3d(slided4.new$dataslide[-which(module_defs_nosph$l_or_r=="l"),,1],radius = 2.5,col=ratecolors)
 shade3d(m1,col="white")
 title3d(main="rates")
 
 
 
 
 
 
 
 
 
 
 
 
 
 library(landvR)
 differences_from_mean <- coordinates.difference(coordinates = GPA_38_crocs$coords[-bilats[,2] ,,],
                                                 reference = GPA_38_crocs$consensus[-bilats[,2] ,],
                                                 type = "vector",
                                                 absolute.distance= FALSE)
 
 get.color.spectrum<-function(data, nbin=100,limits=NULL, colorlist = c("#6e016b","#0c2c84","#225ea8","#005a32","#ffff00","#fe9929","#fc4e2a","red")){
   x=data
   xlims<-NULL
   tol <- 1e-06
   cols1<-colorRampPalette(colorlist)
   cols<-cols1(100)
   if (is.null(limits)==TRUE){
     xlims <- range(x) + c(-tol, tol)}
   else{
     xlims <- limits}
   nbin=nbin
   breaks <- 0:nbin/nbin * (xlims[2] - xlims[1]) + xlims[1]
   whichColor <- function(p, cols, breaks) {
     i <- 1
     while (p >= breaks[i] && p > breaks[i + 1]) i <- i +
         1
     cols[i]
   }
   spec1colors <- sapply(x, whichColor, cols = cols, breaks = breaks)
 }
 
 open3d()
 spheres3d(GPA_38_crocs$consensus[-bilats[,2],],radius = .002)#,col=spec1colors)
 title3d(main="citipati")
 
 #save rgl window parameters
 
 zoom<-par3d()$zoom
 userMatrix<-par3d()$userMatrix
 windowRect<-par3d()$windowRect
 
 folderA<-"~/Google Drive/_UCL/Projects/Crocs 2019/distfrommeanplots/"
 for (i in 1:length(differences_from_mean)){
   datcol1<-get.color.spectrum(differences_from_mean[[i]][,1])
   name<-names(differences_from_mean)[i]
   open3d(zoom = zoom, userMatrix = userMatrix, windowRect=windowRect)
   spheres3d(GPA_38_crocs$consensus[-bilats[,2],],radius = .00075,col=datcol1)
   rgl.snapshot(filename = paste(folderA, names(differences_from_mean)[i],".png",sep=""))
   rgl.close()
 }
 
 
 overall.max.dist<-sapply(differences_from_mean,function(i) apply(i,2,max))
 overall.min.dist<-sapply(differences_from_mean,function(i) apply(i,2,min))
 overall.max.dist<-max(overall.max.dist[1,])
 overall.min.dist<-max(overall.min.dist[1,])
 rangex<-c(overall.min.dist,overall.max.dist)+c(-1e-06,1e-06)
 
 for (i in 1:length(differences_from_mean)){
   datcol1<-get.color.spectrum(differences_from_mean[[i]][,1],limits = rangex)
   name<-names(differences_from_mean)[i]
   open3d(zoom = zoom, userMatrix = userMatrix, windowRect=windowRect)
   spheres3d(dinos.only.f.v.occ.jj.bilateral.GPA.2$consensus[c(1:n.surf.points),],radius = .002,col=datcol1)
   rgl.snapshot(filename = paste(folderA, names(differences_from_mean)[i],"_rescaled",".png",sep=""))
   rgl.close()
 }
 
 
 
 
 
########################
 #geomorph rates!
########################
 
 
 
 reordered.data<-abind(proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==1),,],
                       proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==2),,],
                       proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==3),,],
                       proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==4),,],
                       proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==5),,],
                       proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==6),,],
                       proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==7),,],
                       proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==8),,],
                       proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==9),,],
                       proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==11),,],
                       proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==12),,],
                       proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==13),,],
                       proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==14),,],
                       proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==15),,],
                       proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==16),,],
                       proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==20),,],
                       along = 1)
 module.defs.for.rate.test<-sort(module_defs_right$occipitals_merged_and_pterygoids_merged)
 
 
 evoratesgeomorph_1<-compare.multi.evol.rates(A = reordered.data,
                                              gp =  as.factor(module.defs.for.rate.test),
                                              phy = treelist2[[1]],
                                              Subset = TRUE,
                                              iter = 999) 
 
 evoratesgeomorph_2<-compare.multi.evol.rates(A = reordered.data,
                                              gp =  as.factor(module.defs.for.rate.test),
                                              phy = treelist2[[2]],
                                              Subset = TRUE,
                                              iter = 999) 
 
 evoratesgeomorph_3<-compare.multi.evol.rates(A = reordered.data,
                                              gp =  as.factor(module.defs.for.rate.test),
                                              phy = treelist2[[3]],
                                              Subset = TRUE,
                                              iter = 999) 
 
 evoratesgeomorph_4<-compare.multi.evol.rates(A = reordered.data,
                                              gp =  as.factor(module.defs.for.rate.test),
                                              phy = treelist2[[4]],
                                              Subset = TRUE,
                                              iter = 999) 
 ##############
 #variance
 ################
 premax.vent.points <- proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==1),,] %>% two.d.array(.)
 max.vent.points <- proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==2),,] %>% two.d.array(.)
 max.dorsal.points <- proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==3),,] %>% two.d.array(.)
 premax.dorsal.points <-  proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==4),,] %>% two.d.array(.)
 nasal.points <- proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==5),,] %>% two.d.array(.)
 frontal.points <- proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==6),,] %>% two.d.array(.)
 parietal.points <- proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==7),,] %>% two.d.array(.)
 occiput.points <- proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==8),,] %>% two.d.array(.)
 pterygoid.points <- proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==11),,] %>% two.d.array(.)
 palatine.points <- proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==12),,] %>% two.d.array(.)
 jawjoint.points <- proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==13),,] %>% two.d.array(.)
 jugal.points <- proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==14),,] %>% two.d.array(.)
 squamosal.points <- proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==15),,] %>% two.d.array(.)
 lacrimal.points <- proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==16),,] %>% two.d.array(.)
 postorb.points <- proc_aligned_rightside[which(module_defs_right$occipitals_merged_and_pterygoids_merged==20),,] %>% two.d.array(.)
 
 
 premax.vent.disp<-morphol.disparity(premax.vent.points~1)
 premax.dorsal.disp<-morphol.disparity(premax.dorsal.points~1)
 max.vent.disp<-morphol.disparity(max.vent.points~1)
 max.dorsal.disp<-morphol.disparity(max.dorsal.points~1)
 nasal.disp<-morphol.disparity(nasal.points~1) 
 frontal.disp<-morphol.disparity(frontal.points~1) 
 parietal.disp<-morphol.disparity(parietal.points~1) 
 occiput.disp<-morphol.disparity(occiput.points~1)
 pterygoid.disp<-morphol.disparity(pterygoid.points~1)
 palatine.disp<-morphol.disparity(palatine.points~1)
 jugal.disp<-morphol.disparity(jugal.points~1)
 jawjoint.disp<-morphol.disparity(jawjoint.points~1) 
 postorb.disp<-morphol.disparity(postorb.points~1)
 squamosal.disp<-morphol.disparity(squamosal.points~1)
 lacrimal.disp<-morphol.disparity(lacrimal.points~1) 
 #############
 
 #loading in the results from the 
rate.disp.table<-read_csv("~/Dropbox/Crocs/Rates_Results/croc_geomorph_rates_tableOct2wide.csv")
 library(ggthemes)
 library(ggrepel)
 library(gridExtra)
 
 g1<-ggplot(rate.disp.table,aes(emmli_rho,geomorph_rates))+
   geom_point(size = 1)+labs(x="Within Module Correlation",y="Evolutionary Rate")+
   theme_bw()+
   theme(aspect.ratio=1,text=element_text(family="Arial",size=6))+
   geom_text_repel(aes(emmli_rho,geomorph_rates_trad,label=Module), size=2,family="Arial")
 #geom_text(aes(label=module_name),hjust="inward", vjust=-.20);g1
 gA<-ggplot(rate.disp.table,aes(emmli_rho,disparities/n_points))+
   geom_point(size = 1)+labs(x="Within Module Correlation",y="Disparity")+
   theme_bw()+
   theme(aspect.ratio=1,text=element_text(family="Arial",size=6))+
   geom_text_repel(aes(emmli_rho,disparities/n_points,label=Module),size=2,family="Arial")
 #geom_text(aes(label=module_name),hjust="inward", vjust=-.20);gA
 g2<-ggplot(rate.disp.table,aes(geomorph_rates,disparity/n_landmarks))+
   geom_point(size = 1)+labs(x="Rates",y="Disparity")+
   theme_bw()+ 
   theme(aspect.ratio=1,text=element_text(family="Arial",size=6))+
   geom_text_repel(aes(geomorph_rates,disparity/n_landmarks,label=Module),size=2,family="Arial")
 #geom_text(aes(label=module_name),hjust="inward", vjust=-.20);g2
 
 
 xplot<-grid.arrange(g1,gA, g2, nrow=1)
 
 
 ###########
 
 
 rate.disp.table<-read_csv("~/Dropbox/Crocs/Rates_Results/croc_geomorph_rates_tableOct2.csv")
 
 ratedisp.2 <-  rate.disp.table %>% spread(., key=key, value) 
 ratedisp.2$tree <- as.factor(ratedisp.2$tree)
 
 
 library(ggthemes)
 library(ggrepel)
 library(gridExtra)
 
 g1<-ggplot(ratedisp.2,aes(x = rho,y = rate, color= Module, shape = tree))+
   geom_point(size = 1)+labs(x="Within Module Correlation",y="Evolutionary Rate")+
   theme_bw()+
   theme(aspect.ratio=1,text=element_text(family="Arial",size=6))+
   scale_color_manual(values=unique(ratedisp.2$color))+
   geom_text_repel(aes(rho,rate,label=Module), size=1,family="Arial", color="black")
 #geom_text(aes(label=module_name),hjust="inward", vjust=-.20);g1
 gA<-ggplot(ratedisp.2,aes(rho,disp/n_landmarks, color= Module, shape = tree))+
   geom_point(size = 1)+labs(x="Within Module Correlation",y="Disparity")+
   theme_bw()+
   theme(aspect.ratio=1,text=element_text(family="Arial",size=6))+
   scale_color_manual(values=unique(ratedisp.2$color))+
   geom_text_repel(aes(rho,disp/n_landmarks,label=Module),size=1,family="Arial", color="black")
 #geom_text(aes(label=module_name),hjust="inward", vjust=-.20);gA
 g2<-ggplot(ratedisp.2,aes(rate,disp/n_landmarks, color= Module, shape = tree))+
   geom_point(size = 1)+labs(x="Rates",y="Disparity")+
   theme_bw()+ 
   theme(aspect.ratio=1,text=element_text(family="Arial",size=6))+
   scale_color_manual(values=unique(ratedisp.2$color))+
   geom_text_repel(aes(rate,disp/n_landmarks,label=Module),size=1,family="Arial", color="black")
 #geom_text(aes(label=module_name),hjust="inward", vjust=-.20);g2
 
 
 xplot<-grid.arrange(g1+ theme(legend.position="none"),gA+ theme(legend.position="none"), g2+ theme(legend.position="none"), nrow=1)
 
#ggsave("/Users/rnf/Dropbox/Crocs/R_plots/rate disp and integration plots all.pdf",xplot,device = cairo_pdf, width=17.8,height=6,units="cm")
 
 #HOTSCATs
 
 
 ################################
 #simulate:
 cov.mat<-matrix(0,nrow=length(diag(global$R)),ncol=length(diag(global$R)))
 diag(cov.mat)<-diag(global$R)
 phy<-treelist2[[1]]
 
 iter=100
 k=3
 
 system.time(
   simdat <- sim.char(phy=phy,par=cov.mat, nsim=iter,model="BM") 
 )
 
 expectedmat<-matrix(nrow=(dim(simdat)[2]/k),ncol=iter)
 for (i in 1:iter){ 
   temp<-arrayspecs(simdat[,,i],p=1291,k=3)
   expected<-rowSums(apply(temp,c(1,2),var))
   expectedmat[,i]<-expected
 }
 colnames(expectedmat)<-paste(c(1:iter))
 
 simulatedvariance<-tbl_df(cbind(rates.vector,expectedmat))
 
 simvarclean<-gather(simulatedvariance, "simnumber", "variance", 2:(iter+1))
 
 # Create prediction interval 
 m1<-lm(variance~rates.vector, data=simvarclean)
 newx <- seq(min(simvarclean$rates.vector), max(simvarclean$rates.vector), length.out = 100)
 pred_interval <- predict(m1, newdata=data.frame(rates.vector=newx), interval="prediction", level = 0.95)
 pred_interval <- as.data.frame(pred_interval)
 pred_interval$rates = newx
 
 #ggthemr_reset()
 p <-ggplot(rates.and.vars, aes(x=rates.vector, y=variances, color=as.factor(module.id2))) 
 p<- p + geom_point(size=1.8)+
   scale_color_manual(values=color_table$croc_colors[-c(9,10)],
                      labels = c("Premaxilla (ventral)","Maxilla (ventral)","Maxilla (dorsal)","Premaxilla (dorsal)","Nasal","Frontal","Parietal","Occiput","Pterygoid","Palatine","Jaw Joint","Jugal+Quadratojugal","Squamosal","Lacrimal+Prefrontal","Ectopterygoid","Pterygoid Flange","Postorbital"),
                      name = "Group")
 p + theme_light()+ stat_smooth(method="lm",
                                se=TRUE,
                                fill="grey",
                                formula=y ~ x,
                                colour="red")+ 
   stat_smooth(mapping=aes(x=rates.vector,y=variance),
               linetype=2,
               se=TRUE,
               data=simvarclean,
               method=lm,fill="grey",
               formula=y ~ x, inherit.aes = FALSE)+
   geom_ribbon(mapping=aes(x=rates,ymin = lwr, ymax = upr),
               data=pred_interval,
               fill = "grey",
               alpha = 0.2,
               inherit.aes = FALSE)+
   theme(aspect.ratio = 1,legend.position="bottom")+
   labs(x="Evolutionary Rate (Sigma)",y="Procrustes Variance")
 #
 
 
 model1 <- lm(variances~rates.vector, data= rates.and.vars)
 summary(model1)
 
 
 