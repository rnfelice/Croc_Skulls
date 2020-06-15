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
library(phytools)
library(here)
library(doParallel)
library(hot.dots) #install_github("rnfelice/hot.dots")
library(SURGE) #install_github("rnfelice/SURGE")


#load coordinate data
#sliding landmarks have already been slid to minimize bending energy
n.lm <- 1332 #1332 landmarks total
n.specs <- 43 #number of specimens
n.dim <- 3 #this is 3d data
coords.raw <- read.csv(file=here("Data", "slid_coordinates.csv"),row.names = 1)
#convert from matrix to array
coords.raw <- arrayspecs(coords.raw, p = n.lm, k = n.dim)             
#estimate missing points
coords <- fixLMtps(coords.raw)$out


#import module definitions
module_defs <- read_csv(here("Data","module_data.csv"))
load(file= here("Data","curve_data.Rdata"))

#### Mirror Data ####
bilat.landmarks <- cbind(which(module_defs$l_or_r=="r"),which(module_defs$l_or_r=="l"))

bilat.landmarks <- bilat.landmarks[c(1:41),]

midline<- which(module_defs$l_or_r=="m")#a list of which landmarks fall on the midline
#a list of the CORRESPONDING left side landmarks. (the order of right.lm and left.lm must represent bilateral pairs of landmarks)
surp <- c((max(unlist(curve_list))+1):dim(coords)[1])

all.curves.present<- c(1:length(curve_list))
right.curves <- all.curves.present[-c(4,13,16,20,24,30,32)]#a list of the curves which are going to be mirrored to the other side (excludes midline curves)
right.curve.list<-unlist(curve_list[right.curves])
rightside<-c(bilat.landmarks[,1],right.curve.list,surp)
num.missing<-(length(rightside)-length(bilat.landmarks[,2]))
blanks<-c((dim(coords)[1]+1):(dim(coords)[1]+num.missing))         
leftside<-c(bilat.landmarks[,2],blanks)

module_defs_right <- module_defs %>% filter(., l_or_r != "l")

add_col_or_row = function(x, n = 1, add_col = T, fill = 0)
{
  m1 = matrix(x, ncol = if(add_col) nrow(x) * ncol(x) else nrow(x), byrow = T)
  m2 = matrix(fill, nrow = if(add_col) dim(x)[3] else prod(dim(x)[-1]), 
              ncol = if(add_col) nrow(x) * n else n)
  array(t(cbind(m1, m2)), 
        c(nrow(x) + ((!add_col) * n), ncol(x) + (add_col * n), dim(x)[3]))
}

coords2<-add_col_or_row(coords,n=num.missing,add_col=FALSE,fill=NA)
dimnames(coords2)[3]<-dimnames(coords)[3]
bilats<-cbind(rightside, leftside)
newarray<-paleomorph::mirrorfill(coords2,l1=as.integer(midline),l2=bilats)
dimnames(newarray)[3]<-dimnames(coords)[3]

#### Procrustes Superimposition ####
Y.gpa <- gpagen(newarray)

proc_aligned_rightside<-Y.gpa$coords[-bilats[,2],,]

pca_results <- gm.prcomp(proc_aligned_rightside)

module_defs_right <-  module_defs %>% filter(., l_or_r != "l")

modulecolors<-as.factor(module_defs_right$full_model)
color.palette <- c( "#6a3d9a","dimgrey","#fb9a99",  "gold", "#009E73",  "#D55E00", "#CC79A7", "cyan2",  "#e31a1c", "#0072B2", "#b2df8a", "#E69F00",  "whitesmoke" ,  "deeppink",   "#a6cee3",   "#F0E442","blue","red","brown", "black")
levels(modulecolors)<-color.palette
  
open3d()
spheres3d(proc_aligned_rightside[,,1],col=modulecolors,radius=.005)


tree1 <- read.nexus(here("Trees","Tree_1_FBD_MCC.nex"))
tree2 <- read.nexus(here("Trees","Tree_2_FBD_MCC.nex"))
tree3 <- read.nexus(here("Trees","Tree_3_FBD_MCC.nex"))
tree4 <- read.nexus(here("Trees","Tree_4_FBD_MCC.nex"))

treelist<-list(tree1,tree2,tree3,tree4)
par(mfrow=c(2,2))
plot(treelist[[1]])
plot(treelist[[2]])
plot(treelist[[3]])
plot(treelist[[4]])
par(mfrow=c(1,1))
#change the names of the data to match the tree:


nc1 <- name.check(tree1,two.d.array(proc_aligned_rightside))

treelist2 <- lapply(treelist,drop.tip,nc1$tree_not_data)


nc1 <- name.check(treelist[[1]],two.d.array(proc_aligned_rightside))

treelist2 <- lapply(treelist,drop.tip,nc1$tree_not_data)

for (i in 1:length(treelist2)){
  treelist2[[i]]$root.time <- max(nodeHeights(treelist2[[i]]))
}


#geoscalePhylo(tree=treelist2[[4]],cex.tip=0.7)

####Per-landmark Rate and Variance####
rateslist <- lapply(1:4, function(k) hot.dots::per_lm_rates(proc_aligned_rightside, treelist2[[k]]))
perlmvar <- hot.dots::per_lm_variance(proc_aligned_rightside)






########################
####geomorph rates####
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


####BayesTraits Rates Analysis####

#export trees and phyloPC scores
btfolder <- here("BayesTraits") 

runs <- c("a", "b") # do two runs of each so that you can check convergence at the end

for (i in 1:length(runs)){
  for (j in 1:length(treelist2)){
    
    phypc <- gm.prcomp(A = proc_aligned_rightside, phy = treelist2[[j]], align.to.phy = TRUE)
    
    
    scores <- phypc$x[,c(1:2)]*1000 #keeping 2 pc axes for 95% of variance, multiply by 1000 to make sure tiny values dont underflow in BayesTraits
    
    write.table(scores, file = paste0(btfolder,"/Phylo_PC_SCORES_","tree_",j,"_run_",runs[i],".txt"), quote = FALSE, col.names=FALSE)
    write.nexus(treelist2[[j]], file = paste0(btfolder,"/tree_",j,".nex"))
  }
}

#Run bayestraits

#### Caution: this is slow and will hog resources, only run if it you are sure you want to 
#import bayestraits control file
cmd_script <-  dir(here("BayesTraits"), pattern = "\\.cmd$")

ntrees <- 4
nruns <- 2
n_cores <- detectCores()
registerDoParallel(n_cores-2)
myOS <-  Sys.info()['sysname']

if(myOS == "Windows"){
  foreach(i = 1:length(nruns)) %dopar% {
    for (j in 1:ntrees){
      system2(
        paste0("cd ", btfolder, "/ & BayesTraitsV3.exe ", "tree_", j, ".nex", " ",
               "Phylo_PC_SCORES_","tree_",j,"_run_",runs[i],".txt", 
               " <", cmd_script), wait=FALSE, stdout = FALSE )
    }
  }
} else {
  btfolder <- paste0( "\'",btfolder, "\'") #need to add extra quotes around path name to make sure spaces in folder names dont mess things up on Unix systems
  foreach(i = 1:length(nruns)) %dopar% {
    for (j in 1:ntrees){
      system(
        paste0("cd ", btfolder, " && ./BayesTraitsV3 ", "tree_", j, ".nex", " ",
               "Phylo_PC_SCORES_","tree_",j,"_run_",runs[i],".txt",
               " <", cmd_script), wait=FALSE, ignore.stdout = TRUE )
    }
  }
}

##### Load and Plot BayesTraits Results####


# Load and plot BayesTraits Results ---------------------------------------


###########################################################
#Check convergence of likelihood:
test1 = tracePlots(file=here("BayesTraits","Phylo_PC_SCORES_tree_1_run_a.txt.VarRates"), plot = FALSE)
test2 = tracePlots(file=here("BayesTraits","Phylo_PC_SCORES_tree_1_run_b.txt.VarRates"), plot = FALSE)
my_list_of_chains = mcmc.list(list(test1[,c("Lh", "Lh...Prior", "No.Pram", "Alpha", "Sigma.2")], test2[,c("Lh", "Lh...Prior", "No.Pram", "Alpha", "Sigma.2")]))
gelman.diag(my_list_of_chains)

tree_1_BTraits<-BTRTools::rjpp(rjlog = here("BayesTraits","Phylo_PC_SCORES_tree_1_run_a.txt.VarRates"),
                                          rjtrees = here("BayesTraits","Phylo_PC_SCORES_tree_1_run_a.txt.Output.trees"),
                                          tree = tree_1_wholeskull)

cophylotrees <- cophylo(tree_1_wholeskull,tree_1_BTraits_wholeskull$meantree) ### keeps 2nd tree constant, change orientation of 1st to match
tree_1_wholeskull.a<- cophylotrees$trees[[1]]


{
  
  wholeskulltree.lab <- tree_1_wholeskull.a$tip.label
  wholeskulltree.lab.df <- data.table::data.table(label = wholeskulltree.lab, label2 = sub("_", " ",tree_1_wholeskull.a$tip.label)) #Make species names italic and remove underscore
  treedat_wholeskull_1<-data.frame(tree_1_wholeskull.a$edge,tree_1_BTraits_wholeskull$data$meanRate[-1],tree_1_BTraits_wholeskull$data$pScaled[-1])
  wholeskulltree.lab.df$label2[which(wholeskulltree.lab.df$label2=="Crocodylus acutus1")]<-"Crocodylus acutus"
  colnames(treedat_wholeskull_1)<-c("parent","node","meanrate","percentscaled")
  
  #cphylogram
  p<-ggtree(tree_1_wholeskull.a, aes(color = log(meanrate)), size=1)#, size = 0.5) #%<+% edge +
  p<- p %<+% treedat_wholeskull_1 +
    scale_colour_gradientn(colours = color3(100))+
    theme(legend.position="right")+
    labs(title="Whole Skull Rates- Tree 1",
         color="log(Rate)")
  p.labed<-p%<+% wholeskulltree.lab.df + geom_tiplab(aes(label = label2), size=5, color = "black", fontface="italic")+xlim(-1,250)
  p.labed
}

#ggsave(p.labed, filename="/Users/felice/Dropbox/Crocs/Croc_FBD_Project/Tree_figs/Rates_Trees/Whole Skull Rates- Tree 1.pdf", device="pdf")





# Convergence Test --------------------------------------------------------


conv.test <- search.conv(RR=RRates, y=PCA1$pc.scores[,c(1:4)], min.dim=5, foldername =  tempdir())
