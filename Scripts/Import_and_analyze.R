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
library(doParallel)
library(RRphylo)
library(patchwork)
library(gridExtra)
library(treeio)
library(extrafont)
#devtools::install_github("willgearty/deeptime")
library(deeptime)
#BiocManager::install("ggtree")
library(ggtree)
#devtools::install_github("rnfelice/hot.dots")
library(hot.dots)
#devtools::install_github("rnfelice/SURGE")
library(SURGE) 

library(here)

source(here("Scripts","utility_functions.R"))

# Load Data and Trees--------------------------------------------


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

# Mirror Data
bilat.landmarks <- cbind(which(module_defs$l_or_r=="r"),which(module_defs$l_or_r=="l"))

bilat.landmarks <- bilat.landmarks[c(1:41),]

midline<- which(module_defs$l_or_r=="m")#a list of which landmarks fall on the midline

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

# Procrustes superimposition
Y.gpa <- gpagen(newarray)

proc_aligned_rightside<-Y.gpa$coords[-bilats[,2],,]

pca_results <- gm.prcomp(proc_aligned_rightside)

module_defs_right <-  module_defs %>% filter(., l_or_r != "l")

modulecolors<-as.factor(module_defs_right$full_model)
color.palette <- c( "#6a3d9a","dimgrey","#fb9a99",  "gold", "#009E73",  "#D55E00", "#CC79A7", "cyan2",  "#e31a1c", "#0072B2", "#b2df8a", "#E69F00",  "whitesmoke" ,  "deeppink",   "#a6cee3",   "#F0E442","blue","red","brown", "black")
levels(modulecolors)<-color.palette
  
open3d()
spheres3d(proc_aligned_rightside[,,1],col=modulecolors,radius=.001)


tree1 <- read.nexus(here("Trees","Tree_1_FBD_MeanBL.nex"))
tree2 <- read.nexus(here("Trees","Tree_2_FBD_MeanBL.nex"))
tree3 <- read.nexus(here("Trees","Tree_3_FBD_MeanBL.nex"))
tree4 <- read.nexus(here("Trees","Tree_4_FBD_MeanBL.nex"))

treelist<-list(tree1,tree2,tree3,tree4)
par(mfrow=c(2,2))
plot(treelist[[1]])
plot(treelist[[2]])
plot(treelist[[3]])
plot(treelist[[4]])
par(mfrow=c(1,1))
#change the names of the data to match the tree:



nc1 <- name.check(treelist[[1]],two.d.array(proc_aligned_rightside))

treelist2 <- lapply(treelist,drop.tip,nc1$tree_not_data)

for (i in 1:length(treelist2)){
  treelist2[[i]]$root.time <- max(nodeHeights(treelist2[[i]]))
}


#geoscalePhylo(tree=treelist2[[4]],cex.tip=0.7)

# Test Phylogenetic Signal --------------------------------------------

physiglist <- lapply(1:4, function(k) geomorph::physignal(proc_aligned_rightside, treelist2[[k]]))
physiglist

# Test allometry --------------------------------------------

allometry1 <- geomorph::procD.lm(proc_aligned_rightside~log(Y.gpa$Csize))
summary(allometry1)
plotAllometry(allometry1, size=Y.gpa$Csize,method="RegScore")

allomtery.phylo <- lapply (1:4, function(k) geomorph::procD.pgls(proc_aligned_rightside~Y.gpa$Csize, phy=treelist2[[k]]))
summary(allomtery.phylo[[1]])
summary(allomtery.phylo[[2]])
summary(allomtery.phylo[[3]])
summary(allomtery.phylo[[4]])



# Per-landmark rate and variance --------------------------------------------

rateslist <- lapply(1:4, function(k) hot.dots::per_lm_rates(proc_aligned_rightside, treelist2[[k]]))
perlmvar <- hot.dots::per_lm_variance(proc_aligned_rightside)






# Per-region rate and variance --------------------------------------------




#rate

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




#variance

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



# Compare Variation Between Groups ----------------------------------------
eco <- read_csv(here("Data", "Ecology_data.csv"))

eco <- eco %>% filter(., Tip_label %in% row.names(pca_results$x)) %>% arrange(., Tip_label)
names(eco$extinct)<-eco$Tip_label
brevirostres <- rep("other", 43)
brevirostres[which(eco$Clade %in% c("Alligatoroidea" , "Crocodylidae"))] <- "Brevirostres"
names(brevirostres) <- eco$Tip_label
names(eco$Diet2)<-eco$Tip_label
diet <- eco$Diet2
gdf <- geomorph.data.frame(coords = proc_aligned_rightside[,,names(eco$extinct)], is_extinct = eco$extinct, clade = brevirostres, diet = diet, size = Y.gpa$Csize)


disparity_comp<-morphol.disparity(f1 = coords ~ is_extinct, groups = ~is_extinct, data = gdf)
disparity_comp2<-morphol.disparity(f1 = coords ~ clade, groups = ~clade, data = gdf)


# MANOVA for ecology  -----------------------------------------------------
manova.nonphylo<-procD.lm(f1 = coords ~ diet, data = gdf)

manova.phylo <- lapply(1:4, function(k) geomorph::procD.pgls(f1 = coords ~ diet, phy = treelist2[[k]], data = gdf))
summary(manova.phylo[[1]])
summary(manova.phylo[[2]])
summary(manova.phylo[[3]])
summary(manova.phylo[[4]])

# BayesTraits Rates Analysis ----------------------------------------------


#export trees and phyloPC scores
btfolder <- here("BayesTraits") 

runs <- c("a", "b") # do two runs of each so that you can check convergence at the end

for (i in 1:length(runs)){
  for (j in 1:length(treelist2)){
  
    phypc <- gm.prcomp(A = proc_aligned_rightside, phy = treelist2[[j]], align.to.phy = TRUE)
    
    #keeping 2 pc axes for 95% of variance, multiply by 1000 to make sure tiny values dont underflow in BayesTraits
    
    scores <- phypc$x[,c(1:which(cumsum(phypc$d/sum(phypc$d))>0.95)[1])]*1000 
    
    write.table(scores, file = paste0(btfolder,"/Phylo_PC_SCORES_","tree_",j,"_run_",runs[i],".txt"), quote = FALSE, col.names=FALSE)
    write.nexus(treelist2[[j]], file = paste0(btfolder,"/tree_",j,".nex"))
  }
}

#Run bayestraits

#!!!!! Caution: this is slow and will hog resources, only run if it you are sure you want to 
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
               "Phylo_PC_SCORES_","tree_",j,"_run_",runs[2],".txt",
               " <", cmd_script), wait=FALSE, ignore.stdout = TRUE )
    }
  }
}



# Load and plot BayesTraits Results ---------------------------------------
color3<-colorRampPalette(c("#0c2c84","#225ea8","#31a354","#ffff00","#fe9929","#fc4e2a","red","darkred"))


#Tree 1
##Check convergence of likelihood:
#test1 = tracePlots(file=here("BayesTraits","Phylo_PC_SCORES_tree_1_run_a.txt.VarRates.txt"), plot = FALSE)
#test2 = tracePlots(file=here("BayesTraits","Phylo_PC_SCORES_tree_1_run_b.txt.VarRates.txt"), plot = FALSE)
#my_list_of_chains = mcmc.list(list(test1[,c("Lh", "Lh...Prior", "No.Pram", "Alpha", "Sigma.2")], test2[,c("Lh", "Lh...Prior", "No.Pram", "Alpha", "Sigma.2")]))
#gelman.diag(my_list_of_chains)
#
#tree_1_BTraits<-BTRTools::rjpp(rjlog = here("BayesTraits","Phylo_PC_SCORES_tree_1_run_a.txt.VarRates.txt"),
#                                          rjtrees = here("BayesTraits","Phylo_PC_SCORES_tree_1_run_a.txt.Output.trees"),
#                                          tree = treelist2[[1]]) #this is your time scaled tree that was used to input into bayestraits
#

palatex<-c( "#D55E00",   "#009E73",  "#56B4E9")
tree_1_w_data <- add_rjpp_to_tree(tree_1_BTraits)
tree_1_w_data <- full_join(tree_1_w_data, mutate(eco, label = eco$Tip_label), by = "label")
threshold <- .15 # the minimum posterior probability you want to plot a symbol for 
p<-ggtree(tree_1_w_data, aes(color = log(meanRate)), size=1)+
  scale_colour_gradientn(colours = color3(100))+
  geom_tippoint(aes(fill=Diet2, x=x+4),color="black", shape=23)+
  scale_fill_manual(breaks = c("Carnivore","Omnivore/Herbivore","Piscivore"), values = palatex)+
  #theme(legend.position="top")+
  theme(legend.position=c(.32,.83),legend.direction = "horizontal",legend.box.background = element_rect(colour = "black",size =1))+
  scale_size(range = c(1,2))+ 
  labs(title="Tree 1",
       color="log(Rate)")+
  geom_nodepoint(aes(subset=ppRate>threshold, size = ppRate),color='black',fill="grey", shape=24)+
  scale_size(range = c(1,2))+ 
  geom_tiplab(label= sub("_", " ",tree_1_w_data@phylo$tip.label), size=3, offset=5.5, color = "black", family = "Arial", fontface="italic")+
  coord_cartesian(xlim = c(-230, 90), #you have to fiddle with these values to get your tip labels to show. the first value should be just before your root time, second value pads out space for tip labels
                  ylim = c(-2, 45), #first value makes room for geo timescale, second value is vertical space and should be a few more than your number of tips
                  expand = FALSE) +
 scale_x_continuous(breaks=-periods$max_age[c(1:5)], labels=periods$max_age[c(1:5)]) + 
 theme(panel.grid.major.x = element_line(colour="grey", size=0.5), legend.key.height =unit(.4,"cm"))#should also be modified based on your time scale limits
 p <- revts(p);p
 ptree1 <-  gggeo_scale(p, neg = FALSE, center_end_labels = TRUE, height = unit(1, "line"), size=3)
ptree1

p2<-p+geom_strip(taxa1 ="Crocodylus_mindorensis", taxa2 = "Prodiplocynodon_langi",  label="Crocodyloidea", offset = 60, offset.text = 3, barsize = 2, angle = 35, family = "Arial")+
  geom_strip(taxa1 ="Caiman_yacare", taxa2 =  "Leidyosuchus_canadensis",  label="Alligatoroidea", offset = 60, offset.text = 3, barsize = 2, angle = 35, family = "Arial", color="grey")+
  geom_strip(taxa1 ="Simosuchus_clarki" , taxa2 =  "Araripesuchus_wegeneri" ,  label="Notosuchia", offset = 60, offset.text = 3, barsize = 2, angle = 35, family = "Arial", color = "grey")+
  geom_strip(taxa1 ="Sarcosuchus_imperator", taxa2 =  "Pholidosaurus_sp",  label="Pholidosauridae", offset = 60, offset.text = 3, barsize = 2, angle = 35, family = "Arial")+
  geom_strip(taxa1 ="Pelagosaurus_typus" , taxa2 =  "Cricosaurus"  ,  label="Thalattosuchia", offset = 60, offset.text = 3, barsize = 2, angle = 35, family = "Arial")
  

ptree1b <-  gggeo_scale(p2, neg = FALSE, center_end_labels = TRUE, height = unit(1, "line"), size=3)
ptree1b 

ggsave(filename = here("Fig_Output", "rate_trees1a.pdf"),
       plot = ptree1b, device = cairo_pdf,
       width = 11, height = 11, units = "cm")

#Tree 2
#Check convergence of likelihood:
#test1 = tracePlots(file=here("BayesTraits","Phylo_PC_SCORES_tree_2_run_a.txt.VarRates.txt"), plot = FALSE)
#test2 = tracePlots(file=here("BayesTraits","Phylo_PC_SCORES_tree_2_run_b.txt.VarRates.txt"), plot = FALSE)
#my_list_of_chains = mcmc.list(list(test1[,c("Lh", "Lh...Prior", "No.Pram", "Alpha", "Sigma.2")], test2[,c("Lh", "Lh...Prior", "No.Pram", "Alpha", "Sigma.2")]))
#gelman.diag(my_list_of_chains)
#
#tree_2_BTraits<-BTRTools::rjpp(rjlog = here("BayesTraits","Phylo_PC_SCORES_tree_2_run_a.txt.VarRates.txt"),
#                               rjtrees = here("BayesTraits","Phylo_PC_SCORES_tree_2_run_a.txt.Output.trees"),
#                               tree = treelist2[[2]])

tree_2_w_data <- add_rjpp_to_tree(tree_2_BTraits)
threshold <- .15 # the minimum posterior probability you want to plot a symbol for 
p<-ggtree(tree_2_w_data, aes(color = log(meanRate)), size=1)+
  scale_colour_gradientn(colours = color3(100))+
  #theme(legend.position="top")+
  theme(legend.position=c(.32,.83),legend.direction = "horizontal",legend.box.background = element_rect(colour = "black",size =1))+
  scale_size(range = c(1,2))+ 
  labs(title="Tree 2",
       color="log(Rate)")+
  geom_nodepoint(aes(subset=ppRate>threshold, size = ppRate),color='black',fill="grey", shape=24)+
  geom_tiplab(label= sub("_", " ",tree_2_w_data@phylo$tip.label), size=2, color = "black", family = "Arial", fontface="italic")+
  coord_cartesian(xlim = c(-230, 90), #you have to fiddle with these values to get your tip labels to show. the first value should be just before your root time, second value pads out space for tip labels
                  ylim = c(-2, 45), #first value makes room for geo timescale, second value is vertical space and should be a few more than your number of tips
                  expand = FALSE) +
  scale_x_continuous(breaks=-periods$max_age[c(1:5)], labels=periods$max_age[c(1:5)]) + 
  theme( panel.grid.major.x = element_line(colour="grey", size=0.5))#should also be modified based on your time scale limits
p <- revts(p) 
ptree2 <- gggeo_scale(p, neg = FALSE, center_end_labels = TRUE, height = unit(1, "line"), size=3)
ptree2

#Tree 3
#Check convergence of likelihood:
#test1 = tracePlots(file=here("BayesTraits","Phylo_PC_SCORES_tree_3_run_a.txt.VarRates.txt"), plot = FALSE)
#test2 = tracePlots(file=here("BayesTraits","Phylo_PC_SCORES_tree_3_run_b.txt.VarRates.txt"), plot = FALSE)
#my_list_of_chains = mcmc.list(list(test1[,c("Lh", "Lh...Prior", "No.Pram", "Alpha", "Sigma.2")], test2[,c("Lh", "Lh...Prior", "No.Pram", "Alpha", "Sigma.2")]))
#gelman.diag(my_list_of_chains)
#
#tree_3_BTraits<-BTRTools::rjpp(rjlog = here("BayesTraits","Phylo_PC_SCORES_tree_3_run_a.txt.VarRates.txt"),
#                               rjtrees = here("BayesTraits","Phylo_PC_SCORES_tree_3_run_a.txt.Output.trees"),
#                               tree = treelist2[[3]])

tree_3_w_data <- add_rjpp_to_tree(tree_3_BTraits)
threshold <- .15 # the minimum posterior probability you want to plot a symbol for 
p<-ggtree(tree_3_w_data, aes(color = log(meanRate)), size=1)+
  scale_colour_gradientn(colours = color3(100))+
  #theme(legend.position="top")+
  theme(legend.position=c(.32,.83),legend.direction = "horizontal",legend.box.background = element_rect(colour = "black",size =1))+
  scale_size(range = c(1,2))+ 
  labs(title="Tree 3",
       color="log(Rate)")+
  geom_nodepoint(aes(subset=ppRate>threshold, size = ppRate),color='black',fill="grey", shape=24)+
  geom_tiplab(label= sub("_", " ",tree_3_w_data@phylo$tip.label), size=2, color = "black", family = "Arial", fontface="italic")+
  coord_cartesian(xlim = c(-230, 90), #you have to fiddle with these values to get your tip labels to show. the first value should be just before your root time, second value pads out space for tip labels
                  ylim = c(-2, 45), #first value makes room for geo timescale, second value is vertical space and should be a few more than your number of tips
                  expand = FALSE) +
  scale_x_continuous(breaks=-periods$max_age[c(1:5)], labels=periods$max_age[c(1:5)]) + 
  theme(panel.grid.major.x = element_line(colour="grey", size=0.5))#should also be modified based on your time scale limits
p <- revts(p) 
ptree3 <- gggeo_scale(p, neg = FALSE, center_end_labels = TRUE, height = unit(1, "line"), size=3)


#Tree 4
#Check convergence of likelihood:
#test1 = tracePlots(file=here("BayesTraits","Phylo_PC_SCORES_tree_4_run_a.txt.VarRates.txt"), plot = FALSE)
#test2 = tracePlots(file=here("BayesTraits","Phylo_PC_SCORES_tree_4_run_b.txt.VarRates.txt"), plot = FALSE)
#my_list_of_chains = mcmc.list(list(test1[,c("Lh", "Lh...Prior", "No.Pram", "Alpha", "Sigma.2")], test2[,c("Lh", "Lh...Prior", "No.Pram", "Alpha", "Sigma.2")]))
#gelman.diag(my_list_of_chains)
#
#tree_4_BTraits<-BTRTools::rjpp(rjlog = here("BayesTraits","Phylo_PC_SCORES_tree_4_run_a.txt.VarRates.txt"),
#                               rjtrees = here("BayesTraits","Phylo_PC_SCORES_tree_4_run_a.txt.Output.trees"),
#                               tree = treelist2[[4]])
#


tree_4_w_data <- add_rjpp_to_tree(tree_4_BTraits)
threshold <- .15 # the minimum posterior probability you want to plot a symbol for 
p<-ggtree(tree_4_w_data, aes(color = log(meanRate)), size=1)+
  scale_colour_gradientn(colours = color3(100))+
  #theme(legend.position="top")+
  theme(legend.position=c(.32,.83),legend.direction = "horizontal",legend.box.background = element_rect(colour = "black",size =1))+
  scale_size(range = c(1,2))+ 
  labs(title="Tree 4",
       color="log(Rate)")+
  geom_nodepoint(aes(subset=ppRate>threshold, size = ppRate),color='black',fill="grey", shape=24)+
  geom_tiplab(label= sub("_", " ",tree_4_w_data@phylo$tip.label), size=2, color = "black", family = "Arial", fontface="italic")+
  coord_cartesian(xlim = c(-230, 90), #you have to fiddle with these values to get your tip labels to show. the first value should be just before your root time, second value pads out space for tip labels
                  ylim = c(-2, 45), #first value makes room for geo timescale, second value is vertical space and should be a few more than your number of tips
                  expand = FALSE) +
  scale_x_continuous(breaks=-periods$max_age[c(1:5)], labels=periods$max_age[c(1:5)]) + 
  theme(panel.grid.major.x = element_line(colour="grey", size=0.5))#should also be modified based on your time scale limits
p <- revts(p) 
ptree4 <- gggeo_scale(p, neg = FALSE, center_end_labels = TRUE, height = unit(1, "line"), size=3)


p.combined <- grid.arrange(ptree1, ptree2, ptree3, ptree4, nrow =2)

ggsave(filename = here("Fig_Output", "rate_trees.pdf"),
       plot = p.combined, device = cairo_pdf,
       width = 18, height = 22, units = "cm")



ggsave(filename = here("Fig_Output", "rate_trees1.pdf"),
       plot = ptree1, device = cairo_pdf,
       width = 11, height = 11, units = "cm")

# Convergence Test --------------------------------------------------------

#using first 3 PCs as they represent 5% or greater of cumulative variance

#note: there is a problem with the parallel package in R 4.0.1 on mac, but this fixes it
if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) &&
Sys.info()["sysname"] == "Darwin" && getRversion() == "4.0.1") {
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
}


traits_vector <- rep("nostate", 43)
names(traits_vector)<-rownames(pca_results$x)

longsnout_forms<-c("Tomistoma_schlegelii","Pelagosaurus_typus","Gavialis_gangeticus","Cricosaurus","Mecistops_cataphractus","Pholidosaurus_sp")
traits_vector[which(names(traits_vector)%in%longsnout_forms)]<-"long"

RRates <- RRphylo(treelist2[[1]], pca_results$x[,c(1:3)],clus=.5)#clus .5 will use half the available cores on your system
conv.test <- search.conv(tree=treelist2[[1]], y=pca_results$x,  state = traits_vector, foldername =  here("Convergence_Tests"))

conv.test

# EMMLIi Test of Modularity --------------------------------------------------------


coords.2d <- two.d.array(proc_aligned_rightside)

contrasts_tree1 <- apply(coords.2d, 2, FUN=pic, phy=treelist2[[1]]) %>% arrayspecs(., p=dim(proc_aligned_rightside)[1], k = 3)
contrasts_tree2 <- apply(coords.2d, 2, FUN=pic, phy=treelist2[[2]]) %>% arrayspecs(., p=dim(proc_aligned_rightside)[1], k = 3)
contrasts_tree3 <- apply(coords.2d, 2, FUN=pic, phy=treelist2[[3]]) %>% arrayspecs(., p=dim(proc_aligned_rightside)[1], k = 3)
contrasts_tree4 <- apply(coords.2d, 2, FUN=pic, phy=treelist2[[4]]) %>% arrayspecs(., p=dim(proc_aligned_rightside)[1], k = 3)

corrmat_contrasts_t1<-dotcorr(contrasts_tree1)
corrmat_contrasts_t2<-dotcorr(contrasts_tree2)
corrmat_contrasts_t3<-dotcorr(contrasts_tree3)
corrmat_contrasts_t4<-dotcorr(contrasts_tree4)

EMMLI_PIC_tree1 <- EMMLi(corr = corrmat_contrasts_t1, N_sample = 43, mod = as.data.frame(module_defs_right[,-c(1,2,4)])) 
EMMLI_PIC_tree2 <- EMMLi(corr = corrmat_contrasts_t2, N_sample = 43, mod = as.data.frame(module_defs_right[,-c(1,2,4)])) 
EMMLI_PIC_tree3 <- EMMLi(corr = corrmat_contrasts_t3, N_sample = 43, mod = as.data.frame(module_defs_right[,-c(1,2,4)])) 
EMMLI_PIC_tree4 <- EMMLi(corr = corrmat_contrasts_t4, N_sample = 43, mod = as.data.frame(module_defs_right[,-c(1,2,4)])) 

# end ---------------------------------------------------------------------



