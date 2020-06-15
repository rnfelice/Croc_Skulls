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
library(SURGE)

library(here)


load(here("Data","patched.crocs_sept_23_2019.R"))

module_defs <- read_csv(here("Data","croc_module_defs.csv"))
 

 
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
tree1<-read.nexus("~/Google Drive/NHM/crocs/phylogenetic trees/Tree1.nex")
tree2<-read.nexus("~/Google Drive/NHM/crocs/phylogenetic trees/Tree2.nex")
tree3<-read.nexus("~/Google Drive/NHM/crocs/phylogenetic trees/Tree3.nex")
tree4<-read.nexus("~/Google Drive/NHM/crocs/phylogenetic trees/Tree4.nex")

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
 
 
 