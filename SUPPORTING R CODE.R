################################################################################
########################  SUPPORTING DATA AND CODE   ###########################
################################################################################

# RRmorph - a new R package to map phenotypic evolutionary rates and 
# patterns on 3D meshes

#    Melchionna M, Castiglione S, Girardi G, Serio C, Esposito A, 
#    Mondanaro A, Profico A, Sansalone G, Raia P



################################################################################

# PLEASE, set your working directory

## Install packages / Upload libraries
{
  to_install<-c("RRphylo","ape","inflection","phytools","LambertW",
                "Arothron","Morpho", "Rvcg","rgl","ggplot2","devtools")
  
  if(any(!to_install%in%installed.packages()[,1])){
    if("ggtree"%in%to_install[which(!to_install%in%installed.packages()[,1])]) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("ggtree",force = TRUE)
    }
    
    install.packages(to_install[which(!to_install%in%installed.packages()[,1])])
  }
  
  devtools::install_github(repo = "pasraia/RRmorph")
  
  
  ## Open packages
  sapply(c(to_install,"RRmorph"),require,character.only=TRUE)
  
  
} # Packages


################################################################################
##########################   CASE STUDY 1 - RATE.MAP   #########################
################################################################################

# data: data.frame with information on species and group
# mshapeS: 3D surface of the consensus shape
# pairedLM: landmarks and semilandmarks indices
# refmat: landmarks and semilandmarks configurations used for the interpolation
# refsur: 3D meshes used for the interpolation
# slid: not-superimposed landmarks and semilandmarks configurations
# tree: phylogenetic tree

load("case study 1 ratemap.rda")

species<-data$species
group<-data$group


# symmetrize the semilandmarks sets
slid_s<-symmetrize(slid,pairedLM=pairedLM)

########################     Relative Warp Analysis     ########################  
pca<-relWarps(slid_s,alpha = 0)
dimnames(pca$rotated)[[3]]<-dimnames(slid)[[3]]

slid_s<-symmetrize(slid,pairedLM = pairedLM)

### RRphylo

# compute the average rw scores per species

PC_mean<-apply(pca$bescores,2, function(x) tapply(x, species, mean))

colnames(PC_mean)<-paste("PC",1:ncol(PC_mean),sep="")
colnames(pca$bePCs)<-paste("PC",1:ncol(PC_mean),sep="")

# run RRphylo pn pcs explaining 95% of shape variance
RR<-RRphylo(tree,PC_mean[,1:length(which(pca$Var$cumVar<=0.95))])


# reconstruct the surfaces
ROT_mean<-aggregate(vecx(pca$rotated),by=list(species),mean)
spec_mean<-ROT_mean[,1]
ROT_mean<-vecx(ROT_mean[,-1],revert=TRUE,lmdim=3)
dimnames(ROT_mean)[[3]]<-spec_mean

mshapeS$vb<-t(pca$mshape)

recsur<-list()
for (i in 1:dim(ROT_mean)[3]){
  recsur[[i]] <- mshapeS
  recsur[[i]]$vb[1:3,]<-t(ROT_mean[,,i])
}
names(recsur)<-dimnames(ROT_mean)[[3]]


# run rate.map for each species against the root  
RM<-list()
for(i in 1:length(rownames(PC_mean))){
  RM[[i]]<-RRmorph::rate.map(x=c("121",rownames(PC_mean)[i]),
                             mshape_sur = mshapeS,
                             RR,scores=PC_mean[,1:length(which(pca$Var$cumVar<=0.95))],
                             pcs=pca$bePCs[,1:length(which(pca$Var$cumVar<=0.95))],
                             mshape=pca$mshape,
                             k=4,plot = FALSE,from = NULL,to=NULL)
}
names(RM)<-rownames(PC_mean)


####  area differences interpolation 

rate_val<-sapply(RM, function(x) x$differences[[1]])
rate_val<-rate_val[,match(names(recsur),colnames(rate_val))]

y<-names(refmat)

palr<-c("darkred","red","orange","white","lightblue","blue","darkblue")

inter_val<-list()
for (i in 1:length(y)){
  inter_val[[i]]<-RRmorph::interpolMesh(sur = recsur[[y[i]]],
                                        refsur = refsur[[y[i]]],
                                        refmat = refmat[[y[i]]],values = rate_val[,y[i]],
                                        element = "vertices")
}
names(inter_val)<-y


col_mesh<-list()
for (i in 1:length(y)){
  col_mesh[[i]]<-RRmorph::col2mesh(mesh = refsur[[y[i]]],values = inter_val[[y[i]]],
                                   pal = palr,
                                   from =  -1,
                                   to =  1)
  RRmorph::plotLegend(col_mesh[[i]],inter_val[[y[i]]],main=y[i])
}

names(col_mesh)<-y

# visualize the colored mesh
open3d();shade3d(col_mesh$Homo_sapiens,specular="black")




##############################       FIGURE       ##############################       

# Figure 3a
cols<-c("#ff1e99",# AdaOmo
        "#00b9f6",# Catarrhinae
        "#A1FF0A",# Platyrrhinae
        "#ffa300",# Strepsirrhinae
        "blue",   # Hominids
        "#7209b7")# Tarsioidea

df<-data.frame(cols=rep(1,nrow(data)),group)
for(i in 1:length(cols)){
  df[which(group==unique(group)[i]),1]<-cols[i]
}

plot(pca$bescores[,1],pca$bescores[,2],pch=21,bg=df$cols,col="black",cex=1.4,
     lwd=0.5,asp=1,xlab="RW1",ylab="RW2")
#text(pca$bescores[,1],pca$bescores[,2],labels = species,cex = 0.5)


# Figure 4
tipdata<-as.data.frame(RR$rates[match(RR$tree$tip.label,rownames(RR$rates)),,drop=FALSE])
colnames(tipdata)<-"rate"
tipdata[,1]<-LambertW::Gaussianize(tipdata$rate)[,1]

p <- ggtree(RR$tree,color="gray20", size=0.5,layout = "fan")
p +  geom_tippoint(aes(x=x+2,color=tipdata[match(p$data$label,rownames(tipdata)),]),size=3) +
  scale_color_continuous(low='cyan1', high='magenta2')+ 
  geom_cladelabel(node=178, label="", align=TRUE, offset=4, barsize=3, color="#A1FF0A")+
  geom_cladelabel(node=125, label="", align=TRUE, offset=4, barsize=3, color="#00b9f6")+
  geom_cladelabel(node=166, label="", align=TRUE, offset=7, barsize=3, color="blue")+
  geom_cladelabel(node=238, label="", align=TRUE, offset=4, barsize=3, color="#ff1e99")+
  geom_cladelabel(node=210, label="", align=TRUE, offset=4, barsize=3, color="#ffa300")+
  geom_cladelabel(node=207, label="", align=TRUE, offset=4, barsize=3, color="#7209b7")+
  
  theme_tree(bgcolor="white" ,legend.position="none")


# empty the environment before moving to the next case study
rm(list=ls())

################################################################################
#########################   CASE STUDY 2 - CONV.MAP   ##########################
################################################################################

# data: data.frame with information on species and group
# mshapeS: 3D surface of the consensus shape
# pairedLM: landmarks and semilandmarks indices
# refmat: landmarks and semilandmarks configurations to be used for the interpolation
# refsur: 3D meshes to be used for the interpolation
# slid: not-superimposed landmarks and semilandmarks configurations
# tree: phylogenetic tree

load("case study 2 convmap.rda")

species<-data$species
group<-data$group


# symmetrize the semilandmarks sets
slid_s<-symmetrize(slid,pairedLM=pairedLM)

########################     Relative Warp Analysis     ########################  
pca<-relWarps(slid_s,alpha = 0)
dimnames(pca$rotated)[[3]]<-dimnames(slid)[[3]]

# Figure 3b
cols<-c("#ff1e99",# AdaOmo
        "#00b9f6",# Catarrhinae
        "#A1FF0A",# Platyrrhinae
        "#ffa300",# Strepsirrhinae
        "blue",   # Hominids
        "#7209b7")# Tarsioidea

df<-data.frame(cols=rep(1,nrow(data)),group)
for(i in 1:length(cols)){
  df[which(group==unique(group)[i]),1]<-cols[i]
}

plot(pca$bescores[,1],pca$bescores[,2],pch=21,bg=df$cols,col="black",cex=1.4,
     lwd=0.5,asp=1,xlab="RW1",ylab="RW2")
#text(pca$bescores[,1],pca$bescores[,2],labels = species,cex = 0.5)



### RRphylo

# compute the average rw scores per species

PC_mean<-apply(pca$bescores,2, function(x) tapply(x, species, mean))

colnames(PC_mean)<-paste("PC",1:ncol(PC_mean),sep="")
colnames(pca$bePCs)<-paste("PC",1:ncol(PC_mean),sep="")

# run RRphylo and search.conv on pcs explaining 95% of shape variance 
RR<-RRphylo(tree,PC_mean[,1:length(which(pca$Var$cumVar<=0.95))])
sc<-search.conv(RR,y = PC_mean[,1:length(which(pca$Var$cumVar<=0.95))],min.dim = 5)


# convergent clades
tips(RR$tree,174)
tips(RR$tree,162)


#############################       CONV.MAP       #############################

mshapeS$vb<-t(pca$mshape)

x1<-tips(RR$tree,174)
x2<-tips(RR$tree,162)

CM<-RRmorph::conv.map(x1=x1,x2=x2,scores = PC_mean, pcs = pca$bePCs,mshape = pca$mshape,
                      mshape_sur = mshapeS,
                      refsur = refsur,
                      refmat = refmat,
                      col = "#00b9f6")


# empty the environment before moving to the next case study
rm(list=ls())


################################################################################
################    CASE STUDY 3 - CONV.MAP CRANIUM/ENDOCAST    ################
################################################################################

# data: data.frame with information on species and group
# mshapeS_e: 3D surface of the endocasts consensus shape
# mshapeS_c: 3D surface of the skulls consensus shape
# pairedLM_e: landmarks and semilandmarks indices of endocasts configuration
# pairedLM_c: landmarks and semilandmarks indices of skulls configuration
# refmat: landmarks and semilandmarks configurations to be used for the interpolation
# refsur: 3D meshes to be used for the interpolation
# slid_e: not-superimposed landmarks and semilandmarks endocast configurations
# slid_c: not-superimposed landmarks and semilandmarks skull configurations
# tree: phylogenetic tree

load("case study 3 convmap skull and endocast.rda")

species<-data$species
group<-data$group

# symmetrize the semilandmarks sets
slid_e_s<-symmetrize(slid_e,pairedLM=pairedLM_e)
slid_c_s<-symmetrize(slid_c,pairedLM=pairedLM_c)

########################     Relative Warp Analysis     ########################  

pca_e<-relWarps(slid_e_s,alpha = 0)
dimnames(pca_e$rotated)[[3]]<-dimnames(slid_e)[[3]]

pca_c<-relWarps(slid_c_s,alpha = 0)
dimnames(pca_c$rotated)[[3]]<-dimnames(slid_c)[[3]]


### RRphylo

# compute the average rw scores per species

PC_mean_e<-apply(pca_e$bescores,2, function(x) tapply(x, species, mean))
PC_mean_c<-apply(pca_c$bescores,2, function(x) tapply(x, species, mean))

colnames(PC_mean_e)<-paste("PC",1:ncol(PC_mean_e),sep="")
colnames(pca_e$bePCs)<-paste("PC",1:ncol(PC_mean_e),sep="")

colnames(PC_mean_c)<-paste("PC",1:ncol(PC_mean_c),sep="")
colnames(pca_c$bePCs)<-paste("PC",1:ncol(PC_mean_c),sep="")


# run RRphylo and search.conv on pcs explaining 95% of shape variance
RR_e<-RRphylo(tree,PC_mean_e[,1:length(which(pca_e$Var$cumVar<=0.95))])
sc_e<-search.conv(RR_e,y = PC_mean_e[,1:length(which(pca_e$Var$cumVar<=0.95))],min.dim = 5)

RR_c<-RRphylo(tree,PC_mean_c[,1:length(which(pca_c$Var$cumVar<=0.95))])
sc_c<-search.conv(RR_c,y = PC_mean_c[,1:length(which(pca_c$Var$cumVar<=0.95))],min.dim = 5)

# convergent clades
tips(RR_c$tree,158)
tips(RR_c$tree,167)

#############################       CONV.MAP       #############################

mshapeS_e$vb<-t(pca_e$mshape)
mshapeS_c$vb<-t(pca_c$mshape)

x1=tips(RR_c$tree,158)
x2=tips(RR_c$tree,167)

CM_e<-RRmorph::conv.map(x1,x2,scores = PC_mean_e, pcs = pca_e$bePCs,mshape = pca_e$mshape,
                        mshape_sur = mshapeS_e, 
                        refsur = refsur_e,
                        refmat = refmat_e,
                        col = "#00b9f6")

CM_c<-RRmorph::conv.map(x1,x2,scores = PC_mean_c, pcs = pca_c$bePCs,mshape = pca_c$mshape,
                        mshape_sur = mshapeS_c,
                        refsur = refsur_c,
                        refmat = refmat_c,
                        col = "#00b9f6")

