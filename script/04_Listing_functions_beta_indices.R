#_______________________________________________________________________________
# Title              : 04_Listing_functions_beta_indices.r
# Date               : 02/02/2025
# Object             : Script to create dataset of values from functions that 
#                      calculate beta diversity indices
# Authors            : Jean-Yves Dias
# R version          : 4.5.0
# Github link        : 
#_______________________________________________________________________________

#_______________________________ Packages_______________________________________####

# Function to install packages if not present and/or load them
loadpackages <- function(packages){
  for (pkg in packages){
    if (!requireNamespace(pkg,quietly = T)){
      install.packages(pkg)
      library(pkg,character.only = T)
    }else{
      library(pkg,character.only = T)
    }
  }
}

packages_needed <- c("readr","dplyr","tidyr","stringr")

loadpackages(packages_needed)

#_______________________________________________________________________________####
#________________________________Loading data___________________________________####
library(vegan)
# We choose the mite data from vegan as the dataset for abundance/count data
data("mite")
detach(package:vegan)
# It is possible to work on sipoo data for presence/absence and varespec for abundance/not integer data

# Relative abundances data to allow some functions working
mite_relat <- mite/rowSums(mite)

# Distance matrix to allow some functions working
mite_dist <- as.matrix(dist(t(mite),method = "euclidean",diag = T,upper = T))
rownames(mite_dist) <- rownames(t(mite))
colnames(mite_dist) <- rownames(t(mite))

#___________________________ Alpha diversity indices ___________________________####

# Anderberg ####
PERMANOVA_1 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_1[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 5,transformation = 1)$D[1,1]
}

PERMANOVA_2 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_2[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 11,transformation = 1)$D[1,1]
}

# Aitchison ####
vegan_3 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_3[i] <- vegan::vegdist(mite_beta,method = "robust.aitchison")
}
pctax_3 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_3[i] <- pctax::mat_dist(t((mite_beta)),method = "robust.aitchison")
}
ecodive_3 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_3[i] <- ecodive::aitchison(mite_beta) 
}

# Bhattacharyya ####
ecodive_4 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_4[i] <- ecodive::bhattacharyya(mite_beta,rescale = F)
}
proxy_4 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_4[i] <- proxy::dist(mite_beta,method = "Bhjattacharyya")
}
philentropy_4 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_4[i] <- philentropy::bhattacharyya(as.numeric(mite[i,]),as.numeric(mite[i+1,])) #-4.223818
}
EnvNJ_4 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_4[i] <- EnvNJ::metrics(t(mite_beta),method = "bhattacharyya")
}
Rfast_4 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_4[i] <- Rfast::Dist(mite_beta,method = "bhattacharyya")[1,2]
}

# Binomial ####
vegan_5 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_5[i] <- vegan::vegdist(mite_beta,method = "binomial")
}
pctax_5 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_5[i] <- pctax::mat_dist(t((mite_beta)),method = "binomial")
}
abdiv_5 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_5[i] <- abdiv::binomial_deviance(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
coda.base_5 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  coda.base_5[i] <- coda.base::dist(mite_beta,"binary")
}
NST_5 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_5[i] <- NST::beta.g(mite_beta,dist.method = "binomial")
}

# Binomial co-occurrence assessment ####
tabula_6 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  tabula_6[i] <- tabula::index_binomial(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

# Brainerd-Robinson ####
tabula_7 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  tabula_7[i] <- tabula::index_brainerd(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
brsim_7 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  brsim_7[i] <- brsim::brsim(mite_beta)$BR.similarity.matrix[1,2]
}

# Braun-Blanquet ####
fossil_8 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  fossil_8[i] <- fossil::braun.blanquet(mite[i,],mite[i+1,])
}
proxy_8 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_8[i] <- proxy::dist(mite_beta,method = "Braun-Blanquet")
}

# Bray-Curtis ####
ecodist_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodist_9[i] <- ecodist::distance(mite_beta,method = "bray-curtis") 
}
vegan_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_9[i] <- vegan::vegdist(mite_beta,method = "bray")
}
provenance_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  provenance_9[i] <- provenance::bray.diss(mite[i,],mite[i+1,])
}
otuSummary_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  otuSummary_9[i] <- otuSummary::calc_bc(mite_beta)
}
provenance_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  provenance_9[i] <- provenance::bray.diss(mite[i,],mite[i+1,])
}
analogue_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  analogue_9[i] <- analogue::distance(mite[i,],mite[i+1,],method = "bray")
}
pctax_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_9[i] <- pctax::mat_dist(t((mite_beta)),method = "bray")
}
bioregion_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  bioregion_9[i] <- bioregion::dissimilarity(as.matrix(mite_beta),metric = "Bray")$Bray
}
fAssets_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fAssets_9[i] <- fAssets::braycurtisDist(t(mite_beta))
}
ecodive_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_9[i] <- ecodive::bray(mite_beta,rescale = F)
}
abdiv_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_9[i] <- abdiv::bray_curtis(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
tabula_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  tabula_9[i] <- tabula::index_bray(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
benthos_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  benthos_9[i] <- benthos::bray_curtis(mite[i,],mite[i+1,])
}
chemodiv_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_relat[c(i,i+1),]
  chemodiv_9[i] <- chemodiv::sampDis(mite_beta,type = "BrayCurtis")$BrayCurtis[1,2]
}
wiqid_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_9[i] <- wiqid::distBrayCurtis(mite[i,],mite[i+1,])
}
fossil_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  fossil_9[i] <- fossil::bray.curtis(mite[i,],mite[i+1,])
}
proxy_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_relat[c(i,i+1),]
  proxy_9[i] <- proxy::dist(mite_beta,method = "Bray")
}
labdsv_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_relat[c(i,i+1),]
  labdsv_9[i] <- labdsv::dsvdis(mite_beta, index = "bray/curtis")
}
PERMANOVA_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_relat[c(i,i+1),]
  PERMANOVA_9[i] <- PERMANOVA::DistContinuous(mite_beta,coef = 8)$D[1,2]
}
ClusterR_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_relat[c(i,i+1),]
  ClusterR_9[i] <- ClusterR::distance_matrix(mite_beta,method = "braycurtis")[2,1]
}
provenance_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  provenance_9[i] <- provenance::bray.diss(mite[i,],mite[i+1,])
}
NST_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_relat[c(i,i+1),]
  NST_9[i] <- NST::beta.g(mite_beta,dist.method = "bray")
}
adespatial_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_relat[c(i,i+1),]
  adespatial_9[i] <- adespatial::beta.div(mite_beta,method = "percentdiff",save.D = T)$D
}

# Manque les Sorensen quanti ##

# Canberra ####

dist(x11,method = "canberra") #27.50188

ecodive::canberra(x12,rescale = F) #20.62641

mgc::mgc.distance(x12,method = "canberra") #27.50188

LearnClust::canberradistance(as.numeric(x1),as.numeric(x2)) #1.068272

ChemoSpecUtils::rowDist(as.matrix(x12),method = "canberra") #27.50188

fAssets::canberraDist(t(x12)) #27.50188

abdiv::canberra(as.numeric(x1),as.numeric(x1)) #20.62641

diverse::dis_entities(t(x12),method = "Canberra",category_row = T)[1,2] #27.50188

proxy::dist(x12,method = "Canberra") #27.50188

proxyC::dist(as.matrix(x12),method = "canberra") ##20.62641

BoutrosLab.plotting.general::dist(x12,method = "canberra") #27.50188

phm::canberra(as.numeric(x1),as.numeric(x2)) #20.62641

amap::Dist(x12,method = "canberra") #27.50188

PERMANOVA::DistContinuous(x12,coef = 7)$D #NA

Mercator::binaryDistance(t(x12),metric = "canberra") #13.75094

dynutils::calculate_distance(x12,method = "canberra") #20.62641

philentropy::canberra(as.numeric(x1),as.numeric(x2)) #20.62641

EnvNJ::metrics(t(x12),method = "canberra") #27.50188

Rfast::Dist(x12,method = "canberra") #NA

Rlof::distmc(x12,method = "canberra") #27.50188

ClusterR::distance_matrix(x12,method = "canberra") #31.62625

rdist::rdist(x12,metric = "canberra") #NA

coda.base::dist(x12,"canberra") #27.50188

fda.usc::metric.dist(x12,method = "canberra") #27.50188

flexclust::dist2(x1,x2,method = "canberra") #27.50188

# Canberra 2 ####

NST::beta.g(x12,dist.method = "canberra") #0.6250428 correspond pas
adespatial::beta.div(x12,method = "canberra",save.D = T)$D #0.6250428
vegan::vegdist(x12,method = "canberra") #0.6250428  
# Canberra index is divided by the number of variables in vegdist, but not in dist. So these differ by a constant multiplier, and the alternative in vegdist is in range (0,1).

pctax::mat_dist(t((x12)),method = "canberra") #0.6250428


# Autres
bioregion::dissimilarity(as.matrix(x12),metric = "Brayturn") #0.5293722





