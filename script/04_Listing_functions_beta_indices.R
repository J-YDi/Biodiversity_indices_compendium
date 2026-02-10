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
convert_to_presence_absence <- function(data) {
  # Vérifie si data est un data.frame ou une matrice
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("L'entrée doit être un data.frame ou une matrice.")
  }
  
  # Convertir en présence (1) / absence (0)
  presence_absence <- apply(data, c(1, 2), function(x) {
    if (is.na(x) || x == 0 || x == "") {
      return(0)
    } else {
      return(1)
    }
  })
  
  # Retourner sous forme de data.frame
  return(as.data.frame(presence_absence))
}
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

mite_pa <- convert_to_presence_absence(mite)

#___________________________ Beta diversity indices ___________________________####

# Anderberg ####
PERMANOVA_1_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  PERMANOVA_1_P[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 5,transformation = 1)$D[1,2]
}

# Put all the values in a single dataframe
P1_P <- ls(pattern = "_1_P$")
P1_P <- mget(P1_P)
P1_P <- as.data.frame(P1_P)

# Aitchison ####
vegan_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_2_A[i] <- vegan::vegdist(mite_beta,method = "robust.aitchison")
}
# vegan_2_P <- rep(NA,69)
# for (i in 1:(nrow(mite)-1)){
#   mite_beta <- mite[c(i,i+1),]
#   vegan_2_P[i] <- vegan::vegdist(mite_beta,method = "robust.aitchison",binary = T)
# }
pctax_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_2_A[i] <- pctax::mat_dist(t((mite_beta)),method = "robust.aitchison")
}
ecodive_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_2_A[i] <- ecodive::aitchison(mite_beta) 
}

# Put all the values in a single dataframe
P2_A <- ls(pattern = "_2_A$")
P2_A <- mget(P2_A)
P2_A <- as.data.frame(P2_A)

# Bhattacharyya ####
ecodive_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_3_A[i] <- ecodive::bhattacharyya(mite_beta,rescale = F)
}
proxy_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_3_A[i] <- proxy::dist(mite_beta,method = "Bhjattacharyya")
}
philentropy_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_3_A[i] <- philentropy::bhattacharyya(as.numeric(mite[i,]),as.numeric(mite[i+1,])) #-4.223818
}
EnvNJ_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_3_A[i] <- EnvNJ::metrics(t(mite_beta),method = "bhattacharyya")
}
Rfast_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_3_A[i] <- Rfast::Dist(mite_beta,method = "bhattacharyya")[1,2]
}

# Put all the values in a single dataframe
P3_A <- ls(pattern = "_3_A$")
P3_A <- mget(P3_A)
P3_A <- as.data.frame(P3_A)


# Binomial ####
vegan_4_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_4_A[i] <- vegan::vegdist(mite_beta,method = "binomial")
}
vegan_4_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_4_P[i] <- vegan::vegdist(mite_beta,method = "binomial",binary = T)
}
pctax_4_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_4_A[i] <- pctax::mat_dist(t((mite_beta)),method = "binomial")
}
abdiv_4_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_4_A[i] <- abdiv::binomial_deviance(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
coda.base_4_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  coda.base_4_A[i] <- coda.base::dist(mite_beta,"binary")
}
NST_4_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_4_A[i] <- NST::beta.g(mite_beta,dist.method = "binomial")
}
NST_4_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_4_P[i] <- NST::beta.g(mite_beta,dist.method = "binomial",transform.method = "pa")
}

# Put all the values in a single dataframe
P4_A <- ls(pattern = "_4_A$")
P4_A <- mget(P4_A)
P4_A <- as.data.frame(P4_A)

# Put all the values in a single dataframe
P4_P <- ls(pattern = "_4_P$")
P4_P <- mget(P4_P)
P4_P <- as.data.frame(P4_P)

# Binomial co-occurrence assessment ####
tabula_5_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  tabula_5_A[i] <- tabula::index_binomial(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P5_A <- ls(pattern = "_5_A$")
P5_A <- mget(P5_A)
P5_A <- as.data.frame(P5_A)

# Brainerd-Robinson ####
tabula_6_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  tabula_6_A[i] <- tabula::index_brainerd(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
brsim_6_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  brsim_6_A[i] <- brsim::brsim(mite_beta)$BR.similarity.matrix[1,2]
}

P6_A <- ls(pattern = "_6_A$")
P6_A <- mget(P6_A)
P6_A <- as.data.frame(P6_A)

# Braun-Blanquet ####
fossil_7_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  fossil_7_P[i] <- fossil::braun.blanquet(mite[i,],mite[i+1,])
}
proxy_7_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_7_P[i] <- 1-proxy::dist(mite_beta,method = "Braun-Blanquet")
}

P7_P <- ls(pattern = "_7_P$")
P7_P <- mget(P7_P)
P7_P <- as.data.frame(P7_P)

# Bray-Curtis ####
ecodist_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodist_8_A[i] <- ecodist::distance(mite_beta,method = "bray-curtis") 
}
vegan_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_8_A[i] <- vegan::vegdist(mite_beta,method = "bray")
}
provenance_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  provenance_8_A[i] <- provenance::bray.diss(mite[i,],mite[i+1,])
}
otuSummary_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  otuSummary_8_A[i] <- otuSummary::calc_bc(mite_beta)
}
provenance_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  provenance_8_A[i] <- provenance::bray.diss(mite[i,],mite[i+1,])
}
analogue_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  analogue_8_A[i] <- analogue::distance(mite[i,],mite[i+1,],method = "bray")
}
pctax_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_8_A[i] <- pctax::mat_dist(t((mite_beta)),method = "bray")
}
bioregion_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  bioregion_8_A[i] <- bioregion::dissimilarity(as.matrix(mite_beta),metric = "Bray")$Bray
}
fAssets_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fAssets_8_A[i] <- fAssets::braycurtisDist(t(mite_beta))
}
ecodive_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_8_A[i] <- ecodive::bray(mite_beta,rescale = F)
}
abdiv_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_8_A[i] <- abdiv::bray_curtis(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
tabula_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  tabula_8_A[i] <- 1-tabula::index_bray(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
benthos_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  benthos_8_A[i] <- benthos::bray_curtis(mite[i,],mite[i+1,])
}
chemodiv_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  chemodiv_8_A[i] <- chemodiv::sampDis(mite_beta,type = "BrayCurtis")$BrayCurtis[1,2]
}
wiqid_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_8_A[i] <- wiqid::distBrayCurtis(mite[i,],mite[i+1,])
}
fossil_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  fossil_8_A[i] <- 1-fossil::bray.curtis(mite[i,],mite[i+1,])
}
proxy_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_8_A[i] <- proxy::dist(mite_beta,method = "Bray")
}
labdsv_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  labdsv_8_A[i] <- labdsv::dsvdis(mite_beta, index = "bray/curtis")
}
PERMANOVA_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_8_A[i] <- PERMANOVA::DistContinuous(mite_beta,coef = 8)$D[1,2]
}
ClusterR_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_8_A[i] <- ClusterR::distance_matrix(mite_beta,method = "braycurtis")[2,1]
}
provenance_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  provenance_8_A[i] <- provenance::bray.diss(mite[i,],mite[i+1,])
}
NST_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_8_A[i] <- NST::beta.g(mite_beta,dist.method = "bray")
}
adespatial_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_8_A[i] <- adespatial::beta.div(mite_beta,method = "percentdiff",save.D = T)$D
}
BAT_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_8_A[i] <- BAT::beta(mite_beta,func = "sorensen",abund = T)$Btotal
}
EnvNJ_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_8_A[i] <- EnvNJ::metrics(t(mite_beta),method = "sorensen")
}
philentropy_8_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  philentropy_8_A[i] <- philentropy::sorensen(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P8_A <- ls(pattern = "_8_A$")
P8_A <- mget(P8_A)
P8_A <- as.data.frame(P8_A)

# Canberra ####
stats_9_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  stats_9_A[i] <- dist(mite_beta,method = "canberra")
}

mgc_9_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  mgc_9_A[i] <- mgc::mgc.distance(mite_beta,method = "canberra")[1,2]
}
LearnClust_9_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  LearnClust_9_A[i] <- LearnClust::canberradistance(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
ChemoSpecUtils_9_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ChemoSpecUtils_9_A[i] <- ChemoSpecUtils::rowDist(as.matrix(mite_beta),method = "canberra")
}
fAssets_9_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fAssets_9_A[i] <- fAssets::canberraDist(t(mite_beta))
}

diverse_9_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  diverse_9_A[i] <- diverse::dis_entities(t(mite_beta),method = "Canberra",category_row = T)[1,2]
}
proxy_9_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_9_A[i] <- proxy::dist(mite_beta,method = "Canberra")
}

BoutrosLab.plotting.general_9_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BoutrosLab.plotting.general_9_A[i] <- BoutrosLab.plotting.general::dist(mite_beta,method = "canberra")
}

amap_9_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  amap_9_A[i] <- amap::Dist(mite_beta,method = "canberra")
}
Mercator_9_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Mercator_9_A[i] <- Mercator::binaryDistance(t(mite_beta),metric = "canberra") 
}

EnvNJ_9_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_9_A[i] <- EnvNJ::metrics(t(mite_beta),method = "canberra")[1,2]
}
Rlof_9_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rlof_9_A[i] <- Rlof::distmc(mite_beta,method = "canberra")
}
ClusterR_9_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_9_A[i] <- ClusterR::distance_matrix(mite_beta,method = "canberra")[2,1]
}
coda.base_9_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  coda.base_9_A[i] <- coda.base::dist(mite_beta,"canberra")
}
fda.usc_9_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fda.usc_9_A[i] <- fda.usc::metric.dist(mite_beta,method = "canberra")[1,2]
}
flexclust_9_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  flexclust_9_A[i] <- flexclust::dist2(mite[i,],mite[i+1,],method = "canberra")
}

P9_A <- ls(pattern = "_9_A$")
P9_A <- mget(P9_A)
P9_A <- as.data.frame(P9_A)


# Canberra 2 ####
NST_10_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_10_A[i] <- NST::beta.g(mite_beta,dist.method = "canberra")
}
adespatial_10_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_10_A[i] <- adespatial::beta.div(mite_beta,method = "canberra",save.D = T)$D
}
vegan_10_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_10_A[i] <- vegan::vegdist(mite_beta,method = "canberra")
}
pctax_10_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_10_A[i] <- pctax::mat_dist(t((mite_beta)),method = "canberra")
}

P10_A <- ls(pattern = "_10_A$")
P10_A <- mget(P10_A)
P10_A <- as.data.frame(P10_A)

# Canberra 3 ####
abdiv_11_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_11_A[i] <- abdiv::canberra(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
dynutils_11_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  dynutils_11_A[i] <- dynutils::calculate_distance(mite_beta,method = "canberra")[1,2]
}
ecodive_11_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_11_A[i] <- ecodive::canberra(mite_beta,rescale = F)
}
philentropy_11_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_11_A[i] <- philentropy::canberra(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
phm_11_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  phm_11_A[i] <- phm::canberra(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
proxyC_11_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxyC_11_A[i] <- proxyC::dist(as.matrix(mite_beta),method = "canberra")[1,2]
}

P11_A <- ls(pattern = "_11_A$")
P11_A <- mget(P11_A)
P11_A <- as.data.frame(P11_A)

# Cao ####
abdiv_12_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_12_A[i] <- abdiv::cy_dissimilarity(as.numeric(mite[i,]),as.numeric(mite[i+1,]),base = 10) # base can be modified
}

P12_A <- ls(pattern = "_12_A$")
P12_A <- mget(P12_A)
P12_A <- as.data.frame(P12_A)

# Cao 2 ####
vegan_13_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_13_A[i] <- vegan::vegdist(mite_beta,method = "cao")
}
pctax_13_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_13_A[i] <- pctax::mat_dist(t((mite_beta)),method = "cao")
}
NST_13_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_13_A[i] <- NST::beta.g(mite_beta,dist.method = "cao")
}

P13_A <- ls(pattern = "_13_A$")
P13_A <- mget(P13_A)
P13_A <- as.data.frame(P13_A)


# Chao-Jaccard ####
NST_14_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_14_A[i] <- NST::beta.g(mite_beta,dist.method = "chao")
}
CommEcol_14_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  CommEcol_14_A[i] <- CommEcol::dis.chao(mite_beta,index = "jaccard",version = "rare") #0.05012367
}
adespatial_14_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_14_A[i] <- adespatial::beta.div(mite_beta,method = "ab.jaccard",save.D = T,sqrt.D = F)$D
}
pctax_14_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_14_A[i] <- pctax::mat_dist(t((mite_beta)),method = "chao") #0.05012367
}
vegan_14_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_14_A[i] <- vegan::vegdist(mite_beta,method = "chao")
}
wiqid_14_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_14_A[i] <- wiqid::distChaoJaccCorr(mite[i,],mite[i+1,])
}

P14_A <- ls(pattern = "_14_A$")
P14_A <- mget(P14_A)
P14_A <- as.data.frame(P14_A)


# Chao-Ochiai ####
vegan_15_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_15_A[i] <- vegan::chaodist(mite_beta,method = "1 - sqrt(U*V)")
}

P15_A <- ls(pattern = "_15_A$")
P15_A <- mget(P15_A)
P15_A <- as.data.frame(P15_A)

# Chao-Sorensen ####
vegan_16_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_16_A[i] <- vegan::chaodist(mite_beta,method = "1 - 2*U*V/(U+V)")
}
adespatial_16_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_16_A[i] <- adespatial::beta.div(mite_beta,method = "ab.sorensen",save.D = T,sqrt.D = F)$D
}
CommEcol_16_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  CommEcol_16_A[i] <- CommEcol::dis.chao(mite_beta,index = "sorensen",version = "rare")
}
wiqid_16_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  wiqid_16_A[i] <- wiqid::distChaoSorCorr(mite[i,],mite[i+1,])
}

P16_A <- ls(pattern = "_16_A$")
P16_A <- mget(P16_A)
P16_A <- as.data.frame(P16_A)

# Chao-Sorensen ####

wiqid_17_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  wiqid_17_A[i] <- wiqid::distChaoSorNaive(mite[i,],mite[i+1,])
} 
NST_17_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_17_A[i] <- NST::beta.g(mite_beta,dist.method = "chao.sorensen")
} 
CommEcol_17_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  CommEcol_17_A[i] <- CommEcol::dis.chao(mite_beta,index = "sorensen",version = "probability") #0.02642424
} 

P17_A <- ls(pattern = "_17_A$")
P17_A <- mget(P17_A)
P17_A <- as.data.frame(P17_A)

# Chebyshev ####
ecodive_18_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_18_A[i] <- ecodive::chebyshev(mite_beta,rescale = F)
}
abdiv_18_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_18_A[i] <- abdiv::chebyshev(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}
philentropy_18_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_18_A[i] <- philentropy::chebyshev(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
EnvNJ_18_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_18_A[i] <- EnvNJ::metrics(t(mite_beta),method = "chebyshev")[1,2]
}
LearnClust_18_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  LearnClust_18_A[i] <- LearnClust::chebyshevDistance(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
comparator_18_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  comparator_18_A[i] <- comparator::Chebyshev()(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}
SBCK_18_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  SBCK_18_A[i] <- print(SBCK::chebyshev(as.matrix(mite[i,]),as.matrix(mite[i+1,]))) 
}
Trading_18_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  Trading_18_A[i] <- Trading::Chebyshev_distance(mite[i,],mite[i+1,])
}
beadplexr_18_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  beadplexr_18_A[i] <- beadplexr::dist_chebyshev(mite_beta)
}
ClusterR_18_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_18_A[i] <- ClusterR::distance_matrix(mite_beta,method = "chebyshev")[2,1]
}
rdist_18_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  rdist_18_A[i] <- rdist::rdist(mite_beta,metric = "chebyshev")
}

P18_A <- ls(pattern = "_18_A$")
P18_A <- mget(P18_A)
P18_A <- as.data.frame(P18_A)

# Chi2 ####
vegan_19_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_19_A[i] <- vegan::vegdist(mite_beta,method = "chisq")
}
pctax_19_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_19_A[i] <- pctax::mat_dist(t((mite_beta)),method = "chisq")
}
svs_19_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  svs_19_A[i] <- svs::dist_chisquare(as.matrix(mite_beta))
}
analogue_19_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  analogue_19_A[i] <- analogue::distance(mite[i,],mite[i+1,],method = "chi.square")
}
adespatial_19_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_19_A[i] <- adespatial::beta.div(mite_beta,method = "chisquare",save.D = T)$D
}
SNFtool_19_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_19_A[i] <- adespatial::beta.div(mite_beta,method = "chisquare",save.D = T)$D
}
proxy_19_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_19_A[i] <- proxy::dist(mite_beta,method = "Chi-squared")
}
spaa_19_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  spaa_19_A[i] <- spaa::sp.pair(t(as.matrix(mite_beta)))$chisq
}

P19_A <- ls(pattern = "_19_A$")
P19_A <- mget(P19_A)
P19_A <- as.data.frame(P19_A)

# Squared chi square ####
ecodive_20_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_20_A[i] <- ecodive::squared_chisq(mite_beta,rescale = F)
}
dynutils_20_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  dynutils_20_A[i] <- dynutils::calculate_distance(mite_beta,method = "chisquared")[1,2] #0
}
EnvNJ_20_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_20_A[i] <- EnvNJ::metrics(t(mite_beta),method = "squared_chi")
}
analogue_20_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  analogue_20_A[i] <- analogue::distance(mite[i,],mite[i+1,],method = "SQchi.square")
}
philentropy_20_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_20_A[i] <- philentropy::squared_chi_sq(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P20_A <- ls(pattern = "_20_A$")
P20_A <- mget(P20_A)
P20_A <- as.data.frame(P20_A) 

# Probabilistic Symmetric chi square distance ####
ecodive_21_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_21_A[i] <- ecodive::psym_chisq(mite_beta,rescale = F)
}

P21_A <- ls(pattern = "_21_A$")
P21_A <- mget(P21_A)
P21_A <- as.data.frame(P21_A) 

# Chord ####
vegan_22_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_22_A[i] <- vegan::vegdist(mite_beta,method = "chord")
}
pctax_22_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_22_A[i] <- pctax::mat_dist(t((mite_beta)),method = "chord")
}
analogue_22_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  analogue_22_A[i] <- analogue::distance(mite[i,],mite[i+1,],method = "chord")
}
ecodive_22_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_22_A[i] <- ecodive::chord(mite_beta)
}
adespatial_22_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_22_A[i] <- adespatial::beta.div(mite_beta,method = "chord",save.D = T)$D
}
proxy_22_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_22_A[i] <- proxy::dist(mite_beta,method = "Chord")
}
abdiv_22_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_22_A[i] <- abdiv::chord(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}
wiqid_22_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_22_A[i] <- wiqid::distChord(mite[i,],mite[i+1,])
}

P22_A <- ls(pattern = "_22_A$")
P22_A <- mget(P22_A)
P22_A <- as.data.frame(P22_A) 

# Squared chord distance ####
ecodive_23_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_23_A[i] <- ecodive::squared_chord(mite_beta,rescale = F)
}
analogue_23_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  analogue_23_A[i] <- analogue::distance(mite[i,],mite[i+1,],method = "SQchord") 
}
philentropy_23_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_23_A[i] <- philentropy::squared_chord(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
EnvNJ_23_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_23_A[i] <- EnvNJ::metrics(t(mite_beta),method = "squared_chord")
}

P23_A <- ls(pattern = "_23_A$")
P23_A <- mget(P23_A)
P23_A <- as.data.frame(P23_A) 

# Log chord distance ####
adespatial_24_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_24_A[i] <- adespatial::beta.div(mite_beta,method = "log.chord",save.D = T)$D
}

P24_A <- ls(pattern = "_24_A$")
P24_A <- mget(P24_A)
P24_A <- as.data.frame(P24_A) 

# Clark ####
vegan_25_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_25_A[i] <- vegan::vegdist(mite_beta,method = "clark") 
}
pctax_25_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_25_A[i] <- pctax::mat_dist(t((mite_beta)),method = "clark")
}
adespatial_25_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_25_A[i] <- adespatial::beta.div(mite_beta,method = "divergence",save.D = T)$D
}
abdiv_25_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_25_A[i] <- abdiv::clark_coefficient_of_divergence(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P25_A <- ls(pattern = "_25_A$")
P25_A <- mget(P25_A)
P25_A <- as.data.frame(P25_A) 

# Clark 2 ####
ecodive_26_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_26_A[i] <- ecodive::clark(mite_beta,rescale = F)
}
philentropy_26_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_26_A[i] <- philentropy::clark_sq(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P26_A <- ls(pattern = "_26_A$")
P26_A <- mget(P26_A)
P26_A <- as.data.frame(P26_A) 

# Cosine ####
vegan_27_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_27_A[i] <- vegan::designdist(mite_beta,method = "1-J/sqrt(A*B)",terms = "quadratic")
}
dynutils_27_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  dynutils_27_A[i] <- dynutils::calculate_distance(mite_beta,method = "cosine")[1,2]
}
SemNeT_27_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  SemNeT_27_A[i] <- SemNeT::similarity(t(mite_beta),method = "cosine")[1,2]
}
EnvNJ_27_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_27_A[i] <- EnvNJ::metrics(t(mite_beta),method = "cosine")[1,2]
}
ChemoSpecUtils_27_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ChemoSpecUtils_27_A[i] <- ChemoSpecUtils::rowDist(as.matrix(mite_beta),method = "cosine")
}
Rfast_27_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_27_A[i] <- Rfast::Dist(mite_beta,method = "cosine")[1,2]
}
ClusterR_27_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_27_A[i] <- ClusterR::distance_matrix(mite_beta,method = "cosine")[2,1]
}
svs_27_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  svs_27_A[i] <- svs::dist_cosine(as.matrix(mite_beta))
}
proxyC_27_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  proxyC_27_A[i] <- proxyC::simil(as.matrix(mite[i,]),as.matrix(mite[i+1,]),method = "cosine")
}
abdiv_27_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_27_A[i] <- abdiv::cosine_distance(mite[i,],mite[i+1,])
}
amap_27_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  amap_27_A[i] <- amap::Dist(mite_beta,method = "pearson")
}

P27_A <- ls(pattern = "_27_A$")
P27_A <- mget(P27_A)
P27_A <- as.data.frame(P27_A) 

# Absolute Pearson ####
ChemoSpecUtils_28_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ChemoSpecUtils_28_A[i] <- ChemoSpecUtils::rowDist(as.matrix(mite_beta),method = "abspearson")
}
amap_28_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  amap_28_A[i] <- amap::Dist(mite_beta,method = "abspearson")
}

P28_A <- ls(pattern = "_28_A$")
P28_A <- mget(P28_A)
P28_A <- as.data.frame(P28_A)

# Divergence ####
ecodive_29_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_29_A[i] <- ecodive::divergence(mite_beta,rescale = F) 
}
proxy_29_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_29_A[i] <- proxy::dist(mite_beta,method = "divergence")
}
philentropy_29_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  philentropy_29_A[i] <- philentropy::divergence_sq(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P29_A <- ls(pattern = "_29_A$")
P29_A <- mget(P29_A)
P29_A <- as.data.frame(P29_A)

# Euclidean distance ####
ecodist_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodist_30_A[i] <- ecodist::distance(mite_beta,method = "euclidean")
}
stats_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  stats_30_A[i] <- dist(mite_beta,method = "euclidean")
}
vegan_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_30_A[i] <- vegan::vegdist(mite_beta,method = "euclidean")
}
cluster_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  cluster_30_A[i] <- cluster::daisy(mite_beta,metric = "euclidean")
}
mgc_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  mgc_30_A[i] <- mgc::mgc.distance(mite_beta)[1,2]
}
distances_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  distances_30_A[i] <- distances::distances(mite_beta)[1,2]
}
ldt_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ldt_30_A[i] <- ldt::s.distance(t(mite_beta),distance = "euclidean")
}
ChemoSpecUtils_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ChemoSpecUtils_30_A[i] <- ChemoSpecUtils::rowDist(mite_beta,method = "euclidean")
}
fAssets_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fAssets_30_A[i] <- fAssets::euclideanDist(t(mite_beta))
}
ecodive_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_30_A[i] <- ecodive::euclidean(mite_beta,rescale = F)
}
proxyC_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxyC_30_A[i] <- proxyC::dist(as.matrix(mite_beta),method = "euclidean")[2,1]
}
adespatial_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_30_A[i] <- adespatial::beta.div(mite_beta,method = "euclidean",save.D = T)$D
}
diverse_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  diverse_30_A[i] <- diverse::dis_entities(t(mite_beta),method = "euclidean",category_row = T)[1,2]
}
pctax_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_30_A[i] <- pctax::mat_dist(t((mite_beta)),method = "euclidean")
}
proxy_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_30_A[i] <- proxy::dist(mite_beta,method = "Euclidean") 
}
amap_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  amap_30_A[i] <- amap::Dist(mite_beta,method = "euclidean") 
}
MultivariateAnalysis_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  MultivariateAnalysis_30_A[i] <- MultivariateAnalysis::Distancia(mite_beta,Metodo = 1)[1]$Distancia 
}
Mercator_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Mercator_30_A[i] <- Mercator::binaryDistance(t(mite_beta),metric = "euclid")
}
dynutils_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  dynutils_30_A[i] <- dynutils::calculate_distance(mite_beta,method = "euclidean")[1,2]
}
fda.usc_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fda.usc_30_A[i] <- fda.usc::metric.dist(mite_beta,method = "euclidean")[1,2]
}
EnvNJ_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_30_A[i] <- EnvNJ::metrics(t(mite_beta),method = "euclidean")[1,2]
}
Rfast_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_30_A[i] <- Rfast::Dist(mite_beta,method = "euclidean")[1,2]
}
ClusterR_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_30_A[i] <- ClusterR::distance_matrix(mite_beta,method = "euclidean")[2,1]
}
rdist_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  rdist_30_A[i] <- rdist::rdist(mite_beta,metric = "euclidean")
}
coda.base_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  coda.base_30_A[i] <- coda.base::dist(mite_beta,"euclidean")
}
NST_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_30_A[i] <- NST::beta.g(mite_beta,dist.method = "euclidean")
}
Rlof_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rlof_30_A[i] <- Rlof::distmc(mite_beta,method = "euclidean")
}
BoutrosLab.plotting.general_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BoutrosLab.plotting.general_30_A[i] <- BoutrosLab.plotting.general::dist(mite_beta,method = "euclidean")
}
fossil_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  fossil_30_A[i] <- fossil::euclidean(mite[i,],mite[i+1,])
}
abdiv_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_30_A[i] <- abdiv::euclidean(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
philentropy_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_30_A[i] <- philentropy::euclidean(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}
comparator_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  comparator_30_A[i] <- comparator::Euclidean()(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
TSdist_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  TSdist_30_A[i] <- TSdist::LPDistance(as.numeric(mite[i,]),as.numeric(mite[i+1,]),method = "euclidean")
}
flexclust_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  flexclust_30_A[i] <- flexclust::dist2(mite[i,],mite[i+1,],method = "euclidean")
}
qkerntool_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  qkerntool_30_A[i] <- qkerntool::Eucdist(as.matrix(mite[i,]),as.matrix(mite[i+1,]),sEuclidean = T) 
}
codep_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  codep_30_A[i] <- codep::Euclid(mite[i,],mite[i+1,],squared = F) 
}
ptm_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  ptm_30_A[i] <- ptm::pairwise.dist(as.matrix(mite[i,]),as.matrix(mite[i+1,]),squared = F)
}
nnspat_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  nnspat_30_A[i] <- nnspat::euc.dist(mite[i,],mite[i+1,])
}
CEGO_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  CEGO_30_A[i] <- CEGO::distanceRealEuclidean(mite[i,],mite[i+1,])
}
RnavGraphImageData_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  RnavGraphImageData_30_A[i] <- RnavGraphImageData::L2Distance(as.matrix(t(mite[i,])),as.matrix(t(mite[i+1,]))) 
}
statisfactory_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  statisfactory_30_A[i] <- statisfactory::euclid(mite[i,],mite[i+1,])
}
neighbr_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  neighbr_30_A[i] <- neighbr::distance(as.numeric(mite[i,]),as.numeric(mite[i+1,]),measure = "euclidean")
}
analogue_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  analogue_30_A[i] <- analogue::distance(mite[i,],mite[i+1,],method = "euclidean")
}
LearnClust_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  LearnClust_30_A[i] <- LearnClust::edistance(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
comparator_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  comparator_30_A[i] <- comparator::Euclidean()(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
hmsr_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  hmsr_30_A[i] <- hmsr::euclidean_distance(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P30_A <- ls(pattern = "_30_A$")
P30_A <- mget(P30_A)
P30_A <- as.data.frame(P30_A)

# Anderson's modified Euclidean distance ####
NST_31_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_31_A[i] <- NST::beta.g(mite_beta,dist.method = "mEuclidean")
}

P31_A <- ls(pattern = "_31_A$")
P31_A <- mget(P31_A)
P31_A <- as.data.frame(P31_A)

# Average euclidean distance ####
abdiv_32_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_32_A[i] <- abdiv::rms_distance(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P32_A <- ls(pattern = "_32_A$")
P32_A <- mget(P32_A)
P32_A <- as.data.frame(P32_A)

# Squared euclidean distance ####
ecodive_33_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_33_A[i] <- ecodive::squared_euclidean(mite_beta,rescale = F)
}
MultivariateAnalysis_33_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  MultivariateAnalysis_33_A[i] <- MultivariateAnalysis::Distancia(mite_beta,Metodo = 3)[1]$Distancia
}
neighbr_33_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  neighbr_33_A[i] <- neighbr::distance(as.numeric(mite[i,]),as.numeric(mite[i+1,]),measure = "squared_euclidean")
}
NMFN_33_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  NMFN_33_A[i] <- NMFN::distance2(mite[i,],mite[i+1,])
}
ptm_33_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  ptm_33_A[i] <- ptm::pairwise.dist(as.matrix(mite[i,]),as.matrix(mite[i+1,]),squared = T)
}
codep_33_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  codep_33_A[i] <- codep::Euclid(mite[i,],mite[i+1,],squared = T)
}
M2SMJF_33_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  M2SMJF_33_A[i] <- M2SMJF::dist2eu(mite[i,],mite[i+1,])
}
philentropy_33_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_33_A[i] <- philentropy::squared_euclidean(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
analogue_33_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  analogue_33_A[i] <- analogue::distance(mite[i,],mite[i+1,],method = "SQeuclidean")
}
qkerntool_33_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  qkerntool_33_A[i] <- qkerntool::Eucdist(as.matrix(mite[i,]),as.matrix(mite[i+1,]),sEuclidean = F)
}
laGP_33_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  laGP_33_A[i] <- laGP::distance(mite[i,],mite[i+1,])
}
Rfast_33_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_33_A[i] <- Rfast::Dist(mite_beta,method = "euclidean",square = T)[1,2]
}

P33_A <- ls(pattern = "_33_A$")
P33_A <- mget(P33_A)
P33_A <- as.data.frame(P33_A)

# Gower ####
cluster_34_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  cluster_34_A[i] <- cluster::daisy(mite_beta,metric = "gower")
}
StatMatch_34_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  StatMatch_34_A[i] <- StatMatch::gower.dist(mite[i,],mite[i+1,])#0.7045455
}
vegan_34_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_34_A[i] <- vegan::vegdist(mite_beta,method = "gower")
}
pctax_34_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_34_A[i] <- pctax::mat_dist(t((mite_beta)),method = "gower")
}
ecodive_34_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_34_A[i] <- ecodive::gower(mite_beta,rescale = F)
}
diverse_34_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  diverse_34_A[i] <- diverse::dis_entities(t(mite_beta),method = "Gower",category_row = T)[1,2]
}
shipunov_34_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  shipunov_34_A[i] <- shipunov::Gower.dist(mite[i,],mite[i+1,])
}
proxy_34_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_34_A[i] <- proxy::dist(mite_beta,method = "Gower") 
}
NST_34_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_34_A[i] <- NST::beta.g(mite_beta,dist.method = "gower")
}

P34_A <- ls(pattern = "_34_A$")
P34_A <- mget(P34_A)
P34_A <- as.data.frame(P34_A)

# Gower 2 ####
NST_35_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_35_A[i] <- NST::beta.g(mite_beta,dist.method = "altGower")
}
vegan_35_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_35_A[i] <- vegan::vegdist(mite_beta,method = "altGower")
}
pctax_35_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_35_A[i] <- pctax::mat_dist(t(mite_beta),method = "altGower")
}
philentropy_35_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  philentropy_35_A[i] <- philentropy::gower(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
Rfast_35_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_35_A[i] <- Rfast::Dist(mite_beta,method = "gower")[2,1]
}
gower_35_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  gower_35_A[i] <- gower::gower_dist(mite[i,],mite[i+1,])
}

P35_A <- ls(pattern = "_35_A$")
P35_A <- mget(P35_A)
P35_A <- as.data.frame(P35_A)

# Gower 3 ####
ecodist_36_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodist_36_A[i] <- ecodist::distance(mite_beta,method = "modgower10")
}

P36_A <- ls(pattern = "_36_A$")
P36_A <- mget(P36_A)
P36_A <- as.data.frame(P36_A)

# Anderson's modified Gower distance ####
NST_37_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_37_A[i] <- NST::beta.g(mite_beta,dist.method = "mGower")
}

P37_A <- ls(pattern = "_37_A$")
P37_A <- mget(P37_A)
P37_A <- as.data.frame(P37_A)

# Hamman coefficient ####
ade4_38_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_38_P[i] <- ade4::dist.binary(mite_beta,method = 6)
}
proxy_38_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_38_P[i] <- proxy::dist(mite_beta,method = "Hamman")
}
MultivariateAnalysis_38_P <- rep(NA,69) # A besoin de PA pour etre correct
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  MultivariateAnalysis_38_P[i] <- MultivariateAnalysis::Distancia(mite_beta,Metodo = 19)[1]$Distancia
}
MultBiplotR_38_P <- rep(NA,69) # A besoin de PA pour fonctionner
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  MultBiplotR_38_P[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 9)[1,2]
}
PERMANOVA_38_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_38_P[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 9,transformation = 1)$D[1,2]
}
proxyC_38_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  proxyC_38_P[i] <- 1-proxyC::simil(as.matrix(mite[i,]),as.matrix(mite[i+1,]),method = "hamann") #0.5 ok
}

P38_P <- ls(pattern = "_38_P$")
P38_P <- mget(P38_P)
P38_P <- as.data.frame(P38_P)

# Hamming distance ####
ecodive_39_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_39_P[i] <- ecodive::hamming(mite_beta) 
}
abdiv_39_A <- rep(NA,69) 
for (i in 1:(nrow(mite)-1)){
  abdiv_39_A[i] <- abdiv::hamming(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
Rankcluster_39_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  Rankcluster_39_A[i] <- Rankcluster::distHamming(mite[i,],mite[i+1,])
}
Mercator_39_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Mercator_39_A[i] <- Mercator::binaryDistance(t(mite_beta),metric = "hamming")
}
pegas_39_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pegas_39_A[i] <- pegas::dist.hamming(mite_beta)
}
proxyC_39_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxyC_39_A[i] <- proxyC::dist(as.matrix(mite_beta),method = "hamming")[2,1]
}
bingat_39_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  bingat_39_P[i] <- bingat::calcDistance(mite[i,],mite[i+1,]) 
}
genMCMCDiag_39_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  genMCMCDiag_39_A[i] <- genMCMCDiag::hammingDist(mite[i,],mite[i+1,]) #11 
}

P39_P <- ls(pattern = "_39_P$")
P39_P <- mget(P39_P)
P39_P <- as.data.frame(P39_P)

P39_A <- ls(pattern = "_39_A$")
P39_A <- mget(P39_A)
P39_A <- as.data.frame(P39_A)

# Hamming distance 2 ####
CEGO_40_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  CEGO_40_A[i] <- CEGO::distanceNumericHamming(mite[i,],mite[i+1,])
}
EnvNJ_40_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_40_A[i] <- EnvNJ::metrics(t(mite_beta),method = "hamming")[2,1]
}
ClusterR_40_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_40_A[i] <- ClusterR::distance_matrix(mite_beta,method = "hamming")[2,1]
}
rdist_40_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  rdist_40_A[i] <- rdist::rdist(mite_beta,metric = "hamming")
}
CEGO_40_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  CEGO_40_A[i] <- CEGO::distanceNumericHamming(mite[i,],mite[i+1,])
}
EnsCat_40_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnsCat_40_A[i] <- EnsCat::hammingD(mite_beta)[2,1]
}

P40_A <- ls(pattern = "_40_A$")
P40_A <- mget(P40_A)
P40_A <- as.data.frame(P40_A)

# Hellinger ####
vegan_41_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_41_A[i] <- vegan::vegdist(mite_beta,method = "hellinger")
}
pctax_41_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_41_A[i] <- pctax::mat_dist(t((mite_beta)),method = "hellinger")
}
ecodive_41_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_41_A[i] <- ecodive::hellinger(mite_beta,rescale = T)
}
adespatial_41_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_41_A[i] <- adespatial::beta.div(mite_beta,method = "hellinger",save.D = T)$D
}
proxy_41_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_41_A[i] <- proxy::dist(mite_beta,method = "Hellinger")
}
Rfast_41_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_41_A[i] <- Rfast::Dist(mite_beta,method = "hellinger")[2,1]
}
abdiv_41_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_41_A[i] <- abdiv::hellinger(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P41_A <- ls(pattern = "_41_A$")
P41_A <- mget(P41_A)
P41_A <- as.data.frame(P41_A)

# Jaccard (Abondance) ####
vegan_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_42_A[i] <- vegan::vegdist(mite_beta,method = "jaccard",binary = F) # 0.6936661
}
adespatial_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_42_A[i] <- adespatial::beta.div(mite_beta,method = "ruzicka",save.D = T)$D
}
abdiv_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_42_A[i] <- abdiv::ruzicka(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}
labdsv_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  labdsv_42_A[i] <- labdsv::dsvdis(mite_beta, index = "ruzicka")
}
philentropy_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_42_A[i] <- abdiv::ruzicka(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}
stylo_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  stylo_42_A[i] <- stylo::dist.minmax(mite_beta)
}
adiv_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_42_A[i] <- adiv::distMS(mite_beta)
}
ecodive_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_42_A[i] <- ecodive::soergel(mite_beta,rescale = F) 
}
proxy_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_42_A[i] <- proxy::dist(mite_beta,method = "Soergel")
}
PERMANOVA_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_42_A[i] <- PERMANOVA::DistContinuous(mite_beta,coef = 9)$D[2,1]
}
EnvNJ_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_42_A[i] <- EnvNJ::metrics(t(mite_beta),method = "soergel")
}
Rfast_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_42_A[i] <- Rfast::Dist(mite_beta,method = "soergel")[2,1]
}
NST_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_42_A[i] <- NST::beta.g(mite_beta,dist.method = "jaccard")
}
picante_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  picante_42_A[i] <- 1-picante::species.dist(t(mite_beta),metric = "jaccard")
}
BoutrosLab.plotting.general_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BoutrosLab.plotting.general_42_A[i] <- BoutrosLab.plotting.general::dist(mite_beta,method = "jaccard") 
}
statisfactory_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  statisfactory_42_A[i] <- 1-statisfactory::fuzzyJaccard(mite[i,],mite[i+1,]) 
}
geocmeans_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  geocmeans_42_A[i] <- 1-geocmeans::calc_jaccard_idx(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
adespatial_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_42_A[i] <- adespatial::beta.div.comp(mite_beta,coef = "J",quant = T)$D 
}
pctax_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_42_A[i] <- pctax::mat_dist(t((mite_beta)),method = "jaccard")
}
BAT_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_42_A[i] <- BAT::beta(mite_beta,func = "jaccard",abund = T)$Btotal
}
prabclus_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  prabclus_42_A[i] <- prabclus::jaccard(t(mite_beta)) 
}
Mercator_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Mercator_42_A[i] <- Mercator::binaryDistance(t(mite_beta),metric = "jaccard") #2.266299
}
rdist_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  rdist_42_A[i] <- rdist::rdist(mite_beta,metric = "jaccard") 
}
ConNEct_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  ConNEct_42_A[i] <- ConNEcT::funCorrJacc(as.numeric(mite[i,]),as.numeric(mite[i+1,]))$value #0.9834216
}

P42_A <- ls(pattern = "_42_A$")
P42_A <- mget(P42_A)
P42_A <- as.data.frame(P42_A)

# Jaccard P/A ####
labdsv_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  labdsv_42_P[i] <- labdsv::dsvdis(mite_beta, index = "steinhaus")
}
ecodist_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodist_42_P[i] <- ecodist::distance(mite_beta,method = "jaccard")
}
stats_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  stats_42_P[i] <- dist(mite_beta,method = "binary")
}
vegan_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_42_P[i] <- vegan::vegdist(mite_beta,method = "jaccard",binary = T)
}
neighbr_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  neighbr_42_P[i] <- 1-neighbr::similarity(mite_pa[i,],mite_pa[i+1,],measure = "jaccard") 
}
ChemoSpecUtils_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ChemoSpecUtils_42_P[i] <- ChemoSpecUtils::rowDist(as.matrix(mite_beta),method = "binary") 
}
mgc_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  mgc_42_P[i] <- mgc::mgc.distance(mite_beta,method = "binary")[1,2] 
}
BoutrosLab.plotting.general_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BoutrosLab.plotting.general_42_P[i] <- BoutrosLab.plotting.general::dist(mite_beta,method = "binary")
}
bioregion_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  bioregion_42_P[i] <- bioregion::dissimilarity(as.matrix(mite_beta),metric = "Jaccard")$Jaccard
}
fAssets_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fAssets_42_P[i] <- fAssets::jaccardDist(t(mite_beta))
}
vegan_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_42_P[i] <- 1-vegan::betadiver(mite_beta,"j")
}
ecodive_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_42_P[i] <- ecodive::jaccard(mite_beta)
}
betapart_42_P <- rep(NA,69) # Besoin de PA pour fonctionner
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  betapart_42_P[i] <- betapart::beta.pair(mite_beta,index.family = "jaccard")$beta.jac 
}
adespatial_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_42_P[i] <- adespatial::beta.div.comp(mite_beta,coef = "J",quant = F)$D 
}
tabula_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  tabula_42_P[i] <- 1-tabula::similarity(mite_beta,method = "jaccard")
}
adiv_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_42_P[i] <- adiv::betastatjac(mite_beta)[1]
}
diverse_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  diverse_42_P[i] <-diverse::dis_entities(t(mite_beta),method = "Jaccard",category_row = T)[1,2]
}
BAT_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_42_P[i] <- BAT::beta(mite_beta,func = "jaccard",abund = F)$Btotal
}
proxy_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_42_P[i] <- proxy::dist(mite_beta,method = "Jaccard")
}
PERMANOVA_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_42_P[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 3,transformation = 1)$D[1,2]
}
ClusterR_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_42_P[i] <- ClusterR::distance_matrix(mite_beta,method = "jaccard_coefficient")[2,1]
}
flexclust_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  flexclust_42_P[i] <- flexclust::dist2(mite[i,],mite[i+1,],method = "binary")
}
iTOP_42_P <- rep(NA,69) # Besoin de PA pour fonctionner
for (i in 1:(nrow(mite)-1)){
  iTOP_42_P[i] <- 1-iTOP::jaccard(mite_pa[i,],mite_pa[i+1,])
}
proxyC_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  proxyC_42_P[i] <- 1-proxyC::simil(as.matrix(mite[i,]),as.matrix(mite[i+1,]),method = "jaccard")
}
abdiv_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_42_P[i] <- abdiv::jaccard(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
wiqid_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_42_P[i] <- wiqid::distJaccard(mite[i,],mite[i+1,])
}
fossil_42_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  fossil_42_P[i] <- 1-fossil::jaccard(mite[i,],mite[i+1,])
}

P42_P <- ls(pattern = "_42_P$")
P42_P <- mget(P42_P)
P42_P <- as.data.frame(P42_P)

# Squared-root Jaccard ####
adespatial_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_43_P[i] <- print(adespatial::dist.ldc(mite_beta,method = "jaccard"))
}
MultBiplotR_43_P <- rep(NA,69) # A besoin de PA pour fonctionner
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  MultBiplotR_43_P[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 3)[1,2] 
}
ade4_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_43_P[i] <- ade4::dist.binary(mite_beta,method = 1) 
}

P43_P <- ls(pattern = "_43_P$")
P43_P <- mget(P43_P)
P43_P <- as.data.frame(P43_P)

# Turnover component of Jaccard dissimilarity defined by Baselga ####
vegan_44_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_44_P[i] <- vegan::nestedbetajac(mite_beta)[1]
}
betapart_44_P <- rep(NA,69) # Besoin de PA pour marcher
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  betapart_44_P[i] <- betapart::beta.pair(mite_beta,index.family = "jaccard")$beta.jtu
}
adespatial_44_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_44_P[i] <- adespatial::beta.div.comp(mite_beta,coef = "BJ")$repl 
}
bioregion_44_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  bioregion_44_P[i] <- bioregion::dissimilarity(as.matrix(mite_beta),metric = "Jaccardturn")$Jaccardturn
}
abdiv_44_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
    abdiv_44_P[i] <- abdiv::jaccard_turnover(as.numeric(mite[i,]),as.numeric(mite[i+1,])) #0.2666667
}

P44_P <- ls(pattern = "_44_P$")
P44_P <- mget(P44_P)
P44_P <- as.data.frame(P44_P)

# Nestedness component of Jaccard dissimilarity defined by Baselga ####
abdiv_45_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_45_P[i] <- abdiv::jaccard_nestedness(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
vegan_45_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_45_P[i] <- vegan::nestedbetajac(mite_beta)[2] 
}
betapart_45_P <- rep(NA,69) # Besoin de PA pour fonctionner
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  betapart_45_P[i] <- betapart::beta.pair(mite_beta,index.family = "jaccard")$beta.jne
}
adespatial_45_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_45_P[i] <- adespatial::beta.div.comp(mite_beta,coef = "BJ")$rich
}

P45_P <- ls(pattern = "_45_P$")
P45_P <- mget(P45_P)
P45_P <- as.data.frame(P45_P)

# Extended Jaccard Similarity ####
adiv_46_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_46_A[i] <- adiv::dsimcom(mite_beta,method = "2",type = "similarity",option = "absolute")[1,2] 
}
diverse_46_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  diverse_46_A[i] <- 1-diverse::dis_entities(t(mite_beta),method = "eJaccard",category_row = T)[1,2] 
}
proxy_46_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_46_A[i] <- 1-proxy::dist(mite_beta,method = "eJaccard")
}
EnvNJ_46_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_46_A[i] <- 1-EnvNJ::metrics(t(mite_beta),method = "jaccard")
}
philentropy_46_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_46_A[i] <- 1-philentropy::jaccard(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
proxyC_46_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  proxyC_46_A[i] <- proxyC::simil(as.matrix(mite[i,]),as.matrix(mite[i+1,]),method = "ejaccard") 
}
adespatial_46_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_46_A[i] <- 1-adespatial::beta.div(mite_beta,method = "wishart",save.D = T)$D 
}
wiqid_46_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_46_A[i] <- 1-wiqid::distSimRatio(mite[i,],mite[i+1,])
}
vegan_46_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_46_A[i] <- 1-vegan::designdist(mite_beta,method = "(A+B-2*J)/(A+B-J)",terms = "quadratic")
}

P46_A <- ls(pattern = "_46_A$")
P46_A <- mget(P46_A)
P46_A <- as.data.frame(P46_A)

# ARRET ICI ####
# Replacement index of Jaccard defined by Podani ####
adespatial_52 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_52[i] <- adespatial::beta.div.comp(mite_beta,coef = "J")$repl
}
BAT_52 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_52[i] <- BAT::beta(mite_beta,func = "jaccard",abund = F)$Brepl
}

# Nestedness component of Jaccard dissimilarity defined by Podani ####
BAT_53 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_53[i] <- BAT::beta(mite_beta,func = "jaccard",abund = F)$Brich
}
adespatial_53 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_53[i] <- adespatial::beta.div.comp(mite_beta,coef = "J")$rich
}

# Nestedness component of Jaccard dissimilarity defined by Podani & Schmera ####
adespatial_54 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_54[i] <- adespatial::beta.div.comp(mite_beta,coef = "N")$rich 
}

# Jaccard richness gain ####
BAT_55 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_55[i] <- BAT::beta(mite_beta,func = "jaccard",abund = F)$Bgain
}
adiv_55 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_55[i] <- adiv::betastatjac(mite_beta)[2]
}

# Jaccard richness loss ####
BAT_56 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_56[i] <- BAT::beta(mite_beta,func = "jaccard",abund = F)$Bloss
}
adiv_56 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_56[i] <- adiv::betastatjac(mite_beta)[3]
}

# Jensen-Shannon distance ####
adiv_57 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_57[i] <- ecodive::jensen(mite_beta,rescale = F)
}

# Jensen-Shannon divergence ####
adiv_58 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_58[i] <- ecodive::jsd(mite_beta,rescale = F)
}
EnvNJ_58 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_58[i] <- EnvNJ::metrics(t(mite_beta),method = "jensen-shannon")
}
Rfast_58 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_58[i] <- Rfast::Dist(mite_beta,method = "jensen_shannon")[1,2] #34.28249
}
Compositional_58 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Compositional_58[i] <- Compositional::divergence(mite_beta,type = "jensen_shannon")[1,2]
}
priorsense_58 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  priorsense_58[i] <- priorsense::cjs_dist(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
philentropy_58 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_58[i] <- philentropy::jensen_difference(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}

# Kulczynski ####
proxy_59 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_59[i] <- proxy::dist(mite_beta,method = "Kulczynski2") 
}
abdiv_59 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_59[i] <- abdiv::kulczynski_second(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}
vegan_59 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_59[i] <- vegan::betadiver(mite_beta,"co")
}
abdiv_59 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_59[i] <- abdiv::kulczynski_second(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
prabclus_59 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  prabclus_59[i] <- prabclus::kulczynski(t(mite_beta))[2,1] 
}
fossil_59 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  fossil_59[i] <- fossil::kulczynski(mite[i,],mite[i+1,]) 
}

# Kulczynski 2
vegan_60 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_60[i] <- vegan::vegdist(mite_beta,method = "kulczynski") 
}
pctax_60 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_60[i] <- pctax::mat_dist(t((mite_beta)),method = "kulczynski") 
}
adespatial_60 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_60[i] <- adespatial::beta.div(mite_beta,method = "kulczynski",save.D = T)$D 
}
NST_60 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_60[i] <- NST::beta.g(mite_beta,dist.method = "kulczynski")
}
prabclus_60 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  prabclus_60[i] <- prabclus::qkulczynski(t(mite_beta))[1,2]
}
proxy_60 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_60[i] <- proxy::dist(mite_beta,method = "Kulczynski1") 
}
abdiv_60 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_60[i] <- abdiv::kulczynski_first(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

# Kulczynski 3 ####
philentropy_61 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_61[i] <- philentropy::kulczynski_d(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
EnvNJ_61 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_61[i] <- EnvNJ::metrics(t(mite_beta),method = "kulczynski")
}
Rfast_61 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_61[i] <- Rfast::Dist(mite_beta,method = "kulczynski")[1,2]
}

# beta l  ####
vegan_62 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_62[i] <- vegan::betadiver(mite_beta,"l") 
}
# beta gl  ####
vegan_63 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_63[i] <- vegan::betadiver(mite_beta,"gl") 
}

# beta sim  ####
vegan_64 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_64[i] <- vegan::betadiver(mite_beta,"sim") 
}

# beta z ####
vegan_65 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_65[i] <- vegan::betadiver(mite_beta,"z") 
}

# Lorentzian distance ####
ecodive_66 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_66[i] <- ecodive::lorentzian(mite_beta,rescale = F) 
}
EnvNJ_66 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_66[i] <- EnvNJ::metrics(t(mite_beta),method = "lorentzian") 
}
philentropy_66 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  philentropy_66[i] <- philentropy::lorentzian(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}

# beta m ####
vegan_67 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_67[i] <- vegan::betadiver(mite_beta,"m")
}

# Mahalanobis ####
vegan_68 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_68[i] <- vegan::vegdist(mite_beta,method = "mahalanobis") 
}
FD_68 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  FD_68[i] <- FD::mahaldis(as.matrix(mite_beta)) 
}
pctax_68 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_68[i] <- pctax::mat_dist(t((mite_beta)),method = "mahalanobis")
}
ClusterR_68 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_68[i] <- ClusterR::distance_matrix(mite_beta,method = "mahalanobis")[2,1] 
}
NST_68 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_68[i] <- NST::beta.g(mite_beta,dist.method = "mahalanobis")
}

# Manhattan ####
ecodist_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodist_69[i] <- ecodist::distance(mite_beta,method = "manhattan")
}
stats_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  stats_69[i] <- dist(mite_beta,method = "manhattan")
}
vegan_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_69[i] <- vegan::vegdist(mite_beta,method = "manhattan")
}
proxyC_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxyC_69[i] <- proxyC::dist(as.matrix(mite_beta),method = "manhattan")[1,2]
}
mgc_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  mgc_69[i] <- mgc::mgc.distance(mite_beta,method = "manhattan")[1,2]
}
cluster_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  cluster_69[i] <- cluster::daisy(mite_beta,metric = "manhattan")
}
ChemoSpecUtils_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ChemoSpecUtils_69[i] <- ChemoSpecUtils::rowDist(as.matrix(mite_beta),method = "manhattan")
}
ldt_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ldt_69[i] <- ldt::s.distance(t(mite_beta),distance = "manhattan")
}
fAssets_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fAssets_69[i] <- fAssets::manhattanDist(t(mite_beta))
}
ecodive_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_69[i] <- ecodive::manhattan(mite_beta,rescale = F)
}
adespatial_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_69[i] <- print(adespatial::dist.ldc(mite_beta,method = "manhattan")) 
}
diverse_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  diverse_69[i] <- diverse::dis_entities(t(mite_beta),method = "Manhattan",category_row = T)[1,2]
}
pctax_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_69[i] <- pctax::mat_dist(t((mite_beta)),method = "manhattan")
}
proxy_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_69[i] <- proxy::dist(mite_beta,method = "Manhattan") 
}
amap_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  amap_69[i] <- amap::Dist(mite_beta,method = "manhattan")
}
Mercator_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Mercator_69[i] <- Mercator::binaryDistance(t(mite_beta),metric = "manhattan")
}
dynutils_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  dynutils_69[i] <- dynutils::calculate_distance(mite_beta,method = "manhattan")[1,2]
}
EnvNJ_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_69[i] <- EnvNJ::metrics(t(mite_beta),method = "manhattan")[1,2]
}
Rfast_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_69[i] <- Rfast::Dist(mite_beta,method = "manhattan")[1,2] #95.06
}
comparator_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  comparator_69[i] <- comparator::Manhattan()(as.numeric(mite[i,]),as.numeric(mite[i+1,])) #95.06
}
ClusterR_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_69[i] <- ClusterR::distance_matrix(mite_beta,method = "manhattan")[2,1]
}
rdist_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  rdist_69[i] <- rdist::rdist(mite_beta,metric = "manhattan")
}
BoutrosLab.plotting.general_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BoutrosLab.plotting.general_69[i] <- BoutrosLab.plotting.general::dist(mite_beta,method = "manhattan")
}
coda.base_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  coda.base_69[i] <- coda.base::dist(mite_beta,"manhattan")
}
NST_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_69[i] <- NST::beta.g(mite_beta,dist.method = "manhattan")
}
Rlof_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rlof_69[i] <- Rlof::distmc(mite_beta,method = "manhattan")
}
fda.usc_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fda.usc_69[i] <- fda.usc::metric.dist(mite_beta,method = "manhattan")[1,2]
}
analogue_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  analogue_69[i] <- analogue::distance(mite[i,],mite[i+1,],method = "manhattan") 
}
hmsr_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  hmsr_69[i] <- hmsr::manhattan_distance(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
abdiv_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_69[i] <- abdiv::manhattan(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
fossil_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  fossil_69[i] <- fossil::manhattan(mite[i,],mite[i+1,])
}
philentropy_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_69[i] <- philentropy::manhattan(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
TSdist_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  TSdist_69[i] <- TSdist::LPDistance(as.numeric(mite[i,]),as.numeric(mite[i+1,]),method = "manhattan") #95.06
}
LearnClust_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  LearnClust_69[i] <- LearnClust::mdistance(as.numeric(mite[i,]),as.numeric(mite[i+1,])) #11.08
}
flexclust_69 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  flexclust_69[i] <- flexclust::dist2(mite[i,],mite[i+1,],method = "manhattan")
}

# Manhattan modified ####
NST_70 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_70[i] <- NST::beta.g(mite_beta,dist.method = "mManhattan")
}

# Matusita ####
ecodive_71 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_71[i] <- ecodive::matusita(mite_beta,rescale = F) #6.514023
}

# Minkowski ####-------------------------------------------------------------
stats_72_1 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  stats_72_1[i] <- dist(mite_beta,method = "minkowski",p=1)
}
stats_72_2 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  stats_72_2[i] <- dist(mite_beta,method = "minkowski",p=2)
}
stats_72_3 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  stats_72_3[i] <- dist(mite_beta,method = "minkowski",p=3)
}
flexclust_72_1 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  flexclust_72_1[i] <- flexclust::dist2(mite[i,],mite[i+1,],method = "minkowski",p=1) #40.37368
}
flexclust_72_2 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  flexclust_72_2[i] <- flexclust::dist2(mite[i,],mite[i+1,],method = "minkowski",p=2)
}
flexclust_72_3 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  flexclust_72_3[i] <- flexclust::dist2(mite[i,],mite[i+1,],method = "minkowski",p=3)
}
mgc_72_2 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  mgc_72_2[i] <- mgc::mgc.distance(mite_beta,method = "minkowski")[1,2]
}
ecodive_72_1 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_72_1[i] <- ecodive::minkowski(mite_beta,rescale = F,power = 1)
}
ecodive_72_2 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_72_2[i] <- ecodive::minkowski(mite_beta,rescale = F,power = 2)
}
ecodive_72_3 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_72_3[i] <- ecodive::minkowski(mite_beta,rescale = F,power = 3)
}
abdiv_72_1 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_72_1[i] <- abdiv::minkowski(as.numeric(mite[i,]),as.numeric(mite[i+1,]),p=1)
}
abdiv_72_2 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_72_2[i] <- abdiv::minkowski(as.numeric(mite[i,]),as.numeric(mite[i+1,]),p=2)
}
abdiv_72_3 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_72_3[i] <- abdiv::minkowski(as.numeric(mite[i,]),as.numeric(mite[i+1,]),p=3)
}
proxy_72_1 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_72_1[i] <- proxy::dist(mite_beta,method = "Minkowski",p=1)
}
proxy_72_2 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_72_2[i] <- proxy::dist(mite_beta,method = "Minkowski",p=2)
}
proxy_72_3 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_72_3[i] <- proxy::dist(mite_beta,method = "Minkowski",p=3)
}
PERMANOVA_72_1 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_72_1[i] <- PERMANOVA::DistContinuous(mite_beta,coef = 4,r=1)$D[1,2]
}
PERMANOVA_72_2 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_72_2[i] <- PERMANOVA::DistContinuous(mite_beta,coef = 4,r=2)$D[1,2]
}
PERMANOVA_72_3 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_72_3[i] <- PERMANOVA::DistContinuous(mite_beta,coef = 4,r=3)$D[1,2]
}
dynutils_72_2 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  dynutils_72_2[i] <- dynutils::calculate_distance(mite_beta,method = "minkowski")[1,2]
}
EnvNJ_72_1 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_72_1[i] <- EnvNJ::metrics(t(mite_beta),method = "minkowski",p=1)[1,2]
}
EnvNJ_72_2 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_72_2[i] <- EnvNJ::metrics(t(mite_beta),method = "minkowski",p=2)[1,2]
}
EnvNJ_72_3 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_72_3[i] <- EnvNJ::metrics(t(mite_beta),method = "minkowski",p=3)[1,2]
}
Rfast_72_1 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_72_1[i] <- Rfast::Dist(mite_beta,method = "minkowski",p = 1)[1,2] #95.06
}
Rfast_72_2 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_72_2[i] <- Rfast::Dist(mite_beta,method = "minkowski",p = 2)[1,2] #95.06
}
Rfast_72_3 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_72_3[i] <- Rfast::Dist(mite_beta,method = "minkowski",p = 3)[1,2] #95.06
}
fAssets_72_2 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fAssets_72_2[i] <- fAssets::minkowskiDist(t(mite_beta)) #40.37368
}
proxyC_72_1 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxyC_72_1[i] <- proxyC::dist(as.matrix(mite_beta),method = "minkowski",p=1)[1,2]
}
proxyC_72_2 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxyC_72_2[i] <- proxyC::dist(as.matrix(mite_beta),method = "minkowski",p=2)[1,2]
}
proxyC_72_3 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxyC_72_3[i] <- proxyC::dist(as.matrix(mite_beta),method = "minkowski",p=3)[1,2]
}
ClusterR_72_1 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_72_1[i] <- ClusterR::distance_matrix(mite_beta,method = "minkowski",minkowski_p = 1)[2,1]
}
ClusterR_72_2 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_72_2[i] <- ClusterR::distance_matrix(mite_beta,method = "minkowski",minkowski_p = 2)[2,1]
}
ClusterR_72_3 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_72_3[i] <- ClusterR::distance_matrix(mite_beta,method = "minkowski",minkowski_p = 3)[2,1]
}
rdist_72_1 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  rdist_72_1[i] <- rdist::rdist(mite_beta,metric = "minkowski",p = 1)
}
rdist_72_2 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  rdist_72_2[i] <- rdist::rdist(mite_beta,metric = "minkowski",p = 2)
}
rdist_72_3 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  rdist_72_3[i] <- rdist::rdist(mite_beta,metric = "minkowski",p = 3)
}
coda.base_72_1 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  coda.base_72_1[i] <- coda.base::dist(mite_beta,"minkowski",p=1)
}
coda.base_72_2 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  coda.base_72_2[i] <- coda.base::dist(mite_beta,"minkowski",p=2) 
}
coda.base_72_3 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  coda.base_72_3[i] <- coda.base::dist(mite_beta,"minkowski",p=3)
}
fda.usc_72_1 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fda.usc_72_1[i] <- fda.usc::metric.dist(mite_beta,method = "minkowski",p=1)[1,2]
}
fda.usc_72_2 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fda.usc_72_2[i] <- fda.usc::metric.dist(mite_beta,method = "minkowski",p=2)[1,2]
}
fda.usc_72_3 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fda.usc_72_3[i] <- fda.usc::metric.dist(mite_beta,method = "minkowski",p=3)[1,2]
}
BoutrosLab.plotting.general_72_1 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BoutrosLab.plotting.general_72_1[i] <- BoutrosLab.plotting.general::dist(mite_beta,method = "minkowski",p=1)
}
BoutrosLab.plotting.general_72_2 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BoutrosLab.plotting.general_72_2[i] <- BoutrosLab.plotting.general::dist(mite_beta,method = "minkowski",p=2)
}
BoutrosLab.plotting.general_72_3 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BoutrosLab.plotting.general_72_3[i] <- BoutrosLab.plotting.general::dist(mite_beta,method = "minkowski",p=3)
}
comparator_72_1 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  comparator_72_1[i] <- comparator::Minkowski(p=1)(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
comparator_72_2 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  comparator_72_2[i] <- comparator::Minkowski(p=2)(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
comparator_72_3 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  comparator_72_3[i] <- comparator::Minkowski(p=3)(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
TSdist_72_1 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  TSdist_72_1[i] <- TSdist::MinkowskiDistance(as.numeric(mite[i,]),as.numeric(mite[i+1,]),p=2)
}
TSdist_72_2 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  TSdist_72_2[i] <- TSdist::MinkowskiDistance(as.numeric(mite[i,]),as.numeric(mite[i+1,]),p=2)
}
TSdist_72_3 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  TSdist_72_3[i] <- TSdist::MinkowskiDistance(as.numeric(mite[i,]),as.numeric(mite[i+1,]),p=2)
}
philentropy_72_1 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_72_1[i] <- philentropy::minkowski(as.numeric(mite[i,]),as.numeric(mite[i+1,]),n=1) 
}
philentropy_72_2 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_72_2[i] <- philentropy::minkowski(as.numeric(mite[i,]),as.numeric(mite[i+1,]),n=2) 
}
philentropy_72_3 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_72_3[i] <- philentropy::minkowski(as.numeric(mite[i,]),as.numeric(mite[i+1,]),n=3)
}

# Morisita ####
abdiv_73 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_73[i] <- abdiv::morisita(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
vegan_73 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_73[i] <- vegan::vegdist(mite_beta,method = "morisita")
}
pctax_73 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_73[i] <- pctax::mat_dist(t((mite_beta)),method = "morisita")
}
ecodive_73 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_73[i] <- ecodive::morisita(mite_beta)
}
NST_73 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_73[i] <- NST::beta.g(mite_beta,dist.method = "morisita")
}

# Morisita-Horn ####
vegan_74 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_74[i] <- vegan::vegdist(mite_beta,method = "horn")
}
NST_74 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_74[i] <- NST::beta.g(mite_beta,dist.method = "horn")
}
pctax_74 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_74[i] <- pctax::mat_dist(t((mite_beta)),method = "horn")
}
ecodive_74 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_74[i] <- ecodive::horn(mite_beta,rescale = F)
}
tabula_74 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  tabula_74[i] <- tabula::similarity(mite_beta,method = "morisita")
}
abdiv_74 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_74[i] <- abdiv::horn_morisita(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
wiqid_74 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_74[i] <- wiqid::distMorisitaHorn(mite[i,],mite[i+1,])
}
fossil_74 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  fossil_74[i] <- fossil::morisita.horn(mite[i,],mite[i+1,])
}

# Motyka ####
pctax_75 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_75[i] <- ecodive::motyka(mite_beta,rescale = F)
}
EnvNJ_75 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_75[i] <- EnvNJ::metrics(t(mite_beta),method = "motyka")
}
Rfast_75 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_75[i] <- Rfast::Dist(mite_beta,method = "motyka")[1,2]
}
philentropy_75 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  philentropy_75[i] <- philentropy::motyka(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

# Mountford ####
vegan_76 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_76[i] <- vegan::vegdist(mite_beta,method = "mountford")
}
pctax_76 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_76[i] <- pctax::mat_dist(t((mite_beta)),method = "mountford")
}
proxy_76 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_76[i] <- proxy::dist(mite_beta,method = "Mountford")
}
NST_76 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_76[i] <- NST::beta.g(mite_beta,dist.method = "mountford")
}

# beta me ####
vegan_77 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_77[i] <- vegan::betadiver(mite_beta,"me")
}

# Ochiai ####
adespatial_78 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_78[i] <- adespatial::beta.div(mite_beta,method = "ochiai",save.D = T,sqrt.D = F)$D
}
ade4_78 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_78[i] <- ade4::dist.binary(mite_beta,method = 7)
}
adiv_78 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_78[i] <- adiv::dsimcom(mite_beta,method = "4",type = "similarity",option = "absolute")[1,2]
}
MultBiplotR_78 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  MultBiplotR_78[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 12)[1,2]
}
spaa_78 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  spaa_78[i] <- spaa::sp.pair(t(as.matrix(mite_beta)))$Ochiai
}

# Otsuka-Ochiai ####
ecodive_79 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_79[i] <- ecodive::ochiai(mite_beta)
}
PERMANOVA_79 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_79[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 12,transformation = 1)$D[1,2]
}
proxy_79 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_79[i] <- proxy::dist(mite_beta,method = "Ochiai")
}
vegan_79 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_79[i] <- vegan::designdist(mite_beta,method = "1-J/sqrt(A*B)",terms = "binary")
}
labdsv_79 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  labdsv_79[i] <- labdsv::dsvdis(mite_beta, index = "ochiai")
}
wiqid_79 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_79[i] <- wiqid::distOchiai(mite[i,],mite[i+1,])
}
fossil_79 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  fossil_79[i] <- fossil::ochiai(mite[i,],mite[i+1,])
}

# Phi-squared ####
ade4_80 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_80[i] <- ade4::dist.binary(mite_beta,method = 9)
}
proxy_80 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_80[i] <- proxy::dist(mite_beta,method = "Phi-squared")
}
MultBiplotR_80 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  MultBiplotR_80[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 14)[1,2]
}
PERMANOVA_80 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_80[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 14,transformation = 3)$D[1,2]
}
philentropy_80 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  philentropy_80[i] <- philentropy::pearson_chi_sq(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

# Preston's coefficient of faunal dissimilarity ####
wiqid_81 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_81[i] <- wiqid::distPreston(mite[i,],mite[i+1,])
}

# Raup ####
vegan_82 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_82[i] <- vegan::vegdist(mite_beta,method = "raup")
}
pctax_82 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_82[i] <- pctax::mat_dist(t((mite_beta)),method = "raup")
}
iCAMP_82 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  iCAMP_82[i] <- iCAMP::RC.pc(mite_beta)$index[1,2]
}
NST_82 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_82[i] <- NST::beta.g(mite_beta,dist.method = "raup")
}

# Rogers & Tanimoto ####
ade4_83 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_83[i] <- ade4::dist.binary(mite_beta,method = 4)
}
abdiv_83 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_83[i] <- abdiv::rogers_tanimoto(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
proxy_83 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_83[i] <- proxy::dist(mite_beta,method = "Tanimoto") 
}
neighbr_83 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  neighbr_83[i] <- neighbr::similarity(mite[i,],mite[i+1,],measure = "tanimoto")
}
wiqid_83 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_83[i] <- wiqid::distRogersTanimoto(mite[i,],mite[i+1,])
}
philentropy_83 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_83[i] <- philentropy::tanimoto(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
MultivariateAnalysis_83 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  MultivariateAnalysis_83[i] <- MultivariateAnalysis::Distancia(mite_beta,Metodo = 15)[1]$Distancia 
}
MultBiplotR_83 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  MultBiplotR_83[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 6)[1,2]
}
PERMANOVA_83 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_83[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 6,transformation = 1)$D[1,2]
}
shipunov_83 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  shipunov_83[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 6,transformation = 1)$D[1,2]
}

# Root mean square ####
EnvNJ_84 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_84[i] <- EnvNJ::metrics(t(mite_beta),method = "avg")
}
abdiv_84 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_84[i] <- abdiv::rms_distance(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

# beta e ####
vegan_85 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_85[i] <- vegan::betadiver(mite_beta,"e")
}

# beta l ####
vegan_86 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_86[i] <- vegan::betadiver(mite_beta,"l")
}

# beta r ####
vegan_87 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_87[i] <- vegan::betadiver(mite_beta,"r")
}
tabula_87 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  tabula_87[i] <- tabula::index_routledge1(as.matrix(mite_beta)) 
}

# beta rlb ####
vegan_88 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_88[i] <- vegan::betadiver(mite_beta,"rlb")
}

# Russel-Rao ####
abdiv_89 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_89[i] <- abdiv::russel_rao(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
proxy_89 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_89[i] <- proxy::dist(mite_beta,method = "Russel")
}
MultivariateAnalysis_89 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  MultivariateAnalysis_89[i] <- MultivariateAnalysis::Distancia(mite_beta,Metodo = 16)[1]$Distancia
}
MultBiplotR_89 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  MultBiplotR_89[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 2)
}
PERMANOVA_89 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_89[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 2,transformation = 1)$D[1,2]
}
Mercator_89 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Mercator_89[i] <- Mercator::binaryDistance(t(mite_beta),metric = "russellRao")
}

# Simple match ####
ade4_90 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_90[i] <- ade4::dist.binary(mite_beta,method = 2)
}
neighbr_90 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  neighbr_90[i] <- neighbr::similarity(mite_beta[i,],mite_beta[i+1,],measure = "simple_matching")
}
proxyC_90 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  proxyC_90[i] <- proxyC::simil(as.matrix(mite[i,]),as.matrix(mite[i+1,]),method = "simple matching")
}
shipunov_90 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  shipunov_90[i] <- shipunov::SM.dist(mite_beta)
}
nomclust_90 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  nomclust_90[i] <- nomclust::sm(mite_beta)
}
proxy_90 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_90[i] <- proxy::dist(mite_beta,method = "simple matching")
}
ClusterR_90 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_90[i] <- ClusterR::distance_matrix(mite_beta,method = "simple_matching_coefficient")[2,1]
}
PERMANOVA_90 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_90[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 4,transformation = 1)$D[1,2]
}
wiqid_90 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_90[i] <- wiqid::distMatching(mite[i,],mite[i+1,])
}
MultBiplotR_90 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  MultBiplotR_90[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 4)
}
arules_90 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  arules_90[i] <- arules::dissimilarity(as.matrix(mite_beta),method = "matching")
}

# Sokal & Sneath ####
abdiv_91 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_91[i] <- abdiv::sokal_sneath(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

# Sokal & Sneath 2 ####
adiv_92 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_92[i] <- adiv::dsimcom(mite_beta,method = "1",type = "dissimilarity",option="absolute")
}

# Sokal & Sneath 3 ####
MultivariateAnalysis_93 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  MultivariateAnalysis_93[i] <- MultivariateAnalysis::Distancia(mite_beta,Metodo = 14)[1]$Distancia
}
MultBiplotR_93 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  MultBiplotR_93[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 8)
}
PERMANOVA_93 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_93[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 8,transformation = 3)$D[2,1]
}
Mercator_93 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Mercator_93[i] <- Mercator::binaryDistance(t(mite_beta),metric = "sokalMichener")
}

# Sokal & Sneath S5 ####
ade4_94 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_94[i] <- ade4::dist.binary(mite_beta,method = 3)
}

# Sokal & Sneath S13 ####
ade4_95 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_95[i] <- ade4::dist.binary(mite_beta,method = 8)
}
MultBiplotR_95 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  MultBiplotR_95[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 13)
}
PERMANOVA_95 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_95[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 13,transformation = 1)$D
}

# Sorensen à faire ####
vegan_96 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_96[i] <- vegan::betadiver(mite_beta,"sor")
}
ecodist_96 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodist_96[i] <- ecodist::distance(mite_beta,method = "sorensen")
}
bioregion_96 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  bioregion_96[i] <- bioregion::dissimilarity(as.matrix(mite_beta),metric = "Sorensen")$Sorensen
}
ecodive_96 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_96[i] <- ecodive::sorensen(mite_beta)
}
betapart_96 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  betapart_96[i] <- betapart::beta.pair(mite_beta,index.family = "sorensen")$beta.sor
}
tabula_96 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  tabula_96[i] <- tabula::similarity(mite_beta,method = "sorensen")
}
adiv_96 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_96[i] <- adiv::betastatsor(mite_beta)[1] #0.2
}
adespatial_96 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_96[i] <- adespatial::beta.div.comp(mite_beta,coef = "S",quant = F)$D #0.2
}
BAT_96 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_96[i] <- BAT::beta(mite_beta,func = "sorensen",abund = F)$Btotal #0.2
}
proxy_96 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_96[i] <- proxy::dist(mite_beta,method = "Dice") #0.2
}
labdsv_96 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  labdsv_96[i] <- labdsv::dsvdis(mite_beta, index = "sorensen")
}
PERMANOVA_96 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_96[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 7,transformation = 1)$D[1,2]
}
fAssets_96 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fAssets_96[i] <- fAssets::sorensenDist(t(mite_beta))
}
abdiv_96 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_96[i] <- abdiv::sorenson(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
proxyC_96 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  proxyC_96[i] <- proxyC::simil(as.matrix(mite[i,]),as.matrix(mite[i+1,]),method = "dice")
}
wiqid_96 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_96[i] <- wiqid::distSorensen(mite[i,],mite[i+1,])
}
fossil_96 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  fossil_96[i] <- fossil::sorenson(mite[i,],mite[i+1,])
}
prabclus_96 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  prabclus_96[i] <- prabclus::dicedist(t(mite_beta))
}
spaa_96 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  spaa_96[i] <- spaa::sp.pair(t(as.matrix(mite_beta)))$Dice
}

# Square-root Sorensen ####
abdiv_97 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_97[i] <- ade4::dist.binary(mite_beta,method = 5)
}
adespatial_97 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_97[i] <- print(adespatial::dist.ldc(mite_beta,method = "sorensen")) #0.4472136
}
MultBiplotR_97 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  MultBiplotR_97[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 7)[1,2]
}

# Turnover component of Sorensen dissimilarity ####
BAT_98_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_98_A[i] <- BAT::beta(mite_beta,func = "sorensen",abund = T)$Brepl #0.5275388
}
adespatial_98_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_98_A[i] <- adespatial::beta.div.comp(mite_beta,coef = "BS",quant = T)$repl
}
adespatial_98_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_98_P[i] <- adespatial::beta.div.comp(mite_beta,coef = "BS")$repl # Baselga 0.1538462
}
BAT_98_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_98_P[i] <- BAT::beta(mite_beta,func = "sorensen",abund = F)$Brepl
}
betapart_98 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  betapart_98[i] <- betapart::beta.pair(mite_beta,index.family = "sorensen")$beta.sim 
}
proxy_98 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_98[i] <- proxy::dist(mite_beta,method = "Simpson")
}
vegan_98 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_98[i] <- vegan::nestedbetasor(mite_beta)[1]
}
fossil_98 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fossil_98[i] <- fossil::simpson(mite[i,],mite[i+1,])
}

# Nestedness-resultant component of Sørensen dissimilarity ####
vegan_99 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_99[i] <- vegan::nestedbetasor(mite_beta)[2]
}
betapart_99 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  betapart_99[i] <- betapart::beta.pair(mite_beta,index.family = "sorensen")$beta.sne
}
abdiv_99 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_99 [i] <- abdiv::sorenson_nestedness(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
adespatial_99_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_99_P [i] <- adespatial::beta.div.comp(mite_beta,coef = "BS")$rich
}
adespatial_99_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_99_A[i] <- adespatial::beta.div.comp(mite_beta,coef = "BS",quant=T)$rich
}

# Sorensen richness gain ####
BAT_100_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_100_A[i] <- BAT::beta(mite_beta,func = "sorensen",abund = T)$Bgain
}
BAT_100_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_100_P[i] <- BAT::beta(mite_beta,func = "sorensen",abund = F)$Bgain
}
adiv_100 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_100[i] <- adiv::betastatsor(mite_beta)[3]
}

# Sorensen richness loss ####
BAT_101_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_101_A[i] <- BAT::beta(mite_beta,func = "sorensen",abund = T)$Bloss
}
BAT_101_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_101_P[i] <- BAT::beta(mite_beta,func = "sorensen",abund = F)$Bloss
}
adiv_101 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_101[i] <- adiv::betastatsor(mite_beta)[2]
}

# Legendre replacement index ####
adespatial_102_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_102_P[i] <- adespatial::beta.div.comp(mite_beta,coef = "J",quant = F)$repl
}
adespatial_102_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_102_A[i] <- adespatial::beta.div.comp(mite_beta,coef = "J",quant = T)$repl
}

# Legendre richness difference index ####
adespatial_103_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_103_P[i] <- adespatial::beta.div.comp(mite_beta,coef = "J",quant = F)$rich
}
adespatial_103_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_103_A[i] <- adespatial::beta.div.comp(mite_beta,coef = "J",quant = T)$rich
}
BAT_103_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_103_P[i] <- BAT::beta(mite_beta,func = "sorensen",abund = F)$Brich
}
BAT_103_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_103_A[i] <- BAT::beta(mite_beta,func = "sorensen",abund = T)$Brich
}

# Extended Sorensen Similarity ####
adiv_104 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_104[i] <- adiv::dsimcom(mite_beta,method = "3",type = "dissimilarity",option = "absolute")
}
diverse_104 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  diverse_104[i] <- diverse::dis_entities(t(mite_beta),method = "eDice",category_row = T)[1,2]
}
proxy_104 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_104[i] <- proxy::dist(mite_beta,method = "eDice")
}
philentropy_104 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_104[i] <- philentropy::dice_dist(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}
proxyC_104 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  proxyC_104[i] <- proxyC::simil(as.matrix(mite[i,]),as.matrix(mite[i+1,]),method = "edice")  #0.5108792 
}



vegan::chaodist(mite_betaint,method = "1 - 2*U*V/(U+V)") # 0.1172414 faux
abdiv::jaccard_turnover(as.numeric(x1),as.numeric(x2)) #0.2666667

# Species profile distance ####
adespatial_105 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_105[i] <- adespatial::beta.div(mite_beta,method = "profiles",save.D = T)$D #0.4505909
}
ade4_105 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_105[i] <- ade4::disc(as.data.frame(t(mite_beta)))
}
ClusterR_105 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_105[i] <- ClusterR::distance_matrix(mite_beta,method = "Rao_coefficient")[2,1]
}

# Topsoe ####
ecodive_106 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_106[i] <- ecodive::topsoe(mite_beta,rescale = F)
}
philentropy_106 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  philentropy_106[i] <- philentropy::topsoe(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

# Wave Hedges distance ####
ecodive_107 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_107[i] <- ecodive::wave_hedges(mite_beta,rescale = F) #22.84192
}
proxy_107 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_107[i] <- proxy::dist(mite_beta,method = "Wave")
}
PERMANOVA_107 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_107[i] <- PERMANOVA::DistContinuous(mite_beta,coef = 10)$D[2,1]
}
EnvNJ_107 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_107[i] <- EnvNJ::metrics(t(mite_beta),method = "wavehedges")
}
philentropy_107 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  philentropy_107[i] <- philentropy::wave_hedges(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
Rfast_107 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_107[i] <- Rfast::Dist(mite_beta,method = "wave_hedges")[2,1]
}

# Whittaker's index of association ####
wiqid_108 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  wiqid_108[i] <- wiqid::distWhittaker(mite[i,],mite[i+1,])
}
proxy_108 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_108[i] <- proxy::dist(mite_beta,method = "Whittaker")
}
adespatial_108 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_108[i] <- adespatial::beta.div(mite_beta,method = "whittaker",save.D = T)$D
}

# beta w ####
vegan_109 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_109[i] <- vegan::betadiver(mite_beta,"w")
}
tabula_109 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  tabula_109[i] <- tabula::index_whittaker(as.matrix(mite_beta))
}
# beta wb ####
vegan_110 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_110[i] <- vegan::betadiver(mite_beta,"wb")
}

# beta -3 ####
vegan_111 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_111[i] <- vegan::betadiver(mite_beta,"-3")
}

# beta t ####
vegan_112 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_112[i] <- vegan::betadiver(mite_beta,"t")
}

tabula_XXX <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  tabula_XXX[i] <- tabula::index_wilson(as.matrix(mite_beta)) #ne sait pas a quoi ca correspond
}

# Yule ####
abdiv_113 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_113[i] <- abdiv::yule_dissimilarity(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
proxy_113 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_113[i] <- proxy::dist(mite_beta,method = "Yule")
}
MultBiplotR_113 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  MultBiplotR_113[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 15)[1,2]
}
PERMANOVA_113 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_113[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 15,transformation = 1)$D[1,2]
}

# Yule 2 ####
proxy_114<- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_114[i] <- proxy::dist(mite_beta,method = "Yule2")
}


# beta g ####
vegan_X <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_X[i] <- vegan::betadiver(mite_beta,"g")
}

# beta c ####
vegan_X <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_X[i] <- vegan::betadiver(mite_beta,"c")
}
tabula_X <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  tabula_X[i] <- tabula::index_cody(as.matrix(mite_beta)) 
}

# Sans nom ####
vegan_115 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_115[i] <- vegan::betadiver(mite_beta,"19")
}

# Jeffreys ####
EnvNJ_116 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_116[i] <- EnvNJ::metrics(t(mite_beta),method = "jeffreys")
}

# Mean character difference ####
abdiv_XX <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_XX[i] <- abdiv::mean_character_difference(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

# Modified mean character difference ####
adespatial_XX <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_XX[i] <- print(adespatial::dist.ldc(mite_beta,method = "modmeanchardiff")) 
}
abdiv_XX <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_XX[i] <- abdiv::modified_mean_character_difference(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

# S2 coeff Gower & Legendre ####
ade4_XX <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_XX[i] <- ade4::dist.binary(mite_beta,method = 10)
}

# Harmonic mean ####
Rfast_XX <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_XX[i] <- Rfast::Dist(mite_beta,method = "harmonic_mean")
}

# beta cc ####
vegan_28 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_28[i] <- vegan::betadiver(mite_beta,"cc")
}

# beta -1 ####
vegan_43 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_43[i] <- vegan::betadiver(mite_beta,"-1")
}

# beta -2 ####
vegan_44 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_44[i] <- vegan::betadiver(mite_beta,"-2") 
}

# beta hk ####
vegan_45 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_45[i] <- vegan::betadiver(mite_beta,"hk") #0.2
}

