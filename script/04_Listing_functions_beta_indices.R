#_______________________________________________________________________________
# Title              : 04_Listing_functions_beta_indices.r
# Date               : 12/02/2025
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

# S2 coeff Gower & Legendre ####
ade4_38_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_38_P[i] <- ade4::dist.binary(mite_beta,method = 10)
}
P38_P <- ls(pattern = "_38_P$")
P38_P <- mget(P38_P)
P38_P <- as.data.frame(P38_P)

# Hamman coefficient ####
ade4_39_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_39_P[i] <- ade4::dist.binary(mite_beta,method = 6)
}
proxy_39_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_39_P[i] <- proxy::dist(mite_beta,method = "Hamman")
}
MultivariateAnalysis_39_P <- rep(NA,69) # A besoin de PA pour etre correct
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  MultivariateAnalysis_39_P[i] <- MultivariateAnalysis::Distancia(mite_beta,Metodo = 19)[1]$Distancia
}
MultBiplotR_39_P <- rep(NA,69) # A besoin de PA pour fonctionner
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  MultBiplotR_39_P[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 9)[1,2]
}
PERMANOVA_39_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_39_P[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 9,transformation = 1)$D[1,2]
}
proxyC_39_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  proxyC_39_P[i] <- 1-proxyC::simil(as.matrix(mite[i,]),as.matrix(mite[i+1,]),method = "hamann") #0.5 ok
}

P39_P <- ls(pattern = "_39_P$")
P39_P <- mget(P39_P)
P39_P <- as.data.frame(P39_P)

# Hamming distance ####
ecodive_40_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_40_P[i] <- ecodive::hamming(mite_beta) 
}
abdiv_40_A <- rep(NA,69) 
for (i in 1:(nrow(mite)-1)){
  abdiv_40_A[i] <- abdiv::hamming(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
Rankcluster_40_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  Rankcluster_40_A[i] <- Rankcluster::distHamming(mite[i,],mite[i+1,])
}
Mercator_40_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Mercator_40_A[i] <- Mercator::binaryDistance(t(mite_beta),metric = "hamming")
}
pegas_40_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pegas_40_A[i] <- pegas::dist.hamming(mite_beta)
}
proxyC_40_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxyC_40_A[i] <- proxyC::dist(as.matrix(mite_beta),method = "hamming")[2,1]
}
bingat_40_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  bingat_40_P[i] <- bingat::calcDistance(mite[i,],mite[i+1,]) 
}
genMCMCDiag_40_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  genMCMCDiag_40_A[i] <- genMCMCDiag::hammingDist(mite[i,],mite[i+1,]) #11 
}

P40_P <- ls(pattern = "_40_P$")
P40_P <- mget(P40_P)
P40_P <- as.data.frame(P40_P)

P40_A <- ls(pattern = "_40_A$")
P40_A <- mget(P40_A)
P40_A <- as.data.frame(P40_A)

# Hamming distance 2 ####
CEGO_41_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  CEGO_41_A[i] <- CEGO::distanceNumericHamming(mite[i,],mite[i+1,])
}
EnvNJ_41_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_41_A[i] <- EnvNJ::metrics(t(mite_beta),method = "hamming")[2,1]
}
ClusterR_41_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_41_A[i] <- ClusterR::distance_matrix(mite_beta,method = "hamming")[2,1]
}
rdist_41_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  rdist_41_A[i] <- rdist::rdist(mite_beta,metric = "hamming")
}
CEGO_41_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  CEGO_41_A[i] <- CEGO::distanceNumericHamming(mite[i,],mite[i+1,])
}
EnsCat_41_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnsCat_41_A[i] <- EnsCat::hammingD(mite_beta)[2,1]
}

P41_A <- ls(pattern = "_41_A$")
P41_A <- mget(P41_A)
P41_A <- as.data.frame(P41_A)

# Hellinger ####
vegan_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_42_A[i] <- vegan::vegdist(mite_beta,method = "hellinger")
}
pctax_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_42_A[i] <- pctax::mat_dist(t((mite_beta)),method = "hellinger")
}
ecodive_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_42_A[i] <- ecodive::hellinger(mite_beta,rescale = T)
}
adespatial_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_42_A[i] <- adespatial::beta.div(mite_beta,method = "hellinger",save.D = T)$D
}
proxy_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_42_A[i] <- proxy::dist(mite_beta,method = "Hellinger")
}
Rfast_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_42_A[i] <- Rfast::Dist(mite_beta,method = "hellinger")[2,1]
}
abdiv_42_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_42_A[i] <- abdiv::hellinger(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P42_A <- ls(pattern = "_42_A$")
P42_A <- mget(P42_A)
P42_A <- as.data.frame(P42_A)

# Jaccard (Abondance) ####
vegan_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_43_A[i] <- vegan::vegdist(mite_beta,method = "jaccard",binary = F) # 0.6936661
}
adespatial_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_43_A[i] <- adespatial::beta.div(mite_beta,method = "ruzicka",save.D = T)$D
}
abdiv_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_43_A[i] <- abdiv::ruzicka(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}
labdsv_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  labdsv_43_A[i] <- labdsv::dsvdis(mite_beta, index = "ruzicka")
}
philentropy_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_43_A[i] <- abdiv::ruzicka(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}
stylo_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  stylo_43_A[i] <- stylo::dist.minmax(mite_beta)
}
adiv_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_43_A[i] <- adiv::distMS(mite_beta)
}
ecodive_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_43_A[i] <- ecodive::soergel(mite_beta,rescale = F) 
}
proxy_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_43_A[i] <- proxy::dist(mite_beta,method = "Soergel")
}
PERMANOVA_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_43_A[i] <- PERMANOVA::DistContinuous(mite_beta,coef = 9)$D[2,1]
}
EnvNJ_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_43_A[i] <- EnvNJ::metrics(t(mite_beta),method = "soergel")
}
Rfast_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_43_A[i] <- Rfast::Dist(mite_beta,method = "soergel")[2,1]
}
NST_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_43_A[i] <- NST::beta.g(mite_beta,dist.method = "jaccard")
}
picante_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  picante_43_A[i] <- 1-picante::species.dist(t(mite_beta),metric = "jaccard")
}
BoutrosLab.plotting.general_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BoutrosLab.plotting.general_43_A[i] <- BoutrosLab.plotting.general::dist(mite_beta,method = "jaccard") 
}
statisfactory_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  statisfactory_43_A[i] <- 1-statisfactory::fuzzyJaccard(mite[i,],mite[i+1,]) 
}
geocmeans_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  geocmeans_43_A[i] <- 1-geocmeans::calc_jaccard_idx(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
adespatial_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_43_A[i] <- adespatial::beta.div.comp(mite_beta,coef = "J",quant = T)$D 
}
pctax_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_43_A[i] <- pctax::mat_dist(t((mite_beta)),method = "jaccard")
}
BAT_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_43_A[i] <- BAT::beta(mite_beta,func = "jaccard",abund = T)$Btotal
}
prabclus_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  prabclus_43_A[i] <- prabclus::jaccard(t(mite_beta)) 
}
Mercator_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Mercator_43_A[i] <- Mercator::binaryDistance(t(mite_beta),metric = "jaccard") #2.266299
}
rdist_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  rdist_43_A[i] <- rdist::rdist(mite_beta,metric = "jaccard") 
}
ConNEct_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  ConNEct_43_A[i] <- ConNEcT::funCorrJacc(as.numeric(mite[i,]),as.numeric(mite[i+1,]))$value #0.9834216
}

P43_A <- ls(pattern = "_43_A$")
P43_A <- mget(P43_A)
P43_A <- as.data.frame(P43_A)

# Jaccard P/A ####
labdsv_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  labdsv_43_P[i] <- labdsv::dsvdis(mite_beta, index = "steinhaus")
}
ecodist_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodist_43_P[i] <- ecodist::distance(mite_beta,method = "jaccard")
}
stats_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  stats_43_P[i] <- dist(mite_beta,method = "binary")
}
vegan_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_43_P[i] <- vegan::vegdist(mite_beta,method = "jaccard",binary = T)
}
neighbr_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  neighbr_43_P[i] <- 1-neighbr::similarity(mite_pa[i,],mite_pa[i+1,],measure = "jaccard") 
}
ChemoSpecUtils_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ChemoSpecUtils_43_P[i] <- ChemoSpecUtils::rowDist(as.matrix(mite_beta),method = "binary") 
}
mgc_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  mgc_43_P[i] <- mgc::mgc.distance(mite_beta,method = "binary")[1,2] 
}
BoutrosLab.plotting.general_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BoutrosLab.plotting.general_43_P[i] <- BoutrosLab.plotting.general::dist(mite_beta,method = "binary")
}
bioregion_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  bioregion_43_P[i] <- bioregion::dissimilarity(as.matrix(mite_beta),metric = "Jaccard")$Jaccard
}
fAssets_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fAssets_43_P[i] <- fAssets::jaccardDist(t(mite_beta))
}
vegan_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_43_P[i] <- 1-vegan::betadiver(mite_beta,"j")
}
ecodive_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_43_P[i] <- ecodive::jaccard(mite_beta)
}
betapart_43_P <- rep(NA,69) # Besoin de PA pour fonctionner
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  betapart_43_P[i] <- betapart::beta.pair(mite_beta,index.family = "jaccard")$beta.jac 
}
adespatial_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_43_P[i] <- adespatial::beta.div.comp(mite_beta,coef = "J",quant = F)$D 
}
tabula_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  tabula_43_P[i] <- 1-tabula::similarity(mite_beta,method = "jaccard")
}
adiv_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_43_P[i] <- adiv::betastatjac(mite_beta)[1]
}
diverse_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  diverse_43_P[i] <-diverse::dis_entities(t(mite_beta),method = "Jaccard",category_row = T)[1,2]
}
BAT_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_43_P[i] <- BAT::beta(mite_beta,func = "jaccard",abund = F)$Btotal
}
proxy_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_43_P[i] <- proxy::dist(mite_beta,method = "Jaccard")
}
PERMANOVA_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_43_P[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 3,transformation = 1)$D[1,2]
}
ClusterR_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_43_P[i] <- ClusterR::distance_matrix(mite_beta,method = "jaccard_coefficient")[2,1]
}
flexclust_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  flexclust_43_P[i] <- flexclust::dist2(mite[i,],mite[i+1,],method = "binary")
}
iTOP_43_P <- rep(NA,69) # Besoin de PA pour fonctionner
for (i in 1:(nrow(mite)-1)){
  iTOP_43_P[i] <- 1-iTOP::jaccard(mite_pa[i,],mite_pa[i+1,])
}
proxyC_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  proxyC_43_P[i] <- 1-proxyC::simil(as.matrix(mite[i,]),as.matrix(mite[i+1,]),method = "jaccard")
}
abdiv_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_43_P[i] <- abdiv::jaccard(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
wiqid_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_43_P[i] <- wiqid::distJaccard(mite[i,],mite[i+1,])
}
fossil_43_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  fossil_43_P[i] <- 1-fossil::jaccard(mite[i,],mite[i+1,])
}

P43_P <- ls(pattern = "_43_P$")
P43_P <- mget(P43_P)
P43_P <- as.data.frame(P43_P)

# Squared-root Jaccard ####
adespatial_44_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_44_P[i] <- print(adespatial::dist.ldc(mite_beta,method = "jaccard"))
}
MultBiplotR_44_P <- rep(NA,69) # A besoin de PA pour fonctionner
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  MultBiplotR_44_P[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 3)[1,2] 
}
ade4_44_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_44_P[i] <- ade4::dist.binary(mite_beta,method = 1) 
}

P44_P <- ls(pattern = "_44_P$")
P44_P <- mget(P44_P)
P44_P <- as.data.frame(P44_P)

# Turnover component of Jaccard dissimilarity defined by Baselga ####
vegan_45_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_45_P[i] <- vegan::nestedbetajac(mite_beta)[1]
}
betapart_45_P <- rep(NA,69) # Besoin de PA pour marcher
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  betapart_45_P[i] <- betapart::beta.pair(mite_beta,index.family = "jaccard")$beta.jtu
}
adespatial_45_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_45_P[i] <- adespatial::beta.div.comp(mite_beta,coef = "BJ")$repl 
}
bioregion_45_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  bioregion_45_P[i] <- bioregion::dissimilarity(as.matrix(mite_beta),metric = "Jaccardturn")$Jaccardturn
}
abdiv_45_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
    abdiv_45_P[i] <- abdiv::jaccard_turnover(as.numeric(mite[i,]),as.numeric(mite[i+1,])) #0.2666667
}

P45_P <- ls(pattern = "_45_P$")
P45_P <- mget(P45_P)
P45_P <- as.data.frame(P45_P)

# Nestedness component of Jaccard dissimilarity defined by Baselga ####
abdiv_46_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_46_P[i] <- abdiv::jaccard_nestedness(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
vegan_46_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_46_P[i] <- vegan::nestedbetajac(mite_beta)[2] 
}
betapart_46_P <- rep(NA,69) # Besoin de PA pour fonctionner
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  betapart_46_P[i] <- betapart::beta.pair(mite_beta,index.family = "jaccard")$beta.jne
}
adespatial_46_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_46_P[i] <- adespatial::beta.div.comp(mite_beta,coef = "BJ")$rich
}

P46_P <- ls(pattern = "_46_P$")
P46_P <- mget(P46_P)
P46_P <- as.data.frame(P46_P)

# Extended Jaccard Similarity ####
adiv_47_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_47_A[i] <- adiv::dsimcom(mite_beta,method = "2",type = "similarity",option = "absolute")[1,2] 
}
diverse_47_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  diverse_47_A[i] <- 1-diverse::dis_entities(t(mite_beta),method = "eJaccard",category_row = T)[1,2] 
}
proxy_47_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_47_A[i] <- 1-proxy::dist(mite_beta,method = "eJaccard")
}
EnvNJ_47_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_47_A[i] <- 1-EnvNJ::metrics(t(mite_beta),method = "jaccard")
}
philentropy_47_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_47_A[i] <- 1-philentropy::jaccard(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
proxyC_47_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  proxyC_47_A[i] <- proxyC::simil(as.matrix(mite[i,]),as.matrix(mite[i+1,]),method = "ejaccard") 
}
adespatial_47_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_47_A[i] <- 1-adespatial::beta.div(mite_beta,method = "wishart",save.D = T)$D 
}
wiqid_47_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_47_A[i] <- 1-wiqid::distSimRatio(mite[i,],mite[i+1,])
}
vegan_47_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_47_A[i] <- 1-vegan::designdist(mite_beta,method = "(A+B-2*J)/(A+B-J)",terms = "quadratic")
}

P47_A <- ls(pattern = "_47_A$")
P47_A <- mget(P47_A)
P47_A <- as.data.frame(P47_A)

# Replacement index of Jaccard defined by Podani ####
adespatial_48_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_48_P[i] <- adespatial::beta.div.comp(mite_beta,coef = "J")$repl
}
BAT_48_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_48_P[i] <- BAT::beta(mite_beta,func = "jaccard",abund = F)$Brepl
}

P48_P <- ls(pattern = "_48_P$")
P48_P <- mget(P48_P)
P48_P <- as.data.frame(P48_P)

# Nestedness component of Jaccard dissimilarity defined by Podani ####
BAT_49_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_49_P[i] <- BAT::beta(mite_beta,func = "jaccard",abund = F)$Brich
}
adespatial_49_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_49_P[i] <- adespatial::beta.div.comp(mite_beta,coef = "J")$rich
}

P49_P <- ls(pattern = "_49_P$")
P49_P <- mget(P49_P)
P49_P <- as.data.frame(P49_P)

# Nestedness component of Jaccard dissimilarity defined by Podani & Schmera ####
adespatial_50_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_50_P[i] <- adespatial::beta.div.comp(mite_beta,coef = "N")$rich 
}

P50_P <- ls(pattern = "_50_P$")
P50_P <- mget(P50_P)
P50_P <- as.data.frame(P50_P)

# Jaccard richness gain ####
BAT_51_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_51_P[i] <- BAT::beta(mite_beta,func = "jaccard",abund = F)$Bgain
}
adiv_51_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_51_P[i] <- adiv::betastatjac(mite_beta)[3]
}

P51_P <- ls(pattern = "_51_P$")
P51_P <- mget(P51_P)
P51_P <- as.data.frame(P51_P)

# Jaccard richness loss ####
BAT_52_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_52_P[i] <- BAT::beta(mite_beta,func = "jaccard",abund = F)$Bloss
}
adiv_52_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_52_P[i] <- adiv::betastatjac(mite_beta)[2]
}

P52_P <- ls(pattern = "_52_P$")
P52_P <- mget(P52_P)
P52_P <- as.data.frame(P52_P)

# Jeffreys ####
EnvNJ_53_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_53_A[i] <- EnvNJ::metrics(t(mite_beta),method = "jeffreys")
}

P53_A <- ls(pattern = "_53_A$")
P53_A <- mget(P53_A)
P53_A <- as.data.frame(P53_A)


# Jensen-Shannon distance ####
adiv_54_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_54_A[i] <- ecodive::jensen(mite_beta,rescale = F)
}

P54_A <- ls(pattern = "_54_A$")
P54_A <- mget(P54_A)
P54_A <- as.data.frame(P54_A)

# Jensen-Shannon divergence ####
adiv_55_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_55_A[i] <- ecodive::jsd(mite_beta,rescale = F)
}
EnvNJ_55_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_55_A[i] <- EnvNJ::metrics(t(mite_beta),method = "jensen-shannon")
}
Rfast_55_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_55_A[i] <- Rfast::Dist(mite_beta,method = "jensen_shannon")[1,2] #34.28249
}
Compositional_55_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Compositional_55_A[i] <- Compositional::divergence(mite_beta,type = "jensen_shannon")[1,2]
}
priorsense_55_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  priorsense_55_A[i] <- priorsense::cjs_dist(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
philentropy_55_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_55_A[i] <- philentropy::jensen_difference(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}

P55_A <- ls(pattern = "_55_A$")
P55_A <- mget(P55_A)
P55_A <- as.data.frame(P55_A)

# Kulczynski ####
proxy_56_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_56_P[i] <- proxy::dist(mite_beta,method = "Kulczynski2") 
}
abdiv_56_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_56_P[i] <- abdiv::kulczynski_second(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}
vegan_56_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_56_P[i] <- vegan::betadiver(mite_beta,"co")
}
prabclus_56_P <- rep(NA,69)  # a besoin de PA pour etre correct
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  prabclus_56_P[i] <- prabclus::kulczynski(t(mite_beta))[2,1] 
}
fossil_56_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  fossil_56_P[i] <- 1-fossil::kulczynski(mite[i,],mite[i+1,]) 
}

P56_P <- ls(pattern = "_56_P$")
P56_P <- mget(P56_P)
P56_P <- as.data.frame(P56_P)

# Kulczynski 2
vegan_57_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_57_A[i] <- vegan::vegdist(mite_beta,method = "kulczynski") 
}
pctax_57_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_57_A[i] <- pctax::mat_dist(t((mite_beta)),method = "kulczynski") 
}
adespatial_57_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_57_A[i] <- adespatial::beta.div(mite_beta,method = "kulczynski",save.D = T)$D 
}
NST_57_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_57_A[i] <- NST::beta.g(mite_beta,dist.method = "kulczynski")
}
prabclus_57_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  prabclus_57_A[i] <- prabclus::qkulczynski(t(mite_beta))[1,2]
}
proxy_57_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_57_A[i] <- proxy::dist(mite_beta,method = "Kulczynski1") 
}
abdiv_57_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_57_A[i] <- abdiv::kulczynski_first(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P57_A <- ls(pattern = "_57_A$")
P57_A <- mget(P57_A)
P57_A <- as.data.frame(P57_A)

# Kulczynski 3 ####
philentropy_58_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_58_A[i] <- philentropy::kulczynski_d(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
EnvNJ_58_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_58_A[i] <- EnvNJ::metrics(t(mite_beta),method = "kulczynski")
}
Rfast_58_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_58_A[i] <- Rfast::Dist(mite_beta,method = "kulczynski")[1,2]
}

P58_A <- ls(pattern = "_58_A$")
P58_A <- mget(P58_A)
P58_A <- as.data.frame(P58_A)

# Lorentzian distance ####
ecodive_59_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_59_A[i] <- ecodive::lorentzian(mite_beta,rescale = F) 
}
EnvNJ_59_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_59_A[i] <- EnvNJ::metrics(t(mite_beta),method = "lorentzian") 
}
philentropy_59_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  philentropy_59_A[i] <- philentropy::lorentzian(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}

P59_A <- ls(pattern = "_59_A$")
P59_A <- mget(P59_A)
P59_A <- as.data.frame(P59_A)

# Mahalanobis ####
vegan_60_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_60_P[i] <- vegan::vegdist(mite_beta,method = "mahalanobis") 
}
FD_60_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  FD_60_P[i] <- FD::mahaldis(as.matrix(mite_beta)) 
}
pctax_60_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_60_P[i] <- pctax::mat_dist(t((mite_beta)),method = "mahalanobis")
}
ClusterR_60_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_60_P[i] <- ClusterR::distance_matrix(mite_beta,method = "mahalanobis")[2,1] 
}
NST_60_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_60_P[i] <- NST::beta.g(mite_beta,dist.method = "mahalanobis")
}

P60_P <- ls(pattern = "_60_P$")
P60_P <- mget(P60_P)
P60_P <- as.data.frame(P60_P)

# Manhattan ####
ecodist_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodist_61_A[i] <- ecodist::distance(mite_beta,method = "manhattan")
}
stats_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  stats_61_A[i] <- dist(mite_beta,method = "manhattan")
}
vegan_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_61_A[i] <- vegan::vegdist(mite_beta,method = "manhattan")
}
proxyC_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxyC_61_A[i] <- proxyC::dist(as.matrix(mite_beta),method = "manhattan")[1,2]
}
mgc_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  mgc_61_A[i] <- mgc::mgc.distance(mite_beta,method = "manhattan")[1,2]
}
cluster_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  cluster_61_A[i] <- cluster::daisy(mite_beta,metric = "manhattan")
}
ChemoSpecUtils_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ChemoSpecUtils_61_A[i] <- ChemoSpecUtils::rowDist(as.matrix(mite_beta),method = "manhattan")
}
ldt_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ldt_61_A[i] <- ldt::s.distance(t(mite_beta),distance = "manhattan")
}
fAssets_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fAssets_61_A[i] <- fAssets::manhattanDist(t(mite_beta))
}
ecodive_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_61_A[i] <- ecodive::manhattan(mite_beta,rescale = F)
}
adespatial_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_61_A[i] <- print(adespatial::dist.ldc(mite_beta,method = "manhattan")) 
}
diverse_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  diverse_61_A[i] <- diverse::dis_entities(t(mite_beta),method = "Manhattan",category_row = T)[1,2]
}
pctax_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_61_A[i] <- pctax::mat_dist(t((mite_beta)),method = "manhattan")
}
proxy_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_61_A[i] <- proxy::dist(mite_beta,method = "Manhattan") 
}
amap_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  amap_61_A[i] <- amap::Dist(mite_beta,method = "manhattan")
}
Mercator_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Mercator_61_A[i] <- Mercator::binaryDistance(t(mite_beta),metric = "manhattan")
}
dynutils_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  dynutils_61_A[i] <- dynutils::calculate_distance(mite_beta,method = "manhattan")[1,2]
}
EnvNJ_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_61_A[i] <- EnvNJ::metrics(t(mite_beta),method = "manhattan")[1,2]
}
Rfast_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_61_A[i] <- Rfast::Dist(mite_beta,method = "manhattan")[1,2] #95.06
}
comparator_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  comparator_61_A[i] <- comparator::Manhattan()(as.numeric(mite[i,]),as.numeric(mite[i+1,])) #95.06
}
ClusterR_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_61_A[i] <- ClusterR::distance_matrix(mite_beta,method = "manhattan")[2,1]
}
rdist_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  rdist_61_A[i] <- rdist::rdist(mite_beta,metric = "manhattan")
}
BoutrosLab.plotting.general_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BoutrosLab.plotting.general_61_A[i] <- BoutrosLab.plotting.general::dist(mite_beta,method = "manhattan")
}
coda.base_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  coda.base_61_A[i] <- coda.base::dist(mite_beta,"manhattan")
}
NST_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_61_A[i] <- NST::beta.g(mite_beta,dist.method = "manhattan")
}
Rlof_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rlof_61_A[i] <- Rlof::distmc(mite_beta,method = "manhattan")
}
fda.usc_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fda.usc_61_A[i] <- fda.usc::metric.dist(mite_beta,method = "manhattan")[1,2]
}
analogue_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  analogue_61_A[i] <- analogue::distance(mite[i,],mite[i+1,],method = "manhattan") 
}
hmsr_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  hmsr_61_A[i] <- hmsr::manhattan_distance(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
abdiv_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_61_A[i] <- abdiv::manhattan(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
fossil_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  fossil_61_A[i] <- fossil::manhattan(mite[i,],mite[i+1,])
}
philentropy_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_61_A[i] <- philentropy::manhattan(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
TSdist_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  TSdist_61_A[i] <- TSdist::LPDistance(as.numeric(mite[i,]),as.numeric(mite[i+1,]),method = "manhattan") #95.06
}
LearnClust_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  LearnClust_61_A[i] <- LearnClust::mdistance(as.numeric(mite[i,]),as.numeric(mite[i+1,])) #11.08
}
flexclust_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  flexclust_61_A[i] <- flexclust::dist2(mite[i,],mite[i+1,],method = "manhattan")
}

P61_A <- ls(pattern = "_61_A$")
P61_A <- mget(P61_A)
P61_A <- as.data.frame(P61_A)

# Manhattan modified ####
NST_62_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_62_A[i] <- NST::beta.g(mite_beta,dist.method = "mManhattan")
}

P62_A <- ls(pattern = "_62_A$")
P62_A <- mget(P62_A)
P62_A <- as.data.frame(P62_A)

# Mean character difference ####
abdiv_63_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  abdiv_63_P[i] <- abdiv::mean_character_difference(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P63_P <- ls(pattern = "_63_P$")
P63_P <- mget(P63_P)
P63_P <- as.data.frame(P63_P)

# Modified mean character difference ####
adespatial_64_A <- rep(NA,69) # Donne la meme valeur que P avec des abondances mais si on donne du PA alors la valeur est differente
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_64_A[i] <- print(adespatial::dist.ldc(mite_beta,method = "modmeanchardiff")) 
}
abdiv_64_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_64_P[i] <- abdiv::modified_mean_character_difference(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P64_A <- ls(pattern = "_64_A$")
P64_A <- mget(P64_A)
P64_A <- as.data.frame(P64_A)

P64_P <- ls(pattern = "_64_P$")
P64_P <- mget(P64_P)
P64_P <- as.data.frame(P64_P)

# Matusita ####
ecodive_65_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_65_A[i] <- ecodive::matusita(mite_beta,rescale = F) #6.514023
}

P65_A <- ls(pattern = "_65_A$")
P65_A <- mget(P65_A)
P65_A <- as.data.frame(P65_A)

# Minkowski ####-------------------------------------------------------------
stats_66_1_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  stats_66_1_A[i] <- dist(mite_beta,method = "minkowski",p=1)
}
stats_66_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  stats_66_2_A[i] <- dist(mite_beta,method = "minkowski",p=2)
}
stats_66_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  stats_66_3_A[i] <- dist(mite_beta,method = "minkowski",p=3)
}
flexclust_66_1_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  flexclust_66_1_A[i] <- flexclust::dist2(mite[i,],mite[i+1,],method = "minkowski",p=1) #40.37368
}
flexclust_66_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  flexclust_66_2_A[i] <- flexclust::dist2(mite[i,],mite[i+1,],method = "minkowski",p=2)
}
flexclust_66_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  flexclust_66_3_A[i] <- flexclust::dist2(mite[i,],mite[i+1,],method = "minkowski",p=3)
}
mgc_66_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  mgc_66_2_A[i] <- mgc::mgc.distance(mite_beta,method = "minkowski")[1,2]
}
ecodive_66_1_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_66_1_A[i] <- ecodive::minkowski(mite_beta,rescale = F,power = 1)
}
ecodive_66_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_66_2_A[i] <- ecodive::minkowski(mite_beta,rescale = F,power = 2)
}
ecodive_66_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_66_3_A[i] <- ecodive::minkowski(mite_beta,rescale = F,power = 3)
}
abdiv_66_1_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_66_1_A[i] <- abdiv::minkowski(as.numeric(mite[i,]),as.numeric(mite[i+1,]),p=1)
}
abdiv_66_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_66_2_A[i] <- abdiv::minkowski(as.numeric(mite[i,]),as.numeric(mite[i+1,]),p=2)
}
abdiv_66_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_66_3_A[i] <- abdiv::minkowski(as.numeric(mite[i,]),as.numeric(mite[i+1,]),p=3)
}
proxy_66_1_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_66_1_A[i] <- proxy::dist(mite_beta,method = "Minkowski",p=1)
}
proxy_66_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_66_2_A[i] <- proxy::dist(mite_beta,method = "Minkowski",p=2)
}
proxy_66_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_66_3_A[i] <- proxy::dist(mite_beta,method = "Minkowski",p=3)
}
PERMANOVA_66_1_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_66_1_A[i] <- PERMANOVA::DistContinuous(mite_beta,coef = 4,r=1)$D[1,2]
}
PERMANOVA_66_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_66_2_A[i] <- PERMANOVA::DistContinuous(mite_beta,coef = 4,r=2)$D[1,2]
}
PERMANOVA_66_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_66_3_A[i] <- PERMANOVA::DistContinuous(mite_beta,coef = 4,r=3)$D[1,2]
}
dynutils_66_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  dynutils_66_2_A[i] <- dynutils::calculate_distance(mite_beta,method = "minkowski")[1,2]
}
EnvNJ_66_1_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_66_1_A[i] <- EnvNJ::metrics(t(mite_beta),method = "minkowski",p=1)[1,2]
}
EnvNJ_66_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_66_2_A[i] <- EnvNJ::metrics(t(mite_beta),method = "minkowski",p=2)[1,2]
}
EnvNJ_66_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_66_3_A[i] <- EnvNJ::metrics(t(mite_beta),method = "minkowski",p=3)[1,2]
}
Rfast_66_1_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_66_1_A[i] <- Rfast::Dist(mite_beta,method = "minkowski",p = 1)[1,2] #95.06
}
Rfast_66_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_66_2_A[i] <- Rfast::Dist(mite_beta,method = "minkowski",p = 2)[1,2] #95.06
}
Rfast_66_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_66_3_A[i] <- Rfast::Dist(mite_beta,method = "minkowski",p = 3)[1,2] #95.06
}
fAssets_66_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fAssets_66_2_A[i] <- fAssets::minkowskiDist(t(mite_beta)) #40.37368
}
proxyC_66_1_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxyC_66_1_A[i] <- proxyC::dist(as.matrix(mite_beta),method = "minkowski",p=1)[1,2]
}
proxyC_66_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxyC_66_2_A[i] <- proxyC::dist(as.matrix(mite_beta),method = "minkowski",p=2)[1,2]
}
proxyC_66_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxyC_66_3_A[i] <- proxyC::dist(as.matrix(mite_beta),method = "minkowski",p=3)[1,2]
}
ClusterR_66_1_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_66_1_A[i] <- ClusterR::distance_matrix(mite_beta,method = "minkowski",minkowski_p = 1)[2,1]
}
ClusterR_66_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_66_2_A[i] <- ClusterR::distance_matrix(mite_beta,method = "minkowski",minkowski_p = 2)[2,1]
}
ClusterR_66_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_66_3_A[i] <- ClusterR::distance_matrix(mite_beta,method = "minkowski",minkowski_p = 3)[2,1]
}
rdist_66_1_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  rdist_66_1_A[i] <- rdist::rdist(mite_beta,metric = "minkowski",p = 1)
}
rdist_66_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  rdist_66_2_A[i] <- rdist::rdist(mite_beta,metric = "minkowski",p = 2)
}
rdist_66_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  rdist_66_3_A[i] <- rdist::rdist(mite_beta,metric = "minkowski",p = 3)
}
coda.base_66_1_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  coda.base_66_1_A[i] <- coda.base::dist(mite_beta,"minkowski",p=1)
}
coda.base_66_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  coda.base_66_2_A[i] <- coda.base::dist(mite_beta,"minkowski",p=2) 
}
coda.base_66_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  coda.base_66_3_A[i] <- coda.base::dist(mite_beta,"minkowski",p=3)
}
fda.usc_66_1_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fda.usc_66_1_A[i] <- fda.usc::metric.dist(mite_beta,method = "minkowski",p=1)[1,2]
}
fda.usc_66_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fda.usc_66_2_A[i] <- fda.usc::metric.dist(mite_beta,method = "minkowski",p=2)[1,2]
}
fda.usc_66_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fda.usc_66_3_A[i] <- fda.usc::metric.dist(mite_beta,method = "minkowski",p=3)[1,2]
}
BoutrosLab.plotting.general_66_1_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BoutrosLab.plotting.general_66_1_A[i] <- BoutrosLab.plotting.general::dist(mite_beta,method = "minkowski",p=1)
}
BoutrosLab.plotting.general_66_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BoutrosLab.plotting.general_66_2_A[i] <- BoutrosLab.plotting.general::dist(mite_beta,method = "minkowski",p=2)
}
BoutrosLab.plotting.general_66_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BoutrosLab.plotting.general_66_3_A[i] <- BoutrosLab.plotting.general::dist(mite_beta,method = "minkowski",p=3)
}
comparator_66_1_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  comparator_66_1_A[i] <- comparator::Minkowski(p=1)(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
comparator_66_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  comparator_66_2_A[i] <- comparator::Minkowski(p=2)(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
comparator_66_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  comparator_66_3_A[i] <- comparator::Minkowski(p=3)(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
TSdist_66_1_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  TSdist_66_1_A[i] <- TSdist::MinkowskiDistance(as.numeric(mite[i,]),as.numeric(mite[i+1,]),p=2)
}
TSdist_66_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  TSdist_66_2_A[i] <- TSdist::MinkowskiDistance(as.numeric(mite[i,]),as.numeric(mite[i+1,]),p=2)
}
TSdist_66_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  TSdist_66_3_A[i] <- TSdist::MinkowskiDistance(as.numeric(mite[i,]),as.numeric(mite[i+1,]),p=2)
}
philentropy_66_1_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_66_1_A[i] <- philentropy::minkowski(as.numeric(mite[i,]),as.numeric(mite[i+1,]),n=1) 
}
philentropy_66_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_66_2_A[i] <- philentropy::minkowski(as.numeric(mite[i,]),as.numeric(mite[i+1,]),n=2) 
}
philentropy_66_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_66_3_A[i] <- philentropy::minkowski(as.numeric(mite[i,]),as.numeric(mite[i+1,]),n=3)
}

P66_1_A <- ls(pattern = "_66_1_A$")
P66_1_A <- mget(P66_1_A)
P66_1_A <- as.data.frame(P66_1_A)

P66_2_A <- ls(pattern = "_66_2_A$")
P66_2_A <- mget(P66_2_A)
P66_2_A <- as.data.frame(P66_2_A)

P66_3_A <- ls(pattern = "_66_3_A$")
P66_3_A <- mget(P66_3_A)
P66_3_A <- as.data.frame(P66_3_A)

# Morisita ####
abdiv_67_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_67_A[i] <- abdiv::morisita(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
vegan_67_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_67_A[i] <- vegan::vegdist(mite_beta,method = "morisita")
}
pctax_67_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_67_A[i] <- pctax::mat_dist(t((mite_beta)),method = "morisita")
}
ecodive_67_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_67_A[i] <- ecodive::morisita(mite_beta)
}
NST_67_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_67_A[i] <- NST::beta.g(mite_beta,dist.method = "morisita")
}

P67_A <- ls(pattern = "_67_A$")
P67_A <- mget(P67_A)
P67_A <- as.data.frame(P67_A)

# Morisita-Horn ####
vegan_68_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_68_A[i] <- vegan::vegdist(mite_beta,method = "horn")
}
NST_68_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_68_A[i] <- NST::beta.g(mite_beta,dist.method = "horn")
}
pctax_68_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_68_A[i] <- pctax::mat_dist(t((mite_beta)),method = "horn")
}
ecodive_68_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_68_A[i] <- ecodive::horn(mite_beta,rescale = F)
}
tabula_68_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  tabula_68_A[i] <- 1-tabula::similarity(mite_beta,method = "morisita")
}
abdiv_68_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_68_A[i] <- abdiv::horn_morisita(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
wiqid_68_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_68_A[i] <- wiqid::distMorisitaHorn(mite[i,],mite[i+1,])
}
fossil_68_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  fossil_68_A[i] <- 1-fossil::morisita.horn(mite[i,],mite[i+1,])
}

P68_A <- ls(pattern = "_68_A$")
P68_A <- mget(P68_A)
P68_A <- as.data.frame(P68_A)

# Motyka ####
pctax_69_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_69_A[i] <- ecodive::motyka(mite_beta,rescale = F)
}
EnvNJ_69_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_69_A[i] <- EnvNJ::metrics(t(mite_beta),method = "motyka")
}
Rfast_69_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_69_A[i] <- Rfast::Dist(mite_beta,method = "motyka")[1,2]
}
philentropy_69_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  philentropy_69_A[i] <- philentropy::motyka(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P69_A <- ls(pattern = "_69_A$")
P69_A <- mget(P69_A)
P69_A <- as.data.frame(P69_A)

# Mountford ####
vegan_70_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_70_P[i] <- vegan::vegdist(mite_beta,method = "mountford")
}
pctax_70_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_70_P[i] <- pctax::mat_dist(t((mite_beta)),method = "mountford")
}
proxy_70_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_70_P[i] <- proxy::dist(mite_beta,method = "Mountford")
}
NST_70_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_70_P[i] <- NST::beta.g(mite_beta,dist.method = "mountford")
}

P70_P <- ls(pattern = "_70_P$")
P70_P <- mget(P70_P)
P70_P <- as.data.frame(P70_P)

# Ochiai ####
adespatial_71_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_71_P[i] <- adespatial::beta.div(mite_beta,method = "ochiai",save.D = T,sqrt.D = F)$D
}
ade4_71_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_71_P[i] <- ade4::dist.binary(mite_beta,method = 7)
}
adiv_71_A <- rep(NA,69) # Donne la meme valeur que spaa avec PA alors que spaa avec A
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_71_A[i] <- adiv::dsimcom(mite_beta,method = "4",type = "similarity",option = "absolute")[1,2]
}
MultBiplotR_71_P <- rep(NA,69) # Besoin de PA pour fonctionner
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  MultBiplotR_71_P[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 12)[1,2]
}
spaa_71_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  spaa_71_A[i] <- spaa::sp.pair(t(as.matrix(mite_beta)))$Ochiai
}

P71_P <- ls(pattern = "_71_P$")
P71_P <- mget(P71_P)
P71_P <- as.data.frame(P71_P)

P71_A <- ls(pattern = "_71_A$")
P71_A <- mget(P71_A)
P71_A <- as.data.frame(P71_A)

# Otsuka-Ochiai ####
ecodive_72_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_72_P[i] <- ecodive::ochiai(mite_beta)
}
PERMANOVA_72_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_72_P[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 12,transformation = 1)$D[1,2]
}
proxy_72_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_72_P[i] <- proxy::dist(mite_beta,method = "Ochiai")
}
vegan_72_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_72_P[i] <- vegan::designdist(mite_beta,method = "1-J/sqrt(A*B)",terms = "binary")
}
labdsv_72_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  labdsv_72_P[i] <- labdsv::dsvdis(mite_beta, index = "ochiai")
}
wiqid_72_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_72_P[i] <- wiqid::distOchiai(mite[i,],mite[i+1,])
}
fossil_72_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  fossil_72_P[i] <- 1-fossil::ochiai(mite[i,],mite[i+1,])
}

P72_P <- ls(pattern = "_72_P$")
P72_P <- mget(P72_P)
P72_P <- as.data.frame(P72_P)

# Phi-squared ####
ade4_73_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  ade4_73_P[i] <- ade4::dist.binary(mite_beta,method = 9)
}
proxy_73_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_73_P[i] <- proxy::dist(mite_beta,method = "Phi-squared")
}
MultBiplotR_73_P <- rep(NA,69) # Besoin de PA pour marcher
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  MultBiplotR_73_P[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 14)[1,2]
}
PERMANOVA_73_P <- rep(NA,69) # Besoin de PA pour etre correct
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  PERMANOVA_73_P[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 14,transformation = 3)$D[1,2]
}
philentropy_73_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  philentropy_73_P[i] <- philentropy::pearson_chi_sq(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P73_P <- ls(pattern = "_73_P$")
P73_P <- mget(P73_P)
P73_P <- as.data.frame(P73_P)

# Preston's coefficient of faunal dissimilarity ####
wiqid_74_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_74_P[i] <- wiqid::distPreston(mite[i,],mite[i+1,])
}

P74_P <- ls(pattern = "_74_P$")
P74_P <- mget(P74_P)
P74_P <- as.data.frame(P74_P)

# Raup ####
vegan_75_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_75_A[i] <- vegan::vegdist(mite_beta,method = "raup")
}
pctax_75_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_75_A[i] <- pctax::mat_dist(t((mite_beta)),method = "raup")
}
iCAMP_75_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  iCAMP_75_A[i] <- iCAMP::RC.pc(mite_beta)$index[1,2]
}
NST_75_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_75_A[i] <- NST::beta.g(mite_beta,dist.method = "raup")
}

P75_A <- ls(pattern = "_75_A$")
P75_A <- mget(P75_A)
P75_A <- as.data.frame(P75_A)

# Rogers & Tanimoto ####
ade4_76_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_76_P[i] <- ade4::dist.binary(mite_beta,method = 4)
}
abdiv_76_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_76_P[i] <- abdiv::rogers_tanimoto(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
proxy_76_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_76_P[i] <- proxy::dist(mite_beta,method = "Tanimoto") 
}
neighbr_76_P <- rep(NA,69) # Besoin de PA pour marcher
for (i in 1:(nrow(mite)-1)){
  neighbr_76_P[i] <- 1-neighbr::similarity(mite_pa[i,],mite_pa[i+1,],measure = "tanimoto")
}
wiqid_76_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_76_P[i] <- wiqid::distRogersTanimoto(mite[i,],mite[i+1,])
}
philentropy_76_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_76_P[i] <- philentropy::tanimoto(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
MultivariateAnalysis_76_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  MultivariateAnalysis_76_P[i] <- MultivariateAnalysis::Distancia(mite_beta,Metodo = 15)[1]$Distancia 
}
MultBiplotR_76_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  MultBiplotR_76_P[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 6)[1,2]
}
PERMANOVA_76_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  PERMANOVA_76_P[i] <- 1-PERMANOVA::DistBinary(mite_beta,coefficient = 6,transformation = 1)$D[1,2]
}
shipunov_76_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  shipunov_76_P[i] <- shipunov::SM.dist(mite_beta)
}

P76_P <- ls(pattern = "_76_P$")
P76_P <- mget(P76_P)
P76_P <- as.data.frame(P76_P)

# Root mean square ####
EnvNJ_77_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_77_A[i] <- EnvNJ::metrics(t(mite_beta),method = "avg")
}
abdiv_77_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_77_A[i] <- abdiv::rms_distance(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P77_A <- ls(pattern = "_77_A$")
P77_A <- mget(P77_A)
P77_A <- as.data.frame(P77_A)

# Russel-Rao ####
abdiv_78_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_78_P[i] <- abdiv::russel_rao(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
proxy_78_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_78_P[i] <- proxy::dist(mite_beta,method = "Russel")
}
MultivariateAnalysis_78_P <- rep(NA,69) # Besoin de PA pour etre correct
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  MultivariateAnalysis_78_P[i] <- MultivariateAnalysis::Distancia(mite_beta,Metodo = 16)[1]$Distancia
}
MultBiplotR_78_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  MultBiplotR_78_P[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 2)[1,2]
}
PERMANOVA_78_P <- rep(NA,69) # Besoin de PA pour etre correct
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  PERMANOVA_78_P[i] <- 1-PERMANOVA::DistBinary(mite_beta,coefficient = 2,transformation = 1)$D[1,2]
}
Mercator_78_P <- rep(NA,69) # Besoin de PA pour etre correct
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  Mercator_78_P[i] <- Mercator::binaryDistance(t(mite_beta),metric = "russellRao")
}

P78_P <- ls(pattern = "_78_P$")
P78_P <- mget(P78_P)
P78_P <- as.data.frame(P78_P)

# Simple match ####
ade4_79_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_79_P[i] <- ade4::dist.binary(mite_beta,method = 2)
}
neighbr_79_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  neighbr_79_P[i] <- neighbr::similarity(mite_pa[i,],mite_pa[i+1,],measure = "simple_matching")
}
proxyC_79_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  proxyC_79_P[i] <- proxyC::simil(as.matrix(mite[i,]),as.matrix(mite[i+1,]),method = "simple matching")
}
shipunov_79_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  shipunov_79_P[i] <- shipunov::SM.dist(mite_beta)
}
nomclust_79_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  nomclust_79_P[i] <- nomclust::sm(mite_beta)
}
proxy_79_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_79_P[i] <- proxy::dist(mite_beta,method = "simple matching")
}
ClusterR_79_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_79_P[i] <- ClusterR::distance_matrix(mite_beta,method = "simple_matching_coefficient")[2,1]
}
PERMANOVA_79_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_79_P[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 4,transformation = 1)$D[1,2]
}
wiqid_79_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_79_P[i] <- 1-wiqid::distMatching(mite[i,],mite[i+1,])
}
MultBiplotR_79_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  MultBiplotR_79_P[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 4)[1,2]
}
arules_79_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  arules_79_P[i] <- 1-arules::dissimilarity(as.matrix(mite_beta),method = "matching")
}

P79_P <- ls(pattern = "_79_P$")
P79_P <- mget(P79_P)
P79_P <- as.data.frame(P79_P)

# Sokal & Sneath ####
abdiv_80_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_80_P[i] <- abdiv::sokal_sneath(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P80_P <- ls(pattern = "_80_P$")
P80_P <- mget(P80_P)
P80_P <- as.data.frame(P80_P)

# Sokal & Sneath 2 ####
adiv_81_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_81_A[i] <- adiv::dsimcom(mite_beta,method = "1",type = "dissimilarity",option="absolute")
}

P81_A <- ls(pattern = "_81_A$")
P81_A <- mget(P81_A)
P81_A <- as.data.frame(P81_A)

# Sokal & Sneath 3 ####
MultivariateAnalysis_82_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  MultivariateAnalysis_82_P[i] <- MultivariateAnalysis::Distancia(mite_beta,Metodo = 14)[1]$Distancia
}
MultBiplotR_82_P <- rep(NA,69) # Besoin de PA
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  MultBiplotR_82_P[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 8)[1,2]
}
PERMANOVA_82_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  PERMANOVA_82_P[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 8,transformation = 3)$D[2,1]
}
Mercator_82_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  Mercator_82_P[i] <- Mercator::binaryDistance(t(mite_beta),metric = "sokalMichener")
}

P82_P <- ls(pattern = "_82_P$")
P82_P <- mget(P82_P)
P82_P <- as.data.frame(P82_P)

# Sokal & Sneath S5 ####
ade4_83_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_83_P[i] <- ade4::dist.binary(mite_beta,method = 3)
}

P83_P <- ls(pattern = "_83_P$")
P83_P <- mget(P83_P)
P83_P <- as.data.frame(P83_P)

# Sokal & Sneath S13 ####
ade4_84_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_84_P[i] <- ade4::dist.binary(mite_beta,method = 8)
}
MultBiplotR_84_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  MultBiplotR_84_P[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 13)[1,2]
}
PERMANOVA_84_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  PERMANOVA_84_P[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 13,transformation = 1)$D[1,2]
}

P84_P <- ls(pattern = "_84_P$")
P84_P <- mget(P84_P)
P84_P <- as.data.frame(P84_P)

# Sorensen ####
vegan_85_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_85_P[i] <- 1-vegan::betadiver(mite_beta,"sor")
}
ecodist_85_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodist_85_P[i] <- ecodist::distance(mite_beta,method = "sorensen")
}
bioregion_85_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  bioregion_85_P[i] <- bioregion::dissimilarity(as.matrix(mite_beta),metric = "Sorensen")$Sorensen
}
ecodive_85_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_85_P[i] <- ecodive::sorensen(mite_beta)
}
betapart_85_P <- rep(NA,69) # Besoin de PA
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  betapart_85_P[i] <- betapart::beta.pair(mite_beta,index.family = "sorensen")$beta.sor
}
tabula_85_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  tabula_85_P[i] <- 1-tabula::similarity(mite_beta,method = "sorensen")
}
adiv_85_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_85_P[i] <- adiv::betastatsor(mite_beta)[1] #0.2
}
adespatial_85_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_85_P[i] <- adespatial::beta.div.comp(mite_beta,coef = "S",quant = F)$D #0.2
}
BAT_85_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_85_P[i] <- BAT::beta(mite_beta,func = "sorensen",abund = F)$Btotal #0.2
}
proxy_85_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_85_P[i] <- proxy::dist(mite_beta,method = "Dice") #0.2
}
labdsv_85_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  labdsv_85_P[i] <- labdsv::dsvdis(mite_beta, index = "sorensen")
}
PERMANOVA_85_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_85_P[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 7,transformation = 1)$D[1,2]
}
fAssets_85_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fAssets_85_P[i] <- fAssets::sorensenDist(t(mite_beta))
}
abdiv_85_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_85_P[i] <- abdiv::sorenson(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
proxyC_85_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  proxyC_85_P[i] <- 1-proxyC::simil(as.matrix(mite[i,]),as.matrix(mite[i+1,]),method = "dice")
}
wiqid_85_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_85_P[i] <- wiqid::distSorensen(mite[i,],mite[i+1,])
}
fossil_85_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  fossil_85_P[i] <- 1-fossil::sorenson(mite[i,],mite[i+1,])
}
prabclus_85_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  prabclus_85_P[i] <- prabclus::dicedist(t(mite_beta))[1,2]
}
spaa_85_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  spaa_85_P[i] <- 1-spaa::sp.pair(t(as.matrix(mite_beta)))$Dice
}

P85_P <- ls(pattern = "_85_P$")
P85_P <- mget(P85_P)
P85_P <- as.data.frame(P85_P)

# Square-root Sorensen ####
abdiv_86_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_86_P[i] <- ade4::dist.binary(mite_beta,method = 5)
}
adespatial_86_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_86_P[i] <- print(adespatial::dist.ldc(mite_beta,method = "sorensen")) #0.4472136
}
MultBiplotR_86_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  MultBiplotR_86_P[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 7)[1,2]
}
P86_P <- ls(pattern = "_86_P$")
P86_P <- mget(P86_P)
P86_P <- as.data.frame(P86_P)

# Turnover component of Sorensen dissimilarity ####
adespatial_87_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_87_P[i] <- adespatial::beta.div.comp(mite_beta,coef = "BS")$repl # Baselga 0.1538462
}
BAT_87_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_87_P[i] <- BAT::beta(mite_beta,func = "sorensen",abund = F)$Brepl
}
betapart_87_P <- rep(NA,69) # Besoin de PA pour fonctionner
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  betapart_87_P[i] <- betapart::beta.pair(mite_beta,index.family = "sorensen")$beta.sim 
}
proxy_87_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_87_P[i] <- proxy::dist(mite_beta,method = "Simpson")
}
vegan_87_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_87_P[i] <- vegan::nestedbetasor(mite_beta)[1]
}
fossil_87_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fossil_87_P[i] <- 1-fossil::simpson(mite[i,],mite[i+1,])
}

P87_P <- ls(pattern = "_87_P$")
P87_P <- mget(P87_P)
P87_P <- as.data.frame(P87_P)


# Nestedness-resultant component of Sørensen dissimilarity ####
vegan_88_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_88_P[i] <- vegan::nestedbetasor(mite_beta)[2]
}
betapart_88_P <- rep(NA,69) # Besoin de PA pour fonctionner
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  betapart_88_P[i] <- betapart::beta.pair(mite_beta,index.family = "sorensen")$beta.sne
}
abdiv_88_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_88_P [i] <- abdiv::sorenson_nestedness(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
adespatial_88_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_88_P [i] <- adespatial::beta.div.comp(mite_beta,coef = "BS")$rich
}

P88_P <- ls(pattern = "_88_P$")
P88_P <- mget(P88_P)
P88_P <- as.data.frame(P88_P)

# Sorensen richness gain ####
BAT_89_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_89_P[i] <- BAT::beta(mite_beta,func = "sorensen",abund = F)$Bgain
}
adiv_89_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_89_P[i] <- adiv::betastatsor(mite_beta)[3]
}

P89_P <- ls(pattern = "_89_P$")
P89_P <- mget(P89_P)
P89_P <- as.data.frame(P89_P)

# Sorensen richness loss ####
BAT_90_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_90_P[i] <- BAT::beta(mite_beta,func = "sorensen",abund = F)$Bloss
}
adiv_90_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_90_P[i] <- adiv::betastatsor(mite_beta)[2]
}

P90_P <- ls(pattern = "_90_P$")
P90_P <- mget(P90_P)
P90_P <- as.data.frame(P90_P)

# Legendre replacement index ####
adespatial_91_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_91_P[i] <- adespatial::beta.div.comp(mite_beta,coef = "J",quant = F)$repl
}

P91_P <- ls(pattern = "_91_P$")
P91_P <- mget(P91_P)
P91_P <- as.data.frame(P91_P)

# Legendre richness difference index ####
adespatial_92_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_92_P[i] <- adespatial::beta.div.comp(mite_beta,coef = "J",quant = F)$rich
}
BAT_92_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_92_P[i] <- BAT::beta(mite_beta,func = "sorensen",abund = F)$Brich
}

P92_P <- ls(pattern = "_92_P$")
P92_P <- mget(P92_P)
P92_P <- as.data.frame(P92_P)

# Extended Sorensen Similarity ####
adiv_93_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_93_A[i] <- adiv::dsimcom(mite_beta,method = "3",type = "dissimilarity",option = "absolute")
}
diverse_93_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  diverse_93_A[i] <- diverse::dis_entities(t(mite_beta),method = "eDice",category_row = T)[1,2]
}
proxy_93_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_93_A[i] <- proxy::dist(mite_beta,method = "eDice")
}
philentropy_93_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_93_A[i] <- philentropy::dice_dist(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}
proxyC_93_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  proxyC_93_A[i] <- 1-proxyC::simil(as.matrix(mite[i,]),as.matrix(mite[i+1,]),method = "edice")  #0.5108792 
}

P93_A <- ls(pattern = "_93_A$")
P93_A <- mget(P93_A)
P93_A <- as.data.frame(P93_A)

# Species profile distance ####
adespatial_94_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_94_A[i] <- adespatial::beta.div(mite_beta,method = "profiles",save.D = T)$D #0.4505909
}
ade4_94_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_94_A[i] <- ade4::disc(as.data.frame(t(mite_beta)))
}
ClusterR_94_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_94_A[i] <- ClusterR::distance_matrix(mite_beta,method = "Rao_coefficient")[2,1]
}

P94_A <- ls(pattern = "_94_A$")
P94_A <- mget(P94_A)
P94_A <- as.data.frame(P94_A)

# Topsoe ####
ecodive_95_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_95_A[i] <- ecodive::topsoe(mite_beta,rescale = F)
}
philentropy_95_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  philentropy_95_A[i] <- philentropy::topsoe(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P95_A <- ls(pattern = "_95_A$")
P95_A <- mget(P95_A)
P95_A <- as.data.frame(P95_A)

# Yule ####
abdiv_96_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_96_P[i] <- abdiv::yule_dissimilarity(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
proxy_96_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_96_P[i] <- proxy::dist(mite_beta,method = "Yule")
}
MultBiplotR_96_P <- rep(NA,69) # besoin de PA
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  MultBiplotR_96_P[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 15)[1,2]
}
PERMANOVA_96_P <- rep(NA,69) # besoin de PA pour etre correct
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  PERMANOVA_96_P[i] <- 1-PERMANOVA::DistBinary(mite_beta,coefficient = 15,transformation = 1)$D[1,2]
}

P96_P <- ls(pattern = "_96_P$")
P96_P <- mget(P96_P)
P96_P <- as.data.frame(P96_P)

# Yule 2 ####
proxy_97_P<- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_97_P[i] <- proxy::dist(mite_beta,method = "Yule2")
}

P97_P <- ls(pattern = "_97_P$")
P97_P <- mget(P97_P)
P97_P <- as.data.frame(P97_P)


# Wave Hedges distance ####
ecodive_98_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_98_A[i] <- ecodive::wave_hedges(mite_beta,rescale = F) #22.84192
}
proxy_98_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_98_A[i] <- proxy::dist(mite_beta,method = "Wave")
}
EnvNJ_98_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_98_A[i] <- EnvNJ::metrics(t(mite_beta),method = "wavehedges")
}
philentropy_98_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  philentropy_98_A[i] <- philentropy::wave_hedges(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P98_A <- ls(pattern = "_98_A$")
P98_A <- mget(P98_A)
P98_A <- as.data.frame(P98_A)

# Whittaker's index of association ####
wiqid_99_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  wiqid_99_A[i] <- wiqid::distWhittaker(mite[i,],mite[i+1,])
}
proxy_99_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_99_A[i] <- proxy::dist(mite_beta,method = "Whittaker")
}
adespatial_99_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_99_A[i] <- adespatial::beta.div(mite_beta,method = "whittaker",save.D = T)$D
}

P99_A <- ls(pattern = "_99_A$")
P99_A <- mget(P99_A)
P99_A <- as.data.frame(P99_A)

# beta -3 ####
vegan_100_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_100_P[i] <- vegan::betadiver(mite_beta,"-3")
}
P100_P <- ls(pattern = "_100_P$")
P100_P <- mget(P100_P)
P100_P <- as.data.frame(P100_P)

# beta -1 ####
vegan_101_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_101_P[i] <- vegan::betadiver(mite_beta,"-1")
}
P101_P <- ls(pattern = "_101_P$")
P101_P <- mget(P101_P)
P101_P <- as.data.frame(P101_P)

# beta -2 ####
vegan_102_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_102_P[i] <- vegan::betadiver(mite_beta,"-2") 
}
P102_P <- ls(pattern = "_102_P$")
P102_P <- mget(P102_P)
P102_P <- as.data.frame(P102_P)

# beta c ####
vegan_103_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_103_P[i] <- vegan::betadiver(mite_beta,"c")
}
tabula_103_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  tabula_103_P[i] <- tabula::index_cody(as.matrix(mite_beta)) 
}
P103_P <- ls(pattern = "_103_P$")
P103_P <- mget(P103_P)
P103_P <- as.data.frame(P103_P)

# beta cc ####
vegan_104_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_104_P[i] <- vegan::betadiver(mite_beta,"cc")
}
P104_P <- ls(pattern = "_104_P$")
P104_P <- mget(P104_P)
P104_P <- as.data.frame(P104_P)

# beta e ####
vegan_105_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_105_P[i] <- vegan::betadiver(mite_beta,"e")
}
P105_P <- ls(pattern = "_105_P$")
P105_P <- mget(P105_P)
P105_P <- as.data.frame(P105_P)

# beta g ####
vegan_106_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_106_P[i] <- vegan::betadiver(mite_beta,"g")
}

P106_P <- ls(pattern = "_106_P$")
P106_P <- mget(P106_P)
P106_P <- as.data.frame(P106_P)

# beta gl  ####
vegan_107_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_107_P[i] <- vegan::betadiver(mite_beta,"gl") 
}

P107_P <- ls(pattern = "_107_P$")
P107_P <- mget(P107_P)
P107_P <- as.data.frame(P107_P)

# beta hk ####
vegan_108_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_108_P[i] <- vegan::betadiver(mite_beta,"hk") #0.2
}

P108_P <- ls(pattern = "_108_P$")
P108_P <- mget(P108_P)
P108_P <- as.data.frame(P108_P)

# beta l  ####
vegan_109_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_109_P[i] <- vegan::betadiver(mite_beta,"l") 
}

P109_P <- ls(pattern = "_109_P$")
P109_P <- mget(P109_P)
P109_P <- as.data.frame(P109_P)

# beta i  ####
vegan_110_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_110_P[i] <- vegan::betadiver(mite_beta,"I") 
}

P110_P <- ls(pattern = "_110_P$")
P110_P <- mget(P110_P)
P110_P <- as.data.frame(P110_P)

# beta m ####
vegan_111_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_111_P[i] <- vegan::betadiver(mite_beta,"m")
}
P111_P <- ls(pattern = "_111_P$")
P111_P <- mget(P111_P)
P111_P <- as.data.frame(P111_P)

# beta me ####
vegan_112_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_112_P[i] <- vegan::betadiver(mite_beta,"me")
}
P112_P <- ls(pattern = "_112_P$")
P112_P <- mget(P112_P)
P112_P <- as.data.frame(P112_P)

# beta r ####
vegan_113_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_113_P[i] <- vegan::betadiver(mite_beta,"r")
}
tabula_113_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  tabula_113_P[i] <- tabula::index_routledge1(as.matrix(mite_beta)) 
}

P113_P <- ls(pattern = "_113_P$")
P113_P <- mget(P113_P)
P113_P <- as.data.frame(P113_P)

# beta rlb ####
vegan_114_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_114_P[i] <- vegan::betadiver(mite_beta,"rlb")
}
P114_P <- ls(pattern = "_114_P$")
P114_P <- mget(P114_P)
P114_P <- as.data.frame(P114_P)


# beta sim  ####
vegan_115_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_115_P[i] <- vegan::betadiver(mite_beta,"sim") 
}
P115_P <- ls(pattern = "_115_P$")
P115_P <- mget(P115_P)
P115_P <- as.data.frame(P115_P)

# beta t  ####
vegan_116_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_116_P[i] <- vegan::betadiver(mite_beta,"t") 
}
P116_P <- ls(pattern = "_116_P$")
P116_P <- mget(P116_P)
P116_P <- as.data.frame(P116_P)


# beta w ####
vegan_117_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_117_P[i] <- vegan::betadiver(mite_beta,"w")
}
tabula_117_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  tabula_117_P[i] <- tabula::index_whittaker(as.matrix(mite_beta))
}
P117_P <- ls(pattern = "_117_P$")
P117_P <- mget(P117_P)
P117_P <- as.data.frame(P117_P)


# beta wb ####
vegan_118_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_118_P[i] <- vegan::betadiver(mite_beta,"wb")
}
P118_P <- ls(pattern = "_118_P$")
P118_P <- mget(P118_P)
P118_P <- as.data.frame(P118_P)

# beta z ####
vegan_119_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_119_P[i] <- vegan::betadiver(mite_beta,"z") 
}

P119_P <- ls(pattern = "_119_P$")
P119_P <- mget(P119_P)
P119_P <- as.data.frame(P119_P)

# Sans nom ####
vegan_120_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_120_P[i] <- vegan::betadiver(mite_beta,"19")
}
P120_P <- ls(pattern = "_120_P$")
P120_P <- mget(P120_P)
P120_P <- as.data.frame(P120_P)


# Save all the dataset one by one ####
df_exclude <- c("mite","mite_pa","mite_relat","mite_beta","beta_combined","mite_beta")
beta_data <- Filter(is.data.frame,
                     mget(setdiff(ls(.GlobalEnv), df_exclude), envir = .GlobalEnv)
)

Map(function(df, nom) write.csv(df, paste0("data/beta/",paste0("beta_mite_",nom,"_all"), ".csv"), row.names = FALSE),
    beta_data,
    names(beta_data))

# Combine all the dataset in only one
beta_combined <- do.call(cbind, beta_data)
colnames(beta_combined) <- sub(".*\\.", "", colnames(beta_combined))

to_numeric_df <- function(df) {
  df[] <- lapply(df, function(col) as.numeric(as.character(col)))
  return(df)
}
beta_combined <- to_numeric_df(beta_combined)
write.csv(beta_combined,"data/beta/b_combined_all_BCI.csv",row.names = F)

# Same but with a long version
beta_combined$Sample <- paste0(1:69, "-", 2:70)
beta_combined_long <- pivot_longer(beta_combined,cols = colnames(beta_combined)[1:622],names_to = "package_index",values_to = "value")
beta_combined_long <- beta_combined_long |>
  separate(package_index, into = c("package", "index"), sep = "_",extra = "merge")
write.csv(beta_combined_long,"data/beta/b_combined_long_all_mite.csv",row.names = F)

################################################################################

# # Manipulate datasets
# correspondance_codes_names <- read_delim("data/correspondance_codes_names.csv", 
#                                          delim = ";", escape_double = FALSE, trim_ws = TRUE)
# 
# data <- read_csv("data/beta/b_combined_long_all_mite.csv")
# 
# data <- left_join(data,correspondance_codes_names)
# write.csv(data,"data/beta/b_combined_long_all_withnames_mite.csv",row.names = F)

