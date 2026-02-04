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
mite <- mite/rowSums(mite)

# Distance matrix to allow some functions working
mite_dist <- as.matrix(dist(t(mite),method = "euclidean",diag = T,upper = T))
rownames(mite_dist) <- rownames(t(mite))
colnames(mite_dist) <- rownames(t(mite))

mite_pa <- convert_to_presence_absence(mite)

#___________________________ Beta diversity indices ___________________________####

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
  mite_beta <- mite[c(i,i+1),]
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
  mite_beta <- mite[c(i,i+1),]
  proxy_9[i] <- proxy::dist(mite_beta,method = "Bray")
}
labdsv_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  labdsv_9[i] <- labdsv::dsvdis(mite_beta, index = "bray/curtis")
}
PERMANOVA_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_9[i] <- PERMANOVA::DistContinuous(mite_beta,coef = 8)$D[1,2]
}
ClusterR_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_9[i] <- ClusterR::distance_matrix(mite_beta,method = "braycurtis")[2,1]
}
provenance_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  provenance_9[i] <- provenance::bray.diss(mite[i,],mite[i+1,])
}
NST_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_9[i] <- NST::beta.g(mite_beta,dist.method = "bray")
}
adespatial_9 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_9[i] <- adespatial::beta.div(mite_beta,method = "percentdiff",save.D = T)$D
}

# Manque les Sorensen quanti ##

# Canberra ####
stats_10 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  stats_10[i] <- dist(mite_beta,method = "canberra")
}

mgc_10 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  mgc_10[i] <- mgc::mgc.distance(mite_beta,method = "canberra")[1,2]
}
LearnClust_10 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  LearnClust_10[i] <- LearnClust::canberradistance(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
ChemoSpecUtils_10 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ChemoSpecUtils_10[i] <- ChemoSpecUtils::rowDist(as.matrix(mite_beta),method = "canberra")
}
fAssets_10 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fAssets_10[i] <- fAssets::canberraDist(t(mite_beta))
}

diverse_10 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  diverse_10[i] <- diverse::dis_entities(t(mite_beta),method = "Canberra",category_row = T)[1,2]
}
proxy_10 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_10[i] <- proxy::dist(mite_beta,method = "Canberra")
}

BoutrosLab.plotting.general_10 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BoutrosLab.plotting.general_10[i] <- BoutrosLab.plotting.general::dist(mite_beta,method = "canberra")
}

amap_10 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  amap_10[i] <- amap::Dist(mite_beta,method = "canberra")
}
Mercator_10 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Mercator_10[i] <- Mercator::binaryDistance(t(mite_beta),metric = "canberra") 
}


EnvNJ_10 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_10[i] <- EnvNJ::metrics(t(mite_beta),method = "canberra")[1,2]
}
Rlof_10 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rlof_10[i] <- Rlof::distmc(mite_beta,method = "canberra")
}
ClusterR_10 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_10[i] <- ClusterR::distance_matrix(mite_beta,method = "canberra")[2,1]
}
coda.base_10 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  coda.base_10[i] <- coda.base::dist(mite_beta,"canberra")
}
fda.usc_10 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fda.usc_10[i] <- fda.usc::metric.dist(mite_beta,method = "canberra")[1,2]
}
flexclust_10 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  flexclust_10[i] <- flexclust::dist2(mite[i,],mite[i+1,],method = "canberra")
}


# Canberra 2 ####
NST_11 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_11[i] <- NST::beta.g(mite_beta,dist.method = "canberra")
}
adespatial_11 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_11[i] <- adespatial::beta.div(mite_beta,method = "canberra",save.D = T)$D
}
vegan_11 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_11[i] <- vegan::vegdist(mite_beta,method = "canberra")
}
pctax_11 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_11[i] <- pctax::mat_dist(t((mite_beta)),method = "canberra")
}

# Canberra 3 ####
abdiv_12 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_12[i] <- abdiv::canberra(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
dynutils_12 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  dynutils_12[i] <- dynutils::calculate_distance(mite_beta,method = "canberra")[1,2]
}
ecodive_12 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_12[i] <- ecodive::canberra(mite_beta,rescale = F)
}
philentropy_12 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_12[i] <- philentropy::canberra(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
phm_12 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  phm_12[i] <- phm::canberra(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
proxyC_12 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxyC_12[i] <- proxyC::dist(as.matrix(mite_beta),method = "canberra")[1,2]
}

# Cao ####
abdiv_13 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_13[i] <- abdiv::cy_dissimilarity(as.numeric(mite[i,]),as.numeric(mite[i+1,]),base = 10) # base can be modified
}

# Cao 2 ####
vegan_14 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_14[i] <- vegan::vegdist(mite_beta,method = "cao")
}
pctax_14 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_14[i] <- pctax::mat_dist(t((mite_beta)),method = "cao")
}
NST_14 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_14[i] <- NST::beta.g(mite_beta,dist.method = "cao")#0.5159614
}

# Chao-Jaccard à faire ####

# Chao-Ochiai ####
vegan_16 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_16[i] <- vegan::chaodist(mite_beta,method = "1 - sqrt(U*V)")
}

# Chao-Sorensen à faire ####

# Chebyshev ####
ecodive_19 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_19[i] <- ecodive::chebyshev(mite_beta,rescale = F)
}
abdiv_19 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_19[i] <- abdiv::chebyshev(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}
philentropy_19 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_19[i] <- philentropy::chebyshev(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
EnvNJ_19 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_19[i] <- EnvNJ::metrics(t(mite_beta),method = "chebyshev")[1,2]
}
LearnClust_19 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  LearnClust_19[i] <- LearnClust::chebyshevDistance(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
comparator_19 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  comparator_19[i] <- comparator::Chebyshev()(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}
SBCK_19 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  SBCK_19[i] <- print(SBCK::chebyshev(as.matrix(mite[i,]),as.matrix(mite[i+1,]))) 
}
Trading_19 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  Trading_19[i] <- Trading::Chebyshev_distance(mite[i,],mite[i+1,])
}
beadplexr_19 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  beadplexr_19[i] <- beadplexr::dist_chebyshev(mite_beta)
}
ClusterR_19 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_19[i] <- ClusterR::distance_matrix(mite_beta,method = "chebyshev")[2,1]
}
rdist_19 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  rdist_19[i] <- rdist::rdist(mite_beta,metric = "chebyshev")
}

# Chi2 ####
vegan_20 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_20[i] <- vegan::vegdist(mite_beta,method = "chisq")
}
pctax_20 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_20[i] <- pctax::mat_dist(t((mite_beta)),method = "chisq")
}
svs_20 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  svs_20[i] <- svs::dist_chisquare(as.matrix(mite_beta))
}
analogue_20 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  analogue_20[i] <- analogue::distance(mite[i,],mite[i+1,],method = "chi.square")
}
adespatial_20 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_20[i] <- adespatial::beta.div(mite_beta,method = "chisquare",save.D = T)$D
}
SNFtool_20 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_20[i] <- adespatial::beta.div(mite_beta,method = "chisquare",save.D = T)$D
}
proxy_20 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_20[i] <- proxy::dist(mite_beta,method = "Chi-squared")
}
spaa_20 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  spaa_20[i] <- spaa::sp.pair(t(as.matrix(mite_beta)))$chisq
}

# Squared chi square ####
ecodive_21 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_21[i] <- ecodive::squared_chisq(mite_beta,rescale = F)
}
dynutils_21 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  dynutils_21[i] <- dynutils::calculate_distance(mite_beta,method = "chisquared")[1,2] #0
}
EnvNJ_21 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_21[i] <- EnvNJ::metrics(t(mite_beta),method = "squared_chi")
}
analogue_21 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  analogue_21[i] <- analogue::distance(mite[i,],mite[i+1,],method = "SQchi.square")
}
philentropy_21 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_21[i] <- philentropy::squared_chi_sq(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

# Probabilistic Symmetric chi square distance ####
ecodive_22 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_22[i] <- ecodive::psym_chisq(mite_beta,rescale = F)
}

# Chord ####
vegan_23 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_23[i] <- vegan::vegdist(mite_beta,method = "chord")
}
pctax_23 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_23[i] <- pctax::mat_dist(t((mite_beta)),method = "chord")
}
analogue_23 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  analogue_23[i] <- analogue::distance(mite[i,],mite[i+1,],method = "chord")
}
ecodive_23 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_23[i] <- ecodive::chord(mite_beta)
}
adespatial_23 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_23[i] <- adespatial::beta.div(mite_beta,method = "chord",save.D = T)$D
}
proxy_23 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_23[i] <- proxy::dist(mite_beta,method = "Chord")
}
abdiv_23 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_23[i] <- abdiv::chord(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}
wiqid_23 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_23[i] <- wiqid::distChord(mite[i,],mite[i+1,])
}

# Squared chord distance ####
ecodive_24 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_24[i] <- ecodive::squared_chord(mite_beta,rescale = F)
}
analogue_24 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  analogue_24[i] <- analogue::distance(mite[i,],mite[i+1,],method = "SQchord") 
}
philentropy_24 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_24[i] <- philentropy::squared_chord(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
EnvNJ_24 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_24[i] <- EnvNJ::metrics(t(mite_beta),method = "squared_chord")
}

# Log chord distance ####
adespatial_25 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_25[i] <- adespatial::beta.div(mite_beta,method = "log.chord",save.D = T)$D
}

# Clark ####
vegan_26 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_26[i] <- vegan::vegdist(mite_beta,method = "clark") 
}
pctax_26 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_26[i] <- pctax::mat_dist(t((mite_beta)),method = "clark")
}
adespatial_26 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_26[i] <- adespatial::beta.div(mite_beta,method = "divergence",save.D = T)$D
}
abdiv_26 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  abdiv_26[i] <- abdiv::clark_coefficient_of_divergence(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}



# Clark 2 ####
ecodive_27 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_27[i] <- ecodive::clark(mite_beta,rescale = F)
}
philentropy_27 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_27[i] <- philentropy::clark_sq(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

# beta cc ####
vegan_28 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_28[i] <- vegan::betadiver(mite_beta,"cc")
}

# cosine ####
vegan_29 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_29[i] <- vegan::designdist(mite_beta,method = "1-J/sqrt(A*B)",terms = "quadratic")
}
dynutils_29 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  dynutils_29[i] <- dynutils::calculate_distance(mite_beta,method = "cosine")[1,2]
}
SemNeT_29 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  SemNeT_29[i] <- SemNeT::similarity(t(mite_beta),method = "cosine")[1,2]
}
EnvNJ_29 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_29[i] <- EnvNJ::metrics(t(mite_beta),method = "cosine")[1,2]
}
ChemoSpecUtils_29 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ChemoSpecUtils_29[i] <- ChemoSpecUtils::rowDist(as.matrix(mite_beta),method = "cosine")
}
Rfast_29 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_29[i] <- Rfast::Dist(mite_beta,method = "cosine")[1,2]
}
ClusterR_29 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_29[i] <- ClusterR::distance_matrix(mite_beta,method = "cosine")[2,1]
}
svs_29 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  svs_29[i] <- svs::dist_cosine(as.matrix(mite_beta))
}
proxyC_29 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  proxyC_29[i] <- proxyC::simil(as.matrix(mite[i,]),as.matrix(mite[i+1,]),method = "cosine")
}
abdiv_29 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_29[i] <- abdiv::cosine_distance(mite[i,],mite[i+1,])
}
amap_29 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  amap_29[i] <- amap::Dist(mite_beta,method = "pearson")
}

# Absolute Pearson ####
ChemoSpecUtils_30 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ChemoSpecUtils_30[i] <- ChemoSpecUtils::rowDist(as.matrix(mite_beta),method = "abspearson")
}
amap_30 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  amap_30[i] <- amap::Dist(mite_beta,method = "abspearson")
}

# Divergence ####
ecodive_31 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_31[i] <- ecodive::divergence(mite_beta,rescale = F) 
}
proxy_31 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_31[i] <- proxy::dist(mite_beta,method = "divergence")
}
philentropy_31 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  philentropy_31[i] <- philentropy::divergence_sq(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

# Euclidean distance ####
ecodist_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodist_32[i] <- ecodist::distance(mite_beta,method = "euclidean")
}
stats_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  stats_32[i] <- dist(mite_beta,method = "euclidean")
}
vegan_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_32[i] <- vegan::vegdist(mite_beta,method = "euclidean")
}
cluster_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  cluster_32[i] <- cluster::daisy(mite_beta,metric = "euclidean")
}
mgc_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  mgc_32[i] <- mgc::mgc.distance(mite_beta)[1,2]
}
distances_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  distances_32[i] <- distances::distances(mite_beta)[1,2]
}
ldt_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ldt_32[i] <- ldt::s.distance(t(mite_beta),distance = "euclidean")
}
ChemoSpecUtils_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ChemoSpecUtils_32[i] <- ChemoSpecUtils::rowDist(mite_beta,method = "euclidean")
}
fAssets_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fAssets_32[i] <- fAssets::euclideanDist(t(mite_beta))
}
ecodive_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_32[i] <- ecodive::euclidean(mite_beta,rescale = F)
}
proxyC_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxyC_32[i] <- proxyC::dist(as.matrix(mite_beta),method = "euclidean")[2,1]
}
adespatial_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_32[i] <- adespatial::beta.div(mite_beta,method = "euclidean",save.D = T)$D
}
diverse_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  diverse_32[i] <- diverse::dis_entities(t(mite_beta),method = "euclidean",category_row = T)[1,2]
}
pctax_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_32[i] <- pctax::mat_dist(t((mite_beta)),method = "euclidean")
}
proxy_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_32[i] <- proxy::dist(mite_beta,method = "Euclidean") 
}
amap_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  amap_32[i] <- amap::Dist(mite_beta,method = "euclidean") 
}
MultivariateAnalysis_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  MultivariateAnalysis_32[i] <- MultivariateAnalysis::Distancia(mite_beta,Metodo = 1)[1]$Distancia 
}
Mercator_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Mercator_32[i] <- Mercator::binaryDistance(t(mite_beta),metric = "euclid")
}
dynutils_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  dynutils_32[i] <- dynutils::calculate_distance(mite_beta,method = "euclidean")[1,2]
}
fda.usc_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fda.usc_32[i] <- fda.usc::metric.dist(mite_beta,method = "euclidean")[1,2]
}
EnvNJ_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_32[i] <- EnvNJ::metrics(t(mite_beta),method = "euclidean")[1,2]
}
Rfast_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_32[i] <- Rfast::Dist(mite_beta,method = "euclidean")[1,2]
}
ClusterR_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_32[i] <- ClusterR::distance_matrix(mite_beta,method = "euclidean")[2,1]
}
rdist_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  rdist_32[i] <- rdist::rdist(mite_beta,metric = "euclidean")
}
coda.base_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  coda.base_32[i] <- coda.base::dist(mite_beta,"euclidean")
}
NST_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_32[i] <- NST::beta.g(mite_beta,dist.method = "euclidean")
}
Rlof_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rlof_32[i] <- Rlof::distmc(mite_beta,method = "euclidean")
}
BoutrosLab.plotting.general_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BoutrosLab.plotting.general_32[i] <- BoutrosLab.plotting.general::dist(mite_beta,method = "euclidean")
}
fossil_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  fossil_32[i] <- fossil::euclidean(mite[i,],mite[i+1,])
}
abdiv_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_32[i] <- abdiv::euclidean(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
philentropy_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_32[i] <- philentropy::euclidean(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}
comparator_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  comparator_32[i] <- comparator::Euclidean()(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
TSdist_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  TSdist_32[i] <- TSdist::LPDistance(as.numeric(mite[i,]),as.numeric(mite[i+1,]),method = "euclidean")
}
flexclust_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  flexclust_32[i] <- flexclust::dist2(mite[i,],mite[i+1,],method = "euclidean")
}
qkerntool_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  qkerntool_32[i] <- qkerntool::Eucdist(as.matrix(mite[i,]),as.matrix(mite[i+1,]),sEuclidean = T) 
}
codep_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  codep_32[i] <- codep::Euclid(mite[i,],mite[i+1,],squared = F) 
}
ptm_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  ptm_32[i] <- ptm::pairwise.dist(as.matrix(mite[i,]),as.matrix(mite[i+1,]),squared = F)
}
nnspat_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  nnspat_32[i] <- nnspat::euc.dist(mite[i,],mite[i+1,])
}
CEGO_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  CEGO_32[i] <- CEGO::distanceRealEuclidean(mite[i,],mite[i+1,])
}
RnavGraphImageData_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  RnavGraphImageData_32[i] <- RnavGraphImageData::L2Distance(as.matrix(t(mite[i,])),as.matrix(t(mite[i+1,]))) 
}
statisfactory_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  statisfactory_32[i] <- statisfactory::euclid(mite[i,],mite[i+1,])
}
neighbr_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  neighbr_32[i] <- neighbr::distance(as.numeric(mite[i,]),as.numeric(mite[i+1,]),measure = "euclidean")
}
analogue_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  analogue_32[i] <- analogue::distance(mite[i,],mite[i+1,],method = "euclidean")
}
LearnClust_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  LearnClust_32[i] <- LearnClust::edistance(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
comparator_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  comparator_32[i] <- comparator::Euclidean()(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
hmsr_32 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  hmsr_32[i] <- hmsr::euclidean_distance(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

# Anderson's modified Euclidean distance ####
NST_33 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_33[i] <- NST::beta.g(mite_beta,dist.method = "mEuclidean")
}

# Average euclidean distance ####
abdiv_34 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_34[i] <- abdiv::rms_distance(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

# Squared euclidean distance ####
ecodive_35 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_35[i] <- ecodive::squared_euclidean(mite_beta,rescale = F)
}
MultivariateAnalysis_35 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  MultivariateAnalysis_35[i] <- MultivariateAnalysis::Distancia(mite_beta,Metodo = 3)[1]$Distancia
}
neighbr_35 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  neighbr_35[i] <- neighbr::distance(as.numeric(mite[i,]),as.numeric(mite[i+1,]),measure = "squared_euclidean")
}
NMFN_35 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  NMFN_35[i] <- NMFN::distance2(mite[i,],mite[i+1,])
}
ptm_35 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  ptm_35[i] <- ptm::pairwise.dist(as.matrix(mite[i,]),as.matrix(mite[i+1,]),squared = T)
}
codep_35 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  codep_35[i] <- codep::Euclid(mite[i,],mite[i+1,],squared = T)
}
M2SMJF_35 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  M2SMJF_35[i] <- M2SMJF::dist2eu(mite[i,],mite[i+1,])
}
philentropy_35 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_35[i] <- philentropy::squared_euclidean(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
analogue_35 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  analogue_35[i] <- analogue::distance(mite[i,],mite[i+1],method = "SQeuclidean")
}
qkerntool_35 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  qkerntool_35[i] <- qkerntool::Eucdist(as.matrix(mite[i,]),as.matrix(mite[i+1,]),sEuclidean = F)
}
laGP_35 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  laGP_35[i] <- laGP::distance(mite[i,],mite[i+1,])
}
Rfast_35 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_35[i] <- Rfast::Dist(mite_beta,method = "euclidean",square = T)[1,2]
}

# Gower ####
cluster_36 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  cluster_36[i] <- cluster::daisy(mite_beta,metric = "gower")
}
StatMatch_36 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  StatMatch_36[i] <- StatMatch::gower.dist(mite[i,],mite[i+1,])#0.7045455
}
vegan_36 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_36[i] <- vegan::vegdist(mite_beta,method = "gower")
}
pctax_36 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_36[i] <- pctax::mat_dist(t((mite_beta)),method = "gower")
}
ecodive_36 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_36[i] <- ecodive::gower(mite_beta,rescale = F)
}
diverse_36 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  diverse_36[i] <- diverse::dis_entities(t(mite_beta),method = "Gower",category_row = T)[1,2]
}
shipunov_36 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  shipunov_36[i] <- shipunov::Gower.dist(mite[i,],mite[i+1,])
}
proxy_36 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_36[i] <- proxy::dist(mite_beta,method = "Gower") 
}
NST_36 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_36[i] <- NST::beta.g(mite_beta,dist.method = "gower")
}

# Anderson's modified Gower distance ####
NST_37 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_37[i] <- NST::beta.g(mite_beta,dist.method = "mGower")
}

# Gower 2 ####
NST_38 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_38[i] <- NST::beta.g(mite_beta,dist.method = "altGower")
}
vegan_38 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_38[i] <- vegan::vegdist(mite_beta,method = "altGower")
}
pctax_38 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_38[i] <- pctax::mat_dist(t(mite_beta),method = "altGower")
}
philentropy_38 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  philentropy_38[i] <- philentropy::gower(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
Rfast_38 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_38[i] <- Rfast::Dist(mite_beta,method = "gower")[2,1]
}
gower_38 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  gower_38[i] <- gower::gower_dist(mite[i,],mite[i+1,])
}

# Gower 3 ####
ecodist_39 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodist_39[i] <- ecodist::distance(mite_beta,method = "modgower10")
}

# Hamman coefficient ####
ade4_40 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_40[i] <- ade4::dist.binary(mite_beta,method = 6)
}
proxy_40 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_40[i] <- proxy::dist(mite_beta,method = "Hamman")
}
MultivariateAnalysis_40 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  MultivariateAnalysis_40[i] <- MultivariateAnalysis::Distancia(mite_beta,Metodo = 19)[1]$Distancia
}
MultBiplotR_40 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  MultBiplotR_40[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 9)[1,2]
}
PERMANOVA_40 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_40[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 9,transformation = 1)$D[1,2]
}
proxyC_40 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  proxyC_40[i] <- proxyC::simil(as.matrix(mite[i,]),as.matrix(mite[i+1,]),method = "hamann") #0.5 ok
}

# Hamming distance ####
ecodive_41 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_41[i] <- ecodive::hamming(mite_beta) #11
}
abdiv_41 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_41[i] <- abdiv::hamming(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
Rankcluster_41 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  Rankcluster_41[i] <- Rankcluster::distHamming(mite[i,],mite[i+1,])
}
Mercator_41 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Mercator_41[i] <- Mercator::binaryDistance(t(mite_beta),metric = "hamming")
}
pegas_41 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pegas_41[i] <- pegas::dist.hamming(mite_beta)
}
proxyC_41 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxyC_41[i] <- proxyC::dist(as.matrix(mite_beta),method = "hamming")[2,1]
}
bingat_41 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  bingat_41[i] <- bingat::calcDistance(mite[i,],mite[i+1,]) 
}
genMCMCDiag_41 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  genMCMCDiag_41[i] <- genMCMCDiag::hammingDist(mite[i,],mite[i+1,]) #11 
}

# Hamming distance 2 ####
CEGO_42 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  CEGO_42[i] <- CEGO::distanceNumericHamming(mite[i,],mite[i+1,])
}
EnvNJ_42 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_42[i] <- EnvNJ::metrics(t(mite_beta),method = "hamming")[2,1]
}
ClusterR_42 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_42[i] <- ClusterR::distance_matrix(mite_beta,method = "hamming")[2,1]
}
rdist_42 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  rdist_42[i] <- rdist::rdist(mite_beta,metric = "hamming")
}
CEGO_42 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  CEGO_42[i] <- CEGO::distanceNumericHamming(mite[i,],mite[i+1,])
}
EnsCat_42 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnsCat_42[i] <- EnsCat::hammingD(mite_beta)[2,1]
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

# Hellinger ####
vegan_46 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_46[i] <- vegan::vegdist(mite_beta,method = "hellinger")
}
pctax_46 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_46[i] <- pctax::mat_dist(t((mite_beta)),method = "hellinger")
}
ecodive_46 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_46[i] <- ecodive::hellinger(mite_beta,rescale = T)
}
adespatial_46 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_46[i] <- adespatial::beta.div(mite_beta,method = "hellinger",save.D = T)$D
}
proxy_46 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_46[i] <- proxy::dist(mite_beta,method = "Hellinger")
}
Rfast_46 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_46[i] <- Rfast::Dist(mite_beta,method = "hellinger")[2,1]
}
abdiv_46 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_46[i] <- abdiv::hellinger(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

# Jaccard (Abondance) ####
vegan_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_47A[i] <- vegan::vegdist(mite_beta,method = "jaccard",binary = F) # 0.6936661
}
adespatial_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_47A[i] <- adespatial::beta.div(mite_beta,method = "ruzicka",save.D = T)$D
}
abdiv_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_47A[i] <- abdiv::ruzicka(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}
labdsv_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  labdsv_47A[i] <- labdsv::dsvdis(mite_beta, index = "ruzicka")
}
philentropy_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_47A[i] <- abdiv::ruzicka(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}
stylo_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  stylo_47A[i] <- stylo::dist.minmax(mite_beta)
}
adiv_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_47A[i] <- adiv::distMS(mite_beta)
}
ecodive_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_47A[i] <- ecodive::soergel(mite_beta,rescale = F) 
}
proxy_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_47A[i] <- proxy::dist(mite_beta,method = "Soergel")
}
PERMANOVA_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_47A[i] <- PERMANOVA::DistContinuous(mite_beta,coef = 9)$D[2,1]
}
EnvNJ_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_47A[i] <- EnvNJ::metrics(t(mite_beta),method = "soergel")
}
Rfast_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_47A[i] <- Rfast::Dist(mite_beta,method = "soergel")[2,1]
}
NST_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_47A[i] <- NST::beta.g(mite_beta,dist.method = "jaccard")
}
picante_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  picante_47A[i] <- picante::species.dist(t(mite_beta),metric = "jaccard")
}
BoutrosLab.plotting.general_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BoutrosLab.plotting.general_47A[i] <- BoutrosLab.plotting.general::dist(mite_beta,method = "jaccard") 
}
statisfactory_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  statisfactory_47A[i] <- statisfactory::fuzzyJaccard(mite[i,],mite[i+1,]) 
}
geocmeans_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  geocmeans_47A[i] <- geocmeans::calc_jaccard_idx(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
adespatial_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_47A[i] <- adespatial::beta.div.comp(mite_beta,coef = "J",quant = T)$D 
}
pctax_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_47A[i] <- pctax::mat_dist(t((mite_beta)),method = "jaccard")
}
BAT_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_47A[i] <- BAT::beta(mite_beta,func = "jaccard",abund = T)$Btotal
}
prabclus_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  prabclus_47A[i] <- prabclus::jaccard(t(mite_beta)) 
}

# Jaccard P/A ####
labdsv_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  labdsv_47P[i] <- labdsv::dsvdis(mite_beta, index = "steinhaus")
}
ecodist_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodist_47P[i] <- ecodist::distance(mite_beta,method = "jaccard")
}
stats_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  stats_47P[i] <- dist(mite_beta,method = "binary")
}
vegan_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_47P[i] <- vegan::vegdist(mite_beta,method = "jaccard",binary = T)
}
neighbr_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  neighbr_47P[i] <- neighbr::similarity(mite_pa[i,],mite_pa[i+1,],measure = "jaccard") #0.6666667
}
ChemoSpecUtils_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ChemoSpecUtils_47P[i] <- ChemoSpecUtils::rowDist(as.matrix(mite_beta),method = "binary") #0.33333
}
mgc_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  mgc_47P[i] <- mgc::mgc.distance(mite_beta,method = "binary")[1,2] #0.33333
}
BoutrosLab.plotting.general_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BoutrosLab.plotting.general_47P[i] <- BoutrosLab.plotting.general::dist(mite_beta,method = "binary")
}
bioregion_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  bioregion_47P[i] <- bioregion::dissimilarity(as.matrix(mite_beta),metric = "Jaccard")$Jaccard
}
fAssets_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  fAssets_47P[i] <- fAssets::jaccardDist(t(mite_beta))
}
vegan_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_47P[i] <- vegan::betadiver(mite_beta,"j")
}
ecodive_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_47P[i] <- ecodive::jaccard(mite_beta)
}
betapart_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  betapart_47P[i] <- betapart::beta.pair(mite_beta,index.family = "jaccard")$beta.jac 
}
adespatial_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_47P[i] <- adespatial::beta.div.comp(mite_beta,coef = "J",quant = F)$D 
}
tabula_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  tabula_47P[i] <- tabula::similarity(mite_beta,method = "jaccard")
}
adiv_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_47P[i] <- adiv::betastatjac(mite_beta)[1]
}
diverse_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  diverse_47P[i] <-diverse::dis_entities(t(mite_beta),method = "Jaccard",category_row = T)[1,2]
}
BAT_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_47P[i] <- BAT::beta(mite_beta,func = "jaccard",abund = F)$Btotal
}
proxy_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_47P[i] <- proxy::dist(mite_beta,method = "Jaccard")
}
PERMANOVA_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  PERMANOVA_47P[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 3,transformation = 1)$D[1,2]
}
ClusterR_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ClusterR_47P[i] <- ClusterR::distance_matrix(mite_beta,method = "jaccard_coefficient")[2,1]
}
flexclust_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  flexclust_47P[i] <- flexclust::dist2(mite[i,],mite[i+1,],method = "binary")
}
iTOP_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  iTOP_47P[i] <- iTOP::jaccard(mite_pa[i,],mite_pa[i+1,])
}
proxyC_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  proxyC_47P[i] <- proxyC::simil(as.matrix(mite[i,]),as.matrix(mite[i+1,]),method = "jaccard")
}
abdiv_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_47P[i] <- abdiv::jaccard(as.numeric(mite[i,]),as.numeric(mite[i+1,])) #0.3333333
}
wiqid_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_47P[i] <- wiqid::distJaccard(mite[i,],mite[i+1,]) #0.3333333
}
fossil_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  fossil_47P[i] <- fossil::jaccard(mite[i,],mite[i+1,]) #0.6666667
}









# Autres a voir #
?wiqid::distChaoJaccCorr(x1,x1) #0.05012367
wiqid::distChaoJaccNaive(x1,x2) #0.05148794
ConNEcT::funCorrJacc(as.numeric(x1),as.numeric(x2)) #0.9834216
proxyC::simil(as.matrix(x1),as.matrix(x2),method = "ejaccard") #0.3430744
adespatial::beta.div(mite_beta,method = "ab.jaccard",save.D = T,sqrt.D = F) #0.05012367
print(adespatial::dist.ldc(mite_beta,method = "jaccard")) #0.5773503
print(adespatial::dist.ldc(mite_beta,method = "ab.jaccard")) #0.05012367
ade4::dist.binary(x11,method = 1) #0.5773503
CommEcol::dis.chao(mite_betaint,index = "jaccard",version = "rare") #0.05012367
CommEcol::dis.chao(x11,index = "jaccard",version = "probability") #0.05148794
vegan::chaodist(mite_betaint,method = "1 - U*V/(U+V-U*V)") # 0.2098765 
adiv::dsimcom(mite_beta,method = "2",type = "similarity",option = "absolute") #0.6564616
diverse::dis_entities(t(mite_beta),method = "eJaccard",category_row = T) #0.6569256
?proxy::dist(mite_beta,method = "eJaccard") #0.6569256

?MultBiplotR::BinaryDistances(as.matrix(mite_beta_pa),coefficient = 3) #0.5773503 

Mercator::binaryDistance(t(mite_beta),metric = "jaccard") #2.266299

philentropy::jaccard(as.numeric(x1),as.numeric(x2)) #0.6569256

EnvNJ::metrics(t(mite_beta),method = "jaccard") #0.6569256

rdist::rdist(mite_beta,metric = "jaccard") #0.9393939

spaa::sp.pair(t(as.matrix(mite_beta)))$Jaccard #0.7142857

?NST::beta.g(mite_beta,dist.method = "chao.jaccard") #0.0478141 Chao-Jaccard



# Autres
bioregion::dissimilarity(as.matrix(mite[i,]),metric = "Brayturn") #0.5293722





