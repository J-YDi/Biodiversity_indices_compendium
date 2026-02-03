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
mite <- mite/rowSums(mite)

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

ecodist::distance(x12,method = "euclidean") #40.37368
dist(x12,method = "euclidean") #40.37368
vegan::vegdist(x12,method = "euclidean") # 40.37368

epca::dist.matrix(as.matrix(x1),as.matrix(x2),method = "euclidean") #marche pas contraignant

RobustGaSP::euclidean_distance(as.numeric(x1),as.numeric(x2)) # sortie bizarre

cluster::daisy(x12,metric = "euclidean") # 40.37368

mgc::mgc.distance(x12) #40.37368

distances::distances(x12) #40.37368

analogue::distance(x1,x2,method = "euclidean") # 40.37368

LearnClust::edistance(as.numeric(x1),as.numeric(x2)) #10.96

comparator::Euclidean()(as.numeric(x1),as.numeric(x2)) # 40.37368

SNFtool::dist2(as.numeric(x1),as.numeric(x2)) # non fonctionnel

hmsr::euclidean_distance(as.numeric(x1),as.numeric(x2)) #40.37368

ldt::s.distance(t(x12),distance = "euclidean") # 40.37368

ChemoSpecUtils::rowDist(x12,method = "euclidean") # 40.37368

fAssets::euclideanDist(t(x12)) #40.37368

vegan::vegdist(x12,method = "euclidean",binary = T) # 3.316625

ecodive::euclidean(x12,rescale = F) #40.37368

proxyC::dist(as.matrix(x12),method = "euclidean") #40.37368

adespatial::beta.div(x12,method = "euclidean",save.D = T)$D #40.37368
print(adespatial::dist.ldc(x12,method = "euclidean")) ##40.37368

abdiv::euclidean(as.numeric(x1),as.numeric(x2)) #40.37368

diverse::dis_entities(t(x12),method = "euclidean",category_row = T)[1,2] #40.37368

pctax::mat_dist(t((x12)),method = "euclidean") #40.37368

fossil::euclidean(x1,x2) #40.37368

proxy::dist(x12,method = "Euclidean") #40.37368

amap::Dist(x12,method = "euclidean") #40.37368

MultivariateAnalysis::Distancia(x12,Metodo = 1)[1] #40.37368

Mercator::binaryDistance(t(x12),metric = "euclid") #1

arules::dissimilarity(as.matrix(x12_pa),method = "euclidean") #3.316625 oui avec PA

?dynutils::calculate_distance(x12,method = "euclidean") #40.37368

?philentropy::euclidean(as.numeric(x1),as.numeric(x2)) #40.37368

fda.usc::metric.dist(x12,method = "euclidean") #40.37368

EnvNJ::metrics(t(x12),method = "euclidean") #40.37368

Rfast::Dist(x12,method = "euclidean") #40.37368

comparator::Euclidean()(as.numeric(x1),as.numeric(x2)) #40.37368

ClusterR::distance_matrix(x12,method = "euclidean") #40.37368

rdist::rdist(x12,metric = "euclidean") #40.37368

coda.base::dist(x12,"euclidean") #40.37368

NST::beta.g(x12,dist.method = "euclidean") #40.37368

TSdist::LPDistance(as.numeric(x1),as.numeric(x2),method = "euclidean") #40.37368

Rlof::distmc(x12,method = "euclidean") #40.37368

BoutrosLab.plotting.general::dist(x12,method = "euclidean") #40.37368

flexclust::dist2(x1,x2,method = "euclidean") #40.37368

qkerntool::Eucdist(as.matrix(x1),as.matrix(x2),sEuclidean = T) #40.37368

codep::Euclid(x1,x2,squared = F) #40.37368

ptm::pairwise.dist(as.matrix(x1),as.matrix(x2),squared = F) #40.37368

RprobitB::euc_dist(as.numeric(x1),as.numeric(x2)) #non fonctionnel

freesurferformats::euclidean.dist(as.numeric(x1),as.numeric(x2)) #marche pas

nnspat::euc.dist(x1,x2) #40.37368

CEGO::distanceRealEuclidean(x1,x2) #40.37368

RnavGraphImageData::L2Distance(as.matrix(t(x1)),as.matrix(t(x2))) #40.37368

statisfactory::euclid(x1,x2) #40.37368

neighbr::distance(as.numeric(x1),as.numeric(x2),measure = "euclidean") #40.37368

# Anderson's modified Euclidean distance ####
?NST::beta.g(x12,dist.method = "mEuclidean") #1.223445 modif par Andersen et al. 2006



# Autres
bioregion::dissimilarity(as.matrix(x12),metric = "Brayturn") #0.5293722





