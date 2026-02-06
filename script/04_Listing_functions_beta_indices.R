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

mite <- convert_to_presence_absence(mite)

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
  NST_14[i] <- NST::beta.g(mite_beta,dist.method = "cao")
}

# Chao-Jaccard à faire ####
NST_15 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_15[i] <- NST::beta.g(mite_beta,dist.method = "chao")
}
CommEcol_15 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  CommEcol_15[i] <- CommEcol::dis.chao(mite_beta,index = "jaccard",version = "rare") #0.05012367
}
adespatial_15 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_15[i] <- adespatial::beta.div(mite_beta,method = "ab.jaccard",save.D = T,sqrt.D = F)$D
}
pctax_15 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_15[i] <- pctax::mat_dist(t((mite_beta)),method = "chao") #0.05012367
}
vegan_15 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_15[i] <- vegan::vegdist(mite_beta,method = "chao")
}
wiqid_15 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_15[i] <- wiqid::distChaoJaccCorr(mite[i,],mite[i+1,])
}

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
  mite_beta <- mite[c(i,i+1),]
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
Mercator_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Mercator_47A[i] <- Mercator::binaryDistance(t(mite_beta),metric = "jaccard") #2.266299
}
rdist_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  rdist_47A[i] <- rdist::rdist(mite_beta,metric = "jaccard") 
}
ConNEct_47A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  ConNEct_47A[i] <- ConNEcT::funCorrJacc(as.numeric(mite[i,]),as.numeric(mite[i+1,]))$value #0.9834216
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
  neighbr_47P[i] <- neighbr::similarity(mite[i,],mite[i+1,],measure = "jaccard") 
}
ChemoSpecUtils_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ChemoSpecUtils_47P[i] <- ChemoSpecUtils::rowDist(as.matrix(mite_beta),method = "binary") 
}
mgc_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  mgc_47P[i] <- mgc::mgc.distance(mite_beta,method = "binary")[1,2] 
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
  iTOP_47P[i] <- iTOP::jaccard(mite[i,],mite[i+1,])
}
proxyC_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  proxyC_47P[i] <- proxyC::simil(as.matrix(mite[i,]),as.matrix(mite[i+1,]),method = "jaccard")
}
abdiv_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_47P[i] <- abdiv::jaccard(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
wiqid_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_47P[i] <- wiqid::distJaccard(mite[i,],mite[i+1,])
}
fossil_47P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  fossil_47P[i] <- fossil::jaccard(mite[i,],mite[i+1,])
}

# Squared-root Jaccard ####
adespatial_48 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_48[i] <- print(adespatial::dist.ldc(mite_beta,method = "jaccard"))
}
MultBiplotR_48 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  MultBiplotR_48[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 3)[1,2] 
}
ade4_48 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_48[i] <- ade4::dist.binary(mite_beta,method = 1) 
}

# Turnover component of Jaccard dissimilarity defined by Baselga ####
vegan_49 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_49[i] <- vegan::nestedbetajac(mite_beta)[1]
}
betapart_49 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  betapart_49[i] <- betapart::beta.pair(mite_beta,index.family = "jaccard")$beta.jtu
}
adespatial_49 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_49[i] <- adespatial::beta.div.comp(mite_beta,coef = "BJ")$repl 
}
bioregion_49 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  bioregion_49[i] <- bioregion::dissimilarity(as.matrix(mite_beta),metric = "Jaccardturn")$Jaccardturn
}
abdiv_49 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
    abdiv_49[i] <- abdiv::jaccard_turnover(as.numeric(mite[i,]),as.numeric(mite[i+1,])) #0.2666667
}

# Nestedness component of Jaccard dissimilarity defined by Baselga ####
abdiv_50 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_50[i] <- abdiv::jaccard_nestedness(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
vegan_50 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_50[i] <- vegan::nestedbetajac(mite_beta)[2] 
}
betapart_50 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  betapart_50[i] <- betapart::beta.pair(mite_beta,index.family = "jaccard")$beta.jne
}
adespatial_50 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_50[i] <- adespatial::beta.div.comp(mite_beta,coef = "BJ")$rich
}
# Extended Jaccard Similarity ####
adiv_51 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_51[i] <- adiv::dsimcom(mite_beta,method = "2",type = "similarity",option = "absolute")[1,2] 
}
diverse_51 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  diverse_51[i] <- diverse::dis_entities(t(mite_beta),method = "eJaccard",category_row = T)[1,2] 
}
proxy_51 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_51[i] <- proxy::dist(mite_beta,method = "eJaccard")
}
EnvNJ_51 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_51[i] <- EnvNJ::metrics(t(mite_beta),method = "jaccard")
}
philentropy_51 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_51[i] <- philentropy::jaccard(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
proxyC_51 <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  proxyC_51[i] <- proxyC::simil(as.matrix(mite[i,]),as.matrix(mite[i+1,]),method = "ejaccard") 
}

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

#### Minkowski ####-------------------------------------------------------------
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




