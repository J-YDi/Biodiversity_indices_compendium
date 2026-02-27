#_______________________________________________________________________________
# Title              : 07_Calculate_beta_indices_transformed_data.r
# Date               : 26/02/2025
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
mite <- mite_relat

#___________________________ Beta diversity indices ___________________________####


# Aitchison ####
vegan_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_2_A[i] <- vegan::vegdist(mite_beta,method = "robust.aitchison")
}

# Put all the values in a single dataframe
P2_A <- ls(pattern = "_2_A$")
P2_A <- mget(P2_A)
P2_A <- as.data.frame(P2_A)

# Average ####
EnvNJ_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_3_A[i] <- EnvNJ::metrics(t(mite_beta),method = "avg")
}

# Put all the values in a single dataframe
P3_A <- ls(pattern = "_3_A$")
P3_A <- mget(P3_A)
P3_A <- as.data.frame(P3_A)

# Bhattacharyya ####
ecodive_4_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_4_A[i] <- ecodive::bhattacharyya(mite_beta,rescale = F)
}

# Put all the values in a single dataframe
P4_A <- ls(pattern = "_4_A$")
P4_A <- mget(P4_A)
P4_A <- as.data.frame(P4_A)


# Binomial ####
vegan_5_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_5_A[i] <- vegan::vegdist(mite_beta,method = "binomial")
}

# Put all the values in a single dataframe
P5_A <- ls(pattern = "_5_A$")
P5_A <- mget(P5_A)
P5_A <- as.data.frame(P5_A)


# Binomial co-occurrence assessment ####
tabula_6_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  tabula_6_A[i] <- tabula::index_binomial(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P6_A <- ls(pattern = "_6_A$")
P6_A <- mget(P6_A)
P6_A <- as.data.frame(P6_A)

# Brainerd-Robinson ####
tabula_7_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  tabula_7_A[i] <- 200-tabula::index_brainerd(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P7_A <- ls(pattern = "_7_A$")
P7_A <- mget(P7_A)
P7_A <- as.data.frame(P7_A)

# Bray-Curtis ####
vegan_9_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_9_A[i] <- vegan::vegdist(mite_beta,method = "bray")
}

P9_A <- ls(pattern = "_9_A$")
P9_A <- mget(P9_A)
P9_A <- as.data.frame(P9_A)

# Canberra ####
stats_10_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  stats_10_A[i] <- dist(mite_beta,method = "canberra")
}

P10_A <- ls(pattern = "_10_A$")
P10_A <- mget(P10_A)
P10_A <- as.data.frame(P10_A)


# Canberra 2 ####
vegan_11_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_11_A[i] <- vegan::vegdist(mite_beta,method = "canberra")
}

P11_A <- ls(pattern = "_11_A$")
P11_A <- mget(P11_A)
P11_A <- as.data.frame(P11_A)

# Canberra 3 ####
abdiv_12_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_12_A[i] <- abdiv::canberra(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P12_A <- ls(pattern = "_12_A$")
P12_A <- mget(P12_A)
P12_A <- as.data.frame(P12_A)

# Cao ####
abdiv_13_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_13_A[i] <- abdiv::cy_dissimilarity(as.numeric(mite[i,]),as.numeric(mite[i+1,]),base = 10) # base can be modified
}

P13_A <- ls(pattern = "_13_A$")
P13_A <- mget(P13_A)
P13_A <- as.data.frame(P13_A)

# Cao 2 ####
pctax_14_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_14_A[i] <- pctax::mat_dist(t((mite_beta)),method = "cao")
}

P14_A <- ls(pattern = "_14_A$")
P14_A <- mget(P14_A)
P14_A <- as.data.frame(P14_A)


# Chao-Jaccard ####

vegan_15_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_15_A[i] <- vegan::vegdist(mite_beta,method = "chao")
}


P15_A <- ls(pattern = "_15_A$")
P15_A <- mget(P15_A)
P15_A <- as.data.frame(P15_A)


# Chao-Ochiai ####
# vegan_16_A <- rep(NA,69)
# for (i in 1:(nrow(mite)-1)){
#   mite_beta <- mite[c(i,i+1),]
#   vegan_16_A[i] <- vegan::chaodist(mite_beta,method = "1 - sqrt(U*V)")
# }
# 
# P16_A <- ls(pattern = "_16_A$")
# P16_A <- mget(P16_A)
# P16_A <- as.data.frame(P16_A)

# Chao-Sorensen ####
# vegan_17_A <- rep(NA,69)
# for (i in 1:(nrow(mite)-1)){
#   mite_beta <- mite[c(i,i+1),]
#   vegan_17_A[i] <- vegan::chaodist(mite_beta,method = "1 - 2*U*V/(U+V)")
# }
# 
# P17_A <- ls(pattern = "_17_A$")
# P17_A <- mget(P17_A)
# P17_A <- as.data.frame(P17_A)

# Chao-Sorensen ####

wiqid_18_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  wiqid_18_A[i] <- wiqid::distChaoSorNaive(mite[i,],mite[i+1,])
} 

P18_A <- ls(pattern = "_18_A$")
P18_A <- mget(P18_A)
P18_A <- as.data.frame(P18_A)

# Chebyshev ####
abdiv_19_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_19_A[i] <- abdiv::chebyshev(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}

P19_A <- ls(pattern = "_19_A$")
P19_A <- mget(P19_A)
P19_A <- as.data.frame(P19_A)

# Chi2 ####
vegan_20_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_20_A[i] <- vegan::vegdist(mite_beta,method = "chisq")
}

P20_A <- ls(pattern = "_20_A$")
P20_A <- mget(P20_A)
P20_A <- as.data.frame(P20_A)

# Squared chi square ####
ecodive_21_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_21_A[i] <- ecodive::squared_chisq(mite_beta,rescale = F)
}

P21_A <- ls(pattern = "_21_A$")
P21_A <- mget(P21_A)
P21_A <- as.data.frame(P21_A) 

# Probabilistic Symmetric chi square distance ####
ecodive_22_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_22_A[i] <- ecodive::psym_chisq(mite_beta,rescale = F)
}

P22_A <- ls(pattern = "_22_A$")
P22_A <- mget(P22_A)
P22_A <- as.data.frame(P22_A) 

# Chord ####
vegan_23_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_23_A[i] <- vegan::vegdist(mite_beta,method = "chord")
}

P23_A <- ls(pattern = "_23_A$")
P23_A <- mget(P23_A)
P23_A <- as.data.frame(P23_A) 

# Squared chord distance ####
analogue_24_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  analogue_24_A[i] <- analogue::distance(mite[i,],mite[i+1,],method = "SQchord") 
}

P24_A <- ls(pattern = "_24_A$")
P24_A <- mget(P24_A)
P24_A <- as.data.frame(P24_A) 

# Log chord distance ####
adespatial_25_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_25_A[i] <- adespatial::beta.div(mite_beta,method = "log.chord",save.D = T)$D
}

P25_A <- ls(pattern = "_25_A$")
P25_A <- mget(P25_A)
P25_A <- as.data.frame(P25_A) 

# Clark ####
vegan_26_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_26_A[i] <- vegan::vegdist(mite_beta,method = "clark") 
}
P26_A <- ls(pattern = "_26_A$")
P26_A <- mget(P26_A)
P26_A <- as.data.frame(P26_A) 

# Clark 2 ####
ecodive_27_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_27_A[i] <- ecodive::clark(mite_beta,rescale = F)
}
P27_A <- ls(pattern = "_27_A$")
P27_A <- mget(P27_A)
P27_A <- as.data.frame(P27_A) 

# Cosine ####
vegan_28_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_28_A[i] <- vegan::designdist(mite_beta,method = "1-J/sqrt(A*B)",terms = "quadratic")
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

P29_A <- ls(pattern = "_29_A$")
P29_A <- mget(P29_A)
P29_A <- as.data.frame(P29_A)

# Euclidean distance ####
vegan_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_30_A[i] <- vegan::vegdist(mite_beta,method = "euclidean")
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

P33_A <- ls(pattern = "_33_A$")
P33_A <- mget(P33_A)
P33_A <- as.data.frame(P33_A)

# Gower ####
vegan_34_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_34_A[i] <- vegan::vegdist(mite_beta,method = "gower")
}

P34_A <- ls(pattern = "_34_A$")
P34_A <- mget(P34_A)
P34_A <- as.data.frame(P34_A)

# Gower 2 ####
vegan_35_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_35_A[i] <- vegan::vegdist(mite_beta,method = "altGower")
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

# Hamming distance ####
proxyC_40_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxyC_40_A[i] <- proxyC::dist(as.matrix(mite_beta),method = "hamming")[2,1]
}

P40_A <- ls(pattern = "_40_A$")
P40_A <- mget(P40_A)
P40_A <- as.data.frame(P40_A)

# Hamming distance 2 ####
EnvNJ_41_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_41_A[i] <- EnvNJ::metrics(t(mite_beta),method = "hamming")[2,1]
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

P42_A <- ls(pattern = "_42_A$")
P42_A <- mget(P42_A)
P42_A <- as.data.frame(P42_A)

# Jaccard (Abondance) ####
vegan_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_43_A[i] <- vegan::vegdist(mite_beta,method = "jaccard",binary = F) # 0.6936661
}
P43_A <- ls(pattern = "_43_A$")
P43_A <- mget(P43_A)
P43_A <- as.data.frame(P43_A)

# Extended Jaccard Similarity ####
vegan_47_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_47_A[i] <- vegan::designdist(mite_beta,method = "(A+B-2*J)/(A+B-J)",terms = "quadratic")
}

P47_A <- ls(pattern = "_47_A$")
P47_A <- mget(P47_A)
P47_A <- as.data.frame(P47_A)


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

P55_A <- ls(pattern = "_55_A$")
P55_A <- mget(P55_A)
P55_A <- as.data.frame(P55_A)

# Kulczynski 2
vegan_57_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_57_A[i] <- vegan::vegdist(mite_beta,method = "kulczynski") 
}

P57_A <- ls(pattern = "_57_A$")
P57_A <- mget(P57_A)
P57_A <- as.data.frame(P57_A)

# Kulczynski 3 ####
EnvNJ_58_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_58_A[i] <- EnvNJ::metrics(t(mite_beta),method = "kulczynski")
}


P58_A <- ls(pattern = "_58_A$")
P58_A <- mget(P58_A)
P58_A <- as.data.frame(P58_A)

# Lorentzian distance ####
EnvNJ_59_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_59_A[i] <- EnvNJ::metrics(t(mite_beta),method = "lorentzian") 
}
P59_A <- ls(pattern = "_59_A$")
P59_A <- mget(P59_A)
P59_A <- as.data.frame(P59_A)

# Manhattan ####
vegan_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_61_A[i] <- vegan::vegdist(mite_beta,method = "manhattan")
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

# Modified mean character difference ####
adespatial_64_A <- rep(NA,69) # Donne la meme valeur que P avec des abondances mais si on donne du PA alors la valeur est differente
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_64_A[i] <- print(adespatial::dist.ldc(mite_beta,method = "modmeanchardiff")) 
}

P64_A <- ls(pattern = "_64_A$")
P64_A <- mget(P64_A)
P64_A <- as.data.frame(P64_A)

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
vegan_67_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_67_A[i] <- vegan::vegdist(mite_beta,method = "morisita")
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

P68_A <- ls(pattern = "_68_A$")
P68_A <- mget(P68_A)
P68_A <- as.data.frame(P68_A)

# Motyka ####
pctax_69_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_69_A[i] <- ecodive::motyka(mite_beta,rescale = F)
}

P69_A <- ls(pattern = "_69_A$")
P69_A <- mget(P69_A)
P69_A <- as.data.frame(P69_A)


# Raup ####
vegan_75_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_75_A[i] <- vegan::vegdist(mite_beta,method = "raup")
}
P75_A <- ls(pattern = "_75_A$")
P75_A <- mget(P75_A)
P75_A <- as.data.frame(P75_A)

# Root mean square ####
abdiv_77_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_77_A[i] <- abdiv::rms_distance(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P77_A <- ls(pattern = "_77_A$")
P77_A <- mget(P77_A)
P77_A <- as.data.frame(P77_A)

# Sokal & Sneath 2 ####
adiv_81_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_81_A[i] <- adiv::dsimcom(mite_beta,method = "1",type = "dissimilarity",option="absolute")
}

P81_A <- ls(pattern = "_81_A$")
P81_A <- mget(P81_A)
P81_A <- as.data.frame(P81_A)

# Extended Sorensen Similarity ####
adiv_93_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_93_A[i] <- adiv::dsimcom(mite_beta,method = "3",type = "dissimilarity",option = "absolute")
}

P93_A <- ls(pattern = "_93_A$")
P93_A <- mget(P93_A)
P93_A <- as.data.frame(P93_A)

# Species profile distance ####
ade4_94_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_94_A[i] <- ade4::disc(as.data.frame(t(mite_beta)))
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
P95_A <- ls(pattern = "_95_A$")
P95_A <- mget(P95_A)
P95_A <- as.data.frame(P95_A)

# Wave Hedges distance ####
ecodive_99_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_99_A[i] <- ecodive::wave_hedges(mite_beta,rescale = F) #22.84192
}
P99_A <- ls(pattern = "_99_A$")
P99_A <- mget(P99_A)
P99_A <- as.data.frame(P99_A)

# Whittaker's index of association ####
proxy_100_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_100_A[i] <- proxy::dist(mite_beta,method = "Whittaker")
}

P100_A <- ls(pattern = "_100_A$")
P100_A <- mget(P100_A)
P100_A <- as.data.frame(P100_A)


# Save all the dataset one by one ####
df_exclude <- c("mite","mite_pa","mite_relat","mite_beta","beta_combined","mite_beta","mite_dist","beta_combined_long")
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
write.csv(beta_combined,"data/beta/b_combined_A_mite_relat.csv",row.names = F)

# Same but with a long version
beta_combined$Sample <- paste0(1:69, "-", 2:70)
beta_combined_long <- pivot_longer(beta_combined,cols = colnames(beta_combined)[1:62],names_to = "package_index",values_to = "value")
beta_combined_long <- beta_combined_long |>
  separate(package_index, into = c("package", "index"), sep = "_",extra = "merge")
write.csv(beta_combined_long,"data/beta/b_combined_long_A_mite_relat.csv",row.names = F)

# Indices A relat #####
#__________________________________Loading data_________________________________####

data <- read_csv("data/beta/b_combined_A_mite_relat.csv")
data_long <- read_csv("data/beta/b_combined_long_A_mite_relat.csv")

#______________________Some basic representations of the data___________________####

# # All
# ggplot(data_long)+
#   geom_line(aes(x=Sample,y=value))+
#   facet_wrap(~index,scale = "free_y")

# Calculate some stats to the plot
data_stats <- data_long |> 
  group_by(index) |> 
  summarise(mean_value = mean(value),
            sd_value = sd(value))

# A indices

levels_index <- data_long$index %>%
  unique() %>%
  .[order(as.numeric(str_extract(., "\\d+")))]

data_long$index <- factor(data_long$index, levels = levels_index)

levels_index <- data_stats$index %>%
  unique() %>%
  .[order(as.numeric(str_extract(., "\\d+")))]

data_stats$index <- factor(data_stats$index, levels = levels_index)

data_long <- data_long %>%
  mutate(
    Sample_num = as.numeric(str_extract(Sample, "^[0-9]+")),
    Sample = fct_reorder(Sample, Sample_num)
  )

ggplot(data_long) +
  geom_segment(aes(x = Sample,y=0, yend = value), col = "royalblue",
               linewidth = 1, alpha = 0.4) +
  geom_point(aes(x = Sample, y = value), col = "royalblue", size = 1.7) +
  geom_label(data = filter(data_stats),
             aes(x = 69/2, y = 0.01, label = paste0(round(mean_value, 3)," +/- ",round(sd_value, 3))),
             color = "black", size = 4,alpha=0.5,linewidth=0) +
  facet_wrap(~ index,scale = "free_y",ncol = 8) +
  labs(x = "Sample", y = "Index value") +
  theme(strip.text = element_text(face = "bold", color = "white",
                                  hjust = 0, size = 10),
        strip.background = element_rect(fill = "royalblue"),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.text.x = element_text(size = 3.7,angle=90,hjust = 1,vjust = 0.5))
ggsave('values_beta_A_mite_relat_indices.png', path = "output/fig/beta/indices/", dpi = 900, width = 600, height = 300, units = 'mm')

# ALL VARIABLES SELECT #############____________________________________________
#________________________________________PCA____________________________________####
colnames(data) <- sub("^[^_]*_", "", colnames(data))
data_pca <- data

PCA_results <- PCA(data_pca)
PCA_results_t <- PCA(t(data_pca))

fviz_screeplot(PCA_results) # Screeplot

# PCA viz with colour arrows

PCA <- fviz_pca_var(PCA_results, axes = c(1, 2), repel = T ,col.var = "royalblue",title="", ggtheme = theme_minimal()) +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )
PCA
ggsave('PCA_beta_A_relat.png', path = "output/fig/beta/indices/", dpi = 1200, width = 250, height = 250, units = 'mm')

fviz_contrib(PCA_results, choice = "var", axes = 1)
fviz_contrib(PCA_results, choice = "var", axes = 2)
fviz_contrib(PCA_results, choice = "var", axes = 3)

corrplot(t(PCA_results$var$contrib),
         is.corr = FALSE,
         method = "pie",col = viridis(200),number.cex = 0.5)

# To check for a spatial dissimilarity we check by see the individuals position by region
fviz_pca_ind(PCA_results_t,addEllipses = F,repel = T,col.ind = "royalblue"
             ,title="", ggtheme = theme_minimal(),legend = "none")

#______________________Clusterings______________________________________________####
# 
data_pca_scaled <- scale(data_pca,T,T)


hc <- hclust(dist(t(data_pca_scaled),method = "euclidean"),method = "ward")
plot(hc)

library(cluster)

cluster_quality(data_pca_scaled, return_table = TRUE)
k=4
clusters <- cutree(hc, k = k)
cluster_cols <- c("red", "blue", "darkgreen", "orange")
label_cols <- cluster_cols[clusters]
dend <- as.dendrogram(hc)
dend <- dendrapply(dend, function(n) {
  if (!is.leaf(n)) attr(n, "height") <- log1p(attr(n, "height"))
  n
})

dend <- color_branches(dend,k = k )
dend <- color_labels(dend,k = k )
dend <- set(dend, "branches_lwd", k)
plot(dend, horiz = T,dLeaf = -0.1,axes=T)


# NMDS 
cluster_quality(data_pca_scaled, return_table = TRUE)
nmds <- metaMDS(dist(t(data_pca_scaled)), k = 4, trymax = 999)

scores_df <- as.data.frame(scores(nmds))  # x,y
scores_df$Sample <- rownames(scores_df)

clusters <- cutree(hclust(dist(t(data_pca_scaled)), method = "ward.D2"), k = 4)
scores_df$Cluster <- factor(clusters)

NMDS <- ggplot(scores_df, aes(x = NMDS1, y = NMDS2, color = Cluster)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample), max.overlaps = Inf, box.padding = 0.5) +
  theme_minimal() +
  geom_label(aes(x=-5,y=-8,label = paste("Stress:",round(nmds$stress,4))),
             color = "black",linewidth = 0)+
  theme(legend.position = "none")+
  labs(title = NULL,
       x = "NMDS1", y = "NMDS2")+
  scale_color_discrete(palette = c("red", "blue", "green2", "orange","magenta"))
NMDS
ggsave('NMDS_k4_beta_A_relat.png', path = "output/fig/beta/indices/", dpi = 1200, width = 250, height = 150, units = 'mm')

# Hellinger transformation ####
#________________________________Loading data___________________________________####
library(vegan)
# We choose the mite data from vegan as the dataset for abundance/count data
data("mite")
mite <- decostand(mite, method = "hellinger")
detach(package:vegan)
# It is possible to work on sipoo data for presence/absence and varespec for abundance/not integer data

# Relative abundances data to allow some functions working
mite_relat <- mite

#___________________________ Beta diversity indices ___________________________####

# Aitchison ####
vegan_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_2_A[i] <- vegan::vegdist(mite_beta,method = "robust.aitchison")
}

# Put all the values in a single dataframe
P2_A <- ls(pattern = "_2_A$")
P2_A <- mget(P2_A)
P2_A <- as.data.frame(P2_A)

# Average ####
EnvNJ_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_3_A[i] <- EnvNJ::metrics(t(mite_beta),method = "avg")
}

# Put all the values in a single dataframe
P3_A <- ls(pattern = "_3_A$")
P3_A <- mget(P3_A)
P3_A <- as.data.frame(P3_A)

# Bhattacharyya ####
ecodive_4_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_4_A[i] <- ecodive::bhattacharyya(mite_beta,rescale = F)
}

# Put all the values in a single dataframe
P4_A <- ls(pattern = "_4_A$")
P4_A <- mget(P4_A)
P4_A <- as.data.frame(P4_A)


# Binomial ####
vegan_5_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_5_A[i] <- vegan::vegdist(mite_beta,method = "binomial")
}

# Put all the values in a single dataframe
P5_A <- ls(pattern = "_5_A$")
P5_A <- mget(P5_A)
P5_A <- as.data.frame(P5_A)


# Binomial co-occurrence assessment ####
tabula_6_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  tabula_6_A[i] <- tabula::index_binomial(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P6_A <- ls(pattern = "_6_A$")
P6_A <- mget(P6_A)
P6_A <- as.data.frame(P6_A)

# Brainerd-Robinson ####
tabula_7_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  tabula_7_A[i] <- 200-tabula::index_brainerd(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P7_A <- ls(pattern = "_7_A$")
P7_A <- mget(P7_A)
P7_A <- as.data.frame(P7_A)

# Bray-Curtis ####
vegan_9_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_9_A[i] <- vegan::vegdist(mite_beta,method = "bray")
}

P9_A <- ls(pattern = "_9_A$")
P9_A <- mget(P9_A)
P9_A <- as.data.frame(P9_A)

# Canberra ####
stats_10_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  stats_10_A[i] <- dist(mite_beta,method = "canberra")
}

P10_A <- ls(pattern = "_10_A$")
P10_A <- mget(P10_A)
P10_A <- as.data.frame(P10_A)


# Canberra 2 ####
vegan_11_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_11_A[i] <- vegan::vegdist(mite_beta,method = "canberra")
}

P11_A <- ls(pattern = "_11_A$")
P11_A <- mget(P11_A)
P11_A <- as.data.frame(P11_A)

# Canberra 3 ####
abdiv_12_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_12_A[i] <- abdiv::canberra(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P12_A <- ls(pattern = "_12_A$")
P12_A <- mget(P12_A)
P12_A <- as.data.frame(P12_A)

# Cao ####
abdiv_13_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_13_A[i] <- abdiv::cy_dissimilarity(as.numeric(mite[i,]),as.numeric(mite[i+1,]),base = 10) # base can be modified
}

P13_A <- ls(pattern = "_13_A$")
P13_A <- mget(P13_A)
P13_A <- as.data.frame(P13_A)

# Cao 2 ####
pctax_14_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_14_A[i] <- pctax::mat_dist(t((mite_beta)),method = "cao")
}

P14_A <- ls(pattern = "_14_A$")
P14_A <- mget(P14_A)
P14_A <- as.data.frame(P14_A)


# Chao-Jaccard ####

vegan_15_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_15_A[i] <- vegan::vegdist(mite_beta,method = "chao")
}


P15_A <- ls(pattern = "_15_A$")
P15_A <- mget(P15_A)
P15_A <- as.data.frame(P15_A)


# Chao-Ochiai ####
# vegan_16_A <- rep(NA,69)
# for (i in 1:(nrow(mite)-1)){
#   mite_beta <- mite[c(i,i+1),]
#   vegan_16_A[i] <- vegan::chaodist(mite_beta,method = "1 - sqrt(U*V)")
# }
#
# P16_A <- ls(pattern = "_16_A$")
# P16_A <- mget(P16_A)
# P16_A <- as.data.frame(P16_A)

# Chao-Sorensen ####
# vegan_17_A <- rep(NA,69)
# for (i in 1:(nrow(mite)-1)){
#   mite_beta <- mite[c(i,i+1),]
#   vegan_17_A[i] <- vegan::chaodist(mite_beta,method = "1 - 2*U*V/(U+V)")
# }
# 
# P17_A <- ls(pattern = "_17_A$")
# P17_A <- mget(P17_A)
# P17_A <- as.data.frame(P17_A)

# Chao-Sorensen ####

wiqid_18_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  wiqid_18_A[i] <- wiqid::distChaoSorNaive(mite[i,],mite[i+1,])
} 

P18_A <- ls(pattern = "_18_A$")
P18_A <- mget(P18_A)
P18_A <- as.data.frame(P18_A)

# Chebyshev ####
abdiv_19_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_19_A[i] <- abdiv::chebyshev(as.numeric(mite[i,]),as.numeric(mite[i+1,])) 
}

P19_A <- ls(pattern = "_19_A$")
P19_A <- mget(P19_A)
P19_A <- as.data.frame(P19_A)

# Chi2 ####
vegan_20_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_20_A[i] <- vegan::vegdist(mite_beta,method = "chisq")
}

P20_A <- ls(pattern = "_20_A$")
P20_A <- mget(P20_A)
P20_A <- as.data.frame(P20_A)

# Squared chi square ####
ecodive_21_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_21_A[i] <- ecodive::squared_chisq(mite_beta,rescale = F)
}

P21_A <- ls(pattern = "_21_A$")
P21_A <- mget(P21_A)
P21_A <- as.data.frame(P21_A) 

# Probabilistic Symmetric chi square distance ####
ecodive_22_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_22_A[i] <- ecodive::psym_chisq(mite_beta,rescale = F)
}

P22_A <- ls(pattern = "_22_A$")
P22_A <- mget(P22_A)
P22_A <- as.data.frame(P22_A) 

# Chord ####
vegan_23_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_23_A[i] <- vegan::vegdist(mite_beta,method = "chord")
}

P23_A <- ls(pattern = "_23_A$")
P23_A <- mget(P23_A)
P23_A <- as.data.frame(P23_A) 

# Squared chord distance ####
analogue_24_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  analogue_24_A[i] <- analogue::distance(mite[i,],mite[i+1,],method = "SQchord") 
}

P24_A <- ls(pattern = "_24_A$")
P24_A <- mget(P24_A)
P24_A <- as.data.frame(P24_A) 

# Log chord distance ####
adespatial_25_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_25_A[i] <- adespatial::beta.div(mite_beta,method = "log.chord",save.D = T)$D
}

P25_A <- ls(pattern = "_25_A$")
P25_A <- mget(P25_A)
P25_A <- as.data.frame(P25_A) 

# Clark ####
vegan_26_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_26_A[i] <- vegan::vegdist(mite_beta,method = "clark") 
}
P26_A <- ls(pattern = "_26_A$")
P26_A <- mget(P26_A)
P26_A <- as.data.frame(P26_A) 

# Clark 2 ####
ecodive_27_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_27_A[i] <- ecodive::clark(mite_beta,rescale = F)
}
P27_A <- ls(pattern = "_27_A$")
P27_A <- mget(P27_A)
P27_A <- as.data.frame(P27_A) 

# Cosine ####
vegan_28_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_28_A[i] <- vegan::designdist(mite_beta,method = "1-J/sqrt(A*B)",terms = "quadratic")
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

P29_A <- ls(pattern = "_29_A$")
P29_A <- mget(P29_A)
P29_A <- as.data.frame(P29_A)

# Euclidean distance ####
vegan_30_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_30_A[i] <- vegan::vegdist(mite_beta,method = "euclidean")
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

P33_A <- ls(pattern = "_33_A$")
P33_A <- mget(P33_A)
P33_A <- as.data.frame(P33_A)

# Gower ####
vegan_34_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_34_A[i] <- vegan::vegdist(mite_beta,method = "gower")
}

P34_A <- ls(pattern = "_34_A$")
P34_A <- mget(P34_A)
P34_A <- as.data.frame(P34_A)

# Gower 2 ####
vegan_35_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_35_A[i] <- vegan::vegdist(mite_beta,method = "altGower")
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

# Hamming distance ####
proxyC_40_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxyC_40_A[i] <- proxyC::dist(as.matrix(mite_beta),method = "hamming")[2,1]
}

P40_A <- ls(pattern = "_40_A$")
P40_A <- mget(P40_A)
P40_A <- as.data.frame(P40_A)

# Hamming distance 2 ####
EnvNJ_41_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_41_A[i] <- EnvNJ::metrics(t(mite_beta),method = "hamming")[2,1]
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

P42_A <- ls(pattern = "_42_A$")
P42_A <- mget(P42_A)
P42_A <- as.data.frame(P42_A)

# Jaccard (Abondance) ####
vegan_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_43_A[i] <- vegan::vegdist(mite_beta,method = "jaccard",binary = F) # 0.6936661
}
P43_A <- ls(pattern = "_43_A$")
P43_A <- mget(P43_A)
P43_A <- as.data.frame(P43_A)

# Extended Jaccard Similarity ####
vegan_47_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_47_A[i] <- vegan::designdist(mite_beta,method = "(A+B-2*J)/(A+B-J)",terms = "quadratic")
}

P47_A <- ls(pattern = "_47_A$")
P47_A <- mget(P47_A)
P47_A <- as.data.frame(P47_A)


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

P55_A <- ls(pattern = "_55_A$")
P55_A <- mget(P55_A)
P55_A <- as.data.frame(P55_A)

# Kulczynski 2
vegan_57_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_57_A[i] <- vegan::vegdist(mite_beta,method = "kulczynski") 
}

P57_A <- ls(pattern = "_57_A$")
P57_A <- mget(P57_A)
P57_A <- as.data.frame(P57_A)

# Kulczynski 3 ####
EnvNJ_58_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_58_A[i] <- EnvNJ::metrics(t(mite_beta),method = "kulczynski")
}


P58_A <- ls(pattern = "_58_A$")
P58_A <- mget(P58_A)
P58_A <- as.data.frame(P58_A)

# Lorentzian distance ####
EnvNJ_59_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_59_A[i] <- EnvNJ::metrics(t(mite_beta),method = "lorentzian") 
}
P59_A <- ls(pattern = "_59_A$")
P59_A <- mget(P59_A)
P59_A <- as.data.frame(P59_A)

# Manhattan ####
vegan_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_61_A[i] <- vegan::vegdist(mite_beta,method = "manhattan")
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

# Modified mean character difference ####
adespatial_64_A <- rep(NA,69) # Donne la meme valeur que P avec des abondances mais si on donne du PA alors la valeur est differente
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_64_A[i] <- print(adespatial::dist.ldc(mite_beta,method = "modmeanchardiff")) 
}

P64_A <- ls(pattern = "_64_A$")
P64_A <- mget(P64_A)
P64_A <- as.data.frame(P64_A)

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
vegan_67_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_67_A[i] <- vegan::vegdist(mite_beta,method = "morisita")
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

P68_A <- ls(pattern = "_68_A$")
P68_A <- mget(P68_A)
P68_A <- as.data.frame(P68_A)

# Motyka ####
pctax_69_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_69_A[i] <- ecodive::motyka(mite_beta,rescale = F)
}

P69_A <- ls(pattern = "_69_A$")
P69_A <- mget(P69_A)
P69_A <- as.data.frame(P69_A)


# Raup ####
vegan_75_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_75_A[i] <- vegan::vegdist(mite_beta,method = "raup")
}
P75_A <- ls(pattern = "_75_A$")
P75_A <- mget(P75_A)
P75_A <- as.data.frame(P75_A)

# Root mean square ####
abdiv_77_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_77_A[i] <- abdiv::rms_distance(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

P77_A <- ls(pattern = "_77_A$")
P77_A <- mget(P77_A)
P77_A <- as.data.frame(P77_A)

# Sokal & Sneath 2 ####
adiv_81_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_81_A[i] <- adiv::dsimcom(mite_beta,method = "1",type = "dissimilarity",option="absolute")
}

P81_A <- ls(pattern = "_81_A$")
P81_A <- mget(P81_A)
P81_A <- as.data.frame(P81_A)

# Extended Sorensen Similarity ####
adiv_93_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adiv_93_A[i] <- adiv::dsimcom(mite_beta,method = "3",type = "dissimilarity",option = "absolute")
}

P93_A <- ls(pattern = "_93_A$")
P93_A <- mget(P93_A)
P93_A <- as.data.frame(P93_A)

# Species profile distance ####
ade4_94_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_94_A[i] <- ade4::disc(as.data.frame(t(mite_beta)))
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
P95_A <- ls(pattern = "_95_A$")
P95_A <- mget(P95_A)
P95_A <- as.data.frame(P95_A)

# Wave Hedges distance ####
ecodive_99_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_99_A[i] <- ecodive::wave_hedges(mite_beta,rescale = F) #22.84192
}
P99_A <- ls(pattern = "_99_A$")
P99_A <- mget(P99_A)
P99_A <- as.data.frame(P99_A)

# Whittaker's index of association ####
proxy_100_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_100_A[i] <- proxy::dist(mite_beta,method = "Whittaker")
}

P100_A <- ls(pattern = "_100_A$")
P100_A <- mget(P100_A)
P100_A <- as.data.frame(P100_A)


# Save all the dataset one by one ####
df_exclude <- c("mite","mite_pa","mite_relat","mite_beta","beta_combined","mite_beta","mite_dist","beta_combined_long")
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
write.csv(beta_combined,"data/beta/b_combined_A_mite_hellinger.csv",row.names = F)

# Same but with a long version
beta_combined$Sample <- paste0(1:69, "-", 2:70)
beta_combined_long <- pivot_longer(beta_combined,cols = colnames(beta_combined)[1:62],names_to = "package_index",values_to = "value")
beta_combined_long <- beta_combined_long |>
  separate(package_index, into = c("package", "index"), sep = "_",extra = "merge")
write.csv(beta_combined_long,"data/beta/b_combined_long_A_mite_hellinger.csv",row.names = F)

# Indices A relat #####
#__________________________________Loading data_________________________________####

data <- read_csv("data/beta/b_combined_A_mite_hellinger.csv")
data_long <- read_csv("data/beta/b_combined_long_A_mite_hellinger.csv")

#______________________Some basic representations of the data___________________####

# # All
# ggplot(data_long)+
#   geom_line(aes(x=Sample,y=value))+
#   facet_wrap(~index,scale = "free_y")

# Calculate some stats to the plot
data_stats <- data_long |> 
  group_by(index) |> 
  summarise(mean_value = mean(value),
            sd_value = sd(value))

# A indices

levels_index <- data_long$index %>%
  unique() %>%
  .[order(as.numeric(str_extract(., "\\d+")))]

data_long$index <- factor(data_long$index, levels = levels_index)

levels_index <- data_stats$index %>%
  unique() %>%
  .[order(as.numeric(str_extract(., "\\d+")))]

data_stats$index <- factor(data_stats$index, levels = levels_index)

data_long <- data_long %>%
  mutate(
    Sample_num = as.numeric(str_extract(Sample, "^[0-9]+")),
    Sample = fct_reorder(Sample, Sample_num)
  )

ggplot(data_long) +
  geom_segment(aes(x = Sample,y=0, yend = value), col = "royalblue",
               linewidth = 1, alpha = 0.4) +
  geom_point(aes(x = Sample, y = value), col = "royalblue", size = 1.7) +
  geom_label(data = filter(data_stats),
             aes(x = 69/2, y = 0.01, label = paste0(round(mean_value, 3)," +/- ",round(sd_value, 3))),
             color = "black", size = 4,alpha=0.5,linewidth=0) +
  facet_wrap(~ index,scale = "free_y",ncol = 8) +
  labs(x = "Sample", y = "Index value") +
  theme(strip.text = element_text(face = "bold", color = "white",
                                  hjust = 0, size = 10),
        strip.background = element_rect(fill = "royalblue"),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.text.x = element_text(size = 3.7,angle=90,hjust = 1,vjust = 0.5))
ggsave('values_beta_A_mite_hellinger_indices.png', path = "output/fig/beta/indices/", dpi = 900, width = 600, height = 300, units = 'mm')

# ALL VARIABLES SELECT #############____________________________________________
#________________________________________PCA____________________________________####
colnames(data) <- sub("^[^_]*_", "", colnames(data))
data_pca <- data

PCA_results <- PCA(data_pca)
PCA_results_t <- PCA(t(data_pca))

fviz_screeplot(PCA_results) # Screeplot

# PCA viz with colour arrows

PCA <- fviz_pca_var(PCA_results, axes = c(1, 2), repel = T ,col.var = "royalblue",title="", ggtheme = theme_minimal()) +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )
PCA
ggsave('PCA_beta_A_hellinger.png', path = "output/fig/beta/indices/", dpi = 1200, width = 250, height = 250, units = 'mm')

fviz_contrib(PCA_results, choice = "var", axes = 1)
fviz_contrib(PCA_results, choice = "var", axes = 2)
fviz_contrib(PCA_results, choice = "var", axes = 3)

corrplot(t(PCA_results$var$contrib),
         is.corr = FALSE,
         method = "pie",col = viridis(200),number.cex = 0.5)

# To check for a spatial dissimilarity we check by see the individuals position by region
fviz_pca_ind(PCA_results_t,addEllipses = F,repel = T,col.ind = "royalblue"
             ,title="", ggtheme = theme_minimal(),legend = "none")

#______________________Clusterings______________________________________________####
# 
data_pca_scaled <- scale(data_pca,T,T)


hc <- hclust(dist(t(data_pca_scaled),method = "euclidean"),method = "ward")
plot(hc)

library(cluster)

cluster_quality(data_pca_scaled, return_table = TRUE)
k=3
clusters <- cutree(hc, k = k)
cluster_cols <- c("red", "blue", "darkgreen", "orange")
label_cols <- cluster_cols[clusters]
dend <- as.dendrogram(hc)
dend <- dendrapply(dend, function(n) {
  if (!is.leaf(n)) attr(n, "height") <- log1p(attr(n, "height"))
  n
})

dend <- color_branches(dend,k = k )
dend <- color_labels(dend,k = k )
dend <- set(dend, "branches_lwd", k)
plot(dend, horiz = T,dLeaf = -0.1,axes=T)


# NMDS 
cluster_quality(data_pca_scaled, return_table = TRUE)
nmds <- metaMDS(dist(t(data_pca_scaled)), k = 3, trymax = 999)

scores_df <- as.data.frame(scores(nmds))  # x,y
scores_df$Sample <- rownames(scores_df)

clusters <- cutree(hclust(dist(t(data_pca_scaled)), method = "ward.D2"), k = 3)
scores_df$Cluster <- factor(clusters)

NMDS <- ggplot(scores_df, aes(x = NMDS1, y = NMDS2, color = Cluster)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample), max.overlaps = Inf, box.padding = 0.5) +
  theme_minimal() +
  geom_label(aes(x=-5,y=-8,label = paste("Stress:",round(nmds$stress,4))),
             color = "black",linewidth = 0)+
  theme(legend.position = "none")+
  labs(title = NULL,
       x = "NMDS1", y = "NMDS2")+
  scale_color_discrete(palette = c("red", "blue", "green2", "orange","magenta"))
NMDS
ggsave('NMDS_k4_beta_A_hellinger.png', path = "output/fig/beta/indices/", dpi = 1200, width = 250, height = 150, units = 'mm')

