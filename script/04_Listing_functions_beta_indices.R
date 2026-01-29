#_______________________________________________________________________________
# Title              : 04_Listing_functions_beta_indices.r
# Date               : 29/01/2025
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
coda.base_5 <- rep(NA,69)é
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
tabula::index_brainerd(as.numeric(x1),as.numeric(x1)) #93.6067

brsim::brsim(x12)$BR.similarity.matrix












