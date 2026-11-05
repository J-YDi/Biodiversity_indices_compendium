#_______________________________________________________________________________
# Title              : 10_Manual_verif_Beta.r
# Date               : 05/11/2026
# Object             : Script to manually check values from packages that 
#                      calculate beta diversity indices
# Authors            : Jean-Yves Dias
# R version          : 4.5.0
# Github link        : 
#_______________________________________________________________________________

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

#____________________________ Beta diversity indices ___________________________####

# 1P ####
PERMANOVA_1_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  PERMANOVA_1_P[i] <- 1-PERMANOVA::DistBinary(mite_beta,coefficient = 5,transformation = 1)$D[1,2]
}

anderberg_similarity <- function(x, na.rm = TRUE) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x)) || nrow(x) < 2)
    stop("`x` must be a data frame or matrix with at least 2 samples (rows).")
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  
  # Missing values
  if (any(is.na(x))) {
    if (!na.rm) stop("`x` contains missing values.")
    x[is.na(x)] <- 0
  }
  
  # Presence/absence
  pa <- (x > 0) * 1
  
  # a: shared species for each pair
  a <- pa %*% t(pa)
  
  # b + c: species present in only one of the two samples
  S <- rowSums(pa)
  b <- S - a          # species exclusive to sample i (row)
  c <- t(b)           # species exclusive to sample j (column)
  
  P <- a / (a + 2 * (b + c))
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}
custom_1_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  custom_1_P[i] <- 1-anderberg_similarity(mite_beta)
}
sum(abs(PERMANOVA_1_P)-abs(custom_1_P),na.rm = T) 

# 5_P ####
binomial_dissimilarity <- function(x, na.rm = TRUE) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x)) || nrow(x) < 2)
    stop("`x` must be a data frame or matrix with at least 2 samples (rows).")
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  
  # Missing values
  if (any(is.na(x))) {
    if (!na.rm) stop("`x` contains missing values.")
    x[is.na(x)] <- 0
  }
  
  # Presence/absence
  pa <- (x > 0) * 1
  
  # a: shared species for each pair
  a <- pa %*% t(pa)
  
  # b and c: species exclusive to each sample of the pair
  S <- rowSums(pa)
  b <- S - a          # exclusive to sample i (row)
  c <- t(b)           # exclusive to sample j (column)
  
  P <- log(2) * (b + c)
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}
vegan_5_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_5_P[i] <- vegan::vegdist(mite_beta,method = "binomial",binary = T)
}
NST_5_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_5_P[i] <- NST::beta.g(mite_beta,dist.method = "binomial",transform.method = "pa")
}
custom_5_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_5_P[i] <- binomial_dissimilarity(mite_beta)
}
sum(abs(vegan_5_P)-abs(custom_5_P),na.rm = T) 
sum(abs(NST_5_P)-abs(custom_5_P),na.rm = T) 

# 8_P ####
fossil_8_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  fossil_8_P[i] <- 1-fossil::braun.blanquet(mite[i,],mite[i+1,])
}
proxy_8_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_8_P[i] <- proxy::dist(mite_beta,method = "Braun-Blanquet")
}

braun_blanquet_similarity <- function(x, na.rm = TRUE) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x)) || nrow(x) < 2)
    stop("`x` must be a data frame or matrix with at least 2 samples (rows).")
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  
  # Missing values
  if (any(is.na(x))) {
    if (!na.rm) stop("`x` contains missing values.")
    x[is.na(x)] <- 0
  }
  
  # Presence/absence
  pa <- (x > 0) * 1
  
  # a: shared species for each pair
  a <- pa %*% t(pa)
  
  # b and c: species exclusive to each sample of the pair
  S <- rowSums(pa)
  b <- S - a          # exclusive to sample i (row)
  c <- t(b)           # exclusive to sample j (column)
  
  P <- a / (a + pmax(b, c))
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

custom_8_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_8_P[i] <- 1-braun_blanquet_similarity(mite_beta)
}

sum(abs(fossil_8_P)-abs(custom_8_P),na.rm = T) 
sum(abs(proxy_8_P)-abs(custom_8_P),na.rm = T) 

# 38_P ####
ade4_38_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_38_P[i] <- ade4::dist.binary(mite_beta,method = 10)
}
russell_rao_similarity <- function(x, na.rm = TRUE) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x)) || nrow(x) < 2)
    stop("`x` must be a data frame or matrix with at least 2 samples (rows).")
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  
  # Missing values
  if (any(is.na(x))) {
    if (!na.rm) stop("`x` contains missing values.")
    x[is.na(x)] <- 0
  }
  
  # Presence/absence
  pa <- (x > 0) * 1
  
  # a: shared species for each pair
  a <- pa %*% t(pa)
  
  # b, c and d
  S <- rowSums(pa)
  b <- S - a                  # exclusive to sample i (row)
  c <- t(b)                   # exclusive to sample j (column)
  d <- ncol(pa) - a - b - c   # absent from both samples
  
  P <- a / (a + b + c + d)
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}
custom_38_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_38_P[i] <- sqrt(1-russell_rao_similarity(mite_beta))
}
sum(abs(ade4_38_P)-abs(custom_38_P),na.rm = T) 

# 44_P ####
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

sqrt_jaccard_similarity <- function(x, na.rm = TRUE) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x)) || nrow(x) < 2)
    stop("`x` must be a data frame or matrix with at least 2 samples (rows).")
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  
  # Missing values
  if (any(is.na(x))) {
    if (!na.rm) stop("`x` contains missing values.")
    x[is.na(x)] <- 0
  }
  
  # Presence/absence
  pa <- (x > 0) * 1
  
  # a: shared species for each pair
  a <- pa %*% t(pa)
  
  # b and c: species exclusive to each sample of the pair
  S <- rowSums(pa)
  b <- S - a          # exclusive to sample i (row)
  c <- t(b)           # exclusive to sample j (column)
  
  P <- sqrt(1 - (a / (a + b + c)))
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}
custom_44_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_44_P[i] <- sqrt_jaccard_similarity(mite_beta)
}

sum(abs(ade4_44_P)-abs(custom_44_P),na.rm = T) 
sum(abs(adespatial_44_P)-abs(custom_44_P),na.rm = T) 
sum(abs(MultBiplotR_44_P)-abs(custom_44_P),na.rm = T) 

# 48_P ####
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
jaccard_replacement <- function(x, na.rm = TRUE) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x)) || nrow(x) < 2)
    stop("`x` must be a data frame or matrix with at least 2 samples (rows).")
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  
  # Missing values
  if (any(is.na(x))) {
    if (!na.rm) stop("`x` contains missing values.")
    x[is.na(x)] <- 0
  }
  
  # Presence/absence
  pa <- (x > 0) * 1
  
  # a: shared species for each pair
  a <- pa %*% t(pa)
  
  # b and c: species exclusive to each sample of the pair
  S <- rowSums(pa)
  b <- S - a          # exclusive to sample i (row)
  c <- t(b)           # exclusive to sample j (column)
  
  P <- 2 * pmin(b, c) / (a + b + c)
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

custom_48_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_48_P[i] <- jaccard_replacement(mite_beta)
}

sum(abs(BAT_48_P)-abs(custom_48_P),na.rm = T) 
sum(abs(adespatial_48_P)-abs(custom_48_P),na.rm = T) 

# 49_P ####
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

jaccard_richness_difference <- function(x, na.rm = TRUE) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x)) || nrow(x) < 2)
    stop("`x` must be a data frame or matrix with at least 2 samples (rows).")
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  
  # Missing values
  if (any(is.na(x))) {
    if (!na.rm) stop("`x` contains missing values.")
    x[is.na(x)] <- 0
  }
  
  # Presence/absence
  pa <- (x > 0) * 1
  
  # a: shared species for each pair
  a <- pa %*% t(pa)
  
  # b and c: species exclusive to each sample of the pair
  S <- rowSums(pa)
  b <- S - a          # exclusive to sample i (row)
  c <- t(b)           # exclusive to sample j (column)
  
  P <- abs(b - c) / (a + b + c)
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

custom_49_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_49_P[i] <- jaccard_richness_difference(mite_beta)
}

sum(abs(BAT_49_P)-abs(custom_49_P),na.rm = T) 
sum(abs(adespatial_49_P)-abs(custom_49_P),na.rm = T) 

# 50_P ####
adespatial_50_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_50_P[i] <- adespatial::beta.div.comp(mite_beta,coef = "N")$rich 
}
podani_nestedness <- function(x, na.rm = TRUE) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x)) || nrow(x) < 2)
    stop("`x` must be a data frame or matrix with at least 2 samples (rows).")
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  
  # Missing values
  if (any(is.na(x))) {
    if (!na.rm) stop("`x` contains missing values.")
    x[is.na(x)] <- 0
  }
  
  # Presence/absence
  pa <- (x > 0) * 1
  
  # a: shared species for each pair
  a <- pa %*% t(pa)
  
  # b and c: species exclusive to each sample of the pair
  S <- rowSums(pa)
  b <- S - a          # exclusive to sample i (row)
  c <- t(b)           # exclusive to sample j (column)
  
  N <- (a + abs(b - c)) / (a + b + c)
  
  # Nestedness requires shared species
  N[a == 0] <- 0
  N[(a + b + c) == 0] <- NA_real_
  
  dimnames(N) <- list(rownames(x), rownames(x))
  stats::as.dist(N)
}

custom_50_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_50_P[i] <- podani_nestedness(mite_beta)
}

sum(abs(adespatial_50_P)-abs(custom_50_P),na.rm = T) 

# 51_P ####
BAT_51_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_51_P[i] <- BAT::beta(mite_beta,func = "jaccard",abund = F)$Bgain
}

jaccard_gain <- function(x, na.rm = TRUE) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x)) || nrow(x) < 2)
    stop("`x` must be a data frame or matrix with at least 2 samples (rows).")
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  
  # Missing values
  if (any(is.na(x))) {
    if (!na.rm) stop("`x` contains missing values.")
    x[is.na(x)] <- 0
  }
  
  # Presence/absence
  pa <- (x > 0) * 1
  
  # a: shared species for each pair
  a <- pa %*% t(pa)
  
  # excl[r, s]: species present in sample r but not in sample s
  S <- rowSums(pa)
  excl <- S - a
  
  # Lower triangle (r > s): r is the later sample, s the earlier one
  # Gains = species in the later sample only = excl[r, s]
  G <- excl / (a + excl + t(excl))
  G[(a + excl + t(excl)) == 0] <- NA_real_
  
  dimnames(G) <- list(rownames(x), rownames(x))
  stats::as.dist(G)
}

custom_51_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_51_P[i] <- jaccard_gain(mite_beta)
}

sum(abs(BAT_51_P)-abs(custom_51_P),na.rm = T) 
# 52_P ####

BAT_52_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_52_P[i] <- BAT::beta(mite_beta,func = "jaccard",abund = F)$Bloss
}

jaccard_loss <- function(x, na.rm = TRUE) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x)) || nrow(x) < 2)
    stop("`x` must be a data frame or matrix with at least 2 samples (rows).")
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  
  # Missing values
  if (any(is.na(x))) {
    if (!na.rm) stop("`x` contains missing values.")
    x[is.na(x)] <- 0
  }
  
  # Presence/absence
  pa <- (x > 0) * 1
  
  # a: shared species for each pair
  a <- pa %*% t(pa)
  
  # excl[r, s]: species present in sample r but not in sample s
  S <- rowSums(pa)
  excl <- S - a
  
  # Lower triangle (r > s): r is the later sample, s the earlier one
  # Losses = species in the earlier sample only = excl[s, r] = t(excl)[r, s]
  L <- t(excl) / (a + excl + t(excl))
  L[(a + excl + t(excl)) == 0] <- NA_real_
  
  dimnames(L) <- list(rownames(x), rownames(x))
  stats::as.dist(L)
}

custom_52_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_52_P[i] <- jaccard_loss(mite_beta)
}

sum(abs(BAT_52_P)-abs(custom_52_P),na.rm = T) 

# 57_P ####
proxy_57_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  proxy_57_P[i] <- proxy::dist(mite_beta,method = "Kulczynski1") 
}
abdiv_57_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  abdiv_57_P[i] <- abdiv::kulczynski_first(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}


