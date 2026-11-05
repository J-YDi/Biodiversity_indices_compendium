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

kulczynski1_similarity <- function(x, na.rm = TRUE) {
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
  
  P <- a / (b + c)                     # a / 0 = Inf for identical samples
  P[(a + b + c) == 0] <- NA_real_      # empty pairs
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

custom_57_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_57_P[i] <- kulczynski1_similarity(mite_beta)
}

# 63_P ####
abdiv_63_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  abdiv_63_P[i] <- abdiv::mean_character_difference(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
# abdiv do not convert to P/A -> wrong
mean_character_difference <- function(x, na.rm = TRUE) {
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
  
  P <- 1 - (a + d) / (a + b + c + d)
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

custom_63_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_63_P[i] <- mean_character_difference(mite_beta)
}

abdiv::mean_character_difference(as.numeric(mite[i,]),as.numeric(mite[i+1,]))

# 70_P ####
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
  mite_beta <- mite_pa[c(i,i+1),]
  proxy_70_P[i] <- proxy::dist(mite_beta,method = "Mountford")
}
NST_70_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_70_P[i] <- NST::beta.g(mite_beta,dist.method = "mountford")
}

mountford_similarity <- function(x, na.rm = TRUE) {
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
  
  # Richness of each sample in the pair
  nA <- a + b
  nB <- a + c
  
  P <- 2 * a / (2 * nA * nB - (nA + nB) * a)
  P[a == 0] <- 0                       # no shared species
  P[(a + b + c) == 0] <- NA_real_      # empty pairs
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
} # arrondi

custom_70_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_70_P[i] <- 1-mountford_similarity(mite_beta)
}

sum(abs(vegan_70_P)-abs(custom_70_P),na.rm = T)
sum(abs(pctax_70_P)-abs(custom_70_P),na.rm = T)
sum(abs(proxy_70_P)-abs(custom_70_P),na.rm = T)
sum(abs(NST_70_P)-abs(custom_70_P),na.rm = T)

mountford_dissimilarity <- function(x, na.rm = TRUE) {
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
  n  <- nrow(pa)
  
  # a, b, c for all pairs
  a_mat <- pa %*% t(pa)          # shared species
  S     <- rowSums(pa)
  b_mat <- S - a_mat             # b_mat[i, j]: exclusive to sample i
  
  P <- matrix(NA_real_, n, n)
  
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      a <- a_mat[i, j]
      b <- b_mat[i, j]           # exclusive to sample i
      c <- b_mat[j, i]           # exclusive to sample j
      
      if (a + b + c == 0) next                  # empty pair: NA
      if (a == 0) { P[i, j] <- log(2); next }    # no shared species
      if (b == 0 || c == 0) { P[i, j] <- 0; next } # nested or identical
      
      # exp((a+b)M) + exp((a+c)M) = 1 + exp((a+b+c)M),
      # divided by exp((a+b+c)M)
      g <- function(M) {
        exp(-c * M) + exp(-b * M) - exp(-(a + b + c) * M) - 1
      }
      
      # g > 0 near 0 and g -> -1 for large M: search for the sign change
      upper <- 1
      while (g(upper) > 0) upper <- upper * 2
      
      M <- stats::uniroot(g, lower = 1e-10, upper = upper,
                          tol = .Machine$double.eps^0.5)$root
      
      P[i, j] <- log(2) - M
    }
  }
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
} # exact solution

custom_70_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_70_P[i] <- mountford_dissimilarity(mite_beta)
}

sum(abs(vegan_70_P)-abs(custom_70_P),na.rm = T)
sum(abs(pctax_70_P)-abs(custom_70_P),na.rm = T)
sum(abs(proxy_70_P)-abs(custom_70_P),na.rm = T)
sum(abs(NST_70_P)-abs(custom_70_P),na.rm = T)

# 73_P ####
ade4_73_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  ade4_73_P[i] <- ade4::dist.binary(mite_beta,method = 9)
}
proxy_73_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_73_P[i] <- proxy::simil(mite_beta,method = "phi")
}
MultBiplotR_73_P <- rep(NA,69) # Besoin de PA pour marcher
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  MultBiplotR_73_P[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 14,transformation = "Identity")[1,2]
}
PERMANOVA_73_P <- rep(NA,69) # Besoin de PA pour etre correct
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  PERMANOVA_73_P[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 14,transformation = 1)$D[1,2]
}

phi_similarity <- function(x, na.rm = TRUE) {
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
  
  denom <- sqrt((a + b) * (a + c) * (b + d) * (c + d))
  
  P <- (a * d - b * c) / denom
  P[denom == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}
custom_73_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_73_P[i] <- phi_similarity(mite_beta)
}

sum(abs(ade4_73_P)-abs(custom_73_P),na.rm = T)
sum(abs(MultBiplotR_73_P)-abs(custom_73_P),na.rm = T)
sum(abs(proxy_73_P)-abs(custom_73_P),na.rm = T)
sum(abs(PERMANOVA_73_P)-abs(custom_73_P),na.rm = T)

# 74_P ####
wiqid_74_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  wiqid_74_P[i] <- wiqid::distPreston(mite[i,],mite[i+1,])
}
preston_dissimilarity <- function(x, na.rm = TRUE) {
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
  n  <- nrow(pa)
  
  # a, b, c for all pairs
  a_mat <- pa %*% t(pa)          # shared species
  S     <- rowSums(pa)
  b_mat <- S - a_mat             # b_mat[i, j]: exclusive to sample i
  
  P <- matrix(NA_real_, n, n)
  
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      a <- a_mat[i, j]
      b <- b_mat[i, j]           # exclusive to sample i
      c <- b_mat[j, i]           # exclusive to sample j
      
      if (a + b + c == 0) next                  # empty pair: NA
      if (a == 0) { P[i, j] <- 1; next }         # no shared species
      if (b == 0 || c == 0) { P[i, j] <- 0; next } # nested or identical
      
      u <- (a + b) / (a + b + c)
      v <- (a + c) / (a + b + c)
      
      # Solve u^t + v^t = 1 for t = 1/P (decreasing in t)
      f <- function(t) u^t + v^t - 1
      upper <- 2
      while (f(upper) > 0) upper <- upper * 2
      
      t <- stats::uniroot(f, lower = 1, upper = upper,
                          tol = .Machine$double.eps^0.5)$root
      
      P[i, j] <- 1 / t
    }
  }
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}
custom_74_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_74_P[i] <- preston_dissimilarity(mite_beta)
}

sum(abs(wiqid_74_P)-abs(custom_74_P),na.rm = T)

# 80_P ####
abdiv_80_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_80_P[i] <- abdiv::sokal_sneath(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

sokal_sneath_dissimilarity <- function(x, na.rm = TRUE) {
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
  
  P <- (2 * b + 2 * c) / (a + 2 * b + 2 * c)
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

custom_80_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_80_P[i] <- sokal_sneath_dissimilarity(mite_beta)
}

sum(abs(abdiv_80_P)-abs(custom_80_P),na.rm = T)

# 82_P ####
sokal_sneath_similarity <- function(x, na.rm = TRUE) {
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
  
  P <- (a + d) / (a + (1 / 2) * (b + c) + d)
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

custom_82_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_82_P[i] <- sokal_sneath_similarity(mite_beta)
}

MultivariateAnalysis_82_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  MultivariateAnalysis_82_P[i] <- 1-MultivariateAnalysis::Distancia(mite_beta,Metodo = 14)[1]$Distancia
}
MultBiplotR_82_P <- rep(NA,69) # Besoin de PA
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  MultBiplotR_82_P[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 8,transformation = "Identity")[1,2]
}
PERMANOVA_82_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  PERMANOVA_82_P[i] <- PERMANOVA::DistBinary(mite_beta,coefficient = 8,transformation = 1)$D[2,1]
}
Mercator_82_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  Mercator_82_P[i] <- Mercator::binaryDistance(t(mite_beta),metric = "sokalMichener")
}

sum(abs(MultivariateAnalysis_82_P)-abs(custom_82_P),na.rm = T)
sum(abs(MultBiplotR_82_P)-abs(custom_82_P),na.rm = T)
sum(abs(PERMANOVA_82_P)-abs(custom_82_P),na.rm = T)
sum(abs(Mercator_82_P)-abs(custom_82_P),na.rm = T)

# 83_P ####
ade4_83_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_83_P[i] <- ade4::dist.binary(mite_beta,method = 3)^2
}

custom_1_P

# 84_P ####
ade4_84_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ade4_84_P[i] <- ade4::dist.binary(mite_beta,method = 8)^2
}
MultBiplotR_84_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  MultBiplotR_84_P[i] <- MultBiplotR::BinaryDistances(as.matrix(mite_beta),coefficient = 13)[1,2]^2
}
PERMANOVA_84_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  PERMANOVA_84_P[i] <- 1-PERMANOVA::DistBinary(mite_beta,coefficient = 13,transformation = 1)$D[1,2]
}
sokal_sneath13_similarity <- function(x, na.rm = TRUE) {
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
  
  denom <- sqrt((a + b) * (a + c) * (d + b) * (d + c))
  
  P <- (a * d) / denom
  P[denom == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

custom_84_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_84_P[i] <- 1-sokal_sneath13_similarity(mite_beta)
}

sum(abs(ade4_84_P)-abs(custom_84_P),na.rm = T)
sum(abs(MultBiplotR_84_P)-abs(custom_84_P),na.rm = T)
sum(abs(PERMANOVA_84_P)-abs(custom_84_P),na.rm = T)

# 86_P ####
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
sqrt_sorensen_distance <- function(x, na.rm = TRUE) {
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
  
  P <- sqrt((b + c) / (2 * a + b + c))
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

custom_86_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_86_P[i] <- sqrt_sorensen_distance(mite_beta)
}

sum(abs(adespatial_86_P)-abs(custom_86_P),na.rm = T)
sum(abs(MultBiplotR_86_P)-abs(custom_86_P),na.rm = T)
sum(abs(abdiv_86_P)-abs(custom_86_P),na.rm = T)
