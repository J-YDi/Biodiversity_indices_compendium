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

# 89_P ####

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

sorensen_gain <- function(x, na.rm = TRUE) {
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
  denom <- 2 * a + excl + t(excl)
  G <- excl / denom
  G[denom == 0] <- NA_real_
  
  dimnames(G) <- list(rownames(x), rownames(x))
  stats::as.dist(G)
}

custom_89_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_89_P[i] <- sorensen_gain(mite_beta)
}

sum(abs(BAT_89_P)-abs(custom_89_P),na.rm = T)
sum(abs(adiv_89_P)-abs(custom_89_P),na.rm = T)

# 90_P ####
sorensen_loss <- function(x, na.rm = TRUE) {
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
  denom <- 2 * a + excl + t(excl)
  L <- t(excl) / denom
  L[denom == 0] <- NA_real_
  
  dimnames(L) <- list(rownames(x), rownames(x))
  stats::as.dist(L)
}

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

custom_90_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_90_P[i] <- sorensen_loss(mite_beta)
}

sum(abs(BAT_90_P)-abs(custom_90_P),na.rm = T)
sum(abs(adiv_90_P)-abs(custom_90_P),na.rm = T)

# 91_P ####
sorensen_replacement <- function(x, na.rm = TRUE) {
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
  
  P <- 2 * pmin(b, c) / (2 * a + b + c)
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

adespatial_91_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_91_P[i] <-  adespatial::beta.div.comp(mite_beta,coef = "S")$repl
}
custom_91_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_91_P[i] <- sorensen_replacement(mite_beta)
}

sum(abs(adespatial_91_P)-abs(custom_91_P),na.rm = T)

# 92_P ####
adespatial_92_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_92_P[i] <- adespatial::beta.div.comp(mite_beta,coef = "S",quant = F)$rich
}
BAT_92_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  BAT_92_P[i] <- BAT::beta(mite_beta,func = "sorensen",abund = F)$Brich
}
sorensen_richness_difference <- function(x, na.rm = TRUE) {
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
  
  P <- abs(b - c) / (2 * a + b + c)
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

custom_92_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_92_P[i] <- sorensen_richness_difference(mite_beta)
}

sum(abs(adespatial_92_P)-abs(custom_92_P),na.rm = T)
sum(abs(BAT_92_P)-abs(custom_92_P),na.rm = T)

# 96_P ####

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

PERMANOVA_96_P <- rep(NA,69) # besoin de PA pour etre correct
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  PERMANOVA_96_P[i] <- 1-PERMANOVA::DistBinary(mite_beta,coefficient = 15,transformation = 1)$D[1,2]
}

yule_dissimilarity <- function(x, na.rm = TRUE) {
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
  
  denom <- a * d + b * c
  
  P <- 2 * b * c / denom
  P[denom == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

custom_96_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_96_P[i] <- yule_dissimilarity(mite_beta)
}

sum(abs(abdiv_96_P)-abs(custom_96_P),na.rm = T)
sum(abs(proxy_96_P)-abs(custom_96_P),na.rm = T)
sum(abs(PERMANOVA_96_P)-abs(custom_96_P),na.rm = T)

# 97_P ####
yule_similarity <- function(x, na.rm = TRUE) {
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
  
  sad <- sqrt(a * d)
  sbc <- sqrt(b * c)
  denom <- sad + sbc
  
  P <- (sad - sbc) / denom
  P[denom == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

proxy_97_P<- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_97_P[i] <- proxy::dist(mite_beta,method = "Yule2")
}

custom_97_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_97_P[i] <- yule_similarity(mite_beta)
}

sum(abs(proxy_97_P)-abs(1-custom_97_P),na.rm = T)

# 98_P ####

yule_dissimilarity <- function(x, na.rm = TRUE) {
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
  
  denom <- a*d + b*c
  
  P <- (a * d - b * c) / denom
  P[denom == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

MultBiplotR_98_P <- rep(NA,69) # besoin de PA
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_pa[c(i,i+1),]
  pa <- (as.matrix(mite_beta) > 0) * 1
  MultBiplotR_98_P[i] <- MultBiplotR::BinaryDistances(pa, coefficient = 15, transformation = "Identity")[1,2]
}

custom_98_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_98_P[i] <- yule_dissimilarity(mite_beta)
}

sum(abs(MultBiplotR_98_P)-abs(custom_98_P),na.rm = T)

# 101_P ####
vegan_101_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_101_P[i] <- vegan::betadiver(mite_beta,"-3")
}

williams_beta3 <- function(x, na.rm = TRUE) {
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
  
  P <- pmin(b, c) / (a + b + c)
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

custom_101_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_101_P[i] <- williams_beta3(mite_beta)
}

sum(abs(vegan_101_P)-abs(custom_101_P),na.rm = T)

# 118_P ####
vegan_118_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_118_P[i] <- vegan::betadiver(mite_beta,"w")
}
tabula_118_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  tabula_118_P[i] <- tabula::index_whittaker(as.matrix(mite_beta))
}

whittaker_beta_w <- function(x, na.rm = TRUE) {
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
  
  P <- (a + b + c) / ((2 * a + b + c) / 2)
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

custom_118_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_118_P[i] <- whittaker_beta_w(mite_beta)
}

sum(abs(vegan_118_P)-abs(1-custom_118_P),na.rm = T)
sum(abs(tabula_118_P)-abs(1-custom_118_P),na.rm = T)

# 102_P ####
vegan_102_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_102_P[i] <- vegan::betadiver(mite_beta,"-1")
}

whittaker_beta_w <- function(x, na.rm = TRUE) {
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
  
  P <- (a + b + c) / ((2 * a + b + c) / 2)
  P <- P - 1
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

custom_102_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_102_P[i] <- whittaker_beta_w(mite_beta)
}

sum(abs(vegan_102_P)-abs(custom_102_P),na.rm = T)
sum(abs(vegan_118_P)-abs(custom_118_P),na.rm = T)

# 103_P ####

harrison_beta_2 <- function(x, na.rm = TRUE) {
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
  
  P <- pmin(b, c) / (pmax(b,c) + a)
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

vegan_103_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_103_P[i] <- vegan::betadiver(mite_beta,"-2") 
}

custom_103_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_103_P[i] <- harrison_beta_2(mite_beta)
}

sum(abs(vegan_103_P)-abs(custom_103_P),na.rm = T)

# 104_P ####
vegan_104_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_104_P[i] <- vegan::betadiver(mite_beta,"c")
}
tabula_104_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  tabula_104_P[i] <- tabula::index_cody(as.matrix(mite_beta)) 
}

cody_beta_c <- function(x, na.rm = TRUE) {
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
  
  P <- (b + c )/2
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

custom_104_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_104_P[i] <- cody_beta_c(mite_beta)
}

sum(abs(vegan_104_P)-abs(custom_104_P),na.rm = T)
sum(abs(tabula_104_P)-abs(custom_104_P),na.rm = T)

# 105_P ####
vegan_105_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_105_P[i] <- vegan::betadiver(mite_beta,"cc")
}

colwell_beta_cc <- function(x, na.rm = TRUE) {
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
  
  P <- (b + c )/(a + b + c)
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

custom_105_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_105_P[i] <- colwell_beta_cc(mite_beta)
}

sum(abs(vegan_105_P)-abs(custom_105_P),na.rm = T)

# 106_P ####
vegan_106_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_106_P[i] <- vegan::betadiver(mite_beta,"e")
}

routledge_beta_e <- function(x, na.rm = TRUE) {
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
  
  # x * log(x), with 0 * log(0) = 0
  xlogx <- function(z) ifelse(z > 0, z * log(z), 0)
  
  tot <- 2 * a + b + c
  
  P <- exp(log(tot) - (1 / tot) * 2 * a * log(2) -
             (1 / tot) * (xlogx(a + b) + xlogx(a + c))) - 1
  
  # Remove rounding errors around 0 (identical samples)
  P[abs(P) < 1e-12] <- 0
  P[tot == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

custom_106_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_106_P[i] <- routledge_beta_e(mite_beta)
}

sum(abs(vegan_106_P)-abs(custom_106_P),na.rm = T)

# 107 ####

vegan_107_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_107_P[i] <- vegan::betadiver(mite_beta,"g")
}
sum(abs(vegan_107_P)-abs(custom_105_P),na.rm = T)

# 108 ####
vegan_108_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_108_P[i] <- vegan::betadiver(mite_beta,"gl") 
}
lennon_beta_gl <- function(x, na.rm = TRUE) {
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
  
  P <- 2 * abs(b - c) / (2 * a + b + c)
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}
custom_108_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_108_P[i] <- lennon_beta_gl(mite_beta)
}

sum(abs(vegan_108_P)-abs(custom_108_P),na.rm = T)

# 109 ####

harte_kinzig_beta <- function(x, na.rm = TRUE) {
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
  
  P <- 1 - 2 * a / (2 * a + b + c)
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}
vegan_109_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_109_P[i] <- vegan::betadiver(mite_beta,"hk") #0.2
}

custom_109_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_109_P[i] <- harte_kinzig_beta(mite_beta)
}

sum(abs(vegan_109_P)-abs(custom_109_P),na.rm = T)
sum(abs(custom_102_P)-abs(custom_109_P),na.rm = T)

# 110 ####
vegan_110_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_110_P[i] <- vegan::betadiver(mite_beta,"l") 
}
sum(abs(vegan_110_P)-abs(custom_104_P),na.rm = T)

# 111_P ####
vegan_111_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_111_P[i] <- vegan::betadiver(mite_beta,"I") 
}

routledge_beta_i <- function(x, na.rm = TRUE) {
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
  
  # x * log(x), with 0 * log(0) = 0
  xlogx <- function(z) ifelse(z > 0, z * log(z), 0)
  
  tot <- 2 * a + b + c
  
  P <- log(tot) - (1 / tot) * 2 * a * log(2) -
    (1 / tot) * (xlogx(a + b) + xlogx(a + c))
  
  # Remove rounding errors around 0 (identical samples)
  P[abs(P) < 1e-12] <- 0
  P[tot == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

custom_111_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_111_P[i] <- routledge_beta_i(mite_beta)
}

sum(abs(vegan_111_P)-abs(custom_111_P),na.rm = T)

# 112_P ####
vegan_112_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_112_P[i] <- vegan::betadiver(mite_beta,"m")
}
magurran_beta_m <- function(x, na.rm = TRUE) {
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
  
  P <- (2 * a + b + c) * (1 - a / (a + b + c))
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

custom_112_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_112_P[i] <- magurran_beta_m(mite_beta)
}


sum(abs(vegan_112_P)-abs(custom_112_P),na.rm = T)

# 113_P ####
vegan_113_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_113_P[i] <- vegan::betadiver(mite_beta,"me")
}

sum(abs(vegan_113_P)-abs(custom_102_P),na.rm = T)
sum(abs(vegan_113_P)-abs(custom_118_P),na.rm = T)

# 114_P ####

vegan_114_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_114_P[i] <- vegan::betadiver(mite_beta,"r")
}
tabula_114_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  tabula_114_P[i] <- tabula::index_routledge1(as.matrix(mite_beta)) 
}

routledge_beta_r <- function(x, na.rm = TRUE) {
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
  
  gamma2 <- (a + b + c)^2
  
  P <- gamma2 / (gamma2 - 2 * b * c)
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

custom_114_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_114_P[i] <- routledge_beta_r(mite_beta)
}

sum(abs(vegan_114_P)-abs(custom_114_P),na.rm = T)
sum(abs(tabula_114_P)-abs(custom_114_P),na.rm = T)

# 115_P ####
vegan_115_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_115_P[i] <- vegan::betadiver(mite_beta,"rlb")
}

ruggiero_beta_rlb <- function(x, na.rm = TRUE) {
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

  
  P <- a / (a + c)
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}


custom_115_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_115_P[i] <- ruggiero_beta_rlb(mite_beta)
}

sum(abs(vegan_115_P)-abs(custom_115_P),na.rm = T)

# 116_P ####
vegan_116_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_116_P[i] <- vegan::betadiver(mite_beta,"sim") 
}

lennon_beta_sim <- function(x, na.rm = TRUE) {
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
  
  
  P <- pmin(b,c) / (pmin(b,c) + a)
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}


custom_116_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_116_P[i] <- lennon_beta_sim(mite_beta)
}

sum(abs(vegan_116_P)-abs(custom_116_P),na.rm = T)

# 117_P ####
vegan_117_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_117_P[i] <- vegan::betadiver(mite_beta,"t") 
}

wilson_beta_t <- function(x, na.rm = TRUE) {
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
  
  
  P <- (b + c) / ((2 * a) + b + c)
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}


custom_117_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_117_P[i] <- wilson_beta_t(mite_beta)
}

sum(abs(vegan_117_P)-abs(custom_117_P),na.rm = T)

# 119_P ####

vegan_119_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_119_P[i] <- vegan::betadiver(mite_beta,"wb")
}

weiher_beta_wb <- function(x, na.rm = TRUE) {
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
  
  
  P <- (b + c) 
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}


custom_119_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_119_P[i] <- weiher_beta_wb(mite_beta)
}

sum(abs(vegan_119_P)-abs(custom_119_P),na.rm = T)

# 120_P ####
vegan_120_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_120_P[i] <- vegan::betadiver(mite_beta,"z") 
}

lennon_beta_z <- function(x, na.rm = TRUE) {
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
  
  P <- 1 - log((2 * a + b + c) / (a + b + c)) / log(2)
  P[(a + b + c) == 0] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

custom_120_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_120_P[i] <- lennon_beta_z(mite_beta)
}

sum(abs(vegan_120_P)-abs(custom_120_P),na.rm = T)

# 121_P ####
vegan_121_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_121_P[i] <- vegan::betadiver(mite_beta,"19")
}

williams_beta19 <- function(x, na.rm = TRUE) {
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
  
  N <- a + b + c      # total number of species of the pair
  
  P <- 2 * (b * c + 1) / (N^2 - N)
  P[N <= 1] <- NA_real_
  
  dimnames(P) <- list(rownames(x), rownames(x))
  stats::as.dist(P)
}

custom_121_P <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_121_P[i] <- williams_beta19(mite_beta)
}

sum(abs(vegan_121_P)-abs(custom_121_P),na.rm = T)


# ABONDANCE BASED INDICES ####

# 2_A ####
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
vegan_2_A <- rep(NA, 69)
for (i in 1:(nrow(mite) - 1)) {
  mite_beta <- mite[c(i, i + 1), ]
  vegan_2_A[i] <- vegan::vegdist(mite_beta, "aitchison", pseudocount = 1)
}

aitchison_distance <- function(x, pseudocount = 1, na.rm = TRUE) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x)) || nrow(x) < 2)
    stop("`x` must be a data frame or matrix with at least 2 samples (rows).")
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  if (!is.numeric(pseudocount) || length(pseudocount) != 1 || pseudocount < 0)
    stop("`pseudocount` must be a single non-negative value.")
  
  # Missing values
  if (any(is.na(x))) {
    if (!na.rm) stop("`x` contains missing values.")
    x[is.na(x)] <- 0
  }
  
  # Add pseudocount to all abundances
  x <- x + pseudocount
  
  if (any(x == 0))
    warning("Some samples contain zero abundances: ln(0) is undefined, distances set to NA.")
  
  n_samples <- nrow(x)
  N <- ncol(x)
  
  A <- matrix(NA_real_, n_samples, n_samples)
  
  for (j in 2:n_samples) {
    for (k in 1:(j - 1)) {
      if (any(x[j, ] == 0) || any(x[k, ] == 0)) next
      
      s <- 0
      for (i in 1:N) {
        term_j <- log(x[j, i]) - (1 / N) * sum(log(x[j, ]))
        term_k <- log(x[k, i]) - (1 / N) * sum(log(x[k, ]))
        s <- s + (term_j - term_k)^2
      }
      
      A[j, k] <- sqrt(s)
    }
  }
  
  dimnames(A) <- list(rownames(x), rownames(x))
  stats::as.dist(A)
}

custom_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_2_A[i] <- aitchison_distance(mite_beta)
}

sum(abs(vegan_2_A)-abs(custom_2_A),na.rm = T)
sum(abs(pctax_2_A)-abs(custom_2_A),na.rm = T)
sum(abs(ecodive_2_A)-abs(custom_2_A),na.rm = T)
# robust 
vegan_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_2_A[i] <- vegan::vegdist(mite_beta,method = "robust.aitchison")
}
robust_aitchison_distance <- function(x, na.rm = TRUE) {
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
  
  n_samples <- nrow(x)
  N <- ncol(x)
  
  # Robust centred log-ratio of each sample
  r <- matrix(0, n_samples, N)
  for (j in 1:n_samples) {
    present <- x[j, ] > 0
    if (!any(present)) next
    N_j <- sum(present)
    mean_log_j <- (1 / N_j) * sum(log(x[j, present]))
    r[j, present] <- log(x[j, present]) - mean_log_j
  }
  
  empty <- rowSums(x > 0) == 0
  
  A <- matrix(NA_real_, n_samples, n_samples)
  
  for (j in 2:n_samples) {
    for (k in 1:(j - 1)) {
      if (empty[j] || empty[k]) next
      
      s <- 0
      for (i in 1:N) {
        s <- s + (r[j, i] - r[k, i])^2
      }
      
      A[j, k] <- sqrt(s)
    }
  }
  
  dimnames(A) <- list(rownames(x), rownames(x))
  stats::as.dist(A)
}

custom_2_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_2_A[i] <- robust_aitchison_distance(mite_beta)
}

sum(abs(vegan_2_A)-abs(custom_2_A),na.rm = T)
sum(abs(pctax_2_A)-abs(custom_2_A),na.rm = T)
sum(abs(ecodive_2_A)-abs(custom_2_A),na.rm = T)

# 3_A ####
EnvNJ_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_3_A[i] <- EnvNJ::metrics(t(mite_beta),method = "avg")
}

average_distance <- function(x, na.rm = TRUE) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x)) || nrow(x) < 2)
    stop("`x` must be a data frame or matrix with at least 2 samples (rows).")
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  
  # Missing values
  if (any(is.na(x))) {
    if (!na.rm) stop("`x` contains missing values.")
    x[is.na(x)] <- 0
  }
  
  n_samples <- nrow(x)
  A <- matrix(NA_real_, n_samples, n_samples)
  
  for (j in 2:n_samples) {
    for (k in 1:(j - 1)) {
      diff <- abs(x[j, ] - x[k, ])
      A[j, k] <- (1 / 2) * (sum(diff) + max(diff))
    }
  }
  
  dimnames(A) <- list(rownames(x), rownames(x))
  stats::as.dist(A)
}

custom_3_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_3_A[i] <- average_distance(mite_beta)
}

sum(abs(EnvNJ_3_A)-abs(custom_3_A),na.rm = T)

# 4_A  ####
ecodive_4_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_4_A[i] <- ecodive::bhattacharyya(mite_beta,rescale = F)
}
proxy_4_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_4_A[i] <- proxy::dist(mite_beta,method = "Bhjattacharyya")
}
philentropy_4_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_4_A[i] <- philentropy::bhattacharyya(as.numeric(mite[i,]),as.numeric(mite[i+1,])) #-4.223818
}
EnvNJ_4_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  EnvNJ_4_A[i] <- EnvNJ::metrics(t(mite_beta),method = "bhattacharyya")
}
Rfast_4_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Rfast_4_A[i] <- Rfast::Dist(mite_beta,method = "bhattacharyya")[1,2]
}

bhattacharyya_distance <- function(x, na.rm = TRUE) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x)) || nrow(x) < 2)
    stop("`x` must be a data frame or matrix with at least 2 samples (rows).")
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  
  # Missing values
  if (any(is.na(x))) {
    if (!na.rm) stop("`x` contains missing values.")
    x[is.na(x)] <- 0
  }
  
  n_samples <- nrow(x)
  A <- matrix(NA_real_, n_samples, n_samples)
  
  for (j in 2:n_samples) {
    for (k in 1:(j - 1)) {
      A[j, k] <- - log(sum(sqrt(x[j, ] * x[k, ])))
    }
  }
  
  dimnames(A) <- list(rownames(x), rownames(x))
  stats::as.dist(A)
}

custom_4_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_4_A[i] <- bhattacharyya_distance(mite_beta)
}

sum(abs(ecodive_4_A)-abs(custom_4_A),na.rm = T)
sum(abs(proxy_4_A)-abs(custom_4_A),na.rm = T)
sum(abs(philentropy_4_A)-abs(custom_4_A),na.rm = T)
sum(abs(EnvNJ_4_A)-abs(custom_4_A),na.rm = T)
sum(abs(Rfast_4_A)-abs(custom_4_A),na.rm = T)

# 6_A -> 6_P ####
tabula_6_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  tabula_6_A[i] <- tabula::index_binomial(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

binomial_cooccurrence <- function(x, na.rm = TRUE) {
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
  
  N <- ncol(pa)       # total number of species (a + b + c + d)
  
  p <- (a + b) * (a + c) / N^2
  denom <- sqrt(N * p * (1 - p))
  
  A <- (a - N * p) / denom
  A[denom == 0] <- NA_real_
  
  dimnames(A) <- list(rownames(x), rownames(x))
  stats::as.dist(A)
}

custom_6_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_6_A[i] <- binomial_cooccurrence(mite_beta)
}

sum(abs(tabula_6_A)-abs(custom_6_A),na.rm = T)

# 7_A ####
tabula_7_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  tabula_7_A[i] <- tabula::index_brainerd(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}
brsim_7_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  brsim_7_A[i] <- brsim::brsim(mite_beta)$BR.similarity.matrix[1,2]
}

brainerd_robinson_similarity <- function(x, na.rm = TRUE) {
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
  
  # Percentages (each sample sums to 100)
  empty <- rowSums(x) == 0
  x <- 100 * x / rowSums(x)
  
  n_samples <- nrow(x)
  A <- matrix(NA_real_, n_samples, n_samples)
  
  for (j in 2:n_samples) {
    for (k in 1:(j - 1)) {
      if (empty[j] || empty[k]) next
      diff <- abs(x[j, ] - x[k, ])
      A[j, k] <- 200 - sum(diff)
    }
  }
  
  dimnames(A) <- list(rownames(x), rownames(x))
  stats::as.dist(A)
}

custom_7_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_7_A[i] <- brainerd_robinson_similarity(mite_beta)
}

sum(abs(tabula_7_A)-abs(custom_7_A),na.rm = T)
sum(abs(brsim_7_A)-abs(custom_7_A),na.rm = T)

# 13_A ####
abdiv_13_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_13_A[i] <- abdiv::cy_dissimilarity(as.numeric(mite[i,]),as.numeric(mite[i+1,]),base = 10) # base can be modified
}
cao_original <- function(x, na.rm = TRUE) {
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
  
  n_samples <- nrow(x)
  A <- matrix(NA_real_, n_samples, n_samples)
  
  for (j in 2:n_samples) {
    for (k in 1:(j - 1)) {
      # Species present in at least one of the two samples
      keep <- x[j, ] > 0 | x[k, ] > 0
      N <- sum(keep)
      if (N == 0) next
      
      xj <- x[j, keep]
      xk <- x[k, keep]
      
      # Abundances <= 0.1 (including zeros) replaced by 0.1 (Cao et al. 1997)
      xj[xj <= 0.1] <- 0.1
      xk[xk <= 0.1] <- 0.1
      
      n_i <- xj + xk
      
      A[j, k] <- (1 / N) * sum(log10(n_i / 2) - (xj * log10(xk) + xk * log10(xj)) / n_i)
    }
  }
  
  dimnames(A) <- list(rownames(x), rownames(x))
  stats::as.dist(A)
}

custom_13_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_13_A[i] <- cao_original(mite_beta)
}

sum(abs(abdiv_13_A)-abs(custom_13_A),na.rm = T)

# 14_A ####
vegan_14_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_14_A[i] <- vegan::vegdist(mite_beta,method = "cao")
}
pctax_14_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_14_A[i] <- pctax::mat_dist(t((mite_beta)),method = "cao")
}
NST_14_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_14_A[i] <- NST::beta.g(mite_beta,dist.method = "cao")
}
cao_index <- function(x, na.rm = TRUE) {
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
  
  n_samples <- nrow(x)
  A <- matrix(NA_real_, n_samples, n_samples)
  
  for (j in 2:n_samples) {
    for (k in 1:(j - 1)) {
      # Species present in at least one of the two samples
      keep <- x[j, ] > 0 | x[k, ] > 0
      N <- sum(keep)
      if (N == 0) next
      
      xj <- x[j, keep]
      xk <- x[k, keep]
      
      # Zero abundances replaced by 0.1 (Cao et al. 1997)
      xj[xj == 0] <- 0.1
      xk[xk == 0] <- 0.1
      
      n_i <- xj + xk
      
      A[j, k] <- (1 / N) * sum(log(n_i / 2) - (xj * log(xk) + xk * log(xj)) / n_i)
    }
  }
  
  dimnames(A) <- list(rownames(x), rownames(x))
  stats::as.dist(A)
}

custom_14_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_14_A[i] <- cao_index(mite_beta)
}

sum(abs(vegan_14_A)-abs(custom_14_A),na.rm = T)
sum(abs(pctax_14_A)-abs(custom_14_A),na.rm = T)
sum(abs(NST_14_A)-abs(custom_14_A),na.rm = T)

# 16_A ####
vegan_16_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_16_A[i] <- vegan::chaodist(mite_beta,formula = "1 - sqrt(U*V)")
}

chao_ochiai_distance <- function(x, na.rm = TRUE) {
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
  
  n_samples <- nrow(x)
  A <- matrix(NA_real_, n_samples, n_samples)
  
  for (j in 2:n_samples) {
    for (k in 1:(j - 1)) {
      if (sum(x[j, ]) == 0 || sum(x[k, ]) == 0) next
      
      # Shared species
      shared <- x[j, ] > 0 & x[k, ] > 0
      
      U <- sum(x[j, shared]) / sum(x[j, ])
      V <- sum(x[k, shared]) / sum(x[k, ])
      
      A[j, k] <- 1 - sqrt(U * V)
    }
  }
  
  dimnames(A) <- list(rownames(x), rownames(x))
  stats::as.dist(A)
}

custom_16_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_16_A[i] <- chao_ochiai_distance(mite_beta)
}

sum(abs(vegan_16_A)-abs(custom_16_A),na.rm = T)

# 18_A ####
wiqid_18_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  wiqid_18_A[i] <- wiqid::distChaoSorNaive(mite[i,],mite[i+1,])
} 
NST_18_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_18_A[i] <- NST::beta.g(mite_beta,dist.method = "chao.sorensen")
} 
CommEcol_18_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  CommEcol_18_A[i] <- CommEcol::dis.chao(mite_beta,index = "sorensen",version = "probability") #0.02642424
} 


chao_sorensen <- function(x, na.rm = TRUE) {
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
  
  n_samples <- nrow(x)
  A <- matrix(NA_real_, n_samples, n_samples)
  
  for (j in 2:n_samples) {
    for (k in 1:(j - 1)) {
      if (sum(x[j, ]) == 0 || sum(x[k, ]) == 0) next
      
      # Shared species
      shared <- x[j, ] > 0 & x[k, ] > 0
      
      U <- sum(x[j, shared]) / sum(x[j, ])
      V <- sum(x[k, shared]) / sum(x[k, ])
      
      A[j, k] <- 1 - ((2 * U * V)/(U + V))
    }
  }
  
  dimnames(A) <- list(rownames(x), rownames(x))
  stats::as.dist(A)
}

custom_18_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_18_A[i] <- chao_sorensen(mite_beta)
}

sum(abs(wiqid_18_A)-abs(custom_18_A),na.rm = T)
sum(abs(NST_18_A)-abs(custom_18_A),na.rm = T)
sum(abs(CommEcol_18_A)-abs(custom_18_A),na.rm = T)

# 20_A ####
vegan_20_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  vegan_20_A[i] <- vegan::vegdist(mite_beta,method = "chisq")
}
pctax_20_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  pctax_20_A[i] <- pctax::mat_dist(t((mite_beta)),method = "chisq")
}
svs_20_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  svs_20_A[i] <- svs::dist_chisquare(as.matrix(mite_beta))
}
analogue_20_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  analogue_20_A[i] <- analogue::distance(mite[i,],mite[i+1,],method = "chi.square")
}
adespatial_20_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_20_A[i] <- adespatial::beta.div(mite_beta,method = "chisquare",save.D = T)$D
}
proxy_20_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  proxy_20_A[i] <- proxy::dist(mite_beta,method = "Chi-squared")
}
spaa_20_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  spaa_20_A[i] <- spaa::sp.pair(t(as.matrix(mite_beta)))$chisq
}

chi_square_distance <- function(x, na.rm = TRUE) {
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
  
  # Species absent from all samples are excluded
  x <- x[, colSums(x) > 0, drop = FALSE]
  
  y_plus_plus <- sum(x)       # total of the table
  y_plus_j <- colSums(x)      # total of each species over all samples
  y_i_plus <- rowSums(x)      # total of each sample
  p <- ncol(x)
  
  n_samples <- nrow(x)
  A <- matrix(NA_real_, n_samples, n_samples)
  
  for (j1 in 2:n_samples) {
    for (j2 in 1:(j1 - 1)) {
      if (y_i_plus[j1] == 0 || y_i_plus[j2] == 0) next
      
      s <- 0
      for (j in 1:p) {
        s <- s + (1 / y_plus_j[j]) *
          (x[j1, j] / y_i_plus[j1] - x[j2, j] / y_i_plus[j2])^2
      }
      
      A[j1, j2] <- sqrt(y_plus_plus) * sqrt(s)
    }
  }
  
  dimnames(A) <- list(rownames(x), rownames(x))
  stats::as.dist(A)
}

custom_20_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_20_A[i] <- chi_square_distance(mite_beta)
}

sum(abs(vegan_20_A)-abs(custom_20_A),na.rm = T)
sum(abs(pctax_20_A)-abs(custom_20_A),na.rm = T)
sum(abs(svs_20_A)-abs(custom_20_A),na.rm = T)
sum(abs(analogue_20_A)-abs(custom_20_A),na.rm = T)
sum(abs(adespatial_20_A)-abs(custom_20_A),na.rm = T)
sum(abs(proxy_20_A)-abs(custom_20_A),na.rm = T)
sum(abs(spaa_20_A)-abs(custom_20_A),na.rm = T)

# 22_A ####
ecodive_22_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_22_A[i] <- ecodive::psym_chisq(mite_beta,rescale = F)
}

prob_symmetric_chisq <- function(x, na.rm = TRUE) {
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
  
  n_samples <- nrow(x)
  A <- matrix(NA_real_, n_samples, n_samples)
  
  for (j in 2:n_samples) {
    for (k in 1:(j - 1)) {
      # Species absent from both samples are excluded (0/0)
      keep <- (x[j, ] + x[k, ]) > 0
      
      A[j, k] <- 2 * sum((x[j, keep] - x[k, keep])^2 / (x[j, keep] + x[k, keep]))
    }
  }
  
  dimnames(A) <- list(rownames(x), rownames(x))
  stats::as.dist(A)
}

custom_22_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_22_A[i] <- prob_symmetric_chisq(mite_beta)
}

sum(abs(ecodive_22_A)-abs(custom_22_A),na.rm = T)

# 25_A ####
adespatial_25_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_25_A[i] <- adespatial::beta.div(mite_beta,method = "log.chord",save.D = T)$D
}

log_chord_distance <- function(x, na.rm = TRUE) {
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
  
  # log(y + 1) transformation
  ly <- log(x + 1)
  
  n_samples <- nrow(ly)
  N <- ncol(ly)
  A <- matrix(NA_real_, n_samples, n_samples)
  
  for (j1 in 2:n_samples) {
    for (j2 in 1:(j1 - 1)) {
      norm_1 <- sqrt(sum(ly[j1, ]^2))
      norm_2 <- sqrt(sum(ly[j2, ]^2))
      if (norm_1 == 0 || norm_2 == 0) next
      
      s <- 0
      for (j in 1:N) {
        s <- s + (ly[j1, j] / norm_1 - ly[j2, j] / norm_2)^2
      }
      
      A[j1, j2] <- sqrt(s)
    }
  }
  
  dimnames(A) <- list(rownames(x), rownames(x))
  stats::as.dist(A)
}

custom_25_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_25_A[i] <- log_chord_distance(mite_beta)
}

sum(abs(adespatial_25_A)-abs(custom_25_A),na.rm = T)

# 27_A ####
ecodive_27_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodive_27_A[i] <- ecodive::clark(mite_beta,rescale = F)
}
philentropy_27_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  philentropy_27_A[i] <- philentropy::clark_sq(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

clark_distance <- function(x, na.rm = TRUE) {
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
  
  n_samples <- nrow(x)
  A <- matrix(NA_real_, n_samples, n_samples)
  
  for (j1 in 2:n_samples) {
    for (j2 in 1:(j1 - 1)) {
      # Species present in at least one of the two samples
      keep <- (x[j1, ] + x[j2, ]) > 0
      N <- sum(keep)
      if (N == 0) next
      
      y1 <- x[j1, keep]
      y2 <- x[j2, keep]
      
      A[j1, j2] <- sqrt(sum(((y1 - y2) / (y1 + y2))^2))
    }
  }
  
  dimnames(A) <- list(rownames(x), rownames(x))
  stats::as.dist(A)
}

custom_27_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_27_A[i] <- clark_distance(mite_beta)
}

sum(abs(ecodive_27_A)-abs(custom_27_A),na.rm = T)
sum(abs(philentropy_27_A)-abs(custom_27_A),na.rm = T)

# 29_A ####
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
divergence_distance <- function(x, na.rm = TRUE) {
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
  
  n_samples <- nrow(x)
  A <- matrix(NA_real_, n_samples, n_samples)
  
  for (j in 2:n_samples) {
    for (k in 1:(j - 1)) {
      # Species absent from both samples are excluded (0/0)
      keep <- (x[j, ] + x[k, ]) > 0
      if (!any(keep)) next
      
      xj <- x[j, keep]
      xk <- x[k, keep]
      
      A[j, k] <- 2 * sum((xj - xk)^2 / (xj + xk)^2)
    }
  }
  
  dimnames(A) <- list(rownames(x), rownames(x))
  stats::as.dist(A)
}

custom_29_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_29_A[i] <- divergence_distance(mite_beta)
}

sum(abs(ecodive_29_A)-abs(custom_29_A),na.rm = T)
sum(abs(philentropy_29_A)-abs(custom_29_A),na.rm = T)
sum(abs(proxy_29_A)-abs(custom_27_A),na.rm = T)

# 31_A ####
NST_31_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_31_A[i] <- NST::beta.g(mite_beta,dist.method = "mEuclidean")
}

modified_euclidean <- function(x, na.rm = TRUE) {
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
  
  n_samples <- nrow(x)
  A <- matrix(NA_real_, n_samples, n_samples)
  
  for (j in 2:n_samples) {
    for (k in 1:(j - 1)) {
      # Weights: 0 for joint absences, 1 otherwise
      w <- ifelse(x[j, ] + x[k, ] == 0, 0, 1)
      if (sum(w) == 0) next
      
      A[j, k] <- sqrt(sum((x[j, ] - x[k, ])^2)) / sum(w)
    }
  }
  
  dimnames(A) <- list(rownames(x), rownames(x))
  stats::as.dist(A)
}

custom_31_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_31_A[i] <- modified_euclidean(mite_beta)
}

sum(abs(NST_31_A)-abs(custom_31_A),na.rm = T)

# 32_A ####
abdiv_32_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  abdiv_32_A[i] <- abdiv::rms_distance(as.numeric(mite[i,]),as.numeric(mite[i+1,]))
}

average_euclidean <- function(x, na.rm = TRUE) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x)) || nrow(x) < 2)
    stop("`x` must be a data frame or matrix with at least 2 samples (rows).")
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  
  # Missing values
  if (any(is.na(x))) {
    if (!na.rm) stop("`x` contains missing values.")
    x[is.na(x)] <- 0
  }
  
  n_samples <- nrow(x)
  N <- ncol(x)
  A <- matrix(NA_real_, n_samples, n_samples)
  
  for (j in 2:n_samples) {
    for (k in 1:(j - 1)) {
      A[j, k] <- sqrt(sum((x[j, ] - x[k, ])^2) / N)
    }
  }
  
  dimnames(A) <- list(rownames(x), rownames(x))
  stats::as.dist(A)
}

custom_32_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_32_A[i] <- average_euclidean(mite_beta)
}

sum(abs(abdiv_32_A)-abs(custom_32_A),na.rm = T)

# 35_A ####
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

alt_gower <- function(x, na.rm = TRUE) {
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
  
  n_samples <- nrow(x)
  A <- matrix(NA_real_, n_samples, n_samples)
  
  for (j in 2:n_samples) {
    for (k in 1:(j - 1)) {
      # Weights: 0 for joint absences, 1 otherwise
      omega <- ifelse(x[j, ] + x[k, ] == 0, 0, 1)
      if (sum(omega) == 0) next
      
      A[j, k] <- sum(omega * abs(x[j, ] - x[k, ])) / sum(omega)
    }
  }
  
  dimnames(A) <- list(rownames(x), rownames(x))
  stats::as.dist(A)
}


custom_35_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_35_A[i] <- alt_gower(mite_beta)
}

sum(abs(NST_35_A)-abs(custom_35_A),na.rm = T)
sum(abs(vegan_35_A)-abs(custom_35_A),na.rm = T)
sum(abs(pctax_35_A)-abs(custom_35_A),na.rm = T)
sum(abs(philentropy_35_A)-abs(custom_35_A),na.rm = T)
sum(abs(Rfast_35_A)-abs(custom_35_A),na.rm = T)
sum(abs(gower_35_A)-abs(custom_35_A),na.rm = T)

# 36_A ####
ecodist_36_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  ecodist_36_A[i] <- ecodist::distance(mite_beta,method = "modgower10")
}

modgower10_ecodist <- function(x, na.rm = TRUE) {
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
  
  # log10 transformation: log10(x) if x > 0, 0 otherwise
  y <- ifelse(x > 0, log10(x), 0)
  
  n_samples <- nrow(y)
  A <- matrix(NA_real_, n_samples, n_samples)
  
  for (j in 2:n_samples) {
    for (k in 1:(j - 1)) {
      A[j, k] <- sum(abs(y[j, ] - y[k, ]))
    }
  }
  
  dimnames(A) <- list(rownames(x), rownames(x))
  stats::as.dist(A)
}

custom_36_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_36_A[i] <- modgower10_ecodist(mite_beta)
}

sum(abs(ecodist_36_A)-abs(custom_36_A),na.rm = T)

# 37_A ####
NST_37_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  NST_37_A[i] <- NST::beta.g(mite_beta,dist.method = "mGower")
}
mgower_anderson <- function(x, na.rm = TRUE) {
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
  
  # Anderson et al. (2006) transformation: log10(x) + 1 if x > 0, 0 otherwise
  y <- ifelse(x > 0, log10(x) + 1, 0)
  
  n_samples <- nrow(y)
  A <- matrix(NA_real_, n_samples, n_samples)
  
  for (j in 2:n_samples) {
    for (k in 1:(j - 1)) {
      # Weights: 0 for joint absences (on transformed values), 1 otherwise
      omega <- ifelse(y[j, ] == 0 & y[k, ] == 0, 0, 1)
      if (sum(omega) == 0) next
      
      A[j, k] <- sum(omega * abs(y[j, ] - y[k, ])) / sum(omega)
    }
  }
  
  dimnames(A) <- list(rownames(x), rownames(x))
  stats::as.dist(A)
}

custom_37_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  custom_37_A[i] <- mgower_anderson(mite_beta)
}

sum(abs(NST_37_A)-abs(custom_37_A),na.rm = T)

# 45_A ####
adespatial_45_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  adespatial_45_A[i] <- adespatial::beta.div.comp(mite_beta,coef = "BJ",quant = T)$repl 
}