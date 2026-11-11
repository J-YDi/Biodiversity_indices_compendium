#_______________________________________________________________________________
# Title              : 09_Manual_verif_Alpha.r
# Date               : 11/11/2026
# Object             : Script to manually check values from packages that 
#                      calculate alpha diversity indices
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

# Functions to run custom functions #
richness <- function(x, binary = FALSE) {
  y <- rowSums(x > 0)
  return(y)
}


#___________________________ Alpha diversity indices ___________________________####

#DQ ####

abdiv_DQ <- rep(1,70)
for (i in 1:nrow(mite)){
  abdiv_DQ[i] <- abdiv::kempton_taylor_q(as.numeric(mite[i,]),lower_quantile = 0.5,upper_quantile = 0.9)
}

kempton_q <- function(x, q1 = 0.25, q2 = 0.75) {

  # Keep species present only
  x <- x[!is.na(x)]
  S <- sum(x > 0)
  x <- x[x > 0]
  
  # Abundance classes r and number of species per class nr
  r   <- sort(unique(x))
  nr  <- tabulate(match(x, r), nbins = length(r))
  cum <- cumsum(nr)
  
  # Quantile classes on the cumulative species curve
  i1 <- which(cum >= q1 * S)[1]
  i2 <- which(cum >= q2 * S)[1]
  R1 <- r[i1]
  R2 <- r[i2]
  
  if (R1 == R2) {
    warning("R1 = R2: Q is undefined (ln(R2/R1) = 0).")
    return(NA)
  }
  
  # Sum from R1 to R2 (inclusive, as in the formula)
  inner <- sum(nr[i1:i2])
  
  (0.5 * nr[i1] + inner + 0.5 * nr[i2]) / log(R2 / R1)
}

custom_DQ <- rep(1,70)
for (i in 1:nrow(mite)){
  custom_DQ[i] <- kempton_q(as.numeric(mite[i,]),q1 = 0.5,q2 = 0.9)
}

sum(abs(abdiv_DQ)-abs(custom_DQ),na.rm = T) 

# DBRI ####
tabula_DBRI <- as.numeric(tabula::heterogeneity(mite,method = "brillouin"))
abdiv_DBRI <- apply(mite, 1, abdiv::brillouin_d)
ecodive_DBRI <-ecodive::alpha_div(mite,metric = "brillouin")
wiqid_DBRI <- apply(mite, 1, wiqid::biodBrillouin) #FALSE


brillouin <- function(x, binary = FALSE) {
  y <- (lfactorial(rowSums(x)) - rowSums(lfactorial(x))) / rowSums(x)
  return(y)
}
custom_DBRI <- brillouin(mite)

sum(abs(wiqid_DBRI)-abs(custom_DBRI))
sum(abs(tabula_DBRI)-abs(custom_DBRI))
sum(abs(abdiv_DBRI)-abs(custom_DBRI))
sum(abs(ecodive_DBRI)-abs(custom_DBRI))

# DRAO ####
BAT_DRAO <- BAT::rao(as.matrix(mite))[,1] 
ade4_DRAO <- rep(1,70)
for (i in 1:nrow(mite)){
  ade4_DRAO[i] <- ade4::apqe(as.data.frame(t(mite[i,])),dis = NULL)$results[2,1]
}
SYNCSA_DRAO <- SYNCSA::rao.diversity(mite)$Simpson # package down

rao_q <- function(x, d) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
  }
  
  # Check distance: Euclidean `dist` object
  if (!inherits(d, "dist")) stop("`d` must be a `dist` object.")
  if (!ade4::is.euclid(d)) stop("`d` must be a Euclidean distance.")
  d <- as.matrix(d)^2
  if (nrow(d) != ncol(x)) stop("`d` must have as many species as `x` has columns.")
  
  # Align d on species names when available
  sp <- colnames(x)
  if (!is.null(sp) && !is.null(rownames(d))) {
    if (!all(sp %in% rownames(d)))
      stop("Some species in `x` are missing from `d`.")
    d <- d[sp, sp]
  }
  
  # NA treated as absence
  x[is.na(x)] <- 0
  
  Q <- rep(NA, nrow(x))
  
  for (k in seq_len(nrow(x))) {
    # Species present in sample k
    present <- which(x[k, ] > 0)
    S <- length(present)
    if (S == 0) next
    
    p   <- x[k, present] / sum(x[k, present])
    d_k <- d[present, present, drop = FALSE]
    
    res <- 0
    for (i in seq_len(S - 1)) {
      for (j in seq(from = i + 1, length.out = max(0, S - 1 - i))) {
        res <- res + d_k[i, j] * p[i] * p[j]
      }
    }
    Q[k] <- res
  }
  
  if (any(is.na(Q))) warning("Some samples are empty: Q set to NA.")
  
  names(Q) <- rownames(x)
  Q
}
custom_DRAO <- rao_q(mite,dist(t(mite),method = "euclidean",diag = T,upper = T))

sum(abs(BAT_DRAO)-abs(custom_DRAO))
sum(abs(ade4_DRAO)-abs(custom_DRAO))
sum(abs(SYNCSA_DRAO)-abs(custom_DRAO))


# EBRI ####
tabula_EBRI <- as.numeric(tabula::evenness(mite,method = "brillouin"))

brillouin_evenness <- function(x) {

  # NA treated as absence
  x[is.na(x)] <- 0
  
  E <- rep(NA, nrow(x))
  
  for (k in seq_len(nrow(x))) {
    n <- x[k, x[k, ] > 0]
    S <- length(n)
    if (S < 2) next
    
    N <- sum(n)
    L <- floor(N / S)
    r <- N - S * L
    
    # Brillouin index
    HB <- (lfactorial(N) - sum(lfactorial(n))) / N
    
    # Maximum Brillouin index
    HBmax <- (lfactorial(N) - (S - r) * lfactorial(L) - r * lfactorial(L + 1)) / N
    
    E[k] <- HB / HBmax
  }
  names(E) <- rownames(x)
  E
}

custom_EBRI <- brillouin_evenness(mite)
sum(abs(tabula_EBRI)-abs(custom_EBRI))

# EBU ####
BAT_EBU <- as.numeric(BAT::evenness(mite,func = "bulla"))
microbiome_EBU <- microbiome::evenness(t(mite), index = "bulla")[,1] # OK

bulla_evenness <- function(x) {
  
  # NA treated as absence
  x[is.na(x)] <- 0
  
  E <- rep(NA, nrow(x))
  
  for (k in seq_len(nrow(x))) {
    n <- x[k, x[k, ] > 0]
    S <- length(n)
    if (S < 2) next
    
    p <- n / sum(n)
    
    E[k] <- (sum(pmin(p, 1 / S)) - 1 / S) / (1 - 1 / S)
  }
  
  names(E) <- rownames(x)
  E
}
custom_EBU <- bulla_evenness(mite)

sum(abs(microbiome_EBU)-abs(custom_EBU))
sum(abs(BAT_EBU)-abs(custom_EBU))

#ECAM ####
BAT_ECAM <- as.numeric(BAT::evenness(mite,func = "camargo"))
microbiome_ECAM <- microbiome::evenness(t(mite), index = "camargo")[,1] # OK

camargo_evenness <- function(x) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
  }
  
  # NA treated as absence
  x[is.na(x)] <- 0
  
  E <- rep(NA, nrow(x))
  
  for (k in seq_len(nrow(x))) {
    n <- x[k, x[k, ] > 0]
    S <- length(n)
    if (S == 0) next
    
    p <- n / sum(n)
    
    # Pairwise differences, pairs i < j only
    diffs <- abs(outer(p, p, "-"))
    E[k] <- 1 - sum(diffs[upper.tri(diffs)] / S)
  }
  
  names(E) <- rownames(x)
  E
}

custom_ECAM <- camargo_evenness(mite)
sum(abs(microbiome_ECAM)-abs(custom_ECAM))
sum(abs(BAT_ECAM)-abs(custom_ECAM))

#EEQ ####
eq_evenness <- function(x, na.rm = TRUE) {
  
  E <- rep(NA, nrow(x))
  
  for (k in seq_len(nrow(x))) {
    n <- x[k, ]
    
    # Missing values
    if (any(is.na(n))) {
      if (!na.rm) next
      n <- n[!is.na(n)]
    }
    
    # Species present only
    n <- n[n > 0]
    S <- length(n)
    if (S < 2) next
    
    # Perfectly even sample: infinite slope, E_Q = 1
    if (max(n) == min(n)) {
      E[k] <- 1
      next
    }
    
    # Decreasing ranks (rank 1 = most abundant), ties averaged
    r <- rank(-n, ties.method = "average")
    scaled_rank <- r / max(r)
    log_ab <- log(n)
    
    # Slope z of scaled rank ~ log(abundance)
    z <- cov(log_ab, scaled_rank) / var(log_ab)
    
    E[k] <- -2 / pi * atan(z)
  }
  
  names(E) <- rownames(x)
  E
}
custom_EEQ <- rep(1,70) 
for (i in 1:nrow(mite)){
  custom_EEQ[i] <- eq_evenness(mite[i,]) 
}
codyn_EEQ <- rep(1,70) 
for (i in 1:nrow(mite)){
  mite_long <- tidyr::pivot_longer(mite[i,],cols = colnames(mite)[1]:colnames(mite)[35],names_to = "Esp",values_to = "Value")
  mite_long <- mite_long[1:35,]
  codyn_EEQ[i] <- codyn::community_structure(mite_long,abundance.var = "Value",metric = "EQ")$EQ
}
sum(abs(codyn_EEQ)-abs(custom_EEQ))


# EEVAR  ####

evar_index <- function(x, na.rm = TRUE) {

  E <- rep(NA, nrow(x))
  
  for (k in seq_len(nrow(x))) {
    n <- x[k, ]
    
    # Missing values
    if (any(is.na(n))) {
      if (!na.rm) next
      n <- n[!is.na(n)]
    }
    
    # Species present only
    n <- n[n > 0]
    S <- length(n)
    if (S < 2) next
    
    # Population variance of log abundances divided by S
    log_ab <- log(n)
    var_log <- sum((log_ab - mean(log_ab))^2) / S
    
    E[k] <- 1 - (2 / pi) * atan(var_log)
  }
  
  names(E) <- rownames(x)
  E
}
custom_EEVAR <- rep(1,70) 
for (i in 1:nrow(mite)){
  custom_EEVAR[i] <- evar_index(mite[i,]) # OK S&W Eveness 
}
codyn_EEVAR <- rep(1,70) 
for (i in 1:nrow(mite)){
  mite_long <- tidyr::pivot_longer(mite[i,],cols = colnames(mite)[1]:colnames(mite)[35],names_to = "Esp",values_to = "Value")
  mite_long <- mite_long[1:35,]
  codyn_EEVAR[i] <- codyn::community_structure(mite_long,abundance.var = "Value",metric = "Evar")$Evar 
  
}

adiv_EEVAR <- adiv::specieseve(mite,method = "SmithWilson")[,1]
microbiome_EEVAR <- microbiome::evenness(t(mite), index = "evar")[,1]

sum(abs(microbiome_EEVAR)-abs(custom_EEVAR))
sum(abs(codyn_EEVAR)-abs(custom_EEVAR))
sum(abs(adiv_EEVAR)-abs(custom_EEVAR))


# EHEI ####
adiv_EHEI <- adiv::specieseve(mite,method = "Heip")[,1] # OK
abdiv_EHEI <- rep(1,70) 
for (i in 1:nrow(mite)){
  abdiv_EHEI[i] <- abdiv::heip_e(as.numeric(mite[i,])) # OK
}
heip_evenness <- function(x, na.rm = TRUE) {
  
  E <- rep(NA, nrow(x))
  
  for (k in seq_len(nrow(x))) {
    n <- x[k, ]
    
    # Missing values
    if (any(is.na(n))) {
      if (!na.rm) next
      n <- n[!is.na(n)]
    }
    
    # Species present only
    n <- n[n > 0]
    S <- length(n)
    if (S < 2) next
    
    # Shannon entropy
    p <- n / sum(n)
    H <- -sum(p * log(p))
    
    E[k] <- (exp(H) - 1) / (S - 1)
  }
  
  names(E) <- rownames(x)
  E
}

custom_EHEI <- rep(1,70) 
for (i in 1:nrow(mite)){
  custom_EHEI[i] <- heip_evenness(mite[i,]) 
}

sum(abs(abdiv_EHEI)-abs(custom_EHEI))
sum(abs(adiv_EHEI)-abs(custom_EHEI))

# EMC ####
tabula_EMC <- as.numeric(tabula::evenness(mite,method = "mcintosh"))
abdiv_EMC <- apply(mite,1,abdiv::mcintosh_e)

mcintosh_evenness <- function(x, na.rm = TRUE) {
  E <- rep(NA, nrow(x))
  
  for (k in seq_len(nrow(x))) {
    n <- x[k, ]
    
    # Missing values
    if (any(is.na(n))) {
      if (!na.rm) next
      n <- n[!is.na(n)]
    }
    
    # Species present only
    n <- n[n > 0]
    S <- length(n)
    if (S < 2) next
    
    N <- sum(n)
    
    E[k] <- (N - sqrt(sum(n^2))) / (N - N / sqrt(S))
  }
  
  names(E) <- rownames(x)
  E
}
custom_EMC <- mcintosh_evenness(mite)

sum(abs(tabula_EMC)-abs(custom_EMC))
sum(abs(abdiv_EMC)-abs(custom_EMC))

# EMN ####

microbiome_EMN <- rep(1,70) 
for (i in 1:nrow(mite)){
  microbiome_EMN[i] <- microbiome::dominance(as.numeric(mite[i,]),index = "DMN")$dmn  #OK
}

mcnaughton_dominance <- function(x, na.rm = TRUE) {
  E <- rep(NA, nrow(x))
  
  for (k in seq_len(nrow(x))) {
    n <- x[k, ]
    
    # Missing values
    if (any(is.na(n))) {
      if (!na.rm) next
      n <- n[!is.na(n)]
    }
    
    # Species present only
    n <- n[n > 0]
    S <- length(n)
    if (S == 0) next
    
    N <- sum(n)
    
    # Two most abundant species (second = 0 if only one species)
    n_sorted <- sort(n, decreasing = TRUE)
    Nmax1 <- n_sorted[1]
    Nmax2 <- if (S >= 2) n_sorted[2] else 0
    
    E[k] <- (Nmax1 + Nmax2) / N
  }
  
  
  names(E) <- rownames(x)
  E
}
custom_EMN <- mcnaughton_dominance(mite)

sum(abs(microbiome_EMN)-abs(custom_EMN))

# ESP ####

simpson_evenness <- function(x, na.rm = TRUE) {

  E <- rep(NA, nrow(x))
  
  for (k in seq_len(nrow(x))) {
    n <- x[k, ]
    
    # Missing values
    if (any(is.na(n))) {
      if (!na.rm) next
      n <- n[!is.na(n)]
    }
    
    # Species present only
    n <- n[n > 0]
    S <- length(n)
    if (S == 0) next
    
    # Simpson's index
    p <- n / sum(n)
    D <- sum(p^2)
    
    E[k] <- (1 / D) / S
  }
  
  names(E) <- rownames(x)
  E
}

custom_ESP <- simpson_evenness(mite)
tabula_ESP <- as.numeric(tabula::evenness(mite,method = "simpson")) # Inv Simpson E
sprex_ESP <- rep(1,70) 
for (i in 1:nrow(mite)){
  sprex_ESP[i] <- sprex::diversity(as.numeric(mite[i,]),type = "eveness.simpson") #FALSE
}
adiv_ESP <- adiv::specieseve(mite,method = "Simpson")[,1] # Inv Simpson E
codyn_ESP <- rep(1,70) 
for (i in 1:nrow(mite)){
  mite_long <- tidyr::pivot_longer(mite[i,],cols = colnames(mite)[1]:colnames(mite)[35],names_to = "Esp",values_to = "Value")
  mite_long <- mite_long[1:35,]
  codyn_ESP[i] <- codyn::community_structure(mite_long,abundance.var = "Value",metric = "SimpsonEvenness")$SimpsonEvenness
}

sum(abs(codyn_ESP)-abs(custom_ESP))
sum(abs(tabula_ESP)-abs(custom_ESP))
sum(abs(sprex_ESP)-abs(custom_ESP))
sum(abs(adiv_ESP)-abs(custom_ESP))

# EST ####
strong_dominance <- function(x, na.rm = TRUE) {

  E <- rep(NA, nrow(x))
  
  for (k in seq_len(nrow(x))) {
    n <- x[k, ]
    
    # Missing values
    if (any(is.na(n))) {
      if (!na.rm) next
      n <- n[!is.na(n)]
    }
    
    # Species present only
    n <- n[n > 0]
    S <- length(n)
    if (S < 2) next
    
    N <- sum(n)
    
    # Cumulative abundances, species ranked from highest to lowest
    t_i <- cumsum(sort(n, decreasing = TRUE))
    i <- seq_len(S)
    
    E[k] <- max(t_i / N - i / S)
  }
  
  names(E) <- rownames(x)
  E
}

custom_EST <- strong_dominance(mite)
abdiv_EST <- rep(1,70) 
for (i in 1:nrow(mite)){
  abdiv_EST [i] <- abdiv::strong(as.numeric(mite[i,]))
}

sum(abs(abdiv_EST)-abs(custom_EST))

# TRY ####
rygg <- function(x, binary = FALSE) {
  y <- log(richness(x)) / log(log(rowSums(x)))
  return(y)
}

custom_TRY <- rygg(mite)
benthos_TRY <- rep(1,70)
for (i in 1:nrow(mite)){
  benthos_TRY[i] <- benthos::rygg(taxon = names(mite[i,]), count = mite[i,],adjusted = F)
}

sum(abs(benthos_TRY)-abs(custom_TRY))

# TRYA ####
benthos_TRYA <- rep(1,70)
for (i in 1:nrow(mite)){
  benthos_TRYA[i] <- benthos::rygg(taxon = names(mite[i,]), count = mite[i,],adjusted = T)
}

rygg_adjusted <- function(x, binary = FALSE) {
  y <- log(richness(x)) / log(log(rowSums(x)+1)+1)
  return(y)
}
custom_TRYA <- rygg_adjusted(mite)

sum(abs(benthos_TRYA)-abs(custom_TRYA))

# RCA2 ####
chao2 <- function(x, na.rm = TRUE) {
  # Missing values
  if (any(is.na(x))) {
    if (!na.rm) return(NA)
    x[is.na(x)] <- 0
  }
  
  # Number of samples in which each species occurs
  occ <- colSums(x > 0)
  
  S_obs <- sum(occ > 0)
  Q1 <- sum(occ == 1)
  Q2 <- sum(occ == 2)
  
  if (Q2 == 0) {
    message("Q2 = 0: bias-corrected Chao 2 used.")
    return(S_obs + Q1 * (Q1 - 1) / (2 * (Q2 + 1)))
  }
  
  S_obs + Q1^2 / (2 * Q2)
}

# We cannot calculate this index sample by sample
# fossil_RCA2 <- rep(1,70)
# for (i in 1:nrow(mite)){
#   fossil_RCA2[i] <- fossil::chao2(as.numeric(mite[i,]))
# }
# wiqid_RCA2 <- rep(1,70)
# for (i in 1:nrow(mite)){
#   wiqid_RCA2[i] <- wiqid::richChao2(t(mite[i,]))[1]
# }
# tabula_RCA2 <- rep(1,70)
# for (i in 1:nrow(mite)){
#   tabula_RCA2[i] <- tabula::index_chao2(as.matrix(mite))
# }
# divDyn_RCA2 <- rep(1,70)
# for (i in 1:nrow(mite)){
#   divDyn_RCA2[i] <- divDyn::indices(as.matrix(mite),method = "chao2")
# }
# BAT_RCA2 <- BAT::alpha.accum(mite)[,18]
# 
# custom_RCA2 <- rep(1,70)
# for (i in 1:nrow(mite)){
#   custom_RCA2[i]  <- chao2(mite)
# }

fossil::chao2(mite) #70
wiqid::richChao2(t(mite))[1] #35
tabula::index_chao2(as.matrix(mite)) #35
divDyn::indices(as.matrix(mite),method = "chao2") # Don't respect the purpose
BAT::alpha.accum(mite) #35
chao2(mite) #35
# but according to formula tabula and wiqid compute Chao2 modified
# BAT compute another form of correction factor, either Chao2 and Chao2 modified

# RCA2M ####
chao2_corrected <- function(x, na.rm = TRUE) {

  # Missing values
  if (any(is.na(x))) {
    if (!na.rm) return(NA)
    x[is.na(x)] <- 0
  }
  
  # Number of samples in which each species occurs
  occ <- colSums(x > 0)
  
  S_obs <- sum(occ > 0)
  Q1 <- sum(occ == 1)
  Q2 <- sum(occ == 2)
  t  <- nrow(x)
  
  if (Q2 == 0) {
    message("Q2 = 0: bias-corrected Chao 2 used.")
    return(S_obs + ((t - 1) / t) * Q1 * (Q1 - 1) / (2 * (Q2 + 1)))
  }
  
  S_obs + ((t - 1) / t) * Q1^2 / (2 * Q2)
}

# ICE ####
fossil_RICE <- rep(1,70)
for (i in 1:nrow(mite)){
  fossil_RICE[i] <- fossil::ICE(t(mite[i,])) #44
}
wiqid_RICE <- rep(1,70)
for (i in 1:nrow(mite)){
  wiqid_RICE[i] <- wiqid::richICE(t(mite[i,])) #44
}
tabula_RICE <- rep(1,70)
for (i in 1:nrow(mite)){
  tabula_RICE[i] <- tabula::index_ice(t(mite[i,])) #44
}
BAT_RICE <- BAT::alpha.accum(mite)[,22] #NaN

ice <- function(x, na.rm = TRUE) {
  # Missing values
  if (any(is.na(x))) {
    if (!na.rm) return(NA)
    x[is.na(x)] <- 0
  }
  
  # Number of samples in which each species occurs
  occ <- colSums(x > 0)
  
  # Frequent (> 10 samples) and infrequent (<= 10 samples) species
  infr <- occ > 0 & occ <= 10
  S_freq <- sum(occ > 10)
  S_infr <- sum(infr)
  
  # No infrequent species: nothing to estimate
  if (S_infr == 0) return(S_freq)
  
  Q1 <- sum(occ == 1)
  N_infr <- sum(occ[infr])
  
  # Samples with at least one infrequent species
  m_infr <- sum(rowSums(x[, infr, drop = FALSE] > 0) > 0)
  
  # Sample coverage
  C_ice <- 1 - Q1 / N_infr
  if (C_ice == 0) {
    warning("All infrequent species are uniques (C_ICE = 0)")
    return(NA)
  }
  
  sum_jj <- sum(occ[infr] * (occ[infr] - 1))
  gamma2 <- max((S_infr / C_ice) * (m_infr / (m_infr - 1)) * sum_jj / N_infr^2 - 1, 0)
  
  S_freq + S_infr / C_ice + (Q1 / C_ice) * gamma2
}
custom_RICE <- rep(1,70)
for (i in 1:nrow(mite)){
  custom_RICE[i] <- ice(mite[i,])

}

ice(mite)
wiqid::richICE(t(mite)) # change to Chao2 if C_ICE = 0
fossil::ICE(t(mite)) #wrong formula
tabula::index_ice(t(mite)) # wrong formula
BAT::alpha.accum(mite)[,22] # Don't respect the purpose of the index

# RJA2 ####
SSP_RJA2 <- rep(1,70)
for (i in 1:nrow(mite)){
  SSP_RJA2[i] <- SSP:::assempar(mite[i,],type = "counts",Sest.method = "jack2")$Sest #10
}

wiqid_RJA2 <- rep(1,70)
for (i in 1:nrow(mite)){
  wiqid_RJA2[i] <- wiqid::richJackA2(t(mite[i,])) #13
}
fossil_RJA2 <- rep(1,70)
for (i in 1:nrow(mite)){
  fossil_RJA2[i] <- fossil::jack2(as.numeric(mite[i,]))
}

BiodiversityR_RJA2 <- rep(1,70)
for (i in 1:nrow(mite)){
  BiodiversityR_RJA2[i] <- as.numeric(BiodiversityR::diversityresult(mite[i,],index = "jack2",method = "each site"))
}

BAT_RJA2 <- BAT::alpha.estimate(mite)[,7]

jackknife2 <- function(x, na.rm = TRUE) {

  # Missing values
  if (any(is.na(x))) {
    if (!na.rm) return(NA)
    x[is.na(x)] <- 0
  }
  
  # Number of samples in which each species occurs
  occ <- colSums(x > 0)
  
  S_obs <- sum(occ > 0)
  Q1 <- sum(occ == 1)
  Q2 <- sum(occ == 2)
  
  S_obs + (Q1 * (2 * m - 3) / m - Q2 * (m - 2)^2 / (m * (m - 1)))
}
custom_RJA2 <- rep(1,70)
for (i in 1:nrow(mite)){
  custom_RJA2[i] <- jackknife2(mite[i,])
}

BiodiversityR::diversityresult(mite,index = "jack2",method = "pooled")
BAT::alpha.estimate(mite) # Don't follow index purpose
fossil::jack2(mite)
wiqid::richJackA2(t(mite))
SSP:::assempar(mite,type = "counts",Sest.method = "jack2")$Sest
jackknife2(mite)

# RMM ####
wiqid_RMM <- rep(1,70)
for (i in 1:nrow(mite)){
  wiqid_RMM[i] <- wiqid::richMM(t(mite[i,]))
}

michaelis_menten <- function(x, na.rm = TRUE) {

  # Missing values
  if (any(is.na(x))) {
    if (!na.rm) return(NA)
    x[is.na(x)] <- 0
  }

  # Number of samples in which each species occurs
  occ <- colSums(x > 0)
  occ <- occ[occ > 0]
  
  # Exact expected accumulation curve
  k <- seq_len(m)
  S_k <- vapply(k, function(kk) sum(1 - choose(m - occ, kk) / choose(m, kk)),
                numeric(1))
  
  # Raaijmakers maximum likelihood estimator
  Y <- S_k
  X <- S_k / k
  Sxx <- sum((X - mean(X))^2)
  Syy <- sum((Y - mean(Y))^2)
  Sxy <- sum((X - mean(X)) * (Y - mean(Y)))
  
  B    <- (mean(X) * Syy - mean(Y) * Sxy) / (mean(Y) * Sxx - mean(X) * Sxy)
  Smax <- mean(Y) + B * mean(X)
  
  Smax
}
wiqid::richMM(t(mite))
michaelis_menten(mite)

# RSQ ####
squares_richness <- function(x, na.rm = TRUE) {

  S_est <- rep(NA, nrow(x))
  
  for (k in seq_len(nrow(x))) {
    n <- x[k, ]
    
    # Missing values
    if (any(is.na(n))) {
      if (!na.rm) next
      n <- n[!is.na(n)]
    }
    
    # Species present only
    n <- n[n > 0]
    S_obs <- length(n)
    if (S_obs == 0) next
    
    N  <- sum(n)
    F1 <- sum(n == 1)
    
    # No singleton: no unseen species estimated
    if (F1 == 0) {
      S_est[k] <- S_obs
      next
    }
    
    denom <- N^2 - F1 * S_obs
    if (denom == 0) next
    
    S_est[k] <- S_obs + F1^2 * sum(n^2) / denom
  }
  
  names(S_est) <- rownames(x)
  S_est
}

tabula_RSQ <- rep(1,70)
for (i in 1:nrow(mite)){
  tabula_RSQ[i] <- tabula::index_squares(as.numeric(mite[i,]))
}
divDyn_RSQ <- rep(1,70)
for (i in 1:nrow(mite)){
  divDyn_RSQ[i] <- divDyn::indices(as.matrix(mite[i,]),method = "squares")
}
ecodive_RSQ <- rep(1,70)
for (i in 1:nrow(mite)){
  ecodive_RSQ[i] <- ecodive::squares(mite[i,]) #12.178889
}

custom_RSQ <- rep(1,70)
for (i in 1:nrow(mite)){
  custom_RSQ[i] <- squares_richness(mite[i,])
}

sum(abs(ecodive_RSQ)-abs(custom_RSQ))
sum(abs(divDyn_RSQ)-abs(custom_RSQ))
sum(abs(tabula_RSQ)-abs(custom_RSQ))


# RBOO ####
wiqid_RBOO <- rep(1,70)
for (i in 1:nrow(mite)){
  wiqid_RBOO[i] <- wiqid::richBoot(t(mite[i,])) #44
}
BiodiversityR_RBOO <- rep(1,70)
for (i in 1:nrow(mite)){
  BiodiversityR_RBOO[i] <- as.numeric(BiodiversityR::diversityresult(mite[i,],index = "boot",method = "each site")) #11
}

bootstrap_richness <- function(x, na.rm = TRUE) {

  # Missing values
  if (any(is.na(x))) {
    if (!na.rm) return(NA)
    x[is.na(x)] <- 0
  }
  
  m <- nrow(x)
  
  # Number of samples in which each species occurs
  occ <- colSums(x > 0)
  occ <- occ[occ > 0]
  
  S_obs <- length(occ)
  p <- occ / m
  
  S_obs + sum((1 - p)^m)
}

custom_RBOO <- rep(1,70)
for (i in 1:nrow(mite)){
  custom_RBOO[i] <- bootstrap_richness(mite[i,])
}

bootstrap_richness(mite)
wiqid::richBoot(t(mite))
BiodiversityR::diversityresult(mite,index = "boot",method = "pooled")

# RLMD ####
lognormal_lambda <- function(x, na.rm = TRUE) {
  R <- rep(NA, nrow(x))
  
  for (k in seq_len(nrow(x))) {
    n <- x[k, ]
    
    # Missing values
    if (any(is.na(n))) {
      if (!na.rm) next
      n <- n[!is.na(n)]
    }
    
    # Species present only
    n <- n[n > 0]
    if (length(n) < 3) next
    
    # Fit Preston's model
    fit <- tryCatch(vegan::prestondistr(n), error = function(e) NULL)
    if (is.null(fit)) next
    
    # lambda = S* / sigma = S0 * sqrt(2 * pi)
    S0 <- unname(fit$coefficients["S0"])
    if (!is.finite(S0)) next
    
    R[k] <- S0 * sqrt(2 * pi)
  }

  names(R) <- rownames(x)
  R
}

sads_RLMD <- rep(1,70)
for (i in 1:nrow(mite)){
  sads_RLMD[i] <- sads::fitsad(mite[i, ][mite[i, ] != 0],sad = "lnorm")@coef
}
asbio_RLMD <- numeric(nrow(mite))
for (i in seq_len(nrow(mite))) {
  asbio_RLMD[i] <- tryCatch(
    {
      as.numeric(
        asbio::Preston.dist(as.numeric(mite[i, ]),plot = FALSE)$Est.no.of.spp
      )
    },
    error = function(e) NA
  )
}
vegan_RLMD <- rep(1,70)
for (i in 1:nrow(mite)){
  n <- mite[i, ]
  fit <- vegan::prestondistr(n)
  vegan_RLMD[i] <- vegan::veiledspec(fit)["Extrapolated"] / fit$coefficients["width"]
}

custom_RLMD <- lognormal_lambda(mite)

sum(abs(sads_RLMD)-abs(custom_RLMD))
sum(abs(asbio_RLMD)-abs(custom_RLMD))
sum(abs(vegan_RLMD)-abs(custom_RLMD))
