#_______________________________________________________________________________
# Title              : 09_Manual_verif_Alpha.r
# Date               : 02/11/2026
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

# FUnctions to run custom functions #
check_comm_matrix <- function(x, binary = FALSE) {
  
  # 1. Check if the value is a matrix or a dataframe
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("x must be a matrix or a data.frame")
  }
  # We convert into a matrix to simplify tests
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  
  # 2. Check if the values are numeric
  if (!is.numeric(x)) {
    stop("x must contain numeric values")
  }
  
  # 3. Check for negative values
  if (any(x < 0)) {
    stop("x must contain only positive (or 0) values")
  }
  
  # 4. Check for NAs
  if (any(is.na(x))) {
    stop("x contains NA values")
  }
  
  # 5. Check for empty sites (rows)
  if (any(rowSums(x) == 0)) {
    warning("At least one site (row) contains only zeros")
  }
  
  # 6. Check for empty species (columns)
  if (any(colSums(x) == 0)) {
    warning("At least one species (column) contains only zeros")
  }
  
  # 7. If binary = FALSE (default) but data is binary, warn the user
  if (!binary && all(x %in% c(0, 1))) {
    warning("x appears to be binary (0/1 only). Some abundance-based indices
    may not be appropriate. Set binary = TRUE if this is intentional.")
  }
  
  # Return x invisibly if all checks pass
  invisible(x)
  
}

richness <- function(x, binary = FALSE) {
  check_comm_matrix(x, binary)
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
  if (!is.numeric(x)) stop("`x` must be a numeric vector.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  if (!is.numeric(q1) || !is.numeric(q2) || length(q1) != 1 || length(q2) != 1)
    stop("`q1` and `q2` must be single numeric values.")
  if (q1 <= 0 || q2 >= 1 || q1 >= q2)
    stop("Quantiles must satisfy 0 < q1 < q2 < 1.")
  
  # Keep species present only
  x <- x[!is.na(x)]
  S <- sum(x > 0)
  x <- x[x > 0]
  
  if (S < 2) {
    warning("At least 2 species are needed to compute Q.")
    return(NA_real_)
  }
  
  # Abundance classes r and number of species per class n_r
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
    return(NA_real_)
  }
  
  # Sum from R1 to R2 (inclusive, as in the formula)
  inner <- sum(nr[i1:i2])
  
  (0.5 * nr[i1] + inner + 0.5 * nr[i2]) / log(R2 / R1)
}

custom_DQ <- rep(1,70)
for (i in 1:nrow(mite)){
  custom_DQ[i] <- kempton_q(as.numeric(mite[i,]),q1 = 0.5,q2 = 0.9)
}

# DBRI ####
tabula_DBRI <- as.numeric(tabula::heterogeneity(mite,method = "brillouin"))
abdiv_DBRI <- apply(mite, 1, abdiv::brillouin_d)
ecodive_DBRI <-ecodive::alpha_div(mite,metric = "brillouin")
wiqid_DBRI <- apply(mite, 1, wiqid::biodBrillouin) #FALSE


brillouin <- function(x, binary = FALSE) {
  check_comm_matrix(x, binary)
  if (binary) {
    stop("Brillouin index requires abundance data, not binary data")
  }
  if (any(rowSums(x) == 0)) {
    warning("Brillouin index may not be meaningful when N = 0")
  }
  if (!all(x == floor(x))) {
    warning("Brillouin index requires integer counts. Non-integer values detected.")
  }
  y <- (lfactorial(rowSums(x)) - rowSums(lfactorial(x))) / rowSums(x)
  return(y)
}

custom_DBRI <- brillouin(mite)

sum(abs(wiqid_DBRI)-abs(custom_DBRI))

# DRAO ####
BAT_DRAO <- BAT::rao(as.matrix(mite))[,1] 
ade4_DRAO <- rep(1,70)
for (i in 1:nrow(mite)){
  ade4_DRAO[i] <- ade4::apqe(as.data.frame(t(mite[i,])),dis = NULL)$results[2,1]
}
SYNCSA_DRAO <- SYNCSA::rao.diversity(mite)$Simpson

rao_q <- function(x, d) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
  }
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  
  # Check distance: Euclidean `dist` object, as in ade4::apqe
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
  
  Q <- rep(NA_real_, nrow(x))
  
  for (k in seq_len(nrow(x))) {
    # Species present in sample k
    present <- which(x[k, ] > 0)
    S <- length(present)
    if (S == 0) next
    
    p   <- x[k, present] / sum(x[k, present])
    d_k <- d[present, present, drop = FALSE]
    
    # Q = sum_{i=1}^{S-1} sum_{j=i+1}^{S-1} d_ij p_i p_j
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

# EBRI ####
tabula_EBRI <- as.numeric(tabula::evenness(mite,method = "brillouin"))

brillouin_evenness <- function(x) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
  }
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  if (any(x != round(x), na.rm = TRUE))
    warning("`x` contains non-integer values: Brillouin index expects counts.")
  
  # NA treated as absence
  x[is.na(x)] <- 0
  
  E <- rep(NA_real_, nrow(x))
  
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
  
  if (any(is.na(E)))
    warning("Some samples have fewer than 2 species: E set to NA.")
  
  names(E) <- rownames(x)
  E
}

custom_EBRI <- brillouin_evenness(mite)
abs(tabula_EBRI)-abs(custom_EBRI)

# EBU ####
BAT_EBU <- as.numeric(BAT::evenness(mite,func = "bulla"))
microbiome_EBU <- microbiome::evenness(t(mite), index = "bulla")[,1] # OK

bulla_evenness <- function(x) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
  }
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  
  # NA treated as absence
  x[is.na(x)] <- 0
  
  E <- rep(NA_real_, nrow(x))
  
  for (k in seq_len(nrow(x))) {
    n <- x[k, x[k, ] > 0]
    S <- length(n)
    if (S < 2) next
    
    p <- n / sum(n)
    
    E[k] <- (sum(pmin(p, 1 / S)) - 1 / S) / (1 - 1 / S)
  }
  
  if (any(is.na(E)))
    warning("Some samples have fewer than 2 species: E set to NA.")
  
  names(E) <- rownames(x)
  E
}
custom_EBU <- bulla_evenness(mite)

sum(abs(microbiome_EBU)-abs(custom_EBU))

#ECAM ####
BAT_ECAM <- as.numeric(BAT::evenness(mite,func = "camargo"))
microbiome_ECAM <- microbiome::evenness(t(mite), index = "camargo")[,1] # OK

camargo_evenness <- function(x) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
  }
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  
  # NA treated as absence
  x[is.na(x)] <- 0
  
  E <- rep(NA_real_, nrow(x))
  
  for (k in seq_len(nrow(x))) {
    n <- x[k, x[k, ] > 0]
    S <- length(n)
    if (S == 0) next
    
    p <- n / sum(n)
    
    # Pairwise differences, pairs i < j only
    diffs <- abs(outer(p, p, "-"))
    E[k] <- 1 - sum(diffs[upper.tri(diffs)] / S)
  }
  
  if (any(is.na(E)))
    warning("Some samples are empty: E set to NA.")
  
  names(E) <- rownames(x)
  E
}

custom_ECAM <- camargo_evenness(mite)
sum(abs(microbiome_ECAM)-abs(custom_ECAM))

#EEQ ####
eq_evenness <- function(x, na.rm = TRUE) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
  }
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  
  E <- rep(NA_real_, nrow(x))
  
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
  
  if (any(is.na(E)))
    warning("Some samples have fewer than 2 species or contain NA: E set to NA.")
  
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
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
  }
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  
  E <- rep(NA_real_, nrow(x))
  
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
    
    # Population variance of log abundances (divided by S)
    log_ab <- log(n)
    var_log <- sum((log_ab - mean(log_ab))^2) / S
    
    E[k] <- 1 - (2 / pi) * atan(var_log)
  }
  
  if (any(is.na(E)))
    warning("Some samples have fewer than 2 species or contain NA: E set to NA.")
  
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

# EHEI ####
adiv_EHEI <- adiv::specieseve(mite,method = "Heip")[,1] # OK
abdiv_EHEI <- rep(1,70) 
for (i in 1:nrow(mite)){
  abdiv_EHEI[i] <- abdiv::heip_e(as.numeric(mite[i,])) # OK
}

heip_evenness <- function(x, na.rm = TRUE) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
  }
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  
  E <- rep(NA_real_, nrow(x))
  
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
    
    # Shannon entropy (natural log)
    p <- n / sum(n)
    H <- -sum(p * log(p))
    
    E[k] <- (exp(H) - 1) / (S - 1)
  }
  
  if (any(is.na(E)))
    warning("Some samples have fewer than 2 species or contain NA: E set to NA.")
  
  names(E) <- rownames(x)
  E
}

custom_EHEI <- rep(1,70) 
for (i in 1:nrow(mite)){
  custom_EHEI[i] <- heip_evenness(mite[i,]) 
}

sum(abs(abdiv_EHEI)-abs(custom_EHEI))

# EMC ####
tabula_EMC <- as.numeric(tabula::evenness(mite,method = "mcintosh"))
abdiv_EMC <- apply(mite,1,abdiv::mcintosh_e)

mcintosh_evenness <- function(x, na.rm = TRUE) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
  }
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  
  E <- rep(NA_real_, nrow(x))
  
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
  
  if (any(is.na(E)))
    warning("Some samples have fewer than 2 species or contain NA: E set to NA.")
  
  names(E) <- rownames(x)
  E
}


custom_EMC <- mcintosh_evenness(mite)
sum(abs(tabula_EMC)-abs(custom_EMC))
# EMN ####

microbiome_EMN <- rep(1,70) 
for (i in 1:nrow(mite)){
  microbiome_EMN[i] <- microbiome::dominance(as.numeric(mite[i,]),index = "DMN")$dmn  #OK
}


mcnaughton_dominance <- function(x, na.rm = TRUE) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
  }
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  
  d <- rep(NA_real_, nrow(x))
  
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
    
    d[k] <- (Nmax1 + Nmax2) / N
  }
  
  if (any(is.na(d)))
    warning("Some samples are empty or contain NA: d set to NA.")
  
  names(d) <- rownames(x)
  d
}

custom_EMN <- mcnaughton_dominance(mite)

sum(abs(microbiome_EMN)-abs(custom_EMN))

# ESP ####

simpson_evenness <- function(x, na.rm = TRUE) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
  }
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  
  E <- rep(NA_real_, nrow(x))
  
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
  
  if (any(is.na(E)))
    warning("Some samples are empty or contain NA: E set to NA.")
  
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

# EST ####
strong_dominance <- function(x, na.rm = TRUE) {
  # Convert input to a samples x species matrix
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(dim(x))) {
    x <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
  }
  if (!is.numeric(x)) stop("`x` must contain numeric values only.")
  if (any(x < 0, na.rm = TRUE)) stop("`x` must not contain negative values.")
  
  D <- rep(NA_real_, nrow(x))
  
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
    
    D[k] <- max(t_i / N - i / S)
  }
  
  if (any(is.na(D)))
    warning("Some samples have fewer than 2 species or contain NA: D set to NA.")
  
  names(D) <- rownames(x)
  D
}

custom_EST <- strong_dominance(mite)
abdiv_EST <- rep(1,70) 
for (i in 1:nrow(mite)){
  abdiv_EST [i] <- abdiv::strong(as.numeric(mite[i,]))
}

sum(abs(abdiv_EST)-abs(custom_EST))

# TRY ####
rygg <- function(x, binary = FALSE) {
  check_comm_matrix(x, binary)
  if (binary) {
    stop("Rygg index requires abundance data, not binary data")
  }
  if (any(rowSums(x) == 0)) {
    warning("Rygg index may not be meaningful when N = 0")
  }
  if (any(rowSums(x) <= exp(1))) {
    warning("Rygg index may not be meaningful when N <= e (log(log(N)) <= 0)")
  }
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
  check_comm_matrix(x, binary)
  if (binary) {
    stop("Rygg adjusted index requires abundance data, not binary data")
  }
  if (any(rowSums(x) == 0)) {
    warning("Rygg adjusted index may not be meaningful when N = 0")
  }
  if (any(rowSums(x) <= exp(1))) {
    warning("Rygg adjusted index may not be meaningful when N <= e (log(log(N)) <= 0)")
  }
  y <- log(richness(x)) / log(log(rowSums(x)+1)+1)
  return(y)
}
custom_TRYA <- rygg_adjusted(mite)

sum(abs(benthos_TRYA)-abs(custom_TRYA))
