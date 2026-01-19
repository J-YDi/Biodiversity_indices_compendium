#_______________________________________________________________________________
# Nom               : Testing_package_div_index.r
# Date de modif     : 30/09/2025
# Objet             : Tester les fonctions pour calculer les indices de diversite
# Auteurs           : J-Y. Dias
# Version R         : 4.5.0
#_______________________________________________________________________________

#________________Utiliser les donnees vegan comme test__________________________####
library(vegan)
# Varespec contient 44 abondances d'espece sur 24 echantillons
data(varespec)
# Pour eviter tout potentiel conflit :
detach(package:vegan)


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


x <- varespec[1,]
x_relative <- varespec[1,]/sum(varespec[1,])
x_entier <- round(x,digits = 0)

varespec_long <- tidyr::pivot_longer(varespec,cols = colnames(varespec)[1]:colnames(varespec)[44],names_to = "Esp",values_to = "Value")

varespec_long <- varespec_long[1:44,]

species_counts <- colSums(varespec)

freq_table <- table(species_counts)

freq_df <- data.frame(
  j = as.numeric(names(freq_table)),
  n_j = as.numeric(freq_table)
)

varespec_pa <- convert_to_presence_absence(varespec)

varespec <- t(varespec)
species_occurrences <- rowSums(varespec > 0)  # nombre de sites par espèce
max_occ <- max(species_occurrences)
f <- tabulate(species_occurrences, nbins = max_occ)

f
# Exemple : f[1] = nb d'espèces présentes dans 1 site, etc.

# Étape 2 : Reconstruire une matrice respectant ce vecteur
n_samples <- ncol(varespec)
n_species <- sum(f)

# Créer une matrice vide
new_varespec <- matrix(0, nrow = n_species, ncol = n_samples)

# Remplir la matrice
set.seed(123)
row_idx <- 1
for (i in seq_along(f)) {
  for (j in 1:f[i]) {
    cols <- sample(n_samples, i)
    new_varespec[row_idx, cols] <- 1
    row_idx <- row_idx + 1
  }
}

# Ajouter noms
rownames(new_varespec) <- paste0("sp", 1:n_species)
colnames(new_varespec) <- colnames(varespec)

# Convertir en data.frame
new_varespec <- as.data.frame(new_varespec)


#________________________Indices de diversite alpha_____________________________####

#### Pielou ####----------------------------------------------------------------
OTUtable::pielou(x)
abdiv::pielou_e(as.numeric(x_relative))

tabula::evenness(x,method = "shannon")

BiodiversityR::diversityresult(varespec,y=NULL,index="Jevenness",method = "each site")

install.packages("forestmangr")
forestmangr::species_diversity(varespec_long,species = "Esp",index = "all") # incorrect pour tous les indices

install.packages("adiv")
adiv::specieseve(x,method = "Shannon")

install.packages("pctax")
pctax::a_diversity(round(t(x),digits = 0),method = "pielou")

install.packages("diverse")
diverse::diversity(t(x),type = "evenness",category_row = T)

microbiome::evenness(t(varespec), index = "pielou") # OK

chemodiv::calcDiv(x,type = "PielouEven") # OK

sprex::diversity(as.numeric(x),type = "eveness.pielou")

breakaway::true_shannon_e(x_relative) # Pielou
#### Berger-Parker ####---------------------------------------------------------
abdiv::berger_parker_d(x)

divDyn::indices(as.matrix(x),method = "dominance")

tabula::index_berger(as.numeric(x))
tabula::heterogeneity(x,method = c("berger"))

BiodiversityR::diversityresult(varespec,y=NULL,index="Berger",method = "each site")

diverse::diversity(t(x),type = "berger-parker",category_row = T)$berger.parker.D

install.packages("agricolae")

agricolae::index.bio(x,method = "Berger.Parker")$index

install.packages("divDyn")
divDyn::indices(varespec_long$Esp,varespec_long$Value) # faux

ecodive::alpha_div(round(x,digits = 0),metric="berger")

wiqid::biodBerger(abVec=varespec) # faux

microbiome::dominance(as.numeric(x),index = "DBP") # Berger-Parker

triversity::get_diversity_from_distribution(as.numeric(x_relative),measure = "bergerparker" ) # OK



#### Hulburt ####---------------------------------------------------------------
hulburt_index <- function(x) {
  s <- sum(x, na.rm = TRUE)
  if (s == 0) return(NA_real_)
  sum(sort(x, decreasing = TRUE)[1:2], na.rm = TRUE) / s
}

hulburt_index(as.numeric(x))

#### Gini-Simpson - redite Simpson ####------------------------------------------------------------------

install.packages("catsim")

#divseg::ds_gini(varespec,.cols = dplyr::everything()) # marche pas tres bien

lawstat::gini.index(as.numeric(x)) # faux

CUB::gini(as.numeric(x_relative)) # faux

RoughSets::X.gini(as.numeric(x)) #OK

breakaway::true_gini(x_relative) #OK

concstats::concstats_gini(as.numeric(x_relative)) #faux

catsim::gini(as.numeric(x)) # faux

catsim::sqrtgini(as.numeric(x)) # faux

ade4::divc(as.data.frame(t(varespec)),dis = NULL) # 115

DescTools::DivCoef(as.data.frame(t(varespec)),dis = NULL) # 115

adiv::speciesdiv(x,method = "GiniSimpson") #115

diverse::diversity(t(x),type = "gini-simpson",category_row = T)$gini.simpson #115

ecodive::alpha_div(x,metric = "simpson") #115

agricolae::index.bio(x,method = "Simpson.Div")$index

entropart::GenSimpsonD(as.numeric(x_relative),Correction ="None")

#### Patten ####

# patten_index : calcule l'indice de Patten (R) à partir d'une matrice/data.frame d'abondances
# Utilise vegan::diversity() pour H' (Shannon).

patten_index <- function(comm, base = exp(1)) {
  
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Le package 'vegan' est requis")
  }
  # calcul de H' (Shannon) par échantillon
  H <- vegan::diversity(comm, index = "shannon", base = base)
  
  # Hmax et Hmin (parmi tous les échantillons fournis)
  Hmax <- max(H, na.rm = TRUE)
  Hmin <- min(H, na.rm = TRUE)
  
  # si Hmax == Hmin (pas de variation), on evite division par 0 : on renvoie NA
  denom <- Hmax - Hmin
  if (is.na(denom) || denom == 0) {
    warning("Hmax equals Hmin (no variation in Shannon)")
    R <- rep(NA_real_, length(H))
  } else {
    R <- (Hmax - H) / denom
  }
  
  # retourne un data.frame utile
  out <- R
  return(out)
}

patten_index(varespec)

#### Hurlbert ####--------------------------------------------------------------

benthos::hurlbert(taxon = colnames(x),count = as.numeric(x), n=2) # OK

entropart::Hurlbert(as.numeric(x),k = 2) # OK
entropart::HurlbertD(as.numeric(x),k = 2) # correspond à 1/Hulbert pour faire equivalence de Inv Simpson

tabula::index_hurlbert(as.numeric(x),sample = 2) # OK

benthos::hpie(taxon = colnames(x),count = as.integer(x)) # Hulbert PIE

install.packages("BAT")
mobsim::spec_sample(as.numeric(x),n=2) # OK

mobr::calc_div(as.integer(x),index="PIE",effort = 2) # Hulbert PIE

vegan::simpson.unb(as.integer(x)) # Hulbert PIE 

vegan::rarefy(as.integer(x),sample = 2)

BiodiversityR::diversityresult(as.integer(x),y=NULL,index="simpson.unb",method = "each site") # Hulbert PIE = unbiaised Simpson

divent::div_hurlbert(as.integer(x),k = 2,estimator = "Hurlbert") # ne correspond pas

#### Bulla ####-----------------------------------------------------------------
BAT::evenness(x,func = "bulla") # 0.28191063 OK

bulla_index <- function(comm,
                        margin = 1,
                        na.rm = TRUE,
                        include_zero = FALSE,
                        min_spp = 2L,
                        raw = FALSE,
                        is.proportion = NULL) {
  
  # helper: compute for one numeric vector
  .one <- function(v) {
    if (!is.numeric(v)) stop("Abundances must be numeric.")
    if (na.rm) v <- v[!is.na(v)] else if (any(is.na(v))) return(NA_real_)
    
    if (length(v) == 0) return(NA_real_)
    if (any(v < 0)) stop("Abundances must be non-negative.")
    
    # detect if input are proportions
    if (is.null(is.proportion)) {
      tol <- 1e-8
      is_prop <- all(v <= 1 + tol) && abs(sum(v, na.rm=TRUE) - 1) < 1e-6
    } else is_prop <- is.proportion
    
    if (is_prop) {
      p <- v / sum(v)    # keep stable normalization even if sum ~ 1
    } else {
      total <- sum(v)
      if (total <= 0) return(NA_real_)  # no individuals -> undefined
      p <- v / total
    }
    
    if (!include_zero) {
      p_sub <- p[p > 0]
    } else {
      p_sub <- p
    }
    
    S <- length(p_sub)
    if (S < min_spp) return(NA_real_)
    
    O <- sum(pmin(p_sub, 1 / S))
    if (raw) return(O)
    
    # normalized Bulla index
    denom <- 1 - 1 / S
    # numerical safety (denom > 0 because S >= 2 here)
    E_O <- (O - 1 / S) / denom
    # clamp to [0,1] for floating errors
    E_O <- pmax(0, pmin(1, E_O))
    return(as.numeric(E_O))
  }
  
  # If single vector
  if (is.null(dim(comm))) {
    return(.one(as.numeric(comm)))
  }
  
  # matrix / data.frame
  if (is.data.frame(comm)) comm <- as.matrix(comm)
  if (!is.matrix(comm)) stop("comm must be numeric vector, matrix or data.frame.")
  
  if (!(margin %in% c(1,2))) stop("margin must be 1 (rows) or 2 (columns).")
  apply(comm, margin, function(x) .one(as.numeric(x)))
}

bulla_index(x) #OK

microbiome::evenness(t(varespec), index = "bulla") # OK

#### Camargo ####---------------------------------------------------------------
# Pas possible de determiner la valeur correcte
BAT::evenness(x,func = "camargo") # 0.94280965
microbiome::evenness(t(varespec), index = "camargo") # 0.19641256

camargo_evenness <- function(comm,
                             margin = 1,
                             na.rm = TRUE,
                             include_zero = FALSE) {
  
  # fonction interne (calcul pour un vecteur)
  .one <- function(v) {
    v <- as.numeric(v)
    if (na.rm) v <- v[!is.na(v)] else if (any(is.na(v))) return(NA_real_)
    if (length(v) == 0) return(NA_real_)
    if (any(v < 0)) stop("Les abondances doivent être ≥ 0.")
    
    total <- sum(v)
    if (total <= 0) return(NA_real_)
    p <- v / total
    
    if (!include_zero) p <- p[p > 0]
    S <- length(p)
    if (S < 2) return(NA_real_)
    
    # somme des différences absolues entre toutes les paires
    diff_sum <- 0
    for (i in 1:(S-1)) {
      for (j in (i+1):S) {
        diff_sum <- diff_sum + abs(p[i] - p[j])
      }
    }
    E <- 1 - diff_sum / S
    # borne numérique dans [0,1]
    return(pmax(0, pmin(1, E)))
  }
  
  # si vecteur simple
  if (is.null(dim(comm))) {
    return(.one(comm))
  }
  
  # matrice / data.frame
  if (is.data.frame(comm)) comm <- as.matrix(comm)
  if (!is.matrix(comm)) stop("comm doit être un vecteur, une matrice ou un data.frame.")
  if (!(margin %in% c(1,2))) stop("margin doit être 1 (lignes) ou 2 (colonnes).")
  
  apply(comm, margin, .one)
}


camargo_evenness(x) # 0.19933509 

EcoIndR::DER(varespec) # pas fonctionnel

#### Heip ####------------------------------------------------------------------

adiv::specieseve(x,method = "Heip") # OK
abdiv::heip_e(as.numeric(x)) # OK

OnomasticDiversity::fHeip(varespec_long,k= "Value",n="Value",location = "Site",s = "Ntot") # pas fonctionnel

heip_evenness <- function(x, na.rm = TRUE, include_zero = FALSE) {
  x <- as.numeric(x)
  if (na.rm) x <- x[!is.na(x)] else if (any(is.na(x))) return(NA_real_)
  if (!include_zero) x <- x[x>0]
  S <- length(x)
  if (S < 2) return(NA_real_)
  p <- x / sum(x)
  H <- -sum(p * log(p))
  E <- (exp(H) - 1)/(S - 1)
  E <- pmax(0, pmin(1, E))
  return(E)
}
heip_evenness(x)

#### Herfindahl-Hirschman ####--------------------------------------------------
# Mathematiquement equivalent a Simpson 

EconGeo::herfindahl(varespec) # 0.17828850

concstats::concstats_hhi(as.numeric(x)) # 0.1782885

REAT::herf(x) # 0.1782885

varespec_long$Relative <- as.numeric(varespec_long$Value/sum(varespec_long$Value))
hhi::hhi(as.data.frame(varespec_long),"Relative") # OK mais ne pas transformer en x100 !

politicsR::hh(as.numeric(x_relative)) # DOIT ETRE EN RELATIF sinon Faux

PDtoolkit::hhi(as.numeric(x_relative)) # A relancer

antitrust::HHI(as.numeric(x_relative)) # Faux, multiplicateur

DescTools::Herfindahl(as.numeric(x)) # OK

divseg::ds_hhi(varespec,.cols = dplyr::everything()) # OK

ineq::Herfindahl(x) # OK

triversity::get_diversity_from_distribution(as.numeric(x_relative),measure = "herfindahl" ) # OK

diverse::diversity(t(x),type = "herfindahl-hirschman",category_row = T) # OK

PCRA::divHHI(x_relative) # 1-HHI = 1-Simpson = Gini-Simpson

#### Rényi ####-----------------------------------------------------------------
EntropyEstimation::Renyi.z(as.numeric(x),r=1) # 1.8219059 normalisation probabiliste différente

EntropyEstimation::Renyi.sd(as.numeric(x),r=2) #0.8170176 normalisation probabiliste différente

EntropyEstimation::RenyiEq.z(as.numeric(x),r=2) #0.19166667 normalisation probabiliste différente ne correspond aps

sprex::diversity(as.numeric(x),type = "renyi",q=2) #OK

statcomp::permutation_entropy_Renyi(x,1) # 0.53320883 ne correspond pas

vegan::renyi(x,scales = c(0,1,2),hill = F)  # 3.367296 2.017763 1.724352 correct 

divo::rd(as.matrix(x),alpha = 2) # marche pas bien

seewave::sh(x,alpha = "shannon") # 0.53320891 Shannon en base 2 (bits), pas en nats
seewave::sh(x,alpha = "simpson") # 0.82171151 ne donne pas la forme log
seewave::sh(x,alpha = 2) # 0.4556728 ne donne pas la forme log

BiodiversityR::renyiresult(varespec,y=NULL,method = "each site",scales = c(0,1,2)) # correct

adiv::divparam(x,method = "renyi",q=c(0,1,2))$div # correct

OnomasticDiversity::fHill(varespec_long,k= "Value",n="Notot",location = "Site",lambda = 0) # pas fonctionnel

renyi_entropy <- function(x, alpha = 1, na.rm = TRUE, include_zero = FALSE) {
  # x : vecteur d'abondances
  # alpha : paramètre d'ordre (>0)
  
  # Nettoyage
  x <- as.numeric(x)
  if (na.rm) x <- x[!is.na(x)] else if (any(is.na(x))) return(NA_real_)
  if (!include_zero) x <- x[x > 0]
  N <- sum(x)
  if (N == 0) return(NA_real_)
  
  # Proportions
  p <- x / N
  
  # Cas alpha = 1 (limite -> Shannon)
  if (abs(alpha - 1) < .Machine$double.eps^0.5) {
    return(-sum(p * log(p)))
  }
  
  # Formule générale de Rényi
  H <- (1 / (1 - alpha)) * log(sum(p^alpha))
  return(H)
}

renyi_entropy(x, alpha = 0)   # OK

diverse::diversity(t(x),type = "renyi",category_row = T) # OK

vegan::renyi(x,scales = c(0,1,2),hill = F)  # 3.367296 2.017763 1.724352

adiv::divparam(x,method = "renyi",q=c(0,1,2))# 3.367296 2.017763 1.724352

EntropyEstimation::Renyi.z(as.numeric(x),r=1) # 1.821906 

diverse::diversity(t(varespec),type = "renyi",category_row = T,q=0) # OK

#### NHC : Nee-Harvey-Cotgreave ####--------------------------------------------
nhc_original <- function(x, na.rm = TRUE) {
  # x : vecteur d'abondances
  x <- as.numeric(x)
  
  # Nettoyage des NA
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) < 2) return(NA_real_)
  
  # Ne garder que les abondances > 0
  x <- x[x > 0]
  if (length(x) < 2) return(NA_real_)
  
  # Tri décroissant
  x_sorted <- sort(x, decreasing = TRUE)
  rank <- 1:length(x_sorted)
  
  # Régression linéaire log(abondance) ~ rang
  model <- lm(log(x_sorted) ~ rank)
  
  # Extraction de la pente
  b <- coef(model)[2]
  return(b)
}

nhc_original(x)

#### NHC modifie EQ/E0 Smith & Wilson 1996 ####---------------------------------
e0_index <- function(x, na.rm = TRUE) {
  # x : vecteur d'abondances
  x <- as.numeric(x)
  
  # Nettoyage des NA et des 0
  if (na.rm) x <- x[!is.na(x)]
  x <- x[x > 0]
  
  if (length(x) < 2) return(NA_real_)
  
  # Tri décroissant
  x_sorted <- sort(x, decreasing = TRUE)
  rank <- 1:length(x_sorted)
  n <- length(x_sorted)
  
  # Régression linéaire log(abondance) ~ rang
  model <- lm(log(x_sorted) ~ rank)
  b <- coef(model)[2]
  
  # Transformation normalisée pour obtenir E0 / EQ
  E0 <- -2 / n * atan(b)
  
  return(E0)
}
e0_index(x) 

codyn::community_structure(varespec_long,abundance.var = "Value",metric = "EQ")$EQ

evar_index <- function(x, na.rm = TRUE) {
  x <- as.numeric(x)
  
  # Nettoyage
  if (na.rm) x <- x[!is.na(x)]
  x <- x[x > 0]
  n <- length(x)
  if (n < 2) return(NA_real_)
  
  # Variance populationnelle (divisé par n)
  mean_logx <- mean(log(x))
  var_logx <- sum((log(x) - mean_logx)^2) / n  # populationnelle
  
  # Transformation arctan
  Evar <- 1 - (2/pi) * atan(var_logx)
  return(Evar)
}
evar_index(x) # OK S&W Eveness

codyn::community_structure(varespec_long,abundance.var = "Value",metric = "Evar")$Evar # Ok mais variance divisé par n

adiv::specieseve(x,method = "SmithWilson") # OK S&W Evenness -> divisé par n-1

microbiome::evenness(t(varespec), index = "evar")

#### McNaughton ####------------------------------------------------------------
mcnaughton_index <- function(x, na.rm = TRUE) {
  x <- as.numeric(x)
  if (na.rm) x <- x[!is.na(x)]
  x <- x[x > 0]
  
  if (length(x) == 0) return(NA_real_)
  
  # Total et 2 espèces dominantes
  N <- sum(x)
  top2 <- sort(x, decreasing = TRUE)[1:min(2, length(x))]
  
  D <- sum(top2) / N
  return(D)
}

mcnaughton_index(x)  #OK

microbiome::dominance(as.numeric(x),index = "DMN") # McNaughton

#### Sheldon ####---------------------------------------------------------------
OnomasticDiversity::fSheldon(varespec_long,k= "Value",n="Value",location = "Site",s = "Ntot") # pas fonctionnel

sheldon_index <- function(x, base = exp(1), na.rm = TRUE) {
  x <- as.numeric(x)
  if (na.rm) x <- x[!is.na(x)]
  x <- x[x > 0]
  N <- sum(x)
  if (N == 0 || length(x) < 2) return(NA_real_)
  
  p <- x / N
  H <- -sum(p * log(p, base = base))
  
  S <- length(p)
  Sh <- exp(H) / S
  return(Sh)
}
sheldon_index(x) # OK

chemodiv::calcDiv(x,type = "HillEven",q=1)# OK

adiv::eveparam(x,method = "hill",q=1) # OK

BiodiversityR::diversityresult(varespec,y=NULL,index="Eevenness",method = "each site")

#### Strong Dominance Index ####------------------------------------------------
abdiv::strong(as.numeric(x))


#### Kothe ####-----------------------------------------------------------------
# Rien

kothe <- function(matrix_data) {
  S_i <- apply(matrix_data, 1, function(x) sum(x > 0))
  S_max <- max(S_i)
  (S_max - S_i) / S_max
}
kothe(varespec)

#### Redundancy =? Patten ####--------------------------------------------------
redundancy_index <- function(abundances) {
  # abundances : vecteur d'abondance pour un échantillon
  
  # Calcul des composantes
  S <- sum(abundances > 0)       # nombre d'espèces
  N <- sum(abundances)           # nombre total d'individus
  H <- vegan::diversity(abundances, index = "shannon")  # indice de Shannon H'
  
  # Vérification des conditions minimales
  if (S < 2 || N <= S) return(NA)
  
  # Calcul du numérateur et du dénominateur
  num <- log(S) - H
  denom <- log(S) - (1/N) * (lfactorial(N) - lfactorial(N - S + 1))
  
  R <- num / denom
  return(R)
}
redundancy_index(x)


#### Tsallis ###----------------------------------------------------------------
vegan::tsallis(x,scales = c(0,1,2),hill = F) # 28.0000000  2.0177633  0.8217115 OK

EntropyEstimation::Tsallis.z(x,2) # proche mais different

adiv::divparam(x,method = "tsallis",q=c(0,1,2)) # OK

entropart::Tsallis(as.numeric(x),q=2) # OK

adiv::abgdivparam(comm = x,"speciesab",method = "tsallis",q=2)[1] # OK

divent::ent_tsallis(as.numeric(x),q=1,probability_estimator = "naive",richness_estimator = "naive")$entropy # OK

#### Tsallis evenness ####------------------------------------------------------
adiv::eveparam(x,method = "tsallis",q=c(1,2,3))

#### Blau = Gini-Simpson ####---------------------------------------------------
diverse::diversity(t(x),type = "blau",category_row = T)

#### Rarity ####----------------------------------------------------------------
Rarity::Irr(t(x),t(x),abundance = T) # il faut des poids qui ne s'applique pas


#### Ludwig-Reynold ####
Ludwig_Reynold <- function(abondances) {
  N <- sum(abondances)
  S <- sum(abondances > 0)  # Nombre d'espèces présentes
  p_i <- abondances / N
  p_i <- p_i[p_i > 0]
  
  H_prime <- -sum(p_i * log(p_i))
  
  if (S <= 1) {
    return(NA)  # Évitons une division par zéro
  }
  
  E <- exp(H_prime - 1) / (S - 1)
  return(E)
}


#### Evenness derives ####------------------------------------------------------
tabula::evenness(x,method = "shannon") # = Pielou deja indique
tabula::evenness(x,method = "simpson") # Inv Simpson E
tabula::evenness(x,method = "brillouin") # Brillouin E
tabula::evenness(x,method = "mcintosh") # McIntosh E deja indique
sprex::diversity(as.numeric(x),type = "eveness.simpson") # faux



codyn::community_structure(varespec_long,abundance.var = "Value",metric = "SimpsonEvenness") # Inv Simpson E
adiv::specieseve(x,method = "GiniSimpson") # Gini Simpson E
adiv::specieseve(x,method = "Simpson") # Inv Simpson E
adiv::specieseve(x,method = "Shannon") # Pielou
adiv::specieseve(x,method = "SmithWilson") # Smith & Wilson Evar deja indique

abdiv::simpson_e(x) # Inv Simpson E
microbiome::evenness(t(varespec), index = "simpson") # Inv Simpson E

#### Gini ####------------------------------------------------------------------
DescTools::Gini(as.numeric(x),unbiased = F) # Gini

microbiome::dominance(as.numeric(x),index = "gini") # Gini

giniVarCI::igini(as.numeric(x),bias.correction = F)

acid::gini(as.numeric(x))$Gini 

dplR::gini.coef(as.numeric(x)) 

shipunov::Gini(as.numeric(x)) 

EconGeo::gini(t(varespec)) 

ineq::Gini(x) 

wINEQ::Gini(as.numeric(x)) #
REAT::gini(as.numeric(x)) #

#### Hoover index ####----------------------------------------------------------
EconGeo::hoover_index(t(x))

wINEQ::Hoover(as.numeric(x)) # la même chose mais pas /100

REAT::hoover(as.numeric(x))

#### Good ####

good_index <- function(p, m, n) {
  if (any(p <= 0 | p > 1)) {
    stop("Les probabilités doivent être dans l'intervalle (0, 1].")
  }
  if (m <= 0 || n <= 0) {
    stop("m et n doivent être des entiers naturels non nuls.")
  }
  return(sum(p^m * (-log(p))^n))
}

#### Hurlbert E ####------------------------------------------------------------
indice_hurlbert <- function(abondances) {
  N <- sum(abondances)
  S <- sum(abondances > 0)
  p_i <- abondances / N
  p_i <- p_i[p_i > 0]
  
  H_prime <- -sum(p_i * log(p_i))
  H_max <- log(S)
  
  if ((N - S + 1) <= 0) {
    return(NA)  # Cas non défini
  }
  
  H_min <- log(N) - ((N - S + 1) * log(N - S + 1)) / N
  
  if ((H_max - H_min) == 0) {
    return(NA)  # Évitons une division par zéro
  }
  
  E_hurlbert <- (H_prime - H_min) / (H_max - H_min)
  return(E_hurlbert)
}

#_________________Estimateurs de la richesse specifique_________________________####

#### Abundance-based coverage estimator ACE ####--------------------------------

tabula::index_ace(as.numeric(x_entier)) #12.1788889
fossil::ACE(x_entier,taxa.row = F) # 11.933333

BAT::alpha.estimate(x_entier)[11] #  12.178889

RSE::DetAbu(varespec) # marche pas

RSE::DetInt(varespec) # marche pas

iNEXT::DataInfo(t(varespec),datatype = "abundance") # pas l'air top

SPECIES::ChaoLee1992(freq_matrix) # marche pas bien

vegan::estimateR(x_entier,index = "chao")[4]  #12.1788889

?ecodive::ace(x_entier) #12.17889 

wiqid::richACE(as.numeric(x_entier)) #12.178889

pctax::a_diversity(t(x_entier),method = "ace") # 12.178889

fossil::spp.est(t(varespec),abund = T) #26.837758

sprex::ACE(as.numeric(new_varespec[1,])) # marche pas 


#### Incidence-based coverage estimator ####------------------------------------
# Tab1 non remplie a partir de la
fossil::ICE(t(varespec_pa)) #44

wiqid::richICE(t(varespec_pa)) #44

tabula::index_ice(as.matrix(varespec_pa)) #44

wiqid::richICE(t(varespec_pa)) #44

fossil::spp.est(t(varespec_pa),abund = F) #388.200000

BAT::alpha.accum(x_entier) #NaN

#### Chao1 ####-----------------------------------------------------------------
tabula::index_chao1(as.numeric(x_entier)) # 11.659176

vegan::estimateR(x_entier,index = "chao")[2]  #11.2500000

BiodiversityR::diversityresult(x_entier,index = "chao",method = "each site") #11

ecodive::chao1(x_entier) # 11.666667

microbiome::richness(t(x_entier),index = c("chao1")) #11.25

mobr::calc_chao1(x_entier) #11.659176

pctax::a_diversity(t(x_entier),method = "chao1") #11.25

breakaway::chao1(x_relative) #marche pas

fossil::chao1(as.numeric(x_entier)) # 11.666667

wiqid::richChao1(t(x_entier)) #11.666667

entropart::bcRichness(x_entier,Correction = "Chao1") # 11.7

divent::div_richness(as.numeric(t(x_entier)),estimator = "Chao1") #  11.7

OTUtable::chao1(as.numeric(x_entier)) #  11.666667

print(rareNMtests::chao1(x_entier)) #11

breakaway::chao_shen(x_relative) # marche pas

SpadeR::ChaoSpecies(varespec,datatype = "abundance") # marche pas

iNEXT::ChaoRichness(t(varespec),datatype = "abundance") #29

BAT::alpha.estimate(x_entier) #  11.25 chao 11.6219 chao corrected

BAT::alpha.accum(x_entier) #  11.25 chao 11.6219 chao corrected

vegan::specpool(x_entier,smallsample = F,pool = T) #66

SSP:::assempar(x_entier,type = "counts",Sest.method = "chao") #10

biosampleR::calc_diversity_indices(x) #29

#### Chao2 ####-----------------------------------------------------------------
fossil::chao2(as.numeric(x)) #29

wiqid::richChao2(t(varespec_pa)) #44

install.packages("rareNMtests")

breakaway::chao1(as.matrix(varespec_pa)) # marche pas

rarestR::es(varespec,m=44) # methode non indiquee

SPECIES::chao1984(freq_matrix)  # marche pas bien

tabula::index_chao2(as.matrix(varespec_pa)) #44

divDyn::indices(as.matrix(x),method = "chao2") #29

BAT::alpha.accum(x_entier) #66 chao 2 et corrected

#### Jackknife 1 ####-----------------------------------------------------------

SPECIES::jackknife() # marche pas bien

divent::div_richness(as.numeric(t(x_entier)),estimator = "jackknife") #  13

tabula::jackknife(tabula::heterogeneity(varespec,method = "shannon")) # ne permet pas de le faire sur la rspe

Omisc::jackknife(varespec) # le fait sur les donnees, ne permet pas de le faire sur rspe

arkhe::jackknife(as.numeric(x),do = sum) #ne correspond pas

asbio::pseudo.v(x,mean) # ne correspond pas

BiodiversityR::diversityresult(x_entier,index = "jack1",method = "each site") #11

entropart::bcRichness(x_entier,Correction = "Jackknife") # 13

BAT::alpha.estimate(x_entier) #  13

BAT::alpha.accum(x_entier) # 13

fossil::jack1(as.numeric(x_entier)) #12.977528

vegan::estimateR(x_entier,index = "jack1") # marche pas correspond a chao

wiqid::richJackA1(t(x_entier)) #13

SSP:::assempar(x_entier,type = "counts",Sest.method = "jack1") #10

#### Jackknife 2 ####-----------------------------------------------------------

SSP:::assempar(x_entier,type = "counts",Sest.method = "jack2") #10

wiqid::richJackA2(t(x_entier)) #12

vegan::estimateR(x_entier,index = "jack2") # marche pas correspond a chao

fossil::jack2(as.numeric(x_entier)) #11.965909

BiodiversityR::diversityresult(x_entier,index = "jack2",method = "each site") #11

BAT::alpha.estimate(x_entier) #12

BAT::alpha.accum(x_entier) # 12

install.packages("divDyn")


#### Squares Richness Estimator ####--------------------------------------------
tabula::index_squares(as.numeric(x_entier)) #11.72161

divDyn::indices(as.matrix(x_entier),method = "squares") #11.72161 

ecodive::squares(x_entier) #12.178889

#### Michaelis-Menten, Clench ####----------------------------------------------
wiqid::richMM(t(varespec)) #7.6426636
sprex::Clench(as.numeric(x)) #86

#### Bootstrap ####-------------------------------------------------------------
vegan::specpool2vect(varespec,index="boot") # marche pas
wiqid::richBoot(t(x_entier)) #44
BiodiversityR::diversityresult(x_entier,index = "boot",method = "each site") #11
fossil::bootstrap(varespec,abund = T) # 18, 24

#________________________Indices de diversite béta______________________________####
library(vegan)
# Varespec contient 44 abondances d'espece sur 24 echantillons
data(varespec)
# Pour eviter tout potentiel conflit :
detach(package:vegan)

x1 <- varespec[1,]
x2 <- varespec[2,]

varespec_2 <- varespec * 2

x12 <- varespec[1:2,]

x11 <- rbind(x1,x2)

x12_pa <- ifelse(x12 != 0, 1, 0)
x11_pa <-ifelse(x11 != 0, 1, 0)

x1int <- as.integer(varespec[1,])
x2int <- as.integer(varespec[2,])
x12int <- rbind(x1int,x2int)

#### Euclidean distance ####----------------------------------------------------

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

?NST::beta.g(x12,dist.method = "mEuclidean") #1.223445 modif par Andersen et al. 2006

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

#### Squared euclidean distance ####--------------------------------------------
ecodive::squared_euclidean(x12,rescale = F) #1630.034

neighbr::distance(as.numeric(x1),as.numeric(x2),measure = "squared_euclidean") #1630.034

NMFN::distance2(x1,x2) #1630.034 squared, pas indiqué comme tel

MultivariateAnalysis::Distancia(x12,Metodo = 3)[1] #37.04623

ptm::pairwise.dist(as.matrix(x1),as.matrix(x2),squared = T) #1630.034

codep::Euclid(x1,x2,squared = T) #1630.034

M2SMJF::dist2eu(x1,x2) #squared, ce qui nest pas indique

philentropy::squared_euclidean(as.numeric(x1),as.numeric(x2)) #1630.034

Rfast::Dist(x12,method = "euclidean",square = T) #1630.034

analogue::distance(x1,x2,method = "SQeuclidean") # 1630.034

qkerntool::Eucdist(as.matrix(x1),as.matrix(x2),sEuclidean = F) #1630.034

laGP::distance(x1,x2) # 1630.034

#### Root mean square/average distance ####
abdiv::rms_distance(as.numeric(x11_pa[1,]),as.numeric(x11_pa[2,])) #6.086561

?EnvNJ::metrics(t(x11),method = "avg") #64.07 surement pas normalisé

#### Bray-Curtis distance ####--------------------------------------------------

ecodist::distance(x11,method = "bray-curtis") #0.5310021
ecodist::bcdist(x11) # 0.5310021
vegan::vegdist(x11,method = "bray") #0.5310021

provenance::bray.diss(x1,x1) #0.5310021

otuSummary::calc_bc(x11) #0.5310021

vegan::vegdist(x11,method = "bray",binary = T) #0.2

analogue::distance(x1,x1,method = "bray") #0.5310021

pctax::mat_dist(t((x12)),method = "bray") #0.5310021

bioregion::dissimilarity(as.matrix(x12),metric = "Bray") #0.5310021

bioregion::dissimilarity(as.matrix(x12),metric = "Brayturn") #0.5293722

fAssets::braycurtisDist(t(x12)) # 0.5310021

vegan::designdist(x12,method = "(A+B-2*J)/(A+B)",terms = "minimum") # 0.5310021

ecodive::bray(x12,rescale = F) #0.5310021

abdiv::bray_curtis(as.numeric(x1),as.numeric(x2)) #0.5310021

tabula::index_bray(as.numeric(x1),as.numeric(x1)) #0.4689979

tabula::similarity(x12,method = "bray") #0.4689979

diverse::dis_entities(t(x12),method = "Bray",category_row = T)[1,2] #marche pas

benthos::bray_curtis(x1,x2) #0.5310021

chemodiv::sampDis(x12,type = "BrayCurtis")$BrayCurtis[1,2] ##0.5310021

wiqid::distBrayCurtis(x1,x2) ##0.5310021

fossil::bray.curtis(x1,x1) #0.4689979

proxy::dist(x12,method = "Bray") #0.5310021

labdsv::dsvdis(x12, index = "bray/curtis") #0.5310021

rnndescent::brute_force_knn(varespec,1,metric = "braycurtis")

PERMANOVA::DistContinuous(x12,coef = 8)$D #0.5310021

ClusterR::distance_matrix(x12,method = "braycurtis") #0.5310021

provenance::bray.diss(x1,x2) #0.5310021

NST::beta.g(x12,dist.method = "bray") #0.5310021

#### Manhattan distance ####----------------------------------------------------
ecodist::distance(x12,method = "manhattan") #95.06
dist(x12,method = "manhattan") #95.06
vegan::vegdist(x12,method = "manhattan") #95.06

proxyC::dist(as.matrix(x12),method = "manhattan") #95.06

mgc::mgc.distance(x12,method = "manhattan") #95.06

cluster::daisy(x12,metric = "manhattan") #95.06

analogue::distance(x1,x2,method = "manhattan") #95.06

ChemoSpecUtils::rowDist(as.matrix(x12),method = "manhattan") #95.06

ldt::s.distance(t(x12),distance = "manhattan") #95.06

hmsr::manhattan_distance(as.numeric(x1),as.numeric(x2)) #95.06

vegan::vegdist(x12,method = "manhattan",binary = T) #11

vegan::designdist(x12,method = "A+B-2*J",terms = "minimum") # 95.06

distantia::distance(as.numeric(x1),as.numeric(x2),method = "hellinger") #marche pas

fAssets::manhattanDist(t(x12)) #95.06

ecodive::manhattan(x12,rescale = F) #95.06

print(adespatial::dist.ldc(x12,method = "manhattan")) #95.06

abdiv::manhattan(as.numeric(x1),as.numeric(x2)) #95.06

diverse::dis_entities(t(x12),method = "Manhattan",category_row = T)[1,2] #95.06

pctax::mat_dist(t((x12)),method = "manhattan") #95.06

fossil::manhattan(x1,x2) #95.06

proxy::dist(x12,method = "Manhattan") #95.06

amap::Dist(x12,method = "manhattan") #95.06

Mercator::binaryDistance(t(x12),metric = "manhattan") #1

dynutils::calculate_distance(x12,method = "manhattan") #95.06

philentropy::manhattan(as.numeric(x1),as.numeric(x2)) #95.06

EnvNJ::metrics(t(x12),method = "manhattan") #95.06

Rfast::Dist(x12,method = "manhattan") #95.06

comparator::Manhattan()(as.numeric(x1),as.numeric(x2)) #95.06

ClusterR::distance_matrix(x12,method = "manhattan") #95.06

rdist::rdist(x12,metric = "manhattan") #95.06

BoutrosLab.plotting.general::dist(x12,method = "manhattan") #95.06

TSdist::LPDistance(as.numeric(x1),as.numeric(x2),method = "manhattan") #95.06

coda.base::dist(x12,"manhattan") #95.06

LearnClust::mdistance(as.numeric(x1),as.numeric(x2)) #11.08

NST::beta.g(x12,dist.method = "manhattan") #95.06

Rlof::distmc(x12,method = "manhattan") #95.06

?NST::beta.g(x12,dist.method = "mManhattan") #2.880606 modif par Andersen et al. 2006

fda.usc::metric.dist(x12,method = "manhattan") #95.06

flexclust::dist2(x1,x2,method = "manhattan") #95.06

#### Mahalanobis distance ####--------------------------------------------------
ecodist::distance(x12,method = "mahalanobis") # (squared Mahalanobis distance) # Matrice singuliere ne fonctionne pas
vegan::vegdist(x12,method = "mahalanobis") #1.414214

vegan::vegdist(x12,method = "mahalanobis",binary = T) #5.567764

ade4::dist.quant(x12,method = "3") #2

distances::distances(x12,normalize = "mahalanobize") #marche pas

assertr::maha_dist(x12) #sortie bizarre 0.5

nipnTK::outliersMD(x1,x2) #NA

diverse::dis_entities(t(x12),method = "Mahalanobis",category_row = T)[1,2] # marche pas

FD::mahaldis(as.matrix(x12)) #1.414214

pctax::mat_dist(t((x12)),method = "mahalanobis") #1.414214

StatMatch::mahalanobis.dist(x1,x2)#marche pas

asbio::D.sq(t(x1),t(x2)) # ???

proxy::dist(x12,method = "Mahalanobis") #marche pas

MultivariateAnalysis::Distancia(x12,Metodo = 7)[1] #marche pas

ClusterR::distance_matrix(x12,method = "mahalanobis") #1.414214

NST::beta.g(x12,dist.method = "mahalanobis")#1.414214

fAssets::mahalanobisDist(t(x12)) # marche pas

#### Jaccard ####---------------------------------------------------------------
ecodist::distance(x12,method = "jaccard") # 0.3333333
dist(x12,method = "binary") # 0.3333333 
vegan::vegdist(x12,method = "jaccard",binary = F) # 0.6936661
vegan::vegdist(x12,method = "jaccard",binary = T) # 0.3333333

neighbr::similarity(x12_pa[1,],x12_pa[2,],measure = "jaccard") #0.6666667

ChemoSpecUtils::rowDist(as.matrix(x12),method = "binary") #0.33333

mgc::mgc.distance(x12,method = "binary") #0.33333

jaccard::jaccard(x1int,x2int) #fonctione pas

statisfactory::fuzzyJaccard(x1,x2) #0.3063339

arulesSequences::similarity(as.numeric(x1),as.numeric(x2),method = "jaccard") # marche pas

SemNeT::similarity(t(x12),method = "jaccard") # marche pas

ConNEcT::funClassJacc(as.numeric(x1),as.numeric(x2)) #marche pas
ConNEcT::funCorrJacc(as.numeric(x1),as.numeric(x2)) #0.9834216

picante::species.dist(t(x12),metric = "jaccard") #0.3063339

iTOP::jaccard(x12_pa[1,],x12_pa[2,]) #0.666667

Signac::Jaccard(x1,x2) #marche pas

geocmeans::calc_jaccard_idx(as.numeric(x1),as.numeric(x2)) #0.3063339

OTclust::jaccard(x1,x2) #marche pas

?BoutrosLab.plotting.general::dist(x12,method = "jaccard") # 0.6936661

BoutrosLab.plotting.general::dist(x12,method = "binary") # 0.33

bioregion::dissimilarity(as.matrix(x12),metric = "Jaccard") # 0.33

fAssets::jaccardDist(t(x12)) #0.33

?pctax::mat_dist(t((x12)),method = "jaccard") # # 0.6936661

vegan::designdist(x12,method = "(A+B-2*J)/(A+B-J)",terms = "binary") #0.3333333

1- vegan::betadiver(x12,"j") #0.3333333 # conversion en dissim

vegan::nestedbetajac(x12)[3] #0.333333

proxyC::simil(as.matrix(x1),as.matrix(x2),method = "jaccard") #0.6666
proxyC::simil(as.matrix(x1),as.matrix(x2),method = "ejaccard") #0.3430744
proxyC::simil(as.matrix(x1),as.matrix(x2),method = "fjaccard") #NA

ecodive::jaccard(x12)#0.333333

betapart::beta.pair(x12_pa,index.family = "jaccard")$beta.jac #0.3333333

?adespatial::beta.div(x12,method = "ab.jaccard",save.D = T,sqrt.D = F) #0.05012367

adespatial::beta.div.comp(x12,coef = "J",quant = F)$D #0.3333333
adespatial::beta.div.comp(x12,coef = "J",quant = T)$D #0.6936661
adespatial::beta.div.comp(x12,coef = "BJ")$D #0.3333333
?adespatial::beta.div.comp(x12,coef = "BJ",quant = T) #0.6936661

print(adespatial::dist.ldc(x12,method = "jaccard")) #0.5773503
print(adespatial::dist.ldc(x12,method = "ab.jaccard")) #0.05012367

ade4::dist.binary(x11,method = 1) #0.5773503

abdiv::jaccard(as.numeric(x1),as.numeric(x2)) #0.3333333

tabula::index_jaccard(as.numeric(x1),as.numeric(x2)) #0.6666667

tabula::similarity(x12,method = "jaccard") #0.6666667

adiv::betastatjac(x12)[1] #0.3333333 

diverse::dis_entities(t(x12),method = "Jaccard",category_row = T)[1,2] #0.3333333

diverse::dis_entities(t(x12),method = "fJaccard",category_row = T)[1,2] #marche pas

BAT::beta(x12,func = "jaccard",abund = F)$Btotal #0.3333333
BAT::beta(x12,func = "jaccard",abund = T)$Btotal #0.6936661

?wiqid::distChaoJaccCorr(x1,x1) #0.05012367
wiqid::distChaoJaccNaive(x1,x2) #0.05148794

wiqid::distJaccard(x1,x2) #0.3333333

fossil::jaccard(x1,x2) #0.6666667

proxy::dist(x12,method = "Jaccard") #0.3333333

proxy::dist(x12,method = "fJaccard") # marche pas

CommEcol::dis.chao(x12int,index = "jaccard",version = "rare") #0.05012367
CommEcol::dis.chao(x11,index = "jaccard",version = "probability") #0.05148794
vegan::chaodist(x12int,method = "1 - U*V/(U+V-U*V)") # 0.2098765 

adiv::dsimcom(x12,method = "2",type = "similarity",option = "absolute") #0.6564616
diverse::dis_entities(t(x12),method = "eJaccard",category_row = T) #0.6569256
?proxy::dist(x12,method = "eJaccard") #0.6569256

MultivariateAnalysis::Distancia(x12,Metodo = 12)[1] #NA

?MultBiplotR::BinaryDistances(as.matrix(x12_pa),coefficient = 3) #0.5773503 

PERMANOVA::DistBinary(x12_pa,coefficient = 3,transformation = 1)$D #0.6666667

Mercator::binaryDistance(t(x12),metric = "jaccard") #2.266299

arules::dissimilarity(as.matrix(x12_pa),method = "dice") #0.2

philentropy::jaccard(as.numeric(x1),as.numeric(x2)) #0.6569256

EnvNJ::metrics(t(x12),method = "jaccard") #0.6569256

ClusterR::distance_matrix(x12,method = "jaccard_coefficient") #0.6666667

rdist::rdist(x12,metric = "jaccard") #0.9393939

spaa::sp.pair(t(as.matrix(x12)))$Jaccard #0.7142857

NST::beta.g(x12,dist.method = "jaccard") #0.6936661

?NST::beta.g(x12,dist.method = "chao.jaccard") #0.0478141 Chao-Jaccard

prabclus::jaccard(t(x12_pa)) #0.3333333

flexclust::dist2(x1,x2,method = "binary")#0.3333333

# Sans comprehension du calcul :
adespatial::beta.div(x12,method = "jaccard",save.D = T,sqrt.D = F,)$D #0.7637626 # ne correspond pas



#### Turnover component of Jaccard dissimilarity ####---------------------------
vegan::nestedbetajac(x12)[1] #0.26666667 

betapart::beta.pair(x12_pa,index.family = "jaccard")$beta.jtu #0.2666667

adespatial::beta.div.comp(x12,coef = "BJ")$repl #0.2666667 baselga
adespatial::beta.div.comp(x12,coef = "BJ",quant = T)$repl # baselga 0.6922739

bioregion::dissimilarity(as.matrix(x12),metric = "Jaccardturn") # 0.2666667

abdiv::jaccard_turnover(as.numeric(x1),as.numeric(x2)) #0.2666667

adiv::betastatjac(x12)[2] # Ricotta & Pavoine 2015 0.1212121 

#### Replacement index of Jaccard, Podani family ####
?adespatial::beta.div.comp(x12,coef = "J")$repl #0.2424242
adespatial::beta.div.comp(x12,coef = "J",quant = T)$repl #0.6891419

?BAT::beta(x12,func = "jaccard",abund = F)$Brepl #0.2424242
BAT::beta(x12,func = "jaccard",abund = T)$Brepl # 0.6891419


#### Nestedness component of Jaccard dissimilarity ####-------------------------
vegan::nestedbetajac(x12)[2] # 0.06666667

betapart::beta.pair(x12_pa,index.family = "jaccard")$beta.jne #0.06666667

adespatial::beta.div.comp(x12,coef = "J")$rich # Podani 0.09090909
adespatial::beta.div.comp(x12_pa,coef = "J",quant=T)$rich # 0.004524227 Podani 
adespatial::beta.div.comp(x12,coef = "BJ")$rich #0.06666667 baselga
adespatial::beta.div.comp(x12,coef = "BJ",quant = T)$rich #0.001392223 baselga
abdiv::jaccard_nestedness(as.numeric(x1),as.numeric(x2)) #0.06666667

BAT::beta(x12,func = "jaccard",abund = F)$Brich #0.09090909 podani
BAT::beta(test,func = "jaccard",abund = F)$Bgain #0.1212121
BAT::beta(x12,func = "jaccard",abund = F)$Bloss #0.2121212

BAT::beta(x12_pa,func = "jaccard",abund = T)$Brich #0.004524227 podani
BAT::beta(x12_pa,func = "jaccard",abund = T)$Bgain #0.3490952
BAT::beta(x12_pa,func = "jaccard",abund = T)$Bloss #0.3445709

?adiv::betastatjac(x12)[3]# Ricotta & Pavoine 2015 0.2121212
#### Sorensen ####--------------------------------------------------------------
vegan::betadiver(x11,"sor") # 0.8 # COnversion en dissim

ecodist::distance(x11,method = "sorensen") # 0.2

vegan::designdist(x11,method = "(A+B-2*J)/(A+B)",terms = "binary") #0.2

vegan::chaodist(x12int,method = "1 - 2*U*V/(U+V)") # 0.1172414 faux

vegan::nestedbetasor(x11)[3] #0.2

bioregion::dissimilarity(as.matrix(x11),metric = "Sorensen") #0.2

ecodive::sorensen(x11) #0.2

forestmangr::similarity_matrix(x12,index = "Sorensen") # non fonctionnel

divo::li(t(x12)) #sorties pas adequate

superbiclust::sorensenMat(x1,x2) #non fonctionnel

betapart::beta.pair(x11_pa,index.family = "sorensen")$beta.sor #0.2

adespatial::beta.div(x12,method = "ab.sorensen",save.D = T,sqrt.D = F)$D #0.02570608

adespatial::beta.div.comp(x12,coef = "S",quant = F)$D #0.2

adespatial::beta.div.comp(x12,coef = "S",quant = T)$D #0.5310021

adespatial::beta.div.comp(x11,coef = "BS",quant = F)$D #0.2
adespatial::beta.div.comp(x12,coef = "BS",quant = T)$D #0.5310021

print(adespatial::dist.ldc(x12,method = "sorensen")) #0.4472136
print(adespatial::dist.ldc(x12,method = "ab.sorensen")) #0.02570608

ade4::dist.binary(x12,method = 5) #0.4472136 Sorensen Dice

abdiv::sorenson(as.numeric(x1),as.numeric(x1)) #0.2

tabula::index_sorensen(as.numeric(x1),as.numeric(x1)) #0.8

tabula::similarity(x11,method = "sorensen") #0.8

adiv::betastatsor(x11)[1] #0.2
adiv::dsimcom(x12,method = "3",type = "dissimilarity",option = "absolute") #0.4891208 

diverse::dis_entities(t(x12),method = "Dice",category_row = T)[1,2] # marche pas

diverse::dis_entities(t(x12),method = "eDice",category_row = T)[1,2] #0.4891208

proxyC::simil(as.matrix(x1),as.matrix(x1),method = "dice")#0.8
proxyC::simil(as.matrix(x1),as.matrix(x2),method = "edice")  #0.5108792

BAT::beta(x11,func = "sorensen",abund = F)$Btotal #0.2
BAT::beta(x11,func = "sorensen",abund = T)$Btotal #0.5310021

wiqid::distChaoSorCorr(x1,x2) #0.02570608
wiqid::distChaoSorNaive(x1,x2) #0.02642424

wiqid::distSorensen(x1,x1) #0.2

fossil::sorenson(x1,x1) #0.8

proxy::dist(x11,method = "Dice") #0.2
proxy::dist(x11,method = "eDice") #0.4891208

labdsv::dsvdis(x11, index = "sorensen") #0.2

CommEcol::dis.chao(x11,index = "sorensen",version = "rare") #0.02570608
CommEcol::dis.chao(x12,index = "sorensen",version = "probability") #0.02642424

adespatial::beta.div(x12,method = "sorensen",save.D = T,sqrt.D = F)$D #0.6416889 ne correspond pas

MultivariateAnalysis::Distancia(x12,Metodo = 13)[1] #NA

MultBiplotR::BinaryDistances(as.matrix(x12_pa),coefficient = 7) #0.4472136

PERMANOVA::DistBinary(x11_pa,coefficient = 7,transformation = 1)$D #0.8
PERMANOVA::DistBinary(x12_pa,coefficient = 16,transformation = 1)$D #0.8
PERMANOVA::DistBinary(x12_pa,coefficient = 17,transformation = 1)$D #0.8

philentropy::dice_dist(as.numeric(x1),as.numeric(x1)) #0.4891208

philentropy::sorensen(as.numeric(x1),as.numeric(x1)) #0.5310021

EnvNJ::metrics(t(x11),method = "sorensen") #0.5310021

Rfast::Dist(x12,method = "sorensen") #NA

spaa::sp.pair(t(as.matrix(x12)))$Dice #0.833333

NST::beta.g(x12,dist.method = "chao.sorensen") #0.0244926

prabclus::dicedist(t(x12)) #0.4375 ???

fAssets::sorensenDist(t(x12)) #0.2

#### Gower ##dice_dist()#### Gower ####-----------------------------------------------------------------

?ecodist::distance(x12,method = "gower") # 95.06 # correspond a rien

ecodist::distance(x12,method = "modgower10") # (modified Gower, base 10) # 23.03452

ecodist::distance(x12,method = "modgower2") # (modified Gower, base 2) # 76.51902

gower::gower_dist(x1,x2) #1

cluster::daisy(x12,metric = "gower") #0.7045455

FD::gowdis(as.matrix(x12)) #NA

StatMatch::gower.dist(x1,x2)#0.7045455

analogue::distance(x1,x2,method = "gower") #marche pas
analogue::distance(x1,x2,method = "alt.gower") #marche pas

vegan::vegdist(x12,method = "gower") #0.7045455
#vegan::vegdist(x12,method = "gower",binary = T) #0.7045455

pctax::mat_dist(t((x12)),method = "gower") #0.7045455

vegan::vegdist(x12,method = "altGower") #2.880606

FD::gowdis(as.numeric(x12),w = NULL) #marche pas

pctax::mat_dist(t((x12)),method = "altGower") #2.880606
# There are two versions of Gower distance in vegan ("gower", "altGower") which differ in scaling: "gower" divides all distances by the number of observations (rows) and scales each column to unit range, but "altGower" omits double-zeros and divides by the number of pairs with at least one above-zero value, and does not scale columns
# pctax reprend la fonction vegan donc la distinction aussi

ecodive::gower(x12,rescale = F) #0.7045455

diverse::dis_entities(t(x12),method = "Gower",category_row = T)[1,2] #0.7045455

BAT::gower(x12) # pas adapte

shipunov::Gower.dist(x1,x2) #0.7045455

proxy::dist(x12,method = "Gower") #0.7045455

CommEcol::dis.goodall(x12,p.simi = "gower",approach = "proportion") #marche pas
CommEcol::dis.goodall(x12,p.simi = "gower",approach = "chisquare") # marche pas

MultivariateAnalysis::Distancia(x12,Metodo = 21)[1] #NA
MultivariateAnalysis::Distancia(x12,Metodo = 22)[1] #NA

philentropy::gower(as.numeric(x1),as.numeric(x2)) #2.160455

Rfast::Dist(x12,method = "gower") #2.160455

shipunov::Gower.dist(x12) #0.70455455

NST::beta.g(x12,dist.method = "gower") #0.7045455

NST::beta.g(x12,dist.method = "altGower") #2.880606

NST::beta.g(x12,dist.method = "mGower") #0.5731683 modif par Andersen et al. 2006
#### Canberra ####--------------------------------------------------------------
dist(x11,method = "canberra") #27.50188
vegan::vegdist(x12,method = "canberra") #0.6250428  
# Canberra index is divided by the number of variables in vegdist, but not in dist. So these differ by a constant multiplier, and the alternative in vegdist is in range (0,1).
vegan::vegdist(x11,method = "canberra",binary = T) #0.3333333

pctax::mat_dist(t((x12)),method = "canberra") #0.6250428

ecodive::canberra(x12,rescale = F) #20.62641

mgc::mgc.distance(x12,method = "canberra") #27.50188

LearnClust::canberradistance(as.numeric(x1),as.numeric(x2)) #1.068272

ChemoSpecUtils::rowDist(as.matrix(x12),method = "canberra") #27.50188

adespatial::beta.div(x12,method = "canberra",save.D = T)$D #0.6250428

fAssets::canberraDist(t(x12)) #27.50188

abdiv::canberra(as.numeric(x1),as.numeric(x1)) #20.62641

diverse::dis_entities(t(x12),method = "Canberra",category_row = T)[1,2] #27.50188

proxy::dist(x12,method = "Canberra") #27.50188

proxyC::dist(as.matrix(x12),method = "canberra") ##20.62641

BoutrosLab.plotting.general::dist(x12,method = "canberra") #27.50188

phm::canberra(as.numeric(x1),as.numeric(x2)) #20.62641

amap::Dist(x12,method = "canberra") #27.50188

PERMANOVA::DistContinuous(x12,coef = 7)$D #NA

Mercator::binaryDistance(t(x12),metric = "canberra") #13.75094

dynutils::calculate_distance(x12,method = "canberra") #20.62641

philentropy::canberra(as.numeric(x1),as.numeric(x2)) #20.62641

EnvNJ::metrics(t(x12),method = "canberra") #27.50188

Rfast::Dist(x12,method = "canberra") #NA

Rlof::distmc(x12,method = "canberra") #27.50188

ClusterR::distance_matrix(x12,method = "canberra") #31.62625

rdist::rdist(x12,metric = "canberra") #NA

coda.base::dist(x12,"canberra") #27.50188

NST::beta.g(x12,dist.method = "canberra") #0.6250428 correspond pas

fda.usc::metric.dist(x12,method = "canberra") #27.50188

flexclust::dist2(x1,x2,method = "canberra") #27.50188
#### Minkowski ####-------------------------------------------------------------
dist(x12,method = "minkowski",p=2) #40.37368

flexclust::dist2(x1,x2,method = "minkowski") #40.37368

mgc::mgc.distance(x12,method = "minkowski") #40.37368

ecodive::minkowski(x12,rescale = F,power = 2) # 40.37368

abdiv::minkowski(as.numeric(x1),as.numeric(x2),p=2) #p modifiable -> p=1 = Manhattan 40.37368

diverse::dis_entities(t(x12),method = "Minkowski",category_row = T)[1,2] # marhe pas

proxy::dist(x12,method = "Minkowski",p=2) #p modifiable -> p=1 = Manhattan 40.37368

PERMANOVA::DistContinuous(x12,coef = 4,r=2)$D # 40.37368

dynutils::calculate_distance(x12,method = "minkowski") #40.37368

philentropy::minkowski(as.numeric(x1),as.numeric(x2),n=2) #40.37368

EnvNJ::metrics(t(x12),method = "minkowski") #40.37368

Rfast::Dist(x12,method = "minkowski",p = 2) #95.06

fAssets::minkowskiDist(t(x12)) #40.37368

proxyC::dist(as.matrix(x12),method = "minkowski",p=2) #40.37368

comparator::Minkowski(p=2)(as.numeric(x1),as.numeric(x2)) #40.37368

ClusterR::distance_matrix(x12,method = "minkowski",minkowski_p = 2) #40.37368

rdist::rdist(x12,metric = "minkowski",p = 2) #40.37368

coda.base::dist(x12,"minkowski",p=2) #40.37368

TSdist::MinkowskiDistance(as.numeric(x1),as.numeric(x2),p=2) #40.37368

TSdist::LPDistance(as.numeric(x1),as.numeric(x2),p=2,method = "minkowski") #40.37368

Rlof::distmc(x12,method = "minkowski",p=2) #40.37368

fda.usc::metric.dist(x12,method = "minkowski",p=2) #40.37368

BoutrosLab.plotting.general::dist(x12,method = "minkowski",p=2) #40.37368
#### Clark ####-----------------------------------------------------------------
vegan::vegdist(x12,method = "clark") #0.7262691
?vegan::vegdist(x12,method = "clark",binary = T) #0.5773503

pctax::mat_dist(t((x12)),method = "clark") #0.7262691

ecodive::clark(x12,rescale = F) # 4.172098 idem sans division par n Clark's divergence distance

abdiv::clark_coefficient_of_divergence(as.numeric(x1),as.numeric(x2)) #0.7262691

?philentropy::clark_sq(as.numeric(x1),as.numeric(x2)) #4.17098
#### Kulczynski ####-----------------------------------------------------------------
vegan::vegdist(x12,method = "kulczynski") #0.5309965
vegan::vegdist(x12,method = "kulczynski",binary = T) #0.1976127

pctax::mat_dist(t((x12)),method = "kulczynski") #0.5309965

adespatial::beta.div(x12,method = "kulczynski",save.D = T)$D #0.5309965

print(adespatial::dist.ldc(x12,method = "kulczynski")) #0.5309965

abdiv::kulczynski_first(as.numeric(x1),as.numeric(x2)) # Kulczynski -1

?abdiv::kulczynski_second(as.numeric(x1),as.numeric(x2)) # Kulczynski-Cody = vegdist avec kulczynski et binary = T 0.1976127

abdiv::weighted_kulczynski_second(as.numeric(x1),as.numeric(x2)) #Weighted Kulczynski distance 0.5309965

diverse::dis_entities(t(x12),method = "Kulczynski1",category_row = T)[1,2] # marche pas
diverse::dis_entities(t(x12),method = "Kulczynski2",category_row = T)[1,2] # marche pas

fossil::kulczynski(x1,x2) #0.8023873

proxy::dist(x12,method = "Kulczynski1") #-1
proxy::dist(x12,method = "Kulczynski2") #0.1976127

vegan::betadiver(x12,"co") #0.1976127
abdiv::kulczynski_second(as.numeric(x1),as.numeric(x2)) # Kulczynski-Cody 0.1976127

MultBiplotR::BinaryDistances(as.matrix(x12_pa),coefficient = 1) # marche pas
MultBiplotR::BinaryDistances(as.matrix(x12_pa),coefficient = 10) # marche pas

PERMANOVA::DistBinary(x12_pa,coefficient = 1,transformation = 1)$D #marche pas
PERMANOVA::DistBinary(x12_pa,coefficient = 10,transformation = 1)$D #marche pas

philentropy::kulczynski_d(as.numeric(x1),as.numeric(x2)) #2.2644

EnvNJ::metrics(t(x12),method = "kulczynski") #2.2644

Rfast::Dist(x12,method = "kulczynski") #2.2644

NST::beta.g(x12,dist.method = "kulczynski") #0.5309965

prabclus::kulczynski(t(x12_pa)) #0.1976127

prabclus::qkulczynski(t(x12)) #0.5309965

#### Morisita ####--------------------------------------------------------------

vegan::vegdist(x12int,method = "morisita") # 0.4806241
vegan::vegdist(x12int,method = "morisita",binary = T) #0

pctax::mat_dist(t((x12int)),method = "morisita") #0.4806241

ecodive::morisita(x12int) #0.4806241

abdiv::morisita(as.numeric(x1int),as.numeric(x2int)) #0.4806241

NST::beta.g(x12int,dist.method = "morisita") #0.4806241

#### Morisita-Horn ####---------------------------------------------------------

vegan::vegdist(x12,method = "horn") #0.4886065
NST::beta.g(x12,dist.method = "horn") #0.4886065

?vegan::vegdist(x12,method = "horn",binary = T) #0.2

pctax::mat_dist(t((x12)),method = "horn") #0.4886065

ecodive::horn(x12,rescale = F) #0.4886065

?abdiv::horn_morisita(as.numeric(x1),as.numeric(x2)) #0.4886065

1-tabula::index_morisita(as.numeric(x1),as.numeric(x2)) #0.5113935

tabula::similarity(x12,method = "morisita") #0.5113935

wiqid::distMorisitaHorn(x1,x2) #0.4886065

fossil::morisita.horn(x1,x2) #0.5113935


#### Mountford ####-------------------------------------------------------------

vegan::vegdist(x12,method = "mountford") #  the proper index is defined as the root of the equation above 0.5666415
vegan::vegdist(x12,method = "mountford",binary = T) #0.5666415

pctax::mat_dist(t((x12)),method = "mountford") #0.5666415
diverse::dis_entities(t(x12),method = "Mountford",category_row = T)[1,2]  #marche pas

?proxy::dist(x12,method = "Mountford") #0.852349 # ne correspond pas

NST::beta.g(x12,dist.method = "mountford") #0.5666415

#### Raup ####-----------------------------------------------------------------

?vegan::vegdist(x12,method = "raup")# 0.002291751 #The index uses equal occurrence probabilities for all species, but Raup and Crick originally suggested that sampling probabilities should be proportional to species frequencies 
vegan::vegdist(x12,method = "raup",binary = T)# 0.002291751 #The index uses equal occurrence probabilities for all species, but Raup and Crick originally suggested that sampling probabilities should be proportional to species frequencies 

vegan::designdist(x12,method = "1-phyper(J-1, A, P-A, B)",terms = "binary") #0.002291751
vegan::raupcrick(x12,"r0") #0.003

pctax::mat_dist(t((x12)),method = "raup") #0.002291751^

iCAMP::RC.pc(x12)$index #louche

NST::beta.g(x12,dist.method = "raup") #0.002291751



#### Binomial ####--------------------------------------------------------------

vegan::vegdist(x11,method = "binomial") # 11.30918
vegan::vegdist(x11,method = "binomial",binary = T) #7.624619
pctax::mat_dist(t((x11)),method = "binomial") # 11.30918

abdiv::binomial_deviance(as.numeric(x1),as.numeric(x1)) # 11.30918

coda.base::dist(x12,"binary") #0.333333 ne sait pas a quoi correspond

NST::beta.g(x11,dist.method = "binomial") #11.30918

#### Binomial co-occurence assessment ####--------------------------------------
?tabula::index_binomial(as.numeric(x1),as.numeric(x1)) #1.503646

#### Chao ####------------------------------------------------------------------

vegan::vegdist(x11,method = "chao") #0.05012367

?pctax::mat_dist(t((x11)),method = "chao") #0.05012367

NST::beta.g(x11,dist.method = "chao") #0.05012367


#### Cao ####-------------------------------------------------------------------

vegan::vegdist(x11,method = "cao") # Cao et al. (1997) used log 10 but the current function uses natural logarithms 0.5159614
vegan::vegdist(x11,method = "cao",binary = T) #0.4984741
pctax::mat_dist(t((x11)),method = "cao") #0.5159614
abdiv::cy_dissimilarity(as.numeric(x1),as.numeric(x1),base = 10) #0.2240792 # Par defaut utilise le log 10 mais peut etre modifie 0.2240792
NST::beta.g(x11,dist.method = "cao")#0.5159614

#### Chi2 ####------------------------------------------------------------------

vegan::vegdist(x12,method = "chisq") #1.197533
vegan::vegdist(x12,method = "chisq",binary = T) #3.316625
pctax::mat_dist(t((x12)),method = "chisq") #1.197533

svs::dist_chisquare(as.matrix(x12)) #1.382792

analogue::distance(x1,x2,method = "chi.square") #8.011437

analogue::distance(x1,x2,method = "chi.distance") #NA

adespatial::beta.div(x12,method = "chisquare",save.D = T)$D #1.197533

SNFtool::chiDist2(x12) #NA

GDAtools::dist.chi2(x12) #marche pas

colordistance::chisqDistance(x1,x2) # fonction n'existe plus

print(adespatial::dist.ldc(x12,method = "chisquare")) #1.197533

proxy::dist(x12,method = "Chi-squared") #-612.7185

labdsv::dsvdis(x12, index = "chisq") #Na

Rfast::Dist(x12,method = "chi_square") #Na

spaa::sp.pair(t(as.matrix(x12)))$chisq #9.573626

proxyC::dist(as.matrix(x12),method = "chisquared") #NA

QFASA::chisq.dist(x1,x2) # marche pas

#### Squared chi2 distance ####
ecodive::squared_chisq(x12,rescale = F) #64.18313

diverse::dis_entities(t(x12),method = "Chi-squared",category_row = T)[1,2] # marche pas

dynutils::calculate_distance(x12,method = "chisquared") #0

analogue::distance(x1,x2,method = "SQchi.square") #64.18313

philentropy::squared_chi_sq(as.numeric(x1),as.numeric(x2)) #64.18313

EnvNJ::metrics(t(x12),method = "squared_chi") #64.18313

#### Probabilistic symmetric chi2 distance ####
?ecodive::psym_chisq(x12,rescale = F) #128.3663

#### Soergel distance ####
ecodive::soergel(x12,rescale = F) #0.6936661

diverse::dis_entities(t(x12),method = "Soergel",category_row = T)[1,2] # marche pas

proxy::dist(x12,method = "Soergel") #2.518146

PERMANOVA::DistContinuous(x12,coef = 9)$D #0.6936661

philentropy::soergel(as.numeric(x1),as.numeric(x2)) #0.6936661

EnvNJ::metrics(t(x12),method = "soergel") #0.6936661

Rfast::Dist(x12,method = "soergel") #0.6936661

#### Chord ####-----------------------------------------------------------------

vegan::vegdist(x12,method = "chord") #0.9832398
vegan::vegdist(x12,method = "chord",binary = T) #0.6305668

pctax::mat_dist(t((x12)),method = "chord") #0.9832398

analogue::distance(x1,x2,method = "chord") #6.514023

ecodive::chord(x12) #0.9832398

adespatial::beta.div(x12_pa,method = "chord",save.D = T)$D #0.9832398

print(adespatial::dist.ldc(x12,method = "chord")) #0.9832398

abdiv::chord(as.numeric(x1),as.numeric(x2)) #0.9832398

diverse::dis_entities(t(x12),method = "Chord",category_row = T)[1,2] # marche pas

wiqid::distChord(x1,x2) #0.9832398

proxy::dist(x12,method = "Chord") #0.9832398
#### Squared chord distance ####

ecodive::squared_chord(x12,rescale = F) # 42.4325

analogue::distance(x1,x1,method = "SQchord") #42.4325

philentropy::squared_chord(as.numeric(x1),as.numeric(x2)) #42.4325

EnvNJ::metrics(t(x12),method = "squared_chord") #42.4325

#### Log chord ####
adespatial::beta.div(x12_pa,method = "log.chord",save.D = T)$D #0.6607558

print(adespatial::dist.ldc(x12,method = "log.chord")) #0.6607558

#### Hellinger ####-------------------------------------------------------------

vegan::vegdist(x12,method = "hellinger") #0.6885085
vegan::vegdist(x12,method = "hellinger",binary = T) #0.6305668

pctax::mat_dist(t((x12)),method = "hellinger") #0.6885085

topicmodels::distHellinger(x1,x2) # marche pas

?ecodive::hellinger(x12,rescale = T) # 0.6885085

statip::hellinger(as.numeric(x1),as.numeric(x2)) #non fonctionnel

textmineR::CalcHellingerDist(as.numeric(x1),as.numeric(x2)) #0.48962385

dad::hellinger(x1,x2) #non fonctionnel NA

adespatial::beta.div(x12,method = "hellinger",save.D = T)$D #0.6885085

print(adespatial::dist.ldc(x12,method = "hellinger")) #0.6885085

abdiv::hellinger(as.numeric(x1),as.numeric(x2)) #0.6885085

diverse::dis_entities(t(x12),method = "Hellinger",category_row = T)[1,2] # ne marche pas

proxy::dist(x12,method = "Hellinger") #0.6885085

philentropy::hellinger(as.numeric(x1),as.numeric(x2)) #NA

Rfast::Dist(x12,method = "hellinger") #1152.608
Rfast::Dist(x12,method = "hellinger",square = T) #20.18684


#### Aitchison ####-------------------------------------------------------------

vegan::vegdist(x12,method = "aitchison") # marche pas
pctax::mat_dist(t((x12)),method = "aitchison") # marche pas

?vegan::vegdist(x12,method = "robust.aitchison") #10.3258
vegan::vegdist(x11,method = "robust.aitchison",binary = T) #2.645751
pctax::mat_dist(t((x12)),method = "robust.aitchison") #10.3258

?ecodive::aitchison(x12) # 10.68875

coda.base::dist(x12,"aitchison")#NA

robCompositions::aDist(x1_pos,x2_pos) #NA

QFASA::AIT.dist(x1_pos,x2_pos) #NA

#### Ochiai ####----------------------------------------------------------------

vegan::chaodist(x12int,method = "1 - sqrt(U*V)") # 0.1111111 ok

adespatial::beta.div(x12int,method = "ochiai",save.D = T,sqrt.D = F)$D #0.6408934

adespatial::beta.div(adespatial::beta.div(x12int,method = "ab.ochiai",save.D = T,sqrt.D = F)$D) #0.02555303 ok avec integer

print(adespatial::dist.ldc(x12int,method = "ochiai")) #0.4458781 ok avec présence absence
print(adespatial::dist.ldc(x12int,method = "ab.ochiai")) #0.02555303 ok avec integer

ade4::dist.binary(x12int,method = 7) #0.4458781 ok avec présence absence

adiv::dsimcom(x12int,method = "4",type = "similarity",option = "absolute") #0.4833802

diverse::dis_entities(t(x12int),method = "Ochiai",category_row = T)[1,2] #marche pas

MultivariateAnalysis::Distancia(x12int,Metodo = 17)[1] #NA
MultivariateAnalysis::Distancia(x12int,Metodo = 18)[1] #NA

MultBiplotR::BinaryDistances(as.matrix(x12int),coefficient = 12) #0.4458781 avec présence absence

spaa::sp.pair(t(as.matrix(x12int)))$Ochiai #0.8451543


#### Otsuka-Ochiai ####---------------------------------------------------------
?ecodive::ochiai(x12) #0.1988073

wiqid::distOchiai(x1,x2) #0.1988073

fossil::ochiai(x1,x2) #0.8011927
PERMANOVA::DistBinary(x12_pa,coefficient = 12,transformation = 1)$D #0.8011927
proxy::dist(x12,method = "Ochiai") #0.1988073

labdsv::dsvdis(x12, index = "ochiai") #0.1988073

vegan::designdist(x12,method = "1-J/sqrt(A*B)",terms = "binary") # 0.1988073

#### Ruzicka ####----------------------------------------------------------------
vegan::designdist(x12,method = "(A+B-2*J)/(A+B-J)",terms = "minimum") # 0.6936661

adespatial::beta.div(x12,method = "ruzicka",save.D = T)$D #0.6936661

print(adespatial::dist.ldc(x12,method = "ruzicka")) #0.6936661

abdiv::ruzicka(as.numeric(x1),as.numeric(x2)) #0.6936661

labdsv::dsvdis(x12, index = "ruzicka") #0.6936661

philentropy::ruzicka(as.numeric(x1),as.numeric(x1)) #0.3063339

stylo::dist.minmax(x12) #0.6936661

#### Dissimilarity ratio ####---------------------------------------------------
vegan::designdist(x12,method = "(A+B-2*J)/(A+B-J)",terms = "quadratic") # 0.6569256

#### Cosine complement ####
vegan::designdist(x12,method = "1-J/sqrt(A*B)",terms = "quadratic") #0.4833802

dynutils::calculate_distance(x12,method = "cosine") #0.3456198

SemNeT::similarity(t(x12),method = "cosine") #0.5166198 #ok

EnvNJ::metrics(t(x12),method = "cosine") #0.2764125

?ChemoSpecUtils::rowDist(as.matrix(x12),method = "cosine") #0.4833802

Rfast::Dist(x12,method = "cosine") #0.9667605

?ClusterR::distance_matrix(x12,method = "cosine") #0.4833802

?proxyC::simil(as.matrix(x1),as.matrix(x2),method = "cosine") #0.5166198

svs::dist_cosine(as.matrix(x12)) #0.4833802

?abdiv::cosine_distance(x1,x1) #0.4833802

resemble::f_diss(as.matrix(x12),diss_method = "cosine") #marche pas

#### Koleff beta-3 Williams index ####------------------------------------------
vegan::betadiver(x12,method = "-3") # 0.1212121

#### Koleff betaC Cody ####-----------------------------------------------------
vegan::betadiver(x12,"c") # 5.5

tabula::index_cody(as.matrix(x12)) #11

#### Koleff betasor Simpson dissimilarity turnover sorensen ####----------------

fossil::simpson(x1,x2) #0.8461538
vegan::nestedbetasor(x12)[1] #0.1538462 
betapart::beta.pair(x12_pa,index.family = "sorensen")$beta.sim #0.1538462
proxy::dist(x12,method = "Simpson") #0.1538462
adespatial::beta.div.comp(x11,coef = "BS")$repl # Baselga 0.1538462

adespatial::beta.div.comp(x12,coef = "BS",quant = T)$repl # Baselga 0.5293722

adespatial::beta.div(x12,method = "ab.simpson",save.D = T)$D #0.008127366 correspond a rien
print(adespatial::dist.ldc(x12,method = "ab.simpson")) #0.008127366 correspond a rien

adespatial::beta.div.comp(x12,coef = "J",quant = F)$repl # Podani 0.2424242
adespatial::beta.div.comp(x12,coef = "J",quant =T)$repl # Podani 0.6891419

abdiv::jaccard_turnover(as.numeric(x1),as.numeric(x2)) #0.2666667

adiv::betastatsor(x12)[2] # Ricotta & Pavoine  0.07272727

BAT::beta(x12,func = "sorensen",abund = F)$Brepl #0.1454545 ne correspond pas
BAT::beta(x12,func = "sorensen",abund = T)$Brepl #0.5275388 ne correspond pas

#### Koleff beta-2 Harrison1 ####-----------------------------------------------
vegan::betadiver(x12,"-2") # 0.137931

#### Koleff betaw Whittaker ####------------------------------------------------
vegan::betadiver(x11,"w") # 0.2 ???

?adespatial::beta.div(x11,method = "whittaker",save.D = T)$D #0.5319661

print(adespatial::dist.ldc(x12,method = "whittaker"))#0.5319661

tabula::index_whittaker(as.matrix(x11)) #0.2 ???

?diverse::dis_entities(t(x12),method = "Whittaker",category_row = T)[1,2] # marche pas

?wiqid::distWhittaker(x1,x1) #0.5319661

?proxy::dist(x12,method = "Whittaker") #0.5319661

#### Koleff betar Routledge ####------------------------------------------------
vegan::betadiver(x12,"r") # 0.05421104

tabula::index_routledge1(as.matrix(x12)) # 0.05421104
tabula::index_routledge2(as.matrix(x12)) #NA
tabula::index_routledge3(as.matrix(x12)) #NA

#### Koleff betal Routledge ####------------------------------------------------
vegan::betadiver(x12,"l") #5.5

#### Koleff betae Routledge ####------------------------------------------------
?vegan::betadiver(x12,"e") #0.14699

#### Koleff betacc Colwell & Coddington ####------------------------------------
vegan::betadiver(x12,"cc") #0.33333

#### Koleff betasim Lennon ####-------------------------------------------------
vegan::betadiver(x12,"sim") #0.1538462

#### Koleff betagl Lennon ####--------------------------------------------------
vegan::betadiver(x12,"gl") #0.1090909

#### Koleff betaz Lennon ####---------------------------------------------------
vegan::betadiver(x12,"z") #0.2630344

#### Koleff betat Wilson and Shmida ####----------------------------------------
vegan::betadiver(x12,"t") #0.2

tabula::index_wilson(as.matrix(x12_pa)) #0.4 ne sait pas a quoi ca correspond

#### Koleff betame Mourelle & Ezcurra ####--------------------------------------
vegan::betadiver(x12,"me") #0.2

#### Koleff betawb Wieher & Boylen ####-----------------------------------------
vegan::betadiver(x12,"wb") #11

#### Koleff betam Magurran ####-------------------------------------------------
vegan::betadiver(x12,"m") #18.33333

#### Koleff betag Gaston ####---------------------------------------------------
vegan::betadiver(x12,"g") #0.3333333

#### Koleff betal Lande ####----------------------------------------------------
vegan::betadiver(x12,"l") #5.5

#### Koleff betahk Harte & Kinzig ####------------------------------------------
vegan::betadiver(x12,"hk") #0.2

#### Koleff betarlb Ruggiero ####-----------------------------------------------
1-vegan::betadiver(x12,"rlb") #0.2413793 # Conversion en dissimalirite

#### Koleff beta-1 Harrison 1 ####----------------------------------------------
vegan::betadiver(x12,"-1") #0.2

#### Koleff 19 Sans nom ####----------------------------------------------------
vegan::betadiver(x12,"19") #0.05492424

#### Nestedness-resultant component of Sørensen dissimilarity ####--------------
vegan::nestedbetasor(x12)[2] #0.04615385 

betapart::beta.pair(x12_pa,index.family = "sorensen")$beta.sne #0.04615385

adespatial::beta.div.comp(x12_pa,coef = "S",quant = F)$rich #Podani 0.05454545
adespatial::beta.div.comp(x12_pa,coef = "S",quant = T)$rich #Podani 0.0034633
adespatial::beta.div.comp(x12,coef = "BS")$rich #Baselga 0.04615385
adespatial::beta.div.comp(x12,coef = "BS",quant=T)$rich #Baselga 0.001629925

abdiv::sorenson_nestedness(as.numeric(x1),as.numeric(x2)) #0.04615385 

?adiv::betastatsor(x12)[3] #Ricotta & Pavoine 0.1272727

BAT::beta(x12,func = "sorensen",abund = F)$Brich #0.05454545
BAT::beta(x12,func = "sorensen",abund = T)$Brich #0.0034633

BAT::beta(x12,func = "sorensen",abund = F)$Bgain #0.07272727
?BAT::beta(x12,func = "sorensen",abund = F)$Bloss #0.1272727
BAT::beta(x12_pa,func = "sorensen",abund = T)$Bgain #0.2672327
BAT::beta(x12_pa,func = "sorensen",abund = T)$Bloss #0.2637694

#### Bhattacharyya distance ####
ecodive::bhattacharyya(x11,rescale = F) #-4.223818 
proxy::dist(x11,method = "Bhjattacharyya") #6.514023 le bon car = 0 si identique
philentropy::bhattacharyya(as.numeric(x1),as.numeric(x1)) #-4.223818
EnvNJ::metrics(t(x11),method = "bhattacharyya") #0.2705261
Rfast::Dist(x11,method = "bhattacharyya") #-4.223818
ICGE::dbhatta(x11) #NA


### Autres beta ####------------------------------------------------------------
?ecodist::distance(x12,method = "difference") # 0.12

### Species profile distance ####
adespatial::beta.div(test,method = "profiles",save.D = T)$D #0.4505909

print(adespatial::dist.ldc(x12,method = "profiles")) #0.4505909

#### chebyshev ####-------------------------------------------------------------
ecodive::chebyshev(x11,rescale = F) #33.08

abdiv::chebyshev(as.numeric(x1),as.numeric(x2)) #33.08

philentropy::chebyshev(as.numeric(x1),as.numeric(x2)) #33.08

EnvNJ::metrics(t(x12),method = "chebyshev") #33.08

LearnClust::chebyshevDistance(as.numeric(x1),as.numeric(x2)) #0.12

comparator::Chebyshev()(as.numeric(x1),as.numeric(x2)) #33.08

print(SBCK::chebyshev(as.matrix(x1),as.matrix(x2))) #1

Trading::Chebyshev_distance(x1,x2) #33.08

beadplexr::dist_chebyshev(x12) #33.08

LearnClust::chebyshevDistance(x1,x2) #0.12

ClusterR::distance_matrix(x12,method = "chebyshev") #33.08

rdist::rdist(x12,metric = "chebyshev") #33.08

#### Divergence ####------------------------------------------------------------
?ecodive::divergence(x12,rescale = F) #34.81281

?adespatial::beta.div(x12,method = "divergence",save.D = T)$D #0.7262691 = Clark

print(adespatial::dist.ldc(x12,method = "divergence")) #0.7262691 = Clark

diverse::dis_entities(t(x12),method = "divergence",category_row = T)[1,2] #marche pas

?proxy::dist(x12,method = "divergence") #17.4064

PERMANOVA::DistContinuous(x12,coef = 5)$D #NA

philentropy::divergence_sq(as.numeric(x1),as.numeric(x2)) #34.81281

#### Jenson-Shannon distance ####
ecodive::jensen(x12,rescale = F) #4.330852 distance

priorsense::cjs_dist(as.numeric(x1),as.numeric(x2)) #0.3056077

proxyC::dist(as.matrix(x12),method = "jensen") #NA

philentropy::jensen_shannon(as.numeric(x1),as.numeric(x2),unit = "log10") #18.75628 #jensen shannon divergence

philentropy::jensen_difference(as.numeric(x1),as.numeric(x2)) # 18.75628

EnvNJ::metrics(t(x12),method = "jensen-shannon") #0.2095511 divergence avant racine carré avec les proportions

Rfast::Dist(x12,method = "jensen_shannon") #34.28249

sccore::jsDist(as.matrix(x12)) # marche pas

#### Jenson-Shannon divergence ####
ecodive::jsd(x12,rescale = F) #18.75628 divergence

philentropy::JSD(as.matrix(y12)) #27.05959 divergence log 2

textmineR::CalcJSDivergence(as.numeric(y1),as.numeric(y2)) #0.2094943

philentropy::gJSD(as.numeric(y1),as.numeric(y2)) #marche pas

EnvNJ::metrics(t(y12),method = "jensen_difference") #0.2095511

dlookr::jsd(as.numeric(x1),as.numeric(x2)) # marche pas

Tmisc::jsd(as.matrix(x12),normalizeCounts = T) # marche pas

Compositional::divergence(y12,type = "jensen_shannon") #34.28249

#### Lorentzian distance ####
ecodive::lorentzian(x12,rescale = F) #23.05264

philentropy::lorentzian(as.numeric(x1),as.numeric(x2)) #23.05264

EnvNJ::metrics(t(x12),method = "lorentzian") #23.05264

#### Matusita distance ####
ecodive::matusita(x12,rescale = F) #6.514023

philentropy::matusita(as.numeric(x1),as.numeric(x2)) #NA

Rfast::Dist(x12,method = "jeffries_matusita") #NA

#### Motyka dissimilarity ####
ecodive::motyka(x12,rescale = F) #0.7655011

philentropy::motyka(as.numeric(x1),as.numeric(x2)) #0.7655011

EnvNJ::metrics(t(x12),method = "motyka") #0.7655011

Rfast::Dist(x12,method = "motyka") #0.7655011

#### Jeffreys ####

EnvNJ::metrics(t(x12),method = "jeffreys") #2.341993
#### Topsoe distance ####-------------------------------------------------------
ecodive::topsoe(x11,rescale = F) #37.51256

philentropy::topsoe(as.numeric(x1),as.numeric(x2)) #37.51256

#### Wave Hedges distance ####--------------------------------------------------
ecodive::wave_hedges(x12,rescale = F) #22.84192

diverse::dis_entities(t(x12),method = "Wave",category_row = T)[1,2] # marche pas

proxy::dist(x12,method = "Wave") #1

PERMANOVA::DistContinuous(x12,coef = 10)$D #NA

philentropy::wave_hedges(as.numeric(x1),as.numeric(x1)) #22.84192

EnvNJ::metrics(t(x12),method = "wavehedges") #22.84192

Rfast::Dist(x12,method = "wave_hedges") #NA

#### Hamming distance ####------------------------------------------------------
ecodive::hamming(x11_pa) #11

abdiv::hamming(as.numeric(x1),as.numeric(x2)) #31

CEGO::distanceNumericHamming(x12_pa[1,],x12_pa[2,]) #0.25

Rankcluster::distHamming(x12_pa[1,],x12_pa[2,]) #11

Mercator::binaryDistance(t(x12_pa),metric = "hamming") #11

e1071::hamming.distance(x12_pa) #11

pegas::dist.hamming(x12_pa) #11

tfaddons::metric_hamming_distance(x12_pa[1,],x12_pa[2,]) # marche pas

SID::hammingDist(x12_pa[1,],x12_pa[2,]) # marche pas

EnsCat::hammingD(x12_pa) #0.25

genMCMCDiag::hammingDist(x12_pa[1,],x12_pa[2,]) #11

CEGO::distanceNumericHamming(x12_pa[1,],x12_pa[2,]) #0.25

proxyC::dist(as.matrix(x12_pa),method = "hamming") #11

bingat::calcDistance(x1,x2) #11

EnvNJ::metrics(t(x12_pa),method = "hamming") #0.25 -> 11/44

ClusterR::distance_matrix(x12_pa,method = "hamming") #0.25

rdist::rdist(x12_pa,metric = "hamming") #0.25

#### Wishart ####---------------------------------------------------------------
adespatial::beta.div(x12,method = "wishart",save.D = T)$D #0.6569256

print(adespatial::dist.ldc(x12,method = "wishart"))#0.6569256


wiqid::distSimRatio(x1,x2) # Similarity ratio 0.6569256

#### Percentage difference ####
adespatial::beta.div(x12,method = "percentdiff",save.D = T)$D #0.5310021

print(adespatial::dist.ldc(x12,method = "percentdiff"))

#### Permet de calculer BDtot, SCBD, LCBD, en fonction de l'indice choisi ####
adespatial::beta.div(x12,method = "",save.D = T)

#### Podani & Schmera 2011 Jaccard ####
adespatial::beta.div.comp(test,coef = "N")$repl # 0.2424242
?adespatial::beta.div.comp(x12,coef = "N")$rich # 0.7575758
adespatial::beta.div.comp(x12,coef = "N")$D # 0.3333333

diverse::dis_entities(t(x12),method = "Podani",category_row = T)[1,2] #marche

#### Mean character difference ####
?abdiv::mean_character_difference(as.numeric(x1),as.numeric(x2)) #2.160455

#### Modified mean character difference ####
print(adespatial::dist.ldc(x12,method = "modmeanchardiff")) #2.880606

abdiv::modified_mean_character_difference(as.numeric(x1),as.numeric(x2)) #2.880606

#### Rao dissimilarity coefficient ####
?ade4::disc(as.data.frame(t(test_pa))) #0.4505909

ClusterR::distance_matrix(test_pa,method = "Rao_coefficient") #1

#### Sokal & Sneath ####
adiv::dsimcom(x12_pa,method = "1",type = "dissimilarity") #0.7926071

MultivariateAnalysis::Distancia(x12_pa,Metodo = 14)[1] #0

MultBiplotR::BinaryDistances(as.matrix(x12_pa),coefficient = 8) #0.3779645

PERMANOVA::DistBinary(x12_pa,coefficient = 8,transformation = 3)$D #0.8571429

Mercator::binaryDistance(t(x12_pa),metric = "sokalMichener") #-34.62558

#### Sokal & Sneath 1963 S5####
ade4::dist.binary(x11,method = 3) #0.7071068


#### Sokal & Sneath 1963 S13####
ade4::dist.binary(x12,method = 8) #0.6809189

MultBiplotR::BinaryDistances(as.matrix(x12_pa),coefficient = 13) #0.6809189

PERMANOVA::DistBinary(x12_pa,coefficient = 13,transformation = 1)$D #0.5363494 marche pas

#### Sokal & Sneath ####
abdiv::sokal_sneath(as.numeric(x1),as.numeric(x2)) #0.5

#### Simple match ####
ade4::dist.binary(x12,method = 2) #0.5

neighbr::similarity(x12_pa[1,],x12_pa[2,],measure = "simple_matching") #0.75

diverse::dis_entities(t(x12),method = "simple matching",category_row = T)[1,2] #marche pas

proxyC::simil(as.matrix(x1),as.matrix(x1),method = "simple matching") #0.75

shipunov::SM.dist(x12) #0.7045455

nomclust::sm(x12) #0.7045455

wiqid::distMatching(x1,x2) #0.25

proxy::dist(x12,method = "simple matching") #0.25

MultBiplotR::BinaryDistances(as.matrix(x12_pa),coefficient = 4) #0.5

PERMANOVA::DistBinary(x12_pa,coefficient = 4,transformation = 1)$D #0.75

ClusterR::distance_matrix(x12,method = "simple_matching_coefficient") #0.75

arules::dissimilarity(as.matrix(x12_pa),method = "matching") #0.25 Sokal & Michener 1958

#### Rogers & Tanimoto 1960 ####
?ade4::dist.binary(x12,method = 4) #0.6324555

abdiv::rogers_tanimoto(as.numeric(x1),as.numeric(x1)) #0.4
abdiv::sokal_michener(as.numeric(x1),as.numeric(x2)) #0.4

diverse::dis_entities(t(x12),method = "Tanimoto",category_row = T)[1,2] #marche pas
proxy::dist(x12,method = "Tanimoto") #0.4

neighbr::similarity(x12_pa[1,],x12_pa[2,],measure = "tanimoto") #0.6

wiqid::distRogersTanimoto(x1,x1) #0.4

MultivariateAnalysis::Distancia(x12,Metodo = 15)[1] #0

MultBiplotR::BinaryDistances(as.matrix(x12_pa),coefficient = 6) #0.6324555

PERMANOVA::DistBinary(x11_pa,coefficient = 6,transformation = 1)$D #0.6

philentropy::tanimoto(as.numeric(x1),as.numeric(x2)) #0.6936661

shipunov::SM.dist(x12) #0.70455455

#### Hamann coefficient ####
ade4::dist.binary(test_pa,method = 6) #0.7071068 non

diverse::dis_entities(t(x12),method = "Hamman",category_row = T)[1,2] #marche pas

proxy::dist(test_pa,method = "Hamman") #0.5 ok

MultivariateAnalysis::Distancia(test_pa,Metodo = 19)[1] #0 non

MultBiplotR::BinaryDistances(as.matrix(test_pa),coefficient = 9) #0.7071068 non

PERMANOVA::DistBinary(test_pa,coefficient = 9,transformation = 1)$D #70.25116 ok

proxyC::simil(as.matrix(test_pa[1,]),as.matrix(test_pa[2,]),method = "hamann") #0.5 ok

#### Phi of Pearson ####
ade4::dist.binary(x12,method = 9) #0.7250569 P/A

diverse::dis_entities(t(x12),method = "Phi",category_row = T)[1,2] # marche pas
diverse::dis_entities(t(x12),method = "Pearson",category_row = T)[1,2] # marche pas

proxy::dist(x12_pa,method = "Phi") #0.5257075 P/A et sans sqrt : transformation pas bonne voir permanova

proxy::dist(x12,method = "Phi-squared") #-12.94815 oui car ² = 0.5257075 P/A et sans sqrt : transformation pas bonne

MultBiplotR::BinaryDistances(as.matrix(x12_pa),coefficient = 14) #0.7250569 P/A

PERMANOVA::DistBinary(x12_pa,coefficient = 14,transformation = 3)$D # 0.7250569 P/A

philentropy::pearson_chi_sq(as.numeric(x1),as.numeric(x2)) #687532.4

#### S2 coeff Gower & Legendre ####
ade4::dist.binary(x12,method = 10) #0.7071068

#### Russel-Rao ####
abdiv::russel_rao(as.numeric(x1),as.numeric(x2)) #0.5

diverse::dis_entities(t(x12),method = "Russel",category_row = T)[1,2] #marche pas

proxy::dist(x12,method = "Russel") #0.5

MultivariateAnalysis::Distancia(x12,Metodo = 16)[1] #1

MultBiplotR::BinaryDistances(as.matrix(x12_pa),coefficient = 2) #0.7071068 correspond pas

PERMANOVA::DistBinary(x12_pa,coefficient = 2,transformation = 1)$D #0.50

Mercator::binaryDistance(t(x12),metric = "russellRao") #-18.34711

#### Yule ####
abdiv::yule_dissimilarity(as.numeric(x1),as.numeric(x1)) #0.2074074

diverse::dis_entities(t(x12),method = "Yule",category_row = T)[1,2]
diverse::dis_entities(t(x12),method = "Yule2",category_row = T)[1,2]

proxy::dist(x11,method = "Yule") #0.2074074
proxy::dist(x11,method = "Yule2") #0.5076305

MultivariateAnalysis::Distancia(x12,Metodo = 20) # marche pas

psych::Yule(varespec) #marche pas

MultBiplotR::BinaryDistances(as.matrix(x12_pa),coefficient = 15) #0.45542

PERMANOVA::DistBinary(x12_pa,coefficient = 15,transformation = 1)$D #0.7925926

#### Brainerd-Robinson ####
?tabula::similarity(x11,method = "brainerd") #93.60678
tabula::index_brainerd(as.numeric(x1),as.numeric(x1)) #93.6067

brsim::brsim(x12)$BR.similarity.matrix


#### Autres turnover ####
tabula::turnover(x12,"whittaker")
tabula::turnover(x12,"cody")
tabula::turnover(x12,"routledge1")
tabula::turnover(x12,"routledge2")
tabula::turnover(x12,"routledge3")
tabula::turnover(x12,"wilson")

proxy::dist(x12,method = "Podani") # Ne sait pas a quoi ca correspond

#### Marczewski-Steinhaus ####
adiv::distMS(x12) #0.6936661

labdsv::dsvdis(x12, index = "steinhaus") #0.333333

CommEcol::dis.goodall(x12,p.simi = "steinhaus",approach = "proportion") #non
CommEcol::dis.goodall(x12,p.simi = "steinhaus",approach = "chisquare") #1

#### Preston's coefficient of faunal dissimilarity ####
?wiqid::distPreston(x1,x1) #0.2568745

#### Braun Blanquet ####
?fossil::braun.blanquet(x1,x2) #0.7586207

?proxy::dist(test,method = "Braun-Blanquet") #0.2413793

#### Roberts ####
?labdsv::dsvdis(x11, index = "roberts") #0.6496285

#### Pearson ####
amap::Dist(x12,method = "pearson") #0.4833802

amap::Dist(x11,method = "abspearson") #0.4833802 absolute pearson

amap::Dist(x12,method = "correlation") #0.54278 centered pearson

amap::Dist(x12,method = "abscorrelation") #0.54278 absolute centered pearson

ChemoSpecUtils::rowDist(as.matrix(x12),method = "pearson") #0.4833802

ChemoSpecUtils::rowDist(as.matrix(x12),method = "abspearson") #0.4833802

ChemoSpecUtils::rowDist(as.matrix(x11),method = "correlation") ##0.54278

ChemoSpecUtils::rowDist(as.matrix(x12),method = "abscorrelation") #0.54278 

Mercator::binaryDistance(t(x12),metric = "pearson") #-2254.882

proxyC::simil(as.matrix(x1),as.matrix(x2),method = "correlation") #0.45722

arules::dissimilarity(as.matrix(x12_pa),method = "pearson") #1

arules::dissimilarity(as.matrix(x12_pa),method = "phi") #1

ClusterR::distance_matrix(x12,method = "pearson_correlation") #0.54278

rdist::rdist(x12,metric = "correlation") #0.520951
rdist::rdist(x11,metric = "absolute_correlation") #0.8893536

spaa::sp.pair(t(as.matrix(x12)))$Pearson #0.45722

hyperSpec::pearson.dist(x12) #0.27139

ldt::s.distance(t(x12),distance = "correlation",correlation = "pearson") #0.520951
ldt::s.distance(t(x12),distance = "absCorrelation",correlation = "pearson") #0.8893536

cor(t(x1),t(x2),method = "pearson")  #0.45722

#### Spearman ####
amap::Dist(x12,method = "spearman") #4104
dynutils::calculate_similarity(x12,method = "spearman") #0.1643567

spaa::sp.pair(t(as.matrix(x12)))$Spearman #0.6712866
cor(t(x1),t(x2),method = "spearman") #0.6712866

ldt::s.distance(t(x12),distance = "correlation",correlation = "spearman") #0.3811327
ldt::s.distance(t(x12),distance = "absCorrelation",correlation = "spearman") #0.7047299

ChemoSpecUtils::rowDist(as.matrix(x12),method = "spearman") #4104

fAssets::spearmanDist(t(x12)) #0.3287134

Rankcluster::distSpearman(x1,x2) #1630.034

#### Kendall ####
amap::Dist(x12,method = "kendall") #0.2114165

Rankcluster::distKendall(x1,x2) #marche pas

ChemoSpecUtils::rowDist(as.matrix(x12),method = "kendall") #0.2114165

cor(t(x1),t(x2),method = "kendall") #0.563432

fAssets::kendallDist(t(x12)) #0.436568

analogue::distance(x1,x2,method = "kendall") #95.06

#### Anderberg ####
?MultBiplotR::BinaryDistances(as.matrix(test_pa),coefficient = 5) #0.7071068 sim
MultBiplotR::BinaryDistances(as.matrix(test_pa),coefficient = 11) #0.5125381 sim

PERMANOVA::DistBinary(test_pa,coefficient = 5,transformation = 1)$D #0.5 dsim
?PERMANOVA::DistBinary(test_pa,coefficient = 11,transformation = 1)$D #0.7373047 dsim

### Kullback Leibler 
Rfast::Dist(x12,method = "kullback_leibler") #172.0308 # non coherent avec formule mais correspond avec Jeffreys
proxyC::dist(as.matrix(x12),method = "kullback") #NA
EnvNJ::metrics(t(x12),method = "kullback_leibler") #marche pas
bapred::kldist(x12) #marche pas

#### Harvesine
Rfast::Dist(x12_pa,method = "haversine") #non fonctionnel

#### Harmonic mean
Rfast::Dist(x12,method = "harmonic_mean") #9.510365


bioregion::dissimilarity(as.matrix(x12),metric = "Simpson") #0.1538462
