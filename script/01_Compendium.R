#_______________________________________________________________________________
# Title              : 01_Compendium.r
# Date               : 19/01/2025
# Object             : Script to determine a list of reliable R packages to 
#                      calculate diversity indices
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
# Brillouin ####
tabula_DBRI <- as.numeric(tabula::heterogeneity(mite,method = "brillouin"))
abdiv_DBRI <- apply(mite, 1, abdiv::brillouin_d)
ecodive_DBRI <-ecodive::alpha_div(mite,metric = "brillouin")
wiqid_DBRI <- apply(mite, 1, wiqid::biodBrillouin) #FALSE

# Put all the values in a single dataframe
DBRI<- ls(pattern = "_DBRI$")
DBRI <- mget(DBRI)
DBRI <- as.data.frame(DBRI)

# Cuba ####
# None 
# Fisher alpha ####
vegan_DFIS <- vegan::fisher.alpha(mite)
microbiome_DFIS <- microbiome::diversity(t(mite), index = "fisher")$fisher 
sads_DFIS <- rep(1,70)
for (i in 1:nrow(mite)){
  sads_DFIS[i] <- sads::fitsad(mite[i, ][mite[i, ] != 0],sad = "ls")@coef
}
ecodive_DFIS <- ecodive::fisher(mite)
BiodiversityR_DFIS <- BiodiversityR::diversityresult(mite,y=NULL,index="Logalpha",method = "each site")$Logalpha
preseqR_DFIS <- rep(1,70)
for (i in 1:nrow(mite)){
  mat_need <- as.data.frame(t(as.data.frame(rbind(freq = mite[i, ]/sum(mite_relat[i, ]), numb = mite[i, ]))))
  mat_need$freq <- as.numeric(mat_need$freq)
  mat_need <- mat_need[order(mat_need$freq,decreasing = T), ]
  preseqR_DFIS[i] <- preseqR::fisher.alpha(t(mat_need))
}

# Put all the values in a single dataframe
DFIS<- ls(pattern = "_DFIS$")
DFIS<- mget(DFIS)
DFIS <- as.data.frame(DFIS)

# Gleason ####
gleason <- function(x) {
  S <- sum(x > 0)    # nombre d'espèces présentes
  N <- sum(x)        # total d'individus
  if (N > 1) {
    return(S / log(N))
  } else {
    return(NA)
  }
}

DGLE <- apply(mite,1,gleason)

# Good ####
good_index <- function(p, m, n) {
  if (any(p <= 0 | p > 1)) {
    stop("Les probabilités doivent être dans l'intervalle (0, 1].")
  }
  if (m <= 0 || n <= 0) {
    stop("m et n doivent être des entiers naturels non nuls.")
  }
  return(sum(p^m * (-log(p))^n))
}
DGOO_1_1 <- rep(1,70)
for (i in 1:nrow(mite)){
  DGOO_1_1[i] <- good_index(mite_relat[i, ][mite_relat[i, ] != 0],1,1)
}
DGOO_1_2 <- rep(1,70)
for (i in 1:nrow(mite)){
  DGOO_1_2[i] <- good_index(mite_relat[i, ][mite_relat[i, ] != 0],1,2)
}

DGOO_1_3 <- rep(1,70)
for (i in 1:nrow(mite)){
  DGOO_1_2[i] <- good_index(mite_relat[i, ][mite_relat[i, ] != 0],1,3)
}

DGOO_2_1 <- rep(1,70)
for (i in 1:nrow(mite)){
  DGOO_2_1[i] <- good_index(mite_relat[i, ][mite_relat[i, ] != 0],2,1)
}
DGOO_2_2 <- rep(1,70)
for (i in 1:nrow(mite)){
  DGOO_1_2[i] <- good_index(mite_relat[i, ][mite_relat[i, ] != 0],2,2)
}

DGOO_2_3 <- rep(1,70)
for (i in 1:nrow(mite)){
  DGOO_2_2[i] <- good_index(mite_relat[i, ][mite_relat[i, ] != 0],2,3)
}

DGOO_3_1 <- rep(1,70)
for (i in 1:nrow(mite)){
  DGOO_3_1[i] <- good_index(mite_relat[i, ][mite_relat[i, ] != 0],3,1)
}
DGOO_3_2 <- rep(1,70)
for (i in 1:nrow(mite)){
  DGOO_3_2[i] <- good_index(mite_relat[i, ][mite_relat[i, ] != 0],3,2)
}

DGOO_3_3 <- rep(1,70)
for (i in 1:nrow(mite)){
  DGOO_3_2[i] <- good_index(mite_relat[i, ][mite_relat[i, ] != 0],3,3)
}

# Put all the values in a single dataframe
DGOO<- ls(pattern = "DGOO_")
DGOO<- mget(DGOO)
DGOO <- as.data.frame(DGOO)

# Kothe ####
kothe <- function(matrix_data) {
  S_i <- apply(matrix_data, 1, function(x) sum(x > 0))
  S_max <- max(S_i)
  (S_max - S_i) / S_max
}
DKO <- kothe(mite)

# Log normal lambda ####
sads_DLMD <- rep(1,70)
for (i in 1:nrow(mite)){
  sads_DLMD[i] <- sads::fitsad(mite[i, ][mite[i, ] != 0],sad = "lnorm")@coef
}
asbio_DLMD <- numeric(nrow(mite))
for (i in seq_len(nrow(mite))) {
  asbio_DLMD[i] <- tryCatch(
    {
      as.numeric(
        asbio::Preston.dist(
          as.numeric(mite[i, ]),
          plot = FALSE
        )$Est.no.of.spp
      )
    },
    error = function(e) NA
  )
}
vegan_DLMD <- rep(1,70)
for (i in 1:nrow(mite)){
  vegan_DLMD[i] <- vegan::veiledspec(vegan::prestonfit(mite[i,]))
}
Compositional_DLMD <- rep(1,70)
for (i in 1:nrow(mite)){
  Compositional_DLMD[i] <- as.numeric(Compositional::alfa.tune(mite[i,])[1])
} #FALSE

# Put all the values in a single dataframe
DLMD <- ls(pattern = "_DLMD$")
DLMD <- mget(DLMD)
DLMD <- as.data.frame(DLMD)

# Margalef ####
abdiv_DMG <- apply(mite,1,abdiv::margalef)
agricolae_DMG <- rep(1,70)
for (i in 1:nrow(mite)){
  agricolae_DMG[i] <- agricolae::index.bio(mite[i,],method = "Margalef",nboot=0,console = F)$index
} #FALSE
adiv_DMG <- adiv::speciesdiv(mite,method = "Margalef")[,1]
tabula_DMG <- tabula::richness(mite,method = "margalef")@.Data
ecodive_DMG <- ecodive::margalef(mite)
benthos_DSHA <- rep(1,70)
for (i in 1:nrow(mite)){
  benthos_DSHA[i] <- benthos::margalef(taxon = names(mite[i,]), count = mite[i,])
}

# Put all the values in a single dataframe
DMG <- ls(pattern = "_DMG$")
DMG <- mget(DMG)
DMG <- as.data.frame(DMG)

# McIntosh ####
forestHES_DMC <- forestHES::mcIntosh(mite) #0.6461754 faux
abdiv_DMC <- apply(mite,1,abdiv::mcintosh_d)
tabula_DMC <- rep(1,70)
for (i in 1:nrow(mite)){
  tabula_DMC[i] <- tabula::index_mcintosh(as.numeric(mite[i,]),evenness = F)
}
adiv_DMC <- adiv::speciesdiv(mite,method = "McIntosh")[,1]
agricolae_DMC <- rep(1,70)
for (i in seq_len(nrow(mite))) {
  agricolae_DMC[i] <- tryCatch(
    {
      as.numeric(
        agricolae::index.bio(mite[i,],method = "McIntosh")$index
      )
    },
    error = function(e) NA
  )
}
ecodive_DMC <- ecodive::alpha_div(mite,metric = "mcintosh")

# Put all the values in a single dataframe
DMC <- ls(pattern = "_DMC")
DMC <- mget(DMC)
DMC <- as.data.frame(DMC)

# Menhinick ####
abdiv_DMN <- apply(mite,1,abdiv::menhinick)
tabula_DMN <- tabula::richness(mite,method = "menhinick")@.Data
adiv_DMN <- adiv::speciesdiv(mite,method = "Menhinick")[,1]
ecodive_DMN <- ecodive::menhinick(mite)

# Put all the values in a single dataframe
DMN <- ls(pattern = "_DMN")
DMN <- mget(DMN)
DMN <- as.data.frame(DMN)


# Odum ####
odum <- function(x) {
  S <- sum(x > 0)       # nombre d'espèces présentes
  N <- sum(x)           # total d'individus
  (S / N) * 1000
}
OD <- apply(mite,1,odum)

# Q-statistic ####
abdiv_DQ <- rep(1,70)
for (i in 1:nrow(mite)){
  abdiv_DQ[i] <- abdiv::kempton_taylor_q(as.numeric(mite[i,]),lower_quantile = 0.5,upper_quantile = 0.9)
}

# Put all the values in a single dataframe
DQ <- ls(pattern = "_DQ")
DQ <- mget(DQ)
DQ <- as.data.frame(DQ)

# Rao ####
chemodiv_DRAO <- chemodiv::calcDiv(mite,compDisMat = mite_dist, type = "RaoQ") 
chemodiv_DRAO <- chemodiv_DRAO$RaoQ
diverse_DRAO <- diverse::diversity(t(mite), type = "rao", category_row = TRUE, dis = mite_dist)[[1]] #FALSE
rao_from_abundance <- function(abundance_matrix, dist_method = "euclidean") {
  # Vérification du format
  abundance_matrix <- as.data.frame(abundance_matrix)
  abundance_matrix[] <- lapply(abundance_matrix, as.numeric)  # conversion en numérique
  
  # Vérification des noms de lignes
  if (is.null(rownames(abundance_matrix))) {
    rownames(abundance_matrix) <- paste0("Site_", seq_len(nrow(abundance_matrix)))
  }
  
  # Calcul de la matrice de dissimilarité entre espèces
  dissimilarity_matrix <- as.matrix(dist(t(abundance_matrix), method = dist_method))
  
  # Initialiser le vecteur des résultats
  rao_values <- numeric(nrow(abundance_matrix))
  names(rao_values) <- rownames(abundance_matrix)
  
  # Boucle sur chaque site
  for (i in 1:nrow(abundance_matrix)) {
    abundances <- as.numeric(abundance_matrix[i, ])
    total <- sum(abundances)
    if (total == 0) {
      rao_values[i] <- NA
      next
    }
    
    # Proportions
    p <- abundances / total
    
    # Calcul de Rao via produit matriciel
    rao <- sum(outer(p, p) * dissimilarity_matrix)
    rao_values[i] <- rao
  }
  
  return(rao_values)
}
custom_DRAO <- rao_from_abundance(mite) 
BAT_DRAO <- BAT::rao(as.matrix(mite))[,1] 
ade4_DRAO <- rep(1,70)
for (i in 1:nrow(mite)){
  ade4_DRAO[i] <- ade4::apqe(as.data.frame(t(mite[i,])),dis = NULL)$results[2,1]
}
SYNCSA_DRAO <- SYNCSA::rao.diversity(mite)$Simpson

# Put all the values in a single dataframe
DRAO <- ls(pattern = "_DRAO")
DRAO <- mget(DRAO)
DRAO <- as.data.frame(DRAO)

# Rygg ####
benthos_DRY <- rep(1,70)
for (i in 1:nrow(mite)){
  benthos_DRY[i] <- benthos::rygg(taxon = names(mite[i,]), count = mite[i,],adjusted = F)
}

# Rygg adjusted ####
benthos_DRYA <- rep(1,70)
for (i in 1:nrow(mite)){
  benthos_DRYA[i] <- benthos::rygg(taxon = names(mite[i,]), count = mite[i,],adjusted = T)
}

# Shannon
vegan_DSHA <- vegan::diversity(mite,index = "shannon")
biosampleR_DSHA <- biosampleR::calc_diversity_indices(mite)$shannon
divDyn_DSHA <- divDyn::indices(as.matrix(mite),method = "shannon")[,1]
sprex_DSHA <- sprex::diversity(t(mite),type = "shannon")
EconGeo_DSHA <- EconGeo::entropy(mite) #FALSE
BiodiversityR_DSHA <- BiodiversityR::diversityresult(mite,y=NULL,index="Shannon",method = "each site")$Shannon
tabula_DSHA <- tabula::diversity(mite)$shannon
iNEXT_DSHA <- iNEXT::ChaoShannon(t(mite),datatype = "abundance",transform = F)$Observed
adiv_DSHA <- adiv::speciesdiv(mite,method = "Shannon")[,1]
diverse_DSHA <- diverse::diversity(t(mite),type = "entropy",category_row = T)$entropy
asbio_DSHA <- asbio::alpha.div(mite,"shan")
ecodive_DSHA <- ecodive::shannon(mite)
microbiome_DSHA <- microbiome::diversity(t(mite), index = "shannon")$shannon
chemodiv_DSHA <- chemodiv::calcDiv(mite,type = "Shannon")$Shannon # OK
aqp_DSHA <- apply(mite, 1, aqp::shannonEntropy) # FALSE
wiqid_DSHA <- apply(mite,1,wiqid::biodSimpson) # FALSE
OTUtable_DSHA <- apply(mite, 1, OTUtable::shannon)
aqp_DSHA <- apply(mite, 1, aqp::shannonEntropy) #FALSE
DescTools_DSHA <- apply(mite, 1, DescTools::Entropy) #FALSE
pgirmess_DSHA <- apply(mite, 1, pgirmess::shannon)[1,] #FALSE
divent_DSHA <- apply(mite, 1, function(x) divent::ent_shannon(x, estimator = "naive")$entropy)  # OK
SpiecEasi_DSHA <- apply(mite, 1, SpiecEasi::shannon) #FALSE
abdiv_DSHA <- apply(mite, 1, abdiv::shannon) #FALSE
HardyWeinberg_DSHA <- apply(mite, 1, function(x) HardyWeinberg::shannon(x)$Hp)
MCPAN_DSHA <- apply(mite, 1, function(x) MCPAN::estShannon(x)$estraw)
entropart_DSHA <- apply(mite, 1, entropart::Shannon)
wavethresh_DSHA <- apply(mite, 1, wavethresh::Shannon.entropy)
codyn_DSHA <- as.numeric(t(as.data.frame(apply(mite, 1, function(x) codyn::community_diversity(
  data.frame(Value = x), abundance.var = "Value", metric = "Shannon")
))))
benthos_DSHA <- rep(1,70)
for (i in 1:nrow(mite)){
  benthos_DSHA[i] <- benthos::shannon(taxon = names(mite[i,]), count = mite[i,])
}
triversity_DSHA <- rep(1,70) 
for (i in 1:nrow(mite)){
  triversity_DSHA[i] <- triversity::get_diversity_from_distribution(mite_relat[i,])
} #FALSE

# Put all the values in a single dataframe
DSHA <- ls(pattern = "_DSHA$")
DSHA <- mget(DSHA)
DSHA <- as.data.frame(DSHA)

# Simpson ####
sprex_DSP <- sprex::diversity(t(mite),type = "simpson")
tabula_DSP <- rep(1,70) 
for (i in 1:nrow(mite)){
  tabula_DSP[i] <- tabula::index_simpson(as.numeric(mite[i,]),eveness = F, unbiaised = F, na.rm=F)
}
breakaway_DSP <- rep(1,70) 
for (i in 1:nrow(mite)){
  breakaway_DSP[i] <- breakaway::true_simpson(mite_relat[i,])
}
benthos_DSP <- rep(1,70)
for (i in 1:nrow(mite)){
  benthos_DSP[i] <- benthos::simpson(taxon = names(mite[i,]), count = as.numeric(mite[i,]))
}
diverse_DSP <- diverse::diversity(t(mite),type = "simpson",category_row = T)$simpson.D
# A regler ####
#agricolae_DSP <- rep(1,70)
#for (i in 1:nrow(mite)){
  agricolae_DSP[i] <- agricolae::index.bio(mite[i,],method = "Simpson.Dom")$index
}
#####
wiqid_DSP <- rep(1,70)
for (i in 1:nrow(mite)){
  wiqid_DSP[i] <- wiqid::biodSimpson(abVec=as.numeric(mite[i,]))
}#FALSE
abdiv_DSP <- rep(1,70)
for (i in 1:nrow(mite)){
  abdiv_DSP[i] <- abdiv::dominance(mite[i,])
}
microbiome_DSP <- rep(1,70)
for (i in 1:nrow(mite)){
  microbiome_DSP[i] <- microbiome::dominance(as.numeric(mite[i,]),index = "simpson")$simpson
}

# Put all the values in a single dataframe
DSP <- ls(pattern = "_DSP$")
DSP <- mget(DSP)
DSP <- as.data.frame(DSP)

# Gini-Simpson ####
vegan_D1SP <- vegan::diversity(mite,index = "simpson")
biosampleR_D1SP <- biosampleR::calc_diversity_indices(mite)$simpson
sprex_D1SP <- rep(1,70)
for (i in 1:nrow(mite)){
  sprex_D1SP[i] <- sprex::diversity(as.numeric(mite[i,]),type = "gini.simpson")
}
abdiv_D1SP <- apply(mite,1,abdiv::simpson )
ecodive_D1SP <- ecodive::simpson(mite)
MCPAN_D1SP <- apply(mite,1,MCPAN::Simpson)
entropart_D1SP <- rep(1,70)
for (i in 1:nrow(mite)){
  entropart_D1SP[i] <- entropart::Simpson(as.numeric(mite[i,]),Correction = "None")
}
BiodiversityR_D1SP <- BiodiversityR::diversityresult(mite,y=NULL,index="Simpson",method = "each site")$Simpson
simboot_D1SP <- rep(1,70)
for (i in 1:nrow(mite)){
  simboot_D1SP[i] <- simboot::Simpson(mite[i,])
}
concstats_D1SP <- rep(1,70)
for (i in 1:nrow(mite)){
  concstats_D1SP[i] <- concstats::concstats_simpson(as.numeric(mite_relat[i,]),na.rm = T)
}
untb_D1SP <- rep(1,70)
for (i in 1:nrow(mite)){
  untb_D1SP[i] <- untb::simpson(mite[i,],with.replacement = F)#
}
DescTools_D1SP <- rep(1,70)
for (i in 1:nrow(mite)){
  DescTools_D1SP[i] <- DescTools::Gini(as.numeric(mite[i,]),na.rm = F)
}
adiv_D1SP <- adiv::speciesdiv(mite,method = "GiniSimpson")[,1]
diverse_D1SP <- diverse::diversity(t(mite),type = "gini-simpson",category_row = T)[[1]]
asbio_D1SP <- asbio::alpha.div(mite,"simp")
# A regler ####
#agricolae_D1SP <- rep(1,70)
#for (i in 1:nrow(mite)){
  agricolae_D1SP[i] <- agricolae::index.bio(mite[i,],method = "Simpson.Div")$index
}
#####
microbiome_D1SP <- microbiome::diversity(t(mite), index = "gini_simpson")$gini_simpson
divent_D1SP <- rep(1,70)
for (i in 1:nrow(mite)){
  divent_D1SP[i] <- divent::ent_simpson(as.numeric(mite[i,]),estimator = "naive")$entropy
}

# Put all the values in a single dataframe
D1SP <- ls(pattern = "_D1SP$")
D1SP <- mget(D1SP)
D1SP <- as.data.frame(D1SP)

# Inverse Simpson ####
sprex_DINVSP <- rep(1,70)
for (i in 1:nrow(mite)){
  sprex_DINVSP[i] <- sprex::diversity(as.numeric(mite[i,]),type = "inv.simpson")
}
vegan_DINVSP <- vegan::diversity(mite,index = "invsimpson")
ecodive_DINVSP <- ecodive::inv_simpson(mite)
abdiv_DINVSP <- rep(1,70)
for (i in 1:nrow(mite)){
  abdiv_DINVSP[i] <- abdiv::invsimpson(mite[i,])
}
breakaway_DINVSP <- rep(1,70) 
for (i in 1:nrow(mite)){
  breakaway_DINVSP[i] <- breakaway::true_inverse_simpson(mite_relat[i,])
}
BiodiversityR_DINVSP <- BiodiversityR::diversityresult(mite,y=NULL,index="inverseSimpson",method = "each site")$inverseSimpson
adiv_DINVSP <- adiv::speciesdiv(mite,method = "Simpson")[,1]
asbio_DINVSP <- asbio::alpha.div(mite,"inv.simp")
microbiome_DINVSP <- microbiome::diversity(t(mite), index = "inverse_simpson")[,1] # OK
divseg_DINVSP <- divseg::ds_simpson(mite,.cols = dplyr::everything()) # OK 
chemodiv_DINVSP <- chemodiv::calcDiv(mite,type = "Simpson")[,1] # OK
codyn_DINVSP <- rep(1,70) 
for (i in 1:nrow(mite)){
  mite_long <- tidyr::pivot_longer(mite[i,],cols = colnames(mite)[1]:colnames(mite)[35],names_to = "Esp",values_to = "Value")
  mite_long <- mite_long[1:35,]
  codyn_DINVSP[i] <- codyn::community_diversity(mite_long,abundance.var = "Value",metric = "InverseSimpson")$InverseSimpson # OK
}

# Put all the values in a single dataframe
DINVSP <- ls(pattern = "_DINVSP$")
DINVSP <- mget(DINVSP)
DINVSP <- as.data.frame(DINVSP)

# Hill ####
vegan_QHIL_0 <- vegan::renyi(mite,scales = 0,hill = T)
vegan_QHIL_1 <- vegan::renyi(mite,scales = 1,hill = T)
vegan_QHIL_2 <- vegan::renyi(mite,scales = 2,hill = T)
vegan_QHIL_3 <- vegan::renyi(mite,scales = 3,hill = T)

adiv_QHIL_0 <- as.numeric(adiv::divparam(mite,method = "hill",q=0))
adiv_QHIL_1 <- as.numeric(adiv::divparam(mite,method = "hill",q=1))
adiv_QHIL_2 <- as.numeric(adiv::divparam(mite,method = "hill",q=2))
adiv_QHIL_3 <- as.numeric(adiv::divparam(mite,method = "hill",q=3))


divDyn_QHIL_2 <- divDyn::indices(as.matrix(mite),method = "hill2")[,1] # que pour q = 2

sprex_QHIL_0 <- rep(1,70) 
for (i in 1:nrow(mite)){
  sprex_QHIL_0[i] <- sprex::diversity(as.numeric(mite[i,]),type = "hill",q=0)
}
sprex_QHIL_1 <- rep(1,70) 
for (i in 1:nrow(mite)){
  sprex_QHIL_1[i] <- sprex::diversity(as.numeric(mite[i,]),type = "hill",q=1)
}
sprex_QHIL_2 <- rep(1,70) 
for (i in 1:nrow(mite)){
  sprex_QHIL_2[i] <- sprex::diversity(as.numeric(mite[i,]),type = "hill",q=2)
}
sprex_QHIL_3 <- rep(1,70) 
for (i in 1:nrow(mite)){
  sprex_QHIL_3[i] <- sprex::diversity(as.numeric(mite[i,]),type = "hill",q=3)
}

EntropyEstimation_QHIL_1 <- rep(1,70) 
for (i in 1:nrow(mite)){
  EntropyEstimation_QHIL_1[i] <- EntropyEstimation::Hill.z(as.numeric(mite[i,]),r=1)
}
EntropyEstimation_QHIL_2 <- rep(1,70) 
for (i in 1:nrow(mite)){
  EntropyEstimation_QHIL_2[i] <- EntropyEstimation::Hill.z(as.numeric(mite[i,]),r=2)
}
EntropyEstimation_QHIL_3 <- rep(1,70) 
for (i in 1:nrow(mite)){
  EntropyEstimation_QHIL_3[i] <- EntropyEstimation::Hill.z(as.numeric(mite[i,]),r=3)
}

extremefit_QHIL <- rep(1,70) 
for (i in 1:nrow(mite)){
  mite_long <- tidyr::pivot_longer(mite[i,],cols = colnames(mite)[1]:colnames(mite)[35],names_to = "Esp",values_to = "Value")
  mite_long <- mite_long[1:35,]
  extremefit_QHIL[i] <- extremefit::hill(mite_long$Value,weights = rep(1,nrow(mite_long)),grid = mite_long$Value)$hill
}

BAT_QHIL_0 <- BAT::hill(mite,q=0)[,1]
BAT_QHIL_1 <- BAT::hill(mite,q=1)[,1]
BAT_QHIL_2 <- BAT::hill(mite,q=2)[,1]
BAT_QHIL_3 <- BAT::hill(mite,q=3)[,1]

breakaway_QHIL_0 <- rep(1,70) 
for (i in 1:nrow(mite)){
  breakaway_QHIL_0[i] <- breakaway::true_hill(as.numeric(mite_relat[i,]),q=0)
}
breakaway_QHIL_1 <- rep(1,70) 
for (i in 1:nrow(mite)){
  breakaway_QHIL_1[i] <- breakaway::true_hill(as.numeric(mite_relat[i,]),q=1)
}
breakaway_QHIL_2 <- rep(1,70) 
for (i in 1:nrow(mite)){
  breakaway_QHIL_2[i] <- breakaway::true_hill(as.numeric(mite_relat[i,]),q=2)
}
breakaway_QHIL_3 <- rep(1,70) 
for (i in 1:nrow(mite)){
  breakaway_QHIL_3[i] <- breakaway::true_hill(as.numeric(mite_relat[i,]),q=3)
}

hillR_QHIL_0 <- hillR::hill_taxa(mite,q=0) # OK
hillR_QHIL_1 <- hillR::hill_taxa(mite,q=1) # OK
hillR_QHIL_2 <- hillR::hill_taxa(mite,q=2) # OK
hillR_QHIL_3 <- hillR::hill_taxa(mite,q=3) # OK

hilldiv_QHIL_0 <- hilldiv::hill_div(t(mite),qvalue = 0)
hilldiv_QHIL_1 <- hilldiv::hill_div(t(mite),qvalue = 1)
hilldiv_QHIL_2 <- hilldiv::hill_div(t(mite),qvalue = 2)
hilldiv_QHIL_3 <- hilldiv::hill_div(t(mite),qvalue = 3)

chemodiv_QHIL_0 <- chemodiv::calcDiv(mite,type = "HillDiv",q=0)[,1]
chemodiv_QHIL_1 <- chemodiv::calcDiv(mite,type = "HillDiv",q=1)[,1]
chemodiv_QHIL_2 <- chemodiv::calcDiv(mite,type = "HillDiv",q=2)[,1]
chemodiv_QHIL_3 <- chemodiv::calcDiv(mite,type = "HillDiv",q=3)[,1]

benthos_QHIL_0 <- rep(1,70) 
for (i in 1:nrow(mite)){
  benthos_QHIL_0[i] <- benthos::hill(taxon = colnames(mite),count = as.numeric(mite[i,]),a = 0) # OK
}
benthos_QHIL_1 <- rep(1,70) 
for (i in 1:nrow(mite)){
  benthos_QHIL_1[i] <- benthos::hill(taxon = colnames(mite),count = as.numeric(mite[i,]),a = 1) # OK
}
benthos_QHIL_2 <- rep(1,70) 
for (i in 1:nrow(mite)){
  benthos_QHIL_2[i] <- benthos::hill(taxon = colnames(mite),count = as.numeric(mite[i,]),a = 2) # OK
}
benthos_QHIL_3 <- rep(1,70) 
for (i in 1:nrow(mite)){
  benthos_QHIL_3[i] <- benthos::hill(taxon = colnames(mite),count = as.numeric(mite[i,]),a = 3) # OK
}

entropart_QHIL_0 <- rep(1,70) 
for (i in 1:nrow(mite)){
  entropart_QHIL_0[i] <- entropart::Diversity(as.numeric(mite[i,]),q=0)
}
entropart_QHIL_1 <- rep(1,70) 
for (i in 1:nrow(mite)){
  entropart_QHIL_1[i] <- entropart::Diversity(as.numeric(mite[i,]),q=1)
}
entropart_QHIL_2 <- rep(1,70) 
for (i in 1:nrow(mite)){
  entropart_QHIL_2[i] <- entropart::Diversity(as.numeric(mite[i,]),q=2)
}
entropart_QHIL_3 <- rep(1,70) 
for (i in 1:nrow(mite)){
  entropart_QHIL_3[i] <- entropart::Diversity(as.numeric(mite[i,]),q=3)
}

divent_QHIL_0 <- rep(1,70) 
for (i in 1:nrow(mite)){
  divent_QHIL_0[i] <- divent::div_hill(as.numeric(mite[i,]),q = 0)$diversity
}
divent_QHIL_1 <- rep(1,70) 
for (i in 1:nrow(mite)){
  divent_QHIL_1[i] <- divent::div_hill(as.numeric(mite[i,]),q = 1)$diversity
}
divent_QHIL_2 <- rep(1,70) 
for (i in 1:nrow(mite)){
  divent_QHIL_2[i] <- divent::div_hill(as.numeric(mite[i,]),q = 2)$diversity
}
divent_QHIL_3 <- rep(1,70) 
for (i in 1:nrow(mite)){
  divent_QHIL_3[i] <- divent::div_hill(as.numeric(mite[i,]),q = 3)$diversity
}

# Put all the values in a single dataframe
QHIL  <- ls(pattern = "_QHIL")
QHIL <- mget(QHIL)
QHIL <- as.data.frame(QHIL)


# Autres ####
abdiv::mcintosh_e(x) # ne correspond pas
tabula::index_mcintosh(as.numeric(x),evenness = T) # 0.70951058
adiv::specieseve(x,method = "McIntosh") # 0.70951058
iNEXT::ChaoSimpson(t(varespec),datatype = "abundance",transform = F)$Observed

iNEXT::ChaoSimpson(t(varespec),datatype = "abundance",transform = T)$Observed
tabula::index_simpson(as.numeric(x),eveness = T, unbiaised = T, na.rm=F)



