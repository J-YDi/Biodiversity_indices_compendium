#_______________________________________________________________________________
# Title              : 01_Compendium.r
# Date               : 16/01/2025
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
mite_relat <- mite/rowSums(mite_relat)

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
microbiome_DFIS <- microbiome::diversity(t(mite), index = "fisher") 
sads_DFIS <- rep(1,70)
for (i in 1:nrow(mite)){
  sads_DFIS[i] <- sads::fitsad(mite[i, ][mite[i, ] != 0],sad = "ls")@coef
}
ecodive_DFIS <- ecodive::fisher(mite)
BiodiversityR_DFIS <- BiodiversityR::diversityresult(mite,y=NULL,index="Logalpha",method = "each site")
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
adiv_DMG <- adiv::speciesdiv(mite,method = "Margalef")
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
adiv_DMC <- adiv::speciesdiv(mite,method = "McIntosh")
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
adiv_DMN <- adiv::speciesdiv(mite,method = "Menhinick")
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
  abdiv_Q[i] <- abdiv::kempton_taylor_q(as.numeric(mite[i,]),lower_quantile = 0.5,upper_quantile = 0.9)
}

# Put all the values in a single dataframe
DQ <- ls(pattern = "_DQ")
DQ <- mget(DQ)
DQ <- as.data.frame(DQ)

# Rao ####

varespec_dist <- as.matrix(dist(t(varespec),method = "euclidean",diag = T,upper = T))
chemodiv::calcDiv(varespec,compDisMat = varespec_dist, type = "RaoQ") # correct 69.64516
varespec_dist <- dist(t(varespec), method = "euclidean")
varespec_dist <- as.matrix(varespec_dist)

rownames(varespec_dist) <- rownames(t(varespec))
colnames(varespec_dist) <- rownames(t(varespec))
diverse::diversity(t(varespec), type = "rao", category_row = TRUE, dis = varespec_dist) # ne correspond pas 34.82258
diverse::diversity(t(varespec), type = "rao-stirling", category_row = TRUE, dis = varespec_dist) # ne correspond pas, donne la même chose que Rao 34.82258, rao-striling

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
rao_from_abundance(varespec) # 69.64516

BAT::rao(as.matrix(varespec)) # 0.8217115 Rao quadratic entropy

entropart::Rao(as.matrix(varespec),Tree = varespec_dist) # basé sur un arbre, et donc distance assez problematique

ade4::apqe(as.data.frame(t(x)),dis = NULL)$results[2,1] # 0.8217115 Rao quadratic entropy

hilldiv::index_div(index = "rao") # basé sur un arbre

entropart::Hqz(as.integer(x),Correction = "ChaoShen") # 1.83 pas sur que ca corresponde

SYNCSA::rao.diversity(varespec)$Simpson #0.8217115

divent::ent_rao(as.numeric(x),estimator = "naive") # necessite un arbere



# Autres ####
abdiv::mcintosh_e(x) # ne correspond pas
tabula::index_mcintosh(as.numeric(x),evenness = T) # 0.70951058
adiv::specieseve(x,method = "McIntosh") # 0.70951058









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


corrplot::corrplot(cor(DSHA),method = "color" ,tl.cex = 0.5,type = "full",diag = F,col = c(rep("white",10000000),"blue"))



