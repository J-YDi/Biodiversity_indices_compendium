#_______________________________________________________________________________
# Title              : 01_Listing_functions_alpha_indices.r
# Date               : 27/01/2025
# Object             : Script to create dataset of values from functions that 
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

#___________________________ Alpha diversity indices ___________________________####

# ACE ####
tabula_RACE <- rep(1,70)
for (i in 1:nrow(mite)){
  tabula_RACE[i] <- tabula::index_ace(as.numeric(mite[i,]))
}
fossil_RACE <- rep(1,70)
for (i in 1:nrow(mite)){
  fossil_RACE[i] <- fossil::ACE(mite[i,],taxa.row = F) 
}
BAT_RACE <- rep(1,70)
for (i in 1:nrow(mite)){
  BAT_RACE[i] <- BAT::alpha.estimate(mite[i,])[11]
}
vegan_RACE <- rep(1,70)
for (i in 1:nrow(mite)){
  vegan_RACE[i] <- vegan::estimateR(mite[i,],index = "chao")[4]
}
ecodive_RACE <- ecodive::ace(mite) #12.17889 
wiqid_RACE <- rep(1,70)
for (i in 1:nrow(mite)){
  wiqid_RACE[i] <- wiqid::richACE(as.numeric(mite[i,]))
}
pctax_RACE <- pctax::a_diversity(t(mite),method = "ace")[,1] # 12.178889
sprex_RACE <- rep(1,70)
for (i in 1:nrow(mite)){
  sprex_RACE[i] <- sprex::ACE(as.numeric(mite[i,]))
}
# Put all the values in a single dataframe
RACE <- ls(pattern = "_RACE$")
RACE <- mget(RACE)
RACE <- as.data.frame(RACE)

# Bootstrap ####
wiqid_RBOO <- rep(1,70)
for (i in 1:nrow(mite)){
  wiqid_RBOO[i] <- wiqid::richBoot(t(mite[i,])) #44
}
BiodiversityR_RBOO <- rep(1,70)
for (i in 1:nrow(mite)){
  BiodiversityR_RBOO[i] <- as.numeric(BiodiversityR::diversityresult(mite[i,],index = "boot",method = "each site")) #11
}

# Put all the values in a single dataframe
RBOO <- ls(pattern = "_RBOO$")
RBOO <- mget(RBOO)
RBOO <- as.data.frame(RBOO)

# Chao1 ####
vegan_RC1 <- rep(1,70)
for (i in 1:nrow(mite)){
  vegan_RC1[i] <- vegan::estimateR(mite[i,],index = "chao")[2]
}
BiodiversityR_RC1 <- rep(1,70)
for (i in 1:nrow(mite)){
  BiodiversityR_RC1[i] <- as.numeric(BiodiversityR::diversityresult(mite[i,],index = "chao",method = "each site"))
}
microbiome_RC1 <- microbiome::richness(t(mite),index = c("chao1"))[,1]
pctax_RC1 <- pctax::a_diversity(t(mite),method = "chao1")[,1]
rareNMtests_RC1 <- rep(1,70)
for (i in 1:nrow(mite)){
  rareNMtests_RC1[i] <- as.numeric(print(rareNMtests::chao1(mite[i,])))[2]
}
SpadeR_RC1 <- rep(1,70)
for (i in 1:nrow(mite)){
  SpadeR_RC1[i] <- as.numeric(SpadeR::ChaoSpecies(mite[i,],datatype = "abundance")$Basic_data_information[2,2])
}
iNEXT_RC1 <- iNEXT::ChaoRichness(t(mite),datatype = "abundance")$Estimator

BAT_RC1 <- BAT::alpha.estimate(mite)[,9] 

SSP_RC1 <- rep(1,70)
for (i in 1:nrow(mite)){
  SSP_RC1[i] <- SSP:::assempar(mite[i,],type = "counts",Sest.method = "chao")$Sest 
}
biosampleR_RC1 <- biosampleR::calc_diversity_indices(mite)[6][,1] #29

wiqid_RC1 <- rep(1,70)
for (i in 1:nrow(mite)){
  wiqid_RC1[i] <- as.numeric(wiqid::richChao1(t(mite[i,]),correct = F))[1]
}

# Put all the values in a single dataframe
RC1 <- ls(pattern = "_RC1$")
RC1 <- mget(RC1)
RC1 <- as.data.frame(RC1)

# Chao 1 modified ####

wiqid_RC1M <- rep(1,70)
for (i in 1:nrow(mite)){
  wiqid_RC1M[i] <- as.numeric(wiqid::richChao1(t(mite[i,]),correct = F))[1]
}
ecodive_RC1M <- ecodive::chao1(mite) 
BAT_RC1M <- BAT::alpha.estimate(mite)[,10]
divent_RC1M <- rep(1,70)
for (i in 1:nrow(mite)){
  divent_RC1M[i] <- divent::div_richness(as.numeric(t(mite[i,])),estimator = "Chao1")$diversity
}
entropart_RC1M <- rep(1,70)
for (i in 1:nrow(mite)){
  entropart_RC1M[i] <- entropart::bcRichness(mite[i,],Correction = "Chao1")
}
fossil_RC1M <- rep(1,70)
for (i in 1:nrow(mite)){
  fossil_RC1M[i] <- fossil::chao1(as.numeric(mite[i,])) 
}
mobr_RC1M <- mobr::calc_chao1(mite)
OTUtable_RC1M <- rep(1,70)
for (i in 1:nrow(mite)){
  OTUtable_RC1M[i] <- OTUtable::chao1(as.numeric(mite[i,]))
}
tabula_RC1M <- rep(1,70)
for (i in 1:nrow(mite)){
  tabula_RC1M[i] <- tabula::index_chao1(as.numeric(mite[i,]))
}

# Put all the values in a single dataframe
RC1M <- ls(pattern = "_RC1M$")
RC1M <- mget(RC1M)
RC1M <- as.data.frame(RC1M)

# Chao 2 ####
fossil_RCA2 <- rep(1,70)
for (i in 1:nrow(mite)){
  fossil_RCA2[i] <- fossil::chao2(as.numeric(mite[i,]))
}
wiqid_RCA2 <- rep(1,70)
for (i in 1:nrow(mite)){
  wiqid_RCA2[i] <- wiqid::richChao2(t(mite[i,]))[1]
}
tabula_RCA2 <- rep(1,70)
for (i in 1:nrow(mite)){
  tabula_RCA2[i] <- tabula::index_chao2(as.matrix(mite[i,]))
}
divDyn_RCA2 <- rep(1,70)
for (i in 1:nrow(mite)){
  divDyn_RCA2[i] <- divDyn::indices(as.matrix(mite[i,]),method = "chao2")
}
BAT_RCA2 <- BAT::alpha.accum(mite)[,18]

# Put all the values in a single dataframe
RCA2 <- ls(pattern = "_RCA2$")
RCA2 <- mget(RCA2)
RCA2 <- as.data.frame(RCA2)

# Hurlbert ####
benthos_RHUR_2 <- rep(1,70)
for (i in 1:nrow(mite)){
  benthos_RHUR_2[i] <- benthos::hurlbert(taxon = colnames(mite),count = as.numeric(mite[i,]), n=2) # OK
}
benthos_RHUR_3 <- rep(1,70)
for (i in 1:nrow(mite)){
  benthos_RHUR_3[i] <- benthos::hurlbert(taxon = colnames(mite),count = as.numeric(mite[i,]), n=3) # OK
}

entropart_RHUR_2 <- rep(1,70)
for (i in 1:nrow(mite)){
  entropart_RHUR_2[i] <- entropart::Hurlbert(as.numeric(mite[i,]),k = 2)
}
entropart_RHUR_3 <- rep(1,70)
for (i in 1:nrow(mite)){
  entropart_RHUR_3[i] <- entropart::Hurlbert(as.numeric(mite[i,]),k = 3)
}
tabula_RHUR_2 <- rep(1,70)
for (i in 1:nrow(mite)){
  tabula_RHUR_2[i] <- tabula::index_hurlbert(as.numeric(mite[i,]),sample = 2)
}
tabula_RHUR_3 <- rep(1,70)
for (i in 1:nrow(mite)){
  tabula_RHUR_3[i] <- tabula::index_hurlbert(as.numeric(mite[i,]),sample = 3)
}
mobsim_RHUR_2 <- rep(1,70)
for (i in 1:nrow(mite)){
  mobsim_RHUR_2[i] <- mobsim::spec_sample(as.numeric(mite[i,]),n=2)
}
mobsim_RHUR_3 <- rep(1,70)
for (i in 1:nrow(mite)){
  mobsim_RHUR_3[i] <- mobsim::spec_sample(as.numeric(mite[i,]),n=3)
}
vegan_RHUR_2 <- rep(1,70)
for (i in 1:nrow(mite)){
  vegan_RHUR_2[i] <- vegan::rarefy(mite[i,],sample = 2)
}
vegan_RHUR_3 <- rep(1,70)
for (i in 1:nrow(mite)){
  vegan_RHUR_3[i] <- vegan::rarefy(mite[i,],sample = 3)
}
divent_RHUR_2 <- rep(1,70)
for (i in 1:nrow(mite)){
  divent_RHUR_2[i] <- divent::div_hurlbert(as.numeric(mite[i,]),k = 2,estimator = "Hurlbert")$diversity
}
divent_RHUR_3 <- rep(1,70)
for (i in 1:nrow(mite)){
  divent_RHUR_3[i] <- divent::div_hurlbert(as.numeric(mite[i,]),k = 3,estimator = "Hurlbert")$diversity
}
# Put all the values in a single dataframe
RHUR <- ls(pattern = "_RHUR")
RHUR <- mget(RHUR)
RHUR <- as.data.frame(RHUR)

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
  tabula_RICE[i] <- wiqid::richICE(t(mite[i,])) #44
}
BAT_RICE <- BAT::alpha.accum(mite)[,22] #NaN
# Put all the values in a single dataframe
RICE <- ls(pattern = "_RICE")
RICE <- mget(RICE)
RICE <- as.data.frame(RICE)

# Jacknife 1 ####
divent_RJA1 <- rep(1,70)
for (i in 1:nrow(mite)){
  divent_RJA1[i] <- divent::div_richness(as.numeric(t(mite[i,])),estimator = "jackknife")$diversity
}
tabula_RJA1 <- tabula::jackknife(tabula::heterogeneity(mite,method = "shannon"),f = summary)[,1]
BiodiversityR_RJA1 <- rep(1,70)
for (i in 1:nrow(mite)){
  BiodiversityR_RJA1[i] <- as.numeric(BiodiversityR::diversityresult(mite[i,],index = "jack1",method = "each site"))
}
entropart_RJA1 <- rep(1,70)
for (i in 1:nrow(mite)){
  entropart_RJA1[i] <- entropart::bcRichness(mite[i,],Correction = "Jackknife")
}
BAT_RJA1 <- BAT::alpha.estimate(mite)[,5]
fossil_RJA1 <- rep(1,70)
for (i in 1:nrow(mite)){
  fossil_RJA1[i] <- fossil::jack1(as.numeric(mite[i,]))
}
wiqid_RJA1 <- rep(1,70)
for (i in 1:nrow(mite)){
  wiqid_RJA1[i] <- wiqid::richJackA1(t(mite[i,])) #13
}
SSP_RJA1 <- rep(1,70)
for (i in 1:nrow(mite)){
  SSP_RJA1[i] <- SSP:::assempar(mite[i,],type = "counts",Sest.method = "jack1")$Sest #10
}

# Put all the values in a single dataframe
RJA1 <- ls(pattern = "_RJA1$")
RJA1 <- mget(RJA1)
RJA1 <- as.data.frame(RJA1)

# Jacknife 2 ####
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

# Put all the values in a single dataframe
RJA2 <- ls(pattern = "_RJA2$")
RJA2 <- mget(RJA2)
RJA2 <- as.data.frame(RJA2)

# Michaelis-Menten ####
wiqid_RMM <- rep(1,70)
for (i in 1:nrow(mite)){
  wiqid_RMM[i] <- wiqid::richMM(t(mite[i,]))
}

# Put all the values in a single dataframe
RMM <- ls(pattern = "_RMM$")
RMM <- mget(RMM)
RMM <- as.data.frame(RMM)

# Squares Richness Estimator ####
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

# Put all the values in a single dataframe
RSQ <- ls(pattern = "_RSQ$")
RSQ <- mget(RSQ)
RSQ <- as.data.frame(RSQ)

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

custom_DGLE <- apply(mite,1,gleason)

# Put all the values in a single dataframe
DGLE <- ls(pattern = "_DGLE$")
DGLE <- mget(DGLE)
DGLE <- as.data.frame(DGLE)

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
custom_DGOO_1_1 <- rep(1,70)
for (i in 1:nrow(mite)){
  custom_DGOO_1_1[i] <- good_index(mite_relat[i, ][mite_relat[i, ] != 0],1,1)
}
custom_DGOO_1_2 <- rep(1,70)
for (i in 1:nrow(mite)){
  custom_DGOO_1_2[i] <- good_index(mite_relat[i, ][mite_relat[i, ] != 0],1,2)
}

custom_DGOO_1_3 <- rep(1,70)
for (i in 1:nrow(mite)){
  custom_DGOO_1_2[i] <- good_index(mite_relat[i, ][mite_relat[i, ] != 0],1,3)
}

custom_DGOO_2_1 <- rep(1,70)
for (i in 1:nrow(mite)){
  custom_DGOO_2_1[i] <- good_index(mite_relat[i, ][mite_relat[i, ] != 0],2,1)
}
custom_DGOO_2_2 <- rep(1,70)
for (i in 1:nrow(mite)){
  custom_DGOO_1_2[i] <- good_index(mite_relat[i, ][mite_relat[i, ] != 0],2,2)
}

custom_DGOO_2_3 <- rep(1,70)
for (i in 1:nrow(mite)){
  custom_DGOO_2_2[i] <- good_index(mite_relat[i, ][mite_relat[i, ] != 0],2,3)
}

custom_DGOO_3_1 <- rep(1,70)
for (i in 1:nrow(mite)){
  custom_DGOO_3_1[i] <- good_index(mite_relat[i, ][mite_relat[i, ] != 0],3,1)
}
custom_DGOO_3_2 <- rep(1,70)
for (i in 1:nrow(mite)){
  custom_DGOO_3_2[i] <- good_index(mite_relat[i, ][mite_relat[i, ] != 0],3,2)
}

custom_DGOO_3_3 <- rep(1,70)
for (i in 1:nrow(mite)){
  custom_DGOO_3_2[i] <- good_index(mite_relat[i, ][mite_relat[i, ] != 0],3,3)
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
custom_DKO <- kothe(mite)

# Put all the values in a single dataframe
DKO <- ls(pattern = "_DKO$")
DKO <- mget(DKO)
DKO <- as.data.frame(DKO)

# Log normal lambda ####
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
  vegan_RLMD[i] <- vegan::veiledspec(vegan::prestonfit(mite[i,]))
}


# Put all the values in a single dataframe
RLMD <- ls(pattern = "_RLMD$")
RLMD <- mget(RLMD)
RLMD <- as.data.frame(RLMD)

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
custom_DOD <- apply(mite,1,odum)

# Put all the values in a single dataframe
DOD <- ls(pattern = "_DOD$")
DOD <- mget(DOD)
DOD <- as.data.frame(DOD)

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
aqp_DSHA <- apply(mite_relat, 1, aqp::shannonEntropy) # FALSE
wiqid_DSHA <- apply(mite,1,wiqid::biodShannon) # FALSE
OTUtable_DSHA <- apply(mite, 1, OTUtable::shannon)
DescTools_DSHA <- apply(mite, 1, DescTools::Entropy) #FALSE
pgirmess_DSHA <- apply(mite, 1, pgirmess::shannon)[1,] #FALSE
divent_DSHA <- apply(mite, 1, function(x) divent::ent_shannon(x, estimator = "naive")$entropy)  # OK
SpiecEasi_DSHA <- apply(mite, 1, SpiecEasi::shannon) #FALSE
abdiv_DSHA <- apply(mite, 1, abdiv::shannon) #FALSE
HardyWeinberg_DSHA <- apply(mite, 1, function(x) HardyWeinberg::shannon(x)$Hp)
MCPAN_DSHA <- apply(mite, 1, function(x) MCPAN::estShannon(x)$estraw)
wavethresh_DSHA <- apply(mite, 1, wavethresh::Shannon.entropy)
codyn_DSHA <- as.numeric(t(as.data.frame(apply(mite, 1, function(x) codyn::community_diversity(
  data.frame(Value = x), abundance.var = "Value", metric = "Shannon")
))))
benthos_DSHA <- rep(1,70)
for (i in 1:nrow(mite)){
  benthos_DSHA[i] <- benthos::shannon(taxon = names(mite_relat[i,]), count = mite_relat[i,])
}
triversity_DSHA <- rep(1,70) 
for (i in 1:nrow(mite)){
  triversity_DSHA[i] <- triversity::get_diversity_from_distribution(mite_relat[i,])[6]
} #FALSE
forestmangr_DSHA <- rep(1,70) 
for (i in 1:nrow(mite)){
  mite_long <- tidyr::pivot_longer(mite[i,],cols = colnames(mite)[1]:colnames(mite)[35],names_to = "Esp",values_to = "Value")
  mite_long <- mite_long[1:35,]
  forestmangr_DSHA[i] <- forestmangr::species_diversity(mite_long,species = "Esp",index = "H")
}

entropart_DSHA <- rep(1,70) 
for (i in 1:nrow(mite)){
  entropart_DSHA[i] <- entropart::Shannon(as.numeric(mite[i,]),Correction = "None")
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
agricolae_DSP <- rep(1,70)
for (i in seq_len(nrow(mite))) {
  agricolae_DSP[i] <- tryCatch(
    {agricolae::index.bio(mite[i,],method = "Simpson.Dom")$index},
    error = function(e) NA
  )
}

abdiv_DSP <- rep(1,70)
for (i in 1:nrow(mite)){
  abdiv_DSP[i] <- abdiv::dominance(mite[i,])
}
microbiome_DSP <- rep(1,70)
for (i in 1:nrow(mite)){
  microbiome_DSP[i] <- microbiome::dominance(as.numeric(mite[i,]),index = "simpson")$simpson
}
EconGeo_DSP <- EconGeo::herfindahl(mite) 
concstats_DSP <- rep(1,70)
for (i in 1:nrow(mite)){
  concstats_DSP[i] <- concstats::concstats_hhi(as.numeric(mite[i,])) 
}
REAT_DSP <- rep(1,70)
for (i in 1:nrow(mite)){
  REAT_DSP[i] <- REAT::herf(mite[i,]) 
}
hhi_DSP <- rep(1,70)
for (i in 1:nrow(mite)){
  mite_long <- tidyr::pivot_longer(mite[i,],cols = colnames(mite)[1]:colnames(mite)[35],names_to = "Esp",values_to = "Value")
  mite_long <- mite_long[1:35,]
  mite_long$Relative <- as.numeric(mite_long$Value/sum(mite_long$Value))
  hhi_DSP[i] <- hhi::hhi(as.data.frame(mite_long),"Relative") 
}

politics_DSP <- rep(1,70)
for (i in 1:nrow(mite)){
  politics_DSP[i] <- politicsR::hh(as.numeric(mite_relat[i,])) 
}

PDtoolkit_DSP <- rep(1,70)
for (i in 1:nrow(mite)){
  PDtoolkit_DSP[i] <- PDtoolkit::hhi(as.numeric(mite_relat[i,]))
}
antitrust_DSP <- rep(1,70)
for (i in 1:nrow(mite)){
  antitrust_DSP[i] <- antitrust::HHI(as.numeric(mite_relat[i,])) # Faux, multiplicateur
}
DescTools_DSP <- rep(1,70)
for (i in 1:nrow(mite)){
  DescTools_DSP[i] <- DescTools::Herfindahl(as.numeric(mite[i,])) # OK
}
divseg_DSP <- divseg::ds_hhi(mite,.cols = dplyr::everything()) # OK
ineq_DSP <- rep(1,70)
for (i in 1:nrow(mite)){
  ineq_DSP[i] <- ineq::Herfindahl(mite[i,]) # OK
}

triversity_DSP <- rep(1,70)
for (i in 1:nrow(mite)){
  triversity_DSP[i] <- triversity::get_diversity_from_distribution(as.numeric(mite_relat[i,]),measure = "herfindahl" ) # OK
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
MCPAN_D1SP <- apply(mite_relat,1,MCPAN::Simpson)
entropart_D1SP <- rep(1,70)
for (i in 1:nrow(mite)){
  entropart_D1SP[i] <- entropart::Simpson(as.numeric(mite[i,]),Correction = "None")
}
BiodiversityR_D1SP <- BiodiversityR::diversityresult(mite,y=NULL,index="Simpson",method = "each site")$Simpson
simboot_D1SP <- rep(1,70)
for (i in 1:nrow(mite)){
  simboot_D1SP[i] <- simboot::Simpson(mite_relat[i,])
}
concstats_D1SP <- rep(1,70)
for (i in 1:nrow(mite)){
  concstats_D1SP[i] <- concstats::concstats_simpson(as.numeric(mite_relat[i,]),na.rm = T)
}
untb_D1SP <- rep(1,70)
for (i in 1:nrow(mite)){
  untb_D1SP[i] <- untb::simpson(mite[i,],with.replacement = F)#
}
adiv_D1SP <- adiv::speciesdiv(mite,method = "GiniSimpson")[,1]
diverse_D1SP <- diverse::diversity(t(mite),type = "gini-simpson",category_row = T)[[1]]
asbio_D1SP <- asbio::alpha.div(mite,"simp")

agricolae_D1SP <- rep(1,70)
for (i in seq_len(nrow(mite))) {
  agricolae_D1SP[i] <- tryCatch(
    {agricolae::index.bio(mite[i,],method = "Simpson.Div")$index},
    error = function(e) NA
  )
}
microbiome_D1SP <- microbiome::diversity(t(mite), index = "gini_simpson")$gini_simpson
divent_D1SP <- rep(1,70)
for (i in 1:nrow(mite)){
  divent_D1SP[i] <- divent::ent_simpson(as.numeric(mite[i,]),estimator = "naive")$entropy
}
lawstat_D1SP <- rep(1,70)
for (i in 1:nrow(mite)){
  lawstat_D1SP[i] <- as.numeric(lawstat::gini.index(as.numeric(mite[i,])))
}
CUB_D1SP <- rep(1,70)
for (i in 1:nrow(mite)){
  CUB_D1SP[i] <- CUB::gini(as.numeric(mite_relat[i,]))
}
RoughSets_D1SP <- rep(1,70)
for (i in 1:nrow(mite)){
  RoughSets_D1SP[i] <- RoughSets::X.gini(as.numeric(mite[i,])) #OK
}
breakaway_D1SP <- rep(1,70)
for (i in 1:nrow(mite)){
  breakaway_D1SP[i] <- breakaway::true_gini(mite_relat[i,]) #OK
}
catsim_D1SP <- rep(1,70)
for (i in 1:nrow(mite)){
  catsim_D1SP[i] <- catsim::gini(as.numeric(mite[i,])) # faux
}
ade4_D1SP <- ade4::divc(as.data.frame(t(mite)),dis = NULL)$diversity
DescTools_D1SP <- DescTools::DivCoef(as.data.frame(t(mite)),dis = NULL)$diversity # 115
PCRA_D1SP <- as.numeric(PCRA::divHHI(mite_relat)) # 1-HHI = 1-Simpson = Gini-Simpson

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
wiqid_DINVSP <- rep(1,70)
for (i in 1:nrow(mite)){
  wiqid_DINVSP[i] <- wiqid::biodSimpson(abVec=as.numeric(mite[i,]),correct = F)
}#FALSE
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
  entropart_QHIL_0[i] <- entropart::Diversity(as.numeric(mite[i,]),q=0,Correction = "None")
}
entropart_QHIL_1 <- rep(1,70) 
for (i in 1:nrow(mite)){
  entropart_QHIL_1[i] <- entropart::Diversity(as.numeric(mite[i,]),q=1,Correction = "None")
}
entropart_QHIL_2 <- rep(1,70) 
for (i in 1:nrow(mite)){
  entropart_QHIL_2[i] <- entropart::Diversity(as.numeric(mite[i,]),q=2,Correction = "None")
}
entropart_QHIL_3 <- rep(1,70) 
for (i in 1:nrow(mite)){
  entropart_QHIL_3[i] <- entropart::Diversity(as.numeric(mite[i,]),q=3,Correction = "None")
}

divent_QHIL_0 <- rep(1,70) 
for (i in 1:nrow(mite)){
  divent_QHIL_0[i] <- divent::div_hill(as.numeric(mite[i,]),q = 0,probability_estimator = "naive")$diversity
}
divent_QHIL_1 <- rep(1,70) 
for (i in 1:nrow(mite)){
  divent_QHIL_1[i] <- divent::div_hill(as.numeric(mite[i,]),q = 1,probability_estimator = "naive")$diversity
}
divent_QHIL_2 <- rep(1,70) 
for (i in 1:nrow(mite)){
  divent_QHIL_2[i] <- divent::div_hill(as.numeric(mite[i,]),q = 2,probability_estimator = "naive")$diversity
}
divent_QHIL_3 <- rep(1,70) 
for (i in 1:nrow(mite)){
  divent_QHIL_3[i] <- divent::div_hill(as.numeric(mite[i,]),q = 3,probability_estimator = "naive")$diversity
}

# Put all the values in a single dataframe
QHIL  <- ls(pattern = "_QHIL")
QHIL <- mget(QHIL)
QHIL <- as.data.frame(QHIL)

# Tsallis ####
vegan_QTSA_0 <- vegan::tsallis(mite,scales = 0,hill = F)
vegan_QTSA_1 <- vegan::tsallis(mite,scales = 1,hill = F)
vegan_QTSA_2 <- vegan::tsallis(mite,scales = 2,hill = F)
vegan_QTSA_3 <- vegan::tsallis(mite,scales = 3,hill = F)

EntropyEstimation_QTSA_1 <- rep(1,70) 
for (i in 1:nrow(mite)){
  EntropyEstimation_QTSA_1[i] <- EntropyEstimation::Tsallis.z(mite[i,],1)
}
EntropyEstimation_QTSA_2 <- rep(1,70) 
for (i in 1:nrow(mite)){
  EntropyEstimation_QTSA_2[i] <- EntropyEstimation::Tsallis.z(mite[i,],2)
}
EntropyEstimation_QTSA_3 <- rep(1,70) 
for (i in 1:nrow(mite)){
  EntropyEstimation_QTSA_3[i] <- EntropyEstimation::Tsallis.z(mite[i,],3)
}

adiv_QTSA_0 <- adiv::divparam(mite,method = "tsallis",q=0)
adiv_QTSA_1 <- adiv::divparam(mite,method = "tsallis",q=1)
adiv_QTSA_2 <- adiv::divparam(mite,method = "tsallis",q=2)
adiv_QTSA_3 <- adiv::divparam(mite,method = "tsallis",q=3)

entropart_QTSA_0 <- rep(1,70) 
for (i in 1:nrow(mite)){
  entropart_QTSA_0[i] <- entropart::Tsallis(as.numeric(mite[i,]),q=0,Correction = "None")
}
entropart_QTSA_1 <- rep(1,70) 
for (i in 1:nrow(mite)){
  entropart_QTSA_1[i] <- entropart::Tsallis(as.numeric(mite[i,]),q=1,Correction = "None")
}
entropart_QTSA_2 <- rep(1,70) 
for (i in 1:nrow(mite)){
  entropart_QTSA_2[i] <- entropart::Tsallis(as.numeric(mite[i,]),q=2,Correction = "None")
}
entropart_QTSA_3 <- rep(1,70) 
for (i in 1:nrow(mite)){
  entropart_QTSA_3[i] <- entropart::Tsallis(as.numeric(mite[i,]),q=3,Correction = "None")
}

divent_QTSA_0 <- rep(1,70) 
for (i in 1:nrow(mite)){
  divent_QTSA_0[i] <- divent::ent_tsallis(as.numeric(mite[i,]),q=0,probability_estimator = "naive",richness_estimator = "naive")$entropy # OK
}
divent_QTSA_1 <- rep(1,70) 
for (i in 1:nrow(mite)){
  divent_QTSA_1[i] <- divent::ent_tsallis(as.numeric(mite[i,]),q=1,probability_estimator = "naive",richness_estimator = "naive")$entropy # OK
}
divent_QTSA_2 <- rep(1,70) 
for (i in 1:nrow(mite)){
  divent_QTSA_2[i] <- divent::ent_tsallis(as.numeric(mite[i,]),q=2,probability_estimator = "naive",richness_estimator = "naive")$entropy # OK
}
divent_QTSA_3 <- rep(1,70) 
for (i in 1:nrow(mite)){
  divent_QTSA_3[i] <- divent::ent_tsallis(as.numeric(mite[i,]),q=3,probability_estimator = "naive",richness_estimator = "naive")$entropy # OK
}

# Put all the values in a single dataframe
QTSA <- ls(pattern = "_QTSA")
QTSA <- mget(QTSA)
QTSA <- as.data.frame(QTSA)

# Rényi ####
EntropyEstimation_QREN_1 <- rep(1,70) 
for (i in 1:nrow(mite)){
  EntropyEstimation_QREN_1[i] <- EntropyEstimation::Renyi.z(as.numeric(mite[i,]),r=1)
}
EntropyEstimation_QREN_2 <- rep(1,70) 
for (i in 1:nrow(mite)){
  EntropyEstimation_QREN_2[i] <- EntropyEstimation::Renyi.z(as.numeric(mite[i,]),r=2)
}
EntropyEstimation_QREN_3 <- rep(1,70) 
for (i in 1:nrow(mite)){
  EntropyEstimation_QREN_3[i] <- EntropyEstimation::Renyi.z(as.numeric(mite[i,]),r=3)
}
sprex_QREN_0 <- rep(1,70) 
for (i in 1:nrow(mite)){
  sprex_QREN_0[i] <- sprex::diversity(as.numeric(mite[i,]),type = "renyi",q=0) #OK
}
sprex_QREN_1 <- rep(1,70) 
for (i in 1:nrow(mite)){
  sprex_QREN_1[i] <- sprex::diversity(as.numeric(mite[i,]),type = "renyi",q=1) #OK
}
sprex_QREN_2 <- rep(1,70) 
for (i in 1:nrow(mite)){
  sprex_QREN_2[i] <- sprex::diversity(as.numeric(mite[i,]),type = "renyi",q=2) #OK
}
sprex_QREN_3 <- rep(1,70) 
for (i in 1:nrow(mite)){
  sprex_QREN_3[i] <- sprex::diversity(as.numeric(mite[i,]),type = "renyi",q=3) #OK
}
statcomp_QREN_1 <- rep(1,70) 
for (i in 1:nrow(mite)){
  statcomp_QREN_1[i] <- statcomp::permutation_entropy_Renyi(mite[i,],1) #FALSE
}
statcomp_QREN_2 <- rep(1,70) 
for (i in 1:nrow(mite)){
  statcomp_QREN_2[i] <- statcomp::permutation_entropy_Renyi(mite[i,],2) #FALSE
}
statcomp_QREN_3 <- rep(1,70) 
for (i in 1:nrow(mite)){
  statcomp_QREN_3[i] <- statcomp::permutation_entropy_Renyi(mite[i,],3) #FALSE
}
vegan_QREN_0 <- vegan::renyi(mite,scales = 0,hill = F)
vegan_QREN_1 <- vegan::renyi(mite,scales = 1,hill = F)
vegan_QREN_2 <- vegan::renyi(mite,scales = 2,hill = F)
vegan_QREN_3 <- vegan::renyi(mite,scales = 3,hill = F)
seewave_QREN_2 <- rep(1,70) 
for (i in 1:nrow(mite)){
  seewave_QREN_2[i] <- seewave::sh(mite[i,],alpha = 2)
}
seewave_QREN_3 <- rep(1,70) 
for (i in 1:nrow(mite)){
  seewave_QREN_3[i] <- seewave::sh(mite[i,],alpha = 3)
}
BiodiversityR_QREN_0 <- as.numeric(BiodiversityR::renyiresult(mite,y=NULL,method = "each site",scales = 0))
BiodiversityR_QREN_1 <- as.numeric(BiodiversityR::renyiresult(mite,y=NULL,method = "each site",scales = 1))
BiodiversityR_QREN_2 <- as.numeric(BiodiversityR::renyiresult(mite,y=NULL,method = "each site",scales = 2))
BiodiversityR_QREN_3 <- as.numeric(BiodiversityR::renyiresult(mite,y=NULL,method = "each site",scales = 3))

adiv_QREN_0 <- adiv::divparam(mite,method = "renyi",q=0)
adiv_QREN_1 <- adiv::divparam(mite,method = "renyi",q=1)
adiv_QREN_2 <- adiv::divparam(mite,method = "renyi",q=2)
adiv_QREN_3 <- adiv::divparam(mite,method = "renyi",q=3)

diverse_QREN_0 <- diverse::diversity(t(mite),type = "renyi",category_row = T,q=0)$renyi.entropy
diverse_QREN_1 <- diverse::diversity(t(mite),type = "renyi",category_row = T,q=1)$renyi.entropy
diverse_QREN_2 <- diverse::diversity(t(mite),type = "renyi",category_row = T,q=2)$renyi.entropy
diverse_QREN_3 <- diverse::diversity(t(mite),type = "renyi",category_row = T,q=3)$renyi.entropy

# Put all the values in a single dataframe
QREN <- ls(pattern = "_QREN")
QREN <- mget(QREN)
QREN <- as.data.frame(QREN)

# Berger-Parker ####
abdiv_EBP <- apply(mite,1,abdiv::berger_parker_d)
dyvDyn_EBP <- rep(1,70) 
for (i in 1:nrow(mite)){
  dyvDyn_EBP[i] <- divDyn::indices(as.matrix(mite[i,]),method = "dominance")
}
tabula_EBP <- rep(1,70) 
for (i in 1:nrow(mite)){
  tabula_EBP[i] <- tabula::index_berger(as.numeric(mite[i,]))
}
BiodiversityR_EBP <- BiodiversityR::diversityresult(mite,y=NULL,index="Berger",method = "each site")$Berger
diverse_EBP <- diverse::diversity(t(mite),type = "berger-parker",category_row = T)$berger.parker.D
agricolae_EBP <- rep(1,70) 
for (i in seq_len(nrow(mite))) {
  agricolae_EBP[i] <- tryCatch(
    {
      agricolae::index.bio(mite[i,],method = "Berger.Parker")$index
    },
    error = function(e) NA
  )
}
ecodive_EBP <- ecodive::alpha_div(mite,metric="berger")

wiqid_EBP <- rep(1,70) 
for (i in 1:nrow(mite)){
  wiqid_EBP[i] <- wiqid::biodBerger(abVec=mite[i,]) # FALSE
}
microbiome_EBP <- rep(1,70) 
for (i in 1:nrow(mite)){
  microbiome_EBP[i] <- microbiome::dominance(as.numeric(mite[i,]),index = "DBP")$dbp # Berger-Parker
}
triversity_EBP <- rep(1,70) 
for (i in 1:nrow(mite)){
  triversity_EBP[i] <- triversity::get_diversity_from_distribution(as.numeric(mite_relat[i,]),measure = "bergerparker" ) # OK
}

# Put all the values in a single dataframe
EBP <- ls(pattern = "_EBP")
EBP <- mget(EBP)
EBP <- as.data.frame(EBP)

# Brillouin ####
tabula_EBRI <- as.numeric(tabula::evenness(mite,method = "brillouin"))

# Put all the values in a single dataframe
EBRI <- ls(pattern = "_EBRI")
EBRI <- mget(EBRI)
EBRI <- as.data.frame(EBRI)

# Bulla ####
BAT_EBU <- as.numeric(BAT::evenness(mite,func = "bulla"))
microbiome_EBU <- microbiome::evenness(t(mite), index = "bulla")[,1] # OK

# Put all the values in a single dataframe
EBU <- ls(pattern = "_EBU")
EBU <- mget(EBU)
EBU <- as.data.frame(EBU)

# Camargo ####
BAT_ECAM <- as.numeric(BAT::evenness(mite,func = "camargo"))
microbiome_ECAM <- microbiome::evenness(t(mite), index = "camargo")[,1] # OK

# Put all the values in a single dataframe
ECAM <- ls(pattern = "_ECAM")
ECAM <- mget(ECAM)
ECAM <- as.data.frame(ECAM)

# Gini ####
DescTools_EGIN <- rep(1,70) 
for (i in 1:nrow(mite)){
  DescTools_EGIN[i] <- DescTools::Gini(as.numeric(mite[i,]),unbiased = F) # Gini
}
microbiome_EGIN <- rep(1,70) 
for (i in 1:nrow(mite)){
  microbiome_EGIN[i] <- microbiome::dominance(as.numeric(mite[i,]),index = "gini")$gini
}
giniVarCI_EGIN <- rep(1,70) 
for (i in 1:nrow(mite)){
  giniVarCI_EGIN[i] <- giniVarCI::igini(as.numeric(mite[i,]),bias.correction = F)
}
acid_EGIN <- rep(1,70) 
for (i in 1:nrow(mite)){
  acid_EGIN[i] <- acid::gini(as.numeric(mite[i,]))$Gini 
}
dplR_EGIN <- rep(1,70) 
for (i in 1:nrow(mite)){
  dplR_EGIN[i] <- dplR::gini.coef(as.numeric(mite[i,]))  
}
shipunov_EGIN <- rep(1,70) 
for (i in 1:nrow(mite)){
  shipunov_EGIN[i] <- shipunov::Gini(as.numeric(mite[i,]))  
}
EconGeo_EGIN <- EconGeo::gini(t(mite))$Gini
ineq_EGIN <- rep(1,70) 
for (i in 1:nrow(mite)){
  ineq_EGIN[i] <- ineq::Gini(mite[i,]) 
}
wINEQ_EGIN <- rep(1,70) 
for (i in 1:nrow(mite)){
  wINEQ_EGIN[i] <- wINEQ::Gini(as.numeric(mite[i,]))
}
REAT_EGIN <- rep(1,70) 
for (i in 1:nrow(mite)){
  REAT_EGIN[i] <- REAT::gini(as.numeric(mite[i,]))
}

# Put all the values in a single dataframe
EGIN <- ls(pattern = "_EGIN")
EGIN <- mget(EGIN)
EGIN <- as.data.frame(EGIN)


#Heip ####
adiv_EHEI <- adiv::specieseve(mite,method = "Heip")[,1] # OK
abdiv_EHEI <- rep(1,70) 
for (i in 1:nrow(mite)){
  abdiv_EHEI[i] <- abdiv::heip_e(as.numeric(mite[i,])) # OK
}

# Put all the values in a single dataframe
EHEI <- ls(pattern = "_EHEI")
EHEI <- mget(EHEI)
EHEI <- as.data.frame(EHEI)

# Hulburt ####
hulburt_index <- function(x) {
  s <- sum(x, na.rm = TRUE)
  if (s == 0) return(NA_real_)
  sum(sort(x, decreasing = TRUE)[1:2], na.rm = TRUE) / s
}
custom_EHUL <- rep(1,70) 
for (i in 1:nrow(mite)){
  custom_EHUL[i] <- hulburt_index(as.numeric(mite[i,]))
}

# Put all the values in a single dataframe
EHUL <- ls(pattern = "_EHUL")
EHUL <- mget(EHUL)
EHUL <- as.data.frame(EHUL)

# Hurlbert E ####
# Non calculable car besoin de n : nombre d'individus dans un echantillon standard

# Hurlbert PIE ####
mobr_EHURE <- mobr::calc_div(mite,index="PIE",effort = NA) # Hulbert PIE
vegan_EHURE <- vegan::simpson.unb(mite) # Hulbert PIE 
BiodiversityR_EHURE <- BiodiversityR::diversityresult(mite,y=NULL,index="simpson.unb",method = "each site")[,1] # Hulbert PIE = unbiaised Simpson

benthos_EHURE <- rep(1,70) 
for (i in 1:nrow(mite)){
  benthos_EHURE[i] <- benthos::hpie(taxon = colnames(mite),count = as.integer(mite[i,])) # Hulbert PIE
}

# Put all the values in a single dataframe
EHURE <- ls(pattern = "_EHURE")
EHURE <- mget(EHURE)
EHURE <- as.data.frame(EHURE)

# Ludwig-Reynolds ####
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
custom_ELR <- rep(1,70) 
for (i in 1:nrow(mite)){
  custom_ELR[i] <- Ludwig_Reynold(mite[i,])
}

# Put all the values in a single dataframe
ELR <- ls(pattern = "_ELR")
ELR <- mget(ELR)
ELR <- as.data.frame(ELR)

# McIntosh ####
tabula_EMC <- as.numeric(tabula::evenness(mite,method = "mcintosh"))
abdiv_EMC <- apply(mite,1,abdiv::mcintosh_e)
# Put all the values in a single dataframe
EMC <- ls(pattern = "_EMC")
EMC <- mget(EMC)
EMC <- as.data.frame(EMC)

# McNaughton ####
microbiome_EMN <- rep(1,70) 
for (i in 1:nrow(mite)){
  microbiome_EMN[i] <- microbiome::dominance(as.numeric(mite[i,]),index = "DMN")$dmn  #OK
}
# Put all the values in a single dataframe
EMN <- ls(pattern = "_EMN")
EMN <- mget(EMN)
EMN <- as.data.frame(EMN)

# NHC : Nee-Harvey-Cotgreave 
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
custom_ENHC <- rep(1,70) 
for (i in 1:nrow(mite)){
  custom_ENHC[i] <- nhc_original(mite[i,])
}
# Put all the values in a single dataframe
ENHC <- ls(pattern = "_ENHC")
ENHC <- mget(ENHC)
ENHC <- as.data.frame(ENHC)


# NHC EQ ####
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
custom_EEQ <- rep(1,70) 
for (i in 1:nrow(mite)){
  custom_EEQ[i] <- e0_index(mite[i,]) 
}
codyn_EEQ <- rep(1,70) 
for (i in 1:nrow(mite)){
  mite_long <- tidyr::pivot_longer(mite[i,],cols = colnames(mite)[1]:colnames(mite)[35],names_to = "Esp",values_to = "Value")
  mite_long <- mite_long[1:35,]
  codyn_EEQ[i] <- codyn::community_structure(mite_long,abundance.var = "Value",metric = "EQ")$EQ
}
# Put all the values in a single dataframe
EEQ <- ls(pattern = "_EEQ")
EEQ <- mget(EEQ)
EEQ <- as.data.frame(EEQ)

# NHC Evar ####
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

# Put all the values in a single dataframe
EEVAR <- ls(pattern = "_EEVAR")
EEVAR <- mget(EEVAR)
EEVAR <- as.data.frame(EEVAR)

# Patten ####
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
custom_EPAT <- patten_index(mite)

# Put all the values in a single dataframe
EPAT <- ls(pattern = "_EPAT")
EPAT <- mget(EPAT)
EPAT <- as.data.frame(EPAT)

# Pielou ####
#### Pielou ####----------------------------------------------------------------
OTUtable_EPIE <- rep(1,70) 
for (i in 1:nrow(mite)){
  OTUtable_EPIE[i] <- OTUtable::pielou(mite[i,]) 
}
abdiv_EPIE <- rep(1,70) 
for (i in 1:nrow(mite)){
  abdiv_EPIE[i] <- abdiv::pielou_e(as.numeric(mite_relat[i,])) 
}
tabula_EPIE <- tabula::evenness(mite,method = "shannon")@.Data
BiodiversityR_EPIE <- BiodiversityR::diversityresult(mite,y=NULL,index="Jevenness",method = "each site")$Jevenness
forestmangr_EPIE <- rep(1,70) 
for (i in 1:nrow(mite)){
  mite_long <- tidyr::pivot_longer(mite[i,],cols = colnames(mite)[1]:colnames(mite)[35],names_to = "Esp",values_to = "Value")
  mite_long <- mite_long[1:35,]
  forestmangr_EPIE[i] <- forestmangr::species_diversity(mite_long,species = "Esp",index = "S")
}

adiv_EPIE <- adiv::specieseve(mite,method = "Shannon")[,1]
pctax_EPIE <- pctax::a_diversity(t(mite),method = "pielou")$Pielou_evenness
diverse_EPIE <- diverse::diversity(t(mite),type = "evenness",category_row = T)[,1]
microbiome_EPIE <- microbiome::evenness(t(mite), index = "pielou")[,1]
chemodiv_EPIE <- chemodiv::calcDiv(mite,type = "PielouEven")[,1] # OK
sprex_EPIE <- rep(1,70) 
for (i in 1:nrow(mite)){
  sprex_EPIE[i] <- sprex::diversity(as.numeric(mite[i,]),type = "eveness.pielou")
}
breakaway_EPIE <- rep(1,70) 
for (i in 1:nrow(mite)){
  breakaway_EPIE[i] <- breakaway::true_shannon_e(mite_relat[i,]) # Pielou
}

# Put all the values in a single dataframe
EPIE <- ls(pattern = "_EPIE")
EPIE <- mget(EPIE)
EPIE <- as.data.frame(EPIE)

# Sheldon ####
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
custom_ESHE <- rep(1,70) 
for (i in 1:nrow(mite)){
  custom_ESHE[i] <- sheldon_index(mite[i,])
}
# Put all the values in a single dataframe
ESHE <- ls(pattern = "_ESHE")
ESHE <- mget(ESHE)
ESHE <- as.data.frame(ESHE)

# Inv Simpson Evenness ####
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

# Put all the values in a single dataframe
ESP <- ls(pattern = "_ESP")
ESP <- mget(ESP)
ESP <- as.data.frame(ESP)

# Strong ####
abdiv_EST <- rep(1,70) 
for (i in 1:nrow(mite)){
  abdiv_EST [i] <- abdiv::strong(as.numeric(mite[i,]))
}
# Put all the values in a single dataframe
EST <- ls(pattern = "_EST")
EST <- mget(EST)
EST <- as.data.frame(EST)

# Indicate species richness for next analysis
custom_RICH <- rep(1,70) 
for (i in 1:nrow(mite)){
  custom_RICH [i] <- sum(mite[i,] > 0)
}
custom_RICH <- as.data.frame(custom_RICH)

# Save all the dataset one by one ####
df_exclude <- c("mite","mite_long","mite_dist","mite_relat","mat_need")
alpha_data <- Filter(is.data.frame,
  mget(setdiff(ls(.GlobalEnv), df_exclude), envir = .GlobalEnv)
)

Map(function(df, nom) write.csv(df, paste0("data/alpha/",paste0("a",nom,"_all"), ".csv"), row.names = FALSE),
    alpha_data,
    names(alpha_data))

# Combine all the dataset in only one
alpha_combined <- do.call(cbind, alpha_data)
colnames(alpha_combined) <- sub(".*\\.", "", colnames(alpha_combined))

to_numeric_df <- function(df) {
  df[] <- lapply(df, function(col) as.numeric(as.character(col)))
  return(df)
}
alpha_combined <- to_numeric_df(alpha_combined)
write.csv(alpha_combined,"data/alpha/a_combined_all.csv",row.names = F)

# Same but with a long version
alpha_combined$Sample <- c(1:nrow(alpha_combined))
alpha_combined_long <- pivot_longer(alpha_combined,cols = colnames(alpha_combined)[1:358],names_to = "package_index",values_to = "value")
alpha_combined_long <- alpha_combined_long |>
  separate(package_index, into = c("package", "index"), sep = "_",extra = "merge")
write.csv(alpha_combined_long,"data/alpha/a_combined_long_all.csv",row.names = F)

################################################################################

# Manipulate datasets
correspondance_codes_names <- read_delim("data/correspondance_codes_names.csv", 
                                         delim = ";", escape_double = FALSE, trim_ws = TRUE)

data <- read_csv("data/alpha/a_combined_long_all.csv")

data <- left_join(data,correspondance_codes_names)
write.csv(data,"data/alpha/a_combined_long_all_withnames.csv",row.names = F)

################################################################################

# Prepare data to have only good values per index ####
data <- read_csv("data/alpha/a_combined_long_all.csv")
data$selection <- paste0(data$index,"_",data$package)

dataindex <- filter(data,selection %in% c("D1SP_vegan","DBRI_abdiv","DFIS_vegan",
                                      "DINVSP_vegan","DMC_abdiv","DMG_abdiv",
                                      "DMN_abdiv","DSHA_vegan",
                                      "DSP_abdiv","EBP_abdiv","EBRI_tabula","EBU_BAT",
                                      "EEVAR_adiv","EGIN_microbiome","EHEI_abdiv",
                                      "EHURE_vegan","EMC_tabula","EMN_microbiome",
                                      "EPIE_abdiv","ESP_tabula","EST_abdiv",
                                      "QHIL_0_vegan","QHIL_1_vegan","QHIL_2_vegan","QHIL_3_vegan",
                                      "QREN_0_vegan","QREN_1_vegan","QREN_2_vegan","QREN_3_vegan",
                                      "QTSA_0_vegan","QTSA_1_vegan","QTSA_2_vegan","QTSA_3_vegan",
                                      "RACE_vegan","RC1_vegan","RC1M_entropart","RCA2_tabula",
                                      "RHUR_2_vegan","RHUR_3_vegan","RICE_tabula","RJA1_entropart",
                                      "RJA2_wiqid","RSQ_tabula","RICH_custom"))
dataindex <- dataindex |>
  select(-c(package,selection)) 

dataindex_wide <- dataindex|>
  pivot_wider(names_from = index,values_from = value)

write.csv(dataindex,"data/alpha/alpha_values_mite_long.csv",row.names = F)
write.csv(dataindex_wide,"data/alpha/alpha_values_mite_wide.csv",row.names = F)

