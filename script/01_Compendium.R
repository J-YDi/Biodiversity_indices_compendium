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

#___________________________ Alpha diversity indices ___________________________####

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


# A regler
codyn_DSHA <-apply(mite, 1, function(x) codyn::community_diversity(
  data.frame(Value = x), abundance.var = "Value", metric = "Shannon")
)  
apply(mite, 1, function(x) benthos::shannon(taxon = names(x), count = x))        # ne correspond pas à l'indice
apply(mite, 1, breakaway::sample_shannon)
apply(mite, 1, function(x) triversity::get_diversity_from_distribution(
  x, measure = "entropy")
)  

# Marche pas : 
apply(mite, 1, function(x) divseg::ds_shannon(
  data.frame(x), .cols = dplyr::everything())
)                                                       # marche pas
apply(mite, 1, MCPAN::Shannon)                          # marche pas
apply(mite, 1, function(x) spatialEco::shannons(
  x, counts = TRUE, ens = FALSE, margin = "row")
)                                                       # marche pas
apply(mite, 1, coda4microbiome::shannon)                # a une dépendance mal conçue
apply(mite, 1, function(x) triversity::get_diversity_from_distribution(
  x, measure = "entropy")
) 










