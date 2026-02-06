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

mite_pa <- convert_to_presence_absence(mite)

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



#### Dissimilarity ratio ####---------------------------------------------------
vegan::designdist(x12,method = "(A+B-2*J)/(A+B-J)",terms = "quadratic") # 0.6569256

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

#### Jeffreys ####

EnvNJ::metrics(t(x12),method = "jeffreys") #2.341993

#### Wishart ####---------------------------------------------------------------
adespatial::beta.div(x12,method = "wishart",save.D = T)$D #0.6569256

print(adespatial::dist.ldc(x12,method = "wishart"))#0.6569256

wiqid::distSimRatio(x1,x2) # Similarity ratio 0.6569256

#### Mean character difference ####
?abdiv::mean_character_difference(as.numeric(x1),as.numeric(x2)) #2.160455

#### Modified mean character difference ####
print(adespatial::dist.ldc(x12,method = "modmeanchardiff")) #2.880606

abdiv::modified_mean_character_difference(as.numeric(x1),as.numeric(x2)) #2.880606

#### S2 coeff Gower & Legendre ####
ade4::dist.binary(x12,method = 10) #0.7071068

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

#### Roberts ####
?labdsv::dsvdis(x11, index = "roberts") #0.6496285

### Kullback Leibler 
Rfast::Dist(x12,method = "kullback_leibler") #172.0308 # non coherent avec formule mais correspond avec Jeffreys
proxyC::dist(as.matrix(x12),method = "kullback") #NA
EnvNJ::metrics(t(x12),method = "kullback_leibler") #marche pas
bapred::kldist(x12) #marche pas

#### Harvesine
Rfast::Dist(x12_pa,method = "haversine") #non fonctionnel

#### Harmonic mean
Rfast::Dist(x12,method = "harmonic_mean") #9.510365

