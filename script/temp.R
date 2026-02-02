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
