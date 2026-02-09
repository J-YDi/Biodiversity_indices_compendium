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
adespatial::beta.div(x12,method = "ab.sorensen",save.D = T,sqrt.D = F)$D #0.02570608

adespatial::beta.div.comp(x12,coef = "S",quant = T)$D #0.5310021

BAT::beta(x11,func = "sorensen",abund = T)$Btotal #0.5310021

wiqid::distChaoSorCorr(x1,x2) #0.02570608
wiqid::distChaoSorNaive(x1,x2) #0.02642424

CommEcol::dis.chao(x11,index = "sorensen",version = "rare") #0.02570608
CommEcol::dis.chao(x12,index = "sorensen",version = "probability") #0.02642424

adespatial::beta.div(x12,method = "sorensen",save.D = T,sqrt.D = F)$D #0.6416889 ne correspond pas

philentropy::sorensen(as.numeric(x1),as.numeric(x1)) #0.5310021

EnvNJ::metrics(t(x11),method = "sorensen") #0.5310021

spaa::sp.pair(t(as.matrix(x12)))$Dice #0.833333

NST::beta.g(x12,dist.method = "chao.sorensen") #0.0244926

prabclus::dicedist(t(x12)) #0.4375 ???



#### Koleff betasor Simpson dissimilarity turnover sorensen ####----------------

adespatial::beta.div(x12,method = "ab.simpson",save.D = T)$D #0.008127366 correspond a rien
print(adespatial::dist.ldc(x12,method = "ab.simpson")) #0.008127366 correspond a rien

adespatial::beta.div.comp(x12,coef = "J",quant = F)$repl # Podani 0.2424242
adespatial::beta.div.comp(x12,coef = "J",quant =T)$repl # Podani 0.6891419

adiv::betastatsor(x12)[2] # Ricotta & Pavoine  0.07272727



#### Nestedness-resultant component of Sørensen dissimilarity ####--------------

adespatial::beta.div.comp(x12_pa,coef = "S",quant = F)$rich #Podani 0.05454545
adespatial::beta.div.comp(x12_pa,coef = "S",quant = T)$rich #Podani 0.0034633

BAT::beta(x12,func = "sorensen",abund = F)$Brich #0.05454545
BAT::beta(x12,func = "sorensen",abund = T)$Brich #0.0034633

BAT::beta(x12,func = "sorensen",abund = F)$Bgain #0.07272727
BAT::beta(x12,func = "sorensen",abund = F)$Bloss #0.1272727
adiv::betastatsor(mite_beta)[3] #Ricotta & Pavoine 0.1272727
BAT::beta(x12_pa,func = "sorensen",abund = T)$Bgain #0.2672327
BAT::beta(x12_pa,func = "sorensen",abund = T)$Bloss #0.2637694



