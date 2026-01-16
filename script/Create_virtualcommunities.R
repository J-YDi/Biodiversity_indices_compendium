#_______________________________________________________________________________
# Nom               : Create_virtualcommunities.r
# Date de modif     : 18/09/2025
# Objet             : Creer des communautes virtuelles
# Auteurs           : J-Y. Dias
# Version R         : 4.5.0
#_______________________________________________________________________________

#_______________________________ Packages_______________________________________

# Fonction permettant d'installer les packages si non presents et/ou les charger
# en l etat ne fonctionne pas pour les packages qui s'installent hors CRAN
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
packages_needed <- c("dplyr","ggplot2","tidyr")

loadpackages(packages_needed)

#_______________________Creer un signal saisonnier______________________________#####

# mois <- 1:24 #creer un vecteur indiquant le mois, ici pour 2 ans
# mois_an <- (mois - 1) %% 12 + 1   # recode le mois pour revenir a 1 apres decembre

# Passage en hebdomadaire

hebdo <- 1:(52*4) # 52 semaines dans un an et 4 an
hebdo_an <- (hebdo-1) %% 52 + 1 # idem mois

# Creer un fonction gaussienne pour obtenir une forme en cloche
# mu = centre de la gaussienne, sigma largeur et amp amplitude
gaussian <- function(x, mu, sigma, amp=1) {
  amp * exp(-((x-mu)^2)/(2*sigma^2))
}

# Creer une fonction chi2 pour obtenir une forme "apparition soudaine et decroissance lente"
# df = forme de la distribution et miu = pic souhaitee
chi2_shifted <- function(x, mu, df, amp=1) {
  # Mode de la loi chi² standard
  mode_chi2 <- ifelse(df > 2, df - 2, 0) # pour ne pas avoir une decroissance depuis 0 car loi de chi2 -> df-2
  
  # Decalage pour que le maximum corresponde a mu
  x_shifted <- (x - mu) + mode_chi2
  
  # Densite chi2 avec decalage, bornee a x >= 0
  dens <- ifelse(x_shifted >= 0, dchisq(x_shifted, df=df), 0) # pour eviter les negatifs
  
  # Mise a l’echelle par l’amplitude
  amp * dens
}

# Fonction pour revenir uniquement entre 0 et 1, sera utile pour controler les abondances
scale01 <- function(x) (x - min(x)) / (max(x) - min(x))

#__________________Generer 10 especes avec profils typiques_____________________####

set.seed(123)
# Espece 1 : pic de printemps (mois 4 donc semaine 16) avec allure chi2
esp1 <- rpois(52*4, 1500 * scale01(chi2_shifted(hebdo_an, mu=16, df=2, amp=0.5)))
#       distrib de poisson pour ne pas avoir trop de variation car variance ~= moyenne
#       24 nombre de points et 1500 facteur pour obtenir une abondance (1500 x [0-1]) : a randomiser aussi 

# Espece 2 : pic de printemps (mois 4) avec allure gaussienne
esp2 <- rpois(52*4, 1000 * scale01(gaussian(hebdo_an, mu=16, sigma=3.5)))

# Espece 3 : pic d'ete (mois 8) allure gaussienne
esp3 <- rpois(52*4, 1000 * scale01(gaussian(hebdo_an, mu=33, sigma=1.5)))

# Espece 4 : bimodale (printemps + automne)
esp4 <- rpois(52*4, 800 * scale01(gaussian(hebdo_an, 16, 1.5) + gaussian(hebdo_an, 36, 1)))

# Espece 5 : pic d'hiver
esp5 <- rpois(52*4, 450 * scale01(gaussian(hebdo_an, 50, 3)+gaussian(hebdo_an, 3, 3)))

# Espece 6 : generaliste pas hiver
esp6 <- rpois(52*4, 400 * scale01(gaussian(hebdo_an, 27, 15)))

# Espece 7 : cosmopolite / perenne
esp7 <- rpois(52*4, 600)

# Espece 8 : aleatoire opportuniste
esp8 <- rpois(52*4, runif(24, 100, 1000))

# Espece 9 : rare assez abondante
esp9 <-  runif(52*4, 0, 600) # plus variable que rpois mais peut etre ressemblant avec 8

# Espece 10 : rare peu abondante
esp10 <-  rbinom(52*4, size=1, prob=0.2) * 200 # presence absence et quand present valeur 200

#______________Combiner les especes pour creer la communaute____________________####

comm <- cbind(esp1, esp2, esp3, esp4, esp5, esp6,esp7,esp8,esp9,esp10)
colnames(comm) <- c("Bloom printemps","Pic printemps","Pic ete","Bimodale Spring/Autumn","Hivernale","Generale hors hiver","Cosmopolite","Opportuniste","Rare abondante","Rare peu abondante")

comm_df <- as.data.frame(comm) %>%
  mutate(Semaine = hebdo) %>%
  pivot_longer(-Semaine, names_to="Esp", values_to="Abond") |>
  mutate(Abond = ifelse(Abond < 200, 0, Abond)) # pour representer la valeur minimale

#__________________________Visualisation________________________________________####

ggplot(comm_df, aes(x=Semaine, y=Abond, colour=Esp)) +
  geom_line(size=1) +
  geom_point(size=2) +
  facet_wrap(~Esp, ncol=2) +
  geom_vline(aes(xintercept = 52),colour = "gray56")+  
  geom_vline(aes(xintercept = 52*2),colour = "gray56")+
  geom_vline(aes(xintercept = 52*3),colour = "gray56")+
  labs(title="Dynamiques especes comm virtuelle",
       x="Semaine", y="Abondance") +
  theme(legend.position = "none")

#_____________Meme principe pour generer 100 especes de chaque type_____________####

#### Type 1 ####----------------------------------------------------------------

abd_type1 <- seq(from = 200, to = 3000, by = 100) # gamme abondance max pour l'espece x de type 1
df_type1 <- c(0:3) # gamme allure pour espece x de type 1
amp_type1 <- seq(from = 0.1, to = 2, by = 0.1) # gamme amplitude espece x de type 1
mu_type1 <- seq(from = 14, to = 17, by = 0.5) # semaine pic espece x de type 1

n_esp <- 100

type1 <- replicate(n_esp, {
  abd <- sample(abd_type1, 1) # selectionne une abondance max pour l'espece x de type 1 parmi la gamme
  mu  <- sample(mu_type1, 1) # meme principe pour mu, df, amp
  df  <- sample(df_type1, 1)
  amp <- sample(amp_type1, 1)
  distrib <- abd * scale01(chi2_shifted(hebdo_an, mu = mu, df = df, amp = amp))
  distrib[is.na(distrib)] <- 0
  if(all(distrib == 0)) rep(0, length(hebdo)) else rpois(length(hebdo), distrib)
})

assemb_type1 <- as.data.frame(type1) %>%
  mutate(Semaine = hebdo) %>%
  pivot_longer(-Semaine, names_to="Esp", values_to="Abond") |>
  mutate(Abond = ifelse(Abond < 200, 0, Abond),
         Type = "1")

ggplot(assemb_type1, aes(x=Semaine, y=Abond, colour=Esp)) +
  geom_line(size=1) +
  geom_point(size=2) +
  facet_wrap(~Esp, ncol=10) +
  geom_vline(aes(xintercept = 52),colour = "gray56")+
  geom_vline(aes(xintercept = 52*2),colour = "gray56")+
  geom_vline(aes(xintercept = 52*3),colour = "gray56")+
  labs(title="Assemblage type 1",
       x="Semaine", y="Abondance") +
  theme(legend.position = "none")

#### Type 2 ####----------------------------------------------------------------

abd_type2 <- seq(from = 200, to = 3000, by = 100)
amp_type2 <- seq(from = 0.1, to = 2, by = 0.1)
mu_type2 <- seq(from = 14, to = 17, by = 0.5)
sigma_type2 <- seq(from = 0.5, to = 3, by =0.2)

n_esp <- 100

type2 <- replicate(n_esp, {
  abd <- sample(abd_type2, 1)
  mu  <- sample(mu_type2, 1)
  amp <- sample(amp_type2, 1)
  sigma <- sample(sigma_type2,1)
  distrib <- abd * scale01(gaussian(hebdo_an, mu = mu, amp = amp,sigma = sigma))
  distrib[is.na(distrib)] <- 0
  if(all(distrib == 0)) rep(0, length(hebdo)) else rpois(length(hebdo), distrib)
})

assemb_type2 <- as.data.frame(type2) %>%
  mutate(Semaine = hebdo) %>%
  pivot_longer(-Semaine, names_to="Esp", values_to="Abond") |>
  mutate(Abond = ifelse(Abond < 200, 0, Abond),
         Type = "2")

ggplot(assemb_type2, aes(x=Semaine, y=Abond, colour=Esp)) +
  geom_line(size=1) +
  geom_point(size=2) +
  facet_wrap(~Esp, ncol=10) +
  geom_vline(aes(xintercept = 52),colour = "gray56")+
  geom_vline(aes(xintercept = 52*2),colour = "gray56")+
  geom_vline(aes(xintercept = 52*3),colour = "gray56")+
  labs(title="Assemblage type2",
       x="Semaine", y="Abondance") +
  theme(legend.position = "none")

#### Type 3 ####----------------------------------------------------------------

abd_type3 <- seq(from = 200, to = 3000, by = 100)
amp_type3 <- seq(from = 0.1, to = 2, by = 0.1)
mu_type3 <- seq(from = 27, to = 33, by = 0.5)
sigma_type3 <- seq(from = 0.5, to = 2, by =0.1)

n_esp <- 100

type3 <- replicate(n_esp, {
  abd <- sample(abd_type3, 1)
  mu  <- sample(mu_type3, 1)
  amp <- sample(amp_type3, 1)
  sigma <- sample(sigma_type3,1)
  distrib <- abd * scale01(gaussian(hebdo_an, mu = mu, amp = amp,sigma = sigma))
  distrib[is.na(distrib)] <- 0
  if(all(distrib == 0)) rep(0, length(hebdo)) else rpois(length(hebdo), distrib)
})

assemb_type3 <- as.data.frame(type3) %>%
  mutate(Semaine = hebdo) %>%
  pivot_longer(-Semaine, names_to="Esp", values_to="Abond") |>
  mutate(Abond = ifelse(Abond < 200, 0, Abond),
         Type = "3")

ggplot(assemb_type3, aes(x=Semaine, y=Abond, colour=Esp)) +
  geom_line(size=1) +
  geom_point(size=2) +
  facet_wrap(~Esp, ncol=10) +
  geom_vline(aes(xintercept = 52),colour = "gray56")+
  geom_vline(aes(xintercept = 52*2),colour = "gray56")+
  geom_vline(aes(xintercept = 52*3),colour = "gray56")+
  labs(title="Assemblage type3",
       x="Semaine", y="Abondance") +
  theme(legend.position = "none")

#### Type 4 ####----------------------------------------------------------------

abd_type4 <- seq(from = 200, to = 2000, by = 100)
amp_type4 <- seq(from = 0.1, to = 2, by = 0.1)
mu_type4_1 <- seq(from = 14, to = 17, by = 0.5)
mu_type4_2 <- seq(from = 30, to = 43, by = 0.5)

sigma_type4 <- seq(from = 0.7, to = 1.7, by =0.2)

n_esp <- 100

type4 <- replicate(n_esp, {
  abd <- sample(abd_type4, 1)
  mu_1  <- sample(mu_type4_1, 1)
  mu_2  <- sample(mu_type4_2, 1)
  amp <- sample(amp_type4, 1)
  sigma_1 <- sample(sigma_type4,1)
  sigma_2 <- sample(sigma_type4,1)
  
  distrib <- abd * scale01(gaussian(hebdo_an, mu_1, sigma_1) + gaussian(hebdo_an, mu_2, sigma_2))
  distrib[is.na(distrib)] <- 0
  if(all(distrib == 0)) rep(0, length(hebdo)) else rpois(length(hebdo), distrib)
})

assemb_type4 <- as.data.frame(type4) %>%
  mutate(Semaine = hebdo) %>%
  pivot_longer(-Semaine, names_to="Esp", values_to="Abond") |>
  mutate(Abond = ifelse(Abond < 200, 0, Abond),
         Type = "4")

ggplot(assemb_type4, aes(x=Semaine, y=Abond, colour=Esp)) +
  geom_line(size=1) +
  geom_point(size=2) +
  facet_wrap(~Esp, ncol=10) +
  geom_vline(aes(xintercept = 52),colour = "gray56")+
  geom_vline(aes(xintercept = 52*2),colour = "gray56")+
  geom_vline(aes(xintercept = 52*3),colour = "gray56")+
  labs(title="Assemblage type4",
       x="Semaine", y="Abondance") +
  theme(legend.position = "none")


#### Type 5 ####----------------------------------------------------------------

abd_type5 <- seq(from = 200, to = 1000, by = 100)
amp_type5 <- seq(from = 0.1, to = 2, by = 0.1)
mu_type5_1 <- seq(from = 46, to = 54, by = 0.5)
mu_type5_2 <- seq(from = 1, to = 10, by = 0.5)

sigma_type5 <- seq(from = 3.4, to = 3.7, by =0.1)

n_esp <- 100

type5 <- replicate(n_esp, {
  abd <- sample(abd_type5, 1)
  mu_1  <- sample(mu_type5_1, 1)
  amp <- sample(amp_type5, 1)
  sigma_1 <- sample(sigma_type5,1)
  
  distrib <- abd * scale01(gaussian(hebdo_an, mu_1, sigma_1) + gaussian(hebdo_an, abs(mu_1-47), sigma_1))
  distrib[is.na(distrib)] <- 0
  if(all(distrib == 0)) rep(0, length(hebdo)) else rpois(length(hebdo), distrib)
})

assemb_type5 <- as.data.frame(type5) %>%
  mutate(Semaine = hebdo) %>%
  pivot_longer(-Semaine, names_to="Esp", values_to="Abond") |>
  mutate(Abond = ifelse(Abond < 200, 0, Abond),
         Type = "5")

ggplot(assemb_type5, aes(x=Semaine, y=Abond, colour=Esp)) +
  geom_line(size=1) +
  geom_point(size=2) +
  facet_wrap(~Esp, ncol=10) +
  geom_vline(aes(xintercept = 52),colour = "gray56")+
  geom_vline(aes(xintercept = 52*2),colour = "gray56")+
  geom_vline(aes(xintercept = 52*3),colour = "gray56")+
  labs(title="Assemblage type5",
       x="Semaine", y="Abondance") +
  theme(legend.position = "none")

#### Type 6 ####----------------------------------------------------------------

abd_type6 <- seq(from = 200, to = 1000, by = 100)
amp_type6 <- seq(from = 0.1, to = 2, by = 0.1)
mu_type6 <- seq(from = 27, to = 33, by = 0.5)
sigma_type6 <- seq(from = 6, to = 10, by =0.4)

n_esp <- 100

type6 <- replicate(n_esp, {
  abd <- sample(abd_type6, 1)
  mu  <- sample(mu_type6, 1)
  amp <- sample(amp_type6, 1)
  sigma <- sample(sigma_type6,1)
  distrib <- abd * scale01(gaussian(hebdo_an, mu = mu, amp = amp,sigma = sigma))
  distrib[is.na(distrib)] <- 0
  if(all(distrib == 0)) rep(0, length(hebdo)) else rpois(length(hebdo), distrib)
})

assemb_type6 <- as.data.frame(type6) %>%
  mutate(Semaine = hebdo) %>%
  pivot_longer(-Semaine, names_to="Esp", values_to="Abond") |>
  mutate(Abond = ifelse(Abond < 200, 0, Abond),
         Type = "6")

ggplot(assemb_type6, aes(x=Semaine, y=Abond, colour=Esp)) +
  geom_line(size=1) +
  geom_point(size=2) +
  facet_wrap(~Esp, ncol=10) +
  geom_vline(aes(xintercept = 52),colour = "gray56")+
  geom_vline(aes(xintercept = 52*2),colour = "gray56")+
  geom_vline(aes(xintercept = 52*3),colour = "gray56")+
  labs(title="Assemblage type6",
       x="Semaine", y="Abondance") +
  theme(legend.position = "none")

#### Type 7 ####----------------------------------------------------------------

abd_type7 <- seq(from = 200, to = 1000, by = 100)

n_esp <- 100

type7 <- replicate(n_esp, {
  abd <- sample(abd_type7, 1)
  distrib <-  rpois(hebdo, abd)
  distrib[is.na(distrib)] <- 0
  if(all(distrib == 0)) rep(0, length(hebdo)) else distrib
})

assemb_type7 <- as.data.frame(type7) %>%
  mutate(Semaine = hebdo) %>%
  pivot_longer(-Semaine, names_to="Esp", values_to="Abond") |>
  mutate(Abond = ifelse(Abond < 200, 0, Abond),
         Type = "7")

ggplot(assemb_type7, aes(x=Semaine, y=Abond, colour=Esp)) +
  geom_line(size=1) +
  geom_point(size=2) +
  facet_wrap(~Esp, ncol=10) +
  geom_vline(aes(xintercept = 52),colour = "gray56")+
  geom_vline(aes(xintercept = 52*2),colour = "gray56")+
  geom_vline(aes(xintercept = 52*3),colour = "gray56")+
  labs(title="Assemblage type7",
       x="Semaine", y="Abondance") +
  theme(legend.position = "none")

#### Type 8 ####----------------------------------------------------------------

abd_type8_1 <- seq(from = 200, to = 600, by = 100)
abd_type8_2 <- seq(from = 400, to = 800, by = 100)
abd_type8_3 <- seq(from = 800, to = 1200, by = 100)

n_esp <- 100

type8 <- replicate(n_esp, {
  abd_1 <- sample(abd_type8_1, 1)
  abd_2 <- sample(abd_type8_2, 1)
  abd_3 <- sample(abd_type8_3, 1)
  distrib <-  rpois(52*4, runif(abd_1, abd_2, abd_3))
  distrib[is.na(distrib)] <- 0
  if(all(distrib == 0)) rep(0, length(hebdo)) else distrib
})

assemb_type8 <- as.data.frame(type8) %>%
  mutate(Semaine = hebdo) %>%
  pivot_longer(-Semaine, names_to="Esp", values_to="Abond") |>
  mutate(Abond = ifelse(Abond < 200, 0, Abond),
         Type = "8")

ggplot(assemb_type8, aes(x=Semaine, y=Abond, colour=Esp)) +
  geom_line(size=1) +
  geom_point(size=2) +
  facet_wrap(~Esp, ncol=10) +
  geom_vline(aes(xintercept = 52),colour = "gray56")+
  geom_vline(aes(xintercept = 52*2),colour = "gray56")+
  geom_vline(aes(xintercept = 52*3),colour = "gray56")+
  labs(title="Assemblage type8",
       x="Semaine", y="Abondance") +
  theme(legend.position = "none")

#### Type 9 ####----------------------------------------------------------------

abd_type9 <- seq(from = 200, to = 800, by = 100)

n_esp <- 100

type9 <- replicate(n_esp, {
  abd <- sample(abd_type9, 1)
  distrib <-  runif(52*4, 0, abd)
  distrib[is.na(distrib)] <- 0
  if(all(distrib == 0)) rep(0, length(hebdo)) else distrib
})

assemb_type9 <- as.data.frame(type9) %>%
  mutate(Semaine = hebdo) %>%
  pivot_longer(-Semaine, names_to="Esp", values_to="Abond") |>
  mutate(Abond = ifelse(Abond < 200, 0, Abond),
         Type = "9")

ggplot(assemb_type9, aes(x=Semaine, y=Abond, colour=Esp)) +
  geom_line(size=1) +
  geom_point(size=2) +
  facet_wrap(~Esp, ncol=10) +
  geom_vline(aes(xintercept = 52),colour = "gray56")+
  geom_vline(aes(xintercept = 52*2),colour = "gray56")+
  geom_vline(aes(xintercept = 52*3),colour = "gray56")+
  labs(title="Assemblage type9",
       x="Semaine", y="Abondance") +
  theme(legend.position = "none")

#### Type 10 ####----------------------------------------------------------------

abd_type10 <- seq(from = 200, to = 400, by = 100)
prob_type10 <- seq(from = 0.05, to = 0.15, by = 0.05)

n_esp <- 100

type10 <- replicate(n_esp, {
  abd <- sample(abd_type10, 1)
  prob <- sample(prob_type10,1)
  distrib <-  rbinom(52*4, size=1, prob=prob) * abd
  distrib[is.na(distrib)] <- 0
  if(all(distrib == 0)) rep(0, length(hebdo)) else distrib
})

assemb_type10 <- as.data.frame(type10) %>%
  mutate(Semaine = hebdo) %>%
  pivot_longer(-Semaine, names_to="Esp", values_to="Abond") |>
  mutate(Abond = ifelse(Abond < 200, 0, Abond),
         Type = "10")

ggplot(assemb_type10, aes(x=Semaine, y=Abond, colour=Esp)) +
  geom_line(size=1) +
  geom_point(size=2) +
  facet_wrap(~Esp, ncol=10) +
  geom_vline(aes(xintercept = 52),colour = "gray56")+
  geom_vline(aes(xintercept = 52*2),colour = "gray56")+
  geom_vline(aes(xintercept = 52*3),colour = "gray56")+
  labs(title="Assemblage type10",
       x="Semaine", y="Abondance") +
  theme(legend.position = "none")

#### Communaute complete ####---------------------------------------------------


comm_virtuel <- rbind(assemb_type1,assemb_type2,assemb_type3,assemb_type4,assemb_type5,
                      assemb_type6,assemb_type7,assemb_type8,assemb_type9,assemb_type10)

