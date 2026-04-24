#_______________________________________________________________________________
# Title              : 07_Compare_transformations_beta_indices.r
# Date               : 22/04/2025
# Object             : Script to determine the influence on transformations on 
#                      on beta biodiversity indices
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

#_______________________________________________________________________________####
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
#________________________________Loading data___________________________________####
library(vegan)
# We choose the mite data from vegan as the dataset for abundance/count data
data("mite")

# Relative abundances data to allow some functions working
mite_relat <- mite/rowSums(mite)

# Hellinger transformation
mite_hellinger <- decostand(mite, method = "hellinger")

# Detach the vegan package
detach(package:vegan)

#___________________________ Beta diversity indices ____________________________####

Raw_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Raw_61_A[i] <- vegan::vegdist(mite_beta,method = "manhattan")
}
Relative_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_relat[c(i,i+1),]
  Relative_61_A[i] <- vegan::vegdist(mite_beta,method = "manhattan")
}
Hellinger_61_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_hellinger[c(i,i+1),]
  Hellinger_61_A[i] <- vegan::vegdist(mite_beta,method = "manhattan")
}

Raw_95_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Raw_95_A[i] <- ecodive::topsoe(mite_beta,rescale = F)
}
Relative_95_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_relat[c(i,i+1),]
  Relative_95_A[i] <- ecodive::topsoe(mite_beta,rescale = F)
}
Hellinger_95_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_hellinger[c(i,i+1),]
  Hellinger_95_A[i] <- ecodive::topsoe(mite_beta,rescale = F)
}

Raw_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Raw_43_A[i] <- vegan::vegdist(mite_beta,method = "jaccard",binary = F) 
}
Relative_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_relat[c(i,i+1),]
  Relative_43_A[i] <- vegan::vegdist(mite_beta,method = "jaccard",binary = F) 
}
Hellinger_43_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_hellinger[c(i,i+1),]
  Hellinger_43_A[i] <- vegan::vegdist(mite_beta,method = "jaccard",binary = F) 
}

Raw_100_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Raw_100_A[i] <- proxy::dist(mite_beta,method = "Whittaker")
}
Relative_100_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_relat[c(i,i+1),]
  Relative_100_A[i] <- proxy::dist(mite_beta,method = "Whittaker")
}
Hellinger_100_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_hellinger[c(i,i+1),]
  Hellinger_100_A[i] <- proxy::dist(mite_beta,method = "Whittaker")
}

Raw_9_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Raw_9_A[i] <- vegan::vegdist(mite_beta,method = "bray")
}
Relative_9_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_relat[c(i,i+1),]
  Relative_9_A[i] <- vegan::vegdist(mite_beta,method = "bray")
}
Hellinger_9_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_hellinger[c(i,i+1),]
  Hellinger_9_A[i] <- vegan::vegdist(mite_beta,method = "bray")
}

Raw_59_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Raw_59_A[i] <- EnvNJ::metrics(t(mite_beta),method = "lorentzian") 
}
Relative_59_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_relat[c(i,i+1),]
  Relative_59_A[i] <- EnvNJ::metrics(t(mite_beta),method = "lorentzian") 
}
Hellinger_59_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Hellinger_59_A[i] <- EnvNJ::metrics(t(mite_beta),method = "lorentzian") 
}

Raw_15_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Raw_15_A[i] <- vegan::vegdist(mite_beta,method = "chao")
}
Relative_15_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_relat[c(i,i+1),]
  Relative_15_A[i] <- vegan::vegdist(mite_beta,method = "chao")
}
Hellinger_15_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_hellinger[c(i,i+1),]
  Hellinger_15_A[i] <- vegan::vegdist(mite_beta,method = "chao")
}

Raw_41_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Raw_41_A[i] <- EnvNJ::metrics(t(mite_beta),method = "hamming")[2,1]
}
Relative_41_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_relat[c(i,i+1),]
  Relative_41_A[i] <- EnvNJ::metrics(t(mite_beta),method = "hamming")[2,1]
}
Hellinger_41_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_hellinger[c(i,i+1),]
  Hellinger_41_A[i] <- EnvNJ::metrics(t(mite_beta),method = "hamming")[2,1]
}

Raw_67_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite[c(i,i+1),]
  Raw_67_A[i] <- vegan::vegdist(mite_beta,method = "morisita")
}
Relative_67_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_relat[c(i,i+1),]
  Relative_67_A[i] <- vegan::vegdist(mite_beta,method = "morisita")
}
Hellinger_67_A <- rep(NA,69)
for (i in 1:(nrow(mite)-1)){
  mite_beta <- mite_hellinger[c(i,i+1),]
  Hellinger_67_A[i] <- vegan::vegdist(mite_beta,method = "morisita")
}

Indexes <- ls(pattern = "_A$")
Indexes <- mget(Indexes)
beta_combined <- as.data.frame(Indexes)


beta_combined$Sample <- paste0(1:69, "-", 2:70)
beta_combined_long <- pivot_longer(beta_combined,cols = colnames(beta_combined)[1:27],names_to = "transformation_index",values_to = "value")
beta_combined_long <- beta_combined_long |>
  separate(transformation_index, into = c("Transformation", "index"), sep = "_",extra = "drop")
write.csv(beta_combined_long,"data/beta/b_combined_long_transformed_mite.csv",row.names = F)

# Data viz #
data_long <- read_csv("data/beta/b_combined_long_transformed_mite.csv")

#______________________Some basic representations of the data___________________####

# Calculate some stats to the plot
data_stats <- data_long |> 
  group_by(index,Transformation) |> 
  summarise(mean_value = mean(value),
            sd_value = sd(value))

data_long$index <- paste(data_long$index,"-",data_long$Transformation)
data_stats$index <- paste(data_stats$index,"-",data_stats$Transformation)


levels_index <- c("9 - Raw","15 - Raw","41 - Raw",
                  "9 - Relative","15 - Relative","41 - Relative",
                  "9 - Hellinger","15 - Hellinger","41 - Hellinger",
                  "43 - Raw","59 - Raw","61 - Raw",
                  "43 - Relative","59 - Relative","61 - Relative",
                  "43 - Hellinger","59 - Hellinger","61 - Hellinger",
                  "67 - Raw","95 - Raw","100 - Raw",
                  "67 - Relative","95 - Relative","100 - Relative",
                  "67 - Hellinger","95 - Hellinger","100 - Hellinger")

data_long$index <- factor(data_long$index, levels = levels_index)

data_stats$index <- factor(data_stats$index, levels = levels_index)

data_long <- data_long %>%
  mutate(
    Sample_num = as.numeric(str_extract(Sample, "^[0-9]+")),
    Sample = fct_reorder(Sample, Sample_num)
  )

data_stats <- data_stats %>%
  left_join(
    data_long %>%
      group_by(index) %>%
      summarise(
        y_label = max(value, na.rm = TRUE) * 0.2
      ),
    by = "index"
  )

ggplot(data_long) +
  geom_segment(aes(x = Sample,y=0, yend = value), col = "royalblue",
               linewidth = 1, alpha = 0.4) +
  geom_point(aes(x = Sample, y = value), col = "royalblue", size = 1.7) +
  geom_label(
    data = data_stats,
    aes(
      x = 69/2,
      y = y_label,
      label = paste0(
        round(mean_value, 3), " \u00B1 ", round(sd_value, 3)
      )
    ),
    hjust = 0.5,
    vjust = 1,
    size = 4,
    alpha = 0.5,
    linewidth = 0
  ) +
  facet_wrap(~index,scale = "free_y",ncol = 3) +
  labs(x = "Sample", y = "Index value") +
  theme(strip.text = element_text(face = "bold", color = "white",
                                  hjust = 0, size = 10),
        strip.background = element_rect(fill = "royalblue"),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.text.x = element_text(size = 3.7,angle=90,hjust = 1,vjust = 0.5))



ggsave('values_beta_A_transform_indices.png', path = "output/fig/beta/indices/", dpi = 300, width = 290, height = 350, units = 'mm')
