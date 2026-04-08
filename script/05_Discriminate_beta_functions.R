#_______________________________________________________________________________
# Title              : 05_Discriminate_beta_functions.r
# Date               : 08/04/2026
# Object             : Script to discrimate functions that return a false value
#                      of beta diverisity indices
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

packages_needed <- c("readr","dplyr","tidyr","stringr","dendextend","ggplot2",
                     "tidyverse","dendextend","circlize","corrplot")

loadpackages(packages_needed)

#_______________________________________________________________________________####
#________________________________Loading data___________________________________####

data <- read_csv("data/beta/b_combined_long_all_mite.csv")

data[data$package == "BoutrosLabplottinggeneral", "package"] <- "BoutrosLab.plotting.general"
data[data$package == "fdausc", "package"] <- "fda.usc"
data[data$package == "codabase", "package"] <- "coda.base"

#__________________________Viz false values from packages_______________________####

# On indique pas les fonctions custom car ce n'est pas le but du graphe
data <- filter(data,package != "custom")

count_same_values_by_col <- function(df, digits = NULL, exclude_self = TRUE,treat_NA_as_value = FALSE) {
  out <- lapply(df, function(col) {
    g <- col
    if (!is.null(digits) && is.numeric(g)) g <- round(g, digits)
    
    # taille du groupe (même valeur) pour chaque élément de la colonne
    counts <- ave(seq_along(g), g, FUN = length)
    
    # gestion des NA si souhaitée
    if (treat_NA_as_value) {
      na_idx <- is.na(g)
      if (any(na_idx)) counts[na_idx] <- sum(na_idx)
    }
    
    # on ne compte pas la ligne elle-même si demandé
    if (exclude_self) counts <- ifelse(is.na(counts), NA_integer_, counts - 1L)
    
    as.integer(counts)
  })
  
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  rownames(out) <- rownames(df)
  colnames(out) <- colnames(df)
  out
}

results <- map_dfr(unique(data$index), function(idx) {
  
  mat <- data |>
    filter(index == idx) |>
    pivot_wider(names_from = package, values_from = value) |>
    select(-Sample, -index) |>
    round(3) |>
    as.data.frame()
  
  # nombre de packages ayant effectivement calculé l'indice
  n_packages <- ncol(mat)
  
  mat_t <- t(mat)
  
  if (n_packages == 1) {
    # cas trivial : un seul package → accord = 1
    tibble(
      package = rownames(mat_t),
      Percentage = 1,
      Index = idx
    )
  } else {
    agreement <- count_same_values_by_col(as.data.frame(mat_t)) / nrow(mat_t)
    
    tibble(
      package = rownames(agreement),
      Percentage = rowSums(agreement) / ncol(agreement),
      Index = idx,
      maxi = max(Percentage)
    ) |>
      mutate(
        Percentage = ifelse(is.na(Percentage), 0, Percentage),
        Percentage = Percentage / max(Percentage)
      )
  }
})
results <- results |>
  mutate(Percentage = replace_na(Percentage, 0)) |>
  group_by(Index) |>
  mutate(
    Percentage = ifelse(
      n() > 1 & Percentage == 1 & sum(Percentage == 1) == 1,
      0,
      Percentage
    )
  ) |>
  ungroup()



results_P <- filter(results,endsWith(Index, "P"))

levels_index <- results_P$Index %>%
  unique() %>%
  .[order(as.numeric(str_extract(., "\\d+")))]

results_P$Index <- factor(results_P$Index, levels = levels_index)

line_df_P <- data.frame(
  package = sort(unique(results_P$package)),
  col = rep(c("black", "grey70"), length.out = length(unique(results_P$package)))
)

ggplot(results_P) +
  aes(x = Index, y = package, fill = Percentage) +
  geom_segment(
    data = line_df_P,
    aes(y = package, yend = package, x = -Inf, xend = Inf, color = col),
    inherit.aes = FALSE,
    linewidth = 0.3
  ) +
  scale_color_identity() +
  geom_tile() +
  geom_text(
    aes(label = scales::percent(Percentage, accuracy = 1)), 
    size = 3.5, na.rm = TRUE
  ) +
  scale_fill_gradient(low = "darkorchid1", high = "gold", name = "Agreement (%)",
                      na.value = "darkorchid1") +
  scale_y_discrete(limits = rev(sort(unique(results_P$package)))) +
  labs(x = "", y = "") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    panel.background = NULL
  ) +
  facet_wrap(~Index, scales = "free_x", ncol = 31)+
  theme(strip.text = element_text(face = "bold", color = "white",
                                  size = 10),
        strip.background = element_rect(fill = "indianred1"))
ggsave('heatmap_beta_P_packages_TF.png', path = "output/fig/beta/packages/", dpi = 600, width = 500, height = 300, units = 'mm')

results_A <- filter(results,endsWith(Index, "A"))


results_A <- results_A |>
  separate(Index, into = c("SIndex", "DT"), sep = "_",extra = "merge")
results_A$SIndex <- as.numeric(results_A$SIndex)
results_A$Index <- paste0(results_A$SIndex,"_",results_A$DT)

results_A_35 <- filter(results_A,SIndex <= 35)


levels_index_35 <- results_A_35$Index %>%
  unique() %>%
  .[order(as.numeric(str_extract(., "\\d+")))]

results_A_35$Index <- factor(results_A_35$Index, levels = levels_index_35)

line_df_A_35 <- data.frame(
  package = sort(unique(results_A_35$package)),
  col = rep(c("black", "grey70"), length.out = length(unique(results_A_35$package)))
)

ggplot(results_A_35) +
  aes(x = Index, y = package, fill = Percentage) +
  geom_segment(
    data = line_df_A_35,
    aes(y = package, yend = package, x = -Inf, xend = Inf, color = col),
    inherit.aes = FALSE,
    linewidth = 0.3
  ) +
  scale_color_identity() +
  geom_tile() +
  geom_text(
    aes(label = scales::percent(Percentage, accuracy = 1)), 
    size = 2.5, na.rm = TRUE
  ) +
  scale_fill_gradient(low = "darkorchid1", high = "gold", name = "Agreement (%)",
                      na.value = "darkorchid1") +
  scale_y_discrete(limits = rev(sort(unique(results_A_35$package)))) +
  labs(x = "", y = "") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    panel.background = NULL
  ) +
  facet_wrap(~Index, scales = "free_x", ncol = 33)+
  theme(strip.text = element_text(face = "bold", color = "white",
                                  size = 10),
        strip.background = element_rect(fill = "royalblue"))
ggsave('heatmap_beta_A_P1_packages_TF.png', path = "output/fig/beta/packages/", dpi = 900, width = 500, height = 250, units = 'mm')

results_A_120 <- filter(results_A,SIndex > 35)


levels_index_120 <- results_A_120$Index %>%
  unique() %>%
  .[order(as.numeric(str_extract(., "\\d+")))]

results_A_120$Index <- factor(results_A_120$Index, levels = levels_index_120)

line_df_A_120 <- data.frame(
  package = sort(unique(results_A_120$package)),
  col = rep(c("black", "grey70"), length.out = length(unique(results_A_120$package)))
)

ggplot(results_A_120) +
  aes(x = Index, y = package, fill = Percentage) +
  geom_segment(
    data = line_df_A_120,
    aes(y = package, yend = package, x = -Inf, xend = Inf, color = col),
    inherit.aes = FALSE,
    linewidth = 0.3
  ) +
  scale_color_identity() +
  geom_tile() +
  geom_text(
    aes(label = scales::percent(Percentage, accuracy = 1)), 
    size = 2.5, na.rm = TRUE
  ) +
  scale_fill_gradient(low = "darkorchid1", high = "gold", name = "Agreement (%)",
                      na.value = "darkorchid1") +
  scale_y_discrete(limits = rev(sort(unique(results_A_120$package)))) +
  labs(x = "", y = "") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    panel.background = NULL
  ) +
  facet_wrap(~Index, scales = "free_x", ncol = 33)+
  theme(strip.text = element_text(face = "bold", color = "white",
                                  size = 10),
        strip.background = element_rect(fill = "royalblue"))
ggsave('heatmap_beta_A_P2_packages_TF.png', path = "output/fig/beta/packages/", dpi = 900, width = 500, height = 250, units = 'mm')
