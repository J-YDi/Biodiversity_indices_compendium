#_______________________________________________________________________________
# Title              : 02_Discriminate_alpha_functions.r
# Date               : 25/03/2025
# Object             : Script to discrimate functions that return a false value
#                      of alpha diverisity indices
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
                     "tidyverse","dendextend","circlize","corrplot","ggh4x")

loadpackages(packages_needed)

#_______________________________________________________________________________####
#________________________________Loading data___________________________________####

data <- read_csv("data/alpha/a_combined_long_all.csv")

#__________________________Viz false values from packages_______________________####

# On indique pas les fonctions custom car ce n'est pas le but du graphe
data <- filter(data,package != "custom")

count_same_values_by_col <- function(df, digits = 3, exclude_self = TRUE,treat_NA_as_value = FALSE) {
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



results_E <- filter(results,startsWith(Index, "E"))

line_df_E <- data.frame(
  package = sort(unique(results_E$package)),
  col = rep(c("black", "grey70"), length.out = length(unique(results_E$package)))
)

ggplot(results_E) +
  aes(x = Index, y = package, fill = Percentage) +
  geom_segment(
    data = line_df_E,
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
  scale_y_discrete(limits = rev(sort(unique(results_E$package)))) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    panel.background = NULL
  ) +
  facet_wrap(~Index, scales = "free_x", ncol = 50)+
  theme(strip.text = element_text(face = "bold", color = "white",
                                  size = 10),
        strip.background = element_rect(fill = "darkmagenta"))
ggsave('heatmap_alpha_E_packages_TF.png', path = "output/fig/alpha/packages/", dpi = 600, width = 400, height = 200, units = 'mm')

results_D <- filter(results,startsWith(Index, "D") | startsWith(Index, "Q"))

line_df_D <- data.frame(
  package = sort(unique(results_D$package)),
  col = rep(c("black", "grey70"), length.out = length(unique(results_D$package)))
)

ggplot(results_D) +
  aes(x = Index, y = package, fill = Percentage) +
  geom_segment(
    data = line_df_D,
    aes(y = package, yend = package, x = -Inf, xend = Inf, color = col),
    inherit.aes = FALSE,
    linewidth = 0.3
  ) +
  scale_color_identity() +
  geom_tile() +
  geom_text(
    aes(label = scales::percent(Percentage, accuracy = 1)), 
    size = 3, na.rm = TRUE
  ) +
  scale_fill_gradient(low = "darkorchid1", high = "gold", name = "Agreement (%)",
                      na.value = "darkorchid1") +
  scale_y_discrete(limits = rev(sort(unique(results_D$package)))) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    panel.background = NULL
  ) +
  facet_wrap(~Index, scales = "free_x", ncol = 50)+
  theme(strip.text = element_text(face = "bold", color = "white",
                                  size = 10),
        strip.background = element_rect(fill = "darkblue"))
ggsave('heatmap_alpha_D_packages_TF.png', path = "output/fig/alpha/packages/", dpi = 600, width = 500, height = 300, units = 'mm')

results_R <- filter(results, startsWith(Index, "R") | startsWith(Index, "T"))

# ordonner les packages en ordre alphabétique
results_R$package <- factor(results_R$package,
                            levels = sort(unique(results_R$package)))

# extrait l'ordre des facettes
index_levels <- sort(unique(results_R$Index))

line_df_R <- data.frame(
  package = sort(unique(results_R$package)),
  col = rep(c("black", "grey70"), length.out = length(unique(results_R$package)))
)

strip_colors <- lapply(index_levels, function(x) {
  if (grepl("^T", x)) {
    element_rect(fill = "darkgreen")
  } else {
    element_rect(fill = "darkgoldenrod4")
  }
})

ggplot(results_R) +
  aes(x = Index, y = package, fill = Percentage) +
  geom_segment(
    data = line_df_R,
    aes(y = package, yend = package, x = -Inf, xend = Inf, color = col),
    inherit.aes = FALSE, linewidth = 0.3
  ) +
  scale_color_identity() +
  geom_tile() +
  geom_text(aes(label = scales::percent(Percentage, accuracy = 1)), 
            size = 2.5, na.rm = TRUE) +
  scale_fill_gradient(low = "darkorchid1", high = "gold",
                      name = "Agreement (%)", na.value = "darkorchid1") +
  scale_y_discrete(limits = rev(levels(results_R$package))) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    legend.position = "bottom",
    panel.background = NULL
  ) +
  ggh4x::facet_wrap2(
    vars(Index),
    scales = "free_x",
    ncol = 50,
    strip = strip_themed(
      background_x = strip_colors,
      text_x = element_text(face = "bold", colour = "white", size = 10)
    )
  )
ggsave('heatmap_alpha_RT_packages_TF.png', path = "output/fig/alpha/packages/", dpi = 600, width = 400, height = 200, units = 'mm')


withNAorInf <- unique(select(filter(data,is.na(value) | value == Inf),-Sample))

