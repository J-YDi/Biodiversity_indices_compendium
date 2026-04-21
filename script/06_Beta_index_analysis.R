#_______________________________________________________________________________
# Title              : 06_Beta_index_analysis.r
# Date               : 20/04/2025
# Object             : Script to analyze alpha biodiversity index
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
                     "tidyverse","dendextend","circlize","corrplot","cowplot",
                     "factoextra","FactoMineR","viridis","GGally","cluster",
                     "clValid","vegan","ggrepel","forcats")

loadpackages(packages_needed)

#_____________________________Useful functions__________________________________####
cluster_quality <- function(data, k_min = 2, k_max = 8,
                            indices = c("silhouette", "dunn", "elbow"),
                            hc_method = "ward.D2",dist_meth = "euclidean",
                            return_table = TRUE) {
  
  k_values <- k_min:k_max
  results <- data.frame()
  
  # Distance pour les indices nécessitant dist
  d <- dist(data,dist_meth)
  
  for (k in k_values) {
    # Clustering hiérarchique
    hc <- hclust(d, method = hc_method)
    clusters <- cutree(hc, k = k)
    
    # Silhouette
    if ("silhouette" %in% indices) {
      sil <- silhouette(clusters, d)
      results <- rbind(results,
                       data.frame(k = k,
                                  method = "silhouette",
                                  value = mean(sil[, "sil_width"])))
    }
    
    # Dunn
    if ("dunn" %in% indices) {
      dunn_val <- dunn(d, clusters)
      results <- rbind(results,
                       data.frame(k = k,
                                  method = "dunn",
                                  value = dunn_val))
    }
    
    # Elbow / WCSS
    if ("elbow" %in% indices) {
      wcss <- sum(sapply(unique(clusters), function(cl) {
        pts <- data[clusters == cl, , drop = FALSE]
        sum(rowSums((pts - colMeans(pts))^2))
      }))
      results <- rbind(results,
                       data.frame(k = k,
                                  method = "elbow",
                                  value = wcss))
    }
  }
  
  # Graphiques
  par(mfrow = c(1, length(indices)))
  for (m in indices) {
    sub <- results[results$method == m, ]
    plot(sub$k, sub$value, type = "b",
         xlab = "Nombre de clusters",
         ylab = m,
         main = paste(m, "index"))
  }
  par(mfrow = c(1,1))
  
  if (return_table) return(results)
  invisible(NULL)
}
# Indices PA #####
#__________________________________Loading data_________________________________####

data <- read_csv("data/beta/beta_values_mite_wide_mite_P.csv")
data_long <- read_csv("data/beta/beta_values_mite_long_mite_P.csv")

#______________________Some basic representations of the data___________________####

# # All
# ggplot(data_long)+
#   geom_line(aes(x=Sample,y=value))+
#   facet_wrap(~index,scale = "free_y")

# Calculate some stats to the plot
data_stats <- data_long |> 
  group_by(index) |> 
  summarise(mean_value = mean(value),
            sd_value = sd(value))

# PA indices

levels_index <- data_long$index %>%
  unique() %>%
  .[order(as.numeric(str_extract(., "\\d+")))]

data_long$index <- factor(data_long$index, levels = levels_index)

levels_index <- data_stats$index %>%
  unique() %>%
  .[order(as.numeric(str_extract(., "\\d+")))]

data_stats$index <- factor(data_stats$index, levels = levels_index)

data_long <- data_long %>%
  mutate(
    Sample_num = as.numeric(str_extract(Sample, "^[0-9]+")),
    Sample = fct_reorder(Sample, Sample_num)
  )

data_long_P <- filter(data_long, !index %in% c("87_P","88_P","89_P","90_P","91_P","92_P",
                                               "45_P","46_P","48_P","49_P","50_P","51_P","52_P") )
data_stats_P <- filter(data_stats, !index %in% c("87_P","88_P","89_P","90_P","91_P","92_P",
                                               "45_P","46_P","48_P","49_P","50_P","51_P","52_P") )

ggplot(data_long_P) +
  geom_segment(aes(x = Sample,y=0, yend = value), col = "indianred1",
               linewidth = 1, alpha = 0.4) +
  geom_point(aes(x = Sample, y = value), col = "indianred1", size = 1.7) +
  geom_label(data = filter(data_stats_P, !index %in% c("5_P","64_P","119_P","112_P","63_P","40_P")),
             aes(x = 69/2, y = 0.04, label = paste0(round(mean_value, 3)," +/- ",round(sd_value, 3))),
             color = "black", size = 4,alpha=0.5,linewidth=0) +
  geom_label(data = filter(data_stats_P, index %in% c("5_P","64_P","119_P","112_P","63_P","40_P")),
             aes(x = 69/2, y = 0.9, label = paste0(round(mean_value, 3)," +/- ",round(sd_value, 3))),
             color = "black", size = 4,alpha=0.5,linewidth=0) +
  facet_wrap(~ index,scale = "free_y",ncol = 5) +
  labs(x = "Sample", y = "Index value") +
  theme(strip.text = element_text(face = "bold", color = "white",
                                  hjust = 0, size = 10),
        strip.background = element_rect(fill = "indianred1"),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.text.x = element_text(size = 3.7,angle=90,hjust = 1,vjust = 0.5))
ggsave('values_beta_P_indices.png', path = "output/fig/beta/indices/", dpi = 900, width = 600, height = 300, units = 'mm')

ggplot(filter(data_stats_P, mean_value <= 1 & mean_value >= 0)) +
  geom_point(aes(x = reorder(index, mean_value), y = mean_value)) +
  geom_errorbar(aes(
    x = reorder(index, mean_value),
    ymin = mean_value - sd_value,
    ymax = mean_value + sd_value
  ), width = 0.2)+
  labs(x = "Index", y = "Mean index value") +
  theme(axis.text.x = element_text(size = 12, angle=90,hjust = 1,vjust = 0.5), 
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title = element_text(size = 15))
ggsave('mean_values_beta_P_indices.png', path = "output/fig/beta/indices/", dpi = 900, width = 250, height = 150, units = 'mm')

ggplot(filter(data_stats_P)) +
  geom_point(aes(x = reorder(index, abs(sd_value/mean_value)), y = abs(sd_value/mean_value))) +
  labs(x = "Index", y = "Mean index value") +
  theme(axis.text.x = element_text(size = 12, angle=90,hjust = 1,vjust = 0.5), 
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title = element_text(size = 15))

levels_index <- data_long$index %>%
  unique() %>%
  .[order(as.numeric(str_extract(., "\\d+")))]

data_long$index <- factor(data_long$index, levels = levels_index)

zeros <- ggplot(filter(data_long, value == 0 & !index %in% c("87_P","88_P","89_P","90_P","91_P","92_P",
                                                          "45_P","46_P","48_P","49_P","50_P","51_P","52_P"))) +
  geom_point(aes(x = Sample, y = index),size=10,shape=15) +
  labs(x = "Samples", y = "Index with 0") +
  theme(axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title = element_text(size = 15))
ones <- ggplot(filter(data_long, value == 1 & !index %in% c("87_P","88_P","89_P","90_P","91_P","92_P",
                                                            "45_P","46_P","48_P","49_P","50_P","51_P","52_P"))) +
  geom_point(aes(x = Sample, y = index),size=12,shape=15) +
  labs(x = "Samples", y = "Index with 1") +
  theme(axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title = element_text(size = 15))
plot_grid(zeros,ones,rel_widths = c(3,1))
ggsave('extreme_values_beta_P_indices.png', path = "output/fig/beta/indices/", dpi = 900, width = 220, height = 100, units = 'mm')



data_long_Pdec <- filter(data_long, index %in% c("87_P","88_P","89_P","90_P","91_P","92_P",
                                               "45_P","46_P","48_P","49_P","50_P","51_P","52_P") )
data_stats_Pdec <- filter(data_stats, index %in% c("87_P","88_P","89_P","90_P","91_P","92_P",
                                                 "45_P","46_P","48_P","49_P","50_P","51_P","52_P") )

ggplot(data_long_Pdec) +
  geom_segment(aes(x = Sample,y=0, yend = value), col = "turquoise3",
               linewidth = 1, alpha = 0.4) +
  geom_point(aes(x = Sample, y = value), col = "turquoise3", size = 1.7) +
  geom_label(data = filter(data_stats_Pdec),
             aes(x = 69/2, y = 0.01, label = paste0(round(mean_value, 3)," +/- ",round(sd_value, 3))),
             color = "black", size = 4,alpha=0.5,linewidth=0) +
  facet_wrap(~ index,scale = "free_y",ncol = 3) +
  labs(x = "Sample", y = "Index value") +
  theme(strip.text = element_text(face = "bold", color = "white",
                                  hjust = 0, size = 10),
        strip.background = element_rect(fill = "turquoise3"),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.text.x = element_text(size = 3.7,angle=90,hjust = 1,vjust = 0.5))
ggsave('values_beta_Pdec_indices.png', path = "output/fig/beta/indices/", dpi = 900, width = 300, height = 200, units = 'mm')




# ALL VARIABLES SELECT #############____________________________________________
#________________________________________PCA____________________________________####
data_pca <- select(data,-Sample)
PCA_results <- PCA(data_pca)
PCA_results_t <- PCA(t(data_pca))

fviz_screeplot(PCA_results) # Screeplot

# PCA viz with colour arrows

# Create the color code
group_PC <- names(data_pca)[names(data_pca) %in% c("87_P","88_P","89_P","90_P","91_P","92_P",
                                          "45_P","46_P","48_P","49_P","50_P","51_P","52_P")]
group_P <-  names(data_pca)[!names(data_pca) %in% c("87_P","88_P","89_P","90_P","91_P","92_P",
                                           "45_P","46_P","48_P","49_P","50_P","51_P","52_P")]


# Distinguishing the dimensions
color_vector <- rep("Misrepresented", length(colnames(data_pca)))
names(color_vector) <- colnames(data_pca)
color_vector[group_PC] <- "PC"
color_vector[group_P] <- "P"


# PCA viz with colour arrows

PCA <- fviz_pca_var(PCA_results, axes = c(1, 2), repel = T ,col.var = color_vector, legend = "none", 
                    palette = c("indianred1","turquoise3"),title="", ggtheme = theme_minimal()) +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )
PCA
ggsave('PCA_beta_P_all.png', path = "output/fig/beta/indices/", dpi = 1200, width = 250, height = 250, units = 'mm')

fviz_contrib(PCA_results, choice = "var", axes = 1)
fviz_contrib(PCA_results, choice = "var", axes = 2)

corrplot(t(PCA_results$var$contrib),
         is.corr = FALSE,
         method = "pie",col = viridis(200),number.cex = 0.5)

# To check for a spatial dissimilarity we check by see the individuals position by region
fviz_pca_ind(PCA_results_t,addEllipses = F,repel = T,col.ind = "indianred1"
             ,title="", ggtheme = theme_minimal(),legend = "none")

#______________________Clusterings______________________________________________####
# 
data_pca_scaled <- scale(data_pca,T,T)
hc <- hclust(dist(t(data_pca_scaled),method = "euclidean"),method = "ward")
plot(hc)

library(cluster)

cluster_quality(data_pca_scaled, return_table = TRUE)
k=4
clusters <- cutree(hc, k = k)
cluster_cols <- c("", "blue", "darkgreen", "orange")


label_cols <- cluster_cols[clusters]
dend <- as.dendrogram(hc)
# dend <- dendrapply(dend, function(n) {
#   if (!is.leaf(n)) attr(n, "height") <- log1p(attr(n, "height"))
#   n
# })

dend <- color_branches(dend,k = k, col = c("orchid","olivedrab","thistle4","yellow3") )

dend <- color_labels(dend,col = c(rep("indianred1",9),rep("turquoise3",1),rep("indianred1",13),"turquoise3",
                                  rep("indianred1",1),rep("turquoise3",3),rep("indianred1",5),"turquoise3",rep("indianred1",12),"turquoise3",rep("indianred1",1),
                                  rep("turquoise3",2),"indianred1",
                                  rep("turquoise3",2),"indianred1","turquoise3",rep("indianred1",2)))
dend <- set(dend, "branches_lwd",3 )

plot(dend, horiz = T,dLeaf = -0.3,axes=T)


# NMDS 
cluster_quality(data_pca_scaled, return_table = TRUE)
nmds <- metaMDS(dist(t(data_pca_scaled)), k = 4, trymax = 999)

scores_df <- as.data.frame(scores(nmds))  # x,y
scores_df$Sample <- rownames(scores_df)

clusters <- cutree(hclust(dist(t(data_pca_scaled)), method = "ward.D2"), k = 3)
scores_df$Cluster <- factor(clusters)

NMDS <- ggplot(scores_df, aes(x = NMDS1, y = NMDS2, color = Cluster)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample), max.overlaps = Inf, box.padding = 0.5) +
  theme_minimal() +
  geom_label(aes(x=-5,y=-10,label = paste("Stress:",round(nmds$stress,4))),
             color = "black",linewidth = 0)+
  theme(legend.position = "none")+
  labs(title = NULL,
       x = "NMDS1", y = "NMDS2")+
  scale_color_discrete(palette = c("red", "blue", "green2", "orange","magenta"))
NMDS
ggsave('NMDS_k3_beta_P_all.png', path = "output/fig/beta/indices/", dpi = 1200, width = 250, height = 150, units = 'mm')


# Indices A #####
#__________________________________Loading data_________________________________####

data <- read_csv("data/beta/beta_values_mite_wide_mite_A.csv")
data_long <- read_csv("data/beta/beta_values_mite_long_mite_A.csv")

#______________________Some basic representations of the data___________________####

# # All
# ggplot(data_long)+
#   geom_line(aes(x=Sample,y=value))+
#   facet_wrap(~index,scale = "free_y")

# Calculate some stats to the plot
data_stats <- data_long |> 
  group_by(index) |> 
  summarise(mean_value = mean(value),
            sd_value = sd(value))

# A indices

levels_index <- data_long$index %>%
  unique() %>%
  .[order(as.numeric(str_extract(., "\\d+")))]

data_long$index <- factor(data_long$index, levels = levels_index)

levels_index <- data_stats$index %>%
  unique() %>%
  .[order(as.numeric(str_extract(., "\\d+")))]

data_stats$index <- factor(data_stats$index, levels = levels_index)

data_long_A <- filter(data_long, !index %in% c("87_A","88_A","89_A","90_A","91_A","92_A",
                                               "45_A","46_A","48_A","49_A","50_A","51_A","52_A") )
data_stats_A <- filter(data_stats, !index %in% c("87_A","88_A","89_A","90_A","91_A","92_A",
                                                 "45_A","46_A","48_A","49_A","50_A","51_A","52_A") )


data_stats_A <- data_stats_A %>%
  left_join(
    data_long_A %>%
      group_by(index) %>%
      summarise(
        y_label = max(value, na.rm = TRUE) * 0.4
      ),
    by = "index"
  )

ggplot(data_long_A) +
  geom_segment(
    aes(x = Sample, y = 0, yend = value),
    color = "royalblue",
    linewidth = 1,
    alpha = 0.4
  ) +
  
  geom_point(
    aes(x = Sample, y = value),
    color = "royalblue",
    size = 1.7
  ) +
  
  geom_label(
    data = data_stats_A,
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
  
  facet_wrap(~ index, scales = "free_y", ncol = 5) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  
  labs(
    x = "Sample",
    y = "Index value"
  ) +
  
  theme(
    strip.text = element_text(
      face = "bold",
      color = "white",
      hjust = 0,
      size = 10
    ),
    strip.background = element_rect(
      fill = "royalblue"
    ),
    axis.title = element_text(size = 15),
    axis.text.y = element_text(
      size = 12,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 3.7,
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )

ggsave('values_beta_A_indices.png', path = "output/fig/beta/indices/", dpi = 300, width = 600, height = 300, units = 'mm')

ggplot(filter(data_stats_A, mean_value <= 1 & mean_value >= 0)) +
  geom_point(aes(x = reorder(index, mean_value), y = mean_value)) +
  geom_errorbar(aes(
    x = reorder(index, mean_value),
    ymin = mean_value - sd_value,
    ymax = mean_value + sd_value
  ), width = 0.2)+
  labs(x = "Index", y = "Mean index value") +
  theme(axis.text.x = element_text(size = 12, angle=90,hjust = 1,vjust = 0.5), 
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title = element_text(size = 15))
ggsave('mean_values_beta_A_indices.png', path = "output/fig/beta/indices/", dpi = 300, width = 250, height = 150, units = 'mm')

ggplot(filter(data_stats_A)) +
  geom_point(aes(x = reorder(index, abs(sd_value/mean_value)), y = abs(sd_value/mean_value))) +
  labs(x = "Index", y = "Mean index value") +
  theme(axis.text.x = element_text(size = 12, angle=90,hjust = 1,vjust = 0.5), 
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title = element_text(size = 15))

levels_index <- data_long$index %>%
  unique() %>%
  .[order(as.numeric(str_extract(., "\\d+")))]

data_long$index <- factor(data_long$index, levels = levels_index)

ggplot(filter(data_long, value == 0 & !index %in% c("87_A","88_A","89_A","90_A","91_A","92_A",
                                                             "45_A","46_A","48_A","49_A","50_A","51_A","52_A"))) +
  geom_point(aes(x = Sample, y = index),size=25,shape=15) +
  labs(x = "Samples", y = "Index with 0") +
  theme(axis.text.x = element_text(size = 8), 
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.title = element_text(size = 15))
ggsave('extreme_values_beta_A_indices.png', path = "output/fig/beta/indices/", dpi = 300, width = 200, height = 100, units = 'mm')



data_long_Adec <- filter(data_long, index %in% c("87_A","88_A","89_A","90_A","91_A","92_A",
                                                 "45_A","46_A","48_A","49_A","50_A","51_A","52_A") )
data_stats_Adec <- filter(data_stats, index %in% c("87_A","88_A","89_A","90_A","91_A","92_A",
                                                   "45_A","46_A","48_A","49_A","50_A","51_A","52_A"))
ggplot(data_long_Adec) +
  geom_segment(aes(x = Sample,y=0, yend = value), col = "turquoise3",
               linewidth = 1, alpha = 0.4) +
  geom_point(aes(x = Sample, y = value), col = "turquoise3", size = 1.7) +
  geom_label(data = filter(data_stats_Adec),
             aes(x = 69/2, y = 0.01, label = paste0(round(mean_value, 3)," +/- ",round(sd_value, 3))),
             color = "black", size = 4,alpha=0.5,linewidth=0) +
  facet_wrap(~ index,scale = "free_y",ncol = 3) +
  labs(x = "Sample", y = "Index value") +
  theme(strip.text = element_text(face = "bold", color = "white",
                                  hjust = 0, size = 10),
        strip.background = element_rect(fill = "turquoise3"),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.text.x = element_text(size = 3.7,angle=90,hjust = 1,vjust = 0.5))
ggsave('values_beta_Adec_indices.png', path = "output/fig/beta/indices/", dpi = 900, width = 300, height = 200, units = 'mm')

# ALL VARIABLES SELECT #############____________________________________________
#________________________________________PCA____________________________________####
data_pca <- select(data,-Sample)
PCA_results <- PCA(data_pca)
PCA_results_t <- PCA(t(data_pca))

fviz_screeplot(PCA_results) # Screeplot

# PCA viz with colour arrows

PCA <- fviz_pca_var(PCA_results, axes = c(1, 2), repel = T ,col.var = "royalblue",title="", ggtheme = theme_minimal()) +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )
PCA
ggsave('PCA_beta_A_all.png', path = "output/fig/beta/indices/", dpi = 1200, width = 250, height = 250, units = 'mm')

fviz_contrib(PCA_results, choice = "var", axes = 1)
fviz_contrib(PCA_results, choice = "var", axes = 2)
fviz_contrib(PCA_results, choice = "var", axes = 3)

corrplot(t(PCA_results$var$contrib),
         is.corr = FALSE,
         method = "pie",col = viridis(200),number.cex = 0.5)

# To check for a spatial dissimilarity we check by see the individuals position by region
fviz_pca_ind(PCA_results_t,addEllipses = F,repel = T,col.ind = "royalblue"
             ,title="", ggtheme = theme_minimal(),legend = "none")

#______________________Clusterings______________________________________________####
# 
data_pca_scaled <- scale(data_pca,T,T)


hc <- hclust(dist(t(data_pca_scaled),method = "euclidean"),method = "ward")
plot(hc)

library(cluster)

cluster_quality(data_pca_scaled, return_table = TRUE)
k=3
clusters <- cutree(hc, k = k)
cluster_cols <- c("red", "blue", "darkgreen", "orange")
label_cols <- cluster_cols[clusters]
dend <- as.dendrogram(hc)
dend <- dendrapply(dend, function(n) {
  if (!is.leaf(n)) attr(n, "height") <- log1p(attr(n, "height"))
  n
})

dend <- color_branches(dend,k = k )
dend <- color_labels(dend,k = k )
dend <- set(dend, "branches_lwd", k)
plot(dend, horiz = T,dLeaf = -0.1,axes=T)


# NMDS 
cluster_quality(data_pca_scaled, return_table = TRUE)
nmds <- metaMDS(dist(t(data_pca_scaled)), k = 3, trymax = 999)

scores_df <- as.data.frame(scores(nmds))  # x,y
scores_df$Sample <- rownames(scores_df)

clusters <- cutree(hclust(dist(t(data_pca_scaled)), method = "ward.D2"), k = 3)
scores_df$Cluster <- factor(clusters)

NMDS <- ggplot(scores_df, aes(x = NMDS1, y = NMDS2, color = Cluster)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample), max.overlaps = Inf, box.padding = 0.5) +
  theme_minimal() +
  geom_label(aes(x=-5,y=-8,label = paste("Stress:",round(nmds$stress,4))),
             color = "black",linewidth = 0)+
  theme(legend.position = "none")+
  labs(title = NULL,
       x = "NMDS1", y = "NMDS2")+
  scale_color_discrete(palette = c("red", "blue", "green2", "orange","magenta"))
NMDS
ggsave('NMDS_k3_beta_A_all.png', path = "output/fig/beta/indices/", dpi = 1200, width = 250, height = 150, units = 'mm')

