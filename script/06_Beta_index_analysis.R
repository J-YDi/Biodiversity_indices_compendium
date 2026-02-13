#_______________________________________________________________________________
# Title              : 06_Beta_index_analysis.r
# Date               : 13/02/2025
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
                     "clValid","vegan","ggrepel")

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

ggplot(data_long) +
  geom_segment(aes(x = Sample,y=0, yend = value), col = "indianred1",
               linewidth = 1, alpha = 0.4) +
  geom_point(aes(x = Sample, y = value), col = "indianred1", size = 1.7) +
  geom_label(data = filter(data_stats),
             aes(x = 69/2, y = 0.01, label = paste0(round(mean_value, 3)," +/- ",round(sd_value, 3))),
             color = "black", size = 4,alpha=0.5,linewidth=0) +
  facet_wrap(~ index,scale = "free_y",ncol = 9) +
  labs(x = "Sample", y = "Index value") +
  theme(strip.text = element_text(face = "bold", color = "white",
                                  hjust = 0, size = 10),
        strip.background = element_rect(fill = "indianred1"),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12,face = "bold"),
        axis.text.x = element_text(size = 3.7,angle=90,hjust = 1,vjust = 0.5))
ggsave('values_beta_P_indices.png', path = "output/fig/beta/indices/", dpi = 900, width = 600, height = 300, units = 'mm')

# ALL VARIABLES SELECT #############____________________________________________
#________________________________________PCA____________________________________####
data_pca <- select(data,-Sample)
PCA_results <- PCA(data_pca)
PCA_results_t <- PCA(t(data_pca))

fviz_screeplot(PCA_results) # Screeplot

# PCA viz with colour arrows

PCA <- fviz_pca_var(PCA_results, axes = c(1, 2), repel = T ,col.var = "indianred1",title="", ggtheme = theme_minimal()) +
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
hc <- hclust(dist(t(data_pca),method = "euclidean"),method = "ward")
plot(hc)

library(cluster)

cluster_quality(data_pca, return_table = TRUE)
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
cluster_quality(data_pca, return_table = TRUE)
nmds <- metaMDS(dist(t(data_pca)), k = 3, trymax = 999)

scores_df <- as.data.frame(scores(nmds))  # x,y
scores_df$Sample <- rownames(scores_df)

clusters <- cutree(hclust(dist(t(data_pca)), method = "ward.D2"), k = 3)
scores_df$Cluster <- factor(clusters)

NMDS <- ggplot(scores_df, aes(x = NMDS1, y = NMDS2, color = Cluster)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample), max.overlaps = Inf, box.padding = 0.5) +
  theme_minimal() +
  geom_label(aes(x=40,y=-15,label = paste("Stress:",round(nmds$stress,4))),
             color = "black",linewidth = 0)+
  theme(legend.position = "none")+
  labs(title = NULL,
       x = "NMDS1", y = "NMDS2")+
  scale_color_discrete(palette = c("red", "blue", "green2", "orange","magenta"))
NMDS
ggsave('NMDS_k3_beta_P_all.png', path = "output/fig/beta/indices/", dpi = 1200, width = 250, height = 150, units = 'mm')

nmds <- metaMDS(dist(t(data_pca)), k = 4, trymax = 999)

scores_df <- as.data.frame(scores(nmds))  # x,y
scores_df$Sample <- rownames(scores_df)

clusters <- cutree(hclust(dist(t(data_pca)), method = "ward.D2"), k = 4)
scores_df$Cluster <- factor(clusters)

ggplot(scores_df, aes(x = NMDS1, y = NMDS2, color = Cluster)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample), max.overlaps = Inf, box.padding = 0.5) +
  theme_minimal() +
  geom_label(aes(x=40,y=-15,label = paste("Stress:",round(nmds$stress,4))),
             color = "black",linewidth = 0)+
  theme(legend.position = "none")+
  labs(title = NULL,
       x = "NMDS1", y = "NMDS2",subtitle = paste("Stress:",round(nmds$stress,4)))+
  scale_color_discrete(palette = c("red", "blue", "green2", "orange","magenta"))
ggsave('NMDS_k4_beta_P_all.png', path = "output/fig/beta/indices/", dpi = 1200, width = 250, height = 150, units = 'mm')

nmds <- metaMDS(dist(t(data_pca)), k = 7, trymax = 999)

scores_df <- as.data.frame(scores(nmds))  # x,y
scores_df$Sample <- rownames(scores_df)

clusters <- cutree(hclust(dist(t(data_pca)), method = "ward.D2"), k = 7)
scores_df$Cluster <- factor(clusters)

ggplot(scores_df, aes(x = NMDS1, y = NMDS2, color = Cluster)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = Sample), max.overlaps = Inf, box.padding = 0.5) +
  theme_minimal() +
  geom_label(aes(x=40,y=-15,label = paste("Stress:",round(nmds$stress,4))),
             color = "black",linewidth = 0)+
  theme(legend.position = "none")+
  labs(title = NULL,
       x = "NMDS1", y = "NMDS2",subtitle = paste("Stress:",round(nmds$stress,4)))+
  scale_color_discrete(palette = c("red", "blue", "green2", "orange","magenta","cyan","yellow"))
ggsave('NMDS_k7_beta_P_all.png', path = "output/fig/beta/indices/", dpi = 1200, width = 250, height = 150, units = 'mm')


