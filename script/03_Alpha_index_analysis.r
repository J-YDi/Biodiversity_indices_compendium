#_______________________________________________________________________________
# Title              : 03_Alpha_index_analysis.r
# Date               : 27/01/2025
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
                     "factoextra","FactoMineR","viridis","GGally","cluster")

loadpackages(packages_needed)

#_______________________________________________________________________________####

#__________________________________Loading data_________________________________####

data <- read_csv("data/alpha/alpha_values_mite_wide.csv")
data_long <- read_csv("data/alpha/alpha_values_mite_long.csv")

#______________________Some basic representations of the data___________________####

# All
ggplot(data_long)+
  geom_line(aes(x=Sample,y=value))+
  facet_wrap(~index,scale = "free_y")

# Calculate some stats to the plot
data_stats <- data_long |> 
  group_by(index) |> 
  summarise(mean_value = mean(value),
            sd_value = sd(value))

# Equitability indices

ggplot(filter(data_long, startsWith(index, "E"))) +
  geom_line(aes(x = Sample, y = value), col = "darkmagenta",
            linewidth = 1.5, alpha = 0.4) +
  geom_point(aes(x = Sample, y = value), col = "darkmagenta", size = 1.7) +
  geom_text(data = filter(data_stats, startsWith(index, "E")),
            aes(x = 62, y = 0, label = paste0(round(mean_value, 3)," +/- ",round(sd_value, 3))),
            color = "black", size = 4) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~ index) +
  labs(x = "Sample", y = "Index value") +
  theme(strip.text = element_text(face = "bold", color = "white",
                                  hjust = 0, size = 10),
        strip.background = element_rect(fill = "darkmagenta"),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12,face = "bold"))
ggsave('values_alpha_E_indices.png', path = "output/fig/alpha/indices/", dpi = 900, width = 400, height = 200, units = 'mm')

# Heterogeneity indices

ggplot(filter(data_long, startsWith(index, "D"))) +
  geom_line(aes(x = Sample, y = value), col = "darkblue",
            linewidth = 1.5, alpha = 0.4) +
  geom_point(aes(x = Sample, y = value), col = "darkblue", size = 1.7) +
  geom_text(data = filter(data_stats, startsWith(index, "D")),
            aes(x = 64, y = 0, label = paste0(round(mean_value, 3)," +/- ",round(sd_value, 3))),
            color = "black", size = 4) +
  facet_wrap(~ index,scales="free_y",ncol=3) +
  labs(x = "Sample", y = "Index value") +
  theme(strip.text = element_text(face = "bold", color = "white",
                                  hjust = 0, size = 10),
        strip.background = element_rect(fill = "darkblue"),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12,face = "bold"))
ggsave('values_alpha_D_indices.png', path = "output/fig/alpha/indices/", dpi = 900, width = 400, height = 200, units = 'mm')

# Estimate of species richness indices

A <- ggplot(filter(data_long, startsWith(index, "R") & !index %in% c("RHUR_2","RHUR_3"))) +
  geom_line(aes(x = Sample, y = value), col = "darkgoldenrod4",
            linewidth = 1.5, alpha = 0.4) +
  geom_point(aes(x = Sample, y = value), col = "darkgoldenrod4", size = 1.7) +
  geom_text(data = filter(data_stats, startsWith(index, "R")& !index %in% c("RHUR_2","RHUR_3")),
            aes(x = 61, y = 58, label = paste0(round(mean_value, 3)," +/- ",round(sd_value, 3))),
            color = "black", size = 4) +
  facet_wrap(~ index,ncol=3) +
  labs(x = "Sample", y = "Index value") +
  theme(strip.text = element_text(face = "bold", color = "white",
                                  hjust = 0, size = 10),
        strip.background = element_rect(fill = "darkgoldenrod4"),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12,face = "bold"))
B <- ggplot(filter(data_long, startsWith(index, "R") & index %in% c("RHUR_2","RHUR_3"))) +
  geom_line(aes(x = Sample, y = value), col = "darkgoldenrod4",
            linewidth = 1.5, alpha = 0.4) +
  geom_point(aes(x = Sample, y = value), col = "darkgoldenrod4", size = 1.7) +
  geom_text(data = filter(data_stats, startsWith(index, "R")& index %in% c("RHUR_2","RHUR_3")),
            aes(x = 61, y = 3, label = paste0(round(mean_value, 3)," +/- ",round(sd_value, 3))),
            color = "black", size = 4) +
  scale_y_continuous(limits = c(1, 3)) +
  facet_wrap(~ index,nrow=3) +
  labs(x = "", y = NULL) +
  theme(strip.text = element_text(face = "bold", color = "white",
                                  hjust = 0, size = 10),
        strip.background = element_rect(fill = "darkgoldenrod4"),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12,face = "bold"))
plot_grid(A,B,rel_widths = c(3,1),rel_heights = c(1,2))
ggsave('values_alpha_R_indices.png', path = "output/fig/alpha/indices/", dpi = 900, width = 400, height = 200, units = 'mm')

# Mixte Q indices

# Heterogeneity indices

ggplot(filter(data_long, startsWith(index, "Q"))) +
  geom_line(aes(x = Sample, y = value), col = "darkgreen",
            linewidth = 1.5, alpha = 0.4) +
  geom_point(aes(x = Sample, y = value), col = "darkgreen", size = 1.7) +
  geom_text(data = filter(data_stats, startsWith(index, "Q")),
            aes(x = 61, y = 0, label = paste0(round(mean_value, 3)," +/- ",round(sd_value, 3))),
            color = "black", size = 4) +
  facet_wrap(~ index,scales="free_y",ncol=4) +
  labs(x = "Sample", y = "Index value") +
  theme(strip.text = element_text(face = "bold", color = "white",
                                  hjust = 0, size = 10),
        strip.background = element_rect(fill = "darkgreen"),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(size = 12,face = "bold"))
ggsave('values_alpha_Q_indices.png', path = "output/fig/alpha/indices/", dpi = 900, width = 400, height = 200, units = 'mm')

#________________________________________PCA____________________________________####
data_pca <- select(data,-Sample)
PCA_results <- PCA(data_pca)
PCA_results_t <- PCA(t(data_pca))

fviz_screeplot(PCA_results) # Screeplot

# Create the color code
group_D <- names(data)[startsWith(names(data), "D")]
group_E <- names(data)[startsWith(names(data), "E")]
group_Q <- names(data)[startsWith(names(data), "Q")]
group_R <- names(data)[startsWith(names(data), "R")]

# Distinguishing the dimensions
color_vector <- rep("Misrepresented", length(colnames(data_pca)))
names(color_vector) <- colnames(data_pca)
color_vector[group_D] <- "D"
color_vector[group_E] <- "E"
color_vector[group_Q] <- "Q"
color_vector[group_R] <- "R"

# PCA viz with colour arrows

fviz_pca_var(PCA_results, axes = c(1, 2), repel = T ,col.var = color_vector, legend = "none", 
             palette = c("darkblue","darkmagenta","darkgreen","darkgoldenrod4"),title="", ggtheme = theme_minimal()) +
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

#ggsave('PCA_1_V2.png', path = "output/graphs/networks/", dpi = 600, width = 100, height = 100, units = 'mm')

fviz_contrib(PCA_results, choice = "var", axes = 1)
fviz_contrib(PCA_results, choice = "var", axes = 2)

corrplot(t(PCA_results$var$contrib),
         is.corr = FALSE,
         method = "pie",col = viridis(200),number.cex = 0.5)

# To check for a spatial dissimilarity we check by see the individuals position by region
fviz_pca_ind(PCA_results_t,addEllipses = F,repel = T,col.ind = color_vector,,
             palette = c("darkblue","darkmagenta","darkgreen","darkgoldenrod4"),
             ,title="", ggtheme = theme_minimal(),legend = "none")

#_______________________Correlations____________________________________________####
ggpairs(data_pca)

cor.mtest <- function(mat, method = "pearson") {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      test <- cor.test(mat[, i], mat[, j], method = method)
      p.mat[i, j] <- p.mat[j, i] <- test$p.value
    }
  }
  diag(p.mat) <- 0
  return(p.mat)
}

p.mat <- cor.mtest(data_pca)

corrplot(cor(data_pca),
         method = "shade",col = viridis(200),number.cex = 0.5,order = "alphabet",
         addCoef.col = NULL,tl.col = "black",
         diag = F,type = "full",addgrid.col = "black",addCoefasPercent = T,
         insig = "blank",sig.level = 0.05,p.mat = p.mat
)

#______________________Clusterings______________________________________________####

hc <- hclust(dist(t(data_pca),method = "euclidean"),method = "ward")
plot(hc)

library(cluster)

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

k=4
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

k=5
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
plot(dend, horiz = T,dLeaf = -0.1,axes=T)cl 

# Kmeans

cl <- kmeans(x=data_pca,centers = 2)
plot(data_pca,col = cl$cluster)
     