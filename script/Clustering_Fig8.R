# Clustering_Fig8.R — Fig. 8 (classification des indices alpha, mite) : classification + dendrogramme
# Packages : cluster, clValid, mclust, dendextend.

REPO <- "../Biodiversity_indices_compendium"
dir.create("sorties", showWarnings = FALSE)
suppressMessages({library(cluster); library(clValid); library(mclust); library(dendextend)})

## ---- Parametres ----
f_data   <- file.path(REPO, "data/alpha/alpha_values_mite_wide.csv")
drop_col <- c("Sample", "DSP")            # DSP deja retire dans le script 03 (l. 97)
k_fig    <- 3                              # k retenu (coude ; ARI bootstrap = 0.91 attendu)
marker_cols <- c(EBP = "cyan3", RICH = "pink3", DSHA = "orange3")  # un indice marqueur par groupe -> couleur
                                           # EBP -> dominance ; RICH -> richesse (avec DQ, DFIS) ; DSHA -> heterogeneite + equitabilite
n_boot   <- 200
out_png  <- "sorties/Fig8_Clustering_alpha_wardD2.png"
px       <- 1500

## ---- 1. Classification ----
wcss <- function(X, cl) sum(sapply(unique(cl), function(g){ p <- X[cl == g, , drop = FALSE]; sum(sweep(p, 2, colMeans(p))^2) }))
d <- read.csv(f_data, check.names = FALSE)
x <- as.matrix(d[, setdiff(names(d), drop_col)])

# Regroupement des indices redondants (|r| > 0.9999) : le premier indice represente les suivants
r <- abs(cor(x)); keep <- rep(TRUE, ncol(x)); dup <- list()
for (j in seq_len(ncol(x))) if (keep[j]) {
  dd <- which(r[j, ] > 0.9999 & seq_len(ncol(x)) > j & keep)
  if (length(dd)) { keep[dd] <- FALSE; dup[[colnames(x)[j]]] <- colnames(x)[dd] }
}
xs <- x[, keep]
Xi <- t(scale(xs)); D <- dist(Xi); hc <- hclust(D, "ward.D2")

tab <- t(sapply(2:10, function(k) {
  cl <- cutree(hc, k)
  c(k = k, sil = mean(silhouette(cl, D)[, "sil_width"]), dunn = dunn(D, cl), wcss = wcss(Xi, cl))
}))
d2  <- c(NA, diff(tab[, "wcss"], differences = 2), NA)
kel <- tab[which.max(d2), "k"]

set.seed(1)
ari <- sapply(2:6, function(k) {
  cl0 <- cutree(hc, k)
  mean(replicate(n_boot, {
    i <- sample(nrow(xs), replace = TRUE); xb <- scale(xs[i, ]); xb[is.na(xb)] <- 0
    adjustedRandIndex(cl0, cutree(hclust(dist(t(xb)), "ward.D2"), k))
  }))
})
cat("\n===== alpha :", ncol(x), "->", ncol(xs), "indices ; k coude =", kel, "\n")
print(round(cbind(tab, d2 = d2), 3))
cat("ARI bootstrap k=2..6 :", round(ari, 3), "\n")
cat("Doublons retires :\n"); for (k in names(dup)) cat("  ", k, "=", paste(dup[[k]], collapse = ", "), "\n")
if (kel != k_fig) warning("Le coude (k = ", kel, ") differe du k de la figure (k = ", k_fig, ").")

## ---- 2. Dendrogramme (couleurs selon le contenu des groupes) ----
fam <- function(l) ifelse(grepl("^R", l), "#8B6914", ifelse(grepl("^T", l), "darkgreen",
                   ifelse(grepl("^[DQ]", l), "navy", "darkmagenta")))
cl <- cutree(hc, k_fig); dend <- as.dendrogram(hc); lab <- labels(dend)
grp_order  <- unique(cl[lab])                      # ordre d'apparition des groupes dans le dendrogramme
col_of_grp <- sapply(grp_order, function(g) {
  m <- names(marker_cols)[sapply(names(marker_cols), function(z) cl[z] == g)]
  if (length(m)) marker_cols[[m[1]]] else "grey40"
})
dend <- color_branches(dend, k = k_fig, col = unname(col_of_grp))
dend <- color_labels(dend, col = fam(lab)); dend <- set(dend, "branches_lwd", 3); dend <- set(dend, "labels_cex", 0.95)
png(out_png, width = px, height = px, res = 150); par(mar = c(3, 1, 1, 7))
plot(dend, horiz = TRUE, axes = TRUE); dev.off()
cat("Figure ecrite :", out_png, "\n")
cat("Groupes (k =", k_fig, ") :\n"); for (g in sort(unique(cl))) cat("  groupe", g, ":", paste(names(cl)[cl == g], collapse = ", "), "\n")

saveRDS(list(hc = hc, keep = colnames(xs), dup = dup, tab = tab, kel = kel, ari = ari, cl = cl),
        "sorties/Fig8_Clustering_alpha.rds")
