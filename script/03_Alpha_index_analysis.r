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
                     "tidyverse","dendextend","circlize","corrplot","cowplot")

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
