#_______________________________________________________________________________
# Title              : 08_Bibliometry.r
# Date               : 09/03/2026
# Object             : Script to represent the number of publications involving
#                      biodiversity and/or marine and/or policy
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

packages_needed <- c("readr","dplyr","tidyr","stringr","ggplot2")

loadpackages(packages_needed)

#_______________________________________________________________________________####
#________________________________Loading data___________________________________####

data <- read_delim("data/Scopus_data.csv", 
                          delim = ";", escape_double = FALSE, trim_ws = TRUE)


#________________________________ Data viz______________________________________####
options(scipen = 999)

data$`Query Groupk` <- factor(data$`Query Groupk`, levels = 
                                c(
                                  "Group Biodiversity",
                                  "Group Biodiversity AND Group policy",
                                  "Group Biodiversity AND Group Marine",
                                  "Group Biodiversity AND Group Marine AND Group policy"
                                )
)

ggplot(filter(data,Year<=2025 & Year >= 1970))+
  geom_point(aes(x=Year,y=Article,colour=`Query Groupk`,shape = `Query Groupk`),size=3,alpha = 0.7)+
  geom_line(aes(x=Year,y=Article,colour=`Query Groupk`),size=1.5,alpha = 0.7)+
  scale_y_continuous(n.breaks = 10,breaks = waiver(),labels = scales::label_number(big.mark = " "))+
  scale_x_continuous(n.breaks = 15,breaks = waiver())+
  scale_shape(solid = T)+
  scale_shape_manual(values = c(
    "Group Biodiversity" = 19,
    "Group Biodiversity AND Group policy" = 15,
    "Group Biodiversity AND Group Marine" = 17,
    "Group Biodiversity AND Group Marine AND Group policy" = 18
  )) +
  
  theme_light()+
  theme(
    legend.position = c(0.025, 0.975),
    legend.justification = c("left", "top"),
    axis.text.y = element_text(size = 12,face = "bold"),
    axis.text.x = element_text(size = 12,face = "bold"),
    axis.title = element_text(size = 15,face = "bold"),
  ) +
  labs(x="Year",y="Number of publications",colour = "Key words",shape = "Key words")+
  scale_colour_discrete(palette = c("sienna","darkorange","navy","steelblue"))
ggsave('bibliometry.png', path = "output/fig/bibliometry", dpi = 900, width = 250, height = 150, units = 'mm')

