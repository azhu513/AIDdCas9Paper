## ------------------------------------------------
## Script name: ggplot_themes.R
## Purpose of script: ggplot themes for publication
##                    plots
## Author: Anqi Zhu
## Date Created: 2023-05-26
## ------------------------------------------------
## Notes:
##   theme_mine - modified based on theme_bw
##                larger fonts
## ------------------------------------------------

theme_mine <- function(base_size = 14, base_family = "") {
  # Starts with theme_grey and then modify some parts
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(size = 18),
      strip.text.y = element_text(size = 18),
      axis.text.x = element_text(size=14),
      axis.text.y = element_text(size=14, hjust=1),
      axis.ticks =  element_line(colour = "black"), 
      axis.title.x= element_text(size=18),
      axis.title.y= element_text(size=18, angle=90),
      #legend.position = "none", 
      panel.background = element_blank(), 
      panel.border = element_blank(), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.margin = unit(1.0, "lines"), 
      plot.background = element_blank(), 
      plot.margin = unit(c(0.5,  0.5, 0.5, 0.5), "lines"),
      axis.line = element_line(colour = "black"),
      legend.background=element_rect(fill='white'),
      legend.title=element_text(size=18, colour="black"),
      legend.text=element_text(size=14, colour="black"),
      legend.key = element_rect( fill = 'white'),
      legend.key.size = unit(c(1, 1), "lines")
    )
}
