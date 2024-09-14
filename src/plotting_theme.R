library(ggthemes)

# Set the theme of the figures
theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(
    text = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    axis.ticks.length = unit(.1, "cm"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
)

# Set up palette
tx_palette <- c(
    "lightblue" = "#A3C6E8",
    "darkblue" = "#3979bb",
    "darkred" = "#DF2A55",
    "lightred" = "#E76A85",
    "darkpurple" = "#6A4C93",
    "lightpurple" = "#B89BCF",
    "darkgreen" = "#2A9D8F",
    "lightgreen" = "#83D0C9", 
    "lightgrey" = "grey90", 
    "darkgrey" = "grey40"
)
