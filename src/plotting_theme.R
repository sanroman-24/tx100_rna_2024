library(ggthemes)

# Set the theme of the figures
theme_set(theme_tufte(base_size = 6, base_family = "ArialMT"))
theme_update(
    text = element_text(size = 6),
    axis.text = element_text(size = 6),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    axis.ticks.length = unit(.1, "cm")
)

# TODO
# Set up palette
