# Plot scatter plots for schematic to illustrate I-TED.
# Artificial data

rm(list = ls(all = TRUE))

# LIBRARIES ------------------------------------------------------------------
library(tidyverse)

# PATHS -------------------------------------------------------------------

BASE <- here::here()
FIG_DIR <- file.path(BASE, "analysis", "figures")

# FUNCTIONS ---------------------------------------------------------------
source(file.path(BASE, "src", "plotting_theme.R"))
source(file.path(BASE, "src", "plotting_functions.R"))

generate_plot_data <- function(sd) {
    x <- rnorm(n = 500, 1)
    y <- x + rnorm(500, sd = sd)
    print(cor.test(x, y))
    df <- data.frame(x = x, y = y)

    p <- ggplot(df, aes(x = x, y = y)) +
        geom_point(pch = 21, fill = "grey80") +
        geom_smooth(se = FALSE, method = "lm")

    p <- change_axes(p)
    return(p)
}

set.seed(727)

# about 0.8 cor, so 0.2 dist
p1 <- generate_plot_data(sd = 0.78)
save_ggplot(p1, file.path(FIG_DIR, "SupFig_08_cor"), 40, 40)
# about 0.6 cor, so 0.4 dist
p2 <- generate_plot_data(sd = 1.25)
save_ggplot(p2, file.path(FIG_DIR, "SupFig_06_cor"), 40, 40)
# about 0.7 cor, so 0.3 dist
p3 <- generate_plot_data(sd = 0.97)
save_ggplot(p3, file.path(FIG_DIR, "SupFig_07_cor"), 40, 40)

