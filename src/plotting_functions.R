library(ggplot2)
library(tidyverse)
library(lemon)
library(survival)
library(ggpubr)
library(survminer)
library(gridExtra)

scatter_plot <- function(df, x_str, y_str, y_title, x_title, fill = NA) {
  if (is.na(fill)) {
    p <- ggplot(df, aes_string(x = x_str, y = y_str)) +
      geom_point(pch = 21, size = 2, fill = "grey40") +
      coord_capped_cart(bottom = "none", left = "none") +
      labs(
        y = y_title,
        x = x_title
      )
    return(p)
  }
  p <- ggplot(df, aes_string(x = x_str, y = y_str)) +
    geom_point(pch = 21, size = 2, aes_string(fill)) +
    coord_capped_cart(bottom = "none", left = "none") +
    labs(
      y = y_title,
      x = x_title
    )
  return(p)
}

violin_plot <- function(
    df, x_str, y_str, y_title,
    x_title, labs_x, fun.y, ylim) {
  ggplot(df, aes_string(x = x_str, y = y_str)) +
    geom_violin(alpha = 0.3) +
    geom_boxplot(width = 0.1, alpha = .3, outlier.shape = NA) +
    geom_point(
      size = .2, position = position_dodge2(width = .3),
      pch = 21, fill = "white", alpha = .8
    ) +
    stat_summary(
      geom = "point",
      fun.y = fun.y,
      col = "black",
      size = 1,
      shape = 21,
      fill = "red"
    ) +
    scale_x_discrete(labels = labs_x) +
    ggpubr::stat_compare_means(size = 2) +
    labs(x = x_title, y = y_title) +
    theme(
      axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 5),
      axis.text.y = element_text(size = 5),
      legend.position = "none"
    ) +
    ylim(ylim) +
    lemon::coord_capped_cart(bottom = "both", left = "bottom") +
    theme(panel.border = element_blank(), axis.line = element_line())
}

change_axes <- function(p) {
  p <- p + coord_capped_cart(bottom = "none", left = "none")
  return(p)
}

save_ggplot <- function(p, fp, w, h) {
  mm_to_inches <- 25.3
  w <- w / mm_to_inches
  h <- h / mm_to_inches
  ggsave(paste0(fp, ".svg"), p, width = w, height = h, device = "svg")
  ggsave(paste0(fp, ".png"), p, width = w, height = h, dpi = 400, device = "png")
  ggsave(paste0(fp, ".pdf"), p, width = w, height = h)
}

save_plist <- function(plist, fp, w, h, ncol) {
  mm_to_inches <- 25.3
  w <- w / mm_to_inches
  h <- h / mm_to_inches
  ggsave(paste0(fp, ".svg"), arrangeGrob(grobs = plist, ncol = ncol),
    width = w, height = h, device = "svg"
  )
  ggsave(paste0(fp, ".png"), arrangeGrob(grobs = plist, ncol = ncol),
    width = w, height = h, dpi = 400, device = "png"
  )
  ggsave(paste0(fp, ".pdf"), arrangeGrob(grobs = plist, ncol = ncol),
    width = w, height = h
  )
}

save_baseplot <- function(p, fp, w, h) {
  mm_to_inches <- 25.3
  w <- w / mm_to_inches
  h <- h / mm_to_inches
  pdf(paste0(fp, ".pdf"), height = h, width = w)
  print(p)
  dev.off()
  svg(paste0(fp, ".svg"), height = h, width = w)
  print(p)
  dev.off()
  png(paste0(fp, ".png"),
    height = h, width = w,
    units = "in", res = 400
  )
  print(p)
  dev.off()
}

plot_km <- function(df, km_fit, plt) {
  km_plot <- ggsurvplot(km_fit,
    data = df,
    risk.table = TRUE, pval = TRUE, conf.int = FALSE,
    font.tickslab = c(5), fontsize = 2.5, size = 0.5, pval.size = 2.5,
    ggtheme = theme(
      axis.line = element_line(size = 0.5),
      axis.ticks = element_line(size = 0.5),
      axis.ticks.length = unit(.1, "cm"),
      axis.text.x = element_text(size = 6),
      axis.text.y = element_text(size = 6)
    ),
    palette = plt,
    risk.table.y.text.col = TRUE,
    risk.table.y.text = FALSE
  ) +
    labs(x = "Months")
  return(km_plot)
}

plot_paired_boxplot <- function(df, cond1, cond2, ylab, xlab, ylim) {
  p <- ggpaired(df,
    cond1 = cond1, cond2 = cond2,
    width = 0.5, line.color = "grey", outlier.shape = NA,
    alpha = 0.7, size = 0.3
  ) +
    labs(y = ylab, x = xlab) + stat_compare_means(paired = TRUE, size = 2) +
    theme(legend.position = "none") +
    lemon::coord_capped_cart(bottom = "both", left = "bottom")
  if (!any(is.na(ylim))) {
    p <- p + ylim(ylim)
  }

  return(p)
}

plot_tile <- function(df, x_str, y_str, fill_str, lgd) {
  p <- ggplot(df, aes_string(y = y_str, x = x_str, fill = fill_str)) +
    geom_tile(col = "black", size = .5) +
    labs(x = "", y = "", fill = "-log 10 (pval) * sign(effect)") +
    scale_fill_gradient2(
      low = "blue", high = "red", mid = "white",
      midpoint = 0, limit = c(-3, 3)
    ) +
    theme(
      axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 4),
      axis.text.y = element_text(size = 5)
    ) +
    theme(line = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  if (lgd == "no") {
    p <- p + theme(legend.position = "none")
  }
  return(p)
}

plot_volcano <- function(
    df, x_str, y_str, fill_str,
    lab_str, alpha_str, lgd, xlim) {
  p <- ggplot(df, aes_string(x = x_str, y = y_str, fill = fill_str)) +
    geom_point(aes_string(alpha = alpha_str), pch = 21) +
    ggrepel::geom_text_repel(
      min.segment.length = 0.5,
      aes_string(label = lab_str), size = 3
    ) +
    geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") +
    geom_vline(xintercept = 0, col = "black", linetype = "dashed") +
    scale_alpha_manual(values = c("soft" = 0.2, "dark" = 1), guide = "none") +
    theme(
      axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 5),
      axis.text.y = element_text(size = 5)
    ) +
    theme(panel.border = element_blank(), axis.line = element_line()) +
    lemon::coord_capped_cart(bottom = "both", left = "bottom") +
    xlim(xlim)
  if (lgd == "no") {
    p <- p + theme(legend.position = "none")
  }
  return(p)
}


plot_lolli <- function(df, x_str, y_str, fill_str, fill_manual, xlim) {
  p <- ggplot(df, aes_string(x = x_str, y = y_str, fill = fill_str)) +
    geom_segment(aes_string(
      x = 0, y = y_str,
      xend = x_str, yend = y_str
    ), color = "grey50") +
    geom_vline(aes(xintercept = 0), color = "grey35", size = 1) +
    geom_point(size = 3, pch = 21) +
    scale_fill_manual(values = fill_manual) +
    xlim(xlim) +
    theme(panel.border = element_blank(), axis.line = element_line()) +
    lemon::coord_capped_cart(bottom = "left", left = "both") +
    theme(
      axis.text.x = element_text(hjust = 0.5, vjust = 1, size = 5),
      axis.text.y = element_text(size = 5),
      legend.position = "none"
    )
  return(p)
}

plot_perc_variation <- function(af) {
  af$var <- rownames(af)
  af <- arrange(af, perc_var)
  af$var <- factor(af$var, levels = af$var)
  af <- af[af$var != "Residuals", ]
  af$sign <- af$`Pr(>F)` <= .05
  ggplot(af, aes(y = var, x = perc_var, fill = sign)) +
    geom_col(col = "black") +
    labs(y = "", x = "Explained variance of I-TED (%) ") +
    lemon::coord_capped_cart(bottom = "left", left = "both") +
    scale_fill_manual(values = c("FALSE" = "grey80", "TRUE" = "coral"))
}
