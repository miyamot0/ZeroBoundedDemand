# Shawn Gilroy (2019) - MIT
# Louisiana State University

library(extrafont)
library(jtools)
library(ggplot2)
library(tidyverse)

transMod <- function(x, theta = 0.5) { 
  log((x * theta) + ((theta^2) * (x^2) + 1)^0.5)/log(10)
}

annotation_logticks2 <- function(base = 10, sides = "bl", scaled = TRUE, short = unit(0.1, "cm"),
                                 mid = unit(0.2, "cm"), long = unit(0.3, "cm"), colour = "black",
                                 size = 0.5, linetype = 1, alpha = 1, data =data.frame(x = NA), color = NULL, ...) {
  if (!is.null(color)) {
    colour <- color
  }
  
  layer(
    data = data,
    mapping = NULL,
    stat = StatIdentity,
    geom = GeomLogticks,
    position = PositionIdentity,
    show.legend = FALSE,
    inherit.aes = FALSE,
    params = list(
      base = base,
      sides = sides,
      scaled = scaled,
      short = short,
      mid = mid,
      long = long,
      colour = colour,
      size = size,
      linetype = linetype,
      alpha = alpha,
      ...
    )
  )
}

plottingFrame <- data.frame(
  X     = c(0, seq(0.01, 1000, length.out = 100000)),
  Log10 = c(NA, log10(seq(0.01, 1000, length.out = 100000))),
  IHS   = c(0, transMod(seq(0.01, 1000, length.out = 100000)))
)

par(family="Times")

png(filename = "plots/Figure 1.png",
    width=8,
    height=6,
    units="in",
    res=600)

plottingFrame$mask <- 1

plottingFrame[plottingFrame$X == 0,]$mask <- 0
plottingFrame[plottingFrame$X == 0,]$X <- 0.00001

plottingFrame.melt <- plottingFrame %>%
  gather(Scale, Y, Log10:IHS, -mask)

hackPoint <- data.frame(X = 0.00001, 
                        Y = 0, 
                        Scale = "IHS",
                        mask = 0)

#ggplot2::geom_point(size=3, shape=21, show.legend=T, colour = "black", fill = "white", alpha = .9, stroke = 1) +

plt <- ggplot(plottingFrame.melt, aes(x=X, y=Y, color = Scale)) +
  geom_line() +
  geom_point(hackPoint, mapping = aes(x=X, y=Y, color = Scale), inherit.aes = FALSE) +
  geom_hline(yintercept = 0, lty = 2, ltw = .5) +
  facet_grid(.~mask, scales="free_x", space="free_x") +
  scale_x_log10(breaks=c(0.00001,  0.01, 0.1, 1, 10, 100),
                labels=c("0.00",   0.01, 0.1, 1, 10, 100)) +
  scale_color_manual(values = c(
    "IHS"   = "black",
    "Log10" = "gray"
  )) +
  theme_apa() +
  theme(strip.background = element_blank(),
                 strip.text = element_blank(),
                 plot.title = element_text(hjust = 0.5),
                 legend.position = "bottom",
                 text = element_text(size=16, family="Times")) +
  annotation_logticks2(sides="b", data = data.frame(X= NA, mask = 1)) +
  ylim(c(-2, 2)) +
  labs(x = "Untransformed Value", 
       y = "Transformed Value")

plt <- plt +
  ggplot2::theme(axis.text.x = element_text(size=12, margin = unit(c(0.3,0.3,0.3,0.3), "cm")),
                 axis.text.y = element_text(size=12, margin = unit(c(0.3,0.3,0.3,0.3), "cm")),
                 axis.ticks.length = unit(-0.15, "cm"),
                 axis.title.x = element_text(face = "bold", margin = unit(c(-0.1, 0, 0, 0), "cm")),
                 axis.title.y = element_text(face = "bold", margin = unit(c(0, -0.1, 0, 0), "cm")))

library(grid)
gt = ggplot_gtable(ggplot_build(plt))
gt$widths[5] = 5*gt$widths[5]

plotter <- grid.draw(gt)

print(plotter)

dev.off()
