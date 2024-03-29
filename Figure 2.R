# Shawn Gilroy (2019) - MIT
# Louisiana State University

library(extrafont)
library(jtools)
library(ggplot2)
library(tidyverse)

cexScale <- 1.5

par(family="Times")

png(filename = "plots/Figure 2.png",
    width=8,
    height=3.5,
    units="in",
    res=600)

dataFrameMod <- data.frame(
  x = c(0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,10.0,15.0, 20, 50, 100), 
  y = c(1000,1000,800,800,700,600,500,400,200,100, 50, 10, 0.1)
)

dataFrameMod2 <- data.frame(
  x = c(0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,10.0,15.0, 20, 50, 100), 
  y = c(1000,1000,800,800,700,600,500,400,200,100, 50, 10, 0.01)
)

dataFrameKeep <- data.frame(
  x = c(0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,10.0,15.0, 20, 50, 100), 
  y = c(1000,1000,800,800,700,600,500,400,200,100, 50, 10, 0.0)
)

xAxisLimits <- c(.1, 100)

kMod <-  transMod(max(dataFrameMod$y)) -  transMod(min(dataFrameMod$y))  + 0.5
hurshMod <- nls(transMod(y) ~ transMod(q0) + kMod * (exp(-alpha * q0 * x) - 1),
                start = list(
                  q0 = 1000, 
                  alpha = 0.00006
                ),
                data = dataFrameMod)

kMod2 <-  transMod(max(dataFrameMod2$y)) -  transMod(min(dataFrameMod2$y))  + 0.5
hurshMod2 <- nls(transMod(y) ~ transMod(q0) + kMod2 * (exp(-alpha * q0 * x) - 1),
                 start = list(
                   q0 = 1000, 
                   alpha = 0.00006
                 ),
                 data = dataFrameMod2)

kKeep <-  transMod(max(dataFrameKeep$y)) -  transMod(min(dataFrameKeep$y))  + 0.5
hurshKeep <- nls(transMod(y) ~ transMod(q0) + kKeep * (exp(-alpha * q0 * x) - 1),
                 start = list(
                   q0 = 1000, 
                   alpha = 0.00006
                 ),
                 data = dataFrameKeep)

par(mfrow = c(1, 3),
    family="Times")

plot(dataFrameKeep$x, transMod(dataFrameKeep$y),
     main = "Zero Included",
     ylim = c(0, 4),
     xlim = xAxisLimits,
     ylab = "Transformed Consumption",
     xlab = "",
     xaxt = 'n',
     log = "x",
     cex.main = cexScale,
     cex.axis = cexScale,
     cex.lab = cexScale,
     cex = cexScale)
axis(1, 
     at=c(0.1, 1, 10, 100), 
     labels=c(0.1, 1, 10, 100),
     cex.axis = cexScale)
lines(dataFrameKeep$x, predict(hurshKeep))

text(.1, 0, 
     paste0("Q0: ", round(coef(hurshKeep)["q0"], 2), "\n",
            "Alpha: ", round(coef(hurshKeep)["alpha"], 6)
     ),
     adj = c(0, 0),
     cex = cexScale
)

plot(dataFrameMod$x, transMod(dataFrameMod$y),
     main = "Zero Replaced with 0.1",
     ylim = c(0, 4),
     xlim = xAxisLimits,
     ylab = "",
     xlab = "Price",
     xaxt = 'n',
     log = "x",
     cex.main = cexScale,
     cex.axis = cexScale,
     cex.lab = cexScale,
     cex = cexScale)
axis(1, 
     at=c(0.1, 1, 10, 100), 
     labels=c(0.1, 1, 10, 100),
     cex.axis = cexScale)
lines(dataFrameMod$x, predict(hurshMod))

text(.1, 0, 
     paste0("Q0: ", round(coef(hurshMod)["q0"], 2), "\n",
            "Alpha: ", round(coef(hurshMod)["alpha"], 6)
     ),
     adj = c(0, 0),
     cex = cexScale
)

plot(dataFrameMod2$x, transMod(dataFrameMod2$y),
     main = "Zero Replaced with 0.01",
     ylim = c(0, 4),
     xlim = xAxisLimits,
     ylab = "",
     xlab = "",
     xaxt = 'n',
     log = "x",
     cex.main = cexScale,
     cex.axis = cexScale,
     cex = cexScale)
axis(1, 
     at=c(0.1, 1, 10, 100), 
     labels=c(0.1, 1, 10, 100),
     cex.axis = cexScale)
lines(dataFrameMod2$x, predict(hurshMod2))

text(.1, 0, 
     paste0("Q0: ", round(coef(hurshMod2)["q0"], 2), "\n",
            "Alpha: ", round(coef(hurshMod2)["alpha"], 6)
     ),
     adj = c(0, 0),
     cex = cexScale
)

dev.off()