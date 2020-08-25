# Shawn Gilroy (2019) - MIT
# Louisiana State University

library(beezdemand)
library(broom)
library(ggplot2)
library(grid)
library(gridExtra)
library(kableExtra)
library(knitr)
library(tidyverse)

##

set.seed(65535)

png(filename = "plots/Figure 3.png",
    width=8,
    height=3.5,
    units="in",
    res=600)

nSimulatedSeries <- 2000

capTo1000 <- TRUE

source("helperFx.R")

##

sdindex <- c(2.1978, 1.9243, 1.5804, 1.2465, 0.8104, 0.1751, 0.0380, 0.0270)
x <- c(.1, 1, 3, 10, 30, 100, 300, 1000)

setparams <- vector(length = 4)
setparams <- c(-2.5547,
               .702521,
               1.239893,
               .320221,
               3.096,
               1.438231)
names(setparams) <- c("alphalm", "alphalsd", "q0lm", "q0lsd", "k", "yvalssd")

sim <- SimulateDemand(nruns = nSimulatedSeries,
                      setparams = setparams,
                      sdindex = sdindex,
                      x = x)

sim$simparams <- sim$simparams %>%
  as.data.frame() %>%
  rownames_to_column(var = "id")

passingIds <- beezdemand::CheckUnsystematic(sim$simvalues) %>%
  filter(TotalPass > 2) %>%
  select(id) %>%
  pull()

passingFrame <- sim$simvalues %>%
  filter(id %in% passingIds)

mColNames <- c("id", "Q0", "Alpha",
               "Q0E",
               "PmaxE",
               "Q0.EXPL",
               "Alpha.EXPL",
               "Pmax.EXPL",
               "Q0.EXPD",
               "Alpha.EXPD",
               "Pmax.EXPD",
               "Q0.IHS",
               "Alpha.IHS",
               "Pmax.IHS")

coreFrame <- data.frame(matrix(ncol = length(mColNames), nrow = 0))
colnames(coreFrame) <- mColNames

highestId <- max(passingFrame$id)

for (i in passingIds) {
  test <- subset(passingFrame, id == i)
  noZero.currentData <- subset(test, y > 0)

  if (nrow(noZero.currentData) < 3) {
    message("Skipping, too few")

    next;
  } else {
    message(paste0("Current ID# ", i, "; Max ID# ", highestId))
  }

  steveModFit <- NULL
  koffModFit  <- NULL
  zbeFit2     <- NULL

  try(steveModFit <- nlmrt::wrapnls(
    formula = log(y)/log(10) ~ log(q0)/log(10) + log(q0)/log(10) * (exp(-alpha * q0 * x) - 1),
    start   = list(q0    = 10^setparams["q0lm"] %>% unname(),
                   alpha = 10^setparams["alphalm"] %>% unname()),
    lower   = c(0, 0),
    upper   = c(Inf, Inf),
    data    = noZero.currentData), silent = TRUE)

  try(koffModFit <- optim(par  = c(10^setparams["q0lm"] %>% unname(),
                                   10^setparams["alphalm"] %>% unname()),
                          fn   = unweighted.EXPD.RSS,
                          data = test))

  try(zbeFit2 <- nlmrt::wrapnls(
    formula = log((y * 0.5) + ((0.5^2) * (y^2) + 1)^0.5)/log(10) ~ log((q0 * 0.5) + ((0.5^2) * (q0^2) + 1)^0.5)/log(10) +
      log((q0 * 0.5) + ((0.5^2) * (q0^2) + 1)^0.5)/log(10) * (exp(-alpha * q0 * x) - 1),
    start   = list(q0    = max(test$y),
                   alpha = 0.01),
    lower   = c(0, 0),
    upper   = c(Inf, Inf),
    data    = test), silent = TRUE)


  if (!is.null(koffModFit) & !is.null(steveModFit) & !is.null(zbeFit2)) {
    ### Seed Estimates
    coreFrame[i, "Q0"]     <- 10^sim$simparams[sim$simparams$id == as.numeric(i), "q0lr"]
    coreFrame[i, "Alpha"]  <- 10^sim$simparams[sim$simparams$id == i, "alphalr"]

    ### Empirical Estimates
    coreFrame[i, "Q0E"] <- test[1, "y"]
    coreFrame[i, "id"] <-  test[1, "id"]

    coreFrame[i, "PmaxE"] <- EmpiricalPmax(test)

    coreFrame[i, "Alpha.EXPL"]  <- coef(steveModFit)["alpha"]
    coreFrame[i, "Q0.EXPL"]     <- coef(steveModFit)["q0"]
    coreFrame[i, "Pmax.EXPL"]  <- GetPmaxEXPLobserved(
      coreFrame[i, "Alpha.EXPL"],
      coreFrame[i, "Q0.EXPL"]
    )

    coreFrame[i, "Alpha.EXPD"] <- koffModFit$par[2]
    coreFrame[i, "Q0.EXPD"]    <- koffModFit$par[1]
    coreFrame[i, "Pmax.EXPD"]  <- GetPmaxEXPDobserved(
      coreFrame[i, "Alpha.EXPD"],
      coreFrame[i, "Q0.EXPD"]
    )

    coreFrame[i, "Alpha.IHS"]  <- coef(zbeFit2)["alpha"]
    coreFrame[i, "Q0.IHS"]     <- coef(zbeFit2)["q0"]

    coreFrame[i, "Pmax.IHS"]  <- GetPmaxIHSobserved(
      coef(zbeFit2)["alpha"] %>% unname(),
      coef(zbeFit2)["q0"] %>% unname()
    )

  } else {
    if (is.null(zbeFit2)) {
      message('zbe fail')
    }
  }
}

if (capTo1000 == TRUE) {
  coreFrame <- coreFrame %>%
    head(1000)
}

###
# Just pull from here if just re-drawing
###

################
### Figure 3 ###
################

coreFrame %>%
  select(Q0E:Pmax.IHS) %>%
  select(-c(Alpha.IHS, Alpha.EXPL, Alpha.EXPD)) %>%
  rename(Q0.ZBE = Q0.IHS,
         Pmax.ZBE = Pmax.IHS) %>%
  gather(Parameter, Value) %>%
  mutate(Model     = str_extract(Parameter, "[^.]+$"),
         Parameter = str_extract(Parameter, "[^.]*"),
         Parameter = ifelse(Parameter == "Q0E", "Q0", Parameter),
         Parameter = ifelse(Parameter == "PmaxE", "Pmax", Parameter),
         Model = ifelse(Model == "Q0E" | Model == "PmaxE", "Empirical", Model)) %>%
  ggplot(., aes(Model, Value)) +
  ggtitle(paste0("Koffarnus et al. (2015) Simulation (n=", nrow(coreFrame), ")")) +
  stat_summary(fun.data = calc_stat, geom="boxplot") +
  coord_flip() +
  scale_y_log10(breaks = c(0.1, 1, 10, 100)) +
  scale_x_discrete(limits=c("Empirical",
                            "ZBE",
                            "EXPL",
                            "EXPD")) +
  facet_wrap(~Parameter, scales = "free_y", ncol = 2) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

dev.off()

####################
### Correlations ###
####################

cor.test(coreFrame$Pmax.IHS,  coreFrame$PmaxE,  method = "spearman", exact = FALSE) %>%
  tidy()

cor.test(coreFrame$Pmax.EXPL,  coreFrame$PmaxE, method = "spearman", exact = FALSE) %>%
  tidy()

cor.test(coreFrame$Pmax.EXPD,  coreFrame$PmaxE, method = "spearman", exact = FALSE) %>%
  tidy()

####

cor.test(coreFrame$Q0.IHS,  coreFrame$Q0E, method = "spearman", exact = FALSE) %>%
  tidy()

cor.test(coreFrame$Q0.EXPL,  coreFrame$Q0E, method = "spearman", exact = FALSE) %>%
  tidy()

cor.test(coreFrame$Q0.EXPD,  coreFrame$Q0E, method = "spearman", exact = FALSE) %>%
  tidy()

