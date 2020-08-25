# Shawn Gilroy (2019) - MIT
# Louisiana State University

transMod <- function(x, theta = 0.5) {
  log((x * theta) + ((theta^2) * (x^2) + 1)^0.5)/log(10)
}

fromIhsToNormal <- function(x) {
  (1/10^(1*x))*((10^(2*x))-1)
}

derivIHS <- function(x) {
  log(sqrt(x^2 + 17/16) + x / 4) / log(10)
}

calc_stat <- function(x) {
  coef <- 1.5
  n <- sum(!is.na(x))
  # calculate quantiles
  stats <- quantile(x, probs = c(0.01, 0.25, 0.5, 0.75, 0.99))
  names(stats) <- c("ymin", "lower", "middle", "upper", "ymax")
  return(stats)
}


# Par = Q, A
unweighted.EXPD.RSS <- function(data, par) {
  TSS <- 0

  for (i in 1:nrow(data)) {
    yHat <- par[1] * 10^(log(par[1])/log(10) * (exp(-par[2] * par[1] * data$x[i]) - 1))

    diff <- (data$y[i] - yHat)^2

    TSS <- TSS + diff
  }

  TSS
}

SimulateDemand <- function(nruns = 10, setparams, sdindex, x, outdir = NULL, fn = NULL) {
  RunOneSim <- function(run = 1, setparams, sdindex, x) {
    ## Save inSate seed
    inState <- .Random.seed
    ## Initialize vector for simulation parameters
    simparams <- vector(length = 4)
    ## Store simulation parameters and label
    simparams[1] <- rnorm(1, setparams[1], setparams[2])
    simparams[2] <- rnorm(1, setparams[3], setparams[4])
    simparams[3] <- 10^simparams[1]
    simparams[4] <- 10^simparams[2]
    names(simparams) <- c("alphalr", "q0lr", "alphar", "q0r")
    ## Save outState seed
    outState <- .Random.seed
    ## Calculate consumption values at each price
    sim <- function(x, sdindex) {
      m <- 10^(simparams["q0lr"] + setparams[5] *
                 ((exp((-1 * simparams["alphar"]) * simparams["q0r"] * x)) - 1))
      s <- sdindex * setparams[6]
      rnorm(1, m, s)
    }
    ysim <- mapply(sim, x, sdindex)
    ## Round resulting consumption values, truncating at 0
    ysimr <- pmax(round(ysim, 0), 0)
    ## Return list: unrounded and rounded consumption values,
    ## simulation parameters, inState and outState of seed
    list(dat = data.frame("y"  = ysim, "yr" = ysimr),
         "simparams" = simparams, "inState" = inState, "outState" = outState)
  }

  runs <- seq_len(nruns)
  manysims <- lapply(runs, RunOneSim, setparams, sdindex, x)
  simvalues <- sapply(manysims, "[[", 1)
  simparams <- sapply(manysims, "[[", 2)
  seeds <- sapply(manysims, "[", 3:4)

  simvaluesraw <- simvalues[-2, , drop = FALSE]
  simvalues <- simvalues[-1, , drop = FALSE]
  simvaluesraw <- do.call("cbind", simvaluesraw)
  simvalues <- do.call("cbind", simvalues)
  simvaluesraw <- cbind(x, simvaluesraw)
  simvalues <- cbind(x, simvalues)
  colnames(simvaluesraw) <- c("x", runs)
  colnames(simvalues) <- c("x", runs)

  rownames(simparams) <- c("alphalr", "q0lr", "alphar", "q0r")
  colnames(simparams) <- runs
  simparams <- t(simparams)

  colnames(seeds) <- runs
  seeds <- t(seeds)

  ## Reshape long
  simvalues <- data.frame(simvalues)
  simvalues <- reshape2::melt(simvalues, id.vars = "x")
  colnames(simvalues) <- c("x", "id", "y")
  simvalues <- simvalues[, c("id", "x", "y")]
  simvalues$id <- gsub("X", "", simvalues$id)
  simvalues$id <- as.numeric(simvalues$id)

  simvaluesraw <- data.frame(simvaluesraw)
  simvaluesraw <- reshape2::melt(simvaluesraw, id.vars = "x")
  colnames(simvaluesraw) <- c("x", "id", "y")
  simvaluesraw <- simvaluesraw[, c("id", "x", "y")]
  simvaluesraw$id <- gsub("X", "", simvaluesraw$id)
  simvaluesraw$id <- as.numeric(simvaluesraw$id)

  if (!is.null(outdir) && !is.null(fn)) {
    if (!dir.exists(outdir)) {
      dir.create(outdir)
    }
    saveobjects <- c("simvalues", "simvaluesraw", "simparams", "seeds")
    save(list = saveobjects, file = paste0(outdir, fn, ".RData"))
  }

  invisible(list("simvalues" = simvalues, "simvaluesraw" = simvaluesraw,
                 "simparams" = simparams, "seeds" = seeds))
}

lambertW = function(z,b=0,maxiter=10,eps=.Machine$double.eps,min.imag=1e-9) {
  if (any(round(Re(b)) != b))
    stop("branch number for W must be an integer")
  if (!is.complex(z) && any(z<0)) z=as.complex(z)
  ## series expansion about -1/e
  ##
  ## p = (1 - 2*abs(b)).*sqrt(2*e*z + 2);
  ## w = (11/72)*p;
  ## w = (w - 1/3).*p;
  ## w = (w + 1).*p - 1
  ##
  ## first-order version suffices:
  ##
  w = (1 - 2*abs(b))*sqrt(2*exp(1)*z + 2) - 1
  ## asymptotic expansion at 0 and Inf
  ##
  v = log(z + as.numeric(z==0 & b==0)) + 2*pi*b*1i;
  v = v - log(v + as.numeric(v==0))
  ## choose strategy for initial guess
  ##
  c = abs(z + exp(-1));
  c = (c > 1.45 - 1.1*abs(b));
  c = c | (b*Im(z) > 0) | (!Im(z) & (b == 1))
  w = (1 - c)*w + c*v
  ## Halley iteration
  ##
  for (n in 1:maxiter) {
    p = exp(w)
    t = w*p - z
    f = (w != -1)
    t = f*t/(p*(w + f) - 0.5*(w + 2.0)*t/(w + f))
    w = w - t
    if (abs(Re(t)) < (2.48*eps)*(1.0 + abs(Re(w)))
        && abs(Im(t)) < (2.48*eps)*(1.0 + abs(Im(w))))
      break
  }
  if (n==maxiter) warning(paste("iteration limit (",maxiter,
                                ") reached, result of W may be inaccurate",sep=""))
  if (all(Im(w)<min.imag)) w = as.numeric(w)
  return(w)
}

GetAnalyticPmaxFallback <- function(K_, A_, Q0_) {
  result <- NULL
  try(result <- optimx::optimx(par = c((1/(Q0_ * A_ * K_^1.5)) * (0.083 * K_ + 0.65)),
                               fn = function(par, data) {
                                 abs((log((10^data$K)) * (-data$A * data$Q0 * par[1] * exp(-data$A * data$Q0 * par[1]))) + 1)
                               },
                               data = data.frame(Q0 = Q0_,
                                                 A = A_,
                                                 K = K_),
                               method = c("BFGS"),
                               control=list(maxit=2500)), silent = TRUE)

  if (class(result) == "try-error" | is.null(result)) {
    return(NA)
  }

  return(result$p1)
}

GetAnalyticPmax <- function(Alpha, K, Q0) {
  if (K <= exp(1)/log(10)) {
    return (GetAnalyticPmaxFallback(K, Alpha, Q0));
  } else {
    return (-lambertW(z = -1/log((10^K))) / (Alpha * Q0));
  }
}

EmpiricalPmax <- function(dat) {
  tempObj <- dat
  tempObj$PmaxE <- tempObj$x * tempObj$y

  tempMax <- max(tempObj$PmaxE, na.rm = TRUE)
  retObj <- tail(tempObj[tempObj$PmaxE == tempMax,], 1)

  retObj[1, "x"]
}

GetPmaxIHSderivative <- function(A_, Q0_) {
  result <- NULL
  try(result <- optimx::optimx(par = c(1),
                               fn = function(par, data) {

                                 abs((-(data$A*data$Q0*par[1]*exp(-data$A*data$Q0*par[1])*log((data$Q0*(data$Q0^2 + 4)^(1/2))/4)*(par[1]^2 + 4))/(2*(par[1]^2 + 2))) + 1)
                               },
                               data = data.frame(Q0 = Q0_,
                                                 A = A_),
                               method = c("BFGS"),
                               control=list(maxit=2500)), silent = TRUE)

  if (is.null(result)) {
    return(NA)
  }

  return(result$p1)
}

GetPmaxIHSobserved <- function(A_, Q0_) {
  result <- NULL
  try(result <- optimx::optimx(par = c(1),
                               fn = function(par, data) {
                                 demand <- log((data$Q0 * 0.5) + ((0.5^2) * (data$Q0^2) + 1)^0.5)/log(10) + log((data$Q0 * 0.5) + ((0.5^2) * (data$Q0^2) + 1)^0.5)/log(10) * (exp(-data$A * data$Q0 * par[1]) - 1)
                                 demand.n <- fromIhsToNormal(demand)
                                 
                                 -(demand.n * par[1])
                               },
                               data = data.frame(Q0 = Q0_,
                                                 A = A_),
                               method = c("BFGS"),
                               control=list(maxit=2500)), silent = TRUE)
  
  if (is.null(result)) {
    return(NA)
  }
  
  return(result$p1)
}
