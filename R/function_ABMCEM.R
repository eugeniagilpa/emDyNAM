# Ascent-based MCEM -----------------------------


#' Auxiliar function to compute \enq{\hat\sigma}.
#'
#' @param mt number of observations
#' @param loglikPrev log-likelihood at time \enq{t-1}
#' @param loglikCur log-likelihood at time \enq{t}
#' @param w weights (for importance sampling)
#'
#' @return estimator of \enq{\sigma^2}
#' @export
#'
#' @references Caffo, B. S., Jank, W., & Jones, G. L. (2005). Ascent-based Monte Carlo expectation–maximization. \emph{Journal of the Royal Statistical Society Series B: Statistical Methodology, 67(2)}, 235-251.
#'
sigmaHat <- function(mt, loglikPrev, loglikCur, w) {
  lambda <- loglikCur - loglikPrev
  sumlambdaw <- sum(w * lambda)
  sumlambdaw2 <- sum((w * lambda)^2)
  sumw2lambda <- sum((w^2) * lambda)

  sigma2 <- mt * sumlambdaw^2 * (sumlambdaw2 / sumlambdaw^2 - 2 * sumw2lambda / sumlambdaw + sum(w^2))
  return(sigma2)
}

#' Auxiliar function to compute \enq{\tilde{\Delta Q}}.
#'
#' @param loglikPrev log-likelihood at time \enq{t-1}
#' @param loglikCur log-likelihood at time \enq{t}
#' @param w weights (for importance sampling)
#'
#' @return estimator of \enq{\Delta Q}
#' @export
#'
#' @references Caffo, B. S., Jank, W., & Jones, G. L. (2005). Ascent-based Monte Carlo expectation–maximization. \emph{Journal of the Royal Statistical Society Series B: Statistical Methodology, 67(2)}, 235-251.
#'
deltaQ <- function(loglikPrev, loglikCur, w) {
  loglikRatio <- loglikCur - loglikPrev
  deltaQ <- sum(w * loglikRatio)

  return(deltaQ)
}


# Ascent-Based Markov-Chain Expectation-Maximization algorithm --------------


#' Ascent-Based Markov-Chain Expectation-Maximization algorithm
#'
#' @description This algorithm allows for an implementation of the ABMCEM algorithm with and without constant rates.
#'              For the MCMC steps, three types of steps are done: shortening, augmenting and permutation of the sequence.
#'              \eqn{pShort + pAug + pPerm = 1}
#'
#' @param nmax0 initial value of number of permuations for step 0
#' @param net0 initial network
#' @param net1 final network
#' @param theta0 initial values for rate parameters
#' @param beta0 initial values for effect parameters
#' @param formula formula of the choice model
#' @param formulaRate formula of the rate model. Default value is `NULL`. Leave default for constant rate models (faster algorithm)
#' @param initTime initial time. If formulaRate is not `NULL`, this parameter cannot be null.
#' @param endTime end time. If formulaRate is not `NULL`, this parameter cannot be null.
#' @param num_cores number of cores for paralelization
#' @param errType2 value of power wanted. Default 0.8
#' @param alpha value for lower bound computation. Default 0.9
#' @param gamma value for stopping rule. Default 0.9
#' @param thr threshold for algorithm stop
#' @param maxIter maximum number of iterations for MCMC
#' @param seqIter number of MCMC steps beween samples
#' @param pShort probability of shortening the sequence (MCMC)
#' @param pAug probability of augmenting the sequence (MCMC)
#' @param pPerm probability of permuting two elements of the sequence (MCMC)
#'
#'
#' @return list of
#' \item{logLik}{log likelihood of the sequences at the last step}
#' \item{beta}{estimation of effect parameters}
#' \item{se}{estimation of standard errors}
#' \item{index}{number of iterations until convergence of EM algorithm}
#' \item{diff}{value of stopping rule at last iteration}
#' \item{permut}{last set of permutations}
#' \item{betaCreaDF}{last data frame of goldfish estimators}
#' \item{betaDelDF}{last data frame of goldfish estimators}
#' @export
#'
#' @references
#' Caffo, B. S., Jank, W., & Jones, G. L. (2005). Ascent-based Monte Carlo expectation–maximization.
#' \emph{Journal of the Royal Statistical Society Series B: Statistical Methodology, 67(2)}, 235-251.
#'
#' Stadtfeld, C., and Block, P. (2017). Interactions, Actors, and Time:
#' Dynamic Network Actor Models for Relational Events.
#' \emph{Sociological Science 4 (1)}, 318-52. \doi{10.15195/v4.a14}
#'
#' Koskinen J.(2004). BAYESIAN INFERENCE FOR LONGITUDINAL SOCIAL NETWORKS.
#' \emph{Research Report, number 2004: 4.}.
#'
#' Ruth, W. (2024). A review of Monte Carlo-based versions of the EM algorithm.
#' \emph{arXiv preprint} arXiv:2401.00945.
#'
#' Snijders, T. A., Koskinen, J., & Schweinberger, M. (2010).
#' Maximum likelihood estimation for social network dynamics.
#' \emph{The annals of applied statistics, 4(2)}, 567.
#'
MCEMalgorithm <- function(nmax0, net0, net1, theta0, beta0,
                          fixedparameters = c(NA, NA, NA, 20, 20),
                          formula, formulaRate = NULL, initTime = NULL,
                          endTime = NULL, num_cores = 1, errType2 = 0.8,
                          alpha = 0.9, gamma = 0.9, thr = 1e-3,
                          maxIter = 1000, seqIter = 50,
                          pShort = 0.5, pAug = 0.5, nPT = 1,
                          T0 = 1, nStepExch = 10,
                          maxIterPT = 1000) {
  actDfnodes <- defineNodes(data.frame(label = colnames(net0)))

  nmax <- nmax0

  # Initialization of parameters
  theta <- theta0
  beta <- beta0
  fixedparameters <- data.frame(
    "Crea" = c(NA, NA, NA, NA, -20, -20),
    "Del" = c(NA, NA, NA, NA, 20, 20)
  )
  se <- data.frame(Crea = rep(0, nrow(beta)), Del = rep(0, nrow(beta)))
  row.names(se) <- row.names(beta)
  # Creation of sequence of events from initial data
  sequence <- EMPreprocessing(net0, net1)
  H <- nrow(sequence) # Humming distance

  acceptanceRates <- data.frame(
    "Acc Aug" = 0, "Total Aug" = 0,
    "Acc Short" = 0, "Total Short" = 0,
    "Acc Perm" = 0, "Total Perm" = 0,
    "Temp" = 0,
    "Acc Switch" = 0, "Total Swithc" = 0,
    "Temps Switch" = 0
  )


  # Creation of permutations
  # permut = permute(sequence, nmax = nmax)
  # permut <- MCMC(nmax = nmax, seq, H, actDfnodes, formula, net0, beta, burnIn = TRUE, maxIter, seqIter, pShort, pAug, pPerm)

  cl <- makeCluster(num_cores)
  on.exit(stopCluster(cl))
  clusterExport(cl, c(
    "actDfnodesLab", "tieNames",
    "formula", "net0", "beta", "theta", "initTime",
    "endTime", "k", "T0", "nStepExch", "pAug", "pShort", "splitIndicesPerCore",
    "H", "actDfnodesLab", "actDfnodes", "nAct", "fixedparameters"
  ))
  clusterEvalQ(cl, {
    library(goldfish)
    # library(matrixStats) # CHEcKEAR
    NULL
  })

  clusterExport(cl, list(
    "stepPT", "stepPTMC", "stepMCMC", "getpDoAugment", "getpDoShort",
    "getpDoPerm", "getKelMeMatrix", "getAuxDfE", "stepAugment", "stepShort",
    "stepPerm", "getlogLikelihood", "GatherPreprocessingDF",
    "sampleVec", "getlogLikelihoodMC"
  ))



  seqsPT <- permute(seq, nmax = nPT) # Initial sequences for Parallel Tempering initialization

  dampingIncreaseFactor <- 1
  dampingDecreaseFactor <- 3
  logLikPrevCrea <- -Inf
  logLikPrevDel <- -Inf

  diff <- 1000
  index <- 0
  while (diff > thr) {
    # browser()
    k <- 2

    cat("Index: ", index, "\n")
    cat("Diff: ", diff, "\n")
    if (index > 1000) {
      cat("No convergence\n")
      break
    }

    # Perform Parallel-Tempering
    # for(i in 1:10){
    #   set.seed(i)
    #   print(i)

    seqsEM <- PT_MCMC(nmax, nPT, seqsPT, H, actDfnodes, formula, net0, beta,
      theta, fixedparameters, initTime, endTime,
      burnIn = FALSE,
      burnInIter = 0, maxIterPT, seqIter, T0, nStepExch,
      pAug, pShort, cl, num_cores
    )
    # }
    seqsPT <- unlist(lapply(seqsEM$resstepPT, "[", "newseq"), recursive = FALSE) # update sequences to start PT
    # Compute new estimator from the sequences





    newNRstep <- newtonraphsonStep(
      parameters = beta,
      fixedparameters = fixedparameters,
      formula = formula,
      seqsEM = seqsEM$seqsEM,
      dampingFactorCrea = dampingFactorCrea,
      dampingFactorDel = dampingFactorDel
    )


    splitIndicesPerCore <- splitIndices(length(seqsEM), num_cores)

    logLikCur <- clusterApply(cl, seq_along(splitIndicesPerCore), getlogLikelihoodMC,
      seqsEM = seqsEM$seqsEM, beta = newNRstep$parameters,
      fixedparameters = fixedparameters,
      splitIndicesPerCore = splitIndicesPerCore,
      initTime = initTime, endTime = endTime,
      actDfnodes = actDfnodes, net0 = net0,
      formula = formula, temp = rep(1, length(seqsEM))
    )


    logLikCur <- unlist(logLikCur, recursive = FALSE)

    logLikCurEM <- sapply(logLikCur, sum)

    w <- rep(1, length(logLikPrev)) / length(logLikPrev)
    sigmaHat <- sigmaHat(nmax, logLikPrev, logLikCurEM, w)
    ase <- sqrt(sigmaHat / nmax)
    deltaQ <- deltaQ(logLikPrev, logLikCurEM, w)
    lowerBound <- deltaQ - qnorm(alpha) * ase # 80%



    if (lowerBoud < 0) {
      while (lowerBound < 0) {
        # Estimator is not accepted, new point must be sampled
        # Perform Parallel-Tempering
        nmax <- nmax + nmax / k
        seqsEMaux <- PT_MCMC(nmax / k, nPT, seqsPT, H, actDfnodes, formula, net0, beta,
          theta, fixedparameters, initTime, endTime,
          burnIn = FALSE,
          burnInIter = 500, maxIterPT, seqIter, T0, nStepExch,
          pAug, pShort, cl, num_cores
        )

        seqsEM <- c(seqsEM, seqsEMaux)
        # Compute new estimator from the sequences

        newNRstep <- newtonraphsonStep(
          parameters.old = beta,
          fixedparameters = fixedparameters,
          formula = formula,
          seqsEM = seqsEM$seqsEM, # TO DO: CHANGE THE OUT PUT OF SEQSem
          dampingIncreaseFactor = 2,
          dampingDecreaseFactor = 3
        )


        splitIndicesPerCore <- splitIndices(length(seqsEM), num_cores)

        logLikCur <- clusterApply(cl, seq_along(splitIndicesPerCore), getlogLikelihoodMC,
          seqsEM = seqsEM, beta = newNRstep$parameters,
          fixedparameters = fixedparameters,
          splitIndicesPerCore = splitIndicesPerCore,
          initTime = initTime, endTime = endTime,
          actDfnodes = actDfnodes, net0 = net0,
          formula = formula, temp = rep(1, length(seqsEM))
        )


        logLikCur <- unlist(logLikCur, recursive = FALSE)

        logLikCurEM <- sapply(logLikCur, sum)


        w <- rep(1, length(logLikPrev)) / length(logLikPrev)
        sigmaHat <- sigmaHat(nmax, logLikPrev, logLikCurEM, w)
        ase <- sqrt(sigmaHat / nmax)
        deltaQ <- deltaQ(logLikPrev, logLikCurEM, w)
        lowerBound <- deltaQ - qnorm(alpha) * ase # 80%

        k <- k + 1
      }
    } else {
      # Update of the paremeter, next iteration
      index <- index + 1
      diff <- deltaQ + qnorm(gamma) * ase # 90%


      beta <- newNRstep$parameters
      se <- newNRstep$stdErrors


      # Update on the number of permutations
      m_start <- sigmaHat^2 * (qnorm(alpha) + qnorm(errType2))^2 / deltaQ^2

      if (m_start > nmax) {
        nmax <- ceiling(m_start)
      }

      logLikCurReduce <- Reduce("+", logLikCur)

      if (logLikCurReduce[1] <= logLikPrevCrea) {
        dampingFactorCrea <- dampingFactorCrea * dampingIncreaseFactor
      } else {
        dampingFactorCrea <- max(
          1,
          dampingFactorCrea / dampingDecreaseFactor
        )
      }
      # logLikPrevCrea <- logLikCurReduce[1]

      if (logLikCurReduce[2] <= logLikPrevDel) {
        dampingFactorDel <- minDampingFactorDel * dampingIncreaseFactor
      } else {
        dampingFactorDel <- max(
          1,
          dampingFactorDel / dampingDecreaseFactor
        )
      }

      # logLikPrevDel <- logLikCurReduce[2]
      #
      # logLikPrev <- logLikCurEM
    }
  }


  return(list(
    "logLik" = logLikCur, "beta" = beta, "se" = se, "index" = index, "diff" = diff
  ))
}
