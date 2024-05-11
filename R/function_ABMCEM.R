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
                          fixedparameters = c(NA,NA,NA,20,20),
                          formula,formulaRate = NULL, initTime = NULL,
                          endTime = NULL, num_cores = 1, errType2 = 0.8,
                          alpha = 0.9, gamma = 0.9, thr = 1e-3,
                          maxIter = 10000, seqIter = 50,
                          pShort = 0.5, pAug = 0.5, nPT = 1,
                          T0 = 1, seqIter = 50, nStepExch = 10,
                          maxIterPT = 10000,num_coresPT = 1) {


  actDfnodes <- defineNodes(data.frame(label = colnames(net0)))

  nmax <- nmax0

  # Initialization of parameters
  theta <- theta0
  beta <- beta0
  se <- data.frame(Crea = rep(0, nrow(beta)), Del = rep(0, nrow(beta)))
  row.names(se) <- row.names(beta)
  # Creation of sequence of events from initial data
  sequence <- EMPreprocessing(net0, net1)
  H <- nrow(sequence) # Humming distance
  # Creation of permutations
  # permut = permute(sequence, nmax = nmax)
  # permut <- MCMC(nmax = nmax, seq, H, actDfnodes, formula, net0, beta, burnIn = TRUE, maxIter, seqIter, pShort, pAug, pPerm)

  seqsPTinit = permute(sequence, nmax = nPT) # Initial sequences for Parallel Tempering initialization


  diff <- 1000
  index <- 0
  while (diff > thr) {
    # browser()

    cat("Index: ", index, "\n")
    cat("Diff: ", diff, "\n")
    if (index > 1000) {
      cat("No convergence\n")
      break
    }

    # Perform Parallel-Tempering
    seqsPT = PT_MCMC(nmax, nPT, seqsPTinit, H, actDfnodes, formula, net0, beta,
                     theta, fixedparameters, initTime, endTime, burnIn = TRUE,
                     maxIterPT, seqIter, T0 ,nStepExch, num_coresPT,
                     pAug, pShort)


    # c <- max(logLikPrev)
    # logsumExp <- c + log(sum(exp(logLikPrev - c)))
    # weights <- exp(logLikPrev - logsumExp)

    # betaCreaAux <- rubinsRule(betaCreaDF, seCreaDF, weights)
    # betaDelAux <- rubinsRule(betaDelDF, seDelDF, weights)
    # diff= sqrt(norm(as.matrix(beta$Crea-betaCreaAux$mean))^2+norm(as.matrix(beta$Del-betaDelAux$mean))^2)

    # betaNew$Crea <- betaCreaAux$mean
    # betaNew$Del <- betaDelAux$mean


    if (is.null(formulaRate)) {
      logLikCur <- clusterApply(cl, seq_along(splitIndicesPerCore), logLikelihoodMC,
        permut = permut,
        splitIndicesPerCore = splitIndicesPerCore,
        beta = betaNew, actDfnodes = actDfnodes, net0 = net0, formula = formula,
        theta = theta, initTime = initTime, endTime = endTime,
      )
    } else {
      logLikCur <- clusterApply(cl, seq_along(splitIndicesPerCore), logLikelihoodTimeMC,
        permut = permut,
        splitIndicesPerCore = splitIndicesPerCore,
        beta = betaNew, theta = theta, actDfnodes = actDfnodes, net0 = net0,
        formulaChoice = formula, formulaRate = formulaRate,
        initTime = initTime, endTime = endTime
      )
    }
    logLikCur <- unlist(logLikCur)

    w <- rep(1, length(logLikPrev)) / length(logLikPrev)
    sigmaHat <- sigmaHat(nmax, logLikPrev, logLikCur, w)
    ase <- sqrt(sigmaHat / nmax)
    deltaQ <- deltaQ(logLikPrev, logLikCur, w)
    lowerBound <- deltaQ - qnorm(alpha) * ase # 80%



    if (lowerBoud < 0) {
      # Estimator is not accepted, new point must be sampled
      newpermut <- MCMC(nmax = nmax / 2, permut[[length(permut)]], H, actDfnodes, formula, net0, beta, burnIn = FALSE, maxIter, seqIter, pShort, pAug, pPerm)
      nmax <- nmax + (nmax0 / 2)
      permut <- c(permut, newpermut)

      resPar <- clusterApply(cl, seq_along(splitIndicesPerCore), parametersMC,
        permut = permut,
        splitIndicesPerCore = splitIndicesPerCore, actDfnodes. = actDfnodes,
        net0. = net0, formula. = formula
      )

      resPar <- unlist(resPar, recursive = FALSE)
      betaCreaDF <- t(sapply(resPar, "[[", 1))
      betaDelDF <- t(sapply(resPar, "[[", 2))
      seCreaDF <- t(sapply(resPar, "[[", 3))
      seDelDF <- t(sapply(resPar, "[[", 4))

      if (is.null(formulaRate)) {
        logLikPrev <- clusterApply(cl, seq_along(splitIndicesPerCore), logLikelihoodMC,
          permut = permut,
          splitIndicesPerCore = splitIndicesPerCore,
          beta = beta, actDfnodes = actDfnodes, net0 = net0, formula = formula,
          theta = theta, initTime = initTime, endTime = endTime,
        )
      } else {
        logLikPrev <- clusterApply(cl, seq_along(splitIndicesPerCore), logLikelihoodTimeMC,
          permut = permut,
          splitIndicesPerCore = splitIndicesPerCore,
          beta = beta, theta = theta, actDfnodes = actDfnodes, net0 = net0,
          formulaChoice = formula, formulaRate = formulaRate,
          initTime = initTime, endTime = endTime
        )
      }
      logLikPrev <- unlist(logLikPrev)
    } else {
      # Update of the paremeter, next iteration
      index <- index + 1
      diff <- deltaQ + qnorm(gamma) * ase # 90%


      beta <- betaNew
      se$Crea <- betaCreaAux$se
      se$Del <- betaDelAux$se


      # Update on the permutations -> new MCMC

      m_start <- sigmaHat^2 * (qnorm(alpha) + qnorm(errType2))^2 / deltaQ^2

      if (m_start > nmax) {
        nmax <- m_start
      }

      indexMaxlogLik <- which(logLikPrev == max(logLikPrev))

      permut <- MCMC(nmax = nmax, permut[[indexMaxlogLik]], H, actDfnodes, formula, net0, beta, burnIn = TRUE, maxIter, seqIter, pShort, pAug, pPerm)

      resPar <- clusterApply(cl, seq_along(splitIndicesPerCore), parametersMC,
        permut = permut,
        splitIndicesPerCore = splitIndicesPerCore, actDfnodes. = actDfnodes,
        net0. = net0, formula. = formula
      )

      resPar <- unlist(resPar, recursive = FALSE)
      betaCreaDF <- t(sapply(resPar, "[[", 1))
      betaDelDF <- t(sapply(resPar, "[[", 2))
      seCreaDF <- t(sapply(resPar, "[[", 3))
      seDelDF <- t(sapply(resPar, "[[", 4))

      if (is.null(formulaRate)) {
        logLikPrev <- clusterApply(cl, seq_along(splitIndicesPerCore), logLikelihoodMC,
          permut = permut,
          splitIndicesPerCore = splitIndicesPerCore,
          beta = beta, actDfnodes = actDfnodes, net0 = net0, formula = formula,
          theta = theta, initTime = initTime, endTime = endTime,
        )
      } else {
        logLikPrev <- clusterApply(cl, seq_along(splitIndicesPerCore), logLikelihoodTimeMC,
          permut = permut,
          splitIndicesPerCore = splitIndicesPerCore,
          beta = beta, theta = theta, actDfnodes = actDfnodes, net0 = net0,
          formulaChoice = formula, formulaRate = formulaRate,
          initTime = initTime, endTime = endTime
        )
      }
      logLikPrev <- unlist(logLikPrev)
    }
  }


  return(list(
    "logLik" = logLik, "beta" = beta, "se" = se, "index" = index, "diff" = diff,
    "permut" = permut, "betaCreaDF" = betaCreaDF, "betaDelDF" = betaDelDF
  ))
}
