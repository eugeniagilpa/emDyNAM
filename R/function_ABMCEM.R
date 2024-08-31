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
#' @param thin number of MCMC steps beween samples
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
                          maxIter = 1000, thin = 50,
                          pShort = 0.5, pAug = 0.5, k=5 ,nPT = 1,
                          T0 = 1, nStepExch = 10,
                          maxIterPT = 1000) {
  actDfnodes <- defineNodes(data.frame(label = colnames(net0)))

  nmax <- nmax0

  # Initialization of parameters
  theta <- theta0

  fixedparameters <- data.frame(
    "Crea" = c(NA, NA, NA, NA, -20, -20),
    "Del" = c(NA, NA, NA, NA, 20, 20)
  )
  se <- data.frame(Crea = rep(0, nrow(beta)), Del = rep(0, nrow(beta)))
  row.names(se) <- row.names(beta)
  # Creation of sequence of events from initial data
  seq <- EMPreprocessing(net0, net1)
  H <- nrow(seq) # Humming distance

  if (any(!beta0 == 0)) {
    beta <- beta0
  } else { # all initial parameters are zero

    parametersAux <- parameters(seq,
      actDfnodes. = actDfnodes,
      net0. = net0, formula. = formula
    )

    beta <- data.frame(
      "Crea" = parametersAux$betaCreaDF,
      "Del" = parametersAux$betaDelDF
    )
  }



  # Creation of permutations
  # permut = permute(sequence, nmax = nmax)
  # permut <- MCMC(nmax = nmax, seq, H, actDfnodes, formula, net0, beta, burnIn = TRUE, maxIter, thin, pShort, pAug, pPerm)

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
  dampingDecreaseFactor <- 1
  dampingFactorCrea <- 1
  dampingFactorDel <- 1


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
      burnInIter = 0, maxIterPT, thin, T0, nStepExch,
      pAug, pShort, pPerm, k,cl, num_cores
    )
    # }
    seqsPT <- unlist(seqsEM$resstepPT, recursive = FALSE)
    # unlist(lapply(seqsEM$resstepPT, "[", "newseq"), recursive = FALSE) # update sequences to start PT

    logLikPrevCrea <- sapply(lapply(
      lapply(
        seqsEM$seqsEM, "[[",
        "newlogLikelihoodStats"
      ), "[[",
      "resCrea"
    ), "[[", "logLikelihood")
    logLikPrevDel <- sapply(lapply(
      lapply(
        seqsEM$seqsEM, "[[",
        "newlogLikelihoodStats"
      ), "[[",
      "resDel"
    ), "[[", "logLikelihood")
    logLikPrev <- logLikPrevCrea + logLikPrevDel # vector of likelihoods lambda^{t-1}
    # Compute new estimator from the sequences

    newNRstep <- newtonraphsonStep(
      parameters = beta,
      fixedparameters = fixedparameters,
      formula = formula,
      seqsEM = seqsEM$seqsEM,
      dampingFactorCrea = dampingFactorCrea,
      dampingFactorDel = dampingFactorDel
    )


    splitIndicesPerCore <- splitIndices(length(seqsEM$seqsEM), num_cores)

    logLikCur <- clusterApplyLB(cl, seq_along(splitIndicesPerCore), getlogLikelihoodMC,
      seqsEM = seqsEM$seqsEM, beta = newNRstep$parameters,
      fixedparameters = fixedparameters,
      splitIndicesPerCore = splitIndicesPerCore,
      initTime = initTime, endTime = endTime,
      actDfnodes = actDfnodes, net0 = net0,
      formula = formula, temp = rep(1, length(seqsEM$seqsEM))
    )


    logLikCur <- unlist(logLikCur, recursive = FALSE)

    logLikCurEM <- sapply(logLikCur, sum) # total likelihoods for each sequence (Crea+Del)

    w <- rep(1, length(logLikPrev)) / length(logLikPrev)
    sigmaHat <- sigmaHat(nmax, logLikPrev, logLikCurEM, w)
    ase <- sqrt(sigmaHat / nmax)
    deltaQ <- deltaQ(logLikPrev, logLikCurEM, w)
    lowerBound <- deltaQ - qnorm(alpha) * ase # 80%



    if (lowerBound < 0) {
      while (lowerBound < 0) {
        # Estimator is not accepted, new point must be sampled
        # Perform Parallel-Tempering
        nmax <- ceiling(nmax + nmax / k)
        seqsEMaux <- PT_MCMC(nmax / k, nPT, seqsPT, H, actDfnodes, formula, net0, beta,
          theta, fixedparameters, initTime, endTime,
          burnIn = FALSE,
          burnInIter = 500, maxIterPT, thin, T0, nStepExch,
          pAug, pShort, pPerm, k,cl, num_cores
        )

        # seqsEM <- c(seqsEM, seqsEMaux) not possible!!
        seqsEM$seqsEM <- c(seqsEM$seqsEM, seqsEMaux$seqsEM)
        seqsEM$resstepPT <- c(seqsEM$resstepPT, seqsEMaux$resstepPT)
        seqsEM$acceptSwitch <- rbind(seqsEM$acceptSwitch, seqsEMaux$acceptSwitch)
        seqsEM$acceptDF <- rbind(seqsEM$acceptDF, seqsEMaux$acceptDF)
        seqsEM$mcmcDiagDF <- rbind(seqsEM$mcmcDiagDF, seqsEMaux$mcmcDiagDF)

        # Compute new estimator from the sequences

        seqsPT <- unlist(seqsEM$resstepPT, recursive = FALSE)

        logLikPrevCrea <- c(
          logLikPrevCrea,
          sapply(lapply(
            lapply(
              seqsEMaux$seqsEM, "[[",
              "newlogLikelihoodStats"
            ), "[[",
            "resCrea"
          ), "[[", "logLikelihood")
        )
        logLikPrevDel <- c(
          logLikPrevDel,
          sapply(lapply(
            lapply(
              seqsEMaux$seqsEM, "[[",
              "newlogLikelihoodStats"
            ), "[[",
            "resDel"
          ), "[[", "logLikelihood")
        )


        logLikPrev <- logLikPrevCrea + logLikPrevDel # vector of likelihoods lambda^{t-1}

        # Compute new estimator from the sequences
        newNRstep <- newtonraphsonStep(
          parameters = beta,
          fixedparameters = fixedparameters,
          formula = formula,
          seqsEM = seqsEM$seqsEM,
          dampingFactorCrea,
          dampingFactorDel
        )


        splitIndicesPerCore <- splitIndices(length(seqsEM$seqsEM), num_cores)

        logLikCur <- clusterApplyLB(cl, seq_along(splitIndicesPerCore), getlogLikelihoodMC,
          seqsEM = seqsEM$seqsEM, beta = newNRstep$parameters,
          fixedparameters = fixedparameters,
          splitIndicesPerCore = splitIndicesPerCore,
          initTime = initTime, endTime = endTime,
          actDfnodes = actDfnodes, net0 = net0,
          formula = formula, temp = rep(1, length(seqsEM$seqsEM))
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

      # logLikCurReduce <- Reduce("+", logLikCur)

      # if (logLikCurReduce[1] <= logLikPrevCrea) {
      #   dampingFactorCrea <- dampingFactorCrea * dampingIncreaseFactor
      # } else {
      #   dampingFactorCrea <- max(
      #     1,
      #     dampingFactorCrea / dampingDecreaseFactor
      #   )
      # }
      # # logLikPrevCrea <- logLikCurReduce[1]
      #
      # if (logLikCurReduce[2] <= logLikPrevDel) {
      #   dampingFactorDel <- minDampingFactorDel * dampingIncreaseFactor
      # } else {
      #   dampingFactorDel <- max(
      #     1,
      #     dampingFactorDel / dampingDecreaseFactor
      #   )
      # }
    }
  }

  acceptSwitch <- table(seqsEM$acceptSwitch)
  acceptDF <- table(seqsEM$acceptDF)

  mcmcObj <- tapply(seqsEM$mcmcDiagDF[, -ncol(seqsEM$mcmcDiagDF)],
    seqsEM$mcmcDiagDF$temp, mcmc,
    thin = thin
  )

  geweke <- lapply(mcmcObj, geweke.diag)
  ess <- lapply(mcmcObj, effectiveSize)

  return(list(
    "logLik" = logLikCur, "beta" = beta, "se" = se,
    "index" = index, "diff" = diff, "acceptSwitch" = acceptSwitch,
    "acceptDF" = acceptDF, "geweke" = geweke, "ess" = ess,
    "acceptSwitchDF" = seqsEM$acceptSwitch
  ))
}





#' Ascent-Based Markov-Chain Expectation-Maximization algorithm with parallel tempering
#'
#' @description This algorithm allows for an implementation of the ABMCEM algorithm with constant rates.
#'              For the MCMC steps in each chain, three types of mutations are done: shortening, augmenting and permutation of the sequence.
#'              \eqn{pShort + pAug + pPerm = 1}
#'
#' @param nmax initial value of number of permuations for step 0
#' @param net0 initial network
#' @param net1 final network
#' @param theta0 initial values for rate parameters
#' @param beta0 initial values for effect parameters
#' @param formula formula of the choice model
#' @param initTime initial time.
#' @param endTime end time.
#' @param num_cores number of cores for paralelization
#' @param errType2 value of power wanted. Default 0.8
#' @param alpha value for lower bound computation. Default 0.9
#' @param gamma value for stopping rule. Default 0.9
#' @param thr threshold for algorithm stop
#' @param maxIter maximum number of iterations for Parallel Tempering
#' @param thin number of mutations steps between samples
#' @param pShort probability of shortening the sequence, Parallel Tempering
#' @param pAug probability of augmenting the sequence, Parallel Tempering
#' @param pPerm probability of permuting two elements of the sequence, Parallel Tempering
#' @param k_perm number of permutations after an augmentation, shortening or permutation mutation.
#' @param nPT number of number of parallel chains with different temperatures.
#' @param T0 maximum temperature of chains, Parallel Tempering
#' @param nStepExch number of steps between proposing a swap. Optimal
#'        setting is thin=nStepExhc so the swaps alternate with samples.
#' @param typeTemp character, type of temperature ladder (sequential, exp, geom).
#' @param r float, only used if typeTemp="geom".
#' @param maxIterPT maximum number of iterations, Parallel Tempering. If less than burn-in and sampling iterations required, it is recomputed.
#' @param burnInIter1 burn-in iterations in first iteration of the ABMCEM algorithm.
#' @param burnInIter2 burn-in iterations after first iteration of the ABMCEM algorithm.
#' @param out_file_names names for saving tables with acceptance of mutations,
#'        loglikelihoods and score vectors, and acceptance of swaps.
#'        Do NOT add .txt extension.
#' @param changepPerm allows for adaptative permutation mutation probability.
#'
#'
#' @return list of
#' \item{logLik}{log likelihood of the sequences at the last step}
#' \item{beta}{estimation of choice parameters}
#' \item{stdErrorsChoice}{estimation of standard errors, choice model}
#' \item{theta}{estimation of rate parameters}
#' \item{stdErrorsRate}{estimation of standard errors, rate model}
#' \item{index}{number of iterations until convergence of EM algorithm}
#' \item{diff}{value of stopping rule at last iteration}
#'
#' @export
#'
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
MCEMalgorithmRatePT <- function(nmax, net0, net1, theta0, beta0,
                          fixedparameters = c(NA, NA, NA, 20, 20),
                          formula, initTime = 0,
                          endTime = 1, num_cores = 1, errType2 = 0.8,
                          alpha = 0.9, gamma = 0.9, thr = 1e-3,
                          maxIter = 1000, thin = 50,
                          pShort = 0.35, pAug = 0.35, pPerm = 0.3,
                          k_perm = 5 ,nPT = 1,
                          T0 = 1, nStepExch = 10, typeTemp = "sequential",r=1/2,
                          maxIterPT = 1000, burnInIter1 = 8000,
                          burnInIter2 = 1000,
                          out_file_names = c("out_PT_acceptDF",
                                             "out_PT_mcmcDiagDF",
                                             "out_PT_acceptSwitch"),
                          changePPerm = FALSE) {


  actDfnodes <- defineNodes(data.frame(label = colnames(net0)))


  # Initialization of parameters

  # Creation of sequence of events from initial data
  seq <- EMPreprocessing(net0, net1)
  H <- nrow(seq) # Humming distance

  if (any(!beta0 == 0)) {
    beta <- beta0
  } else { # all initial parameters are zero

    parametersAux <- parameters(seq,
                                actDfnodes. = actDfnodes,
                                net0. = net0, formula. = formula
    )

    beta <- data.frame(
      "Crea" = parametersAux$betaCreaDF,
      "Del" = parametersAux$betaDelDF
    )
  }

  if (any(!theta0 == 0)) {
    theta <- theta0
  } else { # all initial parameters are zero

   theta = data.frame("Crea" = c(log(sum(seq$replace==1) / 1 / ncol(net0) )) ,
                      "Del" = c(log(sum(seq$replace==0) / 1 / ncol(net0))))
  }

  seqsPT <- permute(seq, nmax = nPT) # Initial sequences for Parallel Tempering initialization



  cl2 <- makeCluster(num_cores)
  on.exit(stopCluster(cl2))
  clusterEvalQ(cl2, {
    library(goldfish)
    NULL
  })
  #
  # clusterExport(cl, list(
  #   "stepRatePT", "stepRatePTMC", "stepMCMC", "getpDoAugment", "getpDoShort",
  #   "getpDoPerm", "getKelMeMatrix", "getAuxDfE", "stepAugment", "stepShort",
  #   "stepPerm", "getlogLikelihood", "GatherPreprocessingDF",
  #   "sampleVec", "getlogLikelihoodMC","PT_Rate_MCMC","getlogLikelihoodRate",
  #   "getlogLikelihoodRateMC"
  # ))

  diff <- 1000
  index <- 1
  while (diff > thr) {
    # browser()
    k <- 2

    cat("Index: ", index, "\n")
    cat("Diff: ", diff, "\n")
    if (index > 1000) {
      cat("No convergence\n")
      break
    }

    if(index==1) {
      burnInIter = burnInIter1
    }else{
      burnInIter = burnInIter2
    }

    seqsEM <- PT_Rate_MCMC(nmax = nmax, nPT = nPT, seqsPT = seqsPT, H = H,
                           actDfnodes =  actDfnodes, formula = formula,
                           net0 = net0, beta = beta, theta = theta,
                           fixedparameters = fixedparameters,
                           initTime = initTime, endTime = endTime,
                           burnIn=TRUE, burnInIter = burnInIter,
                           maxIter = maxIterPT, thin = thin, T0 = T0,
                           nStepExch = nStepExch, pAug = pAug,
                           pShort = pShort, pPerm = pPerm, k = k_perm,
                           num_cores = num_cores, typeTemp = typeTemp ,r = r,
                           index = index,out_file_names=out_file_names,
                           changePPerm = changePPerm)


    cat("After MCMC_Rate, writing tables should happen now")
    # if(index==1){
    #   write.table(seqsEM$acceptDF, file = paste("out_PT_acceptDF.txt",sep=""), sep = "\t",
    #               row.names = FALSE, col.names = TRUE,append = FALSE)
    #   write.table(seqsEM$mcmcDiagDF, file = paste("out_PT_mcmcDiagDF.txt",sep=""), sep = "\t",
    #               row.names = FALSE, col.names = TRUE,append = FALSE)
    #   write.table(seqsEM$acceptSwitch, file = paste("out_PT_acceptSwitch.txt",sep=""), sep = "\t",
    #               row.names = FALSE, col.names = TRUE,append = FALSE)
    # }else{
    #   write.table(seqsEM$acceptDF, file = paste("out_PT_acceptDF.txt",sep=""), sep = "\t",
    #               row.names = FALSE, col.names = FALSE,append = TRUE)
    #   write.table(seqsEM$mcmcDiagDF, file = paste("out_PT_mcmcDiagDF.txt",sep=""), sep = "\t",
    #               row.names = FALSE, col.names = FALSE,append = TRUE)
    #   write.table(seqsEM$acceptSwitch, file = paste("out_PT_acceptSwitch.txt",sep=""), sep = "\t",
    #               row.names = FALSE, col.names = FALSE,append = TRUE)
    # }

    seqsPT <- unlist(seqsEM$resstepPT, recursive = FALSE)

    logLikPrevChoiceCrea <- sapply(lapply(
      lapply(seqsEM$seqsEM, "[[","newlogLikelihoodStats"), "[[","resCrea"
               ), "[[", "logLikelihood")
    logLikPrevChoiceDel <- sapply(lapply(
      lapply(seqsEM$seqsEM, "[[", "newlogLikelihoodStats"), "[[", "resDel"
               ), "[[", "logLikelihood")
    logLikPrevChoice <- logLikPrevChoiceCrea + logLikPrevChoiceDel # vector of likelihoods lambda^{t-1}

    logLikPrevRateCrea <- sapply(lapply(
      lapply(seqsEM$seqsEM, "[[","newlogLikelihoodStats"), "[[","resCrea"
    ), "[[", "logLikelihood")
    logLikPrevRateDel <- sapply(lapply(
      lapply(seqsEM$seqsEM, "[[", "newlogLikelihoodStats"), "[[", "resDel"
    ), "[[", "logLikelihood")
    logLikPrevRate <- logLikPrevRateCrea + logLikPrevRateDel # vector of likelihoods lambda^{t-1}

    logLikPrev = logLikPrevChoice + logLikPrevRate


    newNRstepChoice <- newtonraphsonStepChoice(
      parameters = beta,
      fixedparameters = fixedparameters,
      seqsEM = seqsEM$seqsEM
    )

    newNRstepRate<- newtonraphsonStepRate(
      parameters = theta,
      seqsEM = seqsEM$seqsEM
    )

    cat("\n Choice:")
    print(newNRstepChoice$parameters)
    cat("\n Rate:")
    print(newNRstepRate$parameters)
    cat("\n")
    splitIndicesPerCore <- splitIndices(length(seqsEM$seqsEM), num_cores)

    logLikCurChoice <- clusterApplyLB(cl2, seq_along(splitIndicesPerCore),
                                getlogLikelihoodMC,seqsEM = seqsEM$seqsEM,
                                beta = newNRstepChoice$parameters,
                                fixedparameters = fixedparameters,
                                splitIndicesPerCore = splitIndicesPerCore,
                                initTime = initTime, endTime = endTime,
                                actDfnodes = actDfnodes, net0 = net0,
                                formula = formula,
                                temp = rep(1, length(seqsEM$seqsEM))
    )


    logLikCurChoice <- unlist(logLikCurChoice, recursive = FALSE)

    logLikCurRate <- getlogLikelihoodRateMC(seqs = seqsEM$seqsEM,
                                            theta = newNRstepRate$parameters,
                                            initTime = initTime,
                                            endTime = endTime,
                                            actDfnodes = actDfnodes,
                                            temp = rep(1,length(seqsEM$seqsEM)))
    logLikCurRate <- unlist(logLikCurChoice, recursive = FALSE)

    logLikCurEM <- sum(sapply(logLikCurChoice, sum)) + sum(logLikCurRate) # total likelihoods for each sequence (Crea+Del+RATES)

    w <- rep(1, length(logLikPrev)) / length(logLikPrev)
    sigmaHat <- sigmaHat(nmax, logLikPrev, logLikCurEM, w)
    ase <- sqrt(sigmaHat / nmax)
    deltaQ <- deltaQ(logLikPrev, logLikCurEM, w)
    lowerBound <- deltaQ - qnorm(alpha) * ase # 80%



    if (lowerBound < 0) {
      cat("Lower bound less than 0\n")
      while (lowerBound < 0) {
        # Estimator is not accepted, new point must be sampled
        # Perform Parallel-Tempering
        seqsEMaux <- PT_Rate_MCMC(nmax = ceiling(nmax / k), nPT = nPT,
                                  seqsPT = seqsPT,
                                  H = H, actDfnodes =  actDfnodes,
                                  formula = formula, net0 = net0, beta = beta,
                                  theta = theta,
                                  fixedparameters = fixedparameters,
                                  initTime = initTime, endTime = endTime,
                                  burnIn=FALSE, burnInIter = 0,
                                  maxIter = maxIterPT, thin = thin, T0 = T0,
                                  nStepExch = nStepExch, pAug = pAug,
                                  pShort = pShort, pPerm = pPerm, k = k_perm,
                                  num_cores = num_cores, typeTemp = typeTemp ,
                                  r = r, index = 2,out_file_names=out_file_names,
                                  changePPerm = changePPerm)
        nmax <- nmax + ceiling(nmax / k)

        # write.table(seqsEMaux$acceptDF, file = paste("out_acceptDF.txt",sep=""), sep = "\t",
        #             row.names = FALSE, col.names = FALSE,append = TRUE)
        # write.table(seqsEMaux$mcmcDiagDF, file = paste("out_mcmcDiagDF.txt",sep=""), sep = "\t",
        #             row.names = FALSE, col.names = FALSE,append = TRUE)
        # write.table(seqsEM$acceptSwitch, file = paste("out_PT_acceptSwitch.txt",sep=""), sep = "\t",
        #             row.names = FALSE, col.names = FALSE,append = TRUE)


        seqsPT <- unlist(seqsEMaux$resstepPT, recursive = FALSE)

        # seqsEM <- c(seqsEM, seqsEMaux) not possible!!
        seqsEM$seqsEM <- c(seqsEM$seqsEM, seqsEMaux$seqsEM)
        # seqsEM$resstepPT <- c(seqsEM$resstepPT, seqsEMaux$resstepPT)
        # seqsEM$acceptSwitch <- rbind(seqsEM$acceptSwitch, seqsEMaux$acceptSwitch)
        # seqsEM$acceptDF <- rbind(seqsEM$acceptDF, seqsEMaux$acceptDF)
        # seqsEM$mcmcDiagDF <- rbind(seqsEM$mcmcDiagDF, seqsEMaux$mcmcDiagDF)
        seqsEM$resstepPT <- seqsPT

        logLikPrevChoiceCrea <- sapply(lapply(
          lapply(seqsEM$seqsEM, "[[","newlogLikelihoodStats"), "[[","resCrea"
        ), "[[", "logLikelihood")
        logLikPrevChoiceDel <- sapply(lapply(
          lapply(seqsEM$seqsEM, "[[", "newlogLikelihoodStats"), "[[", "resDel"
        ), "[[", "logLikelihood")
        logLikPrevChoice <- logLikPrevChoiceCrea + logLikPrevChoiceDel # vector of likelihoods lambda^{t-1}

        logLikPrevRateCrea <- sapply(lapply(
          lapply(seqsEM$seqsEM, "[[","newlogLikelihoodStats"), "[[","resCrea"
        ), "[[", "logLikelihood")
        logLikPrevRateDel <- sapply(lapply(
          lapply(seqsEM$seqsEM, "[[", "newlogLikelihoodStats"), "[[", "resDel"
        ), "[[", "logLikelihood")
        logLikPrevRate <- logLikPrevRateCrea + logLikPrevRateDel # vector of likelihoods lambda^{t-1}

        logLikPrev = logLikPrevChoice + logLikPrevRate


        newNRstepChoice <- newtonraphsonStepChoice(
          parameters = beta,
          fixedparameters = fixedparameters,
          seqsEM = seqsEM$seqsEM
        )

        newNRstepRate<- newtonraphsonStepRate(
          parameters = theta,
          seqsEM = seqsEM$seqsEM
        )
        cat("\n Choice:")
        print(newNRstepChoice$parameters)
        cat("\n Rate:")
        print(newNRstepRate$parameters)
        cat("\n")


        splitIndicesPerCore <- splitIndices(length(seqsEM$seqsEM), num_cores)

        logLikCurChoice <- clusterApplyLB(cl2, seq_along(splitIndicesPerCore),
                                          getlogLikelihoodMC,seqsEM = seqsEM$seqsEM,
                                          beta = newNRstepChoice$parameters,
                                          fixedparameters = fixedparameters,
                                          splitIndicesPerCore = splitIndicesPerCore,
                                          initTime = initTime, endTime = endTime,
                                          actDfnodes = actDfnodes, net0 = net0,
                                          formula = formula,
                                          temp = rep(1, length(seqsEM$seqsEM))
        )


        logLikCurChoice <- unlist(logLikCurChoice, recursive = FALSE)

        logLikCurRate <- getlogLikelihoodRateMC(seqs = seqsEM$seqsEM,
                                                theta = newNRstepRate$parameters,
                                                initTime = initTime,
                                                endTime = endTime,
                                                actDfnodes = actDfnodes,
                                                temp = rep(1,length(seqsEM$seqsEM)))
        logLikCurRate <- unlist(logLikCurChoice, recursive = FALSE)

        logLikCurEM <- sum(sapply(logLikCurChoice, sum)) + sum(logLikCurRate) # total likelihoods for each sequence (Crea+Del+RATES)

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


      beta <- newNRstepChoice$parameters
      # seChoice <- newNRstepChoice$stdErrors

      theta <- newNRstepRate$parameters
      # seRate <- newNRstepRate$stdErrors


      # Update on the number of permutations
      m_start <- sigmaHat^2 * (qnorm(alpha) + qnorm(errType2))^2 / deltaQ^2

      if (m_start > nmax) {
        nmax <- ceiling(m_start)
      }

    }
  }

  stdErrorsChoice <- data.frame(
    "Crea" = rep(0, length(beta$Crea)),
    "Del" = rep(0, length(beta$Crea))
  )
  stdErrorsChoice$Crea <-
    sqrt(diag(newNRstepChoice$inverseInformationUnfixedCrea))
  stdErrorsChoice$Del <-
    sqrt(diag(newNRstepChoice$inverseInformationUnfixedDel))

  stdErrorsRate <- data.frame(
    "Crea" = rep(0, length(theta$Crea)),
    "Del" = rep(0, length(theta$Crea))
  )
  stdErrorsRate$Crea <-
    sqrt(diag(newNRstepRate$inverseInformationUnfixedCrea))
  stdErrorsRate$Del <-
    sqrt(diag(newNRstepRate$inverseInformationUnfixedDel))


  # acceptSwitch <- table(seqsEM$acceptSwitch)
  # acceptDF <- table(seqsEM$acceptDF)

  # mcmcObj <- tapply(seqsEM$mcmcDiagDF[, -ncol(seqsEM$mcmcDiagDF)],
  #                   seqsEM$mcmcDiagDF$temp, mcmc,
  #                   thin = thin
  # )

  # geweke <- lapply(mcmcObj, geweke.diag)
  # ess <- lapply(mcmcObj, effectiveSize)

  return(list(
    "logLik" = logLikCur, "beta" = beta, "stdErrorsChoice" = stdErrorsChoice,
    "theta" = theta, "stdErrorsRate" = stdErrorsRate, "index" = index,
    "diff" = diff
  ))
}





#' Ascent-Based Markov-Chain Expectation-Maximization algorithm with MCMC and rate
#'
#' @description This algorithm allows for an implementation of the ABMCEM algorithm with constant rates.
#'              For the MCMC steps, three types of mutations are done: shortening, augmenting and permutation of the sequence.
#'              \eqn{pShort + pAug + pPerm = 1}
#'
#' @param nmax initial value of number of permuations for step 0
#' @param net0 initial network
#' @param net1 final network
#' @param theta0 initial values for rate parameters
#' @param beta0 initial values for effect parameters
#' @param formula formula of the choice model
#' @param initTime initial time.
#' @param endTime end time.
#' @param errType2 value of power wanted. Default 0.8
#' @param alpha value for lower bound computation. Default 0.9
#' @param gamma value for stopping rule. Default 0.9
#' @param thr threshold for algorithm stop
#' @param maxIter maximum number of iterations for Parallel Tempering
#' @param thin number of mutations steps between samples
#' @param pShort probability of shortening the sequence, Parallel Tempering
#' @param pAug probability of augmenting the sequence, Parallel Tempering
#' @param pPerm probability of permuting two elements of the sequence, Parallel Tempering
#' @param k_perm number of permutations after an augmentation, shortening or permutation mutation.
#' @param burnInIter1 burn-in iterations in first iteration of the ABMCEM algorithm.
#' @param burnInIter2 burn-in iterations after first iteration of the ABMCEM algorithm.
#' @param out_file_names names for saving tables with acceptance of mutations, and
#'        loglikelihoods and score vectors. Do NOT add .txt extension.
#' @param changepPerm allows for adaptative permutation mutation probability.
#'
#'
#' @return list of
#' \item{logLik}{log likelihood of the sequences at the last step}
#' \item{beta}{estimation of choice parameters}
#' \item{stdErrorsChoice}{estimation of standard errors, choice model}
#' \item{theta}{estimation of rate parameters}
#' \item{stdErrorsRate}{estimation of standard errors, rate model}
#' \item{index}{number of iterations until convergence of EM algorithm}
#' \item{diff}{value of stopping rule at last iteration}
#'
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
MCEMalgorithmRateMCMC <- function(nmax, net0, net1, theta0, beta0,
                              fixedparameters = data.frame("Crea"=c(NA, NA, NA,NA,-20,-20),
                                                           "Del"=c(NA,NA,NA,NA,20,20)),
                              formula, initTime = 0,
                              endTime = 1, num_cores = 1, errType2 = 0.8,
                              alpha = 0.9, gamma = 0.9, thr = 1e-3,
                              maxIter = 1000, thin = 10,
                              pShort = 0.35, pAug = 0.35, pPerm = 0.3,
                              k_perm=5, burnInIter1 = 10000,
                              burnInIter2 = 1000,
                              out_file_names = c("out_PT_acceptDF",
                                                 "out_PT_mcmcDiagDF"),
                              changePPerm = FALSE) {


  actDfnodes <- defineNodes(data.frame(label = colnames(net0)))
# browser()

  # Initialization of parameters


  # Creation of sequence of events from initial data
  seq <- EMPreprocessing(net0, net1)
  seqInit = seq
  H <- nrow(seq) # Humming distance

  if (any(!beta0 == 0)) {
    beta <- beta0
  } else { # all initial parameters are zero

    parametersAux <- parameters(seq,
                                actDfnodes. = actDfnodes,
                                net0. = net0, formula. = formula
    )

    beta <- data.frame(
      "Crea" = parametersAux$betaCreaDF,
      "Del" = parametersAux$betaDelDF
    )
  }

  if (any(!theta0 == 0)) {
    theta <- theta0
  } else { # all initial parameters are zero

    theta = data.frame("Crea" = c(log(sum(seq$replace==1) / 1 / ncol(net0) )) ,
                       "Del" = c(log(sum(seq$replace==0) / 1 / ncol(net0))))
  }

  cl2 <- makeCluster(num_cores)
  on.exit(stopCluster(cl2))
  clusterEvalQ(cl2, {
    library(goldfish)
    NULL
  })

  clusterExport(cl2, list(
    "getlogLikelihood", "GatherPreprocessingDF",
    "sampleVec", "getlogLikelihoodMC","getlogLikelihoodRate",
    "getlogLikelihoodRateMC"
  ))

  diff <- 1000
  index <- 1
  while (diff > thr) {
    k <- 2

    cat("Index: ", index, "\n")
    cat("Diff: ", diff, "\n")
    if (index > 1000) {
      cat("No convergence\n")
      break
    }
    # browser()

    if(index==1) {
      burnInIter = burnInIter1
    }else{
      burnInIter = burnInIter2
    }

    # browser()
    seqsEM <- MCMC_rate(nmax = nmax, seqInit = seqInit, H = H,
                        actDfnodes = actDfnodes, formula = formula,
                        net0 = net0, beta = beta, theta = theta,
                        fixedparameters = fixedparameters,
                        initTime = initTime, endTime = endTime,
                        burnIn = TRUE, burnInIter = burnInIter,
                        maxIter = 100000, thin = thin,
                        pAug = pAug, pShort = pShort,pPerm = pPerm, k = k_perm,
                        index = index,out_file_names=out_file_names,
                        changePPerm = changePPerm
    )


      cat("After MCMC_Rate, writing tables should happen now")

    seqInit <- seqsEM$resstepPT

    logLikPrevChoiceCrea <- sapply(lapply(
      lapply(seqsEM$seqsEM, "[[","newlogLikelihoodStats"), "[[","resCrea"
    ), "[[", "logLikelihood")
    logLikPrevChoiceDel <- sapply(lapply(
      lapply(seqsEM$seqsEM, "[[", "newlogLikelihoodStats"), "[[", "resDel"
    ), "[[", "logLikelihood")

    # if(!is.numeric(logLikPrevChoiceCrea + logLikPrevChoiceDel ) |
    #    is.na(logLikPrevChoiceCrea + logLikPrevChoiceDel )){
    #   save.image(file='myEnvironment.RData')
    # }

    logLikPrevChoice <- logLikPrevChoiceCrea + logLikPrevChoiceDel # vector of likelihoods lambda^{t-1}

    logLikPrevRateCrea <- sapply(lapply(
      lapply(seqsEM$seqsEM, "[[","newloglikRate"), "[[","resCrea"
    ), "[[", "logLikelihood")
    logLikPrevRateDel <- sapply(lapply(
      lapply(seqsEM$seqsEM, "[[", "newloglikRate"), "[[", "resDel"
    ), "[[", "logLikelihood")
    logLikPrevRate <- logLikPrevRateCrea + logLikPrevRateDel # vector of likelihoods lambda^{t-1}

    logLikPrev = logLikPrevChoice + logLikPrevRate


    newNRstepChoice <- newtonraphsonStepChoice(
      parameters = beta,
      fixedparameters = fixedparameters,
      seqsEM = seqsEM$seqsEM
    )

    newNRstepRate<- newtonraphsonStepRate(
      parameters = theta,
      seqsEM = seqsEM$seqsEM
    )

    cat("\n Choice:")
    print(newNRstepChoice$parameters)
    cat("\n Rate:")
    print(newNRstepRate$parameters)
    cat("\n")
    splitIndicesPerCore <- splitIndices(length(seqsEM$seqsEM), num_cores)

    logLikCurChoice <- clusterApplyLB(cl2, seq_along(splitIndicesPerCore),
                                      getlogLikelihoodMC,seqsEM = seqsEM$seqsEM,
                                      beta = newNRstepChoice$parameters,
                                      fixedparameters = fixedparameters,
                                      splitIndicesPerCore = splitIndicesPerCore,
                                      initTime = initTime, endTime = endTime,
                                      actDfnodes = actDfnodes, net0 = net0,
                                      formula = formula,
                                      temp = rep(1, length(seqsEM$seqsEM))
    )


    logLikCurChoice <- unlist(logLikCurChoice, recursive = FALSE)

    logLikCurRate <- getlogLikelihoodRateMC(seqs = seqsEM$seqsEM,
                                            theta = newNRstepRate$parameters,
                                            initTime = initTime,
                                            endTime = endTime,
                                            actDfnodes = actDfnodes,
                                            temp = rep(1,length(seqsEM$seqsEM)))
    logLikCurRate <- unlist(logLikCurRate, recursive = FALSE)

    logLikCurEM <- sum(sapply(logLikCurChoice, sum)) + sum(logLikCurRate) # total likelihoods for each sequence (Crea+Del+RATES)

    w <- rep(1, length(logLikPrev)) / length(logLikPrev)
    sigmaHat <- sigmaHat(nmax, logLikPrev, logLikCurEM, w)
    ase <- sqrt(sigmaHat / nmax)
    deltaQ <- deltaQ(logLikPrev, logLikCurEM, w)
    lowerBound <- deltaQ - qnorm(alpha) * ase # 80%



    if (lowerBound < 0) {
      cat("Lower bound less than 0\n")
      while (lowerBound < 0) {
        cat("\t k=",k)

        # Estimator is not accepted, new point must be sampled
        # Perform Parallel-Tempering

        seqsEMaux <-MCMC_rate(nmax = ceiling(nmax/k), seqInit = seqInit, H = H,
                  actDfnodes = actDfnodes, formula = formula,
                  net0 = net0, beta = beta, theta = theta,
                  fixedparameters = fixedparameters,
                  initTime = initTime, endTime = endTime,
                  burnIn = FALSE, burnInIter = 0,
                  maxIter = 100000, thin = thin,
                  pAug = pAug, pShort = pShort,pPerm = pPerm, k = k_perm,
                  index = 2,out_file_names=out_file_names,
                  changePPerm = changePPerm
        )

        nmax <- nmax + ceiling(nmax / k)

        # write.table(seqsEMaux$acceptDF, file = paste("out_acceptDF.txt",sep=""), sep = "\t",
        #             row.names = FALSE, col.names = FALSE,append = TRUE)
        # write.table(seqsEMaux$mcmcDiagDF, file = paste("out_mcmcDiagDF.txt",sep=""), sep = "\t",
        #             row.names = FALSE, col.names = FALSE,append = TRUE)

        seqInit <- seqsEMaux$resstepPT

        # seqsEM <- c(seqsEM, seqsEMaux) not possible!!
        seqsEM$seqsEM <- c(seqsEM$seqsEM, seqsEMaux$seqsEM)
        seqsEM$resstepPT <- seqsEMaux$resstepPT
        # seqsEM$acceptSwitch <- rbind(seqsEM$acceptSwitch, seqsEMaux$acceptSwitch)
        # seqsEM$acceptDF <- rbind(seqsEM$acceptDF, seqsEMaux$acceptDF)
        # seqsEM$mcmcDiagDF <- rbind(seqsEM$mcmcDiagDF, seqsEMaux$mcmcDiagDF)


        logLikPrevChoiceCrea <- sapply(lapply(
          lapply(seqsEM$seqsEM, "[[","newlogLikelihoodStats"), "[[","resCrea"
        ), "[[", "logLikelihood")
        logLikPrevChoiceDel <- sapply(lapply(
          lapply(seqsEM$seqsEM, "[[", "newlogLikelihoodStats"), "[[", "resDel"
        ), "[[", "logLikelihood")
        logLikPrevChoice <- logLikPrevChoiceCrea + logLikPrevChoiceDel # vector of likelihoods lambda^{t-1}

        logLikPrevRateCrea <- sapply(lapply(
          lapply(seqsEM$seqsEM, "[[","newlogLikelihoodStats"), "[[","resCrea"
        ), "[[", "logLikelihood")
        logLikPrevRateDel <- sapply(lapply(
          lapply(seqsEM$seqsEM, "[[", "newlogLikelihoodStats"), "[[", "resDel"
        ), "[[", "logLikelihood")
        logLikPrevRate <- logLikPrevRateCrea + logLikPrevRateDel # vector of likelihoods lambda^{t-1}

        logLikPrev = logLikPrevChoice + logLikPrevRate


        newNRstepChoice <- newtonraphsonStepChoice(
          parameters = beta,
          fixedparameters = fixedparameters,
          seqsEM = seqsEM$seqsEM
        )

        newNRstepRate<- newtonraphsonStepRate(
          parameters = theta,
          seqsEM = seqsEM$seqsEM
        )

        cat("\n Choice:")
        print(newNRstepChoice$parameters)
        cat("\n Rate:")
        print(newNRstepRate$parameters)
        cat("\n")

        splitIndicesPerCore <- splitIndices(length(seqsEM$seqsEM), num_cores)

        logLikCurChoice <- clusterApplyLB(cl2, seq_along(splitIndicesPerCore),
                                          getlogLikelihoodMC,seqsEM = seqsEM$seqsEM,
                                          beta = newNRstepChoice$parameters,
                                          fixedparameters = fixedparameters,
                                          splitIndicesPerCore = splitIndicesPerCore,
                                          initTime = initTime, endTime = endTime,
                                          actDfnodes = actDfnodes, net0 = net0,
                                          formula = formula,
                                          temp = rep(1, length(seqsEM$seqsEM))
        )


        logLikCurChoice <- unlist(logLikCurChoice, recursive = FALSE)

        logLikCurRate <- getlogLikelihoodRateMC(seqs = seqsEM$seqsEM,
                                                theta = newNRstepRate$parameters,
                                                initTime = initTime,
                                                endTime = endTime,
                                                actDfnodes = actDfnodes,
                                                temp = rep(1,length(seqsEM$seqsEM)))
        logLikCurRate <- unlist(logLikCurRate, recursive = FALSE)

        logLikCurEM <- sum(sapply(logLikCurChoice, sum)) + sum(logLikCurRate) # total likelihoods for each sequence (Crea+Del+RATES)

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


      beta <- newNRstepChoice$parameters
      theta <- newNRstepRate$parameters

      # Update on the number of permutations
      m_start <- sigmaHat^2 * (qnorm(alpha) + qnorm(errType2))^2 / deltaQ^2

      if (m_start > nmax) {
        nmax <- ceiling(m_start)
      }

    }
  }

  stdErrorsChoice <- data.frame(
    "Crea" = rep(0, length(beta$Crea)),
    "Del" = rep(0, length(beta$Crea))
  )
  stdErrorsChoice$Crea <-
    sqrt(diag(newNRstepChoice$inverseInformationUnfixedCrea))
  stdErrorsChoice$Del <-
    sqrt(diag(newNRstepChoice$inverseInformationUnfixedDel))

  stdErrorsRate <- data.frame(
    "Crea" = rep(0, length(theta$Crea)),
    "Del" = rep(0, length(theta$Crea))
  )
  stdErrorsRate$Crea <-
    sqrt(diag(newNRstepRate$inverseInformationUnfixedCrea))
  stdErrorsRate$Del <-
    sqrt(diag(newNRstepRate$inverseInformationUnfixedDel))

  # acceptDF <- table(seqsEM$acceptDF)
  #
  # mcmcObj <- tapply(seqsEM$mcmcDiagDF[, -ncol(seqsEM$mcmcDiagDF)],
  #                   seqsEM$mcmcDiagDF$temp, mcmc,
  #                   thin = thin
  # )

  # geweke <- lapply(mcmcObj, geweke.diag)
  # ess <- lapply(mcmcObj, effectiveSize)

  return(list(
    "logLik" = logLikCur, "beta" = beta, "stdErrorsChoice" = stdErrorsChoice,
    "theta" = theta, "stdErrorsRate" = stdErrorsRate, "index" = index,
    "diff" = diff
  ))
}





