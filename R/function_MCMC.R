# MCMC auxiliar functions --------------------


#' Get probability of doing a augmenting step.
#'
getpDoAugment <- function(gammaEplus, gammaPlus, p, m, me, sender, receiver, typeA, indexSR) {
  if (typeA == "same") {
    if (length(indexSR) < 1) {
      pDoStep <- (gammaEplus[sender, receiver] / gammaPlus) * p * (1 / (m + 1))
    } else {
      pDoStep <- (gammaEplus[sender, receiver] / gammaPlus) * p * (1 / (m - me[sender, receiver] + 1))
    }
  } else {
    pDoStep <- (gammaEplus[sender, receiver] / gammaPlus) * p * (1 / (m - me[sender, receiver] + 1)) * (1 / (m - me[sender, receiver]))
  }

  return(pDoStep)
}

#' Get probability of doing a shortening step.
#'
getpDoShort <- function(gammaEminus, gammaMinus, p, uniqueR, sender, receiver, typeS) {
  if (typeS == "same") {
    pDoStep <- (gammaEminus[sender, receiver] / gammaMinus) * p * (1 / uniqueR)
  } else {
    pDoStep <- (gammaEminus[sender, receiver] / gammaMinus) * (1 - p) * (1 / uniqueR) * (1 / (uniqueR - 1))
  }

  return(pDoStep)
}

#' Get probability of doing a permutation step.
#'
getpDoPerm <- function(probVec, probVec2, e1, e2, sender1, receiver1, sender2, receiver2, seq) {
  pDoStep <- probVec[e1] * probVec2[e2] * (1 / nrow(seq[seq$sender == sender1 & seq$receiver == receiver1, ])) * (1 / nrow(seq[seq$sender == sender2 & seq$receiver == receiver2, ]))
  return(pDoStep)
}


#' Get auxiliary matrix \eqn{K_{el}} and \eqn{m_e} for MCMC steps.
#'
#' @param seq data frame with the sequence of events up to time \eqn{t}.
#' @param actDfnodes vector with the labels of the actors of the network.
#'
#' @return list with matrix \eqn{K_{el}}, \eqn{K_{e,l>1}}, \eqn{m_e} and auxiliary data frame
#' @export
#'
#' @examples seq <- data.frame("sender" = c(1, 2, 3), "receiver" = c(3, 2, 1), "replace" = c(0, 1, 1))
#' getKelMeMatrix(seq, actDfnodes = c(1, 2, 3))
getKelMeMatrix <- function(seq, actDfnodes) {
  seq <- cbind(seq, "row" = as.integer(rownames(seq)))


  auxDf <- lapply(actDfnodes, function(x) lapply(actDfnodes, function(i) subset(seq, subset = (sender == x & receiver == i))))

  for (i in 1:length(actDfnodes)) {
    for (j in 1:length(actDfnodes)) {
      if (nrow(auxDf[[i]][[j]]) > 0) {
        auxDf[[i]][[j]]$row <- as.integer(auxDf[[i]][[j]]$row)
        auxDf[[i]][[j]]$rowDiff <- c(0, diff(auxDf[[i]][[j]]$row))
      }
    }
  }


  Kel_g1 <- matrix(0, nrow = length(actDfnodes), ncol = length(actDfnodes)) # Notation from Johan´s paper, sum for l greater than 1
  Kel_ge1 <- matrix(0, nrow = length(actDfnodes), ncol = length(actDfnodes)) # Notation from Johan´s paper, sum for l greater or equal than 1

  for (i in 1:nrow(Kel_g1)) {
    for (j in 1:ncol(Kel_g1)) {
      if (nrow(auxDf[[i]][[j]]) > 0) {
        Kel_ge1[i, j] <- nrow(auxDf[[i]][[j]][auxDf[[i]][[j]]$rowDiff != 1, ])
        v_rowDiff <- auxDf[[i]][[j]]$rowDiff
        not_1 <- which(v_rowDiff != 1)
        shifted_vector <- c(v_rowDiff[-1], 0)
        Kel_g1[i, j] <- sum(v_rowDiff[-length(v_rowDiff)] != 1 & shifted_vector[-1] == 1)
      }
    }
  }

  me <- t(sapply(auxDf, function(x) sapply(x, nrow)))

  return(list(Kel_g1 = Kel_g1, Kel_ge1 = Kel_ge1, me = me, auxDf = auxDf))
}




#' Get auxiliary data frame for tie \eqn{e}
#'
#' @param auxDf data frame from getKelMeMatrix.
#' @param sender sender node
#' @param receiver receiver node
#'
#' @return list of auxDfE (data frame).
#'
#' @export
#'
getAuxDfE <- function(auxDf, sender, receiver) {
  auxDfE <- auxDf[[sender]][[receiver]]
  indexNo1 <- which(auxDfE$rowDiff != 1)
  auxDfE$run <- rep(0, nrow(auxDfE))
  run <- 1
  if (length(indexNo1) > 0) {
    for (i in 1:(length(indexNo1) - 1)) {
      auxDfE$run[indexNo1[i]:indexNo1[i + 1]] <- run
      run <- run + 1
    }
    auxDfE$run[indexNo1[length(indexNo1)]:nrow(auxDfE)] <- run
  }
  return(auxDfE)
}

#' Step for augmenting a given sequence by adding a tie \eqn{e} two times.
#'
#' @param seq data frame, sequence of events.
#' @param tieNames vector, labels of ties (e.g., "12" if the tie is from sender 1 and receiver 2).
#' @param gammaEplus matrix, probabilities needed for the computation.
#' @param gammaPlus float, sum of gammaEplus.
#' @param m float, sum of me.
#' @param me matrix, probabilities needed for the computation.
#' @param net0 matrix, initial network.
#'
#' @return list of
#' \item{sender}{sender of the new event.}
#' \item{receiver}{receiver of the new event.}
#' \item{place}{positions were the event was included.}
#' \item{typeA}{type of augmentation (two events together or separetly).}
#' \item{pDoStep}{probability of doing the step.}
#'
#' @export
#'
stepAugment <- function(seq, tieNames, gammaEplus, gammaPlus, m, me, net0,pAug) {
  # Choose element to be inserted:
  # browser()
  e <- sampleVec(tieNames, size = 1, prob = as.vector(t(gammaEplus / gammaPlus)))

  # Choose to augment in the same e-run or different ones
  sender <- as.integer(strsplit(e, "-")[[1]][1])
  receiver <- as.integer(strsplit(e, "-")[[1]][2])

  p <- (m - me[sender, receiver] + 1) / gammaEplus[sender, receiver]
  typeA <- sampleVec(c("same", "diff"), size = 1, prob = c(p, 1 - p))
  indexSR <- which(seq$sender == sender & seq$receiver == receiver)



  if (typeA == "same") {
    if (length(indexSR) < 1) {
      place <- sample(m + 1, size = 1)
    } else {
      place <- sampleVec(seq(m + 1)[-indexSR], size = 1)
    }

    if (place == 1) {
      newseq <- rbind(
        c(sender, receiver, 4, place),
        c(sender, receiver, 4, place + 1),
        seq
      )
    } else {
      newseq <- rbind(
        seq[1:(place - 1), ],
        c(sender, receiver, 4, place),
        c(sender, receiver, 4, place + 1),
        seq[(place):m, ]
      )
    }

    nEdges <- which(newseq$sender == sender & newseq$receiver == receiver)
    initAction <- as.integer(seq[which(seq$sender == sender & seq$receiver == receiver), "replace"][1])
    if (is.na(initAction)) {
      initAction <- net0[sender, receiver] + 1
    }
    newseq[nEdges, "replace"] <- seq(initAction, (initAction + length(nEdges) - 1)) %% 2
    newseq$row <- seq(1, nrow(newseq))

    pDoStep <- getpDoAugment(gammaEplus, gammaPlus, p, m, me, sender, receiver, typeA, indexSR)*pAug
  } else if (typeA == "diff") {
    # Choose thow places separated at least from one tie different than e
    if (length(indexSR) < 1) {
      place <- sample(m + 1, size = 2)
    } else {
      place <- sampleVec(seq(m + 1)[-indexSR], size = 2)
    }

    place <- sort(place)


    if (place[1] == 1) {
      if (place[2] < (m + 1)) {
        newseq <- rbind(
          c(sender, receiver, 4, place[1]),
          seq[(place[1]):(place[2] - 1), ],
          c(sender, receiver, 4, place[2]),
          seq[(place[2]):m, ]
        )
      } else {
        newseq <- rbind(
          c(sender, receiver, 4, place[1]),
          seq[(place[1]):(place[2] - 1), ],
          c(sender, receiver, 4, place[2])
        )
      }
    } else {
      if (place[2] < (m + 1)) {
        newseq <- rbind(
          seq[1:(place[1] - 1), ],
          c(sender, receiver, 4, place[1]),
          seq[(place[1]):(place[2] - 1), ],
          c(sender, receiver, 4, place[2]),
          seq[(place[2]):m, ]
        )
      } else {
        newseq <- rbind(
          seq[1:(place[1] - 1), ],
          c(sender, receiver, 4, place[1]),
          seq[(place[1]):(place[2] - 1), ],
          c(sender, receiver, 4, place[2])
        )
      }
    }

    nEdges <- which(newseq$sender == sender & newseq$receiver == receiver)
    initAction <- as.integer(seq[which(seq$sender == sender & seq$receiver == receiver), "replace"][1])
    if (is.na(initAction)) {
      initAction <- net0[sender, receiver] + 1
    }
    newseq[nEdges, "replace"] <- seq(initAction, (initAction + length(nEdges) - 1)) %% 2
    newseq$row <- seq(1, nrow(newseq))

    pDoStep <- getpDoAugment(gammaEplus, gammaPlus, p, m, me, sender, receiver, typeA, indexSR)*pAug
  }

  return(list(
    newseq = newseq, sender = sender, receiver = receiver, place = place, typeA = typeA,
    pDoStep = pDoStep
  ))
}


#' Step for shortening a given sequence by removing a tie \eqn{e} in two positions.
#'
#' @param seq data frame, sequence of events.
#' @param tieNames vector, labels of ties (e.g., "12" if the tie is from sender 1 and receiver 2).
#' @param gammaEminus matrix, probabilities needed for the computation.
#' @param gammaMinus float, sum of gammaEplus.
#' @param m float, sum of me.
#' @param me matrix, probabilities needed for the computation.
#' @param auxDf data frame, from getKelMeMatrix
#'
#' @return list of
#' \item{sender}{sender of the new event.}
#' \item{receiver}{receiver of the new event.}
#' \item{place}{positions were the event was included.}
#' \item{typeS}{type of shortening (two events together or separetly).}
#' \item{pDoStep}{probability of doing the step.}
#'
#' @export
#'
stepShort <- function(seq, tieNames, gammaEminus, gammaMinus, m, me, Kel_g1, auxDf,pShort) {
  # Choose element to be deleted:
  e <- sampleVec(tieNames, size = 1, prob = as.vector(t(gammaEminus / gammaMinus)))

  # Choose to remove from the same e-run or from two distinct e-runs
  sender <- as.integer(strsplit(e, "-")[[1]][1])
  receiver <- as.integer(strsplit(e, "-")[[1]][2])

  p <- Kel_g1[sender, receiver] / gammaEminus[sender, receiver]
  typeS <- sampleVec(c("same", "diff"), size = 1, prob = c(p, 1 - p))

  # Get auxiliar variables with indexes for the runs
  auxDfE <- getAuxDfE(auxDf, sender, receiver)


  if (typeS == "same") {
    # Choose the e-runs of more than 1 element, select one and delete the first 2 elements.
    sampleRun <- sampleVec(which(table(auxDfE$run) > 1), size = 1)
    indexSampleRun <- auxDfE$row[auxDfE$run == sampleRun][1:2]

    newseq <- seq[-(indexSampleRun), ]
    newseq$row <- 1:nrow(newseq)

    place <- list(sampleRun = sampleRun, indexSampleRun = indexSampleRun)

    pDoStep <- getpDoShort(gammaEminus, gammaMinus, p, length(unique(auxDfE$run)), sender, receiver, typeS)*pShort
  } else if (typeS == "diff") {
    # Choose two different e-runs, and delete the first element of each one of them
    sampleRun <- sampleVec(unique(auxDfE$run), size = 2, replace = FALSE)
    indexSampleRun <- c(auxDfE$row[auxDfE$run == sampleRun[1]][1], auxDfE$row[auxDfE$run == sampleRun[2]][1])

    newseq <- seq[-(indexSampleRun), ]
    newseq$row <- 1:nrow(newseq)

    # now we have to fix the replaceOG
    nEdges <- which(newseq$sender == sender & newseq$receiver == receiver)
    initAction <- as.integer(seq[which(seq$sender == sender & seq$receiver == receiver), "replace"][1])
    newseq[nEdges, "replace"] <- seq(initAction, (initAction + length(nEdges) - 1)) %% 2

    place <- list(sampleRun = sampleRun, indexSampleRun = indexSampleRun)
    pDoStep <- getpDoShort(gammaEminus, gammaMinus, p, length(unique(auxDfE$run)), sender, receiver, typeS)*pShort
  }

  return(list(
    newseq = newseq, sender = sender, receiver = receiver, place = place, typeS = typeS,
    pDoStep = pDoStep, auxDfE = auxDfE
  ))
}



#' Step for permuting two elements in a given sequence.
#'
#' @param seq data frame, sequence of events.
#' @param tieNames vector, labels of ties (e.g., "12" if the tie is from sender 1 and receiver 2).
#' @param m float, sum of me.
#' @param me matrix, probabilities needed for the computation.
#'
#' @return list of
#' \item{sender}{sender of the new event.}
#' \item{receiver}{receiver of the new event.}
#' \item{place}{positions were the event was included.}
#' \item{pDoStep}{probability of doing the step.}
#'
#' @export
#'
stepPerm <- function(seq, tieNames, m, me) {
  # Choose elements to be permuted:
  probVec <- as.vector(t(me / m))
  probVec2 <- as.vector(t(me / (m - me)))
  names(probVec) <- tieNames
  names(probVec2) <- tieNames

  e1 <- sampleVec(tieNames, size = 1, prob = probVec)
  e2 <- sampleVec(tieNames[tieNames != e1], size = 1, prob = probVec2[names(probVec) != e1])

  # After getting the edge identifier, select at random e1 and e2 from all the runs
  sender1 <- as.integer(strsplit(e1, "-")[[1]][1])
  receiver1 <- as.integer(strsplit(e1, "-")[[1]][2])
  place1 <- as.integer(sampleVec(as.numeric(seq[seq$sender == sender1 & seq$receiver == receiver1, "row"]), size = 1))
  sender2 <- as.integer(strsplit(e2, "-")[[1]][1])
  receiver2 <- as.integer(strsplit(e2, "-")[[1]][2])
  place2 <- as.integer(sampleVec(seq[seq$sender == sender2 & seq$receiver == receiver2, ]$row, size = 1))
  auxRow <- seq[place1, ]
  newseq <- seq
  newseq[place1, ] <- newseq[place2, ]
  newseq[place2, ] <- auxRow
  rm(auxRow)

  nEdges <- which(newseq$sender == sender1 & newseq$receiver == receiver1)
  initAction <- as.integer(seq[which(seq$sender == sender1 & seq$receiver == receiver1), "replace"][1])
  newseq[nEdges, "replace"] <- seq(initAction, (initAction + length(nEdges) - 1)) %% 2
  nEdges <- which(newseq$sender == sender2 & newseq$receiver == receiver2)
  initAction <- as.integer(seq[which(seq$sender == sender2 & seq$receiver == receiver2), "replace"][1])
  newseq[nEdges, "replace"] <- seq(initAction, (initAction + length(nEdges) - 1)) %% 2

  newseq$row <- 1:nrow(newseq)

  pDoStep <- probVec[e1] * probVec2[e2] * (1 / nrow(seq[seq$sender == sender1 & seq$receiver == receiver1, ])) * (1 / nrow(seq[seq$sender == sender2 & seq$receiver == receiver2, ]))

  return(list(
    newseq = newseq, sender = c(sender1, sender2), receiver = c(receiver1, receiver2), place = c(place1, place2),
    pDoStep = pDoStep
  ))
}



# MCMC -------------------------------------

#' MCMC step
#'
#' @description Performs an step of MCMC given a type:
#' \item Type 1: augmentation of sequence.
#' \item Type 2: shortening of sequence.
#' \item Type 3: permutation of 2 elements of the sequence.
#'
#' @param seq data frame, sequence of events.
#' @param type integer, type of step.
#' @param actDfnodesLab vector, node labs.
#' @param tieNames vector, tie labs (e.g. "12" is a tie from sender 1 to receiver 2).
#' @param formula formula, the one used in goldfish model.
#' @param net0 matrix, initial observed network.
#' @param beta list, estimator of choice parameters of the model (creation and deletion)
#' @param theta dataframe, estimator of rate parameters of the model (creation and deletion) (constant)
#' @param initTime initial time
#' @param endTime end time
#'
#' @return list of
#' \item{newseq}{new sequence.}
#' \item{loglikSeq}{log-likelihood from original sequence.}
#' \item{newloglikSeq}{log-likelihood from new sequence.}
#' \item{step}{depending on the type of step, result of [stepAugment()],[stepShort()] or [stepPerm()].}
#' \item{pUndoStep}{probability of undoing the step of MH acceptance rate.}
#'
#' @export
#'
stepMCMC <- function(seq, type, actDfnodesLab, tieNames, formula, net0, beta, theta, initTime, endTime) {
  getKelMeMatrix <- getKelMeMatrix(seq, actDfnodesLab)
  Kel_g1 <- getKelMeMatrix$Kel_g1
  Kel_ge1 <- getKelMeMatrix$Kel_ge1
  gammaEminus <- choose(Kel_ge1, 2) + Kel_g1
  gammaMinus <- sum(gammaEminus)
  me <- getKelMeMatrix$me
  gammaEplus <- choose(nrow(seq) - me + 2, 2)
  gammaPlus <- sum(gammaEplus)
  m <- nrow(seq)
  auxDf <- getKelMeMatrix$auxDf



  if (type == 1) { # Augmentation
    step <- stepAugment(seq, tieNames, gammaEplus, gammaPlus, m, me, net0)

    getNewKelMeMatrix <- getKelMeMatrix(step$newseq, actDfnodesLab)
    newAuxDfE <- getAuxDfE(getNewKelMeMatrix$auxDf, step$sender, step$receiver)
    auxDfE <- getAuxDfE(auxDf, step$sender, step$receiver)

    newKel_g1 <- getNewKelMeMatrix$Kel_g1
    newKel_ge1 <- getNewKelMeMatrix$Kel_ge1
    newGammaEminus <- choose(newKel_ge1, 2) + newKel_g1
    newGammaMinus <- sum(newGammaEminus)
    newMe <- getNewKelMeMatrix$me
    newGammaEplus <- choose(nrow(newseq) - newMe + 2, 2)
    newGammaPlus <- sum(newGammaEplus)
    newM <- nrow(newseq)

    newp <- newKel_g1[step$sender, step$receiver] / newGammaEminus[step$sender, step$receiver]
    if ((length(unique(auxDfE$run)) == length(unique(newAuxDfE$run))) | step$typeA == "same") {
      pUndoStep <- getpDoShort(newGammaEminus, newGammaMinus, newp, length(unique(newAuxDfE$run)), step$sender, step$receiver, "same")
    } else {
      pUndoStep <- getpDoShort(newGammaEminus, newGammaMinus, newp, length(unique(newAuxDfE$run)), step$sender, step$receiver, "diff")
    }

    loglikSeq <- logLikelihoodMC(indexCore = 1, list(seq[, -which(colnames(seq) == "row")]), beta, theta, initTime, endTime, splitIndicesPerCore = list(1), actDfnodes = actDfnodes, net0 = net0, formula = formula)
    newloglikSeq <- logLikelihoodMC(indexCore = 1, list(step$newseq[, -which(colnames(step$newseq) == "row")]), beta, theta, initTime, endTime, splitIndicesPerCore = list(1), actDfnodes = actDfnodes, net0 = net0, formula = formula)
  } else if (type == 2) { # Shortening
    step <- stepShort(seq, tieNames, gammaEminus, gammaMinus, m, me, Kel_g1, auxDf)
    getNewKelMeMatrix <- getKelMeMatrix(step$newseq, actDfnodesLab)
    newAuxDfE <- getAuxDfE(getNewKelMeMatrix$auxDf, sender, receiver)
    auxDfE <- getAuxDfE(auxDf, sender, receiver)

    newKel_g1 <- getNewKelMeMatrix$Kel_g1
    newKel_ge1 <- getNewKelMeMatrix$Kel_ge1
    newGammaEminus <- choose(newKel_ge1, 2) + newKel_g1
    newGammaMinus <- sum(newGammaEminus)
    newMe <- getNewKelMeMatrix$me
    newGammaEplus <- choose(nrow(newseq) - newMe + 2, 2)
    newGammaPlus <- sum(newGammaEplus)
    newM <- nrow(newseq)

    newp <- (newM - newMe[step$sender, step$receiver] + 1) / newGammaEplus[step$sender, step$receiver]

    newindexSR <- which(step$newseq$sender == step$sender & step$newseq$receiver == step$receiver)

    if (step$typeS == "same") {
      pUndoStep <- getpDoAugment(newGammaEplus, newGammaPlus, newp, newM, newMe, step$sender, step$receiver, "same", newindexSR)
    } else if (step$typeS == "diff") {
      pUndoStep <- getpDoAugment(newGammaEplus, newGammaPlus, newp, newM, newMe, step$sender, step$receiver, "diff", newindexSR)
    }

    loglikSeq <- logLikelihoodMC(indexCore = 1, list(seq[, -which(colnames(seq) == "row")]), beta, theta, initTime, endTime, splitIndicesPerCore = list(1), actDfnodes = actDfnodes, net0 = net0, formula = formula)
    newloglikSeq <- logLikelihoodMC(indexCore = 1, list(step$newseq[, -which(colnames(step$newseq) == "row")]), beta, theta, initTime, endTime, splitIndicesPerCore = list(1), actDfnodes = actDfnodes, net0 = net0, formula = formula)
  } else if (type == 3) { # Permutation
    step <- stepPerm(seq, tieNames, m, me)
    getNewKelMeMatrix <- getKelMeMatrix(step$newseq, actDfnodesLab)
    newAuxDfE <- getAuxDfE(getNewKelMeMatrix$auxDf, sender, receiver)
    auxDfE <- getAuxDfE(auxDf, sender, receiver)

    newKel_g1 <- getNewKelMeMatrix$Kel_g1
    newKel_ge1 <- getNewKelMeMatrix$Kel_ge1
    newGammaEminus <- choose(newKel_ge1, 2) + newKel_g1
    newGammaMinus <- sum(newGammaEminus)
    newMe <- getNewKelMeMatrix$me
    newGammaEplus <- choose(nrow(newseq) - newMe + 2, 2)
    newGammaPlus <- sum(newGammaEplus)
    newM <- nrow(newseq)


    # probVec <- newMe[step$sender[2], step$receiver[2]] / newM
    # probVec2 <- newMe[step$sender[1], step$receiver[1]] / (newM - newMe[step$sender[2], step$receiver[2]])
    # probVec3 <- newMe[step$sender[1], step$receiver[1]] / newM
    # probVec4 <- newMe[step$sender[2], step$receiver[2]] / (newM - newMe[step$sender[1], step$receiver[1]])
    #
    # pUndoStep <- (probVec * probVec2 + probVec3 * probVec4) *
    #   (1 / nrow(newseq[newseq$sender == step$sender[1] & newseq$receiver == step$receiver[1], ])) *
    #   (1 / nrow(newseq[newseq$sender == step$sender[2] & newseq$receiver == step$receiver[2], ]))
    pUndoStep <- step$pDoStep
    loglikSeq <- logLikelihoodMC(indexCore = 1, list(seq[, -which(colnames(seq) == "row")]), beta, theta, initTime, endTime, splitIndicesPerCore = list(1), actDfnodes = actDfnodes, net0 = net0, formula = formula)
    newloglikSeq <- logLikelihoodMC(indexCore = 1, list(step$newseq[, -which(colnames(step$newseq) == "row")]), beta, theta, initTime, endTime, splitIndicesPerCore = list(1), actDfnodes = actDfnodes, net0 = net0, formula = formula)
  }

  return(list(
    "newseq" = newseq, "loglikSeq" = loglikSeq, "newloglikSeq" = newloglikSeq,
    "step" = step, "pUndoStep" = pUndoStep
  ))
}


#' MCMC step reversible jump proposals
#'
#' @description Performs an step of MCMC given a type:
#' \item Type 1: augmentation of sequence.
#' \item Type 2: shortening of sequence.
#' followed by k permutation steps.
#'
#' @param seq data frame, sequence of events.
#' @param type integer, type of step.
#' @param actDfnodesLab vector, node labs.
#' @param tieNames vector, tie labs (e.g. "12" is a tie from sender 1 to receiver 2).
#' @param formula formula, the one used in goldfish model.
#' @param net0 matrix, initial observed network.
#' @param beta list, estimator of choice parameters of the model (creation and deletion)
#' @param theta dataframe, estimator of rate parameters of the model (creation and deletion) (constant)
#' @param initTime initial time
#' @param endTime end time
#' @param k number of permutations after augmentations/shortening
#'
#' @return list of
#' \item{newseq}{new sequence.}
#' \item{loglikSeq}{log-likelihood from original sequence.}
#' \item{newloglikSeq}{log-likelihood from new sequence.}
#' \item{step}{depending on the type of step, result of [stepAugment()],[stepShort()] or [stepPerm()].}
#' \item{pUndoStep}{probability of undoing the step of MH acceptance rate.}
#'
#' @export
#'
stepPT <- function(seq, type, actDfnodesLab, tieNames, formula, net0, beta, theta, initTime, endTime, k = 5, temp = 1,pAug,pShort) {
  getKelMeMatrix <- getKelMeMatrix(seq, actDfnodesLab)
  Kel_g1 <- getKelMeMatrix$Kel_g1
  Kel_ge1 <- getKelMeMatrix$Kel_ge1
  gammaEminus <- choose(Kel_ge1, 2) + Kel_g1
  gammaMinus <- sum(gammaEminus)
  me <- getKelMeMatrix$me
  gammaEplus <- choose(nrow(seq) - me + 2, 2)
  gammaPlus <- sum(gammaEplus)
  m <- nrow(seq)
  auxDf <- getKelMeMatrix$auxDf


  if (type == 1) { # Augmentation
    step <- stepAugment(seq, tieNames, gammaEplus, gammaPlus, m, me, net0,pAug)

    getNewKelMeMatrix <- getKelMeMatrix(step$newseq, actDfnodesLab)
    newAuxDfE <- getAuxDfE(getNewKelMeMatrix$auxDf, step$sender, step$receiver)
    auxDfE <- getAuxDfE(auxDf, step$sender, step$receiver)

    newKel_g1 <- getNewKelMeMatrix$Kel_g1
    newKel_ge1 <- getNewKelMeMatrix$Kel_ge1
    newGammaEminus <- choose(newKel_ge1, 2) + newKel_g1
    newGammaMinus <- sum(newGammaEminus)
    newMe <- getNewKelMeMatrix$me
    newGammaEplus <- choose(nrow(newseq) - newMe + 2, 2)
    newGammaPlus <- sum(newGammaEplus)
    newM <- nrow(newseq)

    newp <- newKel_g1[step$sender, step$receiver] / newGammaEminus[step$sender, step$receiver]
    if ((length(unique(auxDfE$run)) == length(unique(newAuxDfE$run))) | step$typeA == "same") {
      pUndoStep <- getpDoShort(newGammaEminus, newGammaMinus, newp, length(unique(newAuxDfE$run)), step$sender, step$receiver, "same",pShort)
    } else {
      pUndoStep <- getpDoShort(newGammaEminus, newGammaMinus, newp, length(unique(newAuxDfE$run)), step$sender, step$receiver, "diff",pShort)
    }

    for (i in 1:k) {
      step <- stepPerm(step$newseq, tieNames, newM, newMe)
      newMe <- getNewKelMeMatrix$me
      newGammaEplus <- choose(nrow(newseq) - newMe + 2, 2)
      newGammaPlus <- sum(newGammaEplus)
      newM <- nrow(newseq)
    }

    loglikSeq <- logLikelihoodMCTemp(indexCore = 1, list(seq[, -which(colnames(seq) == "row")]), beta, theta, initTime, endTime, splitIndicesPerCore = list(1), actDfnodes = actDfnodes, net0 = net0, formula = formula, temp)
    newloglikSeq <- logLikelihoodMCTemp(indexCore = 1, list(step$newseq[, -which(colnames(step$newseq) == "row")]), beta, theta, initTime, endTime, splitIndicesPerCore = list(1), actDfnodes = actDfnodes, net0 = net0, formula = formula, temp)
  } else if (type == 2) { # Shortening
    step <- stepShort(seq, tieNames, gammaEminus, gammaMinus, m, me, Kel_g1, auxDf,pShort)
    getNewKelMeMatrix <- getKelMeMatrix(step$newseq, actDfnodesLab)
    newAuxDfE <- getAuxDfE(getNewKelMeMatrix$auxDf, sender, receiver)
    auxDfE <- getAuxDfE(auxDf, sender, receiver)

    newKel_g1 <- getNewKelMeMatrix$Kel_g1
    newKel_ge1 <- getNewKelMeMatrix$Kel_ge1
    newGammaEminus <- choose(newKel_ge1, 2) + newKel_g1
    newGammaMinus <- sum(newGammaEminus)
    newMe <- getNewKelMeMatrix$me
    newGammaEplus <- choose(nrow(newseq) - newMe + 2, 2)
    newGammaPlus <- sum(newGammaEplus)
    newM <- nrow(newseq)

    newp <- (newM - newMe[step$sender, step$receiver] + 1) / newGammaEplus[step$sender, step$receiver]

    newindexSR <- which(step$newseq$sender == step$sender & step$newseq$receiver == step$receiver)

    if (step$typeS == "same") {
      pUndoStep <- getpDoAugment(newGammaEplus, newGammaPlus, newp, newM, newMe, step$sender, step$receiver, "same", newindexSR,pAug)
    } else if (step$typeS == "diff") {
      pUndoStep <- getpDoAugment(newGammaEplus, newGammaPlus, newp, newM, newMe, step$sender, step$receiver, "diff", newindexSR,pAug)
    }

    for (i in 1:k) {
      step <- stepPerm(step$newseq, tieNames, newM, newMe)
      newMe <- getNewKelMeMatrix$me
      newGammaEplus <- choose(nrow(newseq) - newMe + 2, 2)
      newGammaPlus <- sum(newGammaEplus)
      newM <- nrow(newseq)
    }

    loglikSeq <- logLikelihoodMCTemp(indexCore = 1, list(seq[, -which(colnames(seq) == "row")]), beta, theta, initTime, endTime, splitIndicesPerCore = list(1), actDfnodes = actDfnodes, net0 = net0, formula = formula, temp)
    newloglikSeq <- logLikelihoodMCTemp(indexCore = 1, list(step$newseq[, -which(colnames(step$newseq) == "row")]), beta, theta, initTime, endTime, splitIndicesPerCore = list(1), actDfnodes = actDfnodes, net0 = net0, formula = formula, temp)
  }


  u <- runif(1, min = 0, max = 1)
  accept <- (u <= exp(newloglikSeq - loglikSeq) * pUndoStep / pDoStep)
  if (!accept) {
    acceptIndex <- acceptIndex + 1
    # if(accept & !acceptIndex %% 50){
    newseq <- seq
    newloglikSeq <- loglikSeq
  }

  return(list(
    "newseq" = newseq, "newloglikSeq" = newloglikSeq,
    "step" = step, "accept" = accept
  ))
}



#' Parallel tempering step (multi-core)
#'
#' @description Performs an step of MCMC given a type:
#' \item Type 1: augmentation of sequence.
#' \item Type 2: shortening of sequence.
#'
#' @param seqs list of data frames, sequences of events for each temperature.
#' @param type vector of integers, type of step for each chain.
#' @param actDfnodesLab vector, node labs.
#' @param tieNames vector, tie labs (e.g. "12" is a tie from sender 1 to receiver 2).
#' @param formula formula, the one used in goldfish model.
#' @param net0 matrix, initial observed network.
#' @param beta list, estimator of choice parameters of the model (creation and deletion)
#' @param theta dataframe, estimator of rate parameters of the model (creation and deletion) (constant)
#' @param initTime initial time
#' @param endTime end time
#' @param T0 maximum temperature
#'
#' @return list of result of stepMCMC for each sequence with different temperatures.
#'
#' @export
#'
stepPTMC <- function(indexCore, splitIndicesPerCore, seqs, type, actDfnodesLab, tieNames, formula, net0, beta, theta, initTime, endTime, k, temp, nStepSwitch,
                     pAug,pShort) {

  indicesCore <- splitIndicesPerCore[[indexCore]]
  stepPT <- vector("list", length(indicesCore))
  for (i in seq_along(indicesCore)) {
    stepPT[[i]] <- lapply(1:nStepSwitch, stepPT(seqs[[indicesCore[[i]]]], type[indicesCore[[i]]], actDfnodesLab, tieNames, formula, net0, beta, theta, initTime, endTime, k, temp[indicesCore[[i]]],
                          pAug,pShort))
  }
  return(stepPT)
}





#' MCMC function (simple version with augmentation/shortening/permutation)
#'
#' @description Computes an MCMC chain and returns nmax sequences from the MCMC chain.
#'
#' @param nmax integer, number of sequences to be returned.
#' @param seq data frame, sequence from which start/continue MCMC chain.
#' @param H integer, hamming distance between net0 and net1.
#' @param actDfnodes object of type ´nodes´ from goldfish.
#' @param formula formula, the one used in goldfish model.
#' @param net0 matrix, initial observed network.
#' @param beta list, estimator of parameters of the model (creation and deletion).
#' @param burnIn boolean, indicates if burn in must be performed.
#' @param maxIter integer, maximum number of steps of the MCMC chain.
#' @param seqIter integer, number of steps between selected sequences.
#' @param pShort float, probability of shortening step.
#' @param pAugment float, probability of augmenting step.
#' @param pPerm float, probability of permutation step.
#'
#' @return permut, list of sequences.
#'
#' @export
#'
MCMC <- function(nmax, seq, H, actDfnodes, formula, net0, beta, theta, initTime, endTime, burnIn = TRUE, maxIter = 10000, seqIter = 50, pShort = 0.35, pAug = 0.35, pPerm = 0.3) {
  # Compute initial quatities:
  # Type 1: augmentation, type 2: shortening, type 3: permutation

  if (nmax > (maxIter - 500) / seqIter) {
    maxIter <- (nmax + 501) * 50
  }

  actDfnodesLab <- actDfnodes$label

  tieNames <- sapply(actDfnodesLab, function(x) sapply(actDfnodesLab, function(i) paste(x, i, sep = "-")))
  tieNames <- as.vector(tieNames)

  permut <- vector(mode = "list", length = nmax)

  acceptIndex <- 0
  for (i in 1:maxIter) {
    if (nrow(seq) == H) {
      type <- sampleVec(c(1, 2, 3), size = 1, prob = c(pAug / (pAug + pPerm), 0, pPerm / (pAug + pPerm)))
    } else {
      type <- sampleVec(c(1, 2, 3), size = 1, prob = c(pAug, pShort, pPerm))
    }

    seq$row <- 1:nrow(seq)
    stepMCMC <- stepMCMC(seq, type, actDfnodesLab, tieNames, formula, net0, beta, theta, initTime, endTime)

    u <- runif(1, min = 0, max = 1)
    accept <- (u <= exp(stepMCMC$newloglikSeq - stepMCMC$loglikSeq) * stepMCMC$pUndoStep / stepMCMC$step$pDoStep)
    if (accept) {
      acceptIndex <- acceptIndex + 1
      # if(accept & !acceptIndex %% 50){
      seq <- stepMCMC$newseq
    }
  }

  if (burnIn) {
    if (i > 500 & !i %% 50) {
      permut <- c(permut, seq)
    }
  } else {
    if (!i %% 50) {
      permut <- c(permut, seq)
    }
  }

  return(permut)
}

#' Parallel tempering MCMC function (augmentation/shortening)
#'
#' @description Computes an MCMC chain and returns nmax sequences from the MCMC chain.
#'
#' @param nmax integer, number of sequences to sampled.
#' @param nPT integer, number of chains in the parallel tempering algorithm.
#' @param seqsInit list of data frame, sequences from which start/continue MCMC chain. They must have all the same length.
#' @param H integer, hamming distance between net0 and net1.
#' @param actDfnodes object of type ´nodes´ from goldfish.
#' @param formula formula, the one used in goldfish model.
#' @param net0 matrix, initial observed network.
#' @param beta list, estimator of parameters of the model (creation and deletion).
#' @param burnIn boolean, indicates if burn in must be performed.
#' @param maxIter integer, maximum number of steps of the MCMC chain.
#' @param seqIter integer, number of steps between selected sequences.
#' @param pShort float, probability of shortening step.
#' @param pAugment float, probability of augmenting step.
#' @param T0 initial maximal tempetature for the parallale tempering scheme/simulated annealing
#' @param nStepExch number of mutation steps before a exchange
#'
#' @return permut, list of sequences.
#'
#' @export
#'
PT_MCMC <- function(nmax, nPT, seqsInit, H, actDfnodes, formula, net0, beta, theta, initTime, endTime, burnIn = TRUE, maxIter = 10000, seqIter = 50, pShort = 0.5, pAug = 0.5, T0 = 100, nStepExch = 10, num_cores = 1) {
  # Compute initial quatities:
  # Type 1: augmentation, type 2: shortening

  if (nmax > (maxIter - 500) / seqIter) {
    maxIter <- (nmax + 501) * 50
  }

  actDfnodesLab <- actDfnodes$label

  tieNames <- sapply(actDfnodesLab, function(x) sapply(actDfnodesLab, function(i) paste(x, i, sep = "-")))
  tieNames <- as.vector(tieNames)


  if(lenght(seqsPT)==1){
    seqsPT = permute(seqsPT, nmax = nPT)
  }
  permute <- list()

  temp <- seq(1, T0, length = length(seqs))

  acceptIndex <- 0
  for (i in 1:maxIter) {
    if (nrow(seqsPT[[1]]) == H) {
      type <- sampleVec(c(1, 2), size = length(seqsPT), prob = c(1, 0), replace = TRUE)
    } else {
      type <- sampleVec(c(1, 2), size = length(seqsPT), prob = c(pAug, pShort), replace = TRUE)
    }

    seqsPT <- lapply(seqsPT, function(x) cbind(x, "row" = 1:nrow(x)))


    splitIndicesPerCore <- splitIndices(length(seqsPT), num_cores)
    cl <- makeCluster(num_cores)
    on.exit(stopCluster(cl))
    clusterExport(cl, c(
      "seqs", "type", "actDfnodesLab", "tieNames",
      "formula", "net0", "beta", "theta", "initTime",
      "endTime", "k", "T0","nStepExch","pAug","pShort"
    ))
    clusterEvalQ(cl, {
      library(goldfish)
      NULL
    })


    clusterExport(cl, list(
      "stepPT", "stepPTMC", "stepMCMC", "getpDoAugment", "getpDoShort", "petpDoPerm",
      "getKelMeMatrix", "getAuxDfE", "stepAugment", "stepShort", "stepPerm",
      "logLikelihoodMCTemp", "GatherPreprocessingDF","logLikChoice","logLikConstantRate"
    ))


    stepPT <- clusterApply(cl, seq_along(splitIndicesPerCore), stepPT,
      splitIndicesPerCore = splitIndicesPerCore,
      seqs = seqsPT, type = type, actDfnodesLab = actDfnodesLab,
      tieNames = tieNames, formula = formula, net0 = net0,
      beta = beta, theta = theta, initTime = initTime,
      endTime = endTime, k = k, temp = temp, nStepExch = nStepExch,
      pAug = pAug, pShort = pShort
    )

    # Every 10 iterations, try to switch
    order = sample(c(0,1),size=1,prob=c(0.5,0.5)) # True=1: Pairs of sequences 1-2, 3-4, ....
                                                  # False=0: Paris of sequences 1, 2-3, 4-5, ...
    if(order){
      if(! nPT %% 2) {
        tempPairs = data.frame(seq(1,nPT,by=2),seq(2,nPT,by=2))
      }else{
        tempPairs = data.frame(seq(1,nPT-1,by=2),seq(2,nPT,by=2))
      }
    }else{
      if(! nPT %% 2) {
        tempPairs = data.frame(seq(2,nPT-1,by=2),seq(3,nPT,by=2))
      }else{
        tempPairs = data.frame(seq(2,nPT,by=2),seq(3,nPT,by=2))
      }
    }

    # Exchange step!
    # Accept or reject the change
    logUnif = log(runif(nrow(tempPairs)))
    for( j in 1:nrow(tempPairs)){
      logMHRatio = min(0,stepPT[[tempPairs[j,1]]]$newloglikSeq*temp[tempPairs[j,1]]/temp[tempPairs[j,2]]+
                stepPT[[tempPairs[j,2]]]$newloglikSeq*temp[tempPairs[j,2]]/temp[tempPairs[j,1]]-
                stepPT[[tempPairs[j,1]]]$newloglikSeq-stepPT[[tempPairs[j,2]]]$newloglikSeq)
      if(logUnif[j]<logMHRatio){
        aux = steptPT[[tempPairs[j,1]]]
        steptPT[[tempPairs[j,1]]] = steptPT[[tempPairs[j,2]]]
        steptPT[[tempPairs[j,2]]] = aux
        rm(aux)
      }
    }


    if (burnIn) { # avoid, after burn in, to have switch and sample in the same step (not really needed)
      if ((i-5) > 500 & !(i-5) %% seqIter) {
        permut <- c(permut, stepPT[[1]]) # get a sample from the sequence with temperature 1
      }
    } else {
      if (!(i-5) %% seqIter) {
        permut <- c(permut, stepPT[[1]])
      }
    }
  }

  return(permut,stepPT)

}




# MCMC_MC = function(indexCore,splitIndicesPerCore,permut = permut,beta=beta, burnIn = TRUE, H = H,actDfnodes=actDfnodes){
#
#   indicesCore = splitIndicesPerCore[[indexCore]]
#   resMCMC = vector("list",length(indicesCore))
#
#   if(burnIn) {burn_in = burnIn(permut[[1]],beta,H,actDfnodes)}
#
#
#   for(i in seq_along(indicesCore)){
#     resMCMC[[i]] = MCMC(permut[[indicesCore[[i]]]],burn_in,H,actDfnodes)
#   }
#
#   return(resMCMC)
# }
