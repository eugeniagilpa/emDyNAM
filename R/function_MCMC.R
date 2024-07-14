# MCMC auxiliar functions --------------------


#' Get probability of doing a augmenting step.
#'
getpDoAugment <- function(gammaEplus, gammaPlus, p, m, me, sender, receiver,
                          typeA, indexSR, pAug) {
  if (typeA == "same") {
    if (length(indexSR) < 1) { #   indexSR <- which(seq$sender == sender & seq$receiver == receiver)
      pDoStep <- (gammaEplus[sender, receiver] / gammaPlus) * p *
        (1 / (m + 1)) * pAug
    } else {
      pDoStep <- (gammaEplus[sender, receiver] / gammaPlus) * p *
        (1 / (m - me[sender, receiver] + 1)) * pAug
    }
  } else {
    pDoStep <- (gammaEplus[sender, receiver] / gammaPlus) * (1 - p) *
      (1 / (m - me[sender, receiver] + 1)) * (1 / (m - me[sender, receiver])) *
      pAug
  }

  return(pDoStep)
}

#' Get probability of doing a shortening step.
#'
getpDoShort <- function(gammaEminus, gammaMinus, p, uniqueR, sender, receiver,
                        typeS, pShort) {
  if (typeS == "same") {
    pDoStep <- (gammaEminus[sender, receiver] / gammaMinus) * p *
      (1 / uniqueR) * pShort
  } else {
    pDoStep <- (gammaEminus[sender, receiver] / gammaMinus) * (1 - p) *
      (1 / uniqueR) * (1 / (uniqueR - 1)) * pShort
  }

  return(pDoStep)
}

#' Get probability of doing a permutation step.
#'
getpDoPerm <- function(probVec, probVec2, e1, e2, sender1, receiver1, sender2,
                       receiver2, me) {


  pDoStep <- (probVec[e1] * probVec2[e2] + probVec[e2] * probVec2[e1]) *
    (1/me[sender1,receiver1]) * (1/me[sender2,receiver2])

  # pDoStep <- probVec[e1] * probVec2[e2] *
  #   (1 / nrow(seq[seq$sender == sender1 & seq$receiver == receiver1, ])) *
  #   (1 / nrow(seq[seq$sender == sender2 & seq$receiver == receiver2, ]))
  return(pDoStep)
}


#' Get auxiliary matrix \eqn{K_{el}} and \eqn{m_e} for MCMC steps.
#'
#' @param seq data frame with the sequence of events up to time \eqn{t}.
#' @param actDfnodes vector with the labels of the actors of the network.
#'
#' @return list with matrix \eqn{K_{el}}, \eqn{K_{e,l>1}}, \eqn{m_e} and
#' auxiliary data frame
#' @export
#'
#' @examples seq <- data.frame(
#'   "sender" = c(1, 2, 3), "receiver" = c(3, 2, 1),
#'   "replace" = c(0, 1, 1)
#' )
#' getKelMeMatrix(seq, actDfnodes = c(1, 2, 3))
getKelMeMatrix <- function(seq, actDfnodesLab) {
  if (!"row" %in% colnames(seq)) {
    seq <- cbind(seq, "row" = 1:nrow(seq))
  }
# browser()

  auxDf <- lapply(actDfnodesLab, function(x) {
    lapply(actDfnodesLab, function(i) {
      subset(seq, subset = (sender == x & receiver == i))
    })
  })

  for (i in 1:length(actDfnodesLab)) {
    for (j in 1:length(actDfnodesLab)) {
      if (nrow(auxDf[[i]][[j]]) > 0) {
        auxDf[[i]][[j]]$row <- as.integer(auxDf[[i]][[j]]$row)
        auxDf[[i]][[j]]$rowDiff <- c(0, diff(auxDf[[i]][[j]]$row))
      }
    }
  }


  Kel_g1 <- matrix(0, nrow = length(actDfnodesLab), ncol = length(actDfnodesLab)) # Notation from Johan´s paper, sum for l greater than 1
  Kel_ge1 <- matrix(0, nrow = length(actDfnodesLab), ncol = length(actDfnodesLab)) # Notation from Johan´s paper, sum for l greater or equal than 1

  for (i in 1:nrow(Kel_g1)) {
    for (j in 1:ncol(Kel_g1)) {
      # browser()
      auxDfE <- getAuxDfE(auxDf, i, j)
      table <- table(auxDfE$run)
      # cat("\n i ",i," j ",j,"\n")
      # print(table)
      Kel_ge1[i, j] <- length(table)
      Kel_g1[i, j] <- sum(table > 1)

      # auxKel = 1
      # Kel_g1[i, j] = 0
      # index = 2
      # while(auxKel > 0){
      #   auxKel = sum(table>=index)
      #   if(auxKel == 0 ) break
      #   auxIndex = which(table>=index)
      #   Kel_g1[i, j] <- Kel_g1[i, j] + sum(table[auxIndex] - (index -1))
      #   index = index + 1
      # }

      # Kel_ge1[i, j] <- sum(table) + Kel_g1[i,j]

      #   if (nrow(auxDf[[i]][[j]]) > 0) {
      #     Kel_ge1[i, j] <- nrow(auxDf[[i]][[j]][auxDf[[i]][[j]]$rowDiff != 1, ])
      #     v_rowDiff <- auxDf[[i]][[j]]$rowDiff
      #     not_1 <- which(v_rowDiff != 1)
      #     shifted_vector <- c(v_rowDiff[-1], 0)
      #     # Kel_g1[i, j] <- sum(v_rowDiff[-length(v_rowDiff)] != 1 & shifted_vector[-1] == 1)
      #   }
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
  # browser()
  auxDfE <- auxDf[[sender]][[receiver]]
  indexNo1 <- which(auxDfE$rowDiff != 1)
  auxDfE$run <- rep(0, nrow(auxDfE))
  run <- 1
  if (length(indexNo1) > 1) {
    for (i in 1:(length(indexNo1) - 1)) {
      auxDfE$run[indexNo1[i]:indexNo1[i + 1]] <- run
      run <- run + 1
    }
    auxDfE$run[indexNo1[length(indexNo1)]:nrow(auxDfE)] <- run
  } else {
    if (length(indexNo1) == 1) {
      auxDfE$run[indexNo1:nrow(auxDfE)] <- 1
    }
  }
  return(auxDfE)
}

#' Step for augmenting a given sequence by adding a tie \eqn{e} two times.
#'
#' @param seq data frame, sequence of events.
#' @param tieNames vector, labels of ties (e.g., "12" if the tie is from
#' sender 1 and receiver 2).
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
stepAugment <- function(seq, tieNames, gammaEplus, gammaPlus, m, me, net0, pAug) {
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
      if (place == (m + 1)) {
        newseq <- rbind(
          seq[1:(place - 1), ],
          c(sender, receiver, 4, place),
          c(sender, receiver, 4, place + 1)
        )
      } else {
        newseq <- rbind(
          seq[1:(place - 1), ],
          c(sender, receiver, 4, place),
          c(sender, receiver, 4, place + 1),
          seq[(place):m, ]
        )
      }
    }

    nEdges <- which(newseq$sender == sender & newseq$receiver == receiver)
    initAction <- as.integer(seq[
      which(seq$sender == sender & seq$receiver == receiver),
      "replace"
    ][1])
    if (is.na(initAction)) {
      initAction <- net0[sender, receiver] + 1
    }
    newseq[nEdges, "replace"] <- seq(initAction, (initAction + length(nEdges) - 1)) %% 2
    newseq$row <- seq(1, nrow(newseq))

    pDoStep <- getpDoAugment(
      gammaEplus, gammaPlus, p, m, me, sender, receiver,
      typeA, indexSR, pAug
    )
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
    initAction <- as.integer(seq[
      which(seq$sender == sender & seq$receiver == receiver),
      "replace"
    ][1])
    if (is.na(initAction)) {
      initAction <- net0[sender, receiver] + 1
    }
    newseq[nEdges, "replace"] <- seq(initAction, (initAction + length(nEdges) - 1)) %% 2
    newseq$row <- seq(1, nrow(newseq))

    pDoStep <- getpDoAugment(
      gammaEplus, gammaPlus, p, m, me, sender, receiver,
      typeA, indexSR, pAug
    )
  }

  return(list(
    newseq = newseq, sender = sender, receiver = receiver, place = place,
    typeA = typeA, pDoStep = pDoStep
  ))
}


#' Step for shortening a given sequence by removing a tie \eqn{e} in
#' two positions.
#'
#' @param seq data frame, sequence of events.
#' @param tieNames vector, labels of ties (e.g., "12" if the tie is from
#' sender 1 and receiver 2).
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
stepShort <- function(seq, tieNames, gammaEminus, gammaMinus, m, me, Kel_g1,
                      auxDf, pShort) {
  # Choose element to be deleted:
  # browser
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

    pDoStep <- getpDoShort(
      gammaEminus, gammaMinus, p,
      length(which(table(auxDfE$run) > 1)), sender, receiver,
      typeS, pShort
    )
  } else if (typeS == "diff") {
    # Choose two different e-runs, and delete the first element of each one of them
    sampleRun <- sampleVec(unique(auxDfE$run), size = 2, replace = FALSE)
    indexSampleRun <- c(
      auxDfE$row[auxDfE$run == sampleRun[1]][1],
      auxDfE$row[auxDfE$run == sampleRun[2]][1]
    )

    newseq <- seq[-(indexSampleRun), ]
    newseq$row <- 1:nrow(newseq)

    # now we have to fix the replaceOG
    nEdges <- which(newseq$sender == sender & newseq$receiver == receiver)
    initAction <- as.integer(seq[which(seq$sender == sender & seq$receiver == receiver), "replace"][1])
    newseq[nEdges, "replace"] <- seq(initAction, (initAction + length(nEdges) - 1)) %% 2

    place <- list(sampleRun = sampleRun, indexSampleRun = indexSampleRun)
    pDoStep <- getpDoShort(
      gammaEminus, gammaMinus, p,
      length(unique(auxDfE$run)), sender, receiver,
      typeS, pShort
    )
  }

  return(list(
    newseq = newseq, sender = sender, receiver = receiver, place = place,
    typeS = typeS, pDoStep = pDoStep, auxDfE = auxDfE
  ))
}



#' Step for permuting two elements in a given sequence.
#'
#' @param seq data frame, sequence of events.
#' @param tieNames vector, labels of ties (e.g., "12" if the tie is from
#' sender 1 and receiver 2).
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
  initAction <- as.integer(seq[
    which(seq$sender == sender1 & seq$receiver == receiver1),
    "replace"
  ][1])
  newseq[nEdges, "replace"] <- seq(initAction, (initAction + length(nEdges) - 1)) %% 2
  nEdges <- which(newseq$sender == sender2 & newseq$receiver == receiver2)
  initAction <- as.integer(seq[
    which(seq$sender == sender2 & seq$receiver == receiver2),
    "replace"
  ][1])
  newseq[nEdges, "replace"] <- seq(initAction, (initAction + length(nEdges) - 1)) %% 2

  newseq$row <- 1:nrow(newseq)

  pDoStep <- getpDoPerm(probVec, probVec2, e1, e2, sender1, receiver1, sender2,
                                    receiver2, me)
  # pDoStep <- probVec[e1] * probVec2[e2] *
  #   (1 / nrow(seq[seq$sender == sender1 & seq$receiver == receiver1, ])) *
  #   (1 / nrow(seq[seq$sender == sender2 & seq$receiver == receiver2, ]))

  return(list(
    newseq = newseq, sender = c(sender1, sender2),
    receiver = c(receiver1, receiver2),
    place = c(place1, place2),
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
#' @param tieNames vector, tie labs (e.g. "12" is a tie from sender 1
#' to receiver 2).
#' @param formula formula, the one used in goldfish model.
#' @param net0 matrix, initial observed network.
#' @param beta list, estimator of choice parameters of the model
#' (creation and deletion)
#' @param theta dataframe, estimator of rate parameters of the model
#' (creation and deletion) (constant)
#' @param initTime initial time
#' @param endTime end time
#'
#' @return list of
#' \item{newseq}{new sequence.}
#' \item{loglikSeq}{log-likelihood from original sequence.}
#' \item{newloglikSeq}{log-likelihood from new sequence.}
#' \item{step}{depending on the type of step, result of [stepAugment()],
#' [stepShort()] or [stepPerm()].}
#' \item{pUndoStep}{probability of undoing the step of MH acceptance rate.}
#'
#' @export
#'
stepMCMC <- function(seq, type, actDfnodesLab, actDfnodes, tieNames, formula,
                     net0, beta, theta, initTime, endTime, pAug, pShort) {
  getKelMeMatrix <- getKelMeMatrix(seq, actDfnodesLab)
  Kel_g1 <- getKelMeMatrix$Kel_g1
  Kel_ge1 <- getKelMeMatrix$Kel_ge1
  gammaEminus <- choose(Kel_ge1, 2) + Kel_g1
  gammaMinus <- sum(gammaEminus)
  me <- getKelMeMatrix$me
  gammaEplus <- choose(nrow(seq) - me + 2, 2)
  diag(gammaEplus) <- 0
  gammaPlus <- sum(gammaEplus)
  m <- nrow(seq)
  auxDf <- getKelMeMatrix$auxDf



  if (type == 1) { # Augmentation
    step <- stepAugment(seq, tieNames, gammaEplus, gammaPlus, m, me, net0, pAug)

    getNewKelMeMatrix <- getKelMeMatrix(step$newseq, actDfnodesLab)
    newAuxDfE <- getAuxDfE(getNewKelMeMatrix$auxDf, step$sender, step$receiver)
    auxDfE <- getAuxDfE(auxDf, step$sender, step$receiver)

    newKel_g1 <- getNewKelMeMatrix$Kel_g1
    newKel_ge1 <- getNewKelMeMatrix$Kel_ge1
    newGammaEminus <- choose(newKel_ge1, 2) + newKel_g1
    newGammaMinus <- sum(newGammaEminus)
    newMe <- getNewKelMeMatrix$me
    newGammaEplus <- choose(nrow(newseq) - newMe + 2, 2)
    diag(newGammaEplus) <- 0
    newGammaPlus <- sum(newGammaEplus)
    newM <- nrow(newseq)

    newp <- newKel_g1[step$sender, step$receiver] /
      newGammaEminus[step$sender, step$receiver]
    if ((length(unique(auxDfE$run)) == length(unique(newAuxDfE$run))) | step$typeA == "same") {
      pUndoStep <- getpDoShort(
        newGammaEminus, newGammaMinus, newp,
        length(which(table(newAuxDfE$run) > 1)), step$sender,
        step$receiver, "same", pShort
      )
    } else {
      pUndoStep <- getpDoShort(
        newGammaEminus, newGammaMinus, newp,
        length(unique(newAuxDfE$run)), step$sender,
        step$receiver, "diff", pShort
      )
    }

    loglikSeq <- logLikelihoodMC(
      indexCore = 1,
      list(seq[, -which(colnames(seq) == "row")]),
      beta, theta, initTime, endTime,
      splitIndicesPerCore = list(1),
      actDfnodes = actDfnodes, net0 = net0,
      formula = formula
    )
    newloglikSeq <- logLikelihoodMC(
      indexCore = 1,
      list(step$newseq[, -which(colnames(step$newseq) == "row")]),
      beta, theta, initTime, endTime,
      splitIndicesPerCore = list(1),
      actDfnodes = actDfnodes, net0 = net0,
      formula = formula
    )
  } else if (type == 2) { # Shortening
    step <- stepShort(
      seq, tieNames, gammaEminus, gammaMinus, m, me, Kel_g1,
      auxDf, pShort
    )
    getNewKelMeMatrix <- getKelMeMatrix(step$newseq, actDfnodesLab)
    newAuxDfE <- getAuxDfE(getNewKelMeMatrix$auxDf, sender, receiver)
    auxDfE <- getAuxDfE(auxDf, sender, receiver)

    newKel_g1 <- getNewKelMeMatrix$Kel_g1
    newKel_ge1 <- getNewKelMeMatrix$Kel_ge1
    newGammaEminus <- choose(newKel_ge1, 2) + newKel_g1
    newGammaMinus <- sum(newGammaEminus)
    newMe <- getNewKelMeMatrix$me
    newGammaEplus <- choose(nrow(newseq) - newMe + 2, 2)
    diag(newGammaEplus) <- 0
    newGammaPlus <- sum(newGammaEplus)
    newM <- nrow(newseq)

    newp <- (newM - newMe[step$sender, step$receiver] + 1) /
      newGammaEplus[step$sender, step$receiver]

    newindexSR <- which(step$newseq$sender == step$sender &
      step$newseq$receiver == step$receiver)

    if (step$typeS == "same") {
      pUndoStep <- getpDoAugment(
        newGammaEplus, newGammaPlus, newp, newM, newMe,
        step$sender, step$receiver, "same", newindexSR,
        pAug
      )
    } else if (step$typeS == "diff") {
      pUndoStep <- getpDoAugment(
        newGammaEplus, newGammaPlus, newp, newM, newMe,
        step$sender, step$receiver, "diff", newindexSR,
        pAug
      )
    }

    loglikSeq <- logLikelihoodMC(
      indexCore = 1,
      list(seq[, -which(colnames(seq) == "row")]),
      beta, theta, initTime, endTime,
      splitIndicesPerCore = list(1),
      actDfnodes = actDfnodes, net0 = net0,
      formula = formula
    )
    newloglikSeq <- logLikelihoodMC(
      indexCore = 1,
      list(step$newseq[, -which(colnames(step$newseq) == "row")]),
      beta, theta, initTime, endTime,
      splitIndicesPerCore = list(1),
      actDfnodes = actDfnodes, net0 = net0,
      formula = formula
    )
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
    diag(newGammaEplus) <- 0
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
    loglikSeq <- logLikelihoodMC(
      indexCore = 1,
      list(seq[, -which(colnames(seq) == "row")]),
      beta, theta, initTime, endTime,
      splitIndicesPerCore = list(1),
      actDfnodes = actDfnodes, net0 = net0,
      formula = formula
    )
    newloglikSeq <- logLikelihoodMC(
      indexCore = 1,
      list(step$newseq[, -which(colnames(step$newseq) == "row")]),
      beta, theta, initTime, endTime,
      splitIndicesPerCore = list(1),
      actDfnodes = actDfnodes, net0 = net0,
      formula = formula
    )
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
#' @param tieNames vector, tie labs (e.g. "12" is a tie from sender 1
#' to receiver 2).
#' @param formula formula, the one used in goldfish model.
#' @param net0 matrix, initial observed network.
#' @param beta list, estimator of choice parameters of the model
#' (creation and deletion)
#' @param theta dataframe, estimator of rate parameters of the model
#' (creation and deletion) (constant)
#' @param initTime initial time
#' @param endTime end time
#' @param k number of permutations after augmentations/shortening
#'
#' @return list of
#' \item{newseq}{new sequence.}
#' \item{loglikSeq}{log-likelihood from original sequence.}
#' \item{newloglikSeq}{log-likelihood from new sequence.}
#' \item{step}{depending on the type of step, result of [stepAugment()],
#' [stepShort()] or [stepPerm()].}
#' \item{pUndoStep}{probability of undoing the step of MH acceptance rate.}
#'
#' @export
#'
stepPT <- function(seq, type, actDfnodesLab, actDfnodes, tieNames, formula, net0,
                   beta, theta, fixedparameters, initTime, endTime,
                   k = 5, temp = 1,
                   pAug, pShort, pPerm, logLikelihoodStats) {
  getKelMeMatrix <- getKelMeMatrix(seq, actDfnodesLab)
  Kel_g1 <- getKelMeMatrix$Kel_g1
  Kel_ge1 <- getKelMeMatrix$Kel_ge1
  gammaEminus <- choose(Kel_ge1, 2) + Kel_g1
  gammaMinus <- sum(gammaEminus)
  if (gammaMinus == 0) {
    type <- 1
  }
  me <- getKelMeMatrix$me
  gammaEplus <- choose(nrow(seq) - me + 2, 2)
  diag(gammaEplus) <- 0
  gammaPlus <- sum(gammaEplus)
  m <- nrow(seq)
  auxDf <- getKelMeMatrix$auxDf

  if (type == 1) { # Augmentation
    step <- stepAugment(seq, tieNames, gammaEplus, gammaPlus, m, me, net0, pAug)

    getNewKelMeMatrix <- getKelMeMatrix(step$newseq, actDfnodesLab)
    newAuxDfE <- getAuxDfE(getNewKelMeMatrix$auxDf, step$sender, step$receiver)
    auxDfE <- getAuxDfE(auxDf, step$sender, step$receiver)

    newKel_g1 <- getNewKelMeMatrix$Kel_g1
    newKel_ge1 <- getNewKelMeMatrix$Kel_ge1
    newGammaEminus <- choose(newKel_ge1, 2) + newKel_g1
    newGammaMinus <- sum(newGammaEminus)
    newMe <- getNewKelMeMatrix$me
    newGammaEplus <- choose(nrow(step$newseq) - newMe + 2, 2)
    diag(newGammaEplus) <- 0
    newGammaPlus <- sum(newGammaEplus)
    newM <- nrow(step$newseq)

    pDoStep <- step$pDoStep
    newp <- newKel_g1[step$sender, step$receiver] /
      newGammaEminus[step$sender, step$receiver]
    if ((length(unique(auxDfE$run)) == length(unique(newAuxDfE$run))) | step$typeA == "same") {
      pUndoStep <- getpDoShort(
        newGammaEminus, newGammaMinus, newp,
        length(which(table(newAuxDfE$run) > 1)), step$sender,
        step$receiver, "same", pShort
      )
    } else {
      pUndoStep <- getpDoShort(
        newGammaEminus, newGammaMinus, newp,
        length(unique(newAuxDfE$run)), step$sender,
        step$receiver, "diff", pShort
      )
    }

    for (i in 1:k) {
      step <- stepPerm(step$newseq, tieNames, newM, newMe)

      getNewKelMeMatrix <- getKelMeMatrix(step$newseq, actDfnodesLab)
      newMe <- getNewKelMeMatrix$me
      newGammaEplus <- choose(nrow(step$newseq) - newMe + 2, 2)
      diag(newGammaEplus) <- 0
      newGammaPlus <- sum(newGammaEplus)
      newM <- nrow(step$newseq)
    }



    newlogLikelihoodStats <- getlogLikelihood(
      step$newseq, actDfnodes, net0, fixedparameters,
      beta, initTime, endTime, formula
    )

    newloglikSeq <- (newlogLikelihoodStats$resCrea$logLikelihood +
      newlogLikelihoodStats$resDel$logLikelihood) / temp

    # loglikSeq <- logLikelihoodMCTemp(
    #   indexCore = 1,
    #   list(seq[, -which(colnames(seq) == "row")]),
    #   beta, theta, initTime, endTime,
    #   splitIndicesPerCore = list(1),
    #   actDfnodes = actDfnodes, net0 = net0,
    #   formula = formula, temp
    # )
    # newloglikSeq <- logLikelihoodMCTemp(
    #   indexCore = 1,
    #   list(step$newseq[, -which(colnames(step$newseq) == "row")]),
    #   beta, theta, initTime, endTime,
    #   splitIndicesPerCore = list(1),
    #   actDfnodes = actDfnodes, net0 = net0,
    #   formula = formula, temp
    # )
  } else if (type == 2) { # Shortening
    step <- stepShort(
      seq, tieNames, gammaEminus, gammaMinus, m, me,
      Kel_g1, auxDf, pShort
    )
    getNewKelMeMatrix <- getKelMeMatrix(step$newseq, actDfnodesLab)
    newAuxDfE <- getAuxDfE(getNewKelMeMatrix$auxDf, step$sender, step$receiver)
    auxDfE <- getAuxDfE(auxDf, step$sender, step$receiver)

    newKel_g1 <- getNewKelMeMatrix$Kel_g1
    newKel_ge1 <- getNewKelMeMatrix$Kel_ge1
    newGammaEminus <- choose(newKel_ge1, 2) + newKel_g1
    newGammaMinus <- sum(newGammaEminus)
    newMe <- getNewKelMeMatrix$me
    newGammaEplus <- choose(nrow(step$newseq) - newMe + 2, 2)
    diag(newGammaEplus) <- 0
    newGammaPlus <- sum(newGammaEplus)
    newM <- nrow(step$newseq)

    newp <- (newM - newMe[step$sender, step$receiver] + 1) /
      newGammaEplus[step$sender, step$receiver]

    newindexSR <- which(step$newseq$sender == step$sender &
      step$newseq$receiver == step$receiver)


    pDoStep <- step$pDoStep
    if (step$typeS == "same") {
      pUndoStep <- getpDoAugment(
        newGammaEplus, newGammaPlus, newp,
        newM, newMe, step$sender, step$receiver,
        "same", newindexSR, pAug
      )
    } else if (step$typeS == "diff") {
      pUndoStep <- getpDoAugment(
        newGammaEplus, newGammaPlus, newp,
        newM, newMe, step$sender, step$receiver,
        "diff", newindexSR, pAug
      )
    }

    for (i in 1:k) {
      step <- stepPerm(step$newseq, tieNames, newM, newMe)

      getNewKelMeMatrix <- getKelMeMatrix(step$newseq, actDfnodesLab)
      newMe <- getNewKelMeMatrix$me
      newGammaEplus <- choose(nrow(step$newseq) - newMe + 2, 2)
      diag(newGammaEplus) <- 0
      newGammaPlus <- sum(newGammaEplus)
      newM <- nrow(step$newseq)
    }



    newlogLikelihoodStats <- getlogLikelihood(
      step$newseq, actDfnodes, net0, fixedparameters,
      beta, initTime, endTime, formula
    )

    newloglikSeq <- (newlogLikelihoodStats$resCrea$logLikelihood +
      newlogLikelihoodStats$resDel$logLikelihood) / temp
  } else if (type == 3) {
    # Only permutations performed
    step <- stepPerm(seq, tieNames, m, me)

    getNewKelMeMatrix <- getKelMeMatrix(step$newseq, actDfnodesLab)
    newMe <- getNewKelMeMatrix$me
    newGammaEplus <- choose(nrow(step$newseq) - newMe + 2, 2)
    diag(newGammaEplus) <- 0
    newGammaPlus <- sum(newGammaEplus)
    newM <- nrow(step$newseq)


    for (i in 1:k) {
      step <- stepPerm(step$newseq, tieNames, newM, newMe)

      getNewKelMeMatrix <- getKelMeMatrix(step$newseq, actDfnodesLab)
      newMe <- getNewKelMeMatrix$me
      newGammaEplus <- choose(nrow(step$newseq) - newMe + 2, 2)
      diag(newGammaEplus) <- 0
      newGammaPlus <- sum(newGammaEplus)
      newM <- nrow(step$newseq)
    }

    pDoStep <- step$pDoStep
    pUndoStep <- step$pDoStep

    newlogLikelihoodStats <- getlogLikelihood(
      step$newseq, actDfnodes, net0, fixedparameters,
      beta, initTime, endTime, formula
    )

    newloglikSeq <- (newlogLikelihoodStats$resCrea$logLikelihood +
      newlogLikelihoodStats$resDel$logLikelihood) / temp
  }


  loglikSeq <- (logLikelihoodStats$resCrea$logLikelihood +
    logLikelihoodStats$resDel$logLikelihood) / temp

  u <- runif(1, min = 0, max = 1)
  accept <- (u <= exp(newloglikSeq - loglikSeq) * pUndoStep / pDoStep)
  acceptanceDF <- data.frame("Type" = type, "Accept" = accept, "Temp" = temp,
                             "logLikStatsCrea" = logLikelihoodStats$resCrea$logLikelihood,
                             "logLikStatsDel" = logLikelihoodStats$resDel$logLikelihood,
                             "newLogLikStatsCrea" = newlogLikelihoodStats$resCrea$logLikelihood,
                             "newLogLikStatsDel" = newlogLikelihoodStats$resDel$logLikelihood)
  if (!accept) {
    newseq <- seq
    newlogLikelihoodStats <- logLikelihoodStats
    newloglikSeq <- loglikSeq

  } else {
    newseq <- step$newseq
  }

  return(list(
    "newseq" = newseq,
    "step" = step, "accept" = accept,
    "newlogLikelihoodStats" = newlogLikelihoodStats,
    "newloglikSeq" = newloglikSeq,
    "temp" = temp,
    "acceptanceDF" = acceptanceDF
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
#' @param tieNames vector, tie labs (e.g. "12" is a tie from sender 1
#' to receiver 2).
#' @param formula formula, the one used in goldfish model.
#' @param net0 matrix, initial observed network.
#' @param beta list, estimator of choice parameters of the model
#' (creation and deletion)
#' @param theta dataframe, estimator of rate parameters of the model
#' (creation and deletion) (constant)
#' @param initTime initial time
#' @param endTime end time
#' @param T0 maximum temperature
#'
#' @return list of result of stepMCMC for each sequence with different temperatures.
#'
#' @export
#'
stepPTMC <- function(indexCore, splitIndicesPerCore, seqs, H, actDfnodesLab,
                     actDfnodes, tieNames, formula, net0, beta, theta,
                     fixedparameters,
                     initTime, endTime, k, temp, nStepExch, pAug, pShort, pPerm) {
  indicesCore <- splitIndicesPerCore[[indexCore]]
  resStepPT <- vector("list", length(indicesCore))
   # browser()

  idunfixedComponentsCrea <- which(is.na(fixedparameters$Crea))
  idunfixedComponentsDel <- which(is.na(fixedparameters$Del))

  for (i in seq_along(indicesCore)) {
     # cat("Index core ",i,"\n" )

    acceptanceDF <- data.frame("Type" = c(), "Accept" = c(), "Temp" = c(),
                               "logLikStatsCrea" = c(), "logLikStatsDel" =c(),
                                "newLogLikStatsCrea" = c(),
                               "newLogLikStatsDel" = c())
    mcmcDiagDF <- data.frame()

    seq <- seqs[[indicesCore[[i]]]]
    logLikelihoodStats <- getlogLikelihood(
      seq, actDfnodes, net0,
      fixedparameters, beta, initTime,
      endTime, formula
    )
    # browser()

    for (t in 1:nStepExch) {
      if (nrow(seq) == H) {
        type <- sampleVec(c(1, 3),
          size = 1,
          prob = c(pAug / (pAug + pPerm), pPerm / (pAug + pPerm)),
          replace = TRUE
        )
      } else {
        type <- sampleVec(c(1, 2, 3),
          size = 1,
          prob = c(pAug, pShort, pPerm),
          replace = TRUE
        )
      }

      # cat("t in nStepExch ", t,"\n")
      aux <- stepPT(
        seq = seq, type = type, actDfnodesLab = actDfnodesLab,
        actDfnodes = actDfnodes, tieNames = tieNames, formula = formula,
        net0 = net0, beta = beta, fixedparameters = fixedparameters,
        theta = theta, initTime = initTime, endTime = endTime, k = k,
        temp = temp[indicesCore[[i]]], pAug = pAug, pShort = pShort,
        pPerm = pPerm,
        logLikelihoodStats = logLikelihoodStats
      )
      # cat(logLikelihoodStats$resCrea$logLikelihood," , ",logLikelihoodStats$resDel$logLikelihood,"\n")

      acceptanceDF <- rbind(acceptanceDF, aux$acceptanceDF)
      scoreVec <- aux$newlogLikelihoodStats$resCrea$finalScore[idunfixedComponentsCrea] +
        aux$newlogLikelihoodStats$resDel$finalScore[idunfixedComponentsDel]
      mcmcDiagDF <- rbind(
        mcmcDiagDF,
        c(
          "log" = aux$newloglikSeq,
          "score" = scoreVec,
          "temp" = temp[indicesCore[[i]]]
        )
      )

      seq <- aux$newseq
      logLikelihoodStats <- aux$newlogLikelihoodStats
    }
    names(mcmcDiagDF) <- c(
      "log",
      paste("score", 1:length(idunfixedComponentsCrea), sep = ""),
      "temp"
    )
    resStepPT[[i]] <- list(
      "aux" = aux, "acceptanceDF" = acceptanceDF,
      "mcmcDiagDF" = mcmcDiagDF
    )
  }

  # add acceptances for augment and shortening
  return(resStepPT)
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
#' @param thin integer, number of steps between selected sequences.
#' @param pShort float, probability of shortening step.
#' @param pAug float, probability of augmenting step.
#' @param pPerm float, probability of permutation step.
#'
#' @return permut, list of sequences.
#'
#' @export
#'
MCMC <- function(nmax, seqInit, H, actDfnodes, formula, net0, beta,
                 theta, fixedparameters, initTime, endTime, burnIn = TRUE,
                 burnInIter = 500, maxIter = 10000, thin = 50,
                 pAug = 0.35, pShort = 0.35, pPerm = 0.3, k = 5
                ) {
    # Compute initial quatities:
    # Type 1: augmentation, type 2: shortening

    if (nmax > (maxIter - burnInIter) / thin) {
      maxIter <- (nmax + burnInIter + 1) * thin
    }

    # browser()
    actDfnodesLab <- actDfnodes$label

    tieNames <- sapply(actDfnodesLab, function(x) {
      sapply(actDfnodesLab, function(i) paste(x, i, sep = "-"))
    })
    tieNames <- as.vector(tieNames)


    seqsEM <- list()

    acceptDF <- data.frame("Type" = c(), "Accept" = c(), "Temp" = c(),
                           "logLikStatsCrea" = c(), "logLikStatsDel" =c(),
                           "newLogLikStatsCrea" = c(),
                           "newLogLikStatsDel" = c())
    mcmcDiagDF <- data.frame()
    # browser()

    emSampled <- 0
    i <- 1

    while (emSampled < nmax) {
      if (!"row" %in% colnames(seqInit)) {
        seqInit <- cbind(seqInit, "row" = 1:nrow(seqInit))
      }
      # browser()
      cat("Iter i ", i , "\n")

      resstepPT <- stepPTMC(1,seqs = list(seqInit),
                            splitIndicesPerCore = list(1), H = H,
                            actDfnodesLab = actDfnodesLab,
                            actDfnodes = actDfnodes,
                            tieNames = tieNames, formula = formula,
                            net0 = net0, beta = beta,
                            theta = theta, fixedparameters = fixedparameters,
                            initTime = initTime, endTime = endTime, k = k,
                            temp = 1, nStepExch = 1, pAug = pAug,
                            pShort = pShort, pPerm = pPerm
      )

      resstepPT <- unlist(resstepPT, recursive = FALSE)
      seqInit = resstepPT$aux$newseq

      acceptDFaux <- resstepPT$acceptanceDF
      acceptDF <- rbind(acceptDF, acceptDFaux)

      if(burnIn){
        if(i> burnInIter){
          mcmcDiagDFaux <- resstepPT$mcmcDiagDF
          mcmcDiagDF <- rbind(mcmcDiagDF, mcmcDiagDFaux)
        }
      }else{
        mcmcDiagDFaux <- resstepPT$mcmcDiagDF
        mcmcDiagDF <- rbind(mcmcDiagDF, mcmcDiagDFaux)
      }


      if (burnIn) { # avoid, after burn in, to have switch and sample in the same step (not really needed)
        if (i > burnInIter & !i %% thin) {
          seqsEM <- c(seqsEM, list(resstepPT$aux)) # get a sample from the sequence with temperature 1
          emSampled = emSampled + 1
        }
      } else {
        if (!i %% thin) {
          seqsEM <- c(seqsEM, list(resstepPT$aux))
          emSampled = emSampled + 1
        }
      }
      i <- i + 1
      if (i > maxIter) break
    }

    return(list(
      "seqsEM" = seqsEM, "resstepPT" = resstepPT$aux$newseq,
      "acceptDF" = acceptDF, "mcmcDiagDF" = mcmcDiagDF
    ))
  }





#' Parallel tempering MCMC function (augmentation/shortening)
#'
#' @description Computes an MCMC chain and returns nmax sequences from the MCMC
#' chain.
#'
#' @param nmax integer, number of sequences to sampled.
#' @param nPT integer, number of chains in the parallel tempering algorithm.
#' @param seqsInit list of data frame, sequences from which start/continue MCMC
#' chain. They must have all the same length.
#' @param H integer, hamming distance between net0 and net1.
#' @param actDfnodes object of type ´nodes´ from goldfish.
#' @param formula formula, the one used in goldfish model.
#' @param net0 matrix, initial observed network.
#' @param beta list, estimator of parameters of the model (creation and deletion).
#' @param burnIn boolean, indicates if burn in must be performed.
#' @param maxIter integer, maximum number of steps of the MCMC chain.
#' @param thin integer, number of steps between selected sequences.
#' @param pShort float, probability of shortening step.
#' @param pAug float, probability of augmenting step.
#' @param T0 initial maximal tempetature for the parallale tempering scheme/
#' simulated annealing
#' @param nStepExch number of mutation steps before a exchange
#'
#' @return permut, list of sequences.
#'
#' @export
#'
PT_MCMC <- function(nmax, nPT, seqsPT, H, actDfnodes, formula, net0, beta,
                    theta, fixedparameters, initTime, endTime, burnIn = TRUE,
                    burnInIter = 500, maxIter = 10000,
                    thin = 50, T0 = 100, nStepExch = 10,
                    pAug = 0.35, pShort = 0.35, pPerm = 0.3, k = 5,
                    cl = cl, num_cores = num_cores) {
  # Compute initial quatities:
  # Type 1: augmentation, type 2: shortening

  if (nmax > (maxIter - burnInIter) / thin) {
    maxIter <- (nmax + burnInIter + 1) * thin
  }

  # browser()
  actDfnodesLab <- actDfnodes$label

  tieNames <- sapply(actDfnodesLab, function(x) {
    sapply(actDfnodesLab, function(i) paste(x, i, sep = "-"))
  })
  tieNames <- as.vector(tieNames)


  if (is.data.frame(seqsPT)) {
    seqsPT <- permute(seqsPT, nmax = nPT)
  }
  seqsEM <- list()

  temp <- seq(1, T0, length = length(seqsPT))

  acceptSwitch <- data.frame("Accept" = c(), "Temp1" = c(), "Temp2" = c())
  acceptDF <- data.frame("Type" = c(), "Accept" = c(), "Temp" = c(),
                         "logLikStatsCrea" = c(),
                         "logLikStatsDel" = c(),
                         "newLogLikStatsCrea" = c(),
                         "newLogLikStatsDel" = c())
  mcmcDiagDF <- data.frame()
  # browser()

  splitIndicesPerCore <- splitIndices(length(seqsPT), num_cores)

  emSampled <- 0
  i <- 1
  # for (i in 1:maxIter) {
  while (emSampled < nmax) {
    if (!"row" %in% colnames(seqsPT[[1]])) {
      seqsPT <- lapply(seqsPT, function(x) cbind(x, "row" = 1:nrow(x)))
    }
     # browser()
    cat("Iter i ", i , "\n")

    resstepPT <- clusterApplyLB(cl, seq_along(splitIndicesPerCore), stepPTMC,
      seqs = seqsPT, splitIndicesPerCore = splitIndicesPerCore, H = H,
      actDfnodesLab = actDfnodesLab, actDfnodes = actDfnodes,
      tieNames = tieNames, formula = formula, net0 = net0, beta = beta,
      theta = theta, fixedparameters = fixedparameters, initTime = initTime,
      endTime = endTime, k = k, temp = temp,
      nStepExch = nStepExch, pAug = pAug, pShort = pShort, pPerm = pPerm
    )

    resstepPT <- unlist(resstepPT, recursive = FALSE)

    seqsPT = unlist(lapply(lapply(resstepPT, "[[", "aux"), "[", "newseq"),
                    recursive=FALSE)

    acceptDFaux <- lapply(resstepPT, function(x) x$acceptanceDF)
    acceptDFaux <- as.data.frame(do.call(rbind, acceptDFaux))
    acceptDF <- rbind(acceptDF, acceptDFaux)

    if(burnIn){
      if(i> burnInIter){
        mcmcDiagDFaux <- lapply(resstepPT, function(x) x$mcmcDiagDF)
        mcmcDiagDFaux <- as.data.frame(do.call(rbind, mcmcDiagDFaux))
        mcmcDiagDF <- rbind(mcmcDiagDF, mcmcDiagDFaux)
      }
    }else{
        mcmcDiagDFaux <- lapply(resstepPT, function(x) x$mcmcDiagDF)
        mcmcDiagDFaux <- as.data.frame(do.call(rbind, mcmcDiagDFaux))
        mcmcDiagDF <- rbind(mcmcDiagDF, mcmcDiagDFaux)
     }

    # browser()

    # Every 10 iterations, try to switch
    order <- sample(c(0, 1), size = 1, prob = c(0.5, 0.5)) # True=1: Pairs of sequences 1-2, 3-4, ....
    # False=0: Paris of sequences 1, 2-3, 4-5, ...
    if (order) {
      if (nPT == 2) {
        tempPairs <- data.frame(c(1), c(2))
      } else {
        if (!nPT %% 2) {
          tempPairs <- data.frame(seq(1, nPT, by = 2), seq(2, nPT, by = 2))
        } else {
          tempPairs <- data.frame(seq(1, nPT - 1, by = 2), seq(2, nPT, by = 2))
        }
      }
    } else {
      if (nPT == 2) {
        tempPairs <- data.frame(c(1), c(2))
      } else {
        if (!nPT %% 2) {
          tempPairs <- data.frame(seq(2, nPT - 1, by = 2), seq(3, nPT, by = 2))
        } else {
          tempPairs <- data.frame(seq(2, nPT, by = 2), seq(3, nPT, by = 2))
        }
      }
    }


    # Exchange step!
    # Accept or reject the change
    logUnif <- log(runif(nrow(tempPairs)))
    for (j in 1:nrow(tempPairs)) {
      logMHRatio <- min(
        0,
        resstepPT[[tempPairs[j, 1]]]$aux$newloglikSeq *
          temp[tempPairs[j, 1]] / temp[tempPairs[j, 2]] +
          resstepPT[[tempPairs[j, 2]]]$aux$newloglikSeq *
            temp[tempPairs[j, 2]] / temp[tempPairs[j, 1]] -
          resstepPT[[tempPairs[j, 1]]]$aux$newloglikSeq -
          resstepPT[[tempPairs[j, 2]]]$aux$newloglikSeq
      )
      if (logUnif[j] < logMHRatio) {
        acceptSwitch <- rbind(acceptSwitch,
                              data.frame("Accept" = TRUE,
                                         "Temp1" = tempPairs[j, 1],
                                         "Temp2" = tempPairs[j, 2]))
        auxPT <- resstepPT[[tempPairs[j, 1]]]
        resstepPT[[tempPairs[j, 1]]] <- resstepPT[[tempPairs[j, 2]]]
        resstepPT[[tempPairs[j, 2]]] <- auxPT
        rm(auxPT)
      } else {
        acceptSwitch <- rbind(acceptSwitch,
                              data.frame("Accept" = FALSE,
                                         "Temp1" = tempPairs[j, 1],
                                         "Temp2" = tempPairs[j, 2]))
      }
    }

  # browser()
    if (burnIn) { # avoid, after burn in, to have switch and sample in the same step (not really needed)
      if ((i - 5) > burnInIter & !(i - 5) %% thin) {
        seqsEM <- c(seqsEM, list(resstepPT[[1]]$aux)) # get a sample from the sequence with temperature 1
        emSampled = emSampled + 1
        }
    } else {
      if (!(i - 5) %% thin) {
        seqsEM <- c(seqsEM, list(resstepPT[[1]]$aux))
        emSampled = emSampled + 1
      }
    }
    i <- i + 1
    if (i > maxIter) break
  }


  return(list(
    "seqsEM" = seqsEM, "resstepPT" = lapply(lapply(resstepPT, "[[", "aux"),
                                            "[", "newseq"),
    "acceptSwitch" = acceptSwitch, "acceptDF" = acceptDF,
    "mcmcDiagDF" = mcmcDiagDF
  ))
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


# MCMC with rates -------

#' stepRatePT (simple version with augmentation/shortening/permutation)
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
#' @param thin integer, number of steps between selected sequences.
#' @param pShort float, probability of shortening step.
#' @param pAug float, probability of augmenting step.
#' @param pPerm float, probability of permutation step.
#'
#' @return permut, list of sequences.
#'
#' @export
#'
stepRatePT <- function(seq, type, actDfnodesLab, actDfnodes, tieNames, formula,
                       net0,beta, theta, fixedparameters, initTime, endTime,
                       k = 5, temp = 1, pAug, pShort, pPerm,
                       logLikelihoodStats, loglikRate) {

   # browser()
  getKelMeMatrix <- getKelMeMatrix(seq, actDfnodesLab)
  Kel_g1 <- getKelMeMatrix$Kel_g1
  Kel_ge1 <- getKelMeMatrix$Kel_ge1
  gammaEminus <- choose(Kel_ge1, 2) + Kel_g1
  gammaMinus <- sum(gammaEminus)
  if (gammaMinus == 0) {
    type <- 1
  }
  me <- getKelMeMatrix$me
  gammaEplus <- choose(nrow(seq) - me + 2, 2)
  diag(gammaEplus) <- 0
  gammaPlus <- sum(gammaEplus)
  m <- nrow(seq)
  auxDf <- getKelMeMatrix$auxDf

  if (type == 1) { # Augmentation
    step <- stepAugment(seq, tieNames, gammaEplus, gammaPlus, m, me, net0, pAug)

    getNewKelMeMatrix <- getKelMeMatrix(step$newseq, actDfnodesLab)
    newAuxDfE <- getAuxDfE(getNewKelMeMatrix$auxDf, step$sender, step$receiver)
    auxDfE <- getAuxDfE(auxDf, step$sender, step$receiver)

    newKel_g1 <- getNewKelMeMatrix$Kel_g1
    newKel_ge1 <- getNewKelMeMatrix$Kel_ge1
    newGammaEminus <- choose(newKel_ge1, 2) + newKel_g1
    newGammaMinus <- sum(newGammaEminus)
    newMe <- getNewKelMeMatrix$me
    newGammaEplus <- choose(nrow(step$newseq) - newMe + 2, 2)
    diag(newGammaEplus) <- 0
    newGammaPlus <- sum(newGammaEplus)
    newM <- nrow(step$newseq)

    pDoStep <- step$pDoStep
    newp <- newKel_g1[step$sender, step$receiver] /
      newGammaEminus[step$sender, step$receiver]
    if ((length(unique(auxDfE$run)) == length(unique(newAuxDfE$run))) | step$typeA == "same") {
      pUndoStep <- getpDoShort(
        newGammaEminus, newGammaMinus, newp,
        length(which(table(newAuxDfE$run) > 1)), step$sender,
        step$receiver, "same", pShort
      )
    } else {
      pUndoStep <- getpDoShort(
        newGammaEminus, newGammaMinus, newp,
        length(unique(newAuxDfE$run)), step$sender,
        step$receiver, "diff", pShort
      )
    }


    if(k>0){
    for (i in 1:k) {
      step <- stepPerm(step$newseq, tieNames, newM, newMe)

      getNewKelMeMatrix <- getKelMeMatrix(step$newseq, actDfnodesLab)
      newMe <- getNewKelMeMatrix$me
      newGammaEplus <- choose(nrow(step$newseq) - newMe + 2, 2)
      diag(newGammaEplus) <- 0
      newGammaPlus <- sum(newGammaEplus)
      newM <- nrow(step$newseq)
    }
    }

    # browser()

    newlogLikelihoodStats <- getlogLikelihood(
      step$newseq, actDfnodes, net0, fixedparameters,
      beta, initTime, endTime, formula
    )

    newloglikRate <- getlogLikelihoodRate(step$newseq,actDfnodes,theta,
                                          initTime,endTime,temp)

    newloglikSeq <- (newlogLikelihoodStats$resCrea$logLikelihood +
                       newlogLikelihoodStats$resDel$logLikelihood) / temp +
                      newloglikRate$resCrea$logLikelihood +
                      newloglikRate$resDel$logLikelihood


  } else if (type == 2) { # Shortening
    step <- stepShort(
      seq, tieNames, gammaEminus, gammaMinus, m, me,
      Kel_g1, auxDf, pShort
    )
    getNewKelMeMatrix <- getKelMeMatrix(step$newseq, actDfnodesLab)
    newAuxDfE <- getAuxDfE(getNewKelMeMatrix$auxDf, step$sender, step$receiver)
    auxDfE <- getAuxDfE(auxDf, step$sender, step$receiver)

    newKel_g1 <- getNewKelMeMatrix$Kel_g1
    newKel_ge1 <- getNewKelMeMatrix$Kel_ge1
    newGammaEminus <- choose(newKel_ge1, 2) + newKel_g1
    newGammaMinus <- sum(newGammaEminus)
    newMe <- getNewKelMeMatrix$me
    newGammaEplus <- choose(nrow(step$newseq) - newMe + 2, 2)
    diag(newGammaEplus) <- 0
    newGammaPlus <- sum(newGammaEplus)
    newM <- nrow(step$newseq)

    newp <- (newM - newMe[step$sender, step$receiver] + 1) /
      newGammaEplus[step$sender, step$receiver]

    newindexSR <- which(step$newseq$sender == step$sender &
                          step$newseq$receiver == step$receiver)


    pDoStep <- step$pDoStep
    if (step$typeS == "same") {
      pUndoStep <- getpDoAugment(
        newGammaEplus, newGammaPlus, newp,
        newM, newMe, step$sender, step$receiver,
        "same", newindexSR, pAug
      )
    } else if (step$typeS == "diff") {
      pUndoStep <- getpDoAugment(
        newGammaEplus, newGammaPlus, newp,
        newM, newMe, step$sender, step$receiver,
        "diff", newindexSR, pAug
      )
    }

    if(k>0){
    for (i in 1:k) {
      step <- stepPerm(step$newseq, tieNames, newM, newMe)

      getNewKelMeMatrix <- getKelMeMatrix(step$newseq, actDfnodesLab)
      newMe <- getNewKelMeMatrix$me
      newGammaEplus <- choose(nrow(step$newseq) - newMe + 2, 2)
      diag(newGammaEplus) <- 0
      newGammaPlus <- sum(newGammaEplus)
      newM <- nrow(step$newseq)
    }
    }


    newlogLikelihoodStats <- getlogLikelihood(
      step$newseq, actDfnodes, net0, fixedparameters,
      beta, initTime, endTime, formula
    )

    newloglikRate <- getlogLikelihoodRate(step$newseq,actDfnodes,theta,
                                          initTime,endTime,temp)

    newloglikSeq <- (newlogLikelihoodStats$resCrea$logLikelihood +
                       newlogLikelihoodStats$resDel$logLikelihood) / temp +
                      newloglikRate$resCrea$logLikelihood +
                      newloglikRate$resDel$logLikelihood

  } else if (type == 3) {
    # Only permutations performed
    step <- stepPerm(seq, tieNames, m, me)

    getNewKelMeMatrix <- getKelMeMatrix(step$newseq, actDfnodesLab)
    newMe <- getNewKelMeMatrix$me
    newGammaEplus <- choose(nrow(step$newseq) - newMe + 2, 2)
    diag(newGammaEplus) <- 0
    newGammaPlus <- sum(newGammaEplus)
    newM <- nrow(step$newseq)


    if(k>0){
    for (i in 1:k) {
      step <- stepPerm(step$newseq, tieNames, newM, newMe)

      getNewKelMeMatrix <- getKelMeMatrix(step$newseq, actDfnodesLab)
      newMe <- getNewKelMeMatrix$me
      newGammaEplus <- choose(nrow(step$newseq) - newMe + 2, 2)
      diag(newGammaEplus) <- 0
      newGammaPlus <- sum(newGammaEplus)
      newM <- nrow(step$newseq)
    }
    }

    pDoStep <- step$pDoStep
    pUndoStep <- step$pDoStep

    newlogLikelihoodStats <- getlogLikelihood(
      step$newseq, actDfnodes, net0, fixedparameters,
      beta, initTime, endTime, formula
    )

    newloglikRate <- getlogLikelihoodRate(step$newseq,actDfnodes,theta,
                                          initTime,endTime,temp)

    newloglikSeq <- (newlogLikelihoodStats$resCrea$logLikelihood +
                       newlogLikelihoodStats$resDel$logLikelihood) / temp +
                      newloglikRate$resCrea$logLikelihood +
                      newloglikRate$resDel$logLikelihood
  }


  loglikSeq <- (logLikelihoodStats$resCrea$logLikelihood +
                  logLikelihoodStats$resDel$logLikelihood) / temp +
                loglikRate$resCrea$logLikelihood +
                loglikRate$resDel$logLikelihood

  u <- runif(1, min = 0, max = 1)
  # browser()
  accept <- (u <= exp(newloglikSeq - loglikSeq) * pUndoStep / pDoStep)
  acceptanceDF <- data.frame("Type" = type, "Accept" = accept, "Temp" = temp,
                             "logLikStatsCreaChoice" = logLikelihoodStats$resCrea$logLikelihood,
                             "logLikStatsDelChoice" = logLikelihoodStats$resDel$logLikelihood,
                             "newLogLikStatsCreaChoice" = newlogLikelihoodStats$resCrea$logLikelihood,
                             "newLogLikStatsDelChoice" = newlogLikelihoodStats$resDel$logLikelihood,
                             "loglikRateCrea" = loglikRate$resCrea$logLikelihood,
                             "loglikRateDel" = loglikRate$resDel$logLikelihood,
                             "newloglikRateCrea" = newloglikRate$resCrea$logLikelihood,
                             "newloglikRateDel" = newloglikRate$resDel$logLikelihood)
  if (!accept) {
    newseq <- seq
    newlogLikelihoodStats <- logLikelihoodStats
    newloglikSeq <- loglikSeq
    newloglikRate <- loglikRate

  } else {
    newseq <- step$newseq
  }

  return(list(
    "newseq" = newseq,
    "step" = step, "accept" = accept,
    "newlogLikelihoodStats" = newlogLikelihoodStats,
    "newloglikRate" = newloglikRate,
    "newloglikSeq" = newloglikSeq,
    "temp" = temp,
    "acceptanceDF" = acceptanceDF
  ))
}





#' stepRatePTMC (simple version with augmentation/shortening/permutation)
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
#' @param thin integer, number of steps between selected sequences.
#' @param pShort float, probability of shortening step.
#' @param pAug float, probability of augmenting step.
#' @param pPerm float, probability of permutation step.
#'
#' @return permut, list of sequences.
#'
#' @export
#'
stepRatePTMC <- function(indexCore, splitIndicesPerCore, seqs, H, actDfnodesLab,
                     actDfnodes, tieNames, formula, net0, beta, theta,
                     fixedparameters, initTime, endTime, k, temp, nStepExch,
                     pAug, pShort, pPerm) {

  indicesCore <- splitIndicesPerCore[[indexCore]]
  resStepPT <- vector("list", length(indicesCore))

  idunfixedComponentsCrea <- which(is.na(fixedparameters$Crea))
  idunfixedComponentsDel <- which(is.na(fixedparameters$Del))

  # browser()

  for (i in seq_along(indicesCore)) {
    # cat("Index core ",i,"\n" )

    acceptanceDF <- data.frame("Type" = c(), "Accept" = c(), "Temp" = c(),
                                "logLikStatsCreaChoice" = c(),
                                "logLikStatsDelChoice" = c(),
                                "newLogLikStatsCreaChoice" = c(),
                                "newLogLikStatsDelChoice" = c(),
                                "loglikRateCrea" = c(),
                                "loglikRateDel" = c(),
                                "newloglikRateCrea" = c(),
                                "newloglikRateDel" = c())
    mcmcDiagDF <- data.frame()

    seq <- seqs[[indicesCore[[i]]]]

    # browser()
    logLikelihoodStats <- getlogLikelihood(
      seq, actDfnodes, net0,
      fixedparameters, beta, initTime,
      endTime, formula
    )

    loglikRate = getlogLikelihoodRate(
      seq, actDfnodes, theta,
      initTime, endTime, temp[indicesCore[[i]]]
    )
    # browser()

    for (t in 1:nStepExch) {
      if (nrow(seq) == H) {
        type <- sampleVec(c(1, 3),
                          size = 1,
                          prob = c(pAug / (pAug + pPerm), pPerm / (pAug + pPerm)),
                          replace = TRUE
        )
      } else {
        type <- sampleVec(c(1, 2, 3),
                          size = 1,
                          prob = c(pAug, pShort, pPerm),
                          replace = TRUE
        )
      }

      # cat("t in nStepExch ", t,"\n")
      aux <- stepRatePT(
        seq = seq, type = type, actDfnodesLab = actDfnodesLab,
        actDfnodes = actDfnodes, tieNames = tieNames, formula = formula,
        net0 = net0, beta = beta, fixedparameters = fixedparameters,
        theta = theta, initTime = initTime, endTime = endTime, k = k,
        temp = temp[indicesCore[[i]]], pAug = pAug, pShort = pShort,
        pPerm = pPerm, logLikelihoodStats = logLikelihoodStats,
        loglikRate = loglikRate
      )

      # print(aux$acceptanceDF)
      acceptanceDF <- rbind(acceptanceDF, aux$acceptanceDF)

      logChoiceCrea = aux$newlogLikelihoodStats$resCrea$logLikelihood
      logChoiceDel =  aux$newlogLikelihoodStats$resDel$logLikelihood
      scoreVecChoiceCrea <- aux$newlogLikelihoodStats$resCrea$finalScore[idunfixedComponentsCrea]
      scoreVecChoiceDel <-  aux$newlogLikelihoodStats$resDel$finalScore[idunfixedComponentsDel]
      scoreRateCrea <- aux$newloglikRate$resCrea$score
      scoreRateDel <- aux$newloglikRate$resDel$score
      logRateCrea <- aux$newloglikRate$resCrea$logLikelihood
      logRateDel <- aux$newloglikRate$resDel$logLikelihood

      mcmcDiagDF <- rbind(
        mcmcDiagDF,
        c("totalLoglik" = aux$newloglikSeq,
          "logChoiceCrea" = logChoiceCrea/temp[indicesCore[[i]]],
          "logChoiceDel" = logChoiceDel/temp[indicesCore[[i]]],
          "logRateCrea" = logRateCrea,
          "logRateDel" = logRateDel,
          "scoreChoiceCrea" = scoreVecChoiceCrea/temp[indicesCore[[i]]],
          "scoreChoiceDel" = scoreVecChoiceDel/temp[indicesCore[[i]]],
          "scoreRateCrea" = scoreRateCrea,
          "scoreRateDel" = scoreRateDel,
          "temp" = temp[indicesCore[[i]]])
      )

      seq <- aux$newseq
      logLikelihoodStats <- aux$newlogLikelihoodStats
      loglikRate <- aux$newloglikRate
    }
    names(mcmcDiagDF) <- c("totalLoglik","logChoiceCrea","logChoiceDel",
                           "logRateCrea","logRateDel",
      paste("scoreChoiceCrea", 1:length(idunfixedComponentsCrea), sep = ""),
      paste("scoreChoiceDel", 1:length(idunfixedComponentsCrea), sep = ""),
      "scoreRateCrea","scoreRateDel",
      "temp"
    )
    resStepPT[[i]] <- list(
      "aux" = aux, "acceptanceDF" = acceptanceDF,
      "mcmcDiagDF" = mcmcDiagDF
    )
  }

  # add acceptances for augment and shortening
  return(resStepPT)
}



#' MCMC_rate (simple version with augmentation/shortening/permutation)
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
#' @param thin integer, number of steps between selected sequences.
#' @param pShort float, probability of shortening step.
#' @param pAug float, probability of augmenting step.
#' @param pPerm float, probability of permutation step.
#'
#' @return permut, list of sequences.
#'
#' @export
#'
MCMC_rate <- function(nmax, seqInit, H, actDfnodes, formula, net0, beta,
                 theta, fixedparameters, initTime, endTime, burnIn = TRUE,
                 burnInIter = 500, maxIter = 10000, thin = 50,
                 pAug = 0.35, pShort = 0.35, pPerm = 0.3, k = 5
) {
  # Compute initial quatities:
  # Type 1: augmentation, type 2: shortening

  if (nmax > (maxIter - burnInIter) / thin) {
    maxIter <- (nmax + burnInIter + 1) * thin
  }

  # browser()
  actDfnodesLab <- actDfnodes$label

  tieNames <- sapply(actDfnodesLab, function(x) {
    sapply(actDfnodesLab, function(i) paste(x, i, sep = "-"))
  })
  tieNames <- as.vector(tieNames)


  seqsEM <- list()

  acceptDF <- data.frame("Type" = c(), "Accept" = c(), "Temp" = c(),
                         "logLikStatsCreaChoice" = c(),
                         "logLikStatsDelChoice" = c(),
                         "newLogLikStatsCreaChoice" = c(),
                         "newLogLikStatsDelChoice" = c(),
                         "loglikRateCrea" = c(),
                         "loglikRateDel" = c(),
                         "newloglikRateCrea" = c(),
                         "newloglikRateDel" = c())
  mcmcDiagDF <- data.frame()
  # browser()

  emSampled <- 0
  i <- 1

  while (emSampled < nmax) {
    if (!"row" %in% colnames(seqInit)) {
      seqInit <- cbind(seqInit, "row" = 1:nrow(seqInit))
    }
    # browser()
    cat("Iter i ", i , "\n")

    resstepPT <- stepRatePTMC(1,seqs = list(seqInit),
                              splitIndicesPerCore = list(1), H = H,
                              actDfnodesLab = actDfnodesLab,
                              actDfnodes = actDfnodes, tieNames = tieNames,
                              formula = formula, net0 = net0, beta = beta,
                              theta = theta, fixedparameters = fixedparameters,
                              initTime = initTime, endTime = endTime, k = k,
                              temp = 1, nStepExch = 1, pAug = pAug,
                              pShort = pShort, pPerm = pPerm
    )

    resstepPT <- unlist(resstepPT, recursive = FALSE)
    seqInit = resstepPT$aux$newseq

    acceptDFaux <- resstepPT$acceptanceDF
    acceptDF <- rbind(acceptDF, acceptDFaux)

    if(burnIn){
      if(i> burnInIter){
        mcmcDiagDFaux <- resstepPT$mcmcDiagDF
        mcmcDiagDF <- rbind(mcmcDiagDF, mcmcDiagDFaux)
      }
    }else{
      mcmcDiagDFaux <- resstepPT$mcmcDiagDF
      mcmcDiagDF <- rbind(mcmcDiagDF, mcmcDiagDFaux)
    }


    if (burnIn) { # avoid, after burn in, to have switch and sample in the same step (not really needed)
      if (i > burnInIter & !i %% thin) {
        seqsEM <- c(seqsEM, list(resstepPT$aux)) # get a sample from the sequence with temperature 1
        emSampled = emSampled + 1
      }
    } else {
      if (!i %% thin) {
        seqsEM <- c(seqsEM, list(resstepPT$aux))
        emSampled = emSampled + 1
      }
    }
    i <- i + 1
    if (i > maxIter) break
  }

  return(list(
    "seqsEM" = seqsEM, "resstepPT" = resstepPT$aux$newseq,
    "acceptDF" = acceptDF, "mcmcDiagDF" = mcmcDiagDF
  ))
}


#' PT_Rate_MCMC (simple version with augmentation/shortening/permutation)
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
#' @param thin integer, number of steps between selected sequences.
#' @param pShort float, probability of shortening step.
#' @param pAug float, probability of augmenting step.
#' @param pPerm float, probability of permutation step.
#'
#' @return permut, list of sequences.
#'
#' @export
#'
PT_Rate_MCMC <- function(nmax, nPT, seqsPT, H, actDfnodes, formula, net0, beta,
                    theta, fixedparameters, initTime, endTime, burnIn = TRUE,
                    burnInIter = 500, maxIter = 10000,
                    thin = 50, T0 = 100, nStepExch = 10,
                    pAug = 0.35, pShort = 0.35, pPerm = 0.3, k = 5,
                    num_cores = num_cores,typeTemp = "sequential",r=1/2) {

  cl <- makeCluster(num_cores)
  on.exit(stopCluster(cl))
  clusterExport(cl, c(
    "formula", "net0", "beta", "theta", "initTime",
    "endTime", "k", "T0", "nStepExch", "pAug", "pShort",
    "H", "actDfnodes", "nAct", "fixedparameters"
  ),envir = environment())


  clusterEvalQ(cl, {
    library(goldfish)
    NULL
  })


  # clusterExport(cl, list(
  #   "stepRatePT", "stepRatePTMC", "getpDoAugment", "getpDoShort",
  #   "getpDoPerm", "getKelMeMatrix", "getAuxDfE", "stepAugment", "stepShort",
  #   "stepPerm", "getlogLikelihood", "GatherPreprocessingDF",
  #   "sampleVec", "getlogLikelihoodMC","getlogLikelihoodRate",
  #   "getlogLikelihoodRateMC"
  # ))



  if (nmax > (maxIter - burnInIter) / thin) {
    maxIter <- (nmax + burnInIter + 1) * thin
  }

  # browser()
  actDfnodesLab <- actDfnodes$label

  tieNames <- sapply(actDfnodesLab, function(x) {
    sapply(actDfnodesLab, function(i) paste(x, i, sep = "-"))
  })
  tieNames <- as.vector(tieNames)


  if (is.data.frame(seqsPT)) {
    seqsPT <- permute(seqsPT, nmax = nPT)
  }
  seqsEM <- list()

  if(typeTemp == "sequential"){
    temp <- seq(1, T0, length = nPT)
  }else if(typeTemp == "exp"){
    temp <- T0^(seq(0,1,length=nPT))
  }else if (typeTemp == "geom"){
    aux = 0
    seqAux = sort(r^(1:(nPT-1)))
    for(i in 1:(nPT-1)) aux=c(aux,sum(seqAux[1:i]))
    temp <- 1 + aux
    rm(aux)
  }


  acceptSwitch <- data.frame("Accept" = c(), "Temp1" = c(), "Temp2" = c())
  acceptDF <- data.frame("Type" = c(), "Accept" = c(), "Temp" = c(),
                         "logLikStatsCreaChoice" = c(),
                         "logLikStatsDelChoice" = c(),
                         "newLogLikStatsCreaChoice" = c(),
                         "newLogLikStatsDelChoice" = c(),
                         "loglikRateCrea" = c(),
                         "loglikRateDel" = c(),
                         "newloglikRateCrea" = c(),
                         "newloglikRateDel" = c())
  mcmcDiagDF <- data.frame()
  # browser()

  if(nStepExch>thin) nStepExch=thin

  splitIndicesPerCore <- splitIndices(length(seqsPT), num_cores)

  emSampled <- 0
  i <- 0
  # for (i in 1:maxIter) {
  while (emSampled < nmax) {
    if (!"row" %in% colnames(seqsPT[[1]])) {
      seqsPT <- lapply(seqsPT, function(x) cbind(x, "row" = 1:nrow(x)))
    }
     # browser()
    # cat("Iter i ", i , "\n")
     # resstepPTaux <- stepRatePTMC(1,seqs = seqsPT,
     #                             splitIndicesPerCore = list(1),
     #                             H = H, actDfnodesLab = actDfnodesLab,
     #                             actDfnodes = actDfnodes, tieNames = tieNames,
     #                             formula = formula, net0 = net0, beta = beta,
     #                             theta = theta, fixedparameters = fixedparameters,
     #                             initTime = initTime, endTime = endTime, k = k,
     #                             temp = temp, nStepExch = 5, pAug = pAug,
     #                             pShort = pShort, pPerm = pPerm)

    resstepPT <- clusterApplyLB(cl, seq_along(splitIndicesPerCore),
                                stepRatePTMC, seqs = seqsPT,
                                splitIndicesPerCore = splitIndicesPerCore,
                                H = H, actDfnodesLab = actDfnodesLab,
                                actDfnodes = actDfnodes, tieNames = tieNames,
                                formula = formula, net0 = net0, beta = beta,
                                theta = theta, fixedparameters = fixedparameters,
                                initTime = initTime, endTime = endTime, k = k,
                                temp = temp, nStepExch = 5, pAug = pAug,
                                pShort = pShort, pPerm = pPerm
    )
    i <- i + 5

    resstepPT <- unlist(resstepPT, recursive = FALSE)

    seqsPT = unlist(lapply(lapply(resstepPT, "[[", "aux"), "[", "newseq"),
                    recursive=FALSE)

    acceptDFaux <- lapply(resstepPT, function(x) x$acceptanceDF)
    acceptDFaux <- as.data.frame(do.call(rbind, acceptDFaux))
    acceptDF <- rbind(acceptDF, acceptDFaux)

    if(burnIn){
      if(i> burnInIter){
        mcmcDiagDFaux <- lapply(resstepPT, function(x) x$mcmcDiagDF)
        mcmcDiagDFaux <- as.data.frame(do.call(rbind, mcmcDiagDFaux))
        mcmcDiagDF <- rbind(mcmcDiagDF, mcmcDiagDFaux)
      }
    }else{
      mcmcDiagDFaux <- lapply(resstepPT, function(x) x$mcmcDiagDF)
      mcmcDiagDFaux <- as.data.frame(do.call(rbind, mcmcDiagDFaux))
      mcmcDiagDF <- rbind(mcmcDiagDF, mcmcDiagDFaux)
    }

    # browser()

    # Every 10 iterations, try to switch
    order <- sample(c(0, 1), size = 1, prob = c(0.5, 0.5)) # True=1: Pairs of sequences 1-2, 3-4, ....
    # False=0: Paris of sequences 1, 2-3, 4-5, ...
    if (order) {
      if (nPT == 2) {
        tempPairs <- data.frame(c(1), c(2))
      } else {
        if (!nPT %% 2) {
          tempPairs <- data.frame(seq(1, nPT, by = 2), seq(2, nPT, by = 2))
        } else {
          tempPairs <- data.frame(seq(1, nPT - 1, by = 2), seq(2, nPT, by = 2))
        }
      }
    } else {
      if (nPT == 2) {
        tempPairs <- data.frame(c(1), c(2))
      } else {
        if (!nPT %% 2) {
          tempPairs <- data.frame(seq(2, nPT - 1, by = 2), seq(3, nPT, by = 2))
        } else {
          tempPairs <- data.frame(seq(2, nPT, by = 2), seq(3, nPT, by = 2))
        }
      }
    }


    # Exchange step!
    # Accept or reject the change
    logUnif <- log(runif(nrow(tempPairs)))
    for (j in 1:nrow(tempPairs)) {
      logMHRatio <- min(
        0,
        resstepPT[[tempPairs[j, 1]]]$aux$newloglikSeq *
          temp[tempPairs[j, 1]] / temp[tempPairs[j, 2]] +
          resstepPT[[tempPairs[j, 2]]]$aux$newloglikSeq *
          temp[tempPairs[j, 2]] / temp[tempPairs[j, 1]] -
          resstepPT[[tempPairs[j, 1]]]$aux$newloglikSeq -
          resstepPT[[tempPairs[j, 2]]]$aux$newloglikSeq
      )
      if (logUnif[j] < logMHRatio) {
        acceptSwitch <- rbind(acceptSwitch,
                              data.frame("Accept" = TRUE,
                                         "Temp1" = tempPairs[j, 1],
                                         "Temp2" = tempPairs[j, 2]))
        auxPT <- resstepPT[[tempPairs[j, 1]]]
        resstepPT[[tempPairs[j, 1]]] <- resstepPT[[tempPairs[j, 2]]]
        resstepPT[[tempPairs[j, 2]]] <- auxPT
        rm(auxPT)
      } else {
        acceptSwitch <- rbind(acceptSwitch,
                              data.frame("Accept" = FALSE,
                                         "Temp1" = tempPairs[j, 1],
                                         "Temp2" = tempPairs[j, 2]))
      }
    }


    resstepPT <- clusterApplyLB(cl, seq_along(splitIndicesPerCore),
                                stepRatePTMC, seqs = seqsPT,
                                splitIndicesPerCore = splitIndicesPerCore,
                                H = H, actDfnodesLab = actDfnodesLab,
                                actDfnodes = actDfnodes, tieNames = tieNames,
                                formula = formula, net0 = net0, beta = beta,
                                theta = theta, fixedparameters = fixedparameters,
                                initTime = initTime, endTime = endTime, k = k,
                                temp = temp, nStepExch = (nStepExch-5), pAug = pAug,
                                pShort = pShort, pPerm = pPerm
    )
    i <- i + nStepExch-5

    resstepPT <- unlist(resstepPT, recursive = FALSE)

    seqsPT = unlist(lapply(lapply(resstepPT, "[[", "aux"), "[", "newseq"),
                    recursive=FALSE)

    acceptDFaux <- lapply(resstepPT, function(x) x$acceptanceDF)
    acceptDFaux <- as.data.frame(do.call(rbind, acceptDFaux))
    acceptDF <- rbind(acceptDF, acceptDFaux)

    if(burnIn){
      if(i> burnInIter){
        mcmcDiagDFaux <- lapply(resstepPT, function(x) x$mcmcDiagDF)
        mcmcDiagDFaux <- as.data.frame(do.call(rbind, mcmcDiagDFaux))
        mcmcDiagDF <- rbind(mcmcDiagDF, mcmcDiagDFaux)
      }
    }else{
      mcmcDiagDFaux <- lapply(resstepPT, function(x) x$mcmcDiagDF)
      mcmcDiagDFaux <- as.data.frame(do.call(rbind, mcmcDiagDFaux))
      mcmcDiagDF <- rbind(mcmcDiagDF, mcmcDiagDFaux)
    }



    # browser()
    if (burnIn) { # avoid, after burn in, to have switch and sample in the same step (not really needed)
      if (i > burnInIter & !i  %% thin) {
        seqsEM <- c(seqsEM, list(resstepPT[[1]]$aux)) # get a sample from the sequence with temperature 1
        emSampled = emSampled + 1
      }
    } else {
      if (!i  %% thin) {
        seqsEM <- c(seqsEM, list(resstepPT[[1]]$aux))
        emSampled = emSampled + 1
      }
    }

    if (i > maxIter) break
  }


  return(list(
    "seqsEM" = seqsEM,
    "resstepPT" = lapply(lapply(resstepPT, "[[", "aux"), "[", "newseq"),
    "acceptSwitch" = acceptSwitch, "acceptDF" = acceptDF,
    "mcmcDiagDF" = mcmcDiagDF
  ))
}




