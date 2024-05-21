test_that("Step Augment works", {
  net0 <- s501
  net1 <- s502
  colnames(net0) <- 1:ncol(net0)
  rownames(net0) <- 1:nrow(net0)
  colnames(net1) <- 1:ncol(net1)
  rownames(net1) <- 1:nrow(net1)

  actDf <- dimnames(net0)[[2]]
  actDfnodes <- defineNodes(data.frame(label = colnames(net0)))

  nAct <- length(actDf)
  endTime <- 50
  nSimCrea <- 59 # from s501 and s502 data, number of ties created between t_0 and t_1
  nSimDel <- 56 # from s501 and s502 data, number of ties deleted between t_0 and t_1
  parmsCrea <- c(log(nSimCrea / 20 / nAct)) # indegree and outdegree of ego
  parmsDel <- c(log(nSimDel / 20 / (nAct - 4))) # indegree and outdegree of ego

  theta <- data.frame("Crea" = parmsCrea, "Del" = parmsDel)
  beta <- data.frame("Crea" = c(1, 1, 0, 0), "Del" = c(-1, -1, 0, 0))
  row.names(beta) <- c("indeg", "outdeg", "recip", "trans")
  formula <- "indeg + outdeg + recip + trans +inertia + tie(net0)"
  sequence <- EMPreprocessing(net0, net1)
  seq <- sequence
  seq$row <- 1:nrow(seq)
  actDfnodesLab <- actDfnodes$label
  tieNames <- sapply(actDfnodesLab, function(x) sapply(actDfnodesLab, function(i) paste(x, i, sep = "-")))


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
  pAug <- 0.5

  vec.lengths <- c()
  for (i in 1:100) {
    set.seed(i)
    step <- stepAugment(seq, tieNames, gammaEplus, gammaPlus, m, me, net0, pAug)
    vec.lengths <- c(vec.lengths, nrow(step$newseq))
    # vec.lengths <- c(vec.lengths, anyNA(step$newseq))
  }

  expect_equal(vec.lengths, rep(nrow(seq) + 2, 100))
})



test_that("Step Shortening works", {
  net0 <- s501
  net1 <- s502
  colnames(net0) <- 1:ncol(net0)
  rownames(net0) <- 1:nrow(net0)
  colnames(net1) <- 1:ncol(net1)
  rownames(net1) <- 1:nrow(net1)

  actDf <- dimnames(net0)[[2]]
  actDfnodes <- defineNodes(data.frame(label = colnames(net0)))

  nAct <- length(actDf)
  endTime <- 50
  nSimCrea <- 59 # from s501 and s502 data, number of ties created between t_0 and t_1
  nSimDel <- 56 # from s501 and s502 data, number of ties deleted between t_0 and t_1
  parmsCrea <- c(log(nSimCrea / 20 / nAct)) # indegree and outdegree of ego
  parmsDel <- c(log(nSimDel / 20 / (nAct - 4))) # indegree and outdegree of ego

  theta <- data.frame("Crea" = parmsCrea, "Del" = parmsDel)
  beta <- data.frame("Crea" = c(1, 1, 0, 0), "Del" = c(-1, -1, 0, 0))
  row.names(beta) <- c("indeg", "outdeg", "recip", "trans")
  formula <- "indeg + outdeg + recip + trans +inertia + tie(net0)"
  sequence <- EMPreprocessing(net0, net1)
  seq <- rbind(sequence, sequence[1:20, ])
  seq$row <- 1:nrow(seq)
  actDfnodesLab <- actDfnodes$label
  tieNames <- sapply(actDfnodesLab, function(x) sapply(actDfnodesLab, function(i) paste(x, i, sep = "-")))


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
  pShort <- 0.5

  vec.lengths <- c()
  for (i in 1:100) {
    set.seed(i)
    step <- stepShort(
      seq, tieNames, gammaEminus, gammaMinus, m, me, Kel_g1,
      auxDf, pShort
    )
    vec.lengths <- c(vec.lengths, nrow(step$newseq))
  }

  expect_equal(vec.lengths, rep(nrow(seq) - 2, 100))
})



test_that("Step Permutation works", {
  net0 <- s501
  net1 <- s502
  colnames(net0) <- 1:ncol(net0)
  rownames(net0) <- 1:nrow(net0)
  colnames(net1) <- 1:ncol(net1)
  rownames(net1) <- 1:nrow(net1)

  actDf <- dimnames(net0)[[2]]
  actDfnodes <- defineNodes(data.frame(label = colnames(net0)))

  nAct <- length(actDf)
  endTime <- 50
  nSimCrea <- 59 # from s501 and s502 data, number of ties created between t_0 and t_1
  nSimDel <- 56 # from s501 and s502 data, number of ties deleted between t_0 and t_1
  parmsCrea <- c(log(nSimCrea / 20 / nAct)) # indegree and outdegree of ego
  parmsDel <- c(log(nSimDel / 20 / (nAct - 4))) # indegree and outdegree of ego

  theta <- data.frame("Crea" = parmsCrea, "Del" = parmsDel)
  beta <- data.frame("Crea" = c(1, 1, 0, 0), "Del" = c(-1, -1, 0, 0))
  row.names(beta) <- c("indeg", "outdeg", "recip", "trans")
  formula <- "indeg + outdeg + recip + trans +inertia + tie(net0)"
  sequence <- EMPreprocessing(net0, net1)
  seq <- sequence
  seq$row <- 1:nrow(seq)
  actDfnodesLab <- actDfnodes$label
  tieNames <- sapply(actDfnodesLab, function(x) sapply(actDfnodesLab, function(i) paste(x, i, sep = "-")))


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

  vec.lengths <- c()
  for (i in 101:200) {
    set.seed(i)
    step <- stepPerm(seq, tieNames, m, me)
    vec.lengths <- c(vec.lengths, nrow(step$newseq))
  }

  expect_equal(vec.lengths, rep(nrow(seq), 100))
})
