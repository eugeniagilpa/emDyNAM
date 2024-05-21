# test_that("Step Augment works", {
#   net0 <- s501
#   net1 <- s502
#   colnames(net0) <- 1:ncol(net0)
#   rownames(net0) <- 1:nrow(net0)
#   colnames(net1) <- 1:ncol(net1)
#   rownames(net1) <- 1:nrow(net1)
#
#   actDf <- dimnames(net0)[[2]]
#   actDfnodes <- defineNodes(data.frame(label = colnames(net0)))
#
#   nAct <- length(actDf)
#   endTime <- 50
#   nSimCrea <- 59 # from s501 and s502 data, number of ties created between t_0 and t_1
#   nSimDel <- 56 # from s501 and s502 data, number of ties deleted between t_0 and t_1
#   parmsCrea <- c(log(nSimCrea / 20 / nAct)) # indegree and outdegree of ego
#   parmsDel <- c(log(nSimDel / 20 / (nAct - 4))) # indegree and outdegree of ego
#
#   theta <- data.frame("Crea" = parmsCrea, "Del" = parmsDel)
#   beta <- data.frame("Crea" = c(1, 1, 0, 0), "Del" = c(-1, -1, 0, 0))
#   row.names(beta) <- c("indeg", "outdeg", "recip", "trans")
#   formula <- "indeg + outdeg + recip + trans +inertia + tie(net0)"
#   sequence <- EMPreprocessing(net0, net1)
#   seq <- sequence
#   seq$row <- 1:nrow(seq)
#   actDfnodesLab <- actDfnodes$label
#   tieNames <- sapply(actDfnodesLab, function(x) sapply(actDfnodesLab, function(i) paste(x, i, sep = "-")))
#
#
#   for (i in 1:100) {
#     set.seed(i)
#
#     for (i in seq_along(indicesCore)) {
#       seq <- seqs[[indicesCore[[i]]]]
#       for (t in 1:nStepSwitch) {
#         if (nrow(seq) == H) {
#           type <- sampleVec(c(1, 2), size = 1, prob = c(1, 0), replace = TRUE)
#         } else {
#           type <- sampleVec(c(1, 2), size = 1, prob = c(pAug, pShort), replace = TRUE)
#         }
#         aux <- stepPT(
#           seq, type, actDfnodesLab, tieNames, formula, net0, beta, theta, initTime, endTime, k, temp[indicesCore[[i]]],
#           pAug, pShort
#         )
#         seq <- aux$newseq
#       }
#       resStepPT[[i]] <- aux
#     }
#   }
#
#   expect_equal(vec.lengths, rep(nrow(seq) + 2, 100))
# })
