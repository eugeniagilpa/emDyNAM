# rm(list=ls())
# # packages ----------------------------------------------------------------
# library(goldfish)
# library(broom)
# library(pixiedust)
# library(survival)
# library(mlogit)
# library(RSiena)
# library(matrixStats)
# library(stats)
# library(foreach)
# library(doParallel)
# library(parallel)


# Permutations -------------------------------------

#' Permutation of a sequence of events
#'
#' @param x original data frame with a column of senders, a column of receivers and a column of replace (replace = 0 means deletion of tie, replace = 1 means creation of tie)
#' @param nmax maximum number of permutations performed. The number of permutations computed will be the min{length(x)!, nmax}
#'
#' @return list of dataframes. Each dataframe is a permutation of the rows of the original list x
#' @export
#'
#' @examples X0 = s501, X1 = s502
#' seq = EMPreprocessing(X0,X1)
#' permute(seq,nmax = 5)
permute = function(x,nmax = 1000){
  n = nrow(x)
  if(factorial(n)<nmax){
    out = lapply(1:factorial(n),function(t) x[sample(n,n),])
  }else{
    out = lapply(1:nmax,function(t) x[sample(n,n),])
  }
  return(out)
}




# EM algorithm --------------------
hardEMAlgorithm = function(net0,net1,theta0,beta0,formula,nmax = 5){
  # net0 : network matrix in initial state
  # net1: network matrix in final state
  # theta0 : initial parameters for sender (list of creation and deletion with
  #          labels Crea and Del)
  # beta0: initial parameters for receiver (list of creation and deletion with
  #        labels Crea and Del)
  # formula: string with formula for the model (only rhs of formula)

  logLik = c()
  dfmlogit <- NULL

  actDfnodes <- defineNodes(data.frame(label=colnames(net0)))
  # Initialization of parameters
  theta = theta0
  beta = beta0

  # Creation of sequence of events from initial data
  sequence = EMPreprocessing(net0,net1)

  # Creation of permutations
  permut = permute(sequence, nmax = nmax)

  for(seq in permut){
      envirPrepro <- new.env()
      seq$time <- seq(1,nrow(seq))
      envirPrepro$seqTime <- seq
      envirPrepro$actDfnodes <- actDfnodes
      envirPrepro$net0 <- net0
      #seqTime = data.frame("time"=timeGenerator(seq,nAct,theta),seq)
      local(
        {
          netEvents = defineNetwork(nodes = actDfnodes,matrix=net0) |>
            linkEvents(changeEvents = seqTime, nodes = actDfnodes)

          depEvents = defineDependentEvents(
            seqTime, nodes = actDfnodes, defaultNetwork = netEvents
          )
        },
        envirPrepro
        )

      listExpandedDF = GatherPreprocessingDF(formula, envir = envirPrepro)
      logLik = c(logLik,logLikelihood(listExpandedDF,beta))
    }

      # For the maximum value, we perform goldfish to estimate the parameters and update beta
      index.max = sample(which(logLik==max(logLik)),1)
      seqMax = permut[[index.max]]
      seqMax$time = seq(1,nrow(seqMax))

      envirMaxCrea <- new.env()
      envirMaxDel <- new.env()
      envirMaxCrea$seqMax <- seqMax
      envirMaxCrea$actDfnodes <- actDfnodes
      envirMaxCrea$net0 <- net0
      local(
        {
          netEvents = defineNetwork(nodes = actDfnodes,matrix=net0) |>
            linkEvents(changeEvents = seqMax, nodes = actDfnodes)

          depEventsCrea <- defineDependentEvents(
            seqMax[seqMax$replace == 1,],
            nodes = actDfnodes,
            defaultNetwork = netEvents
          )
        },
        envirMaxCrea
      )

      formCrea = paste("depEventsCrea ~",formula,sep="")
      modCrea <- estimate(
        as.formula(formCrea),
        model = "DyNAM", subModel = "choice",
        estimationInit = list(
          startTime = 0,
          fixedParameters = c(NA, NA, NA, -20, -20),
          returnIntervalLogL = TRUE,
          returnEventProbabilities = TRUE
        ),
        progress = FALSE,
        envir = envirMaxCrea
      )
      beta$Crea = coef(modCrea)


      envirMaxDel$seqMax <- seqMax
      envirMaxDel$actDfnodes <- actDfnodes
      envirMaxDel$net0 <- net0
      local(
        {
          netEvents = defineNetwork(nodes = actDfnodes,matrix=net0) |>
            linkEvents(changeEvents = seqMax, nodes = actDfnodes)

          depEventsDel <- defineDependentEvents(
            seqMax[seqMax$replace == 0,],
            nodes = actDfnodes,
            defaultNetwork = netEvents
          )
        },
        envirMaxDel
      )

      formDel = paste("depEventsDel ~",formula,sep="")
      modDel <- estimate(
        as.formula(formDel),
        model = "DyNAM", subModel = "choice",
        estimationInit = list(
          startTime = 0,
          fixedParameters = c(NA, NA, NA, 20, 20),
          returnIntervalLogL = TRUE,
          returnEventProbabilities = TRUE
        ),
        progress = FALSE,
        envir = envirMaxDel
      )
      beta$Del = coef(modDel)


  return(list("logLik" = modCrea$logLikelihood+modDel$logLikelihood,"beta" = beta))
}
softEMAlgorithm = function(net0,net1,theta0,beta0,formula,nmax = 10){
  # net0 : network matrix in initial state
  # net1: network matrix in final state
  # theta0 : initial parameters for sender (list of creation and deletion with
  #          labels Crea and Del)
  # beta0: initial parameters for receiver (list of creation and deletion with
  #        labels Crea and Del)
  # formula: string with formula for the model (only rhs of formula)

  actDfnodes <- defineNodes(data.frame(label=colnames(net0)))


  # Initialization of parameters
  theta = theta0
  beta = beta0
  se=data.frame(Crea=rep(0,nrow(beta)),Del=rep(0,nrow(beta)))
  row.names(se) = row.names(beta)
  # Creation of sequence of events from initial data
  sequence = EMPreprocessing(net0,net1)

  # Creation of permutations
  permut = permute(sequence, nmax = nmax)

  betaCreaDF = c()
  betaDelDF = c()
  seCreaDF = c()
  seDelDF = c()


  for(seq in permut){
    # seq=permut[[i]]
    # if(!i%%100) cat("sequence ",i,"\n")
    envirPrepro <- new.env()
    envirPreproCrea <- new.env()
    envirPreproDel <- new.env()

    seq$time <- seq(1,nrow(seq))
    envirPrepro$seqTime <- seq
    envirPrepro$actDfnodes <- actDfnodes
    envirPrepro$net0 <- net0

    envirPreproCrea$seqTime <- seq
    envirPreproCrea$actDfnodes <- actDfnodes
    envirPreproCrea$net0 <- net0

    envirPreproDel$seqTime <- seq
    envirPreproDel$actDfnodes <- actDfnodes
    envirPreproDel$net0 <- net0

    local(
      {
        netEvents = defineNetwork(nodes = actDfnodes,matrix=net0) |>
          linkEvents(changeEvents = seqTime, nodes = actDfnodes)

        depEvents = defineDependentEvents(
          seqTime, nodes = actDfnodes, defaultNetwork = netEvents
        )
      },
      envirPrepro
    )

    listExpandedDF = GatherPreprocessingDF(formula, envir = envirPrepro)

    local(
      {
        netEvents = defineNetwork(nodes = actDfnodes,matrix=net0) |>
          linkEvents(changeEvents = seqTime, nodes = actDfnodes)

        depEventsCrea <- defineDependentEvents(
          seqTime[seqTime$replace == 1,],
          nodes = actDfnodes,
          defaultNetwork = netEvents,
        )
      },
      envirPreproCrea
    )
    local(
      {
        netEvents = defineNetwork(nodes = actDfnodes,matrix=net0) |>
          linkEvents(changeEvents = seqTime, nodes = actDfnodes)

        depEventsDel <- defineDependentEvents(
          seqTime[seqTime$replace == 0,],
          nodes = actDfnodes,
          defaultNetwork = netEvents,
        )
      },
      envirPreproDel
    )

    formCrea = paste("depEventsCrea ~",formula,sep="")
    modCrea <- estimate(
      as.formula(formCrea),
      model = "DyNAM", subModel = "choice",
      estimationInit = list(
        startTime = 0,
        fixedParameters = c(NA, NA, NA, -20,-20),
        returnIntervalLogL = TRUE,
        returnEventProbabilities = TRUE
      ),
      progress = FALSE,
      envir = envirPreproCrea
    )
    betaCreaDF = rbind(betaCreaDF,coef(modCrea))
    seCreaDF = rbind(seCreaDF, modCrea$standardErrors[modCrea$standardErrors>0])


    formDel = paste("depEventsDel ~",formula,sep="")
    modDel <- estimate(
      as.formula(formDel),
      model = "DyNAM", subModel = "choice",
      estimationInit = list(
        startTime = 0,
        fixedParameters = c(NA, NA, NA, 20, 20),
        returnIntervalLogL = TRUE,
        returnEventProbabilities = TRUE
      ),
      progress = FALSE,
      envir = envirPreproDel
    )
    betaDelDF = rbind(betaDelDF,coef(modDel))
    seDelDF = rbind(seDelDF, modDel$standardErrors[modDel$standardErrors>0])
  }


  diff = 1000
  index = 0
  while (diff>1e-4){ # TO DO: change this to nmax or to condition with while.
    logLik = c()
    index = index +1
    cat("Index: ",index,"\n")
    cat("Diff: ", diff,"\n")
    if (index > 1000){
      cat("No convergence\n")
      break
    }

    for(seq in permut){
      # if(!i%%100) cat("Sequence: ", i,"\n")
      envirPrepro <- new.env()
      seq$time <- seq(1,nrow(seq))
      envirPrepro$seqTime <- seq
      envirPrepro$actDfnodes <- actDfnodes
      envirPrepro$net0 <- net0

      local(
        {
          netEvents = defineNetwork(nodes = actDfnodes,matrix=net0) |>
            linkEvents(changeEvents = seqTime, nodes = actDfnodes)

          depEvents = defineDependentEvents(
            seqTime, nodes = actDfnodes, defaultNetwork = netEvents
          )
        },
        envirPrepro
      )
      listExpandedDF = GatherPreprocessingDF(formula, envir = envirPrepro)
      logLik = c(logLik,logLikelihood(listExpandedDF,beta))
    }
    # weights = exp(logLik)/sum(exp(logLik))
    c = max(logLik)
    logsumExp = c + log(sum(exp(logLik-c)))
    weights = exp(logLik - logsumExp)

    betaCreaAux = rubinsRule(betaCreaDF,seCreaDF,weights)
    betaDelAux = rubinsRule(betaDelDF,seDelDF,weights)
    diff= sqrt(norm(as.matrix(beta$Crea-betaCreaAux$mean))^2+norm(as.matrix(beta$Del-betaDelAux$mean))^2)

    beta$Crea = betaCreaAux$mean
    beta$Del = betaDelAux$mean
    se$Crea = betaCreaAux$se
    se$Del = betaDelAux$se
  }

  return(list("logLik" = logLik,"beta" = beta,"se" = se,"betaCrea" = betaCreaDF,"betaDel" = betaDelDF,"index" = index,
              "diff" = diff))
}
softEMAlgorithmMC = function(nmax,net0,net1,theta0,beta0,formula,num_cores=1){
  # net0 : network matrix in initial state
  # net1: network matrix in final state
  # theta0 : initial parameters for sender (list of creation and deletion with
  #          labels Crea and Del)
  # beta0: initial parameters for receiver (list of creation and deletion with
  #        labels Crea and Del)
  # formula: string with formula for the model (only rhs of formula)

  actDfnodes <- defineNodes(data.frame(label=colnames(net0)))


  # Initialization of parameters
  theta = theta0
  beta = beta0
  se=data.frame(Crea=rep(0,nrow(beta)),Del=rep(0,nrow(beta)))
  row.names(se) = row.names(beta)
  # Creation of sequence of events from initial data
  sequence = EMPreprocessing(net0,net1)

  # Creation of permutations
  permut = permute(sequence, nmax = nmax)
  # permut[[nmax+1]]=events[,2:4]

  splitIndicesPerCore = splitIndices(length(permut),num_cores)
  cl =  makeCluster(num_cores)
  on.exit(stopCluster(cl))
  clusterExport(cl, c("net0","actDfnodes","nAct"))
  clusterEvalQ(cl, {library(goldfish)
                    library(matrixStats)
                    NULL})
  clusterExport(cl, list("softEMAlgorithm","EMPreprocessing", "GatherPreprocessingDF",
                         "hardEMAlgorithm","logLikelihood","permute","rubinsRule",
                         "parameters","timeGenerator","logLikelihood"))
  resPar = clusterApply(cl,seq_along(splitIndicesPerCore),parametersMC,permut = permut,
                       splitIndicesPerCore=splitIndicesPerCore,actDfnodes.=actDfnodes,
                       net0.=net0,formula.=formula)

  resPar = unlist(resPar,recursive = FALSE)
  betaCreaDF = t(sapply(resPar,"[[",1))
  betaDelDF = t(sapply(resPar,"[[",2))
  seCreaDF = t(sapply(resPar,"[[",3))
  seDelDF = t(sapply(resPar,"[[",4))

  diff = 1000
  index = 0
  while (diff>1e-3){
    # browser()
    index = index +1
    cat("Index: ",index,"\n")
    cat("Diff: ", diff,"\n")
    if (index > 1000){
      cat("No convergence\n")
      break
    }
    logLik = clusterApply(cl,seq_along(splitIndicesPerCore),logLikelihoodMC,permut = permut,
                          splitIndicesPerCore=splitIndicesPerCore,
                          beta=beta,actDfnodes=actDfnodes,net0=net0,formula=formula)

    logLik = unlist(logLik)

    # weights = exp(logLik)/sum(exp(logLik))
    c = max(logLik)
    logsumExp = c + log(sum(exp(logLik-c)))
    weights = exp(logLik - logsumExp)

    betaCreaAux = rubinsRule(betaCreaDF,seCreaDF,weights)
    betaDelAux = rubinsRule(betaDelDF,seDelDF,weights)
    diff= sqrt(norm(as.matrix(beta$Crea-betaCreaAux$mean))^2+norm(as.matrix(beta$Del-betaDelAux$mean))^2)

    beta$Crea = betaCreaAux$mean
    beta$Del = betaDelAux$mean
    se$Crea = betaCreaAux$se
    se$Del = betaDelAux$se
    }

  aux_permut = permut[[1]]
  og_permut = aux_permut[order(as.numeric(row.names(aux_permut))),]
  aux_names = lapply(lapply(permut,row.names),as.numeric)

  return(list("logLik" = logLik,"beta" = beta,"se" = se,"index" = index,
              "diff" = diff,"og_permut" =og_permut, "permut"=aux_names,"betaCreaDF" = betaCreaDF,"betaDelDF"=betaDelDF))
}
# s = softEMAlgorithmMC(nmax=10,net0,net1,theta0,beta0,formula,num_cores=1)
softEMAlgorithmForEach = function(net0,net1,theta0,beta0,formula,nmax = 10){
  # net0 : network matrix in initial state
  # net1: network matrix in final state
  # theta0 : initial parameters for sender (list of creation and deletion with
  #          labels Crea and Del)
  # beta0: initial parameters for receiver (list of creation and deletion with
  #        labels Crea and Del)
  # formula: string with formula for the model (only rhs of formula)

  actDfnodes <- defineNodes(data.frame(label=colnames(net0)))


  # Initialization of parameters
  theta = theta0
  beta = beta0
  se=data.frame(Crea=rep(0,nrow(beta)),Del=rep(0,nrow(beta)))
  row.names(se) = row.names(beta)
  # Creation of sequence of events from initial data
  sequence = EMPreprocessing(net0,net1)

  # Creation of permutations
  permut = permute(sequence, nmax = nmax)

  betaCreaDF = c()
  betaDelDF = c()
  seCreaDF = c()
  seDelDF = c()


  foreach(i=1:length(permut))%do%{
    seq = permut[[i]]
    # if(!i%%100) cat("sequence ",i,"\n")
    envirPrepro <- new.env()
    envirPreproCrea <- new.env()
    envirPreproDel <- new.env()

    seq$time <- seq(1,nrow(seq))
    envirPrepro$seqTime <- seq
    envirPrepro$actDfnodes <- actDfnodes
    envirPrepro$net0 <- net0

    envirPreproCrea$seqTime <- seq
    envirPreproCrea$actDfnodes <- actDfnodes
    envirPreproCrea$net0 <- net0

    envirPreproDel$seqTime <- seq
    envirPreproDel$actDfnodes <- actDfnodes
    envirPreproDel$net0 <- net0

    local(
      {
        netEvents = defineNetwork(nodes = actDfnodes,matrix=net0) |>
          linkEvents(changeEvents = seqTime, nodes = actDfnodes)

        depEvents = defineDependentEvents(
          seqTime, nodes = actDfnodes, defaultNetwork = netEvents
        )
      },
      envirPrepro
    )

    listExpandedDF = GatherPreprocessingDF(formula, envir = envirPrepro)

    local(
      {
        netEvents = defineNetwork(nodes = actDfnodes,matrix=net0) |>
          linkEvents(changeEvents = seqTime, nodes = actDfnodes)

        depEventsCrea <- defineDependentEvents(
          seqTime[seqTime$replace == 1,],
          nodes = actDfnodes,
          defaultNetwork = netEvents,
        )
      },
      envirPreproCrea
    )
    local(
      {
        netEvents = defineNetwork(nodes = actDfnodes,matrix=net0) |>
          linkEvents(changeEvents = seqTime, nodes = actDfnodes)

        depEventsDel <- defineDependentEvents(
          seqTime[seqTime$replace == 0,],
          nodes = actDfnodes,
          defaultNetwork = netEvents,
        )
      },
      envirPreproDel
    )

    formCrea = paste("depEventsCrea ~",formula,sep="")
    modCrea <- estimate(
      as.formula(formCrea),
      model = "DyNAM", subModel = "choice",
      estimationInit = list(
        startTime = 0,
        fixedParameters = c(NA, NA, NA, -20,-20),
        returnIntervalLogL = TRUE,
        returnEventProbabilities = TRUE
      ),
      progress = FALSE,
      envir = envirPreproCrea
    )
    betaCreaDF = rbind(betaCreaDF,coef(modCrea))
    seCreaDF = rbind(seCreaDF, modCrea$standardErrors[modCrea$standardErrors>0])


    formDel = paste("depEventsDel ~",formula,sep="")
    modDel <- estimate(
      as.formula(formDel),
      model = "DyNAM", subModel = "choice",
      estimationInit = list(
        startTime = 0,
        fixedParameters = c(NA, NA, NA, 20, 20),
        returnIntervalLogL = TRUE,
        returnEventProbabilities = TRUE
      ),
      progress = FALSE,
      envir = envirPreproDel
    )
    betaDelDF = rbind(betaDelDF,coef(modDel))
    seDelDF = rbind(seDelDF, modDel$standardErrors[modDel$standardErrors>0])
  }


  diff = 1000
  index = 0
  while (diff>1e-4){ # TO DO: change this to nmax or to condition with while.
    logLik = c()
    index = index +1
    cat("Index: ",index,"\n")
    cat("Diff: ", diff,"\n")
    if (index > 1000){
      cat("No convergence\n")
      break
    }

    foreach(i=1:length(permut))%do%{
      seq = permut[[i]]
      # if(!i%%100) cat("Sequence: ", i,"\n")
      envirPrepro <- new.env()
      seq$time <- seq(1,nrow(seq))
      envirPrepro$seqTime <- seq
      envirPrepro$actDfnodes <- actDfnodes
      envirPrepro$net0 <- net0

      local(
        {
          netEvents = defineNetwork(nodes = actDfnodes,matrix=net0) |>
            linkEvents(changeEvents = seqTime, nodes = actDfnodes)

          depEvents = defineDependentEvents(
            seqTime, nodes = actDfnodes, defaultNetwork = netEvents
          )
        },
        envirPrepro
      )
      listExpandedDF = GatherPreprocessingDF(formula, envir = envirPrepro)
      logLik = c(logLik,logLikelihood(listExpandedDF,beta))
    }

    c = max(logLik)
    logsumExp = c + log(sum(exp(logLik-c)))
    weights = exp(logLik - logsumExp)
    # weights = exp(logLik)/sum(exp(logLik))


    betaCreaAux = rubinsRule(betaCreaDF,seCreaDF,weights)
    betaDelAux = rubinsRule(betaDelDF,seDelDF,weights)
    diff= sqrt(norm(as.matrix(beta$Crea-betaCreaAux$mean))^2+norm(as.matrix(beta$Del-betaDelAux$mean))^2)

    beta$Crea = betaCreaAux$mean
    beta$Del = betaDelAux$mean
    se$Crea = betaCreaAux$se
    se$Del = betaDelAux$se
  }

  aux_permut = permut[[1]]
  og_permut = aux_permut[order(as.numeric(row.names(aux_permut))),]
  aux_names = lapply(lapply(permut,row.names),as.numeric)

  return(list("logLik" = logLik,"beta" = beta,"se" = se,"betaCrea" = betaCreaDF,"betaDel" = betaDelDF,"index" = index,
              "diff" = diff,"og_permut" =og_permut, "permut"=aux_names,"betaCreaDF" = betaCreaDF,"betaDelDF"=betaDelDF))
}





# Computations --------
# load("simulationData.RData")
# load("/Users/mariaeugeniagilpallares/MEGA/master/HS23/Semester Paper/emDyNAM/simulation_p_01_for_em.RData")
# load("/cluster/scratch/egilpallares/simulation_p_02_for_em_time160.RData")
#
# seqIndex = ceiling(10^seq(1.6,3.4,0.2))
# results = vector("list",length = length(simulation))
#
# chunk <- strtoi(Sys.getenv(c("SLURM_ARRAY_TASK_ID")))

# for(i in 1:length(simulation)){
#   sim = simulation[[i]]
#   net0 = sim$net0
#   net1 = sim$net1
#
#   actDf <-dimnames(net0)[[2]]
#   actDfnodes <- defineNodes(data.frame(label=colnames(net0)))
#
#   nAct <- length(actDf)
#   endTime <- 50
#   nSimCrea <- 59 # from s501 and s502 data, number of ties created between t_0 and t_1
#   nSimDel <- 56 # from s501 and s502 data, number of ties deleted between t_0 and t_1
#   parmsCrea <- c(log(nSimCrea / 20 / nAct )) #indegree and outdegree of ego
#   parmsDel <- c(log(nSimDel / 20 / (nAct-4))) #indegree and outdegree of ego
#
#   theta0=data.frame("Crea"=parmsCrea,"Del"=parmsDel)
#   beta0=data.frame("Crea"=c(0,0,0,0),"Del"=c(0,0,0,0))
#   row.names(beta0) = c("indeg","outdeg","recip","trans")
#   formula = "indeg + outdeg + recip + trans +inertia + tie(net0)"
#
#
#   results[[i]] = softEMAlgorithmMC(nmax = seqIndex[chunk],net0=net0,net1=net1,
#                         theta0=theta0,beta0=beta0,formula = formula,num_cores=50)
#   # results[[i]] = softEMAlgorithmMC(nmax = 2512,net0=net0,net1=net1,
#   #                       theta0=theta0,beta0=beta0,formula = formula,num_cores=50)
#
#
# }

# filename = paste("results_emAlg_p02_300_time160",seqIndex[chunk],".RData",sep="")
# filename = paste("results_emAlg_",2512,".RData",sep="")

# save(results,file = filename)




# # for(i in 1:length(simulation)){
#   sim = simulation[[chunk]]
#   net0 = sim$net0
#   net1 = sim$net1
#
#   actDf <-dimnames(net0)[[2]]
#   actDfnodes <- defineNodes(data.frame(label=colnames(net0)))
#
#   nAct <- length(actDf)
#   endTime <- 50
#   nSimCrea <- 59 # from s501 and s502 data, number of ties created between t_0 and t_1
#   nSimDel <- 56 # from s501 and s502 data, number of ties deleted between t_0 and t_1
#   parmsCrea <- c(log(nSimCrea / 20 / nAct )) #indegree and outdegree of ego
#   parmsDel <- c(log(nSimDel / 20 / (nAct-4))) #indegree and outdegree of ego
#
#   theta0=data.frame("Crea"=parmsCrea,"Del"=parmsDel)
#   beta0=data.frame("Crea"=c(0,0,0,0),"Del"=c(0,0,0,0))
#   row.names(beta0) = c("indeg","outdeg","recip","trans")
#   formula = "indeg + outdeg + recip + trans +inertia + tie(net0)"
#
#
#   results[[chunk]] = lapply(ceiling(10^seq(1,1.1,0.2)),softEMAlgorithmMC,net0=net0,net1=net1,
#                         theta0=theta0,beta0=beta0,formula = formula,num_cores=38)
#
# # }







# net0= initState
# rownames(net0)=colnames(net0)
# net1=Xstate
# actDf <-dimnames(net0)[[2]]
# actDfnodes <- defineNodes(data.frame(label=colnames(net0)))
# nAct <- length(actDf)
# endTime <- 50
# nSimCrea <- 59 # from s501 and s502 data, number of ties created between t_0 and t_1
# nSimDel <- 56 # from s501 and s502 data, number of ties deleted between t_0 and t_1
# parmsCrea <- c(log(nSimCrea / endTime / nAct )) #indegree and outdegree of ego
# parmsDel <- c(log(nSimDel / endTime / (nAct-4))) #indegree and outdegree of ego
# theta0=data.frame("Crea"=parmsCrea,"Del"=parmsDel)
# beta0=data.frame("Crea"=c(0,0,0,0),"Del"=c(0,0,0,0))
# row.names(beta0) = c("indeg","outdeg","recip","trans")
# formula = "indeg + outdeg + recip + trans +inertia + tie(net0)"
# time1 = Sys.time()
# out1 = softEMAlgorithmMC(nmax = 8,net0 = net0,net1 = net1,theta0 = theta0,beta0 = beta0,
#                          formula = formula,num_cores=8)
# time2 = Sys.time() - time1
#

#
#
# time3 = Sys.time()
# out2 = softEMAlgorithmForEach(net0,net1,theta0,beta0,formula,nmax = 50)
# time4 = Sys.time()-time3


# out = vector("list",1000)
# time1 = Sys.time()
# for(j in 1:1000) {
#   cat("Iteration: ",j,"\n")
#   set.seed(j)
#   out[[j]] =softEMAlgorithmMC(net0,net1,theta0,beta0,formula,num_cores=17,nmax = 100)
#    # out[[j]] =softEMAlgorithm(net0,net1,theta0,beta0,formula,nmax = 50)
#
# }
# time2 = Sys.time() - time1
# save(out,time2,file = "out_loop_1000_nmax_100_17cores.RData")



# Test original sequence -----
# out$beta
# coef(modTest2)
# coef(modTest4)
# which.max(out$logLik) # which is the original one!!
# plot(out$logLik) # way smaller than the others!!
#
# save.image(file='testOriginalSequence_100seq.RData')




# load("out_loop_10000_nmax_50_17cores.RData")
#
# betaDF2 = lapply(out,"[[","beta")
# betaDF2
#
# betaDF= c(betaDF,betaDF2)
# # seDF
# creation = lapply(betaDF,"[[",1)
# indegC = sapply(creation,"[[",1)
# outdegC= sapply(creation,"[[",2)
# recipC = sapply(creation,"[[",3)
# deletion = lapply(betaDF,"[[",2)
# indegD = sapply(deletion,"[[",1)
# outdegD= sapply(deletion,"[[",2)
# recipD = sapply(deletion,"[[",3)
# #
# par(mfrow=c(2,3))
# hist(indegC)
# hist(outdegC)
# hist(recipC)
# hist(indegD)
# hist(outdegD)
# hist(recipD)
#
# library(pastecs)
# round(stat.desc(indegC,basic=F),3)
# round(stat.desc(outdegC,basic = F),3)
# round(stat.desc(recipC,basic = F),3)
# round(stat.desc(indegD,basic = F),3)
# round(stat.desc(outdegD,basic = F),3)
# round(stat.desc(recipD,basic = F),3)




