rm(list=ls())
# packages ----------------------------------------------------------------
library(goldfish)
library(broom)
library(pixiedust)
library(survival)
library(mlogit)
library(RSiena)
library(matrixStats)
library(stats)
library(foreach)
library(doParallel)
library(parallel)
  
# Setup of the EM algorithm ---------------
EMPreprocessing = function(X0, X1){
  # Function that preprocess the network matrixes
  # X0 : initial network matrix (row, col names must be actor names)
  # X1 : final network matrix (row, col names must be actor names)
  actors = colnames(X0)  
  sequence = which(X0!=X1,arr.ind = TRUE)
  colnames(sequence)  = c("sender","receiver")
  sequence = as.data.frame(sequence)
  sequence$replace = ifelse(X0[X0!=X1]==1,0,1)
  sequence$sender = actors[sequence$sender]
  sequence$receiver = actors[sequence$receiver]
  rownames(sequence) = 1:nrow(sequence)
  return(sequence)
}

  
# Permutations -------------------------------------
permute = function(x,nmax = 1000){
  # x: original data frame with a column of senders, a column of receivers and a column of replace
  #    (replace = 0 means deletion of tie, replace = 1 means creation of tie)
  # nmax: maximum number of permutations performed. The number of permutations
  #       computed will be the min{length(x)!, nmax}
  # out: output of the function is a list of dataframes. Each dataframe is a permutation
  #      of the rows of the original list x
  n = nrow(x)
  if(factorial(n)<nmax){
    out = lapply(1:factorial(n),function(t) x[sample(n,n),])
  }else{
    out = lapply(1:nmax,function(t) x[sample(n,n),])
  }
  return(out)
}


# Function to create time difference of events ---------
timeGenerator = function(seq,nAct,theta){
  # seq : sequence of events
  # nAct : number of actors
  # theta : list with Crea and Del objects with parameters for the sender
    
  X_aux <- cbind("intercept"=rep(1,nAct))
    
  expXbCrea <- exp(X_aux %*% theta$Crea)
  sumRateCrea <- sum(expXbCrea)
  expXbDel <- exp(X_aux %*% theta$Del)
  sumRateDel <- sum(expXbDel)
    
  time = numeric(nrow(seq))
  if(seq$replace[1]==0){
    time[1] =  rexp(1, sumRateDel)
  }else{
    time[1] = rexp(1, sumRateCrea)
  }
    
  for(i in 2:nrow(seq)){
    if(seq$replace[i]==0){
      time[i] = time[i-1] + rexp(1, sumRateDel)
    }else{
      time[i] = time[i-1] + rexp(1, sumRateCrea)
    }
  }
  return(time)
}
  
 
GatherPreprocessingDF = function(formula, envir = new.env()){
  # formula : string formula with the effects for the model
  # envir : environment were dependentEvents, nodes and net0 are located.
  # return: listExpandedDF
  
  # browser()
  formGather <- paste("depEvents ~", formula, sep = "")
  
  dataProcessed <- GatherPreprocessing(
    as.formula(formGather),
    model = "DyNAM", subModel = "choice",
    progress = FALSE,
    envir = envir
  )
  namesEffects <- dataProcessed$namesEffects
  nEvents <- length(dataProcessed$n_candidates)
  
  expandedDF <- cbind(
    setNames(
      as.data.frame(dataProcessed$stat_all_events),
      namesEffects
    ),
    data.frame(
      event = rep(seq.int(nEvents), dataProcessed$n_candidates),
      selected = sequence(dataProcessed$n_candidates) ==
        rep(dataProcessed$selected, dataProcessed$n_candidates),
      sender = rep(dataProcessed$sender, dataProcessed$n_candidates)
    )
  
  )
  
  nEvents = expandedDF$event[length(expandedDF$event)]
  indexEvents = data.frame("start"=seq(1,length(expandedDF$event),nAct-1),
                           "end"=seq(nAct-1,length(expandedDF$event),nAct-1),
                           "event"=seq(1,nEvents))
  creation_events = 0
  deletion_events = 0
  for(i in 1:nEvents){
  # Select rows of expandedDF corresponding to the event
  # Check if inertia == 0 (creation of ties) or inertia == 1 (delition of ties) for (receiver == TRUE)
   df_aux = expandedDF[seq(indexEvents$start[i],indexEvents$end[i]),]
      index_aux = df_aux[,"selected"]==TRUE
      if (df_aux[index_aux,"inertia_netEvents"]==0){ #creation of ties
         creation_events=c(creation_events,i)
      }else{
         deletion_events=c(deletion_events,i)
       }
    }
    creation_events=creation_events[-1]
    deletion_events=deletion_events[-1]
    
    expandedDFCreation = subset(expandedDF,event %in% creation_events & inertia_netEvents == 0)
    expandedDFDeletion = subset(expandedDF,event %in% deletion_events & inertia_netEvents == 1)
  
    listExpandedDF = list(
    expandedDFCreation = expandedDFCreation[,!names(expandedDFCreation) %in% c("inertia_netEvents", "tie_net0","sender")],
    expandedDFDeletion = expandedDFDeletion[,!names(expandedDFDeletion) %in% c("inertia_netEvents", "tie_net0","sender")])
    
  return(listExpandedDF)
}  
  
# Log-likelihood computation function ----------------

logLikelihood = function(listExpandedDF, beta){
  # listExpandedDF : list returned from function GatherePreprocessingDF
  # beta : choice parameters 
  # return : loglikelihood (creation + deletion of events models)
  
  # Computations of loglikelihood for creation events
  xbCrea = by(listExpandedDF$expandedDFCreation[,!names(listExpandedDF$expandedDFCreation) %in% c("event","selected")],
                  listExpandedDF$expandedDFCreation[,"event"],function(x) as.matrix(x)%*%beta$Crea)
  xb_auxCrea = by(listExpandedDF$expandedDFCreation[listExpandedDF$expandedDFCreation$selected==TRUE,!names(listExpandedDF$expandedDFCreation) %in% c("event","selected")],
              listExpandedDF$expandedDFCreation[listExpandedDF$expandedDFCreation$selected==TRUE,"event"],function(x) as.matrix(x)%*%beta$Crea)
  loglikCrea =  sum(unlist(xb_auxCrea) - sapply(xbCrea,colLogSumExps))

  # Computation of loglikelihood for deletion events
  xbDel = by(listExpandedDF$expandedDFDeletion[,!names(listExpandedDF$expandedDFDeletion) %in% c("event","selected")],
          listExpandedDF$expandedDFDeletion[,"event"],function(x) as.matrix(x)%*%beta$Del)
  xb_auxDel = by(listExpandedDF$expandedDFDeletion[listExpandedDF$expandedDFDeletion$selected==TRUE,!names(listExpandedDF$expandedDFDeletion) %in% c("event","selected")],
              listExpandedDF$expandedDFDeletion[listExpandedDF$expandedDFDeletion$selected==TRUE,"event"],function(x) as.matrix(x)%*%beta$Del)
  loglikDel = sum(unlist(xb_auxDel) - sapply(xbDel,colLogSumExps))
  
  return(loglikCrea+loglikDel)
  
}

logLikelihoodMC = function(indexCore,permut,beta,splitIndicesPerCore,actDfnodes=actDfnodes,net0=net0,formula=formula){
      
  indicesCore = splitIndicesPerCore[[indexCore]]
  resLikelihood = vector("list",length(indicesCore))
  
  for(i in seq_along(indicesCore)){
      seq = permut[[indicesCore[[i]]]]
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
      
      resLikelihood[[i]] = logLikelihood(listExpandedDF,beta)
    }
  return(unlist(resLikelihood))
}

# Rubin's rule for weighted data -------------
rubinsRule = function(x,se,w){
  # x : values of the parameters 
  # se : standard errors of the parameters
  # w : values of the weights
  # returns the value of the averaged mean and the standard error
  x_mean = apply(x,2,weighted.mean,w)
  v_within = apply(se^2,2,weighted.mean,w)
  v_between = apply((sweep(x,2,x_mean))^2,2,weighted.mean,w)
  se_total = sqrt(v_within+(1+1/length(w))*v_between)
  return(list(mean = x_mean, se = se_total))
}



# Parameters ----------

parameters = function(seq,actDfnodes.=actDfnodes,net0.=net0,formula.=formula){
  envirPreproCrea <- new.env()
  envirPreproDel <- new.env()

  seq$time <- seq(1,nrow(seq))
  
  envirPreproCrea$seqTime <- seq
  envirPreproCrea$actDfnodes <- actDfnodes.
  envirPreproCrea$net0 <- net0.
  
  envirPreproDel$seqTime <- seq
  envirPreproDel$actDfnodes <- actDfnodes.
  envirPreproDel$net0 <- net0.
  
  
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
  
  formCrea = paste("depEventsCrea ~",formula.,sep="")
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

  formDel = paste("depEventsDel ~",formula.,sep="")
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
  return(list(betaCreaDF=coef(modCrea),betaDelDF = coef(modDel),
         seCreaDF = modCrea$standardErrors[modCrea$standardErrors>0],
         seDelDF = modDel$standardErrors[modDel$standardErrors>0]))    
}



parametersMC = function(indexCore,splitIndicesPerCore,permut = permut,actDfnodes.=actDfnodes,net0.=net0,formula.=formula){
  
  indicesCore = splitIndicesPerCore[[indexCore]]
  resParameters = vector("list",length(indicesCore))
  for(i in seq_along(indicesCore)){
    resParameters[[i]] = parameters(permut[[indicesCore[[i]]]],actDfnodes.,net0.,formula.)
  }
  
  return(resParameters)
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
    weights = exp(logLik)/sum(logLik)
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







softEMAlgorithmMC = function(net0,net1,theta0,beta0,formula,num_cores=1,nmax = 10){
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
  while (diff>1e-4){ 
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
    
    weights = exp(logLik)/sum(logLik)
    betaCreaAux = rubinsRule(betaCreaDF,seCreaDF,weights)
    betaDelAux = rubinsRule(betaDelDF,seDelDF,weights)
    diff= sqrt(norm(as.matrix(beta$Crea-betaCreaAux$mean))^2+norm(as.matrix(beta$Del-betaDelAux$mean))^2)

    beta$Crea = betaCreaAux$mean
    beta$Del = betaDelAux$mean
    se$Crea = betaCreaAux$se
    se$Del = betaDelAux$se
    }

  return(list("logLik" = logLik,"beta" = beta,"se" = se,"betaCrea" = betaCreaDF,"betaDel" = betaDelDF,"index" = index,
              "diff" = diff,"permut"=permut))
}
   




# load("/Users/mariaeugeniagilpallares/MEGA/master/HS23/Semester Paper/emDyNAM/simulationData.RData")
load("~/euler/simulationData.RData")

net0= initState
rownames(net0)=colnames(net0)
net1=Xstate

actDf <-dimnames(net0)[[2]]
actDfnodes <- defineNodes(data.frame(label=colnames(net0)))

nAct <- length(actDf)
endTime <- 30
nSimCrea <- 59 # from s501 and s502 data, number of ties created between t_0 and t_1
nSimDel <- 56 # from s501 and s502 data, number of ties deleted between t_0 and t_1
parmsCrea <- c(log(nSimCrea / endTime / nAct )) #indegree and outdegree of ego
parmsDel <- c(log(nSimDel / endTime / (nAct-4))) #indegree and outdegree of ego


theta0=data.frame("Crea"=parmsCrea,"Del"=parmsDel)
beta0=data.frame("Crea"=c(0,0,0),"Del"=c(0,0,0))
row.names(beta0) = c("indeg","outdeg","recip")
formula = "indeg + outdeg + recip + inertia + tie(net0)"




# time1 = Sys.time()
# out = softEMAlgorithmMC(net0,net1,theta0,beta0,formula,num_cores=8,nmax = 50)
# time2 = Sys.time() - time1
# 




out = vector("list",10)
time1 = Sys.time()
for(j in 1:100) {
  cat("Iteration: ",j,"\n")
  set.seed(j)
  out[[j]] =softEMAlgorithmMC(net0,net1,theta0,beta0,formula,num_cores=8,nmax = 50)
}
time2 = Sys.time() - time1


save(likelihood,time2,file = "out_loop_nmax_50.RData")




# seDF
# creation = lapply(betaDF,"[[",1)
# indegC = sapply(creation,"[[",1)
# outdegC= sapply(creation,"[[",2)
# recipC = sapply(creation,"[[",3)
# deletion = lapply(betaDF,"[[",2)
# indegD = sapply(deletion,"[[",1)
# outdegD= sapply(deletion,"[[",2)
# recipD = sapply(deletion,"[[",3)
# 
# par(mfrow=c(2,3))
# hist(indegC)
# hist(outdegC)
# hist(recipC)
# hist(indegD)
# hist(outdegD)
# hist(recipD)
# 
# library(pastecs)
# round(stat.desc(indegC),3)
# round(stat.desc(outdegC),3)
# round(stat.desc(recipC),3)
# round(stat.desc(indegD),3)
# round(stat.desc(outdegD),3)
# round(stat.desc(recipD),3)
# 
