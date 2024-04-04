# Ascent-based MCEM -----------------------------
sigmaHat = function(mt,loglikPrev, loglikCur, w){

  lambda = loglikCur - loglikPrev
  sumlambdaw = sum(w*lambda)
  sumlambdaw2 = sum((w*lambda)^2)
  sumw2lambda = sum((w^2)*lambda)

  sigma2 = mt*sumlambdaw^2*(sumlambdaw2/sumlambdaw^2 - 2*sumw2lambda/sumlambdaw + sum(w^2))
  return(sigma2)

}

deltaQ = function(loglikPrev, loglikCur, w){

  loglikRatio = loglikCur - loglikPrev
  deltaQ = sum(w*loglikRatio)

  return(deltaQ)
}


# Ascent-Based Markov-Chain Expectation-Maximization algorithm --------------

MCEMalgorithm = function(nmax0,net0,net1,theta0,beta0,formula,num_cores=1){
  # net0 : network matrix in initial state
  # net1: network matrix in final state
  # theta0 : initial parameters for sender (list of creation and deletion with
  #          labels Crea and Del)
  # beta0: initial parameters for receiver (list of creation and deletion with
  #        labels Crea and Del)
  # formula: string with formula for the model (only rhs of formula)

  actDfnodes <- defineNodes(data.frame(label=colnames(net0)))

  nmax = nmax0

  # Initialization of parameters
  theta = theta0
  beta = beta0
  se=data.frame(Crea=rep(0,nrow(beta)),Del=rep(0,nrow(beta)))
  row.names(se) = row.names(beta)
  # Creation of sequence of events from initial data
  sequence = EMPreprocessing(net0,net1)
  H = length(sequence) # Humming distance
  # Creation of permutations
  #permut = permute(sequence, nmax = nmax)
  permut = MCMC(nmax=nmax,seq,H,actDfnodes,formula,net0,beta,burnIn = TRUE,maxIter=10000,seqIter = 50,pShort=0.35,pAug = 0.35,pPerm = 0.3)

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

  logLikPrev = clusterApply(cl,seq_along(splitIndicesPerCore),logLikelihoodMC,permut = permut,
                            splitIndicesPerCore=splitIndicesPerCore,
                            beta=beta,actDfnodes=actDfnodes,net0=net0,formula=formula)
  logLikPrev = unlist(logLikPrev)

  diff = 1000
  index = 0
  while (diff>1e-3){
    # browser()

    cat("Index: ",index,"\n")
    cat("Diff: ", diff,"\n")
    if (index > 1000){
      cat("No convergence\n")
      break
    }

    c = max(logLikPrev)
    logsumExp = c + log(sum(exp(logLikPrev-c)))
    weights = exp(logLikPrev - logsumExp)

    betaCreaAux = rubinsRule(betaCreaDF,seCreaDF,weights)
    betaDelAux = rubinsRule(betaDelDF,seDelDF,weights)
    #diff= sqrt(norm(as.matrix(beta$Crea-betaCreaAux$mean))^2+norm(as.matrix(beta$Del-betaDelAux$mean))^2)

    betaNew$Crea = betaCreaAux$mean
    betaNew$Del = betaDelAux$mean


    logLikCur = clusterApply(cl,seq_along(splitIndicesPerCore),logLikelihoodMC,permut = permut,
                             splitIndicesPerCore=splitIndicesPerCore,
                             beta=betaNew,actDfnodes=actDfnodes,net0=net0,formula=formula)

    w = rep(1,length(logLikPrev))/length(logLikPrev)
    sigmaHat = sigmaHat(nmax,logLikPrev,logLikCur,w)
    ase = sqrt(sigmaHat/nmax)
    deltaQ = deltaQ(logLikPrev,logLikCur,w)
    lowerBound = deltaQ - 1.281552*ase # 80%



    if(lowerBoud<0){
      # Estimator is not accepted, new point must be sampled
      newpermut = MCMC(nmax=nmax/2,permut[[length(permut)]],H,actDfnodes,formula,net0,beta,burnIn = FALSE,maxIter=10000,seqIter = 50,pShort=0.35,pAug = 0.35,pPerm = 0.3)
      nmax = nmax+(nmax0/2)
      permut = c(permut,newpermut)

      resPar = clusterApply(cl,seq_along(splitIndicesPerCore),parametersMC,permut = permut,
                            splitIndicesPerCore=splitIndicesPerCore,actDfnodes.=actDfnodes,
                            net0.=net0,formula.=formula)

      resPar = unlist(resPar,recursive = FALSE)
      betaCreaDF = t(sapply(resPar,"[[",1))
      betaDelDF = t(sapply(resPar,"[[",2))
      seCreaDF = t(sapply(resPar,"[[",3))
      seDelDF = t(sapply(resPar,"[[",4))

      logLikPrev = clusterApply(cl,seq_along(splitIndicesPerCore),logLikelihoodMC,permut = permut,
                                splitIndicesPerCore=splitIndicesPerCore,
                                beta=beta,actDfnodes=actDfnodes,net0=net0,formula=formula)
      logLikPrev = unlist(logLikPrev)

    }else{
      # Update of the paremeter, next iteration
      index = index +1
      diff = deltaQ + 1.644854*ase #90%


      beta = betaNew
      se$Crea = betaCreaAux$se
      se$Del = betaDelAux$se


      # Update on the permutations -> new MCMC
      errType2 = 0.8
      m_start = sigmaHat^2*(1.281552+qnorm(errType2))^2/deltaQ^2

      if(m_start > nmax){
        nmax = m_start
      }

      indexMaxlogLik = which(logLikPrev==max(logLikPrev))

      permut = MCMC(nmax=nmax,permut[[indexMaxlogLik]],H,actDfnodes,formula,net0,beta,burnIn = TRUE,maxIter=10000,seqIter = 50,pShort=0.35,pAug = 0.35,pPerm = 0.3)

      logLikPrev = clusterApply(cl,seq_along(splitIndicesPerCore),logLikelihoodMC,permut = permut,
                                splitIndicesPerCore=splitIndicesPerCore,
                                beta=beta,actDfnodes=actDfnodes,net0=net0,formula=formula)
      logLikPrev = unlist(logLikPrev)

    }

  }

  aux_permut = permut[[1]]
  og_permut = aux_permut[order(as.numeric(row.names(aux_permut))),]
  aux_names = lapply(lapply(permut,row.names),as.numeric)

  return(list("logLik" = logLik,"beta" = beta,"se" = se,"index" = index,
              "diff" = diff,"og_permut" =og_permut, "permut"=aux_names,"betaCreaDF" = betaCreaDF,"betaDelDF"=betaDelDF))
}


# nmax = 10





