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

# Setup of the EM algorithm ---------------
EMPreprocessing = function(X0,X1){
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

#' Permutaciones
#'
#' @param x vector
#' @param nmax integer
#'
#' @return
#' @export
#'
#' @examples permute(1:5,nmax = 5)
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
      fixedParameters = c(NA, NA, NA, NA,-20,-20),
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
      fixedParameters = c(NA, NA, NA, NA,20, 20),
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

# MCMC auxiliar functions --------------------

getKelMeMatrix = function(seq,actDfnodes){
  seq = cbind(seq,"row" = as.integer(rownames(seq)))


  auxDf = lapply(actDfnodes,function(x) lapply(actDfnodes, function(i) subset(seq,subset=(sender== x & receiver==i))))

  for(i in 1:length(actDfnodes)){
    for(j in 1:length(actDfnodes)){
      if(nrow(auxDf[[i]][[j]])>0){
        auxDf[[i]][[j]]$row = as.integer(auxDf[[i]][[j]]$row)
        auxDf[[i]][[j]]$rowDiff =c(0,diff(auxDf[[i]][[j]]$row))
      }
    }
  }


  Kel_g1 = matrix(0,nrow=length(actDfnodes),ncol=length(actDfnodes)) # Notation from Johan´s paper, sum for l greater than 1
  Kel_ge1 = matrix(0,nrow=length(actDfnodes),ncol=length(actDfnodes)) # Notation from Johan´s paper, sum for l greater or equal than 1

  for(i in 1:nrow(Kel_g1)){
    for(j in 1:ncol(Kel_g1)){
      if(nrow(auxDf[[i]][[j]])>0){
        Kel_ge1[i,j] = nrow(auxDf[[i]][[j]][auxDf[[i]][[j]]$rowDiff!=1,])
        v_rowDiff = auxDf[[i]][[j]]$rowDiff
        not_1 <- which(v_rowDiff != 1)
        shifted_vector <- c(v_rowDiff[-1], 0)
        Kel_g1[i,j] = sum(v_rowDiff[-length(v_rowDiff)] != 1 & shifted_vector[-1] == 1)
      }
    }
  }

  me = t(sapply(auxDf,function(x) sapply(x,nrow)))

  return(list(Kel_g1 = Kel_g1, Kel_ge1 = Kel_ge1, me = me, auxDf = auxDf))
}

getAuxDfE = function(auxDf, sender,receiver){
  auxDfE = auxDf[[sender]][[receiver]]
  indexNo1 = which(auxDfE$rowDiff!=1)
  auxDfE$run = rep(0,nrow(auxDfE))
  run=1
  for(i in 1:(length(indexNo1)-1)){
    auxDfE$run[indexNo1[i]:indexNo1[i+1]]= run
    run=run+1
  }
  auxDfE$run[indexNo1[length(indexNo1)]:nrow(auxDfE)]=run

  return(auxDfE)
  }


stepAugment = function(seq,tieNames,gammaEplus,gammaPlus,m,me,net0){
  # Choose element to be inserted:
  e = sample(tieNames,size=1, prob=as.vector(t(gammaEplus/gammaPlus)))

  # Choose to augment in the same e-run or different ones
  sender = as.integer(strsplit(e,"V")[[1]][2])
  receiver = as.integer(strsplit(e,"V")[[1]][3])

  p = (m-me[sender,receiver]+1)/gammaEplus[sender,receiver]
  typeA = sample(c("same","diff"),size=1,prob=c(p,1-p))
  indexSR = which(seq$sender==paste("V",sender,sep="") & seq$receiver==paste("V",receiver,sep=""))



  if(typeA == "same"){
    if(length(indexSR)<1){
      place = sample(m+1, size=1)
    }else{
      place = sample(seq(m+1)[-indexSR],size=1)
    }


    if(place==1){
      newseq = rbind(c(paste("V",sender,sep=""),paste("V",receiver,sep=""),4,place),
                     c(paste("V",sender,sep=""),paste("V",receiver,sep=""),4,place+1),
                     seq)
    }else{
       newseq = rbind(seq[1:(place-1),],
                   c(paste("V",sender,sep=""),paste("V",receiver,sep=""),4,place),
                   c(paste("V",sender,sep=""),paste("V",receiver,sep=""),4,place+1),
                   seq[(place):m,])
    }

    nEdges = which(newseq$sender==paste("V",sender,sep="") & newseq$receiver==paste("V",receiver,sep=""))
    initAction = as.integer(seq[which(seq$sender==paste("V",sender,sep="") & seq$receiver==paste("V",receiver,sep="")),"replace"][1])
    if(is.na(initAction)) {initAction=net0[sender,receiver]}
    newseq[nEdges,"replace"] = seq(initAction,(initAction+length(nEdges)-1)) %% 2
    newseq$row = seq(1,nrow(newseq))

    pDoStep = (gammaEplus[sender,receiver]/gammaPlus) * p * (1/m) # change

  }else if(typeA == "diff"){
    # Choose thow places separated at least from one tie different than e
    if(length(indexSR)<1){
      place = sample(m+1, size=2)
    }else{
      place = sample(seq(m+1)[-indexSR],size=2)
    }

    place = sort(place)



    if(place[1]==1){
      newseq = rbind(c(paste("V",sender,sep=""),paste("V",receiver,sep=""),4,place[1]+1),
                     seq[(place[1]):(place[2]-1),],
                     c(paste("V",sender,sep=""),paste("V",receiver,sep=""),4,place[2]+2),
                     seq[(place[2]):m,])
    }else{
    newseq = rbind(seq[1:(place[1]-1),],
                   c(paste("V",sender,sep=""),paste("V",receiver,sep=""),4,place[1]+1),
                   seq[(place[1]):(place[2]-1),],
                   c(paste("V",sender,sep=""),paste("V",receiver,sep=""),4,place[2]+2),
                   seq[(place[2]):m,])
    }
    nEdges = which(newseq$sender==paste("V",sender,sep="") & newseq$receiver==paste("V",receiver,sep=""))
    initAction = as.integer(seq[which(seq$sender==paste("V",sender,sep="") & seq$receiver==paste("V",receiver,sep="")),"replace"][1])
    if(is.na(initAction)) {initAction=net0[sender,receiver]}
    newseq[nEdges,"replace"] = seq(initAction,(initAction+length(nEdges)-1)) %% 2
    newseq$row = seq(1,nrow(newseq))

    pDoStep = (gammaEplus[sender,receiver]/gammaPlus) * (1-p) * (1/m)*(1/(m-3)) # to be done

  }

  return(list(newseq=newseq,sender=paste("V",sender,sep=""),receiver=paste("V",receiver,sep=""),place = place,typeA=typeA,
              pDoStep = pDoStep))
}

stepShort = function(seq,tieNames,gammaEminus,gammaMinus,m,me,auxDf,net0){
  # Choose element to be deleted:
  e = sample(tieNames,size=1, prob=as.vector(t(gammaEminus/gammaMinus)))

  # Choose to remove from the same e-run or from two distinct e-runs
  sender = as.integer(strsplit(e,"V")[[1]][2])
  receiver = as.integer(strsplit(e,"V")[[1]][3])
  p = Kel_g1[sender,receiver]/gammaEminus[sender,receiver]
  typeS = sample(c("same","diff"),size=1,prob=c(p,1-p))

  # Get auxiliar variables with indexes for the runs
  auxDfE = getAuxDfE(auxDf, sender,receiver)


  if(typeS=="same"){
    # Choose the e-runs of more than 1 element, select one and delete the first 2 elements.
    sampleRun = sample(which(table(auxDfE$run)>1),size=1)
    indexSampleRun = auxDfE$row[auxDfE$run==sampleRun][1:2]

    newseq = seq[-(indexSampleRun),]
    newseq$row = 1:nrow(newseq)

    place = list(sampleRun=sampleRun, indexSampleRun= indexSampleRun)

    pDoStep = (gammaEminus[sender,receiver]/gammaMinus) * p * (1/length(unique(auxDfE$run)))

  }else if(typeS=="diff"){
    # Choose two different e-runs, and delete the first element of each one of them
    sampleRun = sample(unique(auxDfE$run),size=2,replace=FALSE)
    indexSampleRun = c(auxDfE$row[auxDfE$run==sampleRun[1]][1],auxDfE$row[auxDfE$run==sampleRun[2]][1])

    newseq = seq[-(indexSampleRun),]
    newseq$row = 1:nrow(newseq)

    # now we have to fix the replaceOG
    nEdges = which(newseq$sender==paste("V",sender,sep="") & newseq$receiver==paste("V",receiver,sep=""))
    initAction = as.integer(seq[which(seq$sender==paste("V",sender,sep="") & seq$receiver==paste("V",receiver,sep="")),"replace"][1])
    newseq[nEdges,"replace"] = seq(initAction,(initAction+length(nEdges)-1)) %% 2

    place = list(sampleRun=sampleRun, indexSampleRun= indexSampleRun)
    pDoStep = (gammaEminus[sender,receiver]/gammaMinus) * (1-p) * (1/length(unique(auxDfE$run)))*(1/(length(unique(auxDfE$run))-1))
  }

  return(list(newseq=newseq,sender=paste("V",sender,sep=""),receiver=paste("V",receiver,sep=""),place = place,typeS=typeS,
              pDoStep = pDoStep, auxDfE = auxDfE))
}

stepPerm = function(seq,tieNames,m,me){
  # Choose elements to be permuted:
  probVec = as.vector(t(me/m))
  probVec2 = as.vector(t(me/(m-me)))
  names(probVec) = tieNames
  names(probVec2) = tieNames
  e1 = sample(tieNames,size=1, prob=probVec)
  e2 = sample(tieNames[tieNames!=e1],size=1, prob=probVec2[names(probVec)!=e1])

  # After getting the edge identifier, select at random e1 and e2 from all the runs
  sender1 = paste("V",strsplit(e1,"V")[[1]][2],sep="")
  receiver1 = paste("V",strsplit(e1,"V")[[1]][3],sep="")
  place1 = as.integer(sample(seq[seq$sender==sender1 & seq$receiver==receiver1,]$row,size=1))
  sender2 = paste("V",strsplit(e2,"V")[[1]][2],sep="")
  receiver2 = paste("V",strsplit(e2,"V")[[1]][3],sep="")
  place2 = as.integer(sample(seq[seq$sender==sender2 & seq$receiver==receiver2,]$row,size=1))
  auxRow = seq[place1,]
  newseq = seq
  newseq[place1,] = newseq[place2,]
  newseq[place2,] = auxRow
  rm(auxRow)


  nEdges = which(newseq$sender==sender1 & newseq$receiver==receiver1)
  initAction = as.integer(seq[which(seq$sender==sender1 & seq$receiver==receiver1),"replace"][1])
  newseq[nEdges,"replace"] = seq(initAction,(initAction+length(nEdges)-1)) %% 2
  nEdges = which(newseq$sender==sender2 & newseq$receiver==receiver2)
  initAction = as.integer(seq[which(seq$sender==sender2 & seq$receiver==receiver2),"replace"][1])
  newseq[nEdges,"replace"] = seq(initAction,(initAction+length(nEdges)-1)) %% 2

  newseq$row = 1:nrow(newseq)

  pDoStep = probVec[e1] * probVec2[e2] * (1/nrow(seq[seq$sender==sender1 & seq$receiver==receiver1,])) * (1/nrow(seq[seq$sender==sender2 & seq$receiver==receiver2,]))

  return(list(newseq=newseq,sender=c(sender1,sender2),receiver=c(receiver1,receiver2),place = c(place1,place2),
              pDoStep = pDoStep))
}



# MCMC -------------------------------------


burnIn = function(seq,beta,H,actDfnodes,formula,net0){
  #burn in for MCMC (to be done only once to change all sequences)

}

MCMC = function(seq,burn_in,H,actDfnodes,n=500){
  # Compute initial quatities:
  # Type 1: augmentation, type 2: shortening, type 3: permutation

  pShort = pAug = 0.35 # These could be parameters of the function
  pPerm = 0.3 # This could be a parameter of the function

  actDfnodesLab=actDfnodes$label

  tieNames = sapply(actDfnodesLab,function(x) sapply(actDfnodesLab, function(i) paste(x,i,sep="")))
  tieNames = as.vector(tieNames)


  getKelMeMatrix = getKelMeMatrix(seq,actDfnodesLab)
  Kel_g1 = getKelMeMatrix$Kel_g1
  Kel_ge1 = getKelMeMatrix$Kel_ge1
  gammaEminus = choose(Kel_ge1,2) + Kel_g1
  gammaMinus = sum(gammaEminus)
  me = getKelMeMatrix$me
  gammaEplus = choose(nrow(seq)-me+2,2)
  gammaPlus = sum(gammaEplus)
  m = nrow(seq)
  auxDf = getKelMeMatrix$auxDf

  for(i in 1:1000){
    acceptIndex = 0

   if(length(seq)==H){
      type = sample(c(1,2,3),size=1,prob=c(pAug/(pAug+pPerm),0,pPerm/(pAug+pPerm)))
    }else{
      type = sample(c(1,2,3),size=1,prob=c(pAug,pShort,pPerm))
    }


   if(type == 1){ # Augmentation
      step = stepAugment(seq,tieNames,gammaEplus,gammaPlus,m,me)

      getNewKelMeMatrix = getKelMeMatrix(step$newseq,actDfnodesLab)
      newAuxDfE = getAuxDfE(getNewKelMeMatrix$auxDf, sender,receiver)
      auxDfE = getAuxDfE(auxDf,sender,receiver)

      newKel_g1 = getNewKelMeMatrix$Kel_g1
      newKel_ge1 = getNewKelMeMatrix$Kel_ge1
      newGammaEminus = choose(newKel_ge1,2) + newKel_g1
      newGammaMinus = sum(newGammaEminus)
      newMe = getNewKelMeMatrix$me
      newGammaEplus = choose(nrow(newseq)-newMe+2,2)
      newGammaPlus = sum(newGammaEplus)
      newM = nrow(newseq)

      newp= newKel_g1[sender,receiver]/newGammaEminus[sender,receiver]
      if((length(unique(auxDfE$run)) == length(unique(newAuxDfE$run))) | step$typeA=="same"){
          pUndoStep = (newGammaEminus[sender,receiver]/newGammaMinus) * newp * (1/length(unique(newAuxDfE$run)))
         }else{
          pUndoStep = (newGammaEminus[sender,receiver]/newGammaMinus) * (1-newp) * (1/length(unique(newauxDfE$run)))*(1/(length(unique(newAuxDfE$run))-1))
         }

      loglikSeq = logLikelihoodMC(1,seq,beta,list(1),actDfnodes=actDfnodes,net0=net0,formula=formula)
      newloglikSeq = logLikelihoodMC(indexCore=1,newseq,beta,splitIndicesPerCore=list(1),actDfnodes=actDfnodes,net0=net0,formula=formula)


   }else if(type == 2){ # Shortening
     step = stepShort(seq,tieNames,gammaEminus,gammaMinus,m,me,auxDf)
     getNewKelMeMatrix = getKelMeMatrix(step$newseq,actDfnodesLab)
     newAuxDfE = getAuxDfE(getNewKelMeMatrix$auxDf, sender,receiver)
     auxDfE = getAuxDfE(auxDf,sender,receiver)

     newKel_g1 = getNewKelMeMatrix$Kel_g1
     newKel_ge1 = getNewKelMeMatrix$Kel_ge1
     newGammaEminus = choose(newKel_ge1,2) + newKel_g1
     newGammaMinus = sum(newGammaEminus)
     newMe = getNewKelMeMatrix$me
     newGammaEplus = choose(nrow(newseq)-newMe+2,2)
     newGammaPlus = sum(newGammaEplus)
     newM = nrow(newseq)
     senderNumber = as.integer(strsplit(step$sender,"V")[[1]][2])
     receiverNumber = as.integer(strsplit(step$receiver,"V")[[1]][2])

     newp = (newM-newMe[senderNumber,receiverNumber]+1)/newGammaEplus[senderNumber,receiverNumber]


     if(step$typeS == "same"){
       pUndoStep = (newGammaEplus[senderNumber,receiverNubmer]/newGammaPlus) * newp * (1/(newM-newMe[senderNumber,receiverNumber]+1))

     }else if(step$typeS == "diff"){
       pUndoStep = (newGammaEplus[senderNumber,receiverNumber]/newGammaPlus) * (1-newp) * (1/(newM-newMe[senderNumber,receiverNumber]+1))*(1/(newM-newMe[senderNumber,receiverNumber]))
     }

     loglikSeq = logLikelihoodMC(1,seq,beta,list(1),actDfnodes=actDfnodes,net0=net0,formula=formula)
     newloglikSeq = logLikelihoodMC(indexCore=1,newseq,beta,splitIndicesPerCore=list(1),actDfnodes=actDfnodes,net0=net0,formula=formula)


   }else if(type == 3){ # Permutation
      step = stepPerm(seq,tieNames,m,me)
      getNewKelMeMatrix = getKelMeMatrix(step$newseq,actDfnodesLab)
      newAuxDfE = getAuxDfE(getNewKelMeMatrix$auxDf, sender,receiver)
      auxDfE = getAuxDfE(auxDf,sender,receiver)

      newKel_g1 = getNewKelMeMatrix$Kel_g1
      newKel_ge1 = getNewKelMeMatrix$Kel_ge1
      newGammaEminus = choose(newKel_ge1,2) + newKel_g1
      newGammaMinus = sum(newGammaEminus)
      newMe = getNewKelMeMatrix$me
      newGammaEplus = choose(nrow(newseq)-newMe+2,2)
      newGammaPlus = sum(newGammaEplus)
      newM = nrow(newseq)

      sender1=as.integer(strsplit(step$sender[1],"V")[[1]][2])
      sender2=as.integer(strsplit(step$sender[2],"V")[[1]][2])
      receiver1=as.integer(strsplit(step$receiver[1],"V")[[1]][2])
      receiver2=as.integer(strsplit(step$receiver[2],"V")[[1]][2])

      probVec = newMe[sender2,receiver2]/newM
      probVec2 = newMe[sender1,receiver1]/(newM-newMe[sender2,receiver2])
      probVec3 = newMe[sender1,receiver1]/newM
      probVec4 = newMe[sender2,receiver2]/(newM-newMe[sender1,receiver1])

      pUndoStep = (probVec * probVec2 + robVec3 * probVec4 )*
        (1/nrow(newseq[newseq$sender==paste("V",sender1,sep="") & newseq$receiver==paste("V",receiver1,sep=""),])) *
        (1/nrow(newseq[newseq$sender==paste("V",sender2,sep="") & newseq$receiver==paste("V",receiver2,sep=""),]))
      loglikSeq = logLikelihoodMC(1,seq,beta,list(1),actDfnodes=actDfnodes,net0=net0,formula=formula)
      newloglikSeq = logLikelihoodMC(indexCore=1,newseq,beta,splitIndicesPerCore=list(1),actDfnodes=actDfnodes,net0=net0,formula=formula)
    }


    u = runif(1,min=0,max=1)
    accept =(u<= exp(newloglikSeq-loglikSeq)*pUndoStep/step$pDoStep)
    if(accept){
      accecptIndex = accecptIndex + 1
      #if(accept & !acceptIndex %% 50){
      #Update sequence every 50 steps
      seq = newseq
      getKelMeMatrix = newGetKelMeMatrix
      Kel_g1 = gnewKel_g1
      Kel_ge1 = newKel_ge1
      gammaEminus = newGammaEminus
      gammaMinus = newGammaMinus
      me = newMe
      gammaEplus = newGammaEplus
      gammaPlus = newGammaPlus
      m = newM
      auxDf = getKelMeMatrix$auxDf
    }
  }

}

MCMC_MC = function(indexCore,splitIndicesPerCore,permut = permut,beta=beta, burnIn = TRUE, H = H,actDfnodes=actDfnodes){

  indicesCore = splitIndicesPerCore[[indexCore]]
  resMCMC = vector("list",length(indicesCore))

  if(burnIn) {burn_in = burnIn(permut[[1]],beta,H,actDfnodes)}


  for(i in seq_along(indicesCore)){
    resMCMC[[i]] = MCMC(permut[[indicesCore[[i]]]],burn_in,H,actDfnodes)
  }

  return(resMCMC)
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


MCEMalgorithm = function(nmax,net0,net1,theta0,beta0,formula,num_cores=1){
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
  H = length(sequence) # Humming distance
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
                         "parameters","timeGenerator","logLikelihood","MCMC","burnIn",
                        "getKelMatrix","getmeMatrix")) # to be updated

  permut = clusterApply(cl,seq_along(splitIndicesPerCore),MCMC_MC,permut = permut,
                        beta=beta,splitIndicesPerCore = splitIndicesPerCore, burnIn = TRUE,
                        H = H,actDfnodes=actDfnodes) # to be updated


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
      newpermut = permute(sequence, nmax = nmax/2)
      newpermut = clusterApply(cl,seq_along(splitIndicesPerCore),MCMC_MC,permut = newpermut,
                            beta=beta,splitIndicesPerCore=splitIndicesPerCore, burnIn = FALSE,
                            H = H,actDfnodes=actDfnodes)
      nmax = nmax+(nmax/2)
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
        newpermut = permute(sequence, nmax = m_start - nmax)
        nmax = m_start
        permut = c(permut,newpermut)
      }

       permut = clusterApply(cl,seq_along(splitIndicesPerCore),MCMC_MC,permut = permut,
                             beta=beta, splitIndicesPerCore=splitIndicesPerCore, burnIn = TRUE,
                             H = H,actDfnodes=actDfnodes)


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




