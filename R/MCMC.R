# MCMC auxiliar functions --------------------


#' Get auxiliary matrix \eqn{K_{el}} and \eqn{m_e} for MCMC steps.
#'
#' @param seq data frame with the sequence of events up to time \eqn{t}.
#' @param actDfnodes vector with the labels of the actors of the network.
#'
#' @return list with matrix \eqn{K_{el}}, \eqn{K_{e,l>1}}, \eqn{m_e} and auxiliary data frame
#' @export
#'
#' @examples seq = data.frame("sender" = c(1,2,3), "receiver" = c(3,2,1), "replace" = c(0,1,1))
#' getKelMeMatrix(seq,actDfnodes = c(1,2,3))
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

