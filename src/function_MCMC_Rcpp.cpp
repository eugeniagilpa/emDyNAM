// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <string>
#include <iostream>
using namespace Rcpp;


// [[Rcpp::export]]
double getpDoAugment(NumericMatrix gammaEplus, double gammaPlus, double p, double m, NumericMatrix me, int sender, int receiver,
                     std::string typeA, IntegerVector indexSR, double pAug) {
  double pDoStep;
  if (typeA == "same") {
    if (indexSR.length() < 1) {
      pDoStep = (gammaEplus(sender, receiver) / gammaPlus) * p * (1 / (m + 1)) * pAug;
    } else {
      pDoStep = (gammaEplus(sender, receiver) / gammaPlus) * p * (1 / (m - me(sender, receiver) + 1)) * pAug;
    }
  } else {
    pDoStep = (gammaEplus(sender, receiver) / gammaPlus) * (1 - p) * (1 / (m - me(sender, receiver) + 1)) * (1 / (m - me(sender, receiver))) * pAug;
  }
  return pDoStep;
}

// [[Rcpp::export]]
double getpDoShort(NumericMatrix gammaEminus, double gammaMinus, double p, double uniqueR, int sender, int receiver,
                   std::string typeS, double pShort) {
  double pDoStep;
  if (typeS == "same") {
    pDoStep = (gammaEminus(sender, receiver) / gammaMinus) * p * (1 / uniqueR) * pShort;
  } else {
    pDoStep = (gammaEminus(sender, receiver) / gammaMinus) * (1 - p) * (1 / uniqueR) * (1 / (uniqueR - 1)) * pShort;
  }
  return pDoStep;
}

// [[Rcpp::export]]
double getpDoPerm(NumericVector probVec, NumericVector probVec2, int e1, int e2, int sender1, int receiver1, int sender2, int receiver2,
                  NumericMatrix me) {
  double pDoStep = (probVec[e1] * probVec2[e2] + probVec[e2] * probVec2[e1]) * (1/me(sender1,receiver1)) * (1/me(sender2,receiver2));
  return pDoStep;
}


// [[Rcpp::export]]
DataFrame getAuxDfE(List auxDf, int sender, int receiver) {
  DataFrame auxDfE = auxDf[sender][receiver];
  IntegerVector indexNo1 = which(auxDfE["rowDiff"] != 1);
  auxDfE["run"] = rep(0, auxDfE.nrows());
  int run = 1;
  if (indexNo1.size() > 1) {
    for (int i = 0; i < (indexNo1.size() - 1); ++i) {
      auxDfE["run"][Range(indexNo1[i], indexNo1[i + 1])] = run;
      run += 1;
    }
    auxDfE["run"][Range(indexNo1[indexNo1.size() - 1], auxDfE.nrows())] = run;
  } else if (indexNo1.size() == 1) {
    auxDfE["run"][Range(indexNo1[0], auxDfE.nrows())] = 1;
  }
  return auxDfE;
}



// [[Rcpp::export]]
List getKelMeMatrix(DataFrame seq, IntegerVector actDfnodesLab) {
  if (!seq.containsElementNamed("row")) {
    seq.push_back(seq.nrows(), "row");
  }

  List auxDf(actDfnodesLab.size());
  for (int x : actDfnodesLab) {
    List innerList(actDfnodesLab.size());
    for (int i : actDfnodesLab) {
      DataFrame subsetSeq = subset(seq, seq["sender"] == x & seq["receiver"] == i);
      innerList[i] = subsetSeq;
    }
    auxDf[x] = innerList;
  }

  for (int i = 0; i < actDfnodesLab.size(); ++i) {
    for (int j = 0; j < actDfnodesLab.size(); ++j) {
      DataFrame df = auxDf[i][j];
      if (df.nrows() > 0) {
        df["row"] = as<IntegerVector>(df["row"]);
        IntegerVector rowDiff = diff(as<IntegerVector>(df["row"]));
        rowDiff.push_front(0);
        df["rowDiff"] = rowDiff;
        auxDf[i][j] = df;
      }
    }
  }

  NumericMatrix Kel_g1(actDfnodesLab.size(), actDfnodesLab.size());
  NumericMatrix Kel_ge1(actDfnodesLab.size(), actDfnodesLab.size());

  for (int i = 0; i < Kel_g1.nrow(); ++i) {
    for (int j = 0; j < Kel_g1.ncol(); ++j) {
      DataFrame auxDfE = getAuxDfE(auxDf, i, j);
      IntegerVector runTable = table(auxDfE["run"]);
      Kel_ge1(i, j) = runTable.size();
      Kel_g1(i, j) = sum(runTable > 1);
    }
  }

  NumericMatrix me = transpose(sapply(auxDf, [](List x) { return sapply(x, nrows); }));

  return List::create(Named("Kel_g1") = Kel_g1, Named("Kel_ge1") = Kel_ge1, Named("me") = me, Named("auxDf") = auxDf);
}



// [[Rcpp::export]]
List stepAugment(DataFrame seq, CharacterVector tieNames, NumericMatrix gammaEplus, double gammaPlus,
                 double m, NumericMatrix me, NumericMatrix net0, double pAug) {
  // Choose element to be inserted
  String e = as<String>(tieNames[Rcpp::RcppArmadillo::sample(tieNames.size(), 1, false, gammaEplus / gammaPlus)(0)]);
  IntegerVector e_split = as<IntegerVector>(strsplit(e, "-"));
  int sender = e_split[0];
  int receiver = e_split[1];

  double p = (m - me(sender, receiver) + 1) / gammaEplus(sender, receiver);
  String typeA = as<String>(Rcpp::RcppArmadillo::sample(2, 1, false, NumericVector::create(p, 1 - p))(0) == 0 ? "same" : "diff");
  IntegerVector indexSR = which(seq["sender"] == sender & seq["receiver"] == receiver);

  List newseq;
  double pDoStep;
  if (typeA == "same") {
    int place;
    if (indexSR.size() < 1) {
      place = Rcpp::RcppArmadillo::sample(m + 1, 1, false)(0);
    } else {
      place = Rcpp::RcppArmadillo::sample(seq(m + 1)[seq(m + 1) != indexSR], 1, false)(0);
    }

    // Create new sequence
    // ... (Handle insertion logic similar to R)

    pDoStep = getpDoAugment(gammaEplus, gammaPlus, p, m, me, sender, receiver, typeA, indexSR, pAug);
  } else if (typeA == "diff") {
    // Handle "diff" case
    // ... (Handle insertion logic similar to R)

    pDoStep = getpDoAugment(gammaEplus, gammaPlus, p, m, me, sender, receiver, typeA, indexSR, pAug);
  }

  return List::create(Named("newseq") = newseq, Named("sender") = sender, Named("receiver") = receiver,
                            Named("place") = place, Named("typeA") = typeA, Named("pDoStep") = pDoStep);
}



// [[Rcpp::export]]
List stepShort(DataFrame seq, CharacterVector tieNames, NumericMatrix gammaEminus, double gammaMinus,
               double m, NumericMatrix me, NumericMatrix Kel_g1, DataFrame auxDf, double pShort) {
  // Choose element to be deleted
  String e = as<String>(tieNames[Rcpp::RcppArmadillo::sample(tieNames.size(), 1, false, gammaEminus / gammaMinus)(0)]);
  IntegerVector e_split = as<IntegerVector>(strsplit(e, "-"));
  int sender = e_split[0];
  int receiver = e_split[1];

  double p = Kel_g1(sender, receiver) / gammaEminus(sender, receiver);
  String typeS = as<String>(Rcpp::RcppArmadillo::sample(2, 1, false, NumericVector::create(p, 1 - p))(0) == 0 ? "same" : "diff");

  List newseq;
  double pDoStep;
  if (typeS == "same") {
    // Handle "same" case
    // ... (Handle deletion logic similar to R)

    pDoStep = getpDoShort(gammaEminus, gammaMinus, p, length(which(table(auxDf["run"]) > 1)), sender, receiver, typeS, pShort);
  } else if (typeS == "diff") {
    // Handle "diff" case
    // ... (Handle deletion logic similar to R)

    pDoStep = getpDoShort(gammaEminus, gammaMinus, p, length(unique(auxDf["run"])), sender, receiver, typeS, pShort);
  }

  return List::create(Named("newseq") = newseq, Named("sender") = sender, Named("receiver") = receiver,
                            Named("place") = place, Named("typeS") = typeS, Named("pDoStep") = pDoStep);
}


// [[Rcpp::export]]
List stepPerm(DataFrame seq, CharacterVector tieNames, double m, NumericMatrix me) {
  NumericVector probVec = as<NumericVector>(me / m);
  NumericVector probVec2 = as<NumericVector>(me / (m - me));
  probVec.names() = tieNames;
  probVec2.names() = tieNames;

  String e1 = as<String>(tieNames[Rcpp::RcppArmadillo::sample(tieNames.size(), 1, false, probVec)(0)]);
  String e2 = as<String>(tieNames[tieNames != e1][Rcpp::RcppArmadillo::sample(tieNames[tieNames != e1].size(), 1, false, probVec2[names(probVec) != e1])(0)]);

  IntegerVector e1_split = as<IntegerVector>(strsplit(e1, "-"));
  int sender1 = e1_split[0];
  int receiver1 = e1_split[1];
  int place1 = as<int>(Rcpp::RcppArmadillo::sample(as<NumericVector>(seq[seq["sender"] == sender1 & seq["receiver"] == receiver1]["row"]), 1, false)(0));

  IntegerVector e2_split = as<IntegerVector>(strsplit(e2, "-"));
  int sender2 = e2_split[0];
  int receiver2 = e2_split[1];
  int place2 = as<int>(Rcpp::RcppArmadillo::sample(as<NumericVector>(seq[seq["sender"] == sender2 & seq["receiver"] == receiver2]["row"]), 1, false)(0));

  DataFrame newseq = clone(seq);
  DataFrame auxRow = newseq.row(place1);
  newseq.row(place1) = newseq.row(place2);
  newseq.row(place2) = auxRow;

  // Update "replace" column
  // ... (Handle update logic similar to R)

  double pDoStep = getpDoPerm(probVec, probVec2, e1, e2, sender1, receiver1, sender2, receiver2, me);

  return List::create(Named("newseq") = newseq, Named("sender") = IntegerVector::create(sender1, sender2),
                      Named("receiver") = IntegerVector::create(receiver1, receiver2),
                      Named("place") = IntegerVector::create(place1, place2), Named("pDoStep") = pDoStep);
}


