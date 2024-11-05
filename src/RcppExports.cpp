// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// timesTwo
NumericVector timesTwo(NumericVector& x);
RcppExport SEXP _emDyNAM_timesTwo(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(timesTwo(x));
    return rcpp_result_gen;
END_RCPP
}
// subsetDataFrameTwoConditions
DataFrame subsetDataFrameTwoConditions(DataFrame df, String columnName1, int value1, String columnName2, int value2);
RcppExport SEXP _emDyNAM_subsetDataFrameTwoConditions(SEXP dfSEXP, SEXP columnName1SEXP, SEXP value1SEXP, SEXP columnName2SEXP, SEXP value2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    Rcpp::traits::input_parameter< String >::type columnName1(columnName1SEXP);
    Rcpp::traits::input_parameter< int >::type value1(value1SEXP);
    Rcpp::traits::input_parameter< String >::type columnName2(columnName2SEXP);
    Rcpp::traits::input_parameter< int >::type value2(value2SEXP);
    rcpp_result_gen = Rcpp::wrap(subsetDataFrameTwoConditions(df, columnName1, value1, columnName2, value2));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_emDyNAM_timesTwo", (DL_FUNC) &_emDyNAM_timesTwo, 1},
    {"_emDyNAM_subsetDataFrameTwoConditions", (DL_FUNC) &_emDyNAM_subsetDataFrameTwoConditions, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_emDyNAM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
