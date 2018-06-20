// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// readBCFQuery_
List readBCFQuery_(SEXP fname, SEXP reg);
RcppExport SEXP _boringLD_readBCFQuery_(SEXP fnameSEXP, SEXP regSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type fname(fnameSEXP);
    Rcpp::traits::input_parameter< SEXP >::type reg(regSEXP);
    rcpp_result_gen = Rcpp::wrap(readBCFQuery_(fname, reg));
    return rcpp_result_gen;
END_RCPP
}
// readBCFContigs_
List readBCFContigs_(SEXP fname);
RcppExport SEXP _boringLD_readBCFContigs_(SEXP fnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type fname(fnameSEXP);
    rcpp_result_gen = Rcpp::wrap(readBCFContigs_(fname));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_boringLD_readBCFQuery_", (DL_FUNC) &_boringLD_readBCFQuery_, 2},
    {"_boringLD_readBCFContigs_", (DL_FUNC) &_boringLD_readBCFContigs_, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_boringLD(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}