#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP regsem_rcpp_fit_fun(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP regsem_rcpp_grad_ram(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP regsem_rcpp_RAMmult(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"regsem_rcpp_fit_fun",  (DL_FUNC) &regsem_rcpp_fit_fun,   8},
    {"regsem_rcpp_grad_ram", (DL_FUNC) &regsem_rcpp_grad_ram, 12},
    {"regsem_rcpp_RAMmult",  (DL_FUNC) &regsem_rcpp_RAMmult,   9},
    {NULL, NULL, 0}
};

void R_init_regsem(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
