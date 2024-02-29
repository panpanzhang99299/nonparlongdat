#include <R.h>
#include <Rinternals.h> // for SEXP
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _nonparlongdat_logdmvnorm(SEXP X, SEXP mean, SEXP Sigma);
extern SEXP _nonparlongdat_ideniforj(SEXP Xi, SEXP Xj);
extern SEXP _nonparlongdat_idensigmaforj(SEXP Xi, SEXP Sigma); 
extern SEXP _nonparlongdat_loglikvalid(SEXP Y, SEXP X, SEXP Z, SEXP theta, SEXP Sigma);
extern SEXP _nonparlongdat_logliknonvalidXinvauxinv(SEXP Yval, SEXP Ynonval, SEXP Xval, SEXP Xnonval, SEXP Z, SEXP auxval, SEXP auxnoval, SEXP theta, SEXP Sigma, SEXP H, SEXP auxcont); 
extern SEXP _nonparlongdat_logliknonvalidXvaryauxinv(SEXP Yval, SEXP Ynonval, SEXP Xval, SEXP Xnonval, SEXP Z, SEXP auxval, SEXP auxnoval, SEXP theta, SEXP Sigma, SEXP H);
extern SEXP _nonparlongdat_logliknonvalidXvaryauxvary(SEXP Yval, SEXP Ynonval, SEXP Xval, SEXP Xnonval, SEXP Z, SEXP auxval, SEXP auxnoval, SEXP theta, SEXP Sigma, SEXP H, SEXP auxcont); 

static const R_CallMethodDef R_CallDef[] = {
    {"_nonparlongdat_logdmvnorm",                 (DL_FUNC) &_nonparlongdat_logdmvnorm,                   3},
    {"_nonparlongdat_ideniforj",                  (DL_FUNC) &_nonparlongdat_ideniforj,                    2},
    {"_nonparlongdat_idensigmaforj",              (DL_FUNC) &_nonparlongdat_idensigmaforj,                2},
    {"_nonparlongdat_loglikvalid",                (DL_FUNC) &_nonparlongdat_loglikvalid,                  5},
    {"_nonparlongdat_logliknonvalidXinvauxinv",   (DL_FUNC) &_nonparlongdat_logliknonvalidXinvauxinv,    11},
    {"_nonparlongdat_logliknonvalidXvaryauxinv",  (DL_FUNC) &_nonparlongdat_logliknonvalidXvaryauxinv,   10},
    {"_nonparlongdat_logliknonvalidXvaryauxvary", (DL_FUNC) &_nonparlongdat_logliknonvalidXvaryauxvary,  11},
    {NULL, NULL, 0}
};

void R_init_nonparlongdat(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, R_CallDef, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
