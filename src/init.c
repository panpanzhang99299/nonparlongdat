#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _nonparlongdat_logdmvnorm(SEXP, SEXP, SEXP);
extern SEXP _nonparlongdat_ideniforj(SEXP, SEXP);
extern SEXP _nonparlongdat_idensigmaforj(SEXP, SEXP); 
extern SEXP _nonparlongdat_loglikvalid(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _nonparlongdat_logliknonvalidXinvauxinv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 
extern SEXP _nonparlongdat_logliknonvalidXvaryauxinv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _nonparlongdat_logliknonvalidXvaryauxvary(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP); 

static const R_CallMethodDef CallEntries[] = {
    {"_nonparlongdat_logdmvnorm",                 (DL_FUNC) &_nonparlongdat_logdmvnorm,                   3},
    {"_nonparlongdat_ideniforj",                  (DL_FUNC) &_nonparlongdat_ideniforj,                    2},
    {"_nonparlongdat_idensigmaforj",              (DL_FUNC) &_nonparlongdat_idensigmaforj,                2},
    {"_nonparlongdat_loglikvalid",                (DL_FUNC) &_nonparlongdat_loglikvalid,                  5},
    {"_nonparlongdat_logliknonvalidXinvauxinv",   (DL_FUNC) &_nonparlongdat_logliknonvalidXinvauxinv,    11},
    {"_nonparlongdat_logliknonvalidXvaryauxinv",  (DL_FUNC) &_nonparlongdat_logliknonvalidXvaryauxinv,   10},
    {"_nonparlongdat_logliknonvalidXvaryauxvary", (DL_FUNC) &_nonparlongdat_logliknonvalidXvaryauxvary,  11},
    {NULL, NULL, 0}
};

void R_init_wdnet(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
