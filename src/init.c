#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* prototypes for functions implemented in C++ (use C linkage) */
extern SEXP do_vegdist(SEXP, SEXP);
extern SEXP do_vegdist_n(SEXP, SEXP, SEXP, SEXP);
extern SEXP do_vegdist_p(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"do_vegdist", (DL_FUNC) &do_vegdist, 2},
  {"do_vegdist_n",      (DL_FUNC) &do_vegdist_n,  4},
  {"do_vegdist_p",      (DL_FUNC) &do_vegdist_p,  3},
  {NULL, NULL, 0}
};

void R_init_frescaloDT(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
