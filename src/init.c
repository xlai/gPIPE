#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP gpipe_mincut(SEXP lower_, SEXP upper_, SEXP d_, SEXP tol_, SEXP largest_);

static const R_CallMethodDef call_methods[] = {
  {"gpipe_mincut", (DL_FUNC) &gpipe_mincut, 5},
  {NULL, NULL, 0}
};

void R_init_gPIPE(DllInfo *dll) {
  R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
