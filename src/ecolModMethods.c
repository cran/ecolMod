#include<R_ext/Rdynload.h>
#ifndef R_R_H
# include <R.h>
#endif

#ifndef R_EXT_DYNLOAD_H_
# include <R_ext/Rdynload.h>
#endif

#include <Rinternals.h>
#include <stdlib.h> // for NULL

/* .Fortran calls */
void F77_NAME(initomexdia)(void (* steadyparms)(int *, double *));
void F77_NAME(omexdiamod)(int *, double *, double *, double *, double *, int *);
void F77_NAME(lattice)(int *, int *, int *, int *, double *, int *, int *);

R_FortranMethodDef FEntries[] = {
    {"initomexdia",    (DL_FUNC) &F77_SUB(initomexdia),     1},
    {"omexdiamod",     (DL_FUNC) &F77_SUB(omexdiamod),      4},
    {"lattice",        (DL_FUNC) &F77_SUB(lattice),         7},
    {NULL, NULL, 0}
};

/* Initialization ---------------------------------------------------------- */
void R_init_ecolMod(DllInfo *dll) {

  // register entry points
  R_registerRoutines(dll, NULL, NULL, FEntries, NULL);

  R_useDynamicSymbols(dll, TRUE); // enable dynamic searching
} 
