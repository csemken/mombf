#ifndef NLPMARG_H
#define NLPMARG_H 1

//#include <R.h>
//#include <Rinternals.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <list>
#include "crossprodmat.h"
#include "Polynomial.h"


/*
 * Function Prototypes
 */
extern "C" {
  SEXP nlpMarginalCI(SEXP Ssel, SEXP Snsel, SEXP Sfamily, SEXP SpriorCoef, SEXP SpriorGroup, SEXP Sn, SEXP Sp, SEXP Sy, SEXP Suncens, SEXP Ssumy2, SEXP Sx, SEXP SXtX, SEXP SytX, SEXP Smethod, SEXP Shesstype, SEXP SoptimMethod, SEXP SB, SEXP Salpha, SEXP Slambda, SEXP Stau, SEXP Staugroup, SEXP Staualpha, SEXP Sfixatanhalpha, SEXP Sr, SEXP Sgroups, SEXP Sngroups, SEXP Snvaringroup, SEXP Sconstraints, SEXP Sinvconstraints, SEXP Slogscale);
}

double nlpMarginal(int *sel, int *nsel, int *family, int *prCoef, int *prGroup, struct marginalPars *pars);

#endif /* NLPMARG_H */
