// Copyright (c) 2021 Ozette Technologies
/*
Header file added to flowCore by John Ramey
The `R_zeroin` function used to be and is no longer apart of the R API.
See the following link for more information:
http://developer.r-project.org/blosxom.cgi/R-devel/2012/08/15#n2012-08-15
*/
#ifndef ZEROIN_H_
#define ZEROIN_H_

#include <R_ext/Boolean.h>
#include <R_ext/RS.h>           /* F77_... */
#include <R_ext/BLAS.h>
#include <R.h>
#include <Rinternals.h>




double  R_zeroin(double ax, double bx, double (*f)(double, void *), void *info, double *Tol, int *Maxit);
double R_zeroin2(double ax,double bx,double fa, double fb,double (*f)(double x, void *info),void *info,double *Tol,int *Maxit);
#endif /* ZEROIN_H_ */
