#ifndef QNORM_H
#define QNORM_H 1

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
 


int qnorm_c(double *data, int *rows, int *cols);
int qnorm_robust_c(double *data,double *weights, int *rows, int *cols, int *use_median, int *use_log2, int *weight_scheme);
int qnorm_c_using_target(double *data, int *rows, int *cols, double *target, int *targetrows);
int qnorm_c_determine_target(double *data, int *rows, int *cols, double *target, int *targetrows);
int qnorm_c_within_blocks(double *x, int *rows, int *cols, int *blocks);


int qnorm_c_l(double *data, size_t rows, size_t cols);
int qnorm_c_using_target_l(double *data, size_t rows, size_t cols, double *target, size_t targetrows);
int qnorm_c_determine_target_l(double *data, size_t rows, size_t cols, double *target, size_t targetrows);




SEXP R_qnorm_c(SEXP X, SEXP copy);

SEXP R_qnorm_robust_weights(SEXP X, SEXP remove_extreme, SEXP n_remove);
SEXP R_qnorm_robust_c(SEXP X, SEXP copy, SEXP R_weights, SEXP R_use_median, SEXP R_use_log2, SEXP R_weight_scheme);

SEXP R_qnorm_determine_target(SEXP X, SEXP targetlength);
SEXP R_qnorm_using_target(SEXP X, SEXP target,SEXP copy);
SEXP R_qnorm_within_blocks(SEXP X,SEXP blocks,SEXP copy);

SEXP R_qnorm_c_handleNA(SEXP X, SEXP copy);
SEXP R_qnorm_determine_target_via_subset(SEXP X, SEXP subset, SEXP targetlength);
SEXP R_qnorm_using_target_via_subset(SEXP X,  SEXP subset, SEXP target,SEXP copy);

#endif

