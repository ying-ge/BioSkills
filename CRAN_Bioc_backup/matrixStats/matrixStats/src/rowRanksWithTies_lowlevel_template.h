/***********************************************************************
 TEMPLATE:
  <col|row>Ranks_dbl_ties<Min|Max|Average>(ARGUMENTS_LIST)

 ARGUMENTS_LIST:
  X_C_TYPE *x, R_xlen_t nrow, R_xlen_t ncol, R_xlen_t *rows, R_xlen_t nrows, int rowsHasNA, R_xlen_t *cols, R_xlen_t ncols, int colsHasNA, ANS_C_TYPE *ans

 Arguments:
   The following macros ("arguments") should be defined for the
   template to work as intended.

  - METHOD_NAME: the name of the resulting function
  - X_TYPE: 'i' or 'r'
  - ANS_TYPE: 'i' or 'r'
  - TIESMETHOD: 'a' (average), 'f' (first), 'l' (last), 'r' (random), '0' (min), '1' (max), 'd' (dense)

 Authors:
  Hector Corrada Bravo [HCB]
  Peter Langfelder [PL]
  Henrik Bengtsson [HB]
  Brian Montgomery 
  Jakob Peder Pettersen [JPP]
 ***********************************************************************/
#include <Rinternals.h>

#undef RANK
#if TIESMETHOD == 'a' /* average */
  #define ANS_TYPE 'r'
  #define RANK(firstTie, aboveTie) ((double) (firstTie + aboveTie + 1))/2
#elif TIESMETHOD == '0' /* min */
  #define ANS_TYPE 'i'
  #define RANK(firstTie, aboveTie) firstTie + 1
#elif TIESMETHOD == '1' /* max */
  #define ANS_TYPE 'i'
  #define RANK(firstTie, aboveTie) aboveTie
#else
  #define ANS_TYPE 'i' /* dense and other(RANK not used) */
  #define RANK(firstTie, aboveTie) firstTie + 1
#endif

/* Expand arguments:
    X_TYPE => (X_C_TYPE, X_IN_C, X_ISNAN)
    ANS_TYPE => (ANS_SXP, ANS_NA, ANS_C_TYPE, ANS_IN_C)
 */
#include "000.templates-types.h"

void SHUFFLE_INT(int *array, size_t i, size_t j); /* prototype for use with "random" */

/* Indexing formula to compute the vector index of element j of vector i.
   Should take arguments element, vector, nElements, nVectors. */
#undef ANS_INDEX_OF


void CONCAT_MACROS(METHOD, X_C_SIGNATURE)(X_C_TYPE *x, R_xlen_t nrow, R_xlen_t ncol, 
                 R_xlen_t *rows, R_xlen_t nrows, int rowsHasNA,
                 R_xlen_t *cols, R_xlen_t ncols, int colsHasNA,
                 int byrow, ANS_C_TYPE *ans) {
  ANS_C_TYPE rank;
  X_C_TYPE *values, current, tmp;
  R_xlen_t *colOffset;
  R_xlen_t ii, jj, kk, rowIdx, idx;
  int *I;
  int lastFinite, firstTie, aboveTie, dense_rank_adj;
  int nvalues, nVec;
  int norows, nocols;
  if (cols == NULL) { nocols = 1; } else { nocols = 0; }
  if (rows == NULL) { norows = 1; } else { norows = 0; }

  if (byrow) {
    nvalues = ncols;
    nVec = nrows;
  
    /* Pre-calculate the column offsets */
    colOffset = (R_xlen_t *) R_alloc(ncols, sizeof(R_xlen_t));
      for (jj=0; jj < ncols; jj++) {
        if (nocols) {
        colOffset[jj] = jj * nrow;
        } else {
          if (!colsHasNA) {
            colOffset[jj] = cols[jj] * nrow;
          } else {
            colOffset[jj] = R_INDEX_OP(cols[jj], *, nrow, 1, 0);
          } 
        }
      }

  } else {
    nvalues = nrows;
    nVec = ncols;
  
    /* Pre-calculate the column offsets */
    colOffset = (R_xlen_t *) R_alloc(nrows, sizeof(R_xlen_t));
    
    for (jj=0; jj < nrows; jj++) {
      if (norows) {
        colOffset[jj] = jj;
      } else {
        colOffset[jj] = rows[jj];
      }
    }
  } // if (byrow)
    
  

  values = (X_C_TYPE *) R_alloc(nvalues, sizeof(X_C_TYPE));
  I = (int *) R_alloc(nvalues, sizeof(int));

  for (ii=0; ii < nVec; ii++) {
    if (byrow) {
      if (norows) {
        rowIdx = ii;
      } else {
        rowIdx = rows[ii];
      }
    } else {
      if (nocols) {
        rowIdx = ii * nrow;
      } else {
        if (!colsHasNA) {
          rowIdx = cols[ii] * nrow;
        } else {
          rowIdx = R_INDEX_OP(cols[ii], *, nrow, 1, 0);
        }
      }
    }
    lastFinite = nvalues-1;

    /* Put the NA/NaN elements at the end of the vector and update
       the index vector appropriately.
       This may be a bit faster since it only uses one loop over
       the length of the vector, plus it shortens the sort in case
       there are missing values. /PL (2012-12-14)
    */
    for (jj = 0; jj <= lastFinite; jj++) {
      /*
       * Checking for the colsHasNA when we already have to check colsHasNA || rowsHasNA
       * is indeed useless, but for keeping the code idiomatic, we still do it
       * Hopefully, the compiler will optimize out the unnecessary instructions [JPP].
       */
        if (!rowsHasNA && !colsHasNA) {
          idx = rowIdx + colOffset[jj];
          tmp = x[idx];
        } else {
          idx = R_INDEX_OP(rowIdx, +, colOffset[jj], 1, 1);
          tmp = R_INDEX_GET(x, idx, X_NA, 1);
        }
      
      if (X_ISNAN(tmp)) {
        R_xlen_t lastFinite_idx;
        X_C_TYPE lastFinite_val;
        while (lastFinite > jj) {
          if (!rowsHasNA && !colsHasNA) {
            lastFinite_idx = rowIdx + colOffset[lastFinite];
            lastFinite_val = x[lastFinite_idx];  
          } else {
            lastFinite_idx = R_INDEX_OP(rowIdx, +, colOffset[lastFinite], 1, 1);
            lastFinite_val = R_INDEX_GET(x, lastFinite_idx, X_NA, 1);  
          }
          
          if (!X_ISNAN(lastFinite_val)) {
            break;
          }
          I[lastFinite] = lastFinite;
          lastFinite--;
        }

        I[lastFinite] = jj;
        I[jj] = lastFinite;
        values[ jj ] = lastFinite_val;
        values[ lastFinite ] = tmp;
        lastFinite--;
      } else {
        I[jj] = jj;
        values[ jj ] = tmp;
      }
    } /* for (jj ...) */

   // Diagnostic print-outs
   
    /*
    Rprintf("Swapped vector:\n");
    for (jj=0; jj < nvalues; jj++)
    {
      Rprintf(" %8.4f,", values[jj]);
      if (((jj+1) % 5==0) || (jj==nvalues-1)) Rprintf("\n");
    }
    Rprintf("Index vector:\n");
    for (jj=0; jj<nvalues; jj++)
    {
      Rprintf(" %d,", I[jj]);
      if (((jj+1) % 5==0) || (jj==nvalues-1)) Rprintf("\n");
    }
    */


    // This will sort the data in increasing order and use the I vector to keep track of the original
    // indices. it only makes sense to do sort if there are at least 2 finite values.
    //
    if (lastFinite > 0) X_QSORT_I(values, I, 1, lastFinite + 1);

    // Calculate the ranks.
    firstTie = 0;
    aboveTie = 1;
    dense_rank_adj = 0;
    
    for (jj=0; jj <= lastFinite;) {
      
      if (TIESMETHOD == 'd') {
        dense_rank_adj += (aboveTie - firstTie - 1);
        firstTie = jj - dense_rank_adj;
      } else {
        firstTie = jj;
      }
      
      current = values[jj];
      
      if (X_ISNAN(current)) {
        /*
         * This is really a runtime check of an internal programming error. Preferentially, it should have been
         * caught by testing. The only problem is that the test would be stuck in an infinite loop unless it is
         * checked for at runtime. 
         */
        error("Internal matrixStats programming error, NaN values not handled correctly");
      }
      
      while ((jj <= lastFinite) && (values[jj] == current)) jj++;
      
      if (TIESMETHOD == 'd') {
        aboveTie = jj - dense_rank_adj;
      } else {
        aboveTie = jj;
      }
      // X_QSORT_I is not stable - ties can be permuted.
      // This restores the original order.
      // It might be more efficient to use a stable sort to begin with.
      if (TIESMETHOD == 'f' || TIESMETHOD == 'l') {
        R_qsort_int(I, firstTie + 1, aboveTie); /* Function is 1-based */
      // SHUFFLE_INT randomizes the order.
      } else if (TIESMETHOD == 'r') {
        SHUFFLE_INT(I, firstTie, aboveTie - 1);
      } else {
        // Get appropriate rank for average, min, max, or dense
        rank = RANK(firstTie, aboveTie);
      }
      
      for (kk=firstTie; kk < aboveTie; kk++) {
        if (byrow) {
          switch (TIESMETHOD) {
          /*
           * Note the fallthrough
           */
          case 'f':
          case 'r':
            ans[ii + I[kk]*nrows] = kk + 1;
            break;
          case 'l':
            ans[ii + I[kk]*nrows] = aboveTie - (kk - firstTie);
            break;
          case 'd':
            ans[ii + I[kk + dense_rank_adj]*nrows] = rank;
            break;
          default:
            ans[ii + I[kk]*nrows] = rank;
            break;
          }
        } else {
          switch (TIESMETHOD) {
          case 'f':
          case 'r':
            ans[I[kk] + ii*nrows] = kk + 1;
            break;
          case 'l':
            ans[I[kk] + ii*nrows] = aboveTie - (kk - firstTie);
            break;
          case 'd':
            ans[I[kk + dense_rank_adj] + ii*nrows] = rank;
            break;
          default:
            ans[I[kk] + ii*nrows] = rank;
            break;
          }
        }
      }
    }
    
    
    // At this point jj = lastFinite + 1, no need to re-initialize again.
    for (; jj < nvalues; jj++) {
      if (byrow) {
        ans[ii + I[jj]*nrows] = ANS_NA;
      } else{
        ans[I[jj] + ii*nrows] = ANS_NA;
      }
    }

    /*
    Rprintf("\n");
     */
  }
}



/***************************************************************************
 HISTORY:
 2019-4-23 [BKM]
  o Added new tiesMethods: first, last, random, and dense.
 2015-06-12 [DJ]
  o Supported subsetted computation.
 2014-11-06 [HB]
 o CLEANUP: Moving away from R data types in low-level C functions.
 2013-04-23 [HB]
 o BUG FIX: Ranks did not work for integers with NAs; now using X_ISNAN().
 2013-01-13 [HB]
 o Template cleanup.  Extended template to integer matrices.
 o Added argument 'tiesMethod' to rowRanks().
 2012-12-14 [PL]
 o Added internal support for "min", "max" and "average" ties.  Using
   template to generate the various versions of the functions.
 2013-01-13 [HCB]
 o Created.  Using "max" ties.
 **************************************************************************/
