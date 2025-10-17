/*******************************************************************
   Copyright (C) 2001-7 Leo Breiman, Adele Cutler and Merck & Co., Inc.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the R_Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
*******************************************************************/

/******************************************************************
 * buildtree and findbestsplit routines translated from Leo's
 * original Fortran code.
 *
 *      copyright 1999 by leo Breiman
 *      this is free software and can be used for any purpose.
 *      It comes with no guarantee.
 *
 ******************************************************************/
#include <Rmath.h>
#include <R.h>
#include "rf.h"

/*build one tree, call this multiple times to create a forest in regrf*/
void regTree(double *x, double *y, int mdim, int nsample, int *lDaughter,
             int *rDaughter,
             double *splitCutoff, double *nodeMeans, int *nodestatus, int numNodes,
             int *treeSize, int nthsize, int mtry, int *bestSplitVarInNode, int *cat,
             double *nodeImpurityDecrease, int *varUsed) {
  int i, j, k, m, currentNode, *rowIndices, *nodestart, *numPointsInEachNode;
  int ndstart, ndend, leftChildLastDataPoint, currentNodeCount, splitResult, msplit;
  double label, sumOfSquares, nodeMean, decsplit, ubest, sumOfNodeData;

  nodestart = (int *) R_Calloc(numNodes, int);
  numPointsInEachNode   = (int *) R_Calloc(numNodes, int);

  /* initialize some arrays for the tree */
  zeroInt(nodestatus, numNodes);
  zeroInt(nodestart, numNodes);
  zeroInt(numPointsInEachNode, numNodes);
  zeroDouble(nodeMeans, numNodes);

  rowIndices = (int *) R_Calloc(nsample, int);
  for (i = 1; i <= nsample; ++i){
    rowIndices[i-1] = i;
  }
  currentNode = 0;
  nodestart[0] = 0;
  numPointsInEachNode[0] = nsample;
  nodestatus[0] = NODE_TOSPLIT;

  /* compute mean and sum of squares for Y */
  nodeMean = 0.0;
  sumOfSquares = 0.0;
  for (i = 0; i < nsample; ++i) {
    label = y[rowIndices[i] - 1];
    sumOfSquares += i * (nodeMean - label) * (nodeMean - label) / (i + 1);
    nodeMean = (i * nodeMean + label) / (i + 1);
  }
  nodeMeans[0] = nodeMean;

  /* start main loop */
  for (k = 0; k < numNodes - 2; ++k) {
    if (k > currentNode || currentNode >= numNodes - 2) {
      break;
    }
    /* skip if the node is not to be split */
    if (nodestatus[k] != NODE_TOSPLIT) {
      continue;
    }

#ifdef RF_DEBUG
    Rprintf("regTree: k=%d, av=%f, ss=%f\n", k, av, ss);
#endif

    /* initialize for next call to findbestsplit */
    ndstart = nodestart[k];
    ndend = ndstart + numPointsInEachNode[k] - 1;
    currentNodeCount = numPointsInEachNode[k];
    sumOfNodeData = currentNodeCount * nodeMeans[k];
    splitResult = 0;
    decsplit = 0.0;

#ifdef RF_DEBUG
    Rprintf("before findBestSplit: ndstart=%d, ndend=%d, jstat=%d, decsplit=%f\n",
            ndstart, ndend, jstat, decsplit);
#endif

    findBestSplit(x, rowIndices, y, mdim, nsample, ndstart, ndend, &msplit,
                  &decsplit, &ubest, &leftChildLastDataPoint, &splitResult, mtry, sumOfNodeData,
                  currentNodeCount, cat);
#ifdef RF_DEBUG
    Rprintf(" after findBestSplit: ndstart=%d, ndend=%d, jstat=%d, decsplit=%f, msplit=%d\n",
            ndstart, ndend, jstat, decsplit, msplit);

#endif
    if (splitResult == NODE_TERMINAL) {
      /* Node is terminal: Mark it as such and move on to the next. */
      nodestatus[k] = splitResult;
      continue;
    }
    /* Found the best split. */
    bestSplitVarInNode[k] = msplit;
    varUsed[msplit - 1] = 1;
    splitCutoff[k] = ubest;
    nodeImpurityDecrease[msplit - 1] += decsplit;
    nodestatus[k] = NODE_INTERIOR;

    /* left node is the node after current, right node is 2 after */
    numPointsInEachNode[currentNode + 1] = leftChildLastDataPoint - ndstart + 1;
    numPointsInEachNode[currentNode + 2] = ndend - leftChildLastDataPoint;
    nodestart[currentNode + 1] = ndstart;
    nodestart[currentNode + 2] = leftChildLastDataPoint + 1;

    /* compute mean and sum of squares for the left daughter node */
    nodeMean = 0.0;
    sumOfSquares = 0.0;
    for (j = ndstart; j <= leftChildLastDataPoint; ++j) {
      label = y[rowIndices[j]-1];
      m = j - ndstart;
      sumOfSquares += m * (nodeMean - label) * (nodeMean - label) / (m + 1);
      nodeMean = (m * nodeMean + label) / (m+1);
    }
    nodeMeans[currentNode+1] = nodeMean;
    nodestatus[currentNode+1] = NODE_TOSPLIT;
    if (numPointsInEachNode[currentNode + 1] <= nthsize) {
      nodestatus[currentNode + 1] = NODE_TERMINAL;
    }

    /* compute mean and sum of squares for the right daughter node */
    nodeMean = 0.0;
    sumOfSquares = 0.0;
    for (j = leftChildLastDataPoint + 1; j <= ndend; ++j) {
      label = y[rowIndices[j]-1];
      m = j - (leftChildLastDataPoint + 1);
      sumOfSquares += m * (nodeMean - label) * (nodeMean - label) / (m + 1);
      nodeMean = (m * nodeMean + label) / (m + 1);
    }
    nodeMeans[currentNode + 2] = nodeMean;
    nodestatus[currentNode + 2] = NODE_TOSPLIT;
    if (numPointsInEachNode[currentNode + 2] <= nthsize) {
      nodestatus[currentNode + 2] = NODE_TERMINAL;
    }

    /* map the daughter nodes */
    lDaughter[k] = currentNode + 1 + 1;
    rDaughter[k] = currentNode + 2 + 1;
    /* Augment the tree by two nodes. */
    currentNode += 2;
#ifdef RF_DEBUG
    Rprintf(" after split: ldaughter=%d, rdaughter=%d, ncur=%d\n",
            lDaughter[k], rDaughter[k], ncur);
#endif

  }
  *treeSize = numNodes;
  for (k = numNodes - 1; k >= 0; --k) {
    if (nodestatus[k] == 0) (*treeSize)--;
    if (nodestatus[k] == NODE_TOSPLIT) {
      nodestatus[k] = NODE_TERMINAL;
    }
  }
  R_Free(nodestart);
  R_Free(rowIndices);
  R_Free(numPointsInEachNode);
}

/*
call on every node in a tree. randomly choose some features and 
find which one has best split
*/
void findBestSplit(double *x, int *jdex, double *y, int mdim, int nsample,
                   int ndstart, int ndend, int *bestVarToReturn, double *decsplit,
                   double *bestSplitToReturn, int *ndendl, int *jstat, int mtry,
                   double sumnode, int nodeCount, int *cat) {
  int last, numCategoriesAllVars[MAX_CAT], icat[MAX_CAT], numCategoriesForVar, nl, nr, npopl, npopr, tieVar, tieVal;
  int i, j, kv, l, *varIndices, *ncase;
  double *xt, *ut, *v, *yl, sumcat[MAX_CAT], avcat[MAX_CAT], tavcat[MAX_CAT], valueAtBestSplit;
  double crit, bestSplitForAllVariables, bestSplitWithinVariable, suml, sumr, d, critParent;
  /*R_Calloc is different from calloc, r handles memory allocation instead of os*/
  ut = (double *) R_Calloc(nsample, double);
  xt = (double *) R_Calloc(nsample, double);
  v  = (double *) R_Calloc(nsample, double);
  yl = (double *) R_Calloc(nsample, double);
  varIndices  = (int *) R_Calloc(mdim, int);
  ncase = (int *) R_Calloc(nsample, int);
  zeroDouble(avcat, MAX_CAT);
  zeroDouble(tavcat, MAX_CAT);

  /* START BIG LOOP */
  *bestVarToReturn = -1;
  *decsplit = 0.0;
  bestSplitForAllVariables = 0.0;
  valueAtBestSplit = 0.0;
  for (i=0; i < mdim; ++i)
  {
    varIndices[i] = i;
  }

  last = mdim - 1;
  tieVar = 1;
  /*choose mtry number of variables, choose the one with best split*/
  for (i = 0; i < mtry; ++i) 
  {
    /*sample without replacement: choose random, move to end, 
    do next random choice in range 1 to len - n*/
    bestSplitWithinVariable = 0.0;
    j = (int) (unif_rand() * (last+1));
    kv = varIndices[j];
    swapInt(varIndices[j], varIndices[last]);
    last--;

    numCategoriesForVar = cat[kv];
    if (numCategoriesForVar == 1) {
      /* numeric variable */
      for (j = ndstart; j <= ndend; ++j) {
        xt[j] = x[kv + (jdex[j] - 1) * mdim]; /*indexing to represent 2d in a 1d vector */
        yl[j] = y[jdex[j] - 1];
      }
    } else {
      /* categorical variable */
      zeroInt(numCategoriesAllVars, MAX_CAT);
      zeroDouble(sumcat, MAX_CAT);
      for (j = ndstart; j <= ndend; ++j) {
        l = (int) x[kv + (jdex[j] - 1) * mdim];
        sumcat[l - 1] += y[jdex[j] - 1];
        numCategoriesAllVars[l - 1] ++;
      }
      /* Compute means of Y by category. */
      for (j = 0; j < numCategoriesForVar; ++j) {
        avcat[j] = numCategoriesAllVars[j] ? sumcat[j] / numCategoriesAllVars[j] : 0.0;
      }
      /* Make the category mean the `pseudo' X data. */
      for (j = 0; j < nsample; ++j) {
        xt[j] = avcat[(int) x[kv + (jdex[j] - 1) * mdim] - 1];
        yl[j] = y[jdex[j] - 1];
      }
    }
    /* copy the x data in this node. */
    for (j = ndstart; j <= ndend; ++j) {
      v[j] = xt[j];
    }

    for (j = 1; j <= nsample; ++j) {
      ncase[j - 1] = j;
    }

    R_qsort_I(v, ncase, ndstart + 1, ndend + 1);
    
    if (v[ndstart] >= v[ndend]) {
      continue;
    }
    /* ncase(n)=case number of v nth from bottom */
    /* Start from the right and search to the left. */
    critParent = sumnode * sumnode / nodeCount;
    suml = 0.0;
    sumr = sumnode;
    npopl = 0;
    npopr = nodeCount;
    crit = 0.0;
    tieVal = 1;
    /* Search through the "gaps" in the x-variable. */
    for (j = ndstart; j <= ndend - 1; ++j) {
      d = yl[ncase[j] - 1];
      suml += d;
      sumr -= d;
      npopl++;
      npopr--;
      if (v[j] < v[j+1]) {
        crit = (suml * suml / npopl) + (sumr * sumr / npopr) - critParent;
        if (crit > bestSplitWithinVariable) {
          valueAtBestSplit = (v[j] + v[j+1]) / 2.0;
          bestSplitWithinVariable = crit;
          tieVal = 1;
        }
        if (crit == bestSplitWithinVariable) {
          tieVal++;
          if (unif_rand() < 1.0 / tieVal) {
            valueAtBestSplit = (v[j] + v[j+1]) / 2.0;
            bestSplitWithinVariable = crit;
          }
        }
      }
    }
    if (bestSplitWithinVariable > bestSplitForAllVariables) {
      *bestSplitToReturn = valueAtBestSplit;
      *bestVarToReturn = kv + 1;
      bestSplitForAllVariables = bestSplitWithinVariable;
      for (j = ndstart; j <= ndend; ++j) {
        ut[j] = xt[j];
      }
      if (cat[kv] > 1) {
        for (j = 0; j < cat[kv]; ++j) tavcat[j] = avcat[j];
      }
      tieVar = 1;
    }
    if (bestSplitWithinVariable == bestSplitForAllVariables) {
      tieVar++;
      if (unif_rand() < 1.0 / tieVar) {
        *bestSplitToReturn = valueAtBestSplit;
        *bestVarToReturn = kv + 1;
        bestSplitForAllVariables = bestSplitWithinVariable;
        for (j = ndstart; j <= ndend; ++j) {
          ut[j] = xt[j];
        }
        if (cat[kv] > 1) {
          for (j = 0; j < cat[kv]; ++j) tavcat[j] = avcat[j];
        }
      }
    }

  }
  *decsplit = bestSplitForAllVariables;

  /* If best split can not be found, set to terminal node and return. */
  if (*bestVarToReturn != -1) {
    nl = ndstart;
    for (j = ndstart; j <= ndend; ++j) {
      if (ut[j] <= *bestSplitToReturn) {
        nl++;
        ncase[nl-1] = jdex[j];
      }
    }
    *ndendl = imax2(nl - 1, ndstart);
    nr = *ndendl + 1;
    for (j = ndstart; j <= ndend; ++j) {
      if (ut[j] > *bestSplitToReturn) {
        if (nr >= nsample) break;
        nr++;
        ncase[nr - 1] = jdex[j];
      }
    }
    if (*ndendl >= ndend) *ndendl = ndend - 1;
    for (j = ndstart; j <= ndend; ++j) jdex[j] = ncase[j];

    numCategoriesForVar = cat[*bestVarToReturn - 1];
    if (numCategoriesForVar > 1) {
      for (j = 0; j < numCategoriesForVar; ++j) {
        icat[j] = (tavcat[j] < *bestSplitToReturn) ? 1 : 0;
      }
      *bestSplitToReturn = pack(numCategoriesForVar, icat);
    }
  } else *jstat = NODE_TERMINAL;

  R_Free(ncase);
  R_Free(varIndices);
  R_Free(v);
  R_Free(yl);
  R_Free(xt);
  R_Free(ut);
}

/*====================================================================*/
void predictRegTree(double *x, int nsample, int mdim,
                    int *lDaughter, int *rDaughter, int *nodestatus,
                    double *ypred, double *split, double *nodepred,
                    int *splitVar, int treeSize, int *cat, int maxcat,
                    int *nodex) {
  int i, j, k, m, *cbestsplit;
  double dpack;

  /* decode the categorical splits */
  if (maxcat > 1) {
    cbestsplit = (int *) R_Calloc(maxcat * treeSize, int);
    zeroInt(cbestsplit, maxcat * treeSize);
    for (i = 0; i < treeSize; ++i) {
      if (nodestatus[i] != NODE_TERMINAL && cat[splitVar[i] - 1] > 1) {
        dpack = split[i];
        /* unpack `npack' into bits */
        /* unpack(dpack, maxcat, cbestsplit + i * maxcat); */
        for (j = 0; j < cat[splitVar[i] - 1]; ++j) {
          cbestsplit[j + i*maxcat] = ((unsigned long) dpack & 1) ? 1 : 0;
          dpack = dpack / 2.0 ;
          /* cbestsplit[j + i*maxcat] = npack & 1; */
        }
      }
    }
  }

  for (i = 0; i < nsample; ++i) {
    k = 0;
    while (nodestatus[k] != NODE_TERMINAL) { /* go down the tree */
          m = splitVar[k] - 1;
      if (cat[m] == 1) {
        k = (x[m + i*mdim] <= split[k]) ?
        lDaughter[k] - 1 : rDaughter[k] - 1;
      } else {
        /* Split by a categorical predictor */
        k = cbestsplit[(int) x[m + i * mdim] - 1 + k * maxcat] ?
        lDaughter[k] - 1 : rDaughter[k] - 1;
      }
    }
    /* terminal node: assign prediction and move on to next */
    ypred[i] = nodepred[k];
    nodex[i] = k + 1;
  }
  if (maxcat > 1) R_Free(cbestsplit);
}
