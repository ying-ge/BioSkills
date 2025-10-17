#include <iostream>
#include <algorithm>
#include <Rcpp.h>

typedef struct PeakSet {
  Rcpp::NumericVector chr;
  Rcpp::NumericVector left;
  Rcpp::NumericVector right;
  Rcpp::NumericVector score;
} PeakSet;

int mergeSet(PeakSet dest,PeakSet src,int maxGap) {
  int vLen = src.chr.size();
  int di = 0;

  dest.chr[0] = src.chr[0];
  dest.left[0] = src.left[0];
  dest.right[0] = src.right[0];

  for (int si=1;si<vLen;si++) {
    if (dest.chr[di] == src.chr[si] && src.left[si] - dest.right[di] <= maxGap) {
      // want to merge
      dest.right[di] = std::max(dest.right[di],src.right[si]);
    } else {
      di++;
      dest.chr[di] = src.chr[si];
      dest.left[di] = src.left[si];
      dest.right[di] = src.right[si];
    }
  }
  return di+1;
}

bool validData(Rcpp::DataFrame peaks,bool score=false) {
  return false;
}

// [[Rcpp::export]]
Rcpp::DataFrame mergePeaks(Rcpp::DataFrame data,int maxGap) {
  PeakSet src;
  src.chr = data[0];
  src.left = data[1];
  src.right = data[2];
  
  int sourceLen = src.chr.size();
  int total;
  PeakSet dest;
  dest.chr = Rcpp::NumericVector(sourceLen);
  dest.left = Rcpp::NumericVector(sourceLen);
  dest.right = Rcpp::NumericVector(sourceLen);

  total = mergeSet(dest,src,maxGap);

  Rcpp::NumericVector fChr(total);
  Rcpp::NumericVector fLeft(total);
  Rcpp::NumericVector fRight(total);
  for (int i=0;i<total;i++) {
    fChr[i] = dest.chr[i];
    fLeft[i] = dest.left[i];
    fRight[i] = dest.right[i];
  }
  Rcpp::DataFrame x = Rcpp::DataFrame::create(Rcpp::Named("chr")=fChr,
                                              Rcpp::Named("left")=fLeft,
                                              Rcpp::Named("right")=fRight);
  return x;
}

// [[Rcpp::export]]
Rcpp::List mergeScores(Rcpp::DataFrame sMerged,Rcpp::NumericVector sScore, Rcpp::DataFrame sPeaks, Rcpp::Nullable<Rcpp::LogicalVector> abs = R_NilValue) {
  PeakSet merged;
  PeakSet peaks;
  merged.chr = sMerged[0];
  merged.left = sMerged[1];
  merged.right = sMerged[2];
  peaks.chr = sPeaks[0];
  peaks.left = sPeaks[1];
  peaks.right = sPeaks[2];
  peaks.score = sPeaks[3];
  int vLen = merged.chr.size();
  int peakLen = peaks.chr.size();
  Rcpp::NumericVector score(vLen);
  Rcpp::NumericVector included(vLen);
  Rcpp::NumericVector absScore(vLen);
  int pi = 0;
  bool abso = false;
  if (abs.isNotNull()) {
    Rcpp::LogicalVector absolute(abs);
    abso = absolute[0];
  }

  for (int mi=0;mi<vLen;mi++) {
    while (pi < peakLen &&
           merged.chr[mi] == peaks.chr[pi] &&
           merged.left[mi] <= peaks.left[pi] &&
           merged.right[mi] >= peaks.right[pi]) {
      // peak is contained in region
      if (abso) {
        // have to retain which of the 3 scores is max(abs(...)), but also the original (signed) score.
        double intermediate;
        if (std::abs(score[mi]) > std::abs(peaks.score[pi])) {
          intermediate = score[mi];
        } else if (std::abs(score[mi]) < std::abs(peaks.score[pi])) {
          intermediate = peaks.score[pi];
        } else if (score[mi] > 0 || peaks.score[pi] > 0) {
          intermediate = std::abs(score[mi]); // they're equal in magnitude, and at least one is positive
        } else {
          intermediate = score[mi]; // they're equal in magnitude, and both are negative
        }
        if (std::abs(intermediate) > std::abs(sScore[mi])) {
          score[mi] = intermediate;
        } else if (std::abs(intermediate) < std::abs(sScore[mi])) {
          score[mi] = sScore[mi];
        } else if (intermediate > 0 || sScore[mi] > 0) {
          score[mi] = std::abs(intermediate);
        } else {
          score[mi] = intermediate;
        }
        //score[mi] = std::max(std::abs(sScore[mi]),std::max(std::abs(score[mi]),std::abs(peaks.score[pi])));
      } else {
        score[mi] = std::max(sScore[mi],std::max(score[mi],peaks.score[pi]));
      }
      included[mi] = 1;
      pi++;
    }
    if (!included[mi]) {
      score[mi] = sScore[mi];
    }
  }
  return Rcpp::List::create(Rcpp::Named("score")=score,
                            Rcpp::Named("included")=included);
}
