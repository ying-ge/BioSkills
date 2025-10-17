#include "utils.h" 

#include <vector>
#include <memory>

//' @importFrom Rcpp sourceCpp
//' @useDynLib SingleR
//[[Rcpp::export(rng=false)]]
SEXP prebuild(Rcpp::RObject ref, Rcpp::IntegerVector labels, Rcpp::List markers, bool approximate, int nthreads) {
    singlepp::BasicBuilder builder;
    builder.set_num_threads(nthreads);

    // Use all available markers; assume subsetting was applied on the R side.
    builder.set_top(-1).set_approximate(approximate);
    
    // Setting up the markers.
    singlepp::Markers markers2(markers.size());
    for (size_t m = 0; m < markers.size(); ++m) {
        Rcpp::List curmarkers(markers[m]);
        auto& curmarkers2 = markers2[m];
        curmarkers2.resize(curmarkers.size());

        for (size_t n = 0; n < curmarkers.size(); ++n) {
            Rcpp::IntegerVector seq(curmarkers[n]);
            auto& seq2 =  curmarkers2[n];
            seq2.insert(seq2.end(), seq.begin(), seq.end());
        }
    }

    // Building the indices.
    auto parsed = Rtatami::BoundNumericPointer(ref);
    auto built = builder.run(parsed->ptr.get(), static_cast<const int*>(labels.begin()), std::move(markers2));

    // Moving it into the external pointer.
    return PrebuiltXPtr(new singlepp::BasicBuilder::Prebuilt(std::move(built)), true);
}

//[[Rcpp::export(rng=false)]]
Rcpp::IntegerVector get_subset(SEXP built) {
    PrebuiltXPtr ptr(built);
    return Rcpp::IntegerVector(ptr->subset.begin(), ptr->subset.end());
}

//[[Rcpp::export(rng=false)]]
Rcpp::LogicalVector is_valid_built(SEXP built) {
    return Rf_ScalarLogical(!!R_ExternalPtrAddr(built));
}
