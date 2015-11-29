#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector zahler1_cpp(NumericVector, x, double a, double b) {
    int n = x.length();
    NumericVector out(n);
    for (int i = 0; i < n; ++i) {
        double xi = x[i];
        if (OC > OtoC_min) {
        out[i] = 1;
        } else if (intensity < intensity_max) {
            out[i] = 1;
        } else if (intensity > intensity_min) {
            out[i] = 1;
        } else if (expmass < expmass_max) {
            out[i] = 1;
        } else if (expmass > expmass_min) {
            out[i] = 1;
        } else {
            out[i] = 0;
        }
    }
    return out;
}
