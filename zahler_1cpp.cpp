#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector zahler1_cpp(NumericVector intensity, NumericVector OC, NumericVector expmass, 
double intensity_min, double intensity_max, double expmass_min, double expmass_max,
double OtoC_min, double OtoC_max) {
    
    int n = intensity.length();
    NumericVector out(n);
    
    for (int i = 0; i < n; ++i) {
        double OCi = OC[i];
        double intensityi = intensity[i];
        double expmassi = expmass[i];
        
        if (OCi > OtoC_min) {
            if (OCi < OtoC_max){
                if(intensityi > intensity_min) {
                    if(intensityi < intensity_max) {
                        if (expmassi > expmass_min) { 
                            if (expmassi < expmass_max) {
                                out[i] = 1;
                            }
                        }
                    }
                }
            }
        } else {
            out[i] = 0;
        }
    }
    return out;
}
