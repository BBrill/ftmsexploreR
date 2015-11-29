#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector zahler2_cpp(NumericVector HIon, NumericVector C, NumericVector N,
NumericVector O, double DBEtoC_min, double DBEtoC_max, double OplusNtoC_max) {
    
    int n = HIon.length();
    NumericVector out(n);
    
    for (int i = 0; i < n; ++i) {
        double HIoni = HIon[i];
        double Ci = C[i];
        double Ni = N[i];
        double Oi = O[i];
        
        if ((1 + 0.5 * (2 * Ci - HIoni + Ni))/Ci > DBEtoC_min) {
            if ((1 + 0.5 * (2 * Ci - HIoni + Ni))/Ci < DBEtoC_max){
                if((Oi + Ni)/Ci < OplusNtoC_max) {
                    if(Ci > 5) {
                        out[i] = 1;
                    }
                }
            }
        } else {
            out[i] = 0;
        }
    }
    return out;
}
