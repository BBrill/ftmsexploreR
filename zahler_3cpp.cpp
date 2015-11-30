#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericVector zahler_3cpp(NumericVector HIon, NumericVector C, NumericVector O, NumericVector N,
double AI_max, double AI_min, double NtoC_max) {
    
    int n = C.length();
    NumericVector out(n);
    
    for (int i = 0; i < n; ++i) {
        double HIoni = HIon[i];
        double Ci = C[i];
        double Ni = N[i];
        double Oi = O[i];
        
        if ((Ci - Oi - Ni)==0) {
            out[i] = 1;
        } else if ((1 + Ci - Oi - 0.5 * HIoni)/(Ci - Oi - Ni) > AI_min){
            if((1 + Ci - Oi - 0.5 * HIoni)/(Ci - Oi - Ni) < AI_max) {
                if (Ni/Ci <= NtoC_max) {
                    out[i] = 1; 
                }
            }
        } else {
            out[i] = 0; 
        }
    }   
    return out;
}
    