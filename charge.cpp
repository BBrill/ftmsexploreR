#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericVector charge(NumericVector ExpMass) {
    
    int  n = ExpMass.length();
    NumericVector out(n);
    
    for (int i = 0; i < n; ++i) {
        
        double ExpMassi = ExpMass[i];
        
        if(ExpMassi != 0) {
            out[i] = -1;
        } else {
            out[i] = 0;
        }
    }
    return(out);
}

    
/*
Original R code
chargeR <- ifelse(formcalc_raw$ExpMass!=0,-1,0)       
*/