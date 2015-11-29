library("microbenchmark")
library("Rcpp")

source("zahler_1.R")
sourceCpp("zahler_1cpp.cpp")

formcalc_raw <- read.table("C_634_FORMULAE.dat",skip=15,header=T)
mol_type <- classify_mol(formcalc_raw)

zahlerR <- zahler_1(formcalc_raw,
                    Intensity_Max = 235000000,
                    Intensity = formcalc_raw$Intensity, 
                    OC = formcalc_raw$O/formcalc_raw$C, 
                    ExpMass = formcalc_raw$ExpMass) 
zahlerCpp <- zahler1_cpp(intensity = formcalc_raw$Intensity, 
                         OC = formcalc_raw$O/formcalc_raw$C,
                         expmass = formcalc_raw$ExpMass, 
                         intensity_max = 235000000, 
                         intensity_min = 0, 
                         expmass_min = 149,
                         expmass_max = 800, 
                         OtoC_min = 0, 
                         OtoC_max = 1)

identical(zahlerR, zahlerCpp)
#TRUE

microbenchmark(
    "zahler_1" = zahler_1(formcalc_raw,
                          Intensity_Max = 235000000,
                          Intensity = formcalc_raw$Intensity, 
                          OC = formcalc_raw$O/formcalc_raw$C, 
                          ExpMass = formcalc_raw$ExpMass) , 
    "zahler1_cpp" = zahler1_cpp(intensity = formcalc_raw$Intensity, 
                                OC = formcalc_raw$O/formcalc_raw$C,
                                expmass = formcalc_raw$ExpMass, 
                                intensity_max = 235000000, 
                                intensity_min = 0, 
                                expmass_min = 149,
                                expmass_max = 800, 
                                OtoC_min = 0, 
                                OtoC_max = 1)
)

#44 times faster, but R code was not optimised before and for the instant the 
#Cpp function does not include defensif statements


