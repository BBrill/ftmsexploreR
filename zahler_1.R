#' Compute Zahler 1 values

zahler_1 <- function(data,OtoC_min = 0, OtoC_max = 1, Intensity_Min=0, 
                     Intensity_Max = 235000000, ExpMass_Min = 149, 
                     ExpMass_Max = 800, OC = NULL, 
                     Intensity = NULL) {
  
  if(is.null(Intensity)){
    stop("You should insert Intensity data")
  }
  
  if(is.null(OC)){
    stop("You should insert Oxygen/Carbon data")
  }
  
  zahler1 <- ifelse(OC<OtoC_max,
                    ifelse(OC>OtoC_min,
                           ifelse(formcalc_raw$Intensity<Intensity_Max,
                                  ifelse(formcalc_raw$Intensity>Intensity_Min,
                                         ifelse(formcalc_raw$ExpMass<ExpMass_Max,
                                                ifelse(formcalc_raw$ExpMass>ExpMass_Min,1
                                                       ,0),0),0),0),0),0)
  
  return (zahler1)
  
}