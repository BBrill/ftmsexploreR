#' Compute Zahler 1 values

zahler_1 <- function(data, OtoC_min = 0, OtoC_max = 1, Intensity_Min=0, 
                     Intensity_Max = 235000000, ExpMass_Min = 149, 
                     ExpMass_Max = 800, OC = NULL, 
                     Intensity = NULL, ExpMass = NULL) {
 
#you should oblige OC be the same length of Intensity and Expmass
#otherwise you should oblige data to be of a fixed format with cols of 
    #fixed names named cols
    
  if(is.null(Intensity)){
    stop("You should insert Intensity data")
  }
  
  if(is.null(OC)){
    stop("You should insert Oxygen/Carbon data")
  }
  
  if(is.null(ExpMass)){
      stop("You should insert ExpMass data")
  }
  
  zahler1 <- ifelse(OC<OtoC_max,
                    ifelse(OC>OtoC_min,
                           ifelse(Intensity<Intensity_Max,
                                  ifelse(Intensity>Intensity_Min,
                                         ifelse(ExpMass<ExpMass_Max,
                                                ifelse(ExpMass>ExpMass_Min,1
                                                       ,0),0),0),0),0),0)
  
  return (zahler1)
  
}
