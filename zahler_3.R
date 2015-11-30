#' Zahler 3
zahler_3 <- function(data, AI_max = NULL, AI_min = NULL, NtoC_max = 1) {
  
  if(is.null(AI_max)){
    stop("You should set AI_max")
  }
  
  if(is.null(AI_min)){
    stop("You should set AI_min")
  }
  
  attach(data)
  zahler3 <- ifelse((C-O-N)==0, 1,
                  ifelse(((1+C-O-0.5*HIon)/(C-O-N))<AI_max,
                         ifelse(((1+C-O-0.5*HIon)/(C-O-N))>AI_min,
                                ifelse(N/C<=NtoC_max,1
                                       ,0),0),0))
  detach(data)
  
  return(zahler3)
}
