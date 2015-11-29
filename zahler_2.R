#' zahler_2 function to compute
zahler_2 <- function(data, DBEtoC_min = NULL, DBEtoC_max = NULL, OplusNtoC_max = NULL) {
  if(is.null(DBEtoC_min)){
    stop("You should insert DBEtoC_min")
  }
  
  if(is.null(DBEtoC_max)){
    stop("You should insert DBEtoC_max")
  }
  
  if(is.null(OplusNtoC_max)){
    stop("You should insert OplusNtoC_max")
  }
  
  attach(data)
  zahler2 <- ifelse((1+0.5*(2*C-HIon+N))/C<DBEtoC_max,
                    ifelse((1+0.5*(2*C-HIon+N))/C>DBEtoC_min,
                           ifelse((O+N)/C<OplusNtoC_max,
                                  ifelse(C>5,1
                                         ,0),0),0),0)
  detach(data)
  
  return(zahler2)
}