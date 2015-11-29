#' Function to classify molecules from raw FTMS data

classify_mol <- function(data) {
  attach(data)
  mol_type <- ifelse(C>=1 & HIon>=1 & O>=1 & N==0 & S==0, "CHO",
                     ifelse(C>=1 & HIon>=1 & O>=1 & N==1 & S==0, "CHON",
                            ifelse(C>=1 & HIon>=1 & O>=1 & N==0 & S==1, "CHOS",
                                   ifelse(C>=1 & HIon>=1 & O>=1 & N==1 & S==1, "CHONS",
                                          ifelse(C>=1 & HIon>=1 & O>=1 & N==2 & S==0, "CHON2",
                                                 ifelse(C>=1 & HIon>=1 & O>=1 & N==2 & S==1, "CHON2S", NA)
                                          )
                                   )
                            )
                     )
  )
  detach(data)
  return(mol_type)
}