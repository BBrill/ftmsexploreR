#setwd("D:/Documents/R-FTMS-proj")
formcalc_raw <- read.table("C_634_FORMULAE.dat",skip=15,header=T)
library(extremevalues)
library(plyr)
library(ggplot2)
source("classify_mol.R")
source("zahler_1.R")
source("zahler_2.R")
source("zahler_3.R")

#here adapted to the molecular constraints: Hx-C100-O80-N2-S1_________####
mol_type <- classify_mol(formcalc_raw)

attach(formcalc_raw)
HC <- HIon/C
OC <- O/C
NC <- N/C
detach(formcalc_raw)

#________________________Check high Intensity_____________####
L <- getOutliers(formcalc_raw$Intensity, method="I", distribution="lognormal")
# outlierPlot(formcalc_raw$Intensity, L, mode="qq")
#L$iRight
#L$limit['Right']

#________________________N=0, mass even: N0_mass_even_____####
#set NtoC_max
NtoC_max <- 1
N0_mass_even <- ifelse(NtoC_max==0 , 
                      ifelse(mod(round(ExpMass+1, digits=0), 2)==0,
                             ifelse(N==0, 1, 0),
                             0),
                      1)
#________________________Nregel___________________________####
Nregel <- 1
#________________________Zahler 1_________________________####
zahler1 <- zahler_1(formcalc_raw, Intensity_Max = L$limit['Right'],
                   Intensity = formcalc_raw$Intensity, 
                   OC = formcalc_raw$O/formcalc_raw$C, 
                   ExpMass = formcalc_raw$ExpMass)

#________________________Zahler 2_________________________####
#set DBEtoC_min & DBEtoC_max & OplusNtoC_max
zahler2 <- zahler_2(formcalc_raw, DBEtoC_min = 0, DBEtoC_max = 5, OplusNtoC_max = 3)

#________________________Zahler 3_________________________####
#set AI min & max
zahler3 <- zahler_3(formcalc_raw, AI_max = 1, AI_min = -20)

attach(formcalc_raw)
#________________________H/C charge condition:HCcc________####
#set charge... I write this way to obtain a vector with the right length.
charge <- ifelse(formcalc_raw$ExpMass!=0,-1,0)
HCcc <- ifelse(Nregel*zahler1*zahler2*zahler3==1,
    ifelse(charge==(-1),((HIon+1)/C),((HIon-1)/C))
       ,0)
#________________________H/C filtered: HCf________________####
#set HC_min & HC_max
HC_min <- 0
HC_max <- 2.5  
HCf <- ifelse(HCcc<HC_max,
    ifelse(HCcc>HC_min,HCcc,0)
    ,0)
#________________________O/C filtered: OCf________________####
OCf <- ifelse(Nregel*zahler1*zahler2*zahler3==1,OC,0)
#________________________N/C filtered: NCf________________####
NCf <- ifelse(Nregel*zahler1*zahler2*zahler3==1,NC,0)  
#________________________DBE______________________________####
dbe <- ifelse(HCf>0,
      ifelse(Nregel*zahler1*zahler2*zahler3==1,(1+0.5*(2*C-HIon+N-1)),0)
         ,0)
#________________________AI_______________________________####
AI <- ifelse((C-O-N-S)>0,
     ifelse((1+C-O-0.5*(HIon+1)-S)>0,
    ifelse(HCf>0,
   ifelse(Nregel*zahler1*zahler2*zahler3==1,(1+C-O-0.5*(HIon+1)-S)/(C-O-N-S),0)
         ,0),0),0)
#________________________Xc_______________________________####
#set m & n values
m <- 1
n <- 1
Xc <- ifelse((3*(dbe-(m*O+n*S))-1)/(dbe-(m*O+n*S))<=0,0,
             ifelse(dbe<=(m*O+n*S),0,(3*(dbe-(m*O+n*S))-1)/(dbe-(m*O+n*S))))
#________________________KMD______________________________####
nommass <- ifelse(HCcc==0,0,ifelse(Nregel*zahler1*zahler2*zahler3==1,
                                   floor(formcalc_raw$ExpMass),0))
kmd <- ifelse(HCcc==0,0,ifelse(Nregel*zahler1*zahler2*zahler3==1,
                               formcalc_raw$ExpMass-nommass,0))
#________________________KMD CH2__________________________####
km_CH2 <- (formcalc_raw$ExpMass+1.007825-0.000549)*14/14.01565
nommass_CH2 <- ifelse(Nregel*zahler1*zahler2*zahler3==1,
                      ifelse(HCcc==0,0,ceiling(km_CH2)),0)
kmd_CH2 <- ifelse(Nregel*zahler1*zahler2*zahler3==1,
                  ifelse(HCcc==0,0,(nommass_CH2-km_CH2)),0)
#________________________KMD COO__________________________####
km_COO <- (formcalc_raw$ExpMass+1.007825-0.000549)*44/43.989829
nommass_COO <- ifelse(Nregel*zahler1*zahler2*zahler3==1,
                      ifelse(HCcc==0,0,ceiling(km_COO)),0)
kmd_COO <- ifelse(Nregel*zahler1*zahler2*zahler3==1,
                  ifelse(HCcc==0,0,(nommass_COO-km_COO)),0)
#________________________compile in dataframe_____________####
dfVK_raw <- data.frame(formcalc_raw, mol_type, 
                       HCf, OCf, NCf, 
                       dbe, AI, Xc,
                       nommass_CH2, kmd_CH2,
                       nommass_COO, kmd_COO)
dfVK_filter <- subset(
                subset(
                 subset(
                  subset(dfVK_raw,mol_type!="NA"),
                OCf!=0),
              HCf!=0),
              NCf<1)
detach(formcalc_raw)
#____________________________Graph_1:Mass spectra_____________________####
 
ggplot(formcalc_raw, aes(x=ExpMass, y=Intensity)) + 
  geom_segment(aes(xend=ExpMass), yend=0, colour="royalblue1") +
  scale_x_continuous(name="m/z", limits=c(149, 1000)) +
  scale_y_continuous(name="Intensity", limits=c(0, L$limit['Right'])) +
  geom_point(size=0) +
  theme(axis.text.x = element_text(size=16, colour = "black"), axis.text.y = element_text(size=12, colour="black"), axis.title.x = element_text(size=20), axis.title.y = element_text(size=20))
#____________________________Graph_2:counts of formula________________####
Hf <- dfVK_filter$HIon
Cf <- dfVK_filter$C
Of <- dfVK_filter$O
Nf <- dfVK_filter$N
Sf <- dfVK_filter$S

CHOf <- (Cf>=1 & Hf>=1 & Of>=1 & Nf==0 & Sf==0)
CHONf <- (Cf>=1 & Hf>=1 & Of>=1 & Nf==1 & Sf==0)
CHOSf <- (Cf>=1 & Hf>=1 & Of>=1 & Nf==0 & Sf==1)
CHONSf <- (Cf>=1 & Hf>=1 & Of>=1 & Nf==1 & Sf==1)
CHON2f <- (Cf>=1 & Hf>=1 & Of>=1 & Nf==2 & Sf==0)
CHONS2f <- (Cf>=1 & Hf>=1 & Of>=1 & Nf==2 & Sf==1)

count_mol_type <- ifelse(CHOf==T , "CHO", 
                ifelse(CHONf==T, "CHON",
              ifelse(CHOSf==T, "CHOS",
            ifelse(CHONSf==T, "CHONS",
          ifelse(CHON2f==T, NA,
         ifelse(CHONS2f==T, NA, NA)
        )))))

count_mol_type_N2 <- ifelse(CHOf==T , NA,
                    ifelse(CHONf==T, NA,
                   ifelse(CHOSf==T, NA,
                  ifelse(CHONSf==T, NA,
                 ifelse(CHON2f==T, "CHON",
                ifelse(CHONS2f==T, "CHONS", NA)
               )))))

d <- data.frame(dfVK_filter, count_mol_type, count_mol_type_N2)

summary(d$count_mol_type)->mol_count
summary(d$count_mol_type_N2)->mol_count_N2

nam_mol_count <- names(mol_count)
nam_mol_count_N2 <- names(mol_count_N2)

df_mol_count <- data.frame(nam_mol_count, mol_count)
df_mol_count_N2 <- data.frame(nam_mol_count_N2, mol_count_N2)

mol_merge <- merge(df_mol_count, df_mol_count_N2, 
                   by.x = "nam_mol_count", by.y = "nam_mol_count_N2",
                   all=T)
mol_merge[is.na(mol_merge)] <- 0

te <- c(mol_merge$mol_count, mol_merge$mol_count_N2)
tename <-c(nam_mol_count,nam_mol_count)
ten2 <- c(1,1,1,1,1,2,2,2,2,2)

tedf <- data.frame(tename,te,ten2)
tedf <- ddply(subset(tedf,te!=0), "tename", transform, label_y=cumsum(te))

tecol2 <- c("royalblue1", 
            "darkorange1", "darkorange3",
            "forestgreen",
            "red", "darkred")

ggplot(subset(tedf,tename!="NA's"), aes(x=tename, y=te, fill=factor(ten2))) + 
  geom_bar(stat="identity", fill=tecol2, colour="black", width=0.6) + 
  scale_x_discrete(limits=c("CHO","CHON","CHOS","CHONS"),name="")  +
  scale_y_continuous(name="counts", limits=c(0, 5000), breaks=c(0, 1000, 2000, 3000, 4000)) +
  geom_text(aes(y=label_y, label=te), vjust=1.5, colour="black") +
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_text(face="bold", size=16, colour = c("royalblue1", "darkorange1", "forestgreen", "red")),
        axis.title.y = element_text(size=20), axis.text.y = element_text(size=12, colour="black"))
#____________________________Graph_3:vanKrevelen diagram: H/C vs O/C__####
ggplot(dfVK_filter, aes(x=OCf, y=HCf, color=mol_type, size=Intensity)) +
  scale_size_area(max_size=6) +
  scale_x_continuous(name="O/C") +
  scale_y_continuous(name="H/C") +
  geom_point() + 
  scale_colour_manual(values=c("royalblue1", "darkorange1", "darkorange3",
                               "darkred", "red", "forestgreen"))
#____________________________Graph_4:vanKrevelen diagram: H/C vs m/z__####
ggplot(dfVK_filter, aes(x=ExpMass, y=HCf, color=mol_type, size=Intensity)) +
  scale_size_area(max_size=6) +
  scale_x_continuous(name="m/z") +
  scale_y_continuous(name="H/C") +
  geom_point() + 
  scale_colour_manual(values=c("royalblue1", "darkorange1", "darkorange3",
                               "darkred", "red", "forestgreen"))
#____________________________Graph_5:AI vs C__________________________####
ggplot(dfVK_filter, aes(x=C, y=AI, color=mol_type, size=Intensity)) +
  scale_size_area(max_size=6) +
  scale_x_continuous(name="C atoms") +
  scale_y_continuous(name="AI") +
  geom_point() + 
  scale_colour_manual(values=c("royalblue1", "darkorange1", "darkorange3",
                               "darkred", "red", "forestgreen"))
#____________________________Graph_6:Xc vs C__________________________####
ggplot(dfVK_filter, aes(x=C, y=Xc, color=mol_type, size=Intensity)) +
  scale_size_area(max_size=6) +
  scale_x_continuous(name="m/z") +
  scale_y_continuous(name="Xc", limits=c(2.4,3)) +
  geom_point() + 
  scale_colour_manual(values=c("royalblue1", "darkorange1", "darkorange3",
                               "darkred", "red", "forestgreen"))
#____________________________Graph_7:KMD CH2__________________________####
ggplot(dfVK_filter, aes(x=nommass_CH2, y=kmd_CH2, color=mol_type, size=Intensity)) +
  scale_size_area(max_size=6) +
  scale_x_continuous(name="m/z") +
  scale_y_continuous(name="kmd CH2") +
  geom_point() + 
  scale_colour_manual(values=c("royalblue1", "darkorange1", "darkorange3",
                               "darkred", "red", "forestgreen"))
#____________________________Graph_8:KMD COO__________________________####
ggplot(dfVK_filter, aes(x=nommass_COO, y=kmd_COO, color=mol_type, size=Intensity)) +
  scale_size_area(max_size=6) +
  scale_x_continuous(name="m/z") +
  scale_y_continuous(name="kmd COO") +
  geom_point() + 
  scale_colour_manual(values=c("royalblue1", "darkorange1", "darkorange3",
                               "darkred", "red", "forestgreen"))









