# TB model in R

rm(list=ls())

# get the ODE solver
library(deSolve)

# make the function of the system
model = function (current_timepoint, state_values, parameters)
{
  # create state variables (local variables)
  Vaccine = state_values [1]       
  Makrofag = state_values [2]
  IL12 = state_values [3]
  Th0_IL12 = state_values [4]
  Th0_IL17 = state_values [5]
  Th1 = state_values [6]
  IL2 = state_values [7]
  IFNg = state_values [8]
  IFNgk = state_values [9]
  TNFa = state_values [10]
  Dendrit = state_values [11]
  IL23 = state_values [12]
  Th17 = state_values [13]
  IL17 = state_values [14]
  GranulotcytKnoglemarv = state_values[15]
  
  
  #      ifelse current_timepoint >= stop
  #      ifelse (dVaccine =  initialværdi / (- Kv_m -  Kv_d - Kv_Gk)
  # Evt. skal dødsraterne ganges på med mængden
  
 # stop = 
  
  with ( 
    as.list (parameters),     # variable names within parameters can be used 
    {
      Kv_m <- ifelse(Vaccine0 > 0.02, Kv_m.value, 0.001)
      Kv_d <- ifelse(Vaccine0 > 0.02, Kv_d.value, 0.2)
      #Kv_Gk <- ifelse(Vaccine > 0, Kv_Gk.value, 0)
      #Vaccine <- ifelse(Kv_d > 0, Vaccine,0)
      
      dVaccine = - (Kv_m) - (Kv_d) - Kv_Gk
      dMakrofag = (Kv_m) + TNFa * KTNFa_m + IFNgk * KIFNgk_m - Makrofag*sigma_m
      dIL12 = Makrofag * Km_IL12 - IL12 * KIL12_Th0 - IL12*sigma_IL12
      dTh0_IL12 = Makrofag * Km_Th0 + IL12 * KIL12_Th0 - Th0_IL12 * KTh0_Th1 + IL2 * KIL2_Th0  - Th0_IL12*sigma_Th0_IL12
      dTh0_IL17 = IL23 * KIL23_Th0 - Th0_IL17 * KTh0_Th17 - Th0_IL17*sigma_Th0_IL17
      dTh1 = Th0_IL12 * KTh0_Th1 - Th1*sigma_Th1
      dIL2 = Th1 * KTh1_IL2 - IL2*sigma_IL2
      dIFNg = Th1 * KTh1_IFNg + Th17 * KTh17_IFNg + Th0_IL12 * KTh0_IL12_IFNg + IFNgk * KIFNgk_IFNg - IFNg * KIFNg_IFNgk + Th0_IL17 * KTh0_IL17_IFNg - IFNg*sigma_IFNg
      dIFNgk = - IFNgk * KIFNgk_IFNg + IFNg * KIFNg_IFNgk - IFNgk*sigma_IFNgk
      dTNFa = Makrofag * Km_TNFa - TNFa*sigma_TNFa
      dDendrit = Kv_d - Dendrit*sigma_d
      dIL23 = Dendrit * Kd_IL23 - IL23*sigma_IL23
      dTh17 = Th0_IL17 * KTh0_Th17 - Th17*sigma_Th17
      dIL17 = Th17 * KTh17_IL17 - IL17*sigma_IL17
      dGranulocytKnoglemarv =  Kv_Gk + IL17 * KIL17_Gk - GranulotcytKnoglemarv*sigma_Gk 
      
      
      
      # combine results
      results = c(dVaccine, dMakrofag, dIL12, dTh0_IL12, dTh0_IL17, dTh1, dIL2, dIFNg, dIFNgk, dTNFa, dDendrit,
                  dIL23, dTh17, dIL17, dGranulocytKnoglemarv)
      list (results)
    }
  )
}

# parameters
KTh0_IL17_IFNg.value <- 0.00009
KTh0_IL12_IFNg.value <- 0.00009
Kv_m.value <- 0.0002
Kv_d.value <- 0.000003
Kv_Gk.value <- 0.000009
KTNFa_m.value <- 0.000009
KIFNgk_m.value <- 0.000009
Km_IL12.value <- 0.000009
KIL12_Th0.value <- 0.00009 
Km_Th0.value <- 0.00009
KTh0_Th1.value <- 0.00009 
KTh0_Th17.value <- 0.00009
KIL23_Th0.value <- 0.00009 
KIL2_Th0.value <- 0.000009 
KTh1_IL2.value <- 0.000009
KIFNgk_IFNg.value <- 0.00009 
KTh0_IFNg.value <- .00009 
KTh1_IFNg.value <- .00009
KTh17_IFNg.value <- .00009
KIFNg_IFNgk.value <- .00009
Km_TNFa.value <- .00009
Kd_IL23.value <- .0009
KTh17_IL17.value <- .00009
KIL17_Gk.value <- .00009
sigma_m.value <- 0.00023
sigma_IL12.value <- 0.0008
sigma_Th0_IL12.value <- 0.00023
sigma_Th0_IL17.value <- 0.00023
sigma_Th1.value <- 0.0008 #0.8
sigma_IL2.value <- 0.0008
sigma_IFNg.value <- 0.0008
sigma_IFNgk.value <- 0.0008
sigma_TNFa.value <- 0.0008
sigma_d.value <- 0.00023
sigma_IL23.value <- 0.0008
sigma_Th17.value <- 0.0008
sigma_IL17.value <- 0.0008
sigma_Gk.value <- 0.00023

parameter.list <- c( KTh0_IL17_IFNg = KTh0_IL17_IFNg.value,Kv_m = Kv_m.value, Kv_d = Kv_d.value, Kv_Gk = Kv_Gk.value, KTNFa_m = KTNFa_m.value, 
                    KIFNgk_m = KIFNgk_m.value,Km_IL12 = Km_IL12.value, KIL12_Th0 = KIL12_Th0.value,Km_Th0 = Km_Th0.value, 
                    KTh0_Th1 = KTh0_Th1.value, KTh0_Th17 = KTh0_Th17.value, KIL23_Th0 = KIL23_Th0.value,
                    KIL2_Th0 = KIL2_Th0.value, KTh1_IL2 = KTh1_IL2.value, KIFNgk_IFNg = KIFNgk_IFNg.value, 
                    KTh0_IL12_IFNg = KTh0_IL12_IFNg.value, KTh1_IFNg = KTh1_IFNg.value, KTh17_IFNg = KTh17_IFNg.value, 
                    KIFNg_IFNgk = KIFNg_IFNgk.value, Km_TNFa = Km_TNFa.value, Kd_IL23 = Kd_IL23.value, 
                    KTh17_IL17 = KTh17_IL17.value, KIL17_Gk = KIL17_Gk.value, sigma_m = sigma_m.value,
                    sigma_IL12 = sigma_IL12.value, sigma_Th0_IL12 = sigma_Th0_IL12.value, sigma_Th1 = sigma_Th1.value,
                    sigma_IL2 = sigma_IL2.value, sigma_IFNg = sigma_IFNg.value, sigma_IFNgk = sigma_IFNgk.value,
                    sigma_TNFa = sigma_TNFa.value, sigma_d = sigma_d.value, sigma_IL23 = sigma_IL23.value,
                    sigma_Th17 = sigma_Th17.value, sigma_IL17 = sigma_IL17.value, sigma_Gk = sigma_Gk.value,sigma_Th0_IL17 = sigma_Th0_IL17.value)

# initial values
Vaccine0 = 0.02
Makrofag0 = 5
IL120 = 0
Th0_IL12_0 = 2
Th0_IL17_0 = 2
Th10 = 0
IL20 = 0
IFNg0 = 0
IFNgk0 = 0
TNFa0 = 0
Dendrit0 = .5
IL230 = 0
Th170 = 0
IL170 = 0
GranulotcytKnoglemarv0 = 0

initial.values <- c(Vaccine = Vaccine0, Makrofag = Makrofag0, IL12 = IL120, Th0_IL12 = Th0_IL12_0, Th0_IL17 = Th0_IL17_0, Th1 = Th10, IL2 = IL20, IFNg = IFNg0, 
                    IFNgk = IFNgk0, TNFa = TNFa0, Dendrit = Dendrit0, IL23 = IL230, Th17 = Th170, IL17 = IL170, 
                    GranulotcytKnoglemarv = GranulotcytKnoglemarv0)

# Output timepoints
time.points <- seq(0,2*365,by=1)

# simulate the epidemic
output <- ode(y=initial.values,times = time.points,func = model, parms = parameter.list)

# Plot the result 
plot(Vaccine~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,0.02),xlim=c(0,10),ylab='Amount',xlab='Time (days)')
plot(Makrofag~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,10),ylab='Amount',xlab='Time (days)')
plot(IL12~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,0.1),ylab='Amount',xlab='Time (days)')
plot(Th0_IL12~time,data=output,type='l',lwd=3,lty=2,col='black',xlim=c(0,10),ylim=c(0,4),ylab='Amount',xlab='Time (days)')
plot(Th0_IL17~time,data=output,type='l',lwd=3,lty=2,col='black',xlim=c(0,10),ylim=c(0,4),ylab='Amount',xlab='Time (days)')
plot(Th1~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,0.1),ylab='Amount',xlab='Time (days)')
plot(IL2~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,0.001),ylab='Amount',xlab='Time (days)')
plot(IFNg~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,0.2),ylab='Amount',xlab='Time (days)')
plot(IFNgk~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,0.1),ylab='Amount',xlab='Time (days)')
plot(TNFa~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,0.6),ylab='Amount',xlab='Time (days)')

plot(Dendrit~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,6),xlim=c(0,6),ylab='Amount',xlab='Time (days)')
plot(IL23~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,0.2),ylab='Amount',xlab='Time (days)')
plot(Th17~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,0.1),ylab='Amount',xlab='Time (days)')
plot(IL17~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,0.01),ylab='Amount',xlab='Time (days)')
plot(GranulotcytKnoglemarv~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,0.1),ylab='Amount',xlab='Time (days)')

# Herfra loades data
DATA <- read.table("DataUdenHuller.csv", header=TRUE, sep=";", as.is=TRUE)

VaccineDATA = as.numeric(gsub(",", ".", gsub("\\.", "", DATA$H56.vacc..Dose)))
RestimulerendeVaccineData = as.numeric(gsub(",", ".", gsub("\\.", "", DATA$H56.restim..ug.ml.)))
IL17DATA = as.numeric(gsub(",", ".", gsub("\\.", "", DATA$mean...IL.17)))
INFgDATA = as.numeric(gsub(",", ".", gsub("\\.", "", DATA$mean...IFNg)))
IL2DATA= as.numeric(gsub(",", ".", gsub("\\.", "", DATA$mean...IL.2)))
INFaDATA= as.numeric(gsub(",", ".", gsub("\\.", "", DATA$mean...TNFa)))


plot(VaccineDATA,IL17DATA,xlim=c(0,50))
plot(RestimulerendeVaccineData,IL17DATA,xlim=c(0,100))

CelleData <- read.table("Celler.csv", header=TRUE, sep=";", as.is=TRUE)
plot(CelleData$Time..days,CelleData$Response..SFU.mill.splenocytes.)




dose0CelleData <- CelleData[1:45,c('Time..days.','Response..SFU.mill.splenocytes.')]
dose01CelleData <- CelleData[46:90,c('Time..days.','Response..SFU.mill.splenocytes.')]

plot(dose01CelleData$Time..days,dose01CelleData$Response..SFU.mill.splenocytes.)



