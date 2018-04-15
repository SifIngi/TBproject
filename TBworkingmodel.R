# TB model in R

rm(list=ls())

# get the ODE solver
library(deSolve)

# make the function of the system
model = function (current_timepoint, state_values, parameters)
{
  # create state variables (local variables)
  #Vaccine = state_values [1]       
  Makrofag = state_values [1]
  IL12 = state_values [2]
  Th0 = state_values [3]
  Th1 = state_values [4]
  IL2 = state_values [5]
  IFNg = state_values [6]
  IFNgk = state_values [7]
  TNFa = state_values [8]
  Dendrit = state_values [9]
  IL23 = state_values [10]
  Th17 = state_values [11]
  IL17 = state_values [12]
  GranulotcytKnoglemarv = state_values[13]
  
  with ( 
    as.list (parameters),     # variable names within parameters can be used 
    {
      #dVaccine = - Vaccine * Kv_m - Vaccine * Kv_d - Vaccine * Kv_Gk
      dMakrofag = Vaccine * Kv_m + TNFa * KTNFa_m + IFNgk * KIFNgk_m - Makrofag*sigma_m
      dIL12 = Makrofag * Km_IL12 - IL12 * KIL12_Th0 - IL12*sigma_IL12
      dTh0 = Makrofag * Km_Th0 + IL12 * KIL12_Th0 - Th0 * KTh0_Th1 + IL2 * KIL2_Th0 + IL23 * KIL23_Th0 - Th0 * KTh0_Th17 - Th0*sigma_Th0
      dTh1 = Th0 * KTh0_Th1 - Th1*sigma_Th1
      dIL2 = Th1 * KTh1_IL2 - IL2*sigma_IL2
      dIFNg = Th1 * KTh1_IFNg + Th17 * KTh17_IFNg + Th0 * KTh0_IFNg + IFNgk * KIFNgk_IFNg - IFNg * KIFNg_IFNgk - IFNg*sigma_IFNg
      dIFNgk = - IFNgk * KIFNgk_IFNg + IFNg * KIFNg_IFNgk - IFNgk*sigma_IFNgk
      dTNFa = Makrofag * Km_TNFa - TNFa*sigma_TNFa
      dDendrit = Vaccine * Kv_d - Dendrit*sigma_d
      dIL23 = Dendrit * Kd_IL23 - IL23*sigma_IL23
      dTh17 = Th0 * KTh0_Th17 - Th17*sigma_Th17
      dIL17 = Th17 * KTh17_IL17 - IL17*sigma_IL17
      dGranulocytKnoglemarv = Vaccine * Kv_Gk + IL17 * KIL17_Gk - GranulotcytKnoglemarv*sigma_Gk 
      
      
      
      # combine results
      results = c(dMakrofag, dIL12, dTh0, dTh1, dIL2, dIFNg, dIFNgk, dTNFa, dDendrit,
                  dIL23, dTh17, dIL17, dGranulocytKnoglemarv)
      list (results)
    }
  )
}

# parameters
Vaccine.value <- 1
Kv_m.value <- 0.09
Kv_d.value <- 0.09
Kv_Gk.value <- 0.09
KTNFa_m.value <- 0.09
KIFNgk_m.value <- 0.09
Km_IL12.value <- 0.09
KIL12_Th0.value <- 0.09 
Km_Th0.value <- 0.09
KTh0_Th1.value <- 0.9 
KTh0_Th17.value <- 0.09
KIL23_Th0.value <- 0.090 
KIL2_Th0.value <- 0.09 
KTh1_IL2.value <- 0.09
KIFNgk_IFNg.value <- 0.09 
KTh0_IFNg.value <- .09 
KTh1_IFNg.value <- .09
KTh17_IFNg.value <- .09
KIFNg_IFNgk.value <- .09
Km_TNFa.value <- .09
Kd_IL23.value <- .09
KTh17_IL17.value <- .09
KIL17_Gk.value <- .09
sigma_m.value <- 8
sigma_IL12.value <- 0.8
sigma_Th0.value <- 10
sigma_Th1.value <- 0 #0.8
sigma_IL2.value <- 0.8
sigma_IFNg.value <- 0.8
sigma_IFNgk.value <- 0.8
sigma_TNFa.value <- 0.8
sigma_d.value <- 0.8
sigma_IL23.value <- 0.8
sigma_Th17.value <- 0.8
sigma_IL17.value <- 0.8
sigma_Gk.value <- 0.8

parameter.list <- c(Vaccine = Vaccine.value, Kv_m = Kv_m.value, Kv_d = Kv_d.value, Kv_Gk = Kv_Gk.value, KTNFa_m = KTNFa_m.value, 
                    KIFNgk_m = KIFNgk_m.value,Km_IL12 = Km_IL12.value, KIL12_Th0 = KIL12_Th0.value,Km_Th0 = Km_Th0.value, 
                    KTh0_Th1 = KTh0_Th1.value, KTh0_Th17 = KTh0_Th17.value, KIL23_Th0 = KIL23_Th0.value,
                    KIL2_Th0 = KIL2_Th0.value, KTh1_IL2 = KTh1_IL2.value, KIFNgk_IFNg = KIFNgk_IFNg.value, 
                    KTh0_IFNg = KTh0_IFNg.value, KTh1_IFNg = KTh1_IFNg.value, KTh17_IFNg = KTh17_IFNg.value, 
                    KIFNg_IFNgk = KIFNg_IFNgk.value, Km_TNFa = Km_TNFa.value, Kd_IL23 = Kd_IL23.value, 
                    KTh17_IL17 = KTh17_IL17.value, KIL17_Gk = KIL17_Gk.value, sigma_m = sigma_m.value,
                    sigma_IL12 = sigma_IL12.value, sigma_Th0 = sigma_Th0.value, sigma_Th1 = sigma_Th1.value,
                    sigma_IL2 = sigma_IL2.value, sigma_IFNg = sigma_IFNg.value, sigma_IFNgk = sigma_IFNgk.value,
                    sigma_TNFa = sigma_TNFa.value, sigma_d = sigma_d.value, sigma_IL23 = sigma_IL23.value,
                    sigma_Th17 = sigma_Th17.value, sigma_IL17 = sigma_IL17.value, sigma_Gk = sigma_Gk.value)

# initial values
#Vaccine0 = 10     
Makrofag0 = 0
IL120 = 0
Th00 = 0
Th10 = 0
IL20 = 0
IFNg0 = 0
IFNgk0 = 0
TNFa0 = 0
Dendrit0 = 0
IL230 = 0
Th170 = 0
IL170 = 0
GranulotcytKnoglemarv0 = 0

initial.values <- c(Makrofag = Makrofag0, IL12 = IL120, Th0 = Th00, Th1 = Th10, IL2 = IL20, IFNg = IFNg0, 
                    IFNgk = IFNgk0, TNFa = TNFa0, Dendrit = Dendrit0, IL23 = IL230, Th17 = Th170, IL17 = IL170, 
                    GranulotcytKnoglemarv = GranulotcytKnoglemarv0)

# Output timepoints
time.points <- seq(0,10,by=1)

# simulate the epidemic
output <- ode(y=initial.values,times = time.points,func = model, parms = parameter.list)

# Plot the result 
plot(Vaccine~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,0.20),ylab='Amount',xlab='Time (days)')
plot(Makrofag~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,0.03),ylab='Amount',xlab='Time (days)')
plot(IL12~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,0.004),ylab='Amount',xlab='Time (days)')
plot(Th0~time,data=output,type='l',lwd=3,lty=2,col='black',xlim=c(0,0.10),ylim=c(0,0.00010),ylab='Amount',xlab='Time (days)')
plot(Th1~time,data=output,type='l',lwd=3,lty=2,col='black',xlim=c(0,0.1),ylim=c(0,0.01),ylab='Amount',xlab='Time (days)')
plot(IL2~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(IFNg~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(IFNgk~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(TNFa~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,100000000),ylab='Amount',xlab='Time (days)')
plot(Dendrit~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(IL23~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(Th17~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(IL17~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(GranulotcytKnoglemarv~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')

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



