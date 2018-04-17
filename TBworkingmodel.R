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
      Kv_m <- ifelse(Vaccine0 > 0.02, Kv_m.value, 0.01)
      Kv_d <- ifelse(Vaccine0 > 0.02, Kv_d.value, 0.02)

      dVaccine = - (Vaccine*Kv_m) - (Vaccine*Kv_d) - Vaccine*Kv_Gk
      dMakrofag = (Vaccine*Kv_m) + TNFa * KTNFa_m + IFNgk * KIFNgk_m/(KM_TNFa+IFNgk) - Makrofag*sigma_m
      dIL12 = Makrofag * Km_IL12 - IL12 * KIL12_Th0 - IL12*sigma_IL12
      dTh0_IL12 = Makrofag * Km_Th0 + IL12 * KIL12_Th0 - Th0_IL12 * KTh0_Th1 + IL2 * KIL2_Th0  - Th0_IL12*sigma_Th0_IL12
      dTh0_IL17 = IL23 * KIL23_Th0 - Th0_IL17 * KTh0_Th17 - Th0_IL17*sigma_Th0_IL17
      dTh1 = Th0_IL12 * KTh0_Th1 - Th1*sigma_Th1
      dIL2 = Th1 * KTh1_IL2 - IL2*sigma_IL2
      dIFNg = Th1 * KTh1_IFNg + Th17 * KTh17_IFNg + Th0_IL12 * KTh0_IL12_IFNg + IFNgk * KIFNgk_IFNg - IFNg * KIFNg_IFNgk + Th0_IL17 * KTh0_IL17_IFNg - IFNg*sigma_IFNg
      dIFNgk = - IFNgk * KIFNgk_IFNg + IFNg * KIFNg_IFNgk - IFNgk*sigma_IFNgk
      dTNFa = Makrofag * Km_TNFa - TNFa*sigma_TNFa
      dDendrit = Vaccine*Kv_d - Dendrit*sigma_d
      dIL23 = Dendrit * Kd_IL23 - IL23*sigma_IL23
      dTh17 = Th0_IL17 * KTh0_Th17 - Th17*sigma_Th17
      dIL17 = Th17 * KTh17_IL17 - IL17*sigma_IL17
      dGranulocytKnoglemarv =  Vaccine*Kv_Gk + IL17 * KIL17_Gk - GranulotcytKnoglemarv*sigma_Gk 
      
      
      
      # combine results
      results = c(dVaccine, dMakrofag, dIL12, dTh0_IL12, dTh0_IL17, dTh1, dIL2, dIFNg, dIFNgk, dTNFa, dDendrit,
                  dIL23, dTh17, dIL17, dGranulocytKnoglemarv)
      list (results)
    }
  )
}

# parameters
KTh0_IL17_IFNg.value <- 0.9
KTh0_IL12_IFNg.value <- 0.9
Kv_m.value <- 0.03
Kv_d.value <- 0.02
Kv_Gk.value <- 0.025
KTNFa_m.value <- 0.002
KIFNgk_m.value <- 100
Km_IL12.value <- 0.09
KIL12_Th0.value <- 0.9 
Km_Th0.value <- 0.9
KTh0_Th1.value <- 0.9 
KTh0_Th17.value <- 0.9
KIL23_Th0.value <- 0.9 
KIL2_Th0.value <- 0.09 
KTh1_IL2.value <- 0.09
KIFNgk_IFNg.value <- 0.9 
KTh0_IFNg.value <- .9 
KTh1_IFNg.value <- .9
KTh17_IFNg.value <- .9
KIFNg_IFNgk.value <- .9
Km_TNFa.value <- .9
Kd_IL23.value <- .09
KTh17_IL17.value <- .9
KIL17_Gk.value <- .9
sigma_m.value <- 3
sigma_IL12.value <- 0.08
sigma_Th0_IL12.value <- 0.023
sigma_Th0_IL17.value <- 0.023
sigma_Th1.value <- 0.08 
sigma_IL2.value <- 0.08
sigma_IFNg.value <- 0.08
sigma_IFNgk.value <- 0.08
sigma_TNFa.value <- 0.08
sigma_d.value <- 0.023
sigma_IL23.value <- 0.08
sigma_Th17.value <- 0.08
sigma_IL17.value <- 0.08
sigma_Gk.value <- 0.023
KM_TNFa.value <- 1000

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
                    sigma_Th17 = sigma_Th17.value, sigma_IL17 = sigma_IL17.value, sigma_Gk = sigma_Gk.value,
                    KM_TNFa = KM_TNFa.value, sigma_Th0_IL17 = sigma_Th0_IL17.value)

# initial values
Vaccine0 = 0.02
Makrofag0 = 1
IL120 = 0
Th0_IL12_0 = 0
Th0_IL17_0 = 0
Th10 = 0
IL20 = 0
IFNg0 = 0
IFNgk0 = 0
TNFa0 = 0
Dendrit0 = 1
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
plot(Vaccine~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,0.03),xlim=c(0,50),ylab='Amount',xlab='Time (days)',main="Vaccine")
plot(Makrofag~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,40),ylab='Amount',xlab='Time (days)',main="Makrofag")
plot(IL12~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,4),ylab='Amount',xlab='Time (days)',main="IL12")
plot(Th0_IL12~time,data=output,type='l',lwd=3,lty=2,col='black',xlim=c(0,700),ylim=c(0,100000),ylab='Amount',xlab='Time (days)',main="Th0_IL12")
plot(Th0_IL17~time,data=output,type='l',lwd=3,lty=2,col='black',xlim=c(0,700),ylim=c(0,1.1),ylab='Amount',xlab='Time (days)',main="Th0_IL17")
plot(Th1~time,data=output,type='l',lwd=3,lty=2,col='black',xlim=c(0,700),ylim=c(0,1000000),ylab='Amount',xlab='Time (days)',main="Th1")
plot(IL2~time,data=output,type='l',lwd=3,lty=2,col='black',xlim=c(0,700),ylim=c(0,1000000),ylab='Amount',xlab='Time (days)',main="IL2")
plot(IFNg~time,data=output,type='l',lwd=3,lty=2,col='black',xlim=c(0,700), ylim=c(0,6000000),ylab='Amount',xlab='Time (days)',main="IFNg")
plot(IFNgk~time,data=output,type='l',lwd=3,lty=2,col='black',xlim=c(0,700),ylim=c(0,5000000),ylab='Amount',xlab='Time (days)',main="IFNgk")
plot(TNFa~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,450),ylab='Amount',xlab='Time (days)',main="TNFa")

plot(Dendrit~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1.5),xlim=c(0,100),ylab='Amount',xlab='Time (days)',main="Dendrit")
plot(IL23~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,0.8),ylab='Amount',xlab='Time (days)',main="IL23")
plot(Th17~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,8),ylab='Amount',xlab='Time (days)',main="Th17")
plot(IL17~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,80),ylab='Amount',xlab='Time (days)',main="IL17")
plot(GranulotcytKnoglemarv~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1800),ylab='Amount',xlab='Time (days)',main="GranulotcytKnoglemarv")

# Herfra loades data
# DATA <- read.table("DataUdenHuller.csv", header=TRUE, sep=";", as.is=TRUE)
# 
# VaccineDATA = as.numeric(gsub(",", ".", gsub("\\.", "", DATA$H56.vacc..Dose)))
# RestimulerendeVaccineData = as.numeric(gsub(",", ".", gsub("\\.", "", DATA$H56.restim..ug.ml.)))
# IL17DATA = as.numeric(gsub(",", ".", gsub("\\.", "", DATA$mean...IL.17)))
# INFgDATA = as.numeric(gsub(",", ".", gsub("\\.", "", DATA$mean...IFNg)))
# IL2DATA= as.numeric(gsub(",", ".", gsub("\\.", "", DATA$mean...IL.2)))
# INFaDATA= as.numeric(gsub(",", ".", gsub("\\.", "", DATA$mean...TNFa)))
# 
# 
# plot(VaccineDATA,IL17DATA,xlim=c(0,50))
# plot(RestimulerendeVaccineData,IL17DATA,xlim=c(0,1))
# 
# CelleData <- read.table("Celler.csv", header=TRUE, sep=";", as.is=TRUE)
# plot(CelleData$Time..days,CelleData$Response..SFU.mill.splenocytes.)
# 
# 
# 
# 
# dose0CelleData <- CelleData[1:45,c('Time..days.','Response..SFU.mill.splenocytes.')]
# dose01CelleData <- CelleData[46:90,c('Time..days.','Response..SFU.mill.splenocytes.')]
# 
# plot(dose01CelleData$Time..days,dose01CelleData$Response..SFU.mill.splenocytes.)



