# TB model in R

rm(list=ls())

# Henter ode solveren
library(deSolve)


# Laver funktionen med de 15 parametre
model = function (current_timepoint, state_values, parameters)
{
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
  
  with ( 
    as.list (parameters),     # variable names within parameters can be used 
    {
      
      # De følgende 2 konstanter bestemmes vha. følgende funktion, 
      # som er afhængig af startdosen.
      Kv_m <- 5*(250-Vaccine0)/((250-Vaccine0)+1)-(250-Vaccine0)*0.02
      Kv_d <- 5*Vaccine0/(Vaccine0+1)-Vaccine0*0.02
      
      dVaccine = - Vaccine*Kv_m - Vaccine*Kv_d #- Vaccine*Kv_Gk
    
      # IL12 pathway
      dMakrofag = Vaccine*(Kv_m) + TNFa * KTNFa_m + IFNgk * KIFNgk_m - Makrofag*sigma_m
      dIL12 = Makrofag * Km_IL12 - IL12 * KIL12_Th0 - IL12*sigma_IL12
      dTh0_IL12 = Makrofag * Km_Th0 + IL12 * KIL12_Th0 - Th0_IL12 * KTh0_Th1 + IL2 * KIL2_Th0  - Th0_IL12*sigma_Th0_IL12
      dTh1 = Th0_IL12 * KTh0_Th1 - Th1*sigma_Th1
      dIL2 = Th1 * KTh1_IL2 - IL2*sigma_IL2
      dIFNg = Th1 * KTh1_IFNg + Th17 * KTh17_IFNg + Th0_IL12 * KTh0_IL12_IFNg + IFNgk * KIFNgk_IFNg - IFNg * KIFNg_IFNgk + Th0_IL17 * KTh0_IL17_IFNg - IFNg*sigma_IFNg
      dIFNgk = - IFNgk * KIFNgk_IFNg + IFNg * KIFNg_IFNgk - IFNgk*sigma_IFNgk
      dTNFa = Makrofag * Km_TNFa - TNFa*sigma_TNFa
      
      # IL17 pathway
      dDendrit = Vaccine*Kv_d - Dendrit*sigma_d
      dTh0_IL17 = IL23 * KIL23_Th0 - Th0_IL17 * KTh0_Th17 - Th0_IL17*sigma_Th0_IL17
      dDendrit = Vaccine*Kv_d - Dendrit*sigma_d
      dIL23 = Dendrit * Kd_IL23 - IL23*sigma_IL23
      dTh17 = Th0_IL17 * KTh0_Th17 - Th17*sigma_Th17
      dIL17 = Th17 * KTh17_IL17 - IL17*sigma_IL17
      dGranulocytKnoglemarv =  Vaccine*Kv_Gk + IL17 * KIL17_Gk - GranulotcytKnoglemarv*sigma_Gk 
      
      # sætter resultaterne sammen
      results = c(dVaccine, dMakrofag, dIL12, dTh0_IL12, dTh0_IL17, dTh1, dIL2, dIFNg, dIFNgk, dTNFa, dDendrit,
                  dIL23, dTh17, dIL17, dGranulocytKnoglemarv)
      list (results)
    }
  )
}


# Parameter og deres værdier
KTh0_IL17_IFNg.value <- 0.00009
KTh0_IL12_IFNg.value <- 0.00009
Kv_Gk.value <- 0.00009 
KTNFa_m.value <- 0.00009
KIFNgk_m.value <- 0.00009
Km_IL12.value <- 0.009
KIL12_Th0.value <- 0.0009 
Km_Th0.value <- 0.00009
KTh0_Th1.value <- 0.009 
KTh0_Th17.value <- 0.009
KIL23_Th0.value <- 0.00009 
KIL2_Th0.value <- 0.009 
KTh1_IL2.value <- 0.09
KIFNgk_IFNg.value <- 0.00009 
KTh0_IFNg.value <- .00009 
KTh1_IFNg.value <- .00009
KTh17_IFNg.value <- .009
KIFNg_IFNgk.value <- .04
Km_TNFa.value <- .0009
Kd_IL23.value <- 0.0009
KTh17_IL17.value <- .005
KIL17_Gk.value <- .009
sigma_m.value <- 0.23
sigma_IL12.value <- 0.0023
sigma_Th0_IL12.value <- 0.0023
sigma_Th0_IL17.value <- 0.0023
sigma_Th1.value <- 0.023 #0.8
sigma_IL2.value <- 0.053
sigma_IFNg.value <- 0.023
sigma_IFNgk.value <- 0.023
sigma_TNFa.value <- 0.0023
sigma_d.value <- 0.23
sigma_IL23.value <- 0.0023
sigma_Th17.value <- 0.0023
sigma_IL17.value <- 0.0053
sigma_Gk.value <- 0.023

parameter.list <- c(KTh0_IL17_IFNg = KTh0_IL17_IFNg.value, Kv_Gk = Kv_Gk.value, KTNFa_m = KTNFa_m.value, 
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

# initial værdier
Vaccine0 = 150
Makrofag0 = 5
IL120 = 0
Th0_IL12_0 = 2
Th0_IL17_0 = 2
Th10 = 0
IL20 = 0
IFNg0 = 0
IFNgk0 = 0
TNFa0 = 0
Dendrit0 = 5
IL230 = 0
Th170 = 0
IL170 = 0
GranulotcytKnoglemarv0 = 0

initial.values <- c(Vaccine = Vaccine0, Makrofag = Makrofag0, IL12 = IL120, Th0_IL12 = Th0_IL12_0, Th0_IL17 = Th0_IL17_0, Th1 = Th10, IL2 = IL20, IFNg = IFNg0, 
                    IFNgk = IFNgk0, TNFa = TNFa0, Dendrit = Dendrit0, IL23 = IL230, Th17 = Th170, IL17 = IL170, 
                    GranulotcytKnoglemarv = GranulotcytKnoglemarv0)

# Tiden vi simulerer i
time.points <- seq(0,2*365,by=0.01)

# Simulering af modellen
output <- ode(y=initial.values,times = time.points,func = model, parms = parameter.list)

# Plot af resultaterne
graphics.off()

# Plot af vaccine, makrofag og dendrit
plot(Vaccine~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,150),xlim=c(0,5),ylab='Mængde',xlab='Tid (dage)')
lines(Makrofag~time,data=output,type='l',lwd=3,lty=2,col='red')
lines(Dendrit~time,data=output,type='l',lwd=3,lty=2,col='blue')
legend(3, 150, legend=c("Vaccine", "Makrofag","Dendrit"),
       col=c("black", "red","blue"), lty=1:4, cex=0.8)
title("Dosis = 150")


# Plot af IL12 pathway
plot(IL12~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,5),xlim=c(0,2*365),ylab='Mængde',xlab='Tid (dage)')
lines(Th0_IL12~time,data=output,type='l',lwd=3,lty=2,col='blue')
lines(Th1~time,data=output,type='l',lwd=3,lty=2,col='red')
lines(IL2~time,data=output,type='l',lwd=3,lty=2,col='green')
lines(TNFa~time,data=output,type='l',lwd=3,lty=2,col='orange')
legend(400, 5, legend=c("IL12", "Th0 i IL12 pathway","Th1","IL2","TNF-alpha"),
       col=c("black", "blue","red","green","orange"), lty=1:4, cex=0.8)
title("Dosis = 150")


plot(IFNg~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,0.3),ylab='Mængde',xlab='Tid (dage)')
lines(IFNgk~time,data=output,type='l',lwd=3,lty=2,col='blue')
legend(400, 0.30, legend=c("IFNgamma", "IFNgammakompleks"),
       col=c("black", "blue"), lty=1:4, cex=0.8)
title("Dosis = 150")


# Plot af IL17
plot(IL23~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,3),ylab='Mængde',xlab='Time (dage)')
lines(Th0_IL17~time,data=output,type='l',lwd=3,lty=2,col='blue',ylim=c(0,3))
lines(Th17~time,data=output,type='l',lwd=3,lty=2,col='red',ylim=c(0,0.1))
lines(IL17~time,data=output,type='l',lwd=3,lty=2,col='green',ylim=c(0,0.001))
lines(GranulotcytKnoglemarv~time,data=output,type='l',lwd=3,lty=2,col='orange')
legend(350, 3, legend=c("IL23", "Th0 i IL17 pathway","Th17","IL17","Granulocyt i knoglemarv"),
       col=c("black", "blue","red","green","orange"), lty=1:4, cex=0.8)
title("Dosis = 150")


