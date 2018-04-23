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
  #kM = state_values[16]
  
  
 # stop = 
  
  with ( 
    as.list (parameters),     # variable names within parameters can be used 
    {
      
      
      Kv_m <- 5*(250-Vaccine0)/((250-Vaccine0)+1)-(250-Vaccine0)*0.02
      Kv_d <- 5*Vaccine0/(Vaccine0+1)-Vaccine0*0.02
      #Kv_Gk <- ifelse(- Vaccine*Kv_m - Vaccine*Kv_d - Vaccine*Kv_Gk < 0,Kv_Gk.value,0)
      
      dVaccine = - Vaccine*Kv_m - Vaccine*Kv_d #- Vaccine*Kv_Gk
      
      
      # = ifelse(Vaccine < 0.00, - Vaccine*Kv_m - Vaccine*Kv_d - Vaccine*Kv_Gk, 0)
      
      
      dMakrofag = Vaccine*(Kv_m) + TNFa * KTNFa_m + IFNgk * KIFNgk_m - Makrofag*sigma_m
      dIL12 = Makrofag * Km_IL12 - IL12 * KIL12_Th0 - IL12*sigma_IL12
      
      dTh0_IL12 = Makrofag * Km_Th0 + IL12 * KIL12_Th0 - Th0_IL12 * KTh0_Th1 + IL2 * KIL2_Th0  - Th0_IL12*sigma_Th0_IL12
      dTh1 = Th0_IL12 * KTh0_Th1 - Th1*sigma_Th1
      dIL2 = Th1 * KTh1_IL2 - IL2*sigma_IL2
      dIFNg = Th1 * KTh1_IFNg + Th17 * KTh17_IFNg + Th0_IL12 * KTh0_IL12_IFNg + IFNgk * KIFNgk_IFNg - IFNg * KIFNg_IFNgk + Th0_IL17 * KTh0_IL17_IFNg - IFNg*sigma_IFNg
      dIFNgk = - IFNgk * KIFNgk_IFNg + IFNg * KIFNg_IFNgk - IFNgk*sigma_IFNgk
      dTNFa = Makrofag * Km_TNFa - TNFa*sigma_TNFa
      
      dTh0_IL17 = IL23 * KIL23_Th0 - Th0_IL17 * KTh0_Th17 - Th0_IL17*sigma_Th0_IL17
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


# 0.00000041 på 5 dag
 # 0.00000041
# parameters

KTh0_IL17_IFNg.value <- 0.00009
KTh0_IL12_IFNg.value <- 0.00009
#Kv_m.value <- 0.02
#Kv_d.value <- 0.000003
Kv_Gk.value <- 0.00009 # Ændre denne værdi for dose elimination
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

# initial values
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

# Output timepoints
time.points <- seq(0,2*365,by=0.01)

# simulate the epidemic
output <- ode(y=initial.values,times = time.points,func = model, parms = parameter.list)

# Plot the result 

# For at se mange decimaler kan den kommando nedenunder bruges. 
# Hvortil der ses, at vaccinen assymptotisk nærmer sig 0, 
# selvom den svinger mellem negative og positive værdier

#options("scipen"=100, "digits"=4) 

graphics.off()
# Plot af vaccine, makrofag og dendrit
plot(Vaccine~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,150),xlim=c(0,5),ylab='Amount',xlab='Time (days)')
lines(Makrofag~time,data=output,type='l',lwd=3,lty=2,col='red',ylim=c(0,150),xlim=c(0,5),ylab='Amount',xlab='Time (days)')
lines(Dendrit~time,data=output,type='l',lwd=3,lty=2,col='blue',ylim=c(0,150),xlim=c(0,5),ylab='Amount',xlab='Time (days)')
legend(3, 150, legend=c("Vaccine", "Makrofag","Dendrit"),
       col=c("black", "red","blue"), lty=1:4, cex=0.8)
title("Dose = 150")


# Plots af IL12 pathway
plot(IL12~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,5),xlim=c(0,2*365),ylab='Amount',xlab='Time (days)')
lines(Th0_IL12~time,data=output,type='l',lwd=3,lty=2,col='blue',ylim=c(0,2),xlim=c(0,2*365),ylab='Amount',xlab='Time (days)')
lines(Th1~time,data=output,type='l',lwd=3,lty=2,col='red',ylim=c(0,0.1),ylab='Amount',xlab='Time (days)')
lines(IL2~time,data=output,type='l',lwd=3,lty=2,col='green',ylim=c(0,0.0005),ylab='Amount',xlab='Time (days)')
lines(TNFa~time,data=output,type='l',lwd=3,lty=2,col='orange',ylim=c(0,0.6),ylab='Amount',xlab='Time (days)')
legend(400, 5, legend=c("IL12", "Th0 i IL12 pathway","Th1","IL2","TNF-alpha"),
       col=c("black", "blue","red","green","orange"), lty=1:4, cex=0.8)
title("Dose = 150")


plot(IFNg~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,0.3),ylab='Amount',xlab='Time (days)')
lines(IFNgk~time,data=output,type='l',lwd=3,lty=2,col='blue',ylim=c(0,0.3),ylab='Amount',xlab='Time (days)')
legend(400, 0.30, legend=c("IFNgamma", "IFNgammakompleks"),
       col=c("black", "blue"), lty=1:4, cex=0.8)
title("Dose = 150")


# Plot af IL17
plot(IL23~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,3),ylab='Amount',xlab='Time (days)')
lines(Th0_IL17~time,data=output,type='l',lwd=3,lty=2,col='blue',ylim=c(0,3),ylab='Amount',xlab='Time (days)')
lines(Th17~time,data=output,type='l',lwd=3,lty=2,col='red',ylim=c(0,0.1),ylab='Amount',xlab='Time (days)')
lines(IL17~time,data=output,type='l',lwd=3,lty=2,col='green',ylim=c(0,0.001),ylab='Amount',xlab='Time (days)')
lines(GranulotcytKnoglemarv~time,data=output,type='l',lwd=3,lty=2,col='orange',ylim=c(0,0.001),ylab='Amount',xlab='Time (days)')
legend(350, 3, legend=c("IL23", "Th0 i IL17 pathway","Th17","IL17","Granulocyt i knoglemarv"),
       col=c("black", "blue","red","green","orange"), lty=1:4, cex=0.8)
title("Dose = 150")

SimulatedDataDose100 <- as.data.frame(output)
Vaccinetest <- SimulatedDataDose100$Vaccine



plot(output)

# Sensitivtets analyse starter her
N.iter <- 100
KTh0_IL17_IFNg.list <- runif(N.iter,min=0.7,max=1.2)
KTh0_IL12_IFNg.list <- runif(N.iter,min=0.7,max=1.2)
Kv_m.list <- runif(N.iter,min=0.01,max=0.05)
Kv_d.list <- runif(N.iter,min=0.001,max=0.04) 
Kv_Gk.list <- runif(N.iter,min=0.005,max=0.05) 
KTNFa_m.list <- runif(N.iter,min=0.0005,max=0.004) 
KIFNgk_m.list <- runif(N.iter,min=50,max=150)
Km_IL12.list <- runif(N.iter,min=0.06,max=0.11)
KIL12_Th0.list <- runif(N.iter,min=0.75,max=1.2) 
Km_Th0.list <- runif(N.iter,min=0.75,max=1.2)
KTh0_Th1.list <- runif(N.iter,min=0.75,max=1.2) 
KTh0_Th17.list <- runif(N.iter,min=0.75,max=1.2)
KIL23_Th0.list <- runif(N.iter,min=0.75,max=1.2) 
KIL2_Th0.list <- runif(N.iter,min=0.075,max=0.12) 
KTh1_IL2.list <- runif(N.iter,min=0.075,max=0.12)
KIFNgk_IFNg.list <- runif(N.iter,min=0.75,max=1.2) 
KTh0_IFNg.list <- runif(N.iter,min=0.75,max=1.2)
KTh1_IFNg.list <- runif(N.iter,min=0.75,max=1.2)
KTh17_IFNg.list <- runif(N.iter,min=0.75,max=1.2)
KIFNg_IFNgk.list <- runif(N.iter,min=0.75,max=1.2)
Km_TNFa.list <- runif(N.iter,min=0.75,max=1.2)
Kd_IL23.list <- runif(N.iter,min=0.075,max=0.12)
KTh17_IL17.list <- runif(N.iter,min=0.75,max=1.2)
KIL17_Gk.list <- runif(N.iter,min=0.75,max=1.2)
sigma_m.list <- runif(N.iter,min=1,max=5)
sigma_IL12.list <- runif(N.iter,min=0.06,max=0.1)
sigma_Th0_IL12.list <- runif(N.iter,min=0.003,max=0.043)
sigma_Th0_IL17.list <- runif(N.iter,min=0.003,max=0.043)
sigma_Th1.list <- runif(N.iter,min=0.06,max=0.1)
sigma_IL2.list <- runif(N.iter,min=0.06,max=0.1)
sigma_IFNg.list <- runif(N.iter,min=0.06,max=0.1)
sigma_IFNgk.list <- runif(N.iter,min=0.06,max=0.1)
sigma_TNFa.list <- runif(N.iter,min=0.06,max=0.1)
sigma_d.list <- runif(N.iter,min=0.003,max=0.043)
sigma_IL23.list <- runif(N.iter,min=0.06,max=0.1)
sigma_Th17.list <- runif(N.iter,min=0.06,max=0.1)
sigma_IL17.list <- runif(N.iter,min=0.0006,max=0.1)
sigma_Gk.list <- runif(N.iter,min=0.003,max=0.043)
KM_TNFa.list <- runif(N.iter,min=5,max=90) 
resultsIL17 <- numeric(N.iter)
resultsIL17.t <-  numeric(N.iter) 

for ( i in 1:N.iter){
  # simulate the epidemic
  parameter.list <- c(KTh0_IL17_IFNg = KTh0_IL17_IFNg.list[i],Kv_m = Kv_m.list[i], Kv_d = Kv_d.list[i], Kv_Gk = Kv_Gk.list[i], KTNFa_m = KTNFa_m.list[i], 
                      KIFNgk_m = KIFNgk_m.list[i],Km_IL12 = Km_IL12.list[i], KIL12_Th0 = KIL12_Th0.list[i],Km_Th0 = Km_Th0.list[i], 
                      KTh0_Th1 = KTh0_Th1.list[i], KTh0_Th17 = KTh0_Th17.list[i], KIL23_Th0 = KIL23_Th0.list[i],
                      KIL2_Th0 = KIL2_Th0.list[i], KTh1_IL2 = KTh1_IL2.list[i], KIFNgk_IFNg = KIFNgk_IFNg.list[i], 
                      KTh0_IL12_IFNg = KTh0_IL12_IFNg.list[i], KTh1_IFNg = KTh1_IFNg.list[i], KTh17_IFNg = KTh17_IFNg.list[i], 
                      KIFNg_IFNgk = KIFNg_IFNgk.list[i], Km_TNFa = Km_TNFa.list[i], Kd_IL23 = Kd_IL23.list[i], 
                      KTh17_IL17 = KTh17_IL17.list[i], KIL17_Gk = KIL17_Gk.list[i], sigma_m = sigma_m.list[i],
                      sigma_IL12 = sigma_IL12.list[i], sigma_Th0_IL12 = sigma_Th0_IL12.list[i], sigma_Th1 = sigma_Th1.list[i],
                      sigma_IL2 = sigma_IL2.list[i], sigma_IFNg = sigma_IFNg.list[i], sigma_IFNgk = sigma_IFNgk.list[i],
                      sigma_TNFa = sigma_TNFa.list[i], sigma_d = sigma_d.list[i], sigma_IL23 = sigma_IL23.list[i],
                      sigma_Th17 = sigma_Th17.list[i], sigma_IL17 = sigma_IL17.list[i], sigma_Gk = sigma_Gk.list[i],
                      KM_TNFa = KM_TNFa.list[i], sigma_Th0_IL17 = sigma_Th0_IL17.list[i])
  output <- ode(y=initial.values,times = time.points,func = model, parms = parameter.list)
  
  resultsIL17[i] <-  max(output[,'IL17'])
  resultsIL17.t[i] <- output[which(output[,'IL17'] == resultsIL17[i]),'time']
  
}

#plot(resultsIL17.t,resultsIL17)

resultatIL17 <- data.frame(sigma_IL17 = sigma_IL17.list,KTh17_IL17 = KTh17_IL17.list,resultsIL17, resultsIL17.t)
pairs(resultatIL17)



epi.prcc(resultatIL17[,c(1,2,3)]) # beta/gamma er positiv/negativ korreleret med i.max
epi.prcc(resultatIL17[,c(1,2,4)]) 

# Scatterplots

pairs(m.res1)
#pairs(m.res2)

# PRCC

library(epiR)

epi.prcc(resultatIL17[,c(1,2,3)]) # beta/gamma er positiv/negativ korreleret med i.max

epi.prcc(resultatIL17[,c(1,2,4)]) # beta/gamma er negativ/svagt positiv korreleret med t.max

epi.prcc(resultatIL17[,c(1,2,5)])




m.res1 <- data.frame(Kv_m = Kv_m.list, KTNFa_m = KTNFa_m.list, 
                     KIFNgk_m = KIFNgk_m.list,Km_IL12 = Km_IL12.list, KIL12_Th0 = KIL12_Th0.list,Km_Th0 = Km_Th0.list, 
                     KTh0_Th1 = KTh0_Th1.list, 
                     KIL2_Th0 = KIL2_Th0.list, KTh1_IL2 = KTh1_IL2.list, KIFNgk_IFNg = KIFNgk_IFNg.list, 
                     KTh0_IL12_IFNg = KTh0_IL12_IFNg.list, KTh1_IFNg = KTh1_IFNg.list, 
                     KIFNg_IFNgk = KIFNg_IFNgk.list, Km_TNFa = Km_TNFa.list,  
                     sigma_m = sigma_m.list,
                     sigma_IL12 = sigma_IL12.list, sigma_Th0_IL12 = sigma_Th0_IL12.list, sigma_Th1 = sigma_Th1.list,
                     sigma_IL2 = sigma_IL2.list, sigma_IFNg = sigma_IFNg.list, sigma_IFNgk = sigma_IFNgk.list,
                     sigma_TNFa = sigma_TNFa.list, KM_TNFa = KM_TNFa.list,
                     KTh0_IL17_IFNg = KTh0_IL17_IFNg.list, Kv_d = Kv_d.list, Kv_Gk = Kv_Gk.list,  
                     KTh0_Th17 = KTh0_Th17.list, KIL23_Th0 = KIL23_Th0.list,
                     KTh17_IFNg = KTh17_IFNg.list, 
                     Kd_IL23 = Kd_IL23.list, 
                     KTh17_IL17 = KTh17_IL17.list, KIL17_Gk = KIL17_Gk.list, 
                     sigma_IFNg = sigma_IFNg.list, 
                     sigma_d = sigma_d.list, sigma_IL23 = sigma_IL23.list,
                     sigma_Th17 = sigma_Th17.list, sigma_IL17 = sigma_IL17.list, sigma_Gk = sigma_Gk.list,
                     sigma_Th0_IL17 = sigma_Th0_IL17.list, resultsIL17, resultsIL17.t)


# Scatterplots

pairs(m.res1)
#pairs(m.res2)

# PRCC

library(epiR)

epi.prcc(m.res1[,c(1,2,3)]) # beta/gamma er positiv/negativ korreleret med i.max

epi.prcc(m.res1[,c(1,2,4)]) # beta/gamma er negativ/svagt positiv korreleret med t.max

epi.prcc(m.res1[,c(1,2,5)])



