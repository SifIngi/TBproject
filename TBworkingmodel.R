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
  Fagosom = state_values [3]      
  Fagolysosom = state_values [4] 
  Peptider = state_values [5] 
  MHC2 = state_values [6]
  IL12 = state_values [7]
  Th0 = state_values [8]
  Th1 = state_values [9]
  Il2 = state_values [10]
  Infg = state_values [11]
  InfgR = state_values [12]
  Granuloma = state_values [13]
  Tnfa = state_values [14]
  Dendrocyt = state_values [15]
  IL23 = state_values [16]
  Th17 = state_values [17]
  IL17 = state_values [18]
  GranulopoeticChemicalFactors = state_values[19]
  GranulotcytKnoglemarv = state_values[20]
  
  with ( 
    as.list (parameters),     # variable names within parameters can be used 
    {
      dVaccine = - Vaccine*Kv_m - Vaccine*Kv_d #+*Kgv_v  # Unit for vaccinen er mikrogram/mL, mikrogram/mL/dage,  Kvaccine = mikrogram/mL/dage
      dMakrofag = Vaccine*Kv_m - Makrofag*Km_f + Infg*Kinfgr_m + Tnfa*Ktnfa_m
      dFagosom = Makrofag*Km_f - Fagosom*Kf_fl 
      dFagolysosom = Fagosom*Kf_fl - Fagolysosom*Kfl_il12 - Fagolysosom*Kfl_p
      dPeptider = Fagolysosom*Kfl_p - Peptider*Kp_mhc2
      dMHC2 = Peptider*Kp_mhc2 - MHC2*Kmhc2_th0
      dIL12 = Fagolysosom*Kfl_il12 + InfgR*Kil12_infgr - IL12*Kil12_th0
      dTh0 = MHC2*Kmhc2_th0 + IL12*Kil12_th0 - Th0*Kth0_th1 - Th0*Kth0_il2 + IL23*Kil23_th0 - Th0*Kth0_th17 - Th0*Kth0_infg
      dTh1 = Th0*Kth0_th1 - Th1*Kth1_infg
      dIl2 = Th0*Kth0_il2 - Il2*Kil2_infg
      dInfg = Th1*Kth1_infg + Il2*Kil2_infg + Th17*Kth17_infg + Th0*Kth0_infg - Infg*Kinfg_infgr - Infg*Kinfg_g
      dInfgR = Infg*Kinfg_infgr - IL12*Kil12_infgr - InfgR*Kinfgr_m - InfgR*KinfgR_tnfa
      dGranuloma = Infg*Kinfg_g + Tnfa*Ktnfa_g
      dTnfa = InfgR*KinfgR_tnfa - Tnfa*Ktnfa_g - Tnfa*Ktnfa_m
      dDendrocyt = Vaccine*Kv_d - Dendrocyt*Kd_il23
      dIL23 = Dendrocyt*Kd_il23 - IL23*Kil23_th0
      dTh17 = Th0*Kth0_th17 - Th17*Kth17_infg - Th17*Kth17_il17
      dIL17 = Th17*Kth17_il17 - IL17*Kil17_gcf
      dGranulopoeticChemicalFactors = IL17*Kil17_gcf - GranulopoeticChemicalFactors*Kgcf_gk
      dGranulotcytKnoglemarv = GranulopoeticChemicalFactors*Kgcf_gk # - dGranulotcytKnoglemarv*Kgk_v
      
      
      # combine results
      results = c(dVaccine, dMakrofag, dFagosom, dFagolysosom, dPeptider, dMHC2, dIL12, dTh0, dTh1, dIl2, 
                  dInfg, dInfgR, dGranuloma, dTnfa, dDendrocyt, dIL23, dTh17, dIL17, dGranulopoeticChemicalFactors, 
                  dGranulotcytKnoglemarv)
      list (results)
    }
  )
}

# parameters
Kv_m.value <- 1
Km_f.value <- 1
Kf_fl.value <- 1
Kfl_il12.value <- 1
Kfl_p.value <- 1
Kp_mhc2.value <- 1
Kmhc2_th0.value <- 1
Kil12_infgr.value <- 1
Kil12_th0.value <- 1
Kth0_th1.value <- 1 #
Kth0_il2.value <- 1 #
Kinfgr_m.value <- 1 #
Ktnfa_m.value <- 1 #
Kth1_infg.value <- 1 #
Kil2_infg.value <- 1 #
Kth17_infg.value <- 1 #
Kth0_infg.value <- 1 #
Kinfg_infgr.value <- 1 #
Kinfg_g.value <- 1 #
KinfgR_tnfa.value <- 1 #
Ktnfa_g.value <- 1 #
Kv_d.value <- 1 #
Kd_il23.value <- 1
Kil23_th0.value <- 1
Kth0_th17.value <- 1
Kth17_il17.value <- 1
Kil17_gcf.value <- 1
Kgcf_gk.value <- 1
#Kgk_v.value <- 1

parameter.list <- c(Kv_m = Kv_m.value, Km_f = Km_f.value, Kf_fl = Kf_fl.value, Kfl_il12 = Kfl_il12.value, 
                    Kfl_p = Kfl_p.value, Kp_mhc2 = Kp_mhc2.value,Kmhc2_th0 = Kmhc2_th0.value, 
                    Kil12_infgr = Kil12_infgr.value, Kil12_th0 = Kil12_th0.value, Kth0_th1 = Kth0_th1.value,
                    Kth0_il2 = Kth0_il2.value, Kinfgr_m = Kinfgr_m.value, Ktnfa_m = Ktnfa_m.value, 
                    Kth1_infg = Kth1_infg.value, Kil2_infg = Kil2_infg.value, Kth17_infg = Kth17_infg.value, 
                    Kth0_infg = Kth0_infg.value, Kinfg_infgr = Kinfg_infgr.value, Kinfg_g = Kinfg_g.value, 
                    KinfgR_tnfa = KinfgR_tnfa.value, Ktnfa_g = Ktnfa_g.value, Kv_d = Kv_d.value, 
                    Kd_il23 = Kd_il23.value, Kil23_th0 = Kil23_th0.value, Kth0_th17 = Kth0_th17.value, Kth17_il17 = Kth17_il17.value,
                    Kil17_gcf = Kil17_gcf.value, Kgcf_gk = Kgcf_gk.value)

# initial values
Vaccine0 = 1       
Makrofag0 = 0
Fagosom0 = 0    
Fagolysosom0 = 0 
Peptider0 = 0 
MHC20 = 0
IL120 = 0
Th00 = 0
Th10 = 0
IL20 = 0
Infg0 = 0
InfgR0 = 0
Granuloma0 = 0
Tnfa0 = 0
Dendrocyt0 = 0
IL230 = 0
Th170 = 0
IL170 = 0
GranulopoeticChemicalFactors0 = 0
GranulotcytKnoglemarv0 = 0

initial.values <- c(Vaccine = Vaccine0, Makrofag = Makrofag0, Fagosom = Fagosom0, Fagolysosom = Fagolysosom0, Peptider = Peptider0, 
                    MHC2 = MHC20, IL12 = IL120, Th0 = Th00, Th1 = Th10, IL2 = IL20, Infg = Infg0, 
                    InfgR = InfgR0, Granuloma = Granuloma0, Tnfa = Tnfa0, Dendrocyt = Dendrocyt0, 
                    IL23 = IL230, Th17 = Th170, IL17 = IL170, 
                    GranulopoeticChemicalFactors = GranulopoeticChemicalFactors0, GranulotcytKnoglemarv = GranulotcytKnoglemarv0)

# Output timepoints
time.points <- seq(0,20,by=1)

# simulate the epidemic
output <- ode(y=initial.values,times = time.points,func = model, parms = parameter.list)

# Plot the result 
plot(Vaccine~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(Makrofag~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(Fagosom~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(Fagolysosom~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(Peptider~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(MHC2~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(IL12~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(Th0~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(Th1~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(Il2~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(Infg~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(InfgR~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(Granuloma~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(Tnfa~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(Dendrocyt~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(IL23~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(Th17~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(IL17~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(GranulopoeticChemicalFactors~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')
plot(GranulotcytKnoglemarv~time,data=output,type='l',lwd=3,lty=2,col='black',ylim=c(0,1),ylab='Amount',xlab='Time (days)')

# Hør Niels omkring sidste parameter (rødbrune streg)! 





