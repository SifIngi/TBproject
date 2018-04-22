# plot af forskellige ting

Dose02 <- read.csv(file = "Data02", row.names = 1)
Dose50 <- read.csv(file = "Data50", row.names = 1)
Dose100 <- read.csv(file = "Data100", row.names = 1)
Dose20 <- read.csv(file = "Data20", row.names = 1)
Dose200 <- read.csv(file = "Data200", row.names = 1)

Tid <- Dose02[,1]
DC02 <- Dose02[,2]
DC50 <- Dose50$V2
DC100 <- Dose100$V2
DC20 <- Dose20$V2
DC200 <- Dose200$V2

# Plot af DC ved forskellig dose
plot(Tid,DC100,lty=2, lwd=1,'l',col="red",xlim = c(0,40),ylab = "Mængde Dendridiske celler",xlab = "Tid [Dage]")
lines(Tid,DC50,lty=2, lwd=1,'l',col="blue") # denne her
lines(Tid,DC02,lty=2, lwd=1,'l',col="black")
lines(Tid,DC20,lty=2, lwd=1,'l',col="pink")
lines(Tid,DC200,lty=2, lwd=1,'l',col="green") # denne her
legend(25, 50, legend=c("Dosis 100", "Dosis 50","Dosis 02","Dosis 20","Dosis 200"),
       col=c("red", "blue","black","pink","green"), lty=1:4, cex=0.8)
title("Dendridiske celler")


MC02 <- Dose02[,4]
MC50 <- Dose50$V4
MC100 <- Dose100$V4
MC20 <- Dose20$V4
MC200 <- Dose200$V4

# Plot af makrofager ved forskellig dose
plot(Tid,MC100,lty=2, lwd=1,'l',col="red",xlim = c(0,40),ylim = c(0,150),ylab = "Mængde Makrofager",xlab = "Tid [Dage]")
lines(Tid,MC50,lty=2, lwd=1,'l',col="blue") # denne her
lines(Tid,MC02,lty=2, lwd=1,'l',col="black")
lines(Tid,MC20,lty=2, lwd=1,'l',col="pink")
lines(Tid,MC200,lty=2, lwd=1,'l',col="green") # denne her
legend(25, 150, legend=c("Dosis 100", "Dosis 50","Dosis 02","Dosis 2.0","Dosis 200"),
       col=c("red", "blue","black","pink","green"), lty=1:4, cex=0.8)
title("Makrofager")




Tid <- Dose02[,1]
IFN02 <- Dose02[,3]
IFN50 <- Dose50$V3
IFN100 <- Dose100$V3
IFN20 <- Dose20$V3
IFN200 <- Dose200$V3

# Plot af makrofager ved forskellig dose
plot(Tid,IFN100,lty=2, lwd=1,'l',col="red",xlim = c(0,40),ylim = c(0,0.1),ylab = "Mængde Makrofager",xlab = "Tid [Dage]")
lines(Tid,IFN50,lty=2, lwd=1,'l',col="blue") # denne her
lines(Tid,IFN02,lty=2, lwd=1,'l',col="black")
lines(Tid,IFN20,lty=2, lwd=1,'l',col="pink")
lines(Tid,IFN200,lty=2, lwd=1,'l',col="green") # denne her
legend(25, 150, legend=c("Dosis 100", "Dosis 50","Dosis 0.2","Dosis 20","Dosis 200"),
       col=c("red", "blue","black","pink","green"), lty=1:4, cex=0.8)
title("Makrofager")
