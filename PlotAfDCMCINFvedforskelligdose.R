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


par(mfrow=c(3,2))

graphics.off()
plot(Tid,DC100,lty=2, lwd=1,'l',col="red")
lines(Tid,DC50,lty=2, lwd=1,'l',col="blue")
lines(Tid,DC02,lty=2, lwd=1,'l',col="black")
lines(Tid,DC20,lty=2, lwd=1,'l',col="pink")
lines(Tid,DC200,lty=2, lwd=1,'l',col="green")
legend(10, 60, legend=c("Dose 100", "Dose 50","Dose 02","Dose 20","Dose 200"),
       col=c("red", "blue","black","pink","green"), lty=1:4, cex=0.8)


plot(Tid,DC100,lty=2, lwd=1,'l',col="red")
lines(Tid,DC50,lty=2, lwd=1,'l',col="blue")
lines(Tid,DC02,lty=2, lwd=1,'l',col="black")




plot(Dose02$V1,Dose02$V2,lty=2, lwd=1,'l',col="red")

#legend(title="Tree")

