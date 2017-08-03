library(rgdal)
library(fields)

setwd("")

###  Load data

#empirical species distributions
PresAbs_BR_NH_WH <- read.csv("PresAbs_BR_NH_WH.csv")
PresAbs_NB_NH_WH <- read.csv("PresAbs_NB_NH_WH.csv")
PresAbs_BR_SH_WH <- read.csv("PresAbs_BR_SH_WH.csv")
PresAbs_NB_SH_WH <- read.csv("PresAbs_NB_SH_WH.csv")
PresAbs_res_WH <- read.csv("PresAbs_res_WH.csv")
PresAbs_BR_NH_EH <- read.csv("PresAbs_BR_NH_EH.csv")
PresAbs_NB_NH_EH <- read.csv("PresAbs_NB_NH_EH.csv")
PresAbs_BR_SH_EH <- read.csv("PresAbs_BR_SH_EH.csv")
PresAbs_NB_SH_EH <- read.csv("PresAbs_NB_SH_EH.csv")
PresAbs_res_EH <- read.csv("PresAbs_res_EH.csv")

#model outputs
bestfit_rangesBR_WH <- read.csv("model_outputs/bestfit_rangesBR_WH.csv", header=F)
bestfit_rangesNB_WH <- read.csv("model_outputs/bestfit_rangesNB_WH.csv", header=F)
bestguess_rangesBR_WH <- read.csv("model_outputs/bestguess_rangesBR_WH.csv", header=F)
bestguess_rangesNB_WH <- read.csv("model_outputs/bestguess_rangesNB_WH.csv", header=F)
noenergy_rangesBR_WH <- read.csv("model_outputs/noenergy_rangesBR_WH.csv", header=F)
noenergy_rangesNB_WH <- read.csv("model_outputs/noenergy_rangesNB_WH.csv", header=F)
random_rangesBR_WH <- read.csv("model_outputs/random_rangesBR_WH.csv", header=F)
random_rangesNB_WH <- read.csv("model_outputs/random_rangesNB_WH.csv", header=F)
bestfit_rangesBR_EH <- read.csv("model_outputs/bestfit_rangesBR_EH.csv", header=F)
bestfit_rangesNB_EH <- read.csv("model_outputs/bestfit_rangesNB_EH.csv", header=F)
bestguess_rangesBR_EH <- read.csv("model_outputs/bestguess_rangesBR_EH.csv", header=F)
bestguess_rangesNB_EH <- read.csv("model_outputs/bestguess_rangesNB_EH.csv", header=F)
noenergy_rangesBR_EH <- read.csv("model_outputs/noenergy_rangesBR_EH.csv", header=F)
noenergy_rangesNB_EH <- read.csv("model_outputs/noenergy_rangesNB_EH.csv", header=F)
random_rangesBR_EH <- read.csv("model_outputs/random_rangesBR_EH.csv", header=F)
random_rangesNB_EH <- read.csv("model_outputs/random_rangesNB_EH.csv", header=F)

#environmental data
lonlat.patterns <- read.csv("lonlatMecha.csv", header=F)
envdata <- read.csv("Env_data_6months.csv", header=T)
for(i in 1:dim(lonlat.patterns)[1]){
	lonlat.patterns[i,3] <- envdata$ISEA7_ID[which(envdata$LATITUDE == lonlat.patterns[i,2] & envdata$LONGITUDE == lonlat.patterns[i,1])	]
}
colnames(lonlat.patterns) <- c("Longitude", "Latitude", "HexID")

#global hexagon grid
hexgrid <- readOGR("isea3h7_analyses_clean", verbose=FALSE)
hexgrid <- hexgrid[,c(1,2,15,16)]
hexgridWH <- hexgrid[which(hexgrid@data[,4] <= -30),]
hexgridEH <- hexgrid[which(hexgrid@data[,4] > -30),]
hexgrid2 <- rbind(hexgridWH, hexgridEH)[match(lonlat.patterns$HexID, c(hexgridWH@data[,1], hexgridEH@data[,1])),]


###  Observed Patterns

perm.richness.obs <- c(apply(PresAbs_BR_NH_WH * PresAbs_NB_NH_WH, 1, sum) + apply(PresAbs_BR_SH_WH * PresAbs_NB_SH_WH, 1, sum), apply(PresAbs_BR_NH_EH * PresAbs_NB_NH_EH, 1, sum) + apply(PresAbs_BR_SH_EH * PresAbs_NB_SH_EH, 1, sum))

breeding.richness.obs <- c(apply(PresAbs_BR_NH_WH, 1, sum) + apply(PresAbs_BR_SH_WH, 1, sum), apply(PresAbs_BR_NH_EH, 1, sum) + apply(PresAbs_BR_SH_EH, 1, sum))
breeding.richness.obs = breeding.richness.obs - perm.richness.obs

nonbreeding.richness.obs <- c(apply(PresAbs_NB_NH_WH, 1, sum) + apply(PresAbs_NB_SH_WH, 1, sum), apply(PresAbs_NB_NH_EH, 1, sum) + apply(PresAbs_NB_SH_EH, 1, sum))
nonbreeding.richness.obs = nonbreeding.richness.obs - perm.richness.obs

resident.richness.obs <- c(apply(PresAbs_res_WH, 1, sum), apply(PresAbs_res_EH, 1, sum))
resident.richness.obs = resident.richness.obs + perm.richness.obs

proportion.migrants.obs = (breeding.richness.obs + nonbreeding.richness.obs) / (breeding.richness.obs + nonbreeding.richness.obs + resident.richness.obs)

difference.richness.obs = breeding.richness.obs - nonbreeding.richness.obs




###  Simulated Patterns 

#best-fit
perm.richness.bestfit <- c(apply(bestfit_rangesBR_WH * bestfit_rangesNB_WH, 1, sum), apply(bestfit_rangesBR_EH * bestfit_rangesNB_EH, 1, sum))
breeding.richness.bestfit <- c(apply(bestfit_rangesBR_WH, 1, sum), apply(bestfit_rangesBR_EH, 1, sum))
breeding.richness.bestfit = breeding.richness.bestfit - perm.richness.bestfit
nonbreeding.richness.bestfit <- c(apply(bestfit_rangesNB_WH, 1, sum), apply(bestfit_rangesNB_EH, 1, sum))
nonbreeding.richness.bestfit = nonbreeding.richness.bestfit - perm.richness.bestfit
resident.richness.bestfit <- perm.richness.bestfit
proportion.migrants.bestfit = (breeding.richness.bestfit + nonbreeding.richness.bestfit) / (breeding.richness.bestfit + nonbreeding.richness.bestfit + resident.richness.bestfit)
difference.richness.bestfit = breeding.richness.bestfit - nonbreeding.richness.bestfit

#best-guess
perm.richness.bestguess <- c(apply(bestguess_rangesBR_WH * bestguess_rangesNB_WH, 1, sum), apply(bestguess_rangesBR_EH * bestguess_rangesNB_EH, 1, sum))
breeding.richness.bestguess <- c(apply(bestguess_rangesBR_WH, 1, sum), apply(bestguess_rangesBR_EH, 1, sum))
breeding.richness.bestguess = breeding.richness.bestguess - perm.richness.bestguess
nonbreeding.richness.bestguess <- c(apply(bestguess_rangesNB_WH, 1, sum), apply(bestguess_rangesNB_EH, 1, sum))
nonbreeding.richness.bestguess = nonbreeding.richness.bestguess - perm.richness.bestguess
resident.richness.bestguess <- perm.richness.bestguess
proportion.migrants.bestguess = (breeding.richness.bestguess + nonbreeding.richness.bestguess) / (breeding.richness.bestguess + nonbreeding.richness.bestguess + resident.richness.bestguess)
difference.richness.bestguess = breeding.richness.bestguess - nonbreeding.richness.bestguess

#no energy
perm.richness.noenergy <- c(apply(noenergy_rangesBR_WH * noenergy_rangesNB_WH, 1, sum), apply(noenergy_rangesBR_EH * noenergy_rangesNB_EH, 1, sum))
breeding.richness.noenergy <- c(apply(noenergy_rangesBR_WH, 1, sum), apply(noenergy_rangesBR_EH, 1, sum))
breeding.richness.noenergy = breeding.richness.noenergy - perm.richness.noenergy
nonbreeding.richness.noenergy <- c(apply(noenergy_rangesNB_WH, 1, sum), apply(noenergy_rangesNB_EH, 1, sum))
nonbreeding.richness.noenergy = nonbreeding.richness.noenergy - perm.richness.noenergy
resident.richness.noenergy <- perm.richness.noenergy
proportion.migrants.noenergy = (breeding.richness.noenergy + nonbreeding.richness.noenergy) / (breeding.richness.noenergy + nonbreeding.richness.noenergy + resident.richness.noenergy)
difference.richness.noenergy = breeding.richness.noenergy - nonbreeding.richness.noenergy

#random
perm.richness.random <- c(apply(random_rangesBR_WH * random_rangesNB_WH, 1, sum), apply(random_rangesBR_EH * random_rangesNB_EH, 1, sum))
breeding.richness.random <- c(apply(random_rangesBR_WH, 1, sum), apply(random_rangesBR_EH, 1, sum))
breeding.richness.random = breeding.richness.random - perm.richness.random
nonbreeding.richness.random <- c(apply(random_rangesNB_WH, 1, sum), apply(random_rangesNB_EH, 1, sum))
nonbreeding.richness.random = nonbreeding.richness.random - perm.richness.random
resident.richness.random <- perm.richness.random
proportion.migrants.random = (breeding.richness.random + nonbreeding.richness.random) / (breeding.richness.random + nonbreeding.richness.random + resident.richness.random)
difference.richness.random = breeding.richness.random - nonbreeding.richness.random






## FIGURE 1 ##

hexgridWH2 <- hexgrid2[1:2323,]
energy_richnessT0 <- read.csv("model_outputs/energy_richnessT0.csv")
energy_richnessTi <- read.csv("model_outputs/energy_richnessTi.csv")
energy_richnessTend <- read.csv("model_outputs/energy_richnessTend.csv")

jpeg("Fig1_energyT0.jpg", width=1000, height=1400, quality=100)
energyT0 <- (energy_richnessT0$energyAvailableNS + energy_richnessT0$energyAvailableNW) / 2
rbPal <- colorRampPalette(c("yellow", "dark green"))
datcol <- rbPal(5)[as.numeric(cut(energyT0, breaks=c(-1,100,200,300,400,500)))]
datcol[which(energyT0 ==0)] <- "light grey"
plot(hexgridWH2, col=datcol, border=datcol, bg="white")
dev.off()

jpeg("Fig1_richnessT0.jpg", width=1000, height=1400, quality=100)
richnessT0 <- (energy_richnessT0$richnessBR + energy_richnessT0$richnessNB + energy_richnessT0$richnessRes) 
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(6)[as.numeric(cut(richnessT0, breaks=c(-1,100,200,300,400,500,600)))]
datcol[which(richnessT0 ==0)] <- "light grey"
plot(hexgridWH2, col=datcol, border=datcol, bg="white")
dev.off()




## FIGURE 2 ##

# Left Part 
jpeg("Fig2a.jpg", width=1400, height=2000, quality=100)

par(mfrow=c(5,2), mar=c(4.5,3,2,0.5), mgp=c(3,1,0))

#Breeding richness
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(breeding.richness.obs, breaks=c(-1,25,50,75,100,125,150,175)))]
datcol[which(breeding.richness.obs ==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("A", cex=2.5, side=3, line=-0.5, at=-200)
#mtext("Empirical patterns", cex=1, side=3, line=0.15, at=0)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 150","125–150", "100–125", "75–100", "50–75", "25–50", "< 25", "No species"), fill=c(rev(rbPal(7)),"white"), cex=2.2, border=c(rev(rbPal(7)),"white"))
#mtext("Number of", cex=0.7, side=1, line=-6.4, at=-155)
#mtext("breeding", cex=0.7, side=1, line=-5.8, at=-155)
#mtext("migrants", cex=0.7, side=1, line=-5.2, at=-155)

rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(breeding.richness.bestfit, breaks=c(-1,25,50,75,100,125,150,175)))]
datcol[which(breeding.richness.bestfit==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("B", cex=2.5, side=3, line=-0.5, at=-200)
#mtext("Simulated patterns", cex=1, side=3, line=0.15, at=0)

#Non-breeding richness
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(nonbreeding.richness.obs, breaks=c(-1,25,50,75,100,125,150,175)))]
datcol[which(nonbreeding.richness.obs ==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("E", cex=2.5, side=3, line=-0.5, at=-200)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 150","125–150", "100–125", "75–100", "50–75", "25–50", "< 25", "No species"), fill=c(rev(rbPal(7)),"white"), cex=2.2, border=c(rev(rbPal(7)),"white"))
#mtext("Number of", cex=0.41, side=1, line=-6.4, at=-155)
#mtext("non-breeding", cex=0.41, side=1, line=-5.8, at=-155)
#mtext("migrants", cex=0.41, side=1, line=-5.2, at=-155)

rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(nonbreeding.richness.bestfit, breaks=c(-1,25,50,75,100,125,150,175)))]
datcol[which(nonbreeding.richness.bestfit==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("F", cex=2.5, side=3, line=-0.5, at=-200)

#Resident richness
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(resident.richness.obs, breaks=c(-1,100,200,300,400,500,600,1100)))]
datcol[which(resident.richness.obs ==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("I", cex=2.5, side=3, line=-0.5, at=-200)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 600","500–600", "400–500", "300–400", "200–300", "100–200", "< 100", "No species"), fill=c(rev(rbPal(7)),"white"), cex=2.2, border=c(rev(rbPal(7)),"white"))
#mtext("Number of", cex=0.41, side=1, line=-5.8, at=-155)
#mtext("residents", cex=0.41, side=1, line=-5.2, at=-155)

rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(resident.richness.bestfit, breaks=c(-1,100,200,300,400,500,600,1100)))]
datcol[which(resident.richness.bestfit==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("J", cex=2.5, side=3, line=-0.5, at=-200)

#Difference in richness
rbPal <- colorRampPalette(c("blue", "white", "red3"))
datcol <- rbPal(7)[as.numeric(cut(difference.richness.obs, breaks=c(-160,-110,-60,-10,10,60,110,160)))]
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("M", cex=2.5, side=3, line=-0.5, at=-200)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 110","60–110", "10–60", "-10–10", "-60–-10", "-110–-60", "< -110"), fill=rev(rbPal(7)), cex=2.2, border=rev(rbPal(7)))
#mtext("Difference", cex=0.41, side=1, line=-6, at=-155)
#mtext("in seasonal", cex=0.41, side=1, line=-5.4, at=-155)
#mtext("richness", cex=0.41, side=1, line=-4.8, at=-155)

rbPal <- colorRampPalette(c("blue", "white", "red3"))
datcol <- rbPal(7)[as.numeric(cut(difference.richness.bestfit, breaks=c(-160,-110,-60,-10,10,60,110,160)))]
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("N", cex=2.5, side=3, line=-0.5, at=-200)

#Proportion of migrants
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(proportion.migrants.obs, breaks=c(-0.1,0.15,0.30,0.45,0.60,0.75,0.90,1.1)))]
datcol[which(proportion.migrants.obs==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("Q", cex=2.5, side=3, line=-0.5, at=-200)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 0.9","0.75–0.9", "0.6–0.75", "0.45–0.6", "0.3–0.45", "0.15–0.3", "< 0.15", "No migrants"), fill=c(rev(rbPal(7)),"white"), cex=2.2, border=c(rev(rbPal(7)),"white"))
#mtext("Proportion", cex=0.41, side=1, line=-5.8, at=-155)
#mtext("of migrants", cex=0.41, side=1, line=-5.2, at=-155)

rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(proportion.migrants.bestfit, breaks=c(-0.1,0.15,0.30,0.45,0.60,0.75,0.90,1.1)))]
datcol[which(proportion.migrants.bestfit==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("R", cex=2.5, side=3, line=-0.5, at=-200)

dev.off()



# Right Part 

jpeg("Fig2b.jpg", width=800, height=2000, quality=100)

par(mfrow=c(5,2), mar=c(4.5,5,2,1), mgp=c(3,1,0))

#Breeding richness
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(0,120), ylab="Latitude", xaxt="n", axes=F, xlab="Species richness", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(lonlat.patterns[,2], breeding.richness.obs, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], breeding.richness.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(lonlat.patterns[,2], breeding.richness.bestfit, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], breeding.richness.bestfit, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("C", cex=2.5, side=3, line=-0.5, at=-15)
legend("bottomright", inset=0.1, bg="white", box.col="white", c("Empirical trend", "Simulated trend"), col=c("black","red"), cex=2.2, lwd=2)

plot(breeding.richness.obs, breeding.richness.bestfit, pch=20, xlim=c(0,210), ylim=c(0,210), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical richness", ylab="Simulated richness", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
#mtext("Empirical vs simulated", cex=1, side=3, line=0.15, at=100)
mtext("r=0.867", cex=2, side=1, line=-2.5, at=180)
mtext("D", cex=2.5, side=3, line=-0.5, at=-50)

#Non-breeding richness
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(0,120), ylab="Latitude", xaxt="n", axes=F, xlab="Species richness", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(lonlat.patterns[,2], nonbreeding.richness.obs, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], nonbreeding.richness.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(lonlat.patterns[,2], nonbreeding.richness.bestfit, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], nonbreeding.richness.bestfit, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("G", cex=2.5, side=3, line=-0.5, at=-15)

plot(nonbreeding.richness.obs, nonbreeding.richness.bestfit, pch=20, xlim=c(0,210), ylim=c(0,210), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical richness", ylab="Simulated richness", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
#mtext("Empirical vs simulated", cex=1, side=3, line=0.15, at=100)
mtext("r=0.748", cex=2, side=1, line=-2.5, at=180)
mtext("H", cex=2.5, side=3, line=-0.5, at=-50)

#Residents
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(0,400), ylab="Latitude", xaxt="n", axes=F, xlab="Species richness", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(lonlat.patterns[,2], resident.richness.obs, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], resident.richness.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(lonlat.patterns[,2], resident.richness.bestfit, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], resident.richness.bestfit, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("K", cex=2.5, side=3, line=-0.5, at=-49)

plot(resident.richness.obs, resident.richness.bestfit, pch=20, xlim=c(0,1000), ylim=c(0,800), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical richness", ylab="Simulated richness", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
#mtext("Empirical vs simulated", cex=1, side=3, line=0.15, at=100)
mtext("r=0.696", cex=2, side=3, line=-6, at=180)
mtext("L", cex=2.5, side=3, line=-0.5, at=-230)

#Difference in richness
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(-100,100), ylab="Latitude", xaxt="n", axes=F, xlab="Difference in richness", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(lonlat.patterns[,2], difference.richness.obs, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], difference.richness.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(lonlat.patterns[,2], difference.richness.bestfit, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], difference.richness.bestfit, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("O", cex=2.5, side=3, line=-0.5, at=-125)

plot(difference.richness.obs, difference.richness.bestfit, pch=20, xlim=c(-190,190), ylim=c(-190,190), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical difference", ylab="Simulated difference", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
#mtext("Empirical vs simulated", cex=1, side=3, line=0.15, at=100)
mtext("r=0.871", cex=2, side=3, line=-6, at=-122)
mtext("P", cex=2.5, side=3, line=-0.5, at=-280)

#Proportion of migrants
tokeep = which(is.na(proportion.migrants.obs)==F & is.na(proportion.migrants.bestfit)==F)
prop.obs = proportion.migrants.obs[tokeep]
prop.bestfit = proportion.migrants.bestfit[tokeep]
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(0,1), ylab="Latitude", xaxt="n", axes=F, xlab="Proportion of migrants", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(lonlat.patterns[tokeep,2], prop.obs, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[tokeep,2], prop.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(lonlat.patterns[tokeep,2], prop.bestfit, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[tokeep,2], prop.bestfit, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("S", cex=2.5, side=3, line=-0.5, at=-0.13)

plot(prop.obs, prop.bestfit, pch=20, xlim=c(0,1), ylim=c(0,1), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical proportion", ylab="Simulated proportion", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
#mtext("Empirical vs simulated", cex=1, side=3, line=0.15, at=100)
mtext("r=0.837", cex=2, side=3, line=-6, at=0.13)
mtext("T", cex=2.5, side=3, line=-0.5, at=-0.23)

dev.off()





## FIGURE 3 ##

EMDrandom = c(15070.04877862, 14874.24106341, 14848.6339411, 15128.38885583, 14718.56616819, 14811.93059197, 14286.96650175, 14936.46057012, 14319.97380384, 14619.64875535, 14938.40594741, 15255.37661177, 15122.1796363, 14799.53505901, 14830.38351451, 15042.20544426, 14607.97031363, 15131.36282917, 15161.70910288, 14588.57941737, 15152.97933711,  14628.04664144,  15142.52945065,  15077.37874877, 14887.5091581,14752.09731643,  15096.52487758,  15421.71639309,  14957.63111244, 14812.64824094,15150.262564  ,  14875.1368278 ,  15358.53082714,  14933.89935371, 14466.93908091,  15235.77118226,  15251.83772329,  14892.98154497, 15152.16477366,  14521.77461875,  14917.38667495,  15030.37183644, 14659.95859803,  14854.59014217,  14799.06598113,  14883.25767138, 14902.41998327,  14802.84578178,  15170.21958684,  14752.77211437, 14951.04198296,  14781.80296817,  14969.79792414,  14672.14018681, 15743.14545695,  14862.4939531 ,  15168.63940031,  14773.46395046, 15188.51281814,  14553.97381536,  14212.67829694,  14407.73396212, 14983.63517247,  14793.77282927,  15111.20908505,  14685.18330175, 14931.51372845,  14590.26302761,  15176.13129261,  14855.50145836, 14441.07842778,  14568.75461929,  14901.47193349,  14800.68740383, 14884.89572168,  15048.38188962,  14586.7864773 ,  15213.56061367, 14841.75650165,  15054.55683057,  15161.08539294,  14909.14784689, 14732.01415058,  15590.13914395,  14589.13799519,  14791.08521445, 15657.41720089,  14759.54391251,  15030.44695098,  14640.21472899, 14828.08698223,  15205.59227586,  15018.20965025,  14812.10644452, 14951.05739924,  15268.19428716,  14856.60248414,  15057.60413885, 14003.22709452,  15016.18668042)
par(mfrow=c(1,1), mar=c(3,0.5,0.1,0.1), mgp=c(1.5,0.5,0))
hist(EMDrandom, xlim=c(0,17000), ylim=c(0,45), main="", xlab="Earth Mover's Distance between empirical and simulated patterns       ", ylab="", axes=F, col="grey", border="dark grey")
#axis(side=2)
axis(side=1)
text(x= mean(EMDrandom), y=37, "Null\ndistribution", col="dark grey", cex=1.1)
#hist(energy.efficiency.distrib[,5], xlim=c(0,17000), ylim=c(0,30), main="", xlab="Earth Mover's Distance between empirical and simulated patterns", ylab="", axes=F, col="grey", border="dark grey", add=T)
#text(x= mean(EMDrandom), y=26, "Energy\nefficiency\ndistribution", col="dark grey", cex=1.15)
EMDbestfit = 410.67679711
lines(c(EMDbestfit, EMDbestfit), c(0,20), col="black", lty=1, cex=2) # orange2
text(x= EMDbestfit+60, y=23.5, "Best fit", col="black", pos=3, srt=90, cex=1.1) # orange2
EMDbestguess = 770.3910842757856
lines(c(EMDbestguess, EMDbestguess), c(0,20), col="black", lty=1, cex=2) # blue
text(x= EMDbestguess+210, y=25.5, "Best guess", col="black", pos=3, srt=90, cex=1.1) # blue
EMDnoenergy = 12345.145599306152
lines(c(EMDnoenergy, EMDnoenergy), c(0,20), col="dark grey", lty=1, cex=2) # green3
text(x= EMDnoenergy+150, y=29.2, "No energetic costs", col="dark grey", pos=3, srt=90, cex=1.1) # green3







## FIGURE S4 – Results for an example of a randomized model ##

# Left Part 
jpeg("FigS4a.jpg", width=1400, height=2000, quality=100)

par(mfrow=c(5,2), mar=c(4.5,3,2,0.5), mgp=c(3,1,0))

#Breeding richness
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(breeding.richness.obs, breaks=c(-1,25,50,75,100,125,150,175)))]
datcol[which(breeding.richness.obs ==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("A", cex=2.5, side=3, line=-0.5, at=-200)
#mtext("Empirical patterns", cex=1, side=3, line=0.15, at=0)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 150","125–150", "100–125", "75–100", "50–75", "25–50", "< 25", "No species"), fill=c(rev(rbPal(7)),"white"), cex=2.2, border=c(rev(rbPal(7)),"white"))
#mtext("Number of", cex=0.7, side=1, line=-6.4, at=-155)
#mtext("breeding", cex=0.7, side=1, line=-5.8, at=-155)
#mtext("migrants", cex=0.7, side=1, line=-5.2, at=-155)

rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(breeding.richness.random, breaks=c(-1,25,50,75,100,125,150,700)))]
datcol[which(breeding.richness.random==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("B", cex=2.5, side=3, line=-0.5, at=-200)
#mtext("Simulated patterns", cex=1, side=3, line=0.15, at=0)

#Non-breeding richness
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(nonbreeding.richness.obs, breaks=c(-1,25,50,75,100,125,150,175)))]
datcol[which(nonbreeding.richness.obs ==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("E", cex=2.5, side=3, line=-0.5, at=-200)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 150","125–150", "100–125", "75–100", "50–75", "25–50", "< 25", "No species"), fill=c(rev(rbPal(7)),"white"), cex=2.2, border=c(rev(rbPal(7)),"white"))
#mtext("Number of", cex=0.41, side=1, line=-6.4, at=-155)
#mtext("non-breeding", cex=0.41, side=1, line=-5.8, at=-155)
#mtext("migrants", cex=0.41, side=1, line=-5.2, at=-155)

rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(nonbreeding.richness.random, breaks=c(-1,25,50,75,100,125,150,700)))]
datcol[which(nonbreeding.richness.random==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("F", cex=2.5, side=3, line=-0.5, at=-200)

#Resident richness
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(resident.richness.obs, breaks=c(-1,100,200,300,400,500,600,1100)))]
datcol[which(resident.richness.obs ==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("I", cex=2.5, side=3, line=-0.5, at=-200)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 600","500–600", "400–500", "300–400", "200–300", "100–200", "< 100", "No species"), fill=c(rev(rbPal(7)),"white"), cex=2.2, border=c(rev(rbPal(7)),"white"))
#mtext("Number of", cex=0.41, side=1, line=-5.8, at=-155)
#mtext("residents", cex=0.41, side=1, line=-5.2, at=-155)

rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(resident.richness.random, breaks=c(-1,100,200,300,400,500,600,1100)))]
datcol[which(resident.richness.random==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("J", cex=2.5, side=3, line=-0.5, at=-200)

#Difference in richness
rbPal <- colorRampPalette(c("blue", "white", "red3"))
datcol <- rbPal(7)[as.numeric(cut(difference.richness.obs, breaks=c(-160,-110,-60,-10,10,60,110,160)))]
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("M", cex=2.5, side=3, line=-0.5, at=-200)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 110","60–110", "10–60", "-10–10", "-60–-10", "-110–-60", "< -110"), fill=rev(rbPal(7)), cex=2.2, border=rev(rbPal(7)))
#mtext("Difference", cex=0.41, side=1, line=-6, at=-155)
#mtext("in seasonal", cex=0.41, side=1, line=-5.4, at=-155)
#mtext("richness", cex=0.41, side=1, line=-4.8, at=-155)

rbPal <- colorRampPalette(c("blue", "white", "red3"))
datcol <- rbPal(7)[as.numeric(cut(difference.richness.random, breaks=c(-400,-110,-60,-10,10,60,110,400)))]
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("N", cex=2.5, side=3, line=-0.5, at=-200)

#Proportion of migrants
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(proportion.migrants.obs, breaks=c(-0.1,0.15,0.30,0.45,0.60,0.75,0.90,1.1)))]
datcol[which(proportion.migrants.obs==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("Q", cex=2.5, side=3, line=-0.5, at=-200)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 0.9","0.75–0.9", "0.6–0.75", "0.45–0.6", "0.3–0.45", "0.15–0.3", "< 0.15", "No migrants"), fill=c(rev(rbPal(7)),"white"), cex=2.2, border=c(rev(rbPal(7)),"white"))
#mtext("Proportion", cex=0.41, side=1, line=-5.8, at=-155)
#mtext("of migrants", cex=0.41, side=1, line=-5.2, at=-155)

rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(proportion.migrants.random, breaks=c(-0.1,0.15,0.30,0.45,0.60,0.75,0.90,1.1)))]
datcol[which(proportion.migrants.random==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("R", cex=2.5, side=3, line=-0.5, at=-200)

dev.off()



# Right Part 

jpeg("FigS4b.jpg", width=800, height=2000, quality=100)

par(mfrow=c(5,2), mar=c(4.5,5,2,1), mgp=c(3,1,0))

#Breeding richness
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(0,300), ylab="Latitude", xaxt="n", axes=F, xlab="Species richness", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(lonlat.patterns[,2], breeding.richness.obs, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], breeding.richness.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(lonlat.patterns[,2], breeding.richness.random, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], breeding.richness.random, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("C", cex=2.5, side=3, line=-0.5, at=-39)
legend("bottomright", inset=0.05, bg="white", box.col="white", c("Empirical trend", "Simulated trend"), col=c("black","red"), cex=2.2, lwd=2)

plot(breeding.richness.obs, breeding.richness.random, pch=20, xlim=c(0,500), ylim=c(0,500), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical richness", ylab="Simulated richness", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
#mtext("Empirical vs simulated", cex=1, side=3, line=0.15, at=100)
mtext("r=0.193", cex=2, side=1, line=-2.5, at=430)
mtext("D", cex=2.5, side=3, line=-0.5, at=-115)

#Non-breeding richness
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(0,300), ylab="Latitude", xaxt="n", axes=F, xlab="Species richness", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(lonlat.patterns[,2], nonbreeding.richness.obs, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], nonbreeding.richness.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(lonlat.patterns[,2], nonbreeding.richness.random, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], nonbreeding.richness.random, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("G", cex=2.5, side=3, line=-0.5, at=-39)

plot(nonbreeding.richness.obs, nonbreeding.richness.random, pch=20, xlim=c(0,500), ylim=c(0,500), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical richness", ylab="Simulated richness", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
#mtext("Empirical vs simulated", cex=1, side=3, line=0.15, at=100)
mtext("r=0.086", cex=2, side=1, line=-2.5, at=430)
mtext("H", cex=2.5, side=3, line=-0.5, at=-115)

#Residents
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(0,400), ylab="Latitude", xaxt="n", axes=F, xlab="Species richness", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(lonlat.patterns[,2], resident.richness.obs, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], resident.richness.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(lonlat.patterns[,2], resident.richness.random, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], resident.richness.random, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("K", cex=2.5, side=3, line=-0.5, at=-51)

plot(resident.richness.obs, resident.richness.random, pch=20, xlim=c(0,1000), ylim=c(0,800), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical richness", ylab="Simulated richness", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
#mtext("Empirical vs simulated", cex=1, side=3, line=0.15, at=100)
mtext("r=0.203", cex=2, side=3, line=-6, at=180)
mtext("L", cex=2.5, side=3, line=-0.5, at=-230)

#Difference in richness
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(-100,100), ylab="Latitude", xaxt="n", axes=F, xlab="Difference in richness", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(lonlat.patterns[,2], difference.richness.obs, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], difference.richness.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(lonlat.patterns[,2], difference.richness.random, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], difference.richness.random, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("O", cex=2.5, side=3, line=-0.5, at=-125)

plot(difference.richness.obs, difference.richness.random, pch=20, xlim=c(-200,200), ylim=c(-200,200), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical difference", ylab="Simulated difference", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
#mtext("Empirical vs simulated", cex=1, side=3, line=0.15, at=100)
mtext("r=-0.031", cex=2, side=3, line=-6, at=-140)
mtext("P", cex=2.5, side=3, line=-0.5, at=-290)

#Proportion of migrants
tokeep = which(is.na(proportion.migrants.obs)==F & is.na(proportion.migrants.random)==F)
prop.obs = proportion.migrants.obs[tokeep]
prop.random = proportion.migrants.random[tokeep]
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(0,1), ylab="Latitude", xaxt="n", axes=F, xlab="Proportion of migrants", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(lonlat.patterns[tokeep,2], prop.obs, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[tokeep,2], prop.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(lonlat.patterns[tokeep,2], prop.random, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[tokeep,2], prop.random, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("S", cex=2.5, side=3, line=-0.5, at=-0.13)

plot(prop.obs, prop.random, pch=20, xlim=c(0,1), ylim=c(0,1), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical proportion", ylab="Simulated proportion", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
#mtext("Empirical vs simulated", cex=1, side=3, line=0.15, at=100)
mtext("r=-0.057", cex=2, side=1, line=-2.5, at=0.85)
mtext("T", cex=2.5, side=3, line=-0.5, at=-0.23)

dev.off()







## FIGURE S5 – no energetic costs results ##

# Left Part 
jpeg("FigS5a.jpg", width=1400, height=2000, quality=100)

par(mfrow=c(5,2), mar=c(4.5,3,2,0.5), mgp=c(3,1,0))

#Breeding richness
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(breeding.richness.obs, breaks=c(-1,25,50,75,100,125,150,175)))]
datcol[which(breeding.richness.obs ==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("A", cex=2.5, side=3, line=-0.5, at=-200)
#mtext("Empirical patterns", cex=1, side=3, line=0.15, at=0)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 150","125–150", "100–125", "75–100", "50–75", "25–50", "< 25", "No species"), fill=c(rev(rbPal(7)),"white"), cex=2.2, border=c(rev(rbPal(7)),"white"))
#mtext("Number of", cex=0.7, side=1, line=-6.4, at=-155)
#mtext("breeding", cex=0.7, side=1, line=-5.8, at=-155)
#mtext("migrants", cex=0.7, side=1, line=-5.2, at=-155)

rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(breeding.richness.noenergy, breaks=c(-1,25,50,75,100,125,150,700)))]
datcol[which(breeding.richness.noenergy==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("B", cex=2.5, side=3, line=-0.5, at=-200)
#mtext("Simulated patterns", cex=1, side=3, line=0.15, at=0)

#Non-breeding richness
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(nonbreeding.richness.obs, breaks=c(-1,25,50,75,100,125,150,175)))]
datcol[which(nonbreeding.richness.obs ==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("E", cex=2.5, side=3, line=-0.5, at=-200)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 150","125–150", "100–125", "75–100", "50–75", "25–50", "< 25", "No species"), fill=c(rev(rbPal(7)),"white"), cex=2.2, border=c(rev(rbPal(7)),"white"))
#mtext("Number of", cex=0.41, side=1, line=-6.4, at=-155)
#mtext("non-breeding", cex=0.41, side=1, line=-5.8, at=-155)
#mtext("migrants", cex=0.41, side=1, line=-5.2, at=-155)

rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(nonbreeding.richness.noenergy, breaks=c(-1,25,50,75,100,125,150,700)))]
datcol[which(nonbreeding.richness.noenergy==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("F", cex=2.5, side=3, line=-0.5, at=-200)

#Resident richness
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(resident.richness.obs, breaks=c(-1,100,200,300,400,500,600,1100)))]
datcol[which(resident.richness.obs ==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("I", cex=2.5, side=3, line=-0.5, at=-200)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 600","500–600", "400–500", "300–400", "200–300", "100–200", "< 100", "No species"), fill=c(rev(rbPal(7)),"white"), cex=2.2, border=c(rev(rbPal(7)),"white"))
#mtext("Number of", cex=0.41, side=1, line=-5.8, at=-155)
#mtext("residents", cex=0.41, side=1, line=-5.2, at=-155)

rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(resident.richness.noenergy, breaks=c(-1,100,200,300,400,500,600,1100)))]
datcol[which(resident.richness.noenergy==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("J", cex=2.5, side=3, line=-0.5, at=-200)

#Difference in richness
rbPal <- colorRampPalette(c("blue", "white", "red3"))
datcol <- rbPal(7)[as.numeric(cut(difference.richness.obs, breaks=c(-160,-110,-60,-10,10,60,110,160)))]
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("M", cex=2.5, side=3, line=-0.5, at=-200)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 110","60–110", "10–60", "-10–10", "-60–-10", "-110–-60", "< -110"), fill=rev(rbPal(7)), cex=2.2, border=rev(rbPal(7)))
#mtext("Difference", cex=0.41, side=1, line=-6, at=-155)
#mtext("in seasonal", cex=0.41, side=1, line=-5.4, at=-155)
#mtext("richness", cex=0.41, side=1, line=-4.8, at=-155)

rbPal <- colorRampPalette(c("blue", "white", "red3"))
datcol <- rbPal(7)[as.numeric(cut(difference.richness.noenergy, breaks=c(-400,-110,-60,-10,10,60,110,400)))]
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("N", cex=2.5, side=3, line=-0.5, at=-200)

#Proportion of migrants
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(proportion.migrants.obs, breaks=c(-0.1,0.15,0.30,0.45,0.60,0.75,0.90,1.1)))]
datcol[which(proportion.migrants.obs==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("Q", cex=2.5, side=3, line=-0.5, at=-200)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 0.9","0.75–0.9", "0.6–0.75", "0.45–0.6", "0.3–0.45", "0.15–0.3", "< 0.15", "No migrants"), fill=c(rev(rbPal(7)),"white"), cex=2.2, border=c(rev(rbPal(7)),"white"))
#mtext("Proportion", cex=0.41, side=1, line=-5.8, at=-155)
#mtext("of migrants", cex=0.41, side=1, line=-5.2, at=-155)

rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(proportion.migrants.noenergy, breaks=c(-0.1,0.15,0.30,0.45,0.60,0.75,0.90,1.1)))]
datcol[which(proportion.migrants.noenergy==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("R", cex=2.5, side=3, line=-0.5, at=-200)

dev.off()



# Right Part 

jpeg("FigS5b.jpg", width=800, height=2000, quality=100)

par(mfrow=c(5,2), mar=c(4.5,5,2,1), mgp=c(3,1,0))

#Breeding richness
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(0,400), ylab="Latitude", xaxt="n", axes=F, xlab="Species richness", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(lonlat.patterns[,2], breeding.richness.obs, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], breeding.richness.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(lonlat.patterns[,2], breeding.richness.noenergy, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], breeding.richness.noenergy, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("C", cex=2.5, side=3, line=-0.5, at=-52)
legend("topright", inset=0.05, bg="white", box.col="white", c("Empirical trend", "Simulated trend"), col=c("black","red"), cex=2.2, lwd=2)

plot(breeding.richness.obs, breeding.richness.noenergy, pch=20, xlim=c(0,700), ylim=c(0,700), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical richness", ylab="Simulated richness", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
#mtext("Empirical vs simulated", cex=1, side=3, line=0.15, at=100)
mtext("r=-0.381", cex=2, side=1, line=-2.5, at=600)
mtext("D", cex=2.5, side=3, line=-0.5, at=-170)

#Non-breeding richness
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(0,400), ylab="Latitude", xaxt="n", axes=F, xlab="Species richness", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(lonlat.patterns[,2], nonbreeding.richness.obs, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], nonbreeding.richness.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(lonlat.patterns[,2], nonbreeding.richness.noenergy, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], nonbreeding.richness.noenergy, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("G", cex=2.5, side=3, line=-0.5, at=-52)

plot(nonbreeding.richness.obs, nonbreeding.richness.noenergy, pch=20, xlim=c(0,700), ylim=c(0,700), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical richness", ylab="Simulated richness", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
#mtext("Empirical vs simulated", cex=1, side=3, line=0.15, at=100)
mtext("r=0.239", cex=2, side=1, line=-2.5, at=600)
mtext("H", cex=2.5, side=3, line=-0.5, at=-170)

#Residents
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(0,400), ylab="Latitude", xaxt="n", axes=F, xlab="Species richness", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(lonlat.patterns[,2], resident.richness.obs, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], resident.richness.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(lonlat.patterns[,2], resident.richness.noenergy, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], resident.richness.noenergy, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("K", cex=2.5, side=3, line=-0.5, at=-51)

plot(resident.richness.obs, resident.richness.noenergy, pch=20, xlim=c(0,1000), ylim=c(0,800), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical richness", ylab="Simulated richness", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
#mtext("Empirical vs simulated", cex=1, side=3, line=0.15, at=100)
mtext("r=0.602", cex=2, side=3, line=-6, at=180)
mtext("L", cex=2.5, side=3, line=-0.5, at=-230)

#Difference in richness
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(-200,200), ylab="Latitude", xaxt="n", axes=F, xlab="Difference in richness", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(lonlat.patterns[,2], difference.richness.obs, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], difference.richness.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(lonlat.patterns[,2], difference.richness.noenergy, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], difference.richness.noenergy, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("O", cex=2.5, side=3, line=-0.5, at=-250)

plot(difference.richness.obs, difference.richness.noenergy, pch=20, xlim=c(-400,400), ylim=c(-400,400), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical difference", ylab="Simulated difference", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
#mtext("Empirical vs simulated", cex=1, side=3, line=0.15, at=100)
mtext("r=-0.79", cex=2, side=3, line=-6, at=-275)
mtext("P", cex=2.5, side=3, line=-0.5, at=-575)

#Proportion of migrants
tokeep = which(is.na(proportion.migrants.obs)==F & is.na(proportion.migrants.noenergy)==F)
prop.obs = proportion.migrants.obs[tokeep]
prop.noenergy = proportion.migrants.noenergy[tokeep]
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(0,1), ylab="Latitude", xaxt="n", axes=F, xlab="Proportion of migrants", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(lonlat.patterns[tokeep,2], prop.obs, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[tokeep,2], prop.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(lonlat.patterns[tokeep,2], prop.noenergy, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[tokeep,2], prop.noenergy, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("S", cex=2.5, side=3, line=-0.5, at=-0.13)

plot(prop.obs, prop.noenergy, pch=20, xlim=c(0,1), ylim=c(0,1), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical proportion", ylab="Simulated proportion", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
#mtext("Empirical vs simulated", cex=1, side=3, line=0.15, at=100)
mtext("r=0.552", cex=2, side=1, line=-2.5, at=0.85)
mtext("T", cex=2.5, side=3, line=-0.5, at=-0.23)

dev.off()








## FIGURE S6 –  Best fit vs Best guess parameter space ##



bestguess <- c(0.0000645, 13.4, 0.195, 81.36)
bestfit <- c(6.55E-05, 24.93301313, 0.413860822, 83.66956208)

par(mfrow=c(3,2), mar=c(3.5,3.5,1,1), mgp=c(2,1,0))

plot(bestguess[1], bestguess[2], xlim=c(0,0.0002), ylim=c(0,80), xaxt="n", axes=F, xlab="α (migration cost)", ylab="β (thermoregulation cost)", pch=3, col="blue", cex=1.5)
axis(side=1, labels=T)
axis(side=2, labels=T)
points(bestfit[1], bestfit[2], pch=3, col="orange2", cex=1.5)
points(0, 80, pch=3, col="green3", cex=1.5)
#points(energy.efficiency.distrib[,1], energy.efficiency.distrib[,3], pch=20, col="green3", cex=1)

plot(bestguess[3], bestguess[2], xlim=c(0,2.25), ylim=c(0,80), xaxt="n", axes=F, xlab="γ (reproduction cost)", ylab="β (thermoregulation cost)", pch=3, col="blue", cex=1.5)
axis(side=1, labels=T)
axis(side=2, labels=T)
points(bestfit[3], bestfit[2], pch=3, col="orange2", cex=1.5)
points(0, 80, pch=3, col="green3", cex=1.5)
#points(energy.efficiency.distrib[,2], energy.efficiency.distrib[,3], pch=20, col="green3", cex=1)

plot(bestguess[1], bestguess[4], xlim=c(0,0.0002), ylim=c(50,110), xaxt="n", axes=F, xlab="α (migration cost)", ylab="μ (scaling NDVI)", pch=3, col="blue", cex=1.5)
axis(side=1, labels=T)
axis(side=2, labels=T)
points(bestfit[1], bestfit[4], pch=3, col="orange2", cex=1.5)
points(0, bestguess[4], pch=3, col="green3", cex=1.5)
#points(energy.efficiency.distrib[,1], energy.efficiency.distrib[,4], pch=20, col="green3", cex=1)

plot(bestguess[3], bestguess[4], xlim=c(0,2.25), ylim=c(50,110), xaxt="n", axes=F, xlab="γ (reproduction cost)", ylab="μ (scaling NDVI)", pch=3, col="blue", cex=1.5)
axis(side=1, labels=T)
axis(side=2, labels=T)
points(bestfit[3], bestfit[4], pch=3, col="orange2", cex=1.5)
points(0, bestguess[4], pch=3, col="green3", cex=1.5)
#points(energy.efficiency.distrib[,2], energy.efficiency.distrib[,4], pch=20, col="green3", cex=1)

plot(bestguess[1], bestguess[3], xlim=c(0,0.0002), ylim=c(0,2.25), xaxt="n", axes=F, xlab="α (migration cost)", ylab="γ (reproduction cost)", pch=3, col="blue", cex=1.5)
axis(side=1, labels=T)
axis(side=2, labels=T)
points(bestfit[1], bestfit[3], pch=3, col="orange2", cex=1.5)
points(0, 0, pch=3, col="green3", cex=1.5)
#points(energy.efficiency.distrib[,1], energy.efficiency.distrib[,2], pch=20, col="green3", cex=1)

plot(bestguess[2], bestguess[4], xlim=c(0,80), ylim=c(50,110), xaxt="n", axes=F, xlab="β (thermoregulation cost)", ylab="μ (scaling NDVI)", pch=3, col="blue", cex=1.5)
axis(side=1, labels=T)
axis(side=2, labels=T)
points(bestfit[2], bestfit[4], pch=3, col="orange2", cex=1.5)
points(80, bestguess[4], pch=3, col="green3", cex=1.5)
#points(energy.efficiency.distrib[,3], energy.efficiency.distrib[,4], pch=20, col="green3", cex=1)







## FIGURE S7 – Results best-guess model ##

# Left Part 
jpeg("FigS7a.jpg", width=1400, height=2000, quality=100)

par(mfrow=c(5,2), mar=c(4.5,3,2,0.5), mgp=c(3,1,0))

#Breeding richness
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(breeding.richness.obs, breaks=c(-1,25,50,75,100,125,150,175)))]
datcol[which(breeding.richness.obs ==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("A", cex=2.5, side=3, line=-0.5, at=-200)
#mtext("Empirical patterns", cex=1, side=3, line=0.15, at=0)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 150","125–150", "100–125", "75–100", "50–75", "25–50", "< 25", "No species"), fill=c(rev(rbPal(7)),"white"), cex=2.2, border=c(rev(rbPal(7)),"white"))
#mtext("Number of", cex=0.7, side=1, line=-6.4, at=-155)
#mtext("breeding", cex=0.7, side=1, line=-5.8, at=-155)
#mtext("migrants", cex=0.7, side=1, line=-5.2, at=-155)

rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(breeding.richness.bestguess, breaks=c(-1,25,50,75,100,125,150,175)))]
datcol[which(breeding.richness.bestguess==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("B", cex=2.5, side=3, line=-0.5, at=-200)
#mtext("Simulated patterns", cex=1, side=3, line=0.15, at=0)

#Non-breeding richness
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(nonbreeding.richness.obs, breaks=c(-1,25,50,75,100,125,150,175)))]
datcol[which(nonbreeding.richness.obs ==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("E", cex=2.5, side=3, line=-0.5, at=-200)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 150","125–150", "100–125", "75–100", "50–75", "25–50", "< 25", "No species"), fill=c(rev(rbPal(7)),"white"), cex=2.2, border=c(rev(rbPal(7)),"white"))
#mtext("Number of", cex=0.41, side=1, line=-6.4, at=-155)
#mtext("non-breeding", cex=0.41, side=1, line=-5.8, at=-155)
#mtext("migrants", cex=0.41, side=1, line=-5.2, at=-155)

rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(nonbreeding.richness.bestguess, breaks=c(-1,25,50,75,100,125,150,175)))]
datcol[which(nonbreeding.richness.bestguess==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("F", cex=2.5, side=3, line=-0.5, at=-200)

#Resident richness
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(resident.richness.obs, breaks=c(-1,100,200,300,400,500,600,1100)))]
datcol[which(resident.richness.obs ==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("I", cex=2.5, side=3, line=-0.5, at=-200)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 600","500–600", "400–500", "300–400", "200–300", "100–200", "< 100", "No species"), fill=c(rev(rbPal(7)),"white"), cex=2.2, border=c(rev(rbPal(7)),"white"))
#mtext("Number of", cex=0.41, side=1, line=-5.8, at=-155)
#mtext("residents", cex=0.41, side=1, line=-5.2, at=-155)

rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(resident.richness.bestguess, breaks=c(-1,100,200,300,400,500,600,1100)))]
datcol[which(resident.richness.bestguess==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("J", cex=2.5, side=3, line=-0.5, at=-200)

#Difference in richness
rbPal <- colorRampPalette(c("blue", "white", "red3"))
datcol <- rbPal(7)[as.numeric(cut(difference.richness.obs, breaks=c(-160,-110,-60,-10,10,60,110,160)))]
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("M", cex=2.5, side=3, line=-0.5, at=-200)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 110","60–110", "10–60", "-10–10", "-60–-10", "-110–-60", "< -110"), fill=rev(rbPal(7)), cex=2.2, border=rev(rbPal(7)))
#mtext("Difference", cex=0.41, side=1, line=-6, at=-155)
#mtext("in seasonal", cex=0.41, side=1, line=-5.4, at=-155)
#mtext("richness", cex=0.41, side=1, line=-4.8, at=-155)

rbPal <- colorRampPalette(c("blue", "white", "red3"))
datcol <- rbPal(7)[as.numeric(cut(difference.richness.bestguess, breaks=c(-160,-110,-60,-10,10,60,110,160)))]
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("N", cex=2.5, side=3, line=-0.5, at=-200)

#Proportion of migrants
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(proportion.migrants.obs, breaks=c(-0.1,0.15,0.30,0.45,0.60,0.75,0.90,1.1)))]
datcol[which(proportion.migrants.obs==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("Q", cex=2.5, side=3, line=-0.5, at=-200)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 0.9","0.75–0.9", "0.6–0.75", "0.45–0.6", "0.3–0.45", "0.15–0.3", "< 0.15", "No migrants"), fill=c(rev(rbPal(7)),"white"), cex=2.2, border=c(rev(rbPal(7)),"white"))
#mtext("Proportion", cex=0.41, side=1, line=-5.8, at=-155)
#mtext("of migrants", cex=0.41, side=1, line=-5.2, at=-155)

rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(proportion.migrants.bestguess, breaks=c(-0.1,0.15,0.30,0.45,0.60,0.75,0.90,1.1)))]
datcol[which(proportion.migrants.bestguess==0)] <- "white"
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("R", cex=2.5, side=3, line=-0.5, at=-200)

dev.off()



# Right Part 

jpeg("FigS7b.jpg", width=800, height=2000, quality=100)

par(mfrow=c(5,2), mar=c(4.5,5,2,1), mgp=c(3,1,0))

#Breeding richness
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(0,120), ylab="Latitude", xaxt="n", axes=F, xlab="Species richness", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(lonlat.patterns[,2], breeding.richness.obs, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], breeding.richness.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(lonlat.patterns[,2], breeding.richness.bestguess, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], breeding.richness.bestguess, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("C", cex=2.5, side=3, line=-0.5, at=-15)
legend("bottomright", inset=0.1, bg="white", box.col="white", c("Empirical trend", "Simulated trend"), col=c("black","red"), cex=2.2, lwd=2)

plot(breeding.richness.obs, breeding.richness.bestguess, pch=20, xlim=c(0,210), ylim=c(0,210), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical richness", ylab="Simulated richness", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
#mtext("Empirical vs simulated", cex=1, side=3, line=0.15, at=100)
mtext("r=0.753", cex=2, side=1, line=-2.5, at=180)
mtext("D", cex=2.5, side=3, line=-0.5, at=-50)

#Non-breeding richness
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(0,120), ylab="Latitude", xaxt="n", axes=F, xlab="Species richness", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(lonlat.patterns[,2], nonbreeding.richness.obs, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], nonbreeding.richness.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(lonlat.patterns[,2], nonbreeding.richness.bestguess, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], nonbreeding.richness.bestguess, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("G", cex=2.5, side=3, line=-0.5, at=-15)

plot(nonbreeding.richness.obs, nonbreeding.richness.bestguess, pch=20, xlim=c(0,210), ylim=c(0,210), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical richness", ylab="Simulated richness", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
#mtext("Empirical vs simulated", cex=1, side=3, line=0.15, at=100)
mtext("r=0.403", cex=2, side=1, line=-2.5, at=180)
mtext("H", cex=2.5, side=3, line=-0.5, at=-50)

#Residents
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(0,400), ylab="Latitude", xaxt="n", axes=F, xlab="Species richness", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(lonlat.patterns[,2], resident.richness.obs, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], resident.richness.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(lonlat.patterns[,2], resident.richness.bestguess, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], resident.richness.bestguess, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("K", cex=2.5, side=3, line=-0.5, at=-49)

plot(resident.richness.obs, resident.richness.bestguess, pch=20, xlim=c(0,1000), ylim=c(0,800), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical richness", ylab="Simulated richness", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
#mtext("Empirical vs simulated", cex=1, side=3, line=0.15, at=100)
mtext("r=0.732", cex=2, side=3, line=-6, at=180)
mtext("L", cex=2.5, side=3, line=-0.5, at=-230)

#Difference in richness
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(-100,100), ylab="Latitude", xaxt="n", axes=F, xlab="Difference in richness", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(lonlat.patterns[,2], difference.richness.obs, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], difference.richness.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(lonlat.patterns[,2], difference.richness.bestguess, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[,2], difference.richness.bestguess, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("O", cex=2.5, side=3, line=-0.5, at=-125)

plot(difference.richness.obs, difference.richness.bestguess, pch=20, xlim=c(-190,190), ylim=c(-190,190), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical difference", ylab="Simulated difference", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
#mtext("Empirical vs simulated", cex=1, side=3, line=0.15, at=100)
mtext("r=0.734", cex=2, side=3, line=-6, at=-122)
mtext("P", cex=2.5, side=3, line=-0.5, at=-280)

#Proportion of migrants
tokeep = which(is.na(proportion.migrants.obs)==F & is.na(proportion.migrants.bestguess)==F)
prop.obs = proportion.migrants.obs[tokeep]
prop.bestguess = proportion.migrants.bestguess[tokeep]
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(0,1), ylab="Latitude", xaxt="n", axes=F, xlab="Proportion of migrants", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(lonlat.patterns[tokeep,2], prop.obs, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[tokeep,2], prop.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(lonlat.patterns[tokeep,2], prop.bestguess, kernel="normal", bandwidth=10)$y, ksmooth(lonlat.patterns[tokeep,2], prop.bestguess, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("S", cex=2.5, side=3, line=-0.5, at=-0.13)

plot(prop.obs, prop.bestguess, pch=20, xlim=c(0,1), ylim=c(0,1), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical proportion", ylab="Simulated proportion", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
#mtext("Empirical vs simulated", cex=1, side=3, line=0.15, at=100)
mtext("r=0.77", cex=2, side=3, line=-8, at=0.13)
mtext("T", cex=2.5, side=3, line=-0.5, at=-0.23)

dev.off()






## FIGURE S8 – Residuals ##

jpeg("FigS8.jpg", width=700, height=2000, quality=100)

par(mfrow=c(5,1), mar=c(1,1,1,1), mgp=c(3,1,0))

#Breeding richness residuals
breeding.richness.residuals <- breeding.richness.obs - breeding.richness.bestfit
rbPal <- colorRampPalette(c("blue", "white", "red3"))
datcol <- rbPal(7)[as.numeric(cut(breeding.richness.residuals, breaks=c(-160,-70,-40,-10,10,40,70,160)))]
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("A", cex=2.5, side=3, line=-3.5, at=-183)
mtext("Residuals for\nrichness in breeding migrants", cex=2.5, side=1, line=-3, at=20)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 70","40–70", "10–40", "-10–10", "-40–-10", "-70–-40", "< -70"), fill=rev(rbPal(7)), cex=2.6, border=rev(rbPal(7)))

#Non-breeding richness residuals
nonbreeding.richness.residuals <- nonbreeding.richness.obs - nonbreeding.richness.bestfit
rbPal <- colorRampPalette(c("blue", "white", "red3"))
datcol <- rbPal(7)[as.numeric(cut(nonbreeding.richness.residuals, breaks=c(-160,-70,-40,-10,10,40,70,160)))]
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("B", cex=2.5, side=3, line=-3.5, at=-183)
mtext("Residuals for\nrichness in non-breeding migrants", cex=2.5, side=1, line=-3, at=20)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 70","40–70", "10–40", "-10–10", "-40–-10", "-70–-40", "< -70"), fill=rev(rbPal(7)), cex=2.6, border=rev(rbPal(7)))

#Resident richness residuals
resident.richness.residuals <- resident.richness.obs - resident.richness.bestfit
rbPal <- colorRampPalette(c("blue", "white", "red3"))
datcol <- rbPal(7)[as.numeric(cut(resident.richness.residuals, breaks=c(-900,-300,-175,-50,50,175,300,900)))]
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("C", cex=2.5, side=3, line=-3.5, at=-183)
mtext("Residuals for\nrichness in residents", cex=2.5, side=1, line=-3, at=20)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 300","175–300", "50–175", "-50–50", "-175–-50", "-300–-175", "< -300"), fill=rev(rbPal(7)), cex=2.6, border=rev(rbPal(7)))

##Difference in richness
difference.richness.residuals <- difference.richness.obs - difference.richness.bestfit
rbPal <- colorRampPalette(c("blue", "white", "red3"))
datcol <- rbPal(7)[as.numeric(cut(difference.richness.residuals, breaks=c(-160,-70,-40,-10,10,40,70,160)))]
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("D", cex=2.5, side=3, line=-3.5, at=-183)
mtext("Residuals for\nseasonal difference in richness", cex=2.5, side=1, line=-3, at=20)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 70","40–70", "10–40", "-10–10", "-40–-10", "-70–-40", "< -70"), fill=rev(rbPal(7)), cex=2.6, border=rev(rbPal(7)))

##Proportion of migrants
proportion.migrants.residuals <- proportion.migrants.obs - proportion.migrants.bestfit
rbPal <- colorRampPalette(c("blue", "white", "red3"))
datcol <- rbPal(7)[as.numeric(cut(proportion.migrants.residuals, breaks=c(-1,-0.35,-0.2,-0.05,0.05,0.2,0.35,1)))]
plot(hexgrid2, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("E", cex=2.5, side=3, line=-3.5, at=-183)
mtext("Residuals for\nproportion of migrants", cex=2.5, side=1, line=-3, at=20)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 0.35","0.2–0.35", "0.05–0.2", "-0.05–0.05", "-0.2–-0.05", "-0.35–-0.2", "< -0.35"), fill=rev(rbPal(7)), cex=2.6, border=rev(rbPal(7)))

dev.off()





## FIGURE S9 – Along the simulation ##

lat_BR_WH <- vector()
lat_NB_WH <- vector()
for(i in 1:ncol(bestfit_rangesBR_WH)){
	lat_BR_WH[i] <- mean(lonlat.patterns[which(bestfit_rangesBR_WH[,i]==1),2])
	lat_NB_WH[i] <- mean(lonlat.patterns[which(bestfit_rangesNB_WH[,i]==1),2])
}
lat_BR_EH <- vector()
lat_NB_EH <- vector()
for(i in 1:ncol(bestfit_rangesBR_EH)){
	lat_BR_EH[i] <- mean(lonlat.patterns[which(bestfit_rangesBR_EH[,i]==1)+2323,2])
	lat_NB_EH[i] <- mean(lonlat.patterns[which(bestfit_rangesNB_EH[,i]==1)+2323,2])
}

iterationsWH <- (1:ncol(bestfit_rangesBR_WH)) / ncol(bestfit_rangesBR_WH)
iterationsEH <- (1:ncol(bestfit_rangesBR_EH)) / ncol(bestfit_rangesBR_EH)
iterations <- c(iterationsWH,iterationsEH)

migra.dist = c(migra.dist_WH, migra.dist_EH)
lat_BR <- c(lat_BR_WH, lat_BR_EH)
lat_NB <- c(lat_NB_WH, lat_NB_EH)
lat_res <- lat_BR[which(migra.dist==0)]
lat_BR <- lat_BR[which(migra.dist>0)]
lat_NB <- lat_NB[which(migra.dist>0)]
iterations_migr <- iterations[which(migra.dist>0)]
iterations_res <- iterations[which(migra.dist==0)]
migra.dist <- migra.dist[which(migra.dist>0)]

par(mfrow=c(2,1), mar=c(3,3,0.5,0.1), mgp=c(1.6,0.5,0))

plot(iterations_res, lat_res, pch=20, xlim=c(0,1), ylim=c(-60,80), col="green2", xaxt="n", axes=F, xlab="Simulation iteration", ylab="Latitude", cex.lab=1, cex=0.3)
axis(side=1)
axis(side=2)
legend("topleft", inset=0.01, bg="white", box.col="white", c("resident", "breeding", "non-breeding"), col=c("green2","red","blue"), cex=0.6, pch=20)
points(iterations_migr, lat_NB, pch=20, col="blue", cex=0.3)
points(iterations_migr, lat_BR, pch=20, col="red", cex=0.3)
mtext("A", cex=1.3, side=3, line=-1.2, at=0.5)

plot(iterations_migr, migra.dist, pch=20, xlim=c(0,1), ylim=c(0,14000), col="black", xaxt="n", axes=F, xlab="Simulation iteration", ylab="Migration distance", cex.lab=1, cex=0.5)
axis(side=1)
axis(side=2)
mtext("B", cex=1.3, side=3, line=-1.2, at=0.5)















