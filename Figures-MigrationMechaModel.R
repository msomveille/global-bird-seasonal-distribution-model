
####  R script to plot figures from outputs of the Seasonally Explicit Distributions Simulator (SEDS) model described in Somveille et al (2018) Nature Ecology & Evolution (doi.org/10.1038/s41559-018-0556-9)


##  Load libraries

library(rgdal)
library(fields)

setwd("global-bird-seasonal-distribution-model")  # set working directory


##  Load the observed presence/absence of species across the global grid of hexagons  ##

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


##  Load global grid of hexagons  ##

hexgrid <- readOGR("Hex_grid","isea3h7", verbose=FALSE)


##  Load simulated data outputed by the model  ##

load("Model-outputs-patterns.RData")		# Spatial patterns outputed by the model

# Western Hemisphere
bestfit_rangesBR_WH <- read.csv("bestfit_rangesBR_WH.csv", header=F)			# simulated breeding ranges outputed by the best-fit model
bestfit_rangesNB_WH <- read.csv("bestfit_rangesNB_WH.csv", header=F)			# simulated non-breeding ranges outputed by the best-fit model
bestguess_rangesBR_WH <- read.csv("bestguess_rangesBR_WH.csv", header=F)		# simulated breeding ranges outputed by the best guess model
bestguess_rangesNB_WH <- read.csv("bestguess_rangesNB_WH.csv", header=F)		# simulated non-breeding ranges outputed by the best guess model

# Eastern Hemisphere
bestfit_rangesBR_EH <- read.csv("bestfit_rangesBR_EH.csv", header=F)		# simulated breeding ranges outputed by the best-fit model
bestfit_rangesNB_EH <- read.csv("bestfit_rangesNB_EH.csv", header=F)			# simulated non-breeding ranges outputed by the best-fit model
bestguess_rangesBR_EH <- read.csv("bestguess_rangesBR_EH.csv", header=F)		# simulated breeding ranges outputed by the best guess model
bestguess_rangesNB_EH <- read.csv("bestguess_rangesNB_EH.csv", header=F)		# simulated non-breeding ranges outputed by the best guess model



##  FIGURE 2: contrast between empirical patterns in the global spatial distribution of terrestrial birds across seasons and the same patterns simulated through the overall best-fit model  ##

# Left-hand side part of the figure
 
jpeg("Fig2a.jpg", width=1400, height=2000, quality=100)

par(mfrow=c(5,2), mar=c(4.5,3,2,0.5), mgp=c(3,1,0))

##  Breeding richness  ##

# Plot empirical pattern
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(breeding.richness.obs, breaks=c(-1,25,50,75,100,125,150,175)))]
datcol[which(breeding.richness.obs ==0)] <- "white"
plot(hexgrid, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("A", cex=2.5, side=3, line=-0.5, at=-200)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 150","125–150", "100–125", "75–100", "50–75", "25–50", "< 25", "No species"), fill=c(rev(rbPal(7)),"white"), cex=2.2, border=c(rev(rbPal(7)),"white"))

# Plot simulated pattern
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(breeding.richness.bestfit, breaks=c(-1,25,50,75,100,125,150,175)))]
datcol[which(breeding.richness.bestfit==0)] <- "white"
plot(hexgrid, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("B", cex=2.5, side=3, line=-0.5, at=-200)


##  Non-breeding richness  ##

# Plot empirical pattern
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(nonbreeding.richness.obs, breaks=c(-1,25,50,75,100,125,150,175)))]
datcol[which(nonbreeding.richness.obs ==0)] <- "white"
plot(hexgrid, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("E", cex=2.5, side=3, line=-0.5, at=-200)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 150","125–150", "100–125", "75–100", "50–75", "25–50", "< 25", "No species"), fill=c(rev(rbPal(7)),"white"), cex=2.2, border=c(rev(rbPal(7)),"white"))

# Plot simulated pattern
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(nonbreeding.richness.bestfit, breaks=c(-1,25,50,75,100,125,150,175)))]
datcol[which(nonbreeding.richness.bestfit==0)] <- "white"
plot(hexgrid, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("F", cex=2.5, side=3, line=-0.5, at=-200)


##  Resident richness  ##

# Plot empirical pattern
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(resident.richness.obs, breaks=c(-1,100,200,300,400,500,600,1100)))]
datcol[which(resident.richness.obs ==0)] <- "white"
plot(hexgrid, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("I", cex=2.5, side=3, line=-0.5, at=-200)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 600","500–600", "400–500", "300–400", "200–300", "100–200", "< 100", "No species"), fill=c(rev(rbPal(7)),"white"), cex=2.2, border=c(rev(rbPal(7)),"white"))

# Plot simulated pattern
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(resident.richness.bestfit, breaks=c(-1,100,200,300,400,500,600,1100)))]
datcol[which(resident.richness.bestfit==0)] <- "white"
plot(hexgrid, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("J", cex=2.5, side=3, line=-0.5, at=-200)


##  Seasonal difference in richness  ##

# Plot empirical pattern
rbPal <- colorRampPalette(c("blue", "white", "red3"))
datcol <- rbPal(7)[as.numeric(cut(difference.richness.obs, breaks=c(-160,-110,-60,-10,10,60,110,160)))]
plot(hexgrid, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("M", cex=2.5, side=3, line=-0.5, at=-200)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 110","60–110", "10–60", "-10–10", "-60–-10", "-110–-60", "< -110"), fill=rev(rbPal(7)), cex=2.2, border=rev(rbPal(7)))

# Plot simulated pattern
rbPal <- colorRampPalette(c("blue", "white", "red3"))
datcol <- rbPal(7)[as.numeric(cut(difference.richness.bestfit, breaks=c(-160,-110,-60,-10,10,60,110,160)))]
plot(hexgrid, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("N", cex=2.5, side=3, line=-0.5, at=-200)


##  Proportion of migrants  ##

# Plot empirical pattern
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(proportion.migrants.obs, breaks=c(-0.1,0.15,0.30,0.45,0.60,0.75,0.90,1.1)))]
datcol[which(proportion.migrants.obs==0)] <- "white"
plot(hexgrid, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("Q", cex=2.5, side=3, line=-0.5, at=-200)
legend("bottomleft", inset=0, bg="grey", box.col="grey", c("> 0.9","0.75–0.9", "0.6–0.75", "0.45–0.6", "0.3–0.45", "0.15–0.3", "< 0.15", "No migrants"), fill=c(rev(rbPal(7)),"white"), cex=2.2, border=c(rev(rbPal(7)),"white"))

# Plot simulated pattern
rbPal <- colorRampPalette(c("yellow", "red3"))
datcol <- rbPal(7)[as.numeric(cut(proportion.migrants.bestfit, breaks=c(-0.1,0.15,0.30,0.45,0.60,0.75,0.90,1.1)))]
datcol[which(proportion.migrants.bestfit==0)] <- "white"
plot(hexgrid, col=datcol, border=datcol, bg="grey", ylim=c(-90,90))
mtext("R", cex=2.5, side=3, line=-0.5, at=-200)

dev.off()



# Right-hand side part of the figure

jpeg("Fig2b.jpg", width=800, height=2000, quality=100)

par(mfrow=c(5,2), mar=c(4.5,5,2,1), mgp=c(3,1,0))

##  Breeding richness  ##

# Plot latitudinal trends
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(0,120), ylab="Latitude", xaxt="n", axes=F, xlab="Species richness", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(hexgrid@data[,2], breeding.richness.obs, kernel="normal", bandwidth=10)$y, ksmooth(hexgrid@data[,2], breeding.richness.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(hexgrid@data[,2], breeding.richness.bestfit, kernel="normal", bandwidth=10)$y, ksmooth(hexgrid@data[,2], breeding.richness.bestfit, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("C", cex=2.5, side=3, line=-0.5, at=-15)
legend("bottomright", inset=0.1, bg="white", box.col="white", c("Empirical trend", "Simulated trend"), col=c("black","red"), cex=2.2, lwd=2)

# Plot relationship between empirical and simulated data
plot(breeding.richness.obs, breeding.richness.bestfit, pch=20, xlim=c(0,210), ylim=c(0,210), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical richness", ylab="Simulated richness", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
mtext("r=0.867", cex=2, side=1, line=-2.5, at=180)
mtext("D", cex=2.5, side=3, line=-0.5, at=-50)


##  Non-breeding richness  ##

# Plot latitudinal trends
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(0,120), ylab="Latitude", xaxt="n", axes=F, xlab="Species richness", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(hexgrid@data[,2], nonbreeding.richness.obs, kernel="normal", bandwidth=10)$y, ksmooth(hexgrid@data[,2], nonbreeding.richness.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(hexgrid@data[,2], nonbreeding.richness.bestfit, kernel="normal", bandwidth=10)$y, ksmooth(hexgrid@data[,2], nonbreeding.richness.bestfit, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("G", cex=2.5, side=3, line=-0.5, at=-15)

# Plot relationship between empirical and simulated data
plot(nonbreeding.richness.obs, nonbreeding.richness.bestfit, pch=20, xlim=c(0,210), ylim=c(0,210), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical richness", ylab="Simulated richness", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
mtext("r=0.748", cex=2, side=1, line=-2.5, at=180)
mtext("H", cex=2.5, side=3, line=-0.5, at=-50)


##  Resident richness  ##

# Plot latitudinal trends
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(0,400), ylab="Latitude", xaxt="n", axes=F, xlab="Species richness", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(hexgrid@data[,2], resident.richness.obs, kernel="normal", bandwidth=10)$y, ksmooth(hexgrid@data[,2], resident.richness.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(hexgrid@data[,2], resident.richness.bestfit, kernel="normal", bandwidth=10)$y, ksmooth(hexgrid@data[,2], resident.richness.bestfit, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("K", cex=2.5, side=3, line=-0.5, at=-49)

# Plot relationship between empirical and simulated data
plot(resident.richness.obs, resident.richness.bestfit, pch=20, xlim=c(0,1000), ylim=c(0,800), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical richness", ylab="Simulated richness", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
mtext("r=0.696", cex=2, side=3, line=-6, at=180)
mtext("L", cex=2.5, side=3, line=-0.5, at=-230)


##  Seasonal difference in richness  ##

# Plot latitudinal trends
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(-100,100), ylab="Latitude", xaxt="n", axes=F, xlab="Difference in richness", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(hexgrid@data[,2], difference.richness.obs, kernel="normal", bandwidth=10)$y, ksmooth(hexgrid@data[,2], difference.richness.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(hexgrid@data[,2], difference.richness.bestfit, kernel="normal", bandwidth=10)$y, ksmooth(hexgrid@data[,2], difference.richness.bestfit, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("O", cex=2.5, side=3, line=-0.5, at=-125)

# Plot relationship between empirical and simulated data
plot(difference.richness.obs, difference.richness.bestfit, pch=20, xlim=c(-190,190), ylim=c(-190,190), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical difference", ylab="Simulated difference", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
mtext("r=0.871", cex=2, side=3, line=-6, at=-122)
mtext("P", cex=2.5, side=3, line=-0.5, at=-280)


##  Proportion of migrants  ##

# Plot latitudinal trends
tokeep = which(is.na(proportion.migrants.obs)==F & is.na(proportion.migrants.bestfit)==F)
prop.obs = proportion.migrants.obs[tokeep]
prop.bestfit = proportion.migrants.bestfit[tokeep]
plot(NULL, pch=20, col="white", ylim=c(-90,90), xlim=c(0,1), ylab="Latitude", xaxt="n", axes=F, xlab="Proportion of migrants", cex.lab=2.5)
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
lines(ksmooth(hexgrid@data[tokeep,2], prop.obs, kernel="normal", bandwidth=10)$y, ksmooth(hexgrid@data[tokeep,2], prop.obs, kernel="normal", bandwidth=10)$x, lwd=2, col="black")
lines(ksmooth(hexgrid@data[tokeep,2], prop.bestfit, kernel="normal", bandwidth=10)$y, ksmooth(hexgrid@data[tokeep,2], prop.bestfit, kernel="normal", bandwidth=10)$x, lwd=2, col="red")
mtext("S", cex=2.5, side=3, line=-0.5, at=-0.13)

# Plot relationship between empirical and simulated data
plot(prop.obs, prop.bestfit, pch=20, xlim=c(0,1), ylim=c(0,1), cex=1.5, xaxt="n", axes=F, col="dark grey", xlab="Empirical proportion", ylab="Simulated proportion", cex.lab=2.5)
abline(0,1, col="red")
axis(side=1, labels=T, cex.axis=2)
axis(side=2, labels=T, cex.axis=2)
mtext("r=0.837", cex=2, side=3, line=-6, at=0.13)
mtext("T", cex=2.5, side=3, line=-0.5, at=-0.23)

dev.off()





##  FIGURE 3: Latitudinal trends for empirical and simulated patterns in the global spatial distribution of terrestrial birds across seasons  ##

par(mfrow=c(5,1), mar=c(2,4.5,2,0.5), mgp=c(1.7,0.5,0), xpd=T)

# Breeding richness
plot(NULL, pch=20, col="white", xlim=c(-65,90), ylim=c(0,235), xlab="", xaxt="n", axes=F, ylab="Richness in\nbreeding migrants", cex.lab=1.15, main="")
axis(side=1, labels=T, at=c(-60,-30,0,30,60,90), cex.axis=1)
axis(side=2, labels=T, at=c(0,60,120,180), cex.axis= 0.92)
lines(ksmooth(lonlat.patterns[,2], breeding.richness.noenergy, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], breeding.richness.noenergy, kernel="normal", bandwidth=10)$y, lwd=1.3, col="pink")			# Latitudinal trends generated by the best-fit model without energetic costs
lines(ksmooth(lonlat.patterns[,2], breeding.richness.nomigra, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], breeding.richness.nomigra, kernel="normal", bandwidth=10)$y, lwd=1.3, col="dark grey")		# Latitudinal trends generated by the best-fit model without migration cost
lines(ksmooth(lonlat.patterns[,2], breeding.richness.norepro, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], breeding.richness.norepro, kernel="normal", bandwidth=10)$y, lwd=1.3, col="cyan2")			# Latitudinal trends generated by the best-fit model without reproduction cost
lines(ksmooth(lonlat.patterns[,2], breeding.richness.bestguess, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], breeding.richness.bestguess, kernel="normal", bandwidth=10)$y, lwd=1.3, col="blue")			# Latitudinal trends generated by the best guess model
lines(ksmooth(lonlat.patterns[,2], breeding.richness.nothermo, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], breeding.richness.nothermo, kernel="normal", bandwidth=10)$y, lwd=1.3, col="chartreuse3")	# Latitudinal trends generated by the best-fit model without thermoregulation cost
lines(ksmooth(lonlat.patterns[,2], breeding.richness.obs, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], breeding.richness.obs, kernel="normal", bandwidth=10)$y, lwd=1.3, lty=3, col="black")					# Latitudinal trends for the empirical pattern
lines(ksmooth(lonlat.patterns[,2], breeding.richness.bestfit, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], breeding.richness.bestfit, kernel="normal", bandwidth=10)$y, lwd=1.3, col="orange")		# Latitudinal trends generated by the overall best-fit model
mtext("a", cex=1.2, at=-107, line=0.7)

# Non-breeding richness
plot(NULL, pch=20, col="white", xlim=c(-65,90), ylim=c(0,200), xlab="", xaxt="n", axes=F, ylab="Richness in\nnon-breeding migrants", cex.lab=1.15, main="")
axis(side=1, labels=T, at=c(-60,-30,0,30,60,90), cex.axis=1)
axis(side=2, labels=T, at=c(0,70,140,210), cex.axis= 0.92)
lines(ksmooth(lonlat.patterns[,2], nonbreeding.richness.noenergy, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], nonbreeding.richness.noenergy, kernel="normal", bandwidth=10)$y, lwd=1.3, col="pink")
lines(ksmooth(lonlat.patterns[,2], nonbreeding.richness.nomigra, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], nonbreeding.richness.nomigra, kernel="normal", bandwidth=10)$y, lwd=1.3, col="dark grey")
lines(ksmooth(lonlat.patterns[,2], nonbreeding.richness.norepro, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], nonbreeding.richness.norepro, kernel="normal", bandwidth=10)$y, lwd=1.3, col="cyan2")
lines(ksmooth(lonlat.patterns[,2], nonbreeding.richness.bestguess, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], nonbreeding.richness.bestguess, kernel="normal", bandwidth=10)$y, lwd=1.3, col="blue")
lines(ksmooth(lonlat.patterns[,2], nonbreeding.richness.nothermo, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], nonbreeding.richness.nothermo, kernel="normal", bandwidth=10)$y, lwd=1.3, col="chartreuse3")
lines(ksmooth(lonlat.patterns[,2], nonbreeding.richness.obs, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], nonbreeding.richness.obs, kernel="normal", bandwidth=10)$y, lwd=1.3, lty=3, col="black")
lines(ksmooth(lonlat.patterns[,2], nonbreeding.richness.bestfit, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], nonbreeding.richness.bestfit, kernel="normal", bandwidth=10)$y, lwd=1.3, col="orange")
mtext("b", cex=1.2, at=-107, line=0.7)

# Resident richness
plot(NULL, pch=20, col="white", xlim=c(-65,90), ylim=c(0,400), xlab="", xaxt="n", axes=F, ylab="Richness in\nresidents", cex.lab=1.15, main="")
axis(side=1, labels=T, at=c(-60,-30,0,30,60,90), cex.axis=1)
axis(side=2, labels=T, at=c(0,130,260,390), cex.axis=0.86)
legend("topright", inset=c(0,-0.40), bg="white", box.col="white", c("Empirical", "Best fit", "Best guess", "Thermo zero", "Repro zero", "Migra zero", "Costs zero"), col=c("black","orange", "blue", "chartreuse3", "cyan2", "dark grey", "pink"), cex=0.75, lwd=1.15, lty=c(3,1,1,1,1,1,1))
lines(ksmooth(lonlat.patterns[,2], resident.richness.noenergy, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], resident.richness.noenergy, kernel="normal", bandwidth=10)$y, lwd=1.3, col="pink")
lines(ksmooth(lonlat.patterns[,2], resident.richness.nomigra, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], resident.richness.nomigra, kernel="normal", bandwidth=10)$y, lwd=1.3, col="dark grey")
lines(ksmooth(lonlat.patterns[,2], resident.richness.norepro, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], resident.richness.norepro, kernel="normal", bandwidth=10)$y, lwd=1.3, col="cyan2")
lines(ksmooth(lonlat.patterns[,2], resident.richness.bestguess, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], resident.richness.bestguess, kernel="normal", bandwidth=10)$y, lwd=1.3, col="blue")
lines(ksmooth(lonlat.patterns[,2], resident.richness.nothermo, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], resident.richness.nothermo, kernel="normal", bandwidth=10)$y, lwd=1.3, col="chartreuse3")
lines(ksmooth(lonlat.patterns[,2], resident.richness.obs, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], resident.richness.obs, kernel="normal", bandwidth=10)$y, lwd=1.3, lty=3, col="black")
lines(ksmooth(lonlat.patterns[,2], resident.richness.bestfit, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], resident.richness.bestfit, kernel="normal", bandwidth=10)$y, lwd=1.3, col="orange")
mtext("c", cex=1.2, at=-107, line=0.7)

# Seasonal difference in richness
plot(NULL, pch=20, col="white", xlim=c(-65,90), ylim=c(-140,100), xlab="", xaxt="n", axes=F, ylab="Seasonal difference\nin richness", cex.lab=1.15, main="")
axis(side=1, labels=T, at=c(-60,-30,0,30,60,90), cex.axis=1)
axis(side=2, labels=T, at=c(-100,-50, 0, 50, 100), cex.axis=0.95)
lines(ksmooth(lonlat.patterns[,2], difference.richness.noenergy, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], difference.richness.noenergy, kernel="normal", bandwidth=10)$y, lwd=1.3, col="pink")
lines(ksmooth(lonlat.patterns[,2], difference.richness.nomigra, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], difference.richness.nomigra, kernel="normal", bandwidth=10)$y, lwd=1.3, col="dark grey")
lines(ksmooth(lonlat.patterns[,2], difference.richness.norepro, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], difference.richness.norepro, kernel="normal", bandwidth=10)$y, lwd=1.3, col="cyan2")
lines(ksmooth(lonlat.patterns[,2], difference.richness.bestguess, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], difference.richness.bestguess, kernel="normal", bandwidth=10)$y, lwd=1.3, col="blue")
lines(ksmooth(lonlat.patterns[,2], difference.richness.nothermo, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], difference.richness.nothermo, kernel="normal", bandwidth=10)$y, lwd=1.3, col="chartreuse3")
lines(ksmooth(lonlat.patterns[,2], difference.richness.obs, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], difference.richness.obs, kernel="normal", bandwidth=10)$y, lwd=1.3, lty=3, col="black")
lines(ksmooth(lonlat.patterns[,2], difference.richness.bestfit, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[,2], difference.richness.bestfit, kernel="normal", bandwidth=10)$y, lwd=1.3, col="orange")
mtext("d", cex=1.2, at=-107, line=1.8)

# Proportion of migrants
plot(NULL, pch=20, col="white", xlim=c(-65,90), ylim=c(0,1), xlab="", xaxt="n", axes=F, ylab="Proportion\nof migrants", cex.lab=1.15, main="")
axis(side=1, labels=T, at=c(-60,-30,0,30,60,90), cex.axis=1)
axis(side=2, labels=T, at=c(0,0.5,1), cex.axis=1)
tokeep = which(is.na(proportion.migrants.noenergy)==F)
prop.noenergy = proportion.migrants.nomigra[tokeep]
lines(ksmooth(lonlat.patterns[tokeep,2], prop.noenergy, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[tokeep,2], prop.noenergy, kernel="normal", bandwidth=10)$y, lwd=1.3, col="pink")
tokeep = which(is.na(proportion.migrants.nomigra)==F)
prop.nomigra = proportion.migrants.nomigra[tokeep]
lines(ksmooth(lonlat.patterns[tokeep,2], prop.nomigra, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[tokeep,2], prop.nomigra, kernel="normal", bandwidth=10)$y, lwd=1.3, col="dark grey")
tokeep = which(is.na(proportion.migrants.norepro)==F)
prop.norepro = proportion.migrants.norepro[tokeep]
lines(ksmooth(lonlat.patterns[tokeep,2], prop.norepro, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[tokeep,2], prop.norepro, kernel="normal", bandwidth=10)$y, lwd=1.3, col="cyan2")
tokeep = which(is.na(proportion.migrants.bestguess)==F)
prop.bestguess = proportion.migrants.bestguess[tokeep]
lines(ksmooth(lonlat.patterns[tokeep,2], prop.bestguess, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[tokeep,2], prop.bestguess, kernel="normal", bandwidth=10)$y, lwd=1.3, col="blue")
tokeep = which(is.na(proportion.migrants.nothermo)==F)
prop.nothermo = proportion.migrants.nothermo[tokeep]
lines(ksmooth(lonlat.patterns[tokeep,2], prop.nothermo, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[tokeep,2], prop.nothermo, kernel="normal", bandwidth=10)$y, lwd=1.3, col="chartreuse3")
tokeep = which(is.na(proportion.migrants.obs)==F)
prop.obs = proportion.migrants.obs[tokeep]
lines(ksmooth(lonlat.patterns[tokeep,2], prop.obs, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[tokeep,2], prop.obs, kernel="normal", bandwidth=10)$y, lwd=1.3, lty=3, col="black")
tokeep = which(is.na(proportion.migrants.bestfit)==F)
prop.bestfit = proportion.migrants.bestfit[tokeep]
lines(ksmooth(lonlat.patterns[tokeep,2], prop.bestfit, kernel="normal", bandwidth=10)$x, ksmooth(lonlat.patterns[tokeep,2], prop.bestfit, kernel="normal", bandwidth=10)$y, lwd=1.3, col="orange")
mtext("e", cex=1.2, at=-107, line=0.5)
mtext("Latitude", side=1, cex=0.8, at=10, line=1.2)





##  FIGURE 4: Predictive value of the overall best-fit model contrasted with the best guess model and the null models  ##

par(mfrow=c(1,1), mar=c(3,0.5,0.1,0.1), mgp=c(1.5,0.5,0))

# Earth Mover's Distance (EMD) generated by randomly selected runs of the null model that does not integrate energetic considerations
EMDrandom = c(15070.04877862, 14874.24106341, 14848.6339411, 15128.38885583, 14718.56616819, 14811.93059197, 14286.96650175, 14936.46057012, 14319.97380384, 14619.64875535, 14938.40594741, 15255.37661177, 15122.1796363, 14799.53505901, 14830.38351451, 15042.20544426, 14607.97031363, 15131.36282917, 15161.70910288, 14588.57941737, 15152.97933711,  14628.04664144,  15142.52945065,  15077.37874877, 14887.5091581,14752.09731643,  15096.52487758,  15421.71639309,  14957.63111244, 14812.64824094,15150.262564  ,  14875.1368278 ,  15358.53082714,  14933.89935371, 14466.93908091,  15235.77118226,  15251.83772329,  14892.98154497, 15152.16477366,  14521.77461875,  14917.38667495,  15030.37183644, 14659.95859803,  14854.59014217,  14799.06598113,  14883.25767138, 14902.41998327,  14802.84578178,  15170.21958684,  14752.77211437, 14951.04198296,  14781.80296817,  14969.79792414,  14672.14018681, 15743.14545695,  14862.4939531 ,  15168.63940031,  14773.46395046, 15188.51281814,  14553.97381536,  14212.67829694,  14407.73396212, 14983.63517247,  14793.77282927,  15111.20908505,  14685.18330175, 14931.51372845,  14590.26302761,  15176.13129261,  14855.50145836, 14441.07842778,  14568.75461929,  14901.47193349,  14800.68740383, 14884.89572168,  15048.38188962,  14586.7864773 ,  15213.56061367, 14841.75650165,  15054.55683057,  15161.08539294,  14909.14784689, 14732.01415058,  15590.13914395,  14589.13799519,  14791.08521445, 15657.41720089,  14759.54391251,  15030.44695098,  14640.21472899, 14828.08698223,  15205.59227586,  15018.20965025,  14812.10644452, 14951.05739924,  15268.19428716,  14856.60248414,  15057.60413885, 14003.22709452,  15016.18668042,  14868.88530811,  15258.92538405, 15227.158873  ,  15301.46651266,  14874.12790364,  14651.22009083)
hist(EMDrandom, xlim=c(0,17000), ylim=c(0,45), main="", xlab="Earth Mover's Distance between empirical and simulated patterns       ", ylab="", axes=F, col="grey", border="dark grey")
axis(side=1)
text(x= mean(EMDrandom), y=37, "Null\ndistribution", col="dark grey", cex=1.1)

# Earth Mover's Distance (EMD) generated by the best-fit model without thermoregulation cost
EMDnothermo = 477.951472836
lines(c(EMDnothermo, EMDnothermo), c(0,23), col="dark grey", lty=2, cex=2) # orange2
text(x= EMDnothermo+100, y=27.7, "No thermo cost", col="dark grey", pos=3, srt=90, cex=1.1)

# Earth Mover's Distance (EMD) generated by the best-fit model without reproduction cost
EMDnorepro = 1341.17693267
lines(c(EMDnorepro, EMDnorepro), c(0,10), col="dark grey", lty=2, cex=2) # orange2
text(x= EMDnorepro+100, y=17, "No reproduction cost", col="dark grey", pos=3, srt=90, cex=1.1)

# Earth Mover's Distance (EMD) generated by the best-fit model without migration cost
EMDnomigra = 11959.4451043
lines(c(EMDnomigra, EMDnomigra), c(0,10), col="dark grey", lty=2, cex=2) # orange2
text(x= EMDnomigra+100, y=16, "No migration cost", col="dark grey", pos=3, srt=90, cex=1.1) # orange2

# Earth Mover's Distance (EMD) generated by the overall best-fit model
EMDbestfit = 410.67679711
lines(c(EMDbestfit, EMDbestfit), c(0,15), col="black", lty=1, cex=2) # orange2
text(x= EMDbestfit+60, y=17, "Best fit", col="black", pos=3, srt=90, cex=1.1) # orange2

# Earth Mover's Distance (EMD) generated by the best guess model
EMDbestguess = 558
lines(c(EMDbestguess, EMDbestguess), c(0,5), col="black", lty=1, cex=2) # blue
text(x= EMDbestguess+190, y=8.5, "Best guess", col="black", pos=3, srt=90, cex=1.1) # blue

# Earth Mover's Distance (EMD) generated by the best-fit model without energetic costs
EMDnoenergy = 12271.9526771
lines(c(EMDnoenergy, EMDnoenergy), c(0,10), col="dark grey", lty=1, cex=2) # green3
text(x= EMDnoenergy+140, y=16.4, "No energetic costs", col="dark grey", pos=3, srt=90, cex=1.1) # green3




##  FIGURE S11: Species generated along the simulation corresponding to the best-fit model  ##

# Compute the mean latitude and longitude of the simulated ranges

# Western Hemisphere
lat_BR_WH <- vector()
lat_NB_WH <- vector()
lon_BR_WH <- vector()
lon_NB_WH <- vector()
for(i in 1:ncol(bestfit_rangesBR_WH)){
	lat_BR_WH[i] <- mean(hexgrid@data[which(bestfit_rangesBR_WH[,i]==1),2])		# mean latitude of the simulated breeding ranges
	lat_NB_WH[i] <- mean(hexgrid@data[which(bestfit_rangesNB_WH[,i]==1),2])		# mean latitude of the simulated non-breeding ranges
	lon_BR_WH[i] <- mean(hexgrid@data[which(bestfit_rangesBR_WH[,i]==1),1])		# mean longitude of the simulated breeding ranges
	lon_NB_WH[i] <- mean(hexgrid@data[which(bestfit_rangesNB_WH[,i]==1),1])		# mean longitude of the simulated non-breeding ranges
}

# Eastern Hemisphere
lat_BR_EH <- vector()
lat_NB_EH <- vector()
lon_BR_EH <- vector()
lon_NB_EH <- vector()
for(i in 1:ncol(bestfit_rangesBR_EH)){
	lat_BR_EH[i] <- mean(hexgrid@data[which(bestfit_rangesBR_EH[,i]==1)+2323,2])
	lat_NB_EH[i] <- mean(hexgrid@data[which(bestfit_rangesNB_EH[,i]==1)+2323,2])
	lon_BR_EH[i] <- mean(hexgrid@data[which(bestfit_rangesBR_EH[,i]==1)+2323,1])
	lon_NB_EH[i] <- mean(hexgrid@data[which(bestfit_rangesNB_EH[,i]==1)+2323,1])
}


# Rescale the number of simulation iterations in order to combine Western and Eastern Hemispheres
iterationsWH <- (1:ncol(bestfit_rangesBR_WH)) / ncol(bestfit_rangesBR_WH)
iterationsEH <- (1:ncol(bestfit_rangesBR_EH)) / ncol(bestfit_rangesBR_EH)
iterations <- c(iterationsWH,iterationsEH)


# Compute the migration distances (i.e. distance separating the simulated breeding and non-breeding ranges)
migra.dist_WH <- vector()
for(i in 1:length(lat_BR_WH)){
	migra.dist_WH[i] <- distance(deg2rad(lon_BR_WH[i]), deg2rad(lat_BR_WH[i]), deg2rad(lon_NB_WH[i]), deg2rad(lat_NB_WH[i]))
}
migra.dist_EH <- vector()
for(i in 1:length(lat_BR_EH)){
	migra.dist_EH[i] <- distance(deg2rad(lon_BR_EH[i]), deg2rad(lat_BR_EH[i]), deg2rad(lon_NB_EH[i]), deg2rad(lat_NB_EH[i]))
}
migra.dist = c(migra.dist_WH, migra.dist_EH)
lat_BR <- c(lat_BR_WH, lat_BR_EH)
lat_NB <- c(lat_NB_WH, lat_NB_EH)
lat_res <- lat_BR[which(migra.dist==0)]
lat_BR <- lat_BR[which(migra.dist>0)]
lat_NB <- lat_NB[which(migra.dist>0)]
iterations_migr <- iterations[which(migra.dist>0)]
iterations_res <- iterations[which(migra.dist==0)]
migra.dist <- migra.dist[which(migra.dist>0)]


# Plot 
par(mfrow=c(2,1), mar=c(3,3,0.5,0.1), mgp=c(1.6,0.5,0))

plot(iterations_res, lat_res, pch=20, xlim=c(0,1), ylim=c(-60,80), col="green2", xaxt="n", axes=F, xlab="Simulation iteration", ylab="Latitude", cex.lab=1, cex=0.3)
axis(side=1)
axis(side=2)
#legend("topleft", inset=0.01, bg="white", box.col="white", c("resident", "breeding", "non-breeding"), col=c("green2","red","blue"), cex=0.6, pch=20)
points(iterations_migr, lat_NB, pch=20, col="blue", cex=0.3)
points(iterations_migr, lat_BR, pch=20, col="red", cex=0.3)
mtext("A", cex=1.3, side=3, line=-1.2, at=0.5)

plot(iterations_migr, migra.dist, pch=20, xlim=c(0,1), ylim=c(0,14000), col="black", xaxt="n", axes=F, xlab="Simulation iteration", ylab="Migration distance", cex.lab=1, cex=0.5)
axis(side=1)
axis(side=2)
mtext("B", cex=1.3, side=3, line=-1.2, at=0.5)



