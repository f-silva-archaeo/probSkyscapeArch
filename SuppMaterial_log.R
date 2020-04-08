#### A probabilistic framework and significance test for the analysis of structural orientations in skyscape archaeology
#### Fabio Silva ::: Bournemouth University ::: 2020 ::: fsilva@bournemouth.ac.uk

###########################################################################
#### Supporting Material

setwd("~/Dropbox/Archaeology/Projects/probSkyArch/2019 11 Paper/Revision")
devtools::install_github('f-silva-archaeo/skyscapeR') ## requires development version of skyscapeR
source('src.R')



# Figure S1 ---------------------------------------------------------------
par(mfrow=c(1,3))

hor.flat <- createHor(c(0,360), c(0,0), 0.5, c(35,0,1000))
eq <- sky.objects(0,epoch=2000, col='blue', lty=1, lwd=2)
plot(hor.flat, xlim=c(70,110), ylim=c(-2,10), obj=eq, refraction=F, show.az=T)
points(90,0, cex=3)
text(70, 9, labels='Latitude 35ºN \n Azimuth 90º', font=2, pos=4, cex=1.5)

hor <- createHor(c(0,45,80,90,100,135,360), c(0,0,0,6.66,0,0,0), 0.5, c(35,0,1000), smooth=T, .scale=1000)
plot(hor, xlim=c(70,110), ylim=c(-2,10), obj=eq, refraction=F, show.az=T)
points(93,hor2alt(hor,93), cex=3)
text(70, 9, labels='Latitude 35ºN \n Azimuth 93º', font=2, pos=4, cex=1.5)

hor <- createHor(c(0,45,80,90,100,135,360), c(0,0,0,6.66,0,0,0), 0.5, c(65,0,1000), smooth=T, .scale=1000)
plot(hor, xlim=c(70,110), ylim=c(-2,10), obj=eq, refraction=F, show.az=T)
points(95.9,hor2alt(hor,95.9), cex=3)
text(70, 9, labels='Latitude 65ºN \n Azimuth 96º', font=2, pos=4, cex=1.5)

dev.print(png, filename = 'FigS1.png', res=300, width=10, height=3, units='in')



# Figure S2 ----------------------------------------------------------------
lat <- 40
hor <- createHor(c(0,360), c(0,0), 0.5, c(lat,0,1000), 'Flat')

par(mfrow=c(2,2))
dec1 <- coordtrans(pdf='normal', az=90, unc=5, hor=hor, refraction=F)
plot(dec1, show.az = F, col=MESS::col.alpha('blue',0.4))
text(par('usr')[1],0.93*par('usr')[4], 'A', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Az = 90º', font=2, cex=1.1)

dec2 <- coordtrans(pdf='normal', az=120, unc=5, hor=hor, refraction=F)
plot(dec2, show.az = F, col=MESS::col.alpha('blue',0.4), xlim=az2dec(120, hor, 0, refraction=F)+c(-15,15))
lines(dec1$data[[1]]$dec$x + az2dec(120, hor, 0, refraction=F), dec1$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'B', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Az = 120º', font=2, cex=1.1)

dec3 <- coordtrans(pdf='normal', az=160, unc=5, hor=hor, refraction=F)
plot(dec3, show.az = F, col=MESS::col.alpha('blue',0.4), xlim=az2dec(160, hor, 0, refraction=F)+c(-15,15))
lines(dec1$data[[1]]$dec$x + az2dec(160, hor, 0, refraction=F), dec1$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'C', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Az = 160º', font=2, cex=1.1)

dec4 <- coordtrans(pdf='normal', az=180, unc=5, hor=hor, refraction=F)
plot(dec4, show.az = F, col=MESS::col.alpha('blue',0.4), xlim=az2dec(180, hor, 0, refraction=F)+c(-15,15))
lines(dec1$data[[1]]$dec$x + az2dec(180, hor, 0, refraction=F), dec1$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'D', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Az = 180º', font=2, cex=1.1)

dev.print(png, filename = 'FigS2.png', res=300, width=8, height=8, units='in')



# Figure S3 ----------------------------------------------------------------
lat <- 40
hor <- createHor(c(0,360), c(0,0), 0.5, c(lat,0,1000), 'Flat')

par(mfrow=c(2,2))
dec1 <- coordtrans(pdf='uniform', az=90, unc=5, hor=hor, refraction=F)
plot(dec1, show.az = F, col=MESS::col.alpha('blue',0.4))
text(par('usr')[1],0.93*par('usr')[4], 'A', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Az = 90º', font=2, cex=1.1)

dec2 <- coordtrans(pdf='uniform', az=120, unc=5, hor=hor, refraction=F)
plot(dec2, show.az = F, col=MESS::col.alpha('blue',0.4), xlim=az2dec(120, hor, 0, refraction=F)+c(-15,15))
lines(dec1$data[[1]]$dec$x + az2dec(120, hor, 0, refraction=F), dec1$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'B', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Az = 120º', font=2, cex=1.1)

dec3 <- coordtrans(pdf='uniform', az=160, unc=5, hor=hor, refraction=F)
plot(dec3, show.az = F, col=MESS::col.alpha('blue',0.4), xlim=az2dec(160, hor, 0, refraction=F)+c(-15,15))
lines(dec1$data[[1]]$dec$x + az2dec(160, hor, 0, refraction=F), dec1$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'C', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Az = 160º', font=2, cex=1.1)

dec4 <- coordtrans(pdf='uniform', az=180, unc=5, hor=hor, refraction=F)
plot(dec4, show.az = F, col=MESS::col.alpha('blue',0.4), xlim=az2dec(180, hor, 0, refraction=F)+c(-15,15))
lines(dec1$data[[1]]$dec$x + az2dec(180, hor, 0, refraction=F), dec1$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'D', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Az = 180º', font=2, cex=1.1)

dev.print(png, filename = 'FigS3.png', res=300, width=8, height=8, units='in')



# Figure S4 ----------------------------------------------------------------
par(mfrow=c(2,2))

lat <- 40
hor <- createHor(c(0,90-7,97,270,360),c(15,15,0,0,15), 0.5, c(lat,0,1000), 'Downwards slope', smooth=T, .scale=100)
dec.down <- coordtrans(pdf='normal', az=90, unc=2, hor=hor, refraction=F)
hor <- createHor(c(0,360),rep(hor2alt(hor,90),2), 0.5, c(lat,0,1000))
dec.flat <- coordtrans(pdf='normal', az=90, unc=2, hor=hor, refraction=F)
plot(dec.down, show.az=F, col=MESS::col.alpha('blue',0.4), xlim=c(-7,17))
lines(dec.flat$data[[1]]$dec$x, dec.flat$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'A', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Downwards slope', font=2, cex=1.1)

hor <- createHor(c(0,90-7,97,360),c(0,0,15,15), 0.5, c(lat,0,1000), 'Upwards slope', smooth=T, .scale=100)
dec.up <- coordtrans(pdf='normal', az=90, unc=2, hor=hor, refraction=F)
hor <- createHor(c(0,360),rep(hor2alt(hor,90), 2), 0.5, c(lat,0,1000))
dec.flat <- coordtrans(pdf='normal', az=90, unc=2, hor=hor, refraction=F)
plot(dec.up, show.az=F, col=MESS::col.alpha('blue',0.4), xlim=c(-7,17))
lines(dec.flat$data[[1]]$dec$x, dec.flat$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'B', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Upwards slope', font=2, cex=1.1)

hor <- createHor(c(0,70,90,110,360),c(0,0,15,0,0), 0.5, c(lat,0,1000), 'Peak', smooth=T, .scale=100)
dec.peak <- coordtrans(pdf='normal', az=90, unc=2, hor=hor, refraction=F)
hor <- createHor(c(0,360),rep(hor2alt(hor,90),2), 0.5, c(lat,0))
dec.flat <- coordtrans(pdf='normal', az=90, unc=2, hor=hor, refraction=F)
plot(dec.peak, show.az=F, col=MESS::col.alpha('blue',0.4), xlim=c(-7,17))
lines(dec.flat$data[[1]]$dec$x, dec.flat$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'C', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Peak', font=2, cex=1.1)

hor <- createHor(c(0,70,90,110,360),c(15,15,0,15,15), 0.5, c(lat,0,1000), 'Notch', smooth=T, .scale=100)
dec.notch <- coordtrans(pdf='normal', az=90, unc=2, hor=hor, refraction=F)
hor <- createHor(c(0,360),rep(hor2alt(hor,90),2), 0.5, c(lat,0,1000))
dec.flat <- coordtrans(pdf='normal', az=90, unc=2, hor=hor, refraction=F)
plot(dec.notch, show.az=F, col=MESS::col.alpha('blue',0.4), xlim=c(-7,17))
lines(dec.flat$data[[1]]$dec$x, dec.flat$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'D', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Notch', font=2, cex=1.1)

dev.print(png, filename = 'FigS4.png', res=300, width=8, height=8, units='in')



# Figure S5 ----------------------------------------------------------------
par(mfrow=c(2,2))

lat <- 40
hor <- createHor(c(0,90-7,97,270,360),c(15,15,0,0,15), 0.5, c(lat,0,1000), 'Downwards slope', smooth=T, .scale=100)
dec.down <- coordtrans(pdf='uniform', az=90, unc=2, hor=hor, refraction=F)
hor <- createHor(c(0,360),rep(hor2alt(hor,90),2), 0.5, c(lat,0,1000))
dec.flat <- coordtrans(pdf='uniform', az=90, unc=2, hor=hor, refraction=F)
plot(dec.down, show.az=F, col=MESS::col.alpha('blue',0.4), xlim=c(-7,17))
lines(dec.flat$data[[1]]$dec$x, dec.flat$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'A', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Downwards slope', font=2, cex=1.1)

hor <- createHor(c(0,90-7,97,360),c(0,0,15,15), 0.5, c(lat,0,1000), 'Upwards slope', smooth=T, .scale=100)
dec.up <- coordtrans(pdf='uniform', az=90, unc=2, hor=hor, refraction=F)
hor <- createHor(c(0,360),rep(hor2alt(hor,90), 2), 0.5, c(lat,0,1000))
dec.flat <- coordtrans(pdf='uniform', az=90, unc=2, hor=hor, refraction=F)
plot(dec.up, show.az=F, col=MESS::col.alpha('blue',0.4), xlim=c(-7,17))
lines(dec.flat$data[[1]]$dec$x, dec.flat$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'B', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Upwards slope', font=2, cex=1.1)

hor <- createHor(c(0,70,90,110,360),c(0,0,15,0,0), 0.5, c(lat,0,1000), 'Peak', smooth=T, .scale=100)
dec.peak <- coordtrans(pdf='uniform', az=90, unc=2, hor=hor, refraction=F)
hor <- createHor(c(0,360),rep(hor2alt(hor,90),2), 0.5, c(lat,0))
dec.flat <- coordtrans(pdf='uniform', az=90, unc=2, hor=hor, refraction=F)
plot(dec.peak, show.az=F, col=MESS::col.alpha('blue',0.4), xlim=c(-7,17))
lines(dec.flat$data[[1]]$dec$x, dec.flat$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'C', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Peak', font=2, cex=1.1)

hor <- createHor(c(0,70,90,110,360),c(15,15,0,15,15), 0.5, c(lat,0,1000), 'Notch', smooth=T, .scale=100)
dec.notch <- coordtrans(pdf='uniform', az=90, unc=2, hor=hor, refraction=F)
hor <- createHor(c(0,360),rep(hor2alt(hor,90),2), 0.5, c(lat,0,1000))
dec.flat <- coordtrans(pdf='uniform', az=90, unc=2, hor=hor, refraction=F)
plot(dec.notch, show.az=F, col=MESS::col.alpha('blue',0.4), xlim=c(-7,17))
lines(dec.flat$data[[1]]$dec$x, dec.flat$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'D', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Notch', font=2, cex=1.1)

dev.print(png, filename = 'FigS5.png', res=300, width=8, height=8, units='in')



# Figure S6 ----------------------------------------------------------------
par(mfrow=c(2,2))
lat <- 40

hor <- createHor(c(0,70,90,110,360),c(0,0,15,0,0), 0.1, c(lat,0,1000), 'Peak', smooth=T, .scale=100)
dec.0.1 <- coordtrans(pdf='normal', az=90, unc=2, hor=hor, refraction=F)
hor <- createHor(c(0,70,90,110,360),c(0,0,15,0,0), 0.5, c(lat,0,1000), 'Peak', smooth=T, .scale=100)
dec.0.5 <- coordtrans(pdf='normal', az=90, unc=2, hor=hor, refraction=F)
plot(dec.0.1, show.az=F, col=MESS::col.alpha('blue',0.4), xlim=c(-7,17))
lines(dec.0.5$data[[1]]$dec$x, dec.0.5$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'A', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Unc = 0.1º', font=2, cex=1.1)

hor <- createHor(c(0,70,90,110,360),c(0,0,15,0,0), 0.25, c(lat,0,1000), 'Peak', smooth=T, .scale=100)
dec.0.25 <- coordtrans(pdf='normal', az=90, unc=2, hor=hor, refraction=F)
hor <- createHor(c(0,70,90,110,360),c(0,0,15,0,0), 0.5, c(lat,0,1000), 'Peak', smooth=T, .scale=100)
dec.0.5 <- coordtrans(pdf='normal', az=90, unc=2, hor=hor, refraction=F)
plot(dec.0.25, show.az=F, col=MESS::col.alpha('blue',0.4), xlim=c(-7,17))
lines(dec.0.5$data[[1]]$dec$x, dec.0.5$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'B', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Unc = 0.25º', font=2, cex=1.1)

hor <- createHor(c(0,70,90,110,360),c(0,0,15,0,0), 1, c(lat,0,1000), 'Peak', smooth=T, .scale=100)
dec.1 <- coordtrans(pdf='normal', az=90, unc=2, hor=hor, refraction=F)
hor <- createHor(c(0,70,90,110,360),c(0,0,15,0,0), 0.5, c(lat,0,1000), 'Peak', smooth=T, .scale=100)
dec.0.5 <- coordtrans(pdf='normal', az=90, unc=2, hor=hor, refraction=F)
plot(dec.1, show.az=F, col=MESS::col.alpha('blue',0.4), xlim=c(-7,17))
lines(dec.0.5$data[[1]]$dec$x, dec.0.5$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'C', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Unc = 1º', font=2, cex=1.1)

hor <- createHor(c(0,70,90,110,360),c(0,0,15,0,0), 2, c(lat,0,1000), 'Peak', smooth=T, .scale=100)
dec.2 <- coordtrans(pdf='normal', az=90, unc=2, hor=hor, refraction=F)
hor <- createHor(c(0,70,90,110,360),c(0,0,15,0,0), 0.5, c(lat,0,1000), 'Peak', smooth=T, .scale=100)
dec.0.5 <- coordtrans(pdf='normal', az=90, unc=2, hor=hor, refraction=F)
plot(dec.2, show.az=F, col=MESS::col.alpha('blue',0.4), xlim=c(-7,17))
lines(dec.0.5$data[[1]]$dec$x, dec.0.5$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'D', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Unc = 2º', font=2, cex=1.1)

dev.print(png, filename = 'FigS6.png', res=300, width=8, height=8, units='in')



# Figure S7 ----------------------------------------------------------------
par(mfrow=c(2,2))
lat <- 40

hor <- createHor(c(0,70,90,110,360),c(0,0,15,0,0), 0.1, c(lat,0,1000), 'Peak', smooth=T, .scale=100)
dec.0.1 <- coordtrans(pdf='uniform', az=90, unc=2, hor=hor, refraction=F)
hor <- createHor(c(0,70,90,110,360),c(0,0,15,0,0), 0.5, c(lat,0,1000), 'Peak', smooth=T, .scale=100)
dec.0.5 <- coordtrans(pdf='uniform', az=90, unc=2, hor=hor, refraction=F)
plot(dec.0.1, show.az=F, col=MESS::col.alpha('blue',0.4), xlim=c(-7,17))
lines(dec.0.5$data[[1]]$dec$x, dec.0.5$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'A', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Unc = 0.1º', font=2, cex=1.1)

hor <- createHor(c(0,70,90,110,360),c(0,0,15,0,0), 0.25, c(lat,0,1000), 'Peak', smooth=T, .scale=100)
dec.0.25 <- coordtrans(pdf='uniform', az=90, unc=2, hor=hor, refraction=F)
hor <- createHor(c(0,70,90,110,360),c(0,0,15,0,0), 0.5, c(lat,0,1000), 'Peak', smooth=T, .scale=100)
dec.0.5 <- coordtrans(pdf='uniform', az=90, unc=2, hor=hor, refraction=F)
plot(dec.0.25, show.az=F, col=MESS::col.alpha('blue',0.4), xlim=c(-7,17))
lines(dec.0.5$data[[1]]$dec$x, dec.0.5$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'B', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Unc = 0.25º', font=2, cex=1.1)

hor <- createHor(c(0,70,90,110,360),c(0,0,15,0,0), 1, c(lat,0,1000), 'Peak', smooth=T, .scale=100)
dec.1 <- coordtrans(pdf='uniform', az=90, unc=2, hor=hor, refraction=F)
hor <- createHor(c(0,70,90,110,360),c(0,0,15,0,0), 0.5, c(lat,0,1000), 'Peak', smooth=T, .scale=100)
dec.0.5 <- coordtrans(pdf='uniform', az=90, unc=2, hor=hor, refraction=F)
plot(dec.1, show.az=F, col=MESS::col.alpha('blue',0.4), xlim=c(-7,17))
lines(dec.0.5$data[[1]]$dec$x, dec.0.5$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'C', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Unc = 1º', font=2, cex=1.1)

hor <- createHor(c(0,70,90,110,360),c(0,0,15,0,0), 2, c(lat,0,1000), 'Peak', smooth=T, .scale=100)
dec.2 <- coordtrans(pdf='uniform', az=90, unc=2, hor=hor, refraction=F)
hor <- createHor(c(0,70,90,110,360),c(0,0,15,0,0), 0.5, c(lat,0,1000), 'Peak', smooth=T, .scale=100)
dec.0.5 <- coordtrans(pdf='uniform', az=90, unc=2, hor=hor, refraction=F)
plot(dec.2, show.az=F, col=MESS::col.alpha('blue',0.4), xlim=c(-7,17))
lines(dec.0.5$data[[1]]$dec$x, dec.0.5$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'D', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Unc = 2º', font=2, cex=1.1)

dev.print(png, filename = 'FigS7.png', res=300, width=8, height=8, units='in')


# Figure S8 ----------------------------------------------------------------
data <- read.csv('data/Belmonteetal2010_Egypt.csv', stringsAsFactors = F)
ind <- which(substr(data$Temple,1,6)=='Serdab'); data <- data[-c(ind),] ## remove Djoser Serdab left and right eyes
ind <- which(is.na(data$Horizon.Altitude)); data$Horizon.Altitude[ind] <- 0 # data[-ind,] ## remove missing horizon data

horFlat <- list()
for (i in 1:NROW(data)) {
  horFlat[[i]] <- createHor(c(0,180,360), rep(data$Horizon.Altitude[i],3), .5, c(data$Latitude[i], data$Longitude[i],1000), data$Temple)
}

decs <- coordtrans(pdf = 'normal', az=data$Azimuth, unc=2, hor = horFlat, name=data$Temple)
ss <- spd(decs)
plot(ss)

cgram <- read.csv('data/Belmonteetal2010_Egypt_Curvigram.csv')
cgram$Freq <- cgram$Freq/max(cgram$Freq)*max(ss$data$y)
lines(cgram$Dec, cgram$Freq, col='red', lwd=2, lty=3)

families <- read.csv('data/Belmonteetal2010_Egypt_Families.csv', header=T, stringsAsFactors = F)
families$Freq <- families$Freq/max(families$Freq)*max(ss$data$y)
text(families$Dec, families$Freq -.4, labels=families$Family, cex=1.5, vfont=c('serif','bold'))
legend(-60, par('usr')[4], lwd=2, lty=c(1,3), box.lwd=0, bg=NULL, col=c('blue','red'), legend=c('SPD', 'Curvigram'), cex=0.8)

dev.print(png, filename = 'FigS8.png', res=300, width=8, height=6, units='in')
