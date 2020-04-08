#### A probabilistic framework and significance test for the analysis of structural orientations in skyscape archaeology
#### Fabio Silva ::: Bournemouth University ::: 2020 ::: fsilva@bournemouth.ac.uk

###########################################################################
#### Main Text

devtools::install_github('f-silva-archaeo/skyscapeR') ## requires latest version of skyscapeR
source('src.R')


# Figure 1 ----------------------------------------------------------------
xx <- seq(0,180,0.1)
nn <- dnorm(xx, 90,5)
uu <- dunif(xx, 90-10,90+10)

par(mfrow=c(2,1), mar= c(4,2,1,1))
xp <- c(xx, rev(xx)); yp <- c(nn, rep(0,length(xx)))
plot(-999,-999, xlab='Azimuth (º)', axes=F, xlim=c(70,110), ylim=c(0,.1), xaxs='i', yaxs='i')
polygon(xp, yp, col=MESS::col.alpha('blue',0.4), border='blue', lwd=2)
box(); axis(1, at=seq(0,360,1), labels=NA, col='grey')
axis(1, at=seq(0,360,5), labels=NA, col='grey5')
axis(1, at=seq(0,360,10))
text(90, 0.087, labels='90º ± 10º', font=2, cex=1.5)
text(par('usr')[1],0.93*par('usr')[4], 'A', pos=4, font=2, cex=1.7)

xp <- c(xx, rev(xx)); yp <- c(uu, rep(0,length(xx)))
plot(-999,-999, xlab='Azimuth (º)', axes=F, xlim=c(70,110), ylim=c(0,.1), xaxs='i', yaxs='i')
polygon(xp, yp, col=MESS::col.alpha('blue',0.4), border='blue', lwd=2)
box(); axis(1, at=seq(0,360,1), labels=NA, col='grey')
axis(1, at=seq(0,360,5), labels=NA, col='grey5')
axis(1, at=seq(0,360,10))
text(90, 0.087, labels='90º ± 10º', font=2, cex=1.5)
text(par('usr')[1],0.93*par('usr')[4], 'B', pos=4, font=2, cex=1.7)
dev.print(png, filename = 'Fig1.png', res=300, width=8, height=8, units='in')



# Figure 2 ----------------------------------------------------------------
par(mfrow=c(1,1), mar= c(4,1,1,1))
hor <- createHor(c(0,45,80,90,100,135,360), c(0,0,0,6.66,0,0,0), 0.5, c(40,0,1000), smooth=T, .scale=1000)
grid <- sky.objects(seq(-20,20,1),epoch=2000, col='blue', lty=2)
plot(hor, xlim=c(70,110), ylim=c(-3,7), obj=grid, show.az = T)
axis(1, at=seq(0,360,1), labels=NA, col='grey')
axis(1, at=seq(0,360,5), labels=NA, col='grey5')
axis(1, at=seq(0,360,10), labels=NA)

points(90,-2, pch=19); arrows(90-10,-2,90+10,-2, code=3, angle=90, length=0.03)
lines(c(80,80),c(-2,hor2alt(hor,80)+.5), lty=3)
lines(c(100,100),c(-2,hor2alt(hor,100)+.5), lty=3)
lines(c(90,90),c(-2,hor2alt(hor,90)+.5), lty=3)
text(90, -2.4, labels='90º ± 10º', font=2, cex=1.5)

xx <- seq(80,100,0.5)
yy <- hor2alt(hor, xx)
xp <- c(xx, rev(xx))
yp <- c(yy+.5,yy-.5)
polygon(xp, yp, col=MESS::col.alpha('grey',0.4), border='grey', lwd=1.2)
dev.print(png, filename = 'Fig2.png', res=300, width=8, height=4, units='in')
-diff(az2dec(c(80,90,100), hor, refraction=F)) ## declination range on left and right side of 90º



# Figure 3 ----------------------------------------------------------------
az <- 90
unc <- 5
dec.pdf <- coordtrans('normal', az, unc, hor)
plot(dec.pdf)
dev.print(png, filename = 'Fig3.png', res=300, width=6, height=6, units='in')



# Figure 4 ----------------------------------------------------------------
par(mfrow=c(2,2))

lat <- 40
hor <- createHor(c(0,360), c(0,0), 0.5, c(lat,0,1000), 'Flat')
dec3 <- coordtrans(pdf='normal', az=160, unc=5, hor=hor, refraction=F)
plot(dec3, show.az = F, col=MESS::col.alpha('blue',0.4), xlim=az2dec(160, hor, 0, refraction=F)+c(-15,15))
dec1 <- coordtrans(pdf='normal', az=90, unc=5, hor=hor, refraction=F)
lines(dec1$data[[1]]$dec$x + az2dec(160, hor, 0, refraction=F), dec1$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'A', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Az = 160º', font=2, cex=1.5)

dec4 <- coordtrans(pdf='normal', az=180, unc=5, hor=hor, refraction=F)
plot(dec4, show.az = F, col=MESS::col.alpha('blue',0.4), xlim=az2dec(180, hor, 0, refraction=F)+c(-15,15))
dec1 <- coordtrans(pdf='normal', az=90, unc=5, hor=hor, refraction=F)
lines(dec1$data[[1]]$dec$x + az2dec(180, hor, 0, refraction=F), dec1$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'B', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Az = 180º', font=2, cex=1.5)

hor <- createHor(c(0,90-7,97,270,360),c(15,15,0,0,15), 0.5, c(lat,0,1000), 'Downwards slope', smooth=T, .scale=100)
dec.down <- coordtrans(pdf='normal', az=90, unc=2, hor=hor, refraction=F)
hor <- createHor(c(0,360),rep(hor2alt(hor,90),2), 0.5, c(lat,0,1000))
dec.flat <- coordtrans(pdf='normal', az=90, unc=2, hor=hor, refraction=F)
plot(dec.down, show.az=F, col=MESS::col.alpha('blue',0.4), xlim=c(-7,17))
lines(dec.flat$data[[1]]$dec$x, dec.flat$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'C', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Downwards slope', font=2, cex=1.5)

hor <- createHor(c(0,70,90,110,360),c(0,0,15,0,0), 0.5, c(lat,0,1000), 'Peak', smooth=T, .scale=100)
dec.peak <- coordtrans(pdf='normal', az=90, unc=2, hor=hor, refraction=F)
hor <- createHor(c(0,360),rep(hor2alt(hor,90),2), 0.5, c(lat,0))
dec.flat <- coordtrans(pdf='normal', az=90, unc=2, hor=hor, refraction=F)
plot(dec.peak, show.az=F, col=MESS::col.alpha('blue',0.4), xlim=c(-7,17))
lines(dec.flat$data[[1]]$dec$x, dec.flat$data[[1]]$dec$y, lwd=2, lty=3, col='red')
text(par('usr')[1],0.93*par('usr')[4], 'D', pos=4, font=2, cex=1.7)
text(mean(c(par('usr')[1],par('usr')[2])),0.93*par('usr')[4], 'Peak', font=2, cex=1.5)

dev.print(png, filename = 'Fig4.png', res=300, width=8, height=8, units='in')



# Figure 5 ----------------------------------------------------------------
hor.Peak <- createHor(c(0,85,90,95,360),c(0,0,5,0,0), 0.5, c(40,0,1000), 'Peak', smooth=T, .scale=300)
# az <- rnorm(5, 90, 5)
az <- c(87.95794,92.26147,82.47100,93.85006,92.86030) ## this was the result of the random pick used in the paper

decs <- coordtrans(pdf='normal', az=az, unc=2, hor=hor.Peak, refraction=F)
ss <- spd(decs)
par(mfrow=c(1,1))
plot(ss)

discrete.decs <- az2dec(az, hor.Peak, hor2alt(hor.Peak,az), refraction=F)
cgram <- density(discrete.decs, bw=az2dec(90,hor.Peak, F) - az2dec(92,hor.Peak, F))
cgram$y <- cgram$y/max(cgram$y)*max(ss$data$y)
lines(cgram$x, cgram$y, col='#7b3294', lwd=2, lty=2)

kde <- density(discrete.decs, bw='nrd')
kde$y <- kde$y/max(kde$y)*max(ss$data$y)
lines(kde$x, kde$y, col='#008837', lwd=2, lty=3)

legend(5.7, par('usr')[4], lwd=1.5, box.lwd=0, bg=NULL, col=c('blue','#7b3294','#008837'), lty=1:3, legend=c('SPD (new method)', 'Curvigram', 'KDE'), cex=0.8)

dev.print(png, filename = 'Fig5.png', res=300, width=6, height=6, units='in')

ss$data$x[which.max(ss$data$y)]
cgram$x[which.max(cgram$y)]



# Figure 6 ----------------------------------------------------------------
xx <- seq(-45,360+45,0.1)
uu <- dunif(xx, 0, 360)

hor.flat <- createHor(c(0,360), c(0,0), 0.5, c(40,0,1000))
null <- coordtrans('uniform', 180, 180, hor.flat, refraction=F)

par(mfrow=c(2,1), mar= c(4,1,1,1))
xp <- c(xx, rev(xx)); yp <- c(uu, rep(0,length(xx)))
plot(-999,-999, xlab='Azimuth (º)', axes=F, xlim=c(-20,360+20), ylim=c(0,.005), xaxs='i', yaxs='i')
polygon(xp, yp, col=MESS::col.alpha('blue',0.4), border='blue', lwd=2)
axis(1, at=seq(0,360,10), labels=NA, col.ticks='grey5', lwd.ticks = 0.5)
axis(1, at=seq(0,360,45)); box()
text(par('usr')[1], par('usr')[4]*0.92, labels='A', font=2, pos=4, cex=1.7)

xp <- c(null$data[[1]]$dec$x, rev(null$data[[1]]$dec$x)); yp <- c(null$data[[1]]$dec$y, rep(0,length(null$data[[1]]$dec$y)))
plot(-999,-999, xlab='Declination (º)', axes=F, xlim=az2dec(0, hor.flat, 0, refraction=F)*c(-1,1)+c(-3.896,3.896), ylim=c(0,.05), xaxs='i', yaxs='i')
polygon(xp, yp, col=MESS::col.alpha('blue',0.4), border='blue', lwd=2)
axis(1, at=seq(-60,60,5), labels=NA, col.ticks='grey5', lwd.ticks = 0.5)
axis(1, at=seq(-60,60,10)); box()
text(par('usr')[1], par('usr')[4]*0.92, labels='B', font=2, pos=4, cex=1.7)

dev.print(png, filename = 'Fig6.png', res=300, width=8, height=8, units='in')



# Figure 7 ----------------------------------------------------------------
# sg <- sigTestDec(decs, 10, ncores=1, save.sim=T) # ten random picks
# sg100 <- sigTestDec(decs, 100) # hundred random picks
# sg1k <- sigTestDec(decs, 1000) # thousand random picks
# sg5k <- sigTestDec(decs, nsims=5000)
# save(list=c('sg', 'sg100', 'sg1k','sg5k'), file='data/method_sigTest.RData')
load('data/method_sigTest.RData')

par(mfrow=c(5,1))
# empirical dataset
ss <- spd(decs)
plot(ss, xlim = c(-60,60), ylim=c(0,2))
text(par('usr')[1], par('usr')[4]*0.85, labels='A: empirical SPD', font=2, pos=4, cex=1.7)

# single random pick
decs1 <- coordtrans('normal', sg$data$simulated$Az[1,], sg$data$simulated$Unc[1,], hor.Peak, refraction=T)
ss1 <- spd(decs1)
plot(ss, xlim = c(-60,60), ylim=c(0,2))
lines(ss1$data$x, ss1$data$y, col='grey', lwd=2)
text(par('usr')[1], par('usr')[4]*0.85, labels='B: single simulated SPD', font=2, pos=4, cex=1.7)

# ten random picks
plot(ss, xlim = c(-60,60), ylim=c(0,2))
lines(ss1$data$x, ss1$data$y, col='grey', lwd=2)
for (i in 2:10) {
  decs1 <- coordtrans('normal', sg$data$simulated$Az[i,], sg$data$simulated$Unc[i,], hor.Peak, refraction=T, verbose=F)
  ss1 <- spd(decs1)
  lines(ss1$data$x, ss1$data$y, col='grey', lwd=2)
}
text(par('usr')[1], par('usr')[4]*0.85, labels='C: ten simulated SPDs', font=2, pos=4, cex=1.7)

# hundred random picks
plot(ss, xlim = c(-60,60), ylim=c(0,2))
xp <- c(sg100$data$null.hyp$x, rev(sg100$data$null.hyp$x))
yp <- c(sg100$data$null.hyp$CE.upper, sg100$data$null.hyp$CE.lower)
polygon(xp, yp, border='grey', col=MESS::col.alpha('grey',.5))
text(par('usr')[1], par('usr')[4]*0.85, labels='D: Confidence Envelope from 100 simulated SPDs', font=2, pos=4, cex=1.7)

# five thousand and p-values
plot(ss, xlim = c(-60,60), ylim=c(-0.3,2))
xp <- c(sg5k$data$null.hyp$x, rev(sg5k$data$null.hyp$x))
yp <- c(sg5k$data$null.hyp$CE.upper, sg5k$data$null.hyp$CE.lower)
polygon(xp, yp, border='grey', col=MESS::col.alpha('grey',.5))
text(par('usr')[1], par('usr')[4]*0.85, labels='E: Confidence Envelope and p-values from 5000 simulated SPDs', font=2, pos=4, cex=1.7)
abline(h=0)
aux <- as.matrix(sg5k$metadata$local.pval[sg5k$metadata$local.pval$type=='+',2:4])
for (i in 1:NROW(aux)) {
  if (sg5k$metadata$local.pval[i,1] == '+') { col <- pal[2] } else { col <- pal[1] }
  xp <- c(aux[i,1], aux[i,1], aux[i,2], aux[i,2])
  yp <- c(-0.3, 0, 0, -0.3)
  polygon(xp, yp, col=MESS::col.alpha(col, .4), border=NA)
  text(mean(aux[i,1:2]), -0.3/2, labels=pval2stars(aux[i,3]), cex=1.5)
}
text(par('usr')[2], par('usr')[4]*0.75, labels=paste0('global p-value = ', sg5k$metadata$global.pval), font=2, pos=2, cex=1.2)
dev.print(png, filename = 'Fig7.png', res=300, width=8, height=10, units='in')
par(mfrow=c(1,1))



# Figure 8 ----------------------------------------------------------------
data <- read.csv('data/Belmonteetal2010_Egypt.csv', stringsAsFactors = F)
ind <- which(substr(data$Temple,1,6)=='Serdab'); data <- data[-c(ind),] ## remove Djoser Serdab left and right eyes
ind <- which(is.na(data$Horizon.Altitude)); data$Horizon.Altitude[ind] <- 0 # data[-ind,] ## remove missing horizon data

horFlat <- list()
for (i in 1:NROW(data)) {
  horFlat[[i]] <- createHor(c(0,180,360), rep(data$Horizon.Altitude[i],3), .5, c(data$Latitude[i], data$Longitude[i],1000), data$Temple)
}

decs <- coordtrans(pdf = 'normal', az=data$Azimuth, unc=2, hor = horFlat, name=data$Temple)
ss <- spd(decs)

# sg1 <- sigTestDec(decs, nsims=5000, ncores=11)
# save(sg1, file='data/Belmonte_sigTest.RData')
load('data/Belmonte_sigTest.RData')

par(mar=c(5,4,1,1))
plot(sg1, show.local=T)
lines(rep(eq(max(data$Latitude)),2), c(0,12), lty=1, lwd=1.5)
lines(rep(dS(-1500, max(data$Latitude)),2), c(0,12), lty=2, lwd=1.5)
lines(rep(-10,2), c(0,12), lty=3, lwd=1.5)
legend('topleft', lwd=1.5, box.lwd=0, bg=NULL, lty=1:3, legend=c('Equinox/East', 'December solstice', 'Peret/Shomu start c. 1500 BC'), cex=0.8)
dev.print(png, filename = 'Fig8_rraw.png', res=300, width=8, height=6, units='in')



# Table 1 ---------------------------------------------------------------
sg1$metadata$local.pval[sg1$metadata$local.pval$type == '+',]



# Figure 9 ---------------------------------------------------------------
ind.OK <- which(data$Dynasty == '3rd' | data$Dynasty == '4th' | data$Dynasty == '4th ?' | data$Dynasty == '5th' | data$Dynasty == '6th')
decs.OK <- coordtrans(pdf = 'normal', az=data$Azimuth[ind.OK], unc=2, hor = horFlat, name=data$Temple[ind.OK])
ss.OK <- spd(decs.OK)

ind.MK <- which(data$Dynasty == '11th' | data$Dynasty == '12th')
decs.MK <- coordtrans(pdf = 'normal', az=data$Azimuth[ind.MK], unc=2, hor = horFlat, name=data$Temple[ind.MK])
ss.MK <- spd(decs.MK)

ind.NK <- which(data$Dynasty == '18th' | data$Dynasty == '18th-19th' | data$Dynasty == '19th')
decs.NK <- coordtrans(pdf = 'normal', az=data$Azimuth[ind.NK], unc=2, hor = horFlat, name=data$Temple[ind.NK])
ss.NK <- spd(decs.NK)

ind.Ptolemaic <- which(data$Dynasty == 'Ptolemaic' | data$Dynasty == 'Ptolemaic ?' | data$Dynasty == 'Ptolemiac')
decs.Ptolemaic <- coordtrans(pdf = 'normal', az=data$Azimuth[ind.Ptolemaic], unc=2, hor = horFlat, name=data$Temple[ind.Ptolemaic])
ss.Ptolemaic <- spd(decs.Ptolemaic)

# sg.OK <- sigTestDec(decs.OK, nsims=5000, ncores=19)
# sg.MK <- sigTestDec(decs.MK, nsims=5000, ncores=19)
# sg.NK <- sigTestDec(decs.NK, nsims=5000, ncores=19)
# sg.Ptolemaic <- sigTestDec(decs.Ptolemaic, nsims=5000, ncores=19)
# save(list=c('sg.OK', 'sg.MK', 'sg.NK','sg.Ptolemaic'), file='data/Belmonte_chrono_sigTest.RData')
load('data/Belmonte_chrono_sigTest.RData')

par(mfrow=c(2,2), mar=c(4,4,1,1))
plot(sg.OK, xlim=c(-62,62), show.local=T, show.pval = F)
plot(sg.MK, xlim=c(-62,62), show.local=T, show.pval = F)
plot(sg.NK, xlim=c(-62,62), show.local=T, show.pval = F)
plot(sg.Ptolemaic, xlim=c(-62,62), show.local=T, show.pval = F)
dev.print(png, filename = 'Fig9_raw.png', res=300, width=11, height=8, units='in')



# Table 2 ---------------------------------------------------------------
sg.OK$metadata$local.pval[sg.OK$metadata$local.pval$type == '+',]
sg.MK$metadata$local.pval[sg.MK$metadata$local.pval$type == '+',]
sg.NK$metadata$local.pval[sg.NK$metadata$local.pval$type == '+',]
sg.Ptolemaic$metadata$local.pval[sg.Ptolemaic$metadata$local.pval$type == '+',]



# Figure 10 ---------------------------------------------------------------
data <- read.csv('data/Ruggles1999_RSC.csv', stringsAsFactors = F)

half.width <- function(x) {
  if (sum(is.na(x)) < 4) {
    aux <- c(sum(range(x, na.rm=T))/2, abs(diff(range(x, na.rm=T)))/2)
  } else { aux <- c(NA,NA) }
  return(aux)
}

ind <- which(data$Name=='Sunhoney')
hw <- half.width(data[ind,7:11])

xx <- seq(0, 360, 0.1)

modA <- dnorm(xx, data$CL_Rec_C[ind], .5)
modB <- dnorm(xx, hw[1], hw[2]/2)
modC <- dunif(xx, hw[1]-hw[2], hw[1]+hw[2])

par(mfrow=c(3,1), mar=c(4,1,1,1))
xp <- c(xx, rev(xx)); yp <- c(modA, rep(0,length(xx)))
plot(-999,-999, xlab='', axes=F, xlim=c(213,253), ylim=c(0,1), xaxs='i', yaxs='i')
polygon(xp, yp, col=MESS::col.alpha('blue',0.4), border='blue', lwd=2)
box()
text(par('usr')[1],0.9*par('usr')[4], 'Model A', pos=4, font=2, cex=1.7)

xp <- c(xx, rev(xx)); yp <- c(modB, rep(0,length(xx)))
plot(-999,-999, xlab='', axes=F, xlim=c(213,253), ylim=c(0,.1), xaxs='i', yaxs='i')
polygon(xp, yp, col=MESS::col.alpha('blue',0.4), border='blue', lwd=2)
box()
text(par('usr')[1],0.9*par('usr')[4], 'Model B', pos=4, font=2, cex=1.7)

xp <- c(xx, rev(xx)); yp <- c(modC, rep(0,length(xx)))
plot(-999,-999, xlab='Azimuth (º)', axes=F, xlim=c(213,253), ylim=c(0,.07), xaxs='i', yaxs='i')
polygon(xp, yp, col=MESS::col.alpha('blue',0.4), border='blue', lwd=2)
box(); axis(1, at=seq(0,360,1), labels=NA, col='grey')
axis(1, at=seq(0,360,5), labels=NA, col='grey5')
axis(1, at=seq(0,360,10))
text(par('usr')[1],0.9*par('usr')[4], 'Model C', pos=4, font=2, cex=1.7)

dev.print(png, filename = 'Fig10_raw.png', res=300, width=6, height=8, units='in')



# Figure 11 ---------------------------------------------------------------
load('data/Ruggles1999_RSC_HWT.RData')
decs.modA <- coordtrans('normal', az = data$CL_Rec_C, unc = .5, hor = horHWT, name=data$Name)
spd.modA <- spd(decs.modA)

hw <- matrix(NA,ncol=2,nrow=NROW(data))
for (i in 1:NROW(data)) {
  hw[i,] <- half.width(data[i,7:11])
}

decs.modB <- coordtrans(pdf='normal', az=hw[,1], unc=hw[,2]/2, hor = horHWT, name=data$Name)
spd.modB <- spd(decs.modB)

decs.modC <- coordtrans(pdf = 'uniform', az=hw[,1], unc=hw[,2], hor = horHWT, name=data$Name)
spd.modC <- spd(decs.modC)

par(mfrow=c(1,1))
plot(spd.modC, ylim = c(0,7.5))
lines(spd.modB$data, col='#7b3294', lwd=2, lty=2)
lines(spd.modA$data, col='#008837', lwd=2, lty=3)
legend(-15, par('usr')[4], lwd=2, lty=c(3,2,1), box.lwd=0, bg=NULL, col=c('#008837','#7b3294','blue'), legend=c('Model A', 'Model B', 'Model C'), cex=1)
dev.print(png, filename = 'Fig11.png', res=300, width=8, height=6, units='in')



# Figure 12---------------------------------------------------------------
# sg.modA <- sigTestDec(decs.modA, nsims=5000, ncores=11)
# sg.modB <- sigTestDec(decs.modB, nsims=5000, ncores=11)
# sg.modC <- sigTestDec(decs.modC, nsims=5000, ncores=11)
# save(list=c('sg.modA','sg.modB','sg.modC'), file='data/Ruggles_sigTest.RData')
load('data/Ruggles_sigTest.RData')

par(mfrow=c(1,3), mar=c(4,4,2,1))
plot(sg.modA, show.pval=T, show.local=T, xlim=c(-35,-5))
mtext('Model A', side=3, font=2, cex=1.5)
plot(sg.modB, show.pval=T, show.local=T, xlim=c(-35,-5))
mtext('Model B', side=3, font=2, cex=1.5)
plot(sg.modC, show.pval=T, show.local=T, xlim=c(-35,-5))
mtext('Model C', side=3, font=2, cex=1.5)
dev.print(png, filename = 'Fig12.png', res=300, width=14, height=4, units='in')



# Table 3---------------------------------------------------------------
sg.modA$metadata$local.pval[sg.modA$metadata$local.pval$type == '+',]
sg.modB$metadata$local.pval[sg.modB$metadata$local.pval$type == '+',]
sg.modC$metadata$local.pval[sg.modC$metadata$local.pval$type == '+',]



# Figure 13 ---------------------------------------------------------------
par(mfrow=c(1,1), mar=c(4,1,1,1))
plot(-999,-999, xlim=c(-40,-15), ylim=c(0,9.99), xlab='Declination (º)', ylab='', axes=F, xaxs='i', yaxs='i'); box()
axis(1, at=pretty(seq(par('usr')[1],par('usr')[2])))
axis(1, at=seq(-90,90,1), lwd=0.5, labels=F)

text(-38, 8.325, labels='Model A', font=2, cex=1.5)
aux <- sg.modA$metadata$local.pval
ind <- which(aux$type=='+')
for (i in ind) {
  polygon(c(aux$startDec[i], aux$endDec[i], aux$endDec[i], aux$startDec[i]), c(6.66,6.66,9.99,9.99), col=MESS::col.alpha('#008837',0.4), border=NA)
  text(mean(as.numeric(aux[i,2:3])), 8.325, labels=pval2stars(aux[i,4]), cex=1.5)
}
lines(spd.modA$data$x, spd.modA$data$y/max(spd.modA$data$y)*3+6.66 , lwd=1.5, lty=1)

text(-38, 4.995, labels='Model B', font=2, cex=1.5)
aux <- sg.modB$metadata$local.pval
ind <- which(aux$type=='+')
for (i in ind) {
  polygon(c(aux$startDec[i], aux$endDec[i], aux$endDec[i], aux$startDec[i]), c(3.33,3.33,6.66,6.66), col=MESS::col.alpha('#008837',0.4), border=NA)
  text(mean(as.numeric(aux[i,2:3])), 4.995, labels=pval2stars(aux[i,4]), cex=1.5)
}
lines(spd.modB$data$x, spd.modB$data$y/max(spd.modB$data$y)*3+3.33 , lwd=1.5, lty=1)

text(-38, 1.665, labels='Model C', font=2, cex=1.5)
aux <- sg.modC$metadata$local.pval
ind <- which(aux$type=='+')
for (i in ind) {
  polygon(c(aux$startDec[i], aux$endDec[i], aux$endDec[i], aux$startDec[i]), c(0,0,3.33,3.33), col=MESS::col.alpha('#008837',0.4), border=NA)
  text(mean(as.numeric(aux[i,2:3])), 1.665, labels=pval2stars(aux[i,4]), cex=1.5)
}
lines(spd.modC$data$x, spd.modC$data$y/max(spd.modC$data$y)*3 , lwd=1.5, lty=1)

abline(h=c(3.33,6.66))

decs <- c(sMjLX(-2500, max(data$Latitude)), smnLX(-1750, min(data$Latitude)))
polygon(c(decs, rev(decs)), c(-5,-5,20,20), col=MESS::col.alpha('orange',0.3), border=NA)
decs <- c(dS(-2600), dS(-1750))
lines(rep(mean(decs),2), c(-5,10), lwd=2, lty=2, col='black')
dev.print(png, filename = 'Fig13.png', res=300, width=8, height=4, units='in')


