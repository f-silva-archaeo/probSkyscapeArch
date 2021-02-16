library(skyscapeR)
library(foreach)
library(RColorBrewer)
pal <- brewer.pal(5, 'PRGn')[c(1,5)]


# Updates to skyscapeR library ----------------------------------------
eq2hor_new <- function(ra, dec, loc, refraction=T, atm=1013.25, temp=15, out='both') {
  jd <- swephR::swe_julday(2000,1,1,12,1)
  xx <- swephR::swe_azalt(jd, 1, c(loc[2],loc[1],loc[3]), atm, temp, c(ra,dec))$xaz
  
  if (refraction) { alt <- xx[3] } else { alt <- xx[2] }
  az <- xx[1]-180
  if (az < 0) { az <- az + 360 }
  if (az > 360) { az <- az - 360 }
  
  if (out=='both') { aux <- c(az, alt) }
  if (out=='alt') { aux <- alt }
  if (out=='az') { aux <- az }
  return(aux)
}

veq2hor <- Vectorize(eq2hor_new, 'ra', SIMPLIFY=T)

az2dec <- function(az, loc, alt, refraction=T, atm=1013.25, temp=15){
  jd <- swephR::swe_julday(2000,1,1,12,1)
  prec <- max(nchar(sub('.*\\.', '', as.character(az))))
  
  if (missing(alt)) {
    if (class(loc) == 'skyscapeR.horizon') { 
      alt <- hor2alt(loc, az) 
    } else if (class(loc) == 'list') { 
      alt <- c()
      for (i in 1:NROW(loc)) {
        alt[i] <- hor2alt(loc[[i]], az)
      }
    } else { stop('Altitude missing.') } }
  
  if (length(alt) == 1) { alt <- rep(alt, NROW(az)) }
  
  if (class(loc) == 'skyscapeR.horizon') { georefs <- loc$metadata$georef; georefs <- matrix(georefs, ncol=3, nrow=NROW(az), byrow=T) }
  if (class(loc) == 'list') { 
    georefs <- matrix(NA, ncol=3, nrow=NROW(loc))
    for (i in 1:NROW(loc)) { georefs[i,] <- loc[[i]]$metadata$georef } }
  if (class(loc) == 'numeric') {
    georefs <- matrix(NA, ncol=3, nrow=NROW(loc))
    for (i in 1:NROW(loc)) { georefs[i,] <- c(loc[i],0,0) }
  }
  if (class(loc) == 'matrix') {
    georefs <- matrix(NA, ncol=3, nrow=NROW(loc))
    for (i in 1:NROW(loc)) {
      if (dim(loc)[2]==2) { georefs[i,] <- c(loc[i,1], loc[i,2], 0) }
      if (dim(loc)[2]==3) { georefs[i,] <- loc[i,] }
    }
  }
  
  
  dec <- c()
  for (i in 1:NROW(az)) {
    georef <- georefs[i,]
    if (refraction) {
      alt[i] <- alt[i] - (swephR::swe_refrac_extended(alt[i], georef[3], atm, temp, 0, 0)$return - alt[i])
    }
    dec[i] <- round( swephR::swe_azalt_rev(jd, 1, c(georef[2],georef[1],georef[3]), c(az[i]-180, alt[i]))$xout[2], prec)
  }
  
  return(dec)
}


createHor <- function(az, alt, alt.unc=0.5, loc, name='', smooth=F, .scale=1000) {
  # return result
  hor <- c()
  hor$metadata$name <- name
  if (length(loc)<3) { cat('No site elevation given. Assuming 0 metres above sea level.'); loc <- c(loc,0) }
  if (length(loc)<2) { stop('Site latitude, longitude and elevation needed.') }
  hor$metadata$georef <- loc; names(hor$metadata$georef) <- c('Lat','Lon','Elev'); dim(hor$metadata$georef) <- c(1,3)
  if (length(alt.unc) > 1 & length(alt.unc) < length(az)) { stop('Altitude uncertainty must be either a single value or have the same length as "az"') }
  if (length(alt.unc) == 1 ) { alt.unc <- rep(alt.unc, length(az)) }
  
  az <- c(az-360, az, az+360); alt <- rep(alt, 3); alt.unc <- rep(alt.unc, 3)
  ind <- which(duplicated(az)); if (length(ind)>0) { az <- az[-ind]; alt <- alt[-ind]; alt.unc <- alt.unc[-ind] }
  xx <- seq(0-90, 360+90, 0.01)
  alt <- approx(az, alt, xx)$y
  alt.unc <- approx(az, alt.unc, xx)$y
  hor$data <- data.frame(az = xx, alt = alt, alt.unc = alt.unc)
  
  if (smooth) {
    yy <- zoo::rollmean(alt, .scale, fill=0)
    hor$data$alt <- yy
    
    yy <- zoo::rollmean(alt.unc, .scale, fill=0)
    hor$data$alt.unc <- yy
  }
  class(hor) <- "skyscapeR.horizon"
  return(hor)
}


hor2alt <- function (hor, az) {
  hh <- splinefun(hor$data$az, hor$data$alt)
  alt <- round(hh(az), 2)
  return(alt)
}
#############################################################################################################


## replaces R function (twice as fast)
dnorm <- function(x, mean, sd) { return(exp(-(x-mean)^2/(2*sd^2)) / sqrt(2*pi*sd^2)) }


# Coordinate-Transformation of Azimuthal Distribution ---------------------
coordtrans_unc <- function(pdf='normal', az, unc, hor, name='', refraction=T, atm=1013.25, temp=15, verbose=T, .cutoff=1e-4, .prec_az=0.01, .prec_dec=0.1) {
  
  ## probability distribution
  if (class(pdf) == 'character') {
    
    code <- switch(pdf, N = 0, Normal = 0, n = 0, normal = 0, U = 1, Uniform = 1, u = 1, uniform = 1, -1)
    
    if (code == 0) {
      if (verbose) { cat('Normal probability distribution(s) chosen\n') }
      az.pdf <- function(x,i) dnorm(x, az[i], unc[i])
    } else 
      
      if (code == 1) {
        if (verbose) { cat('Uniform probability distribution(s) chosen\n') }
        az.pdf <- function(x,i) dunif(x, az[i]-unc[i], az[i]+unc[i])
      } else { stop('Probability density function not recognized. Please use only normal or uniform.')}
    
    if (missing(az) | missing (unc)) {
      stop('Missing azimuths and/or uncertainties.')
    }
    
  } else if (class(pdf) == 'data.frame' & sum(colnames(pdf) == c('az','dens'))) {
    if (verbose) { stop('Custom probability distribution found. This feature is not yet functional.') }
    
  } else if (class(pdf) == 'list' & class(pdf[[1]]) == 'data.frame' & sum(colnames(pdf[[1]]) == c('az','dens'))) { 
    if (verbose) { stop('Custom probability distribution(s) found. This feature is not yet functional.') }
    
  } else {
    stop('Probability density function not recognized. Please use only normal or uniform.')# or a list of data.frames with columns "az" and "dens".') 
  }
  
  ## init parameters
  n <- length(az)
  if (n > 1 & length(unc) == 1) { 
    if (verbose) { cat('Single uncertainty value found. Using it for all measurements.\n') }
    unc <- rep(unc, length(az)) 
  }
  
  if (class(hor) == "skyscapeR.horizon") { 
    if (verbose) { cat('Single horizon profile found. Using it for all measurements.\n') }
    hor <- rep(list(hor),n)
  }
  
  allaz <- seq(0-90, 360+90, .prec_az)
  xx <- seq(-90, 360+90, .1)
  
  out <- list()
  if (verbose & length(az)>1) { pb <- txtProgressBar(max = length(az), style = 3) }
  
  for (k in 1:n) {
    ind <- which(az.pdf(allaz,k) >= .cutoff); azs <- allaz[ind]
    
    hh <- hor[[k]]$data
    loc <- hor[[k]]$metadata$georef
    fhor <- approx(hh$az, hh$alt, xout=azs)$y
    fhor.unc <- approx(hh$az, hh$alt.unc, xout=azs)$y
    
    aux1 <- az2dec(azs, hor[[k]], fhor-4*fhor.unc, refraction=refraction, atm=atm, temp=temp)
    aux2 <- az2dec(azs, hor[[k]], fhor+4*fhor.unc, refraction=refraction, atm=atm, temp=temp)
    
    min.dec <- min(aux1, aux2, na.rm = T); max.dec <- max(aux1, aux2, na.rm = T)
    
    dec <- seq(min.dec, max.dec, by = .prec_dec)
    mean.dec <- (aux2+aux1)/2
    sd.dec <- (aux2-aux1)/8
    
    dens <- rep(NA, length(dec))
    if (verbose & n==1) { pb <- txtProgressBar(max = length(dec), style = 3) }
    for (i in 1:length(dec)) {
      d <- dnorm(dec[i], mean.dec, sd.dec)  
      f <- approxfun(azs, d, yleft=0, yright=0)
      
      # do convolution manually
      conv <- f(xx)*(az.pdf(xx,k))
      dens[i] <- MESS::auc(xx, conv)
      
      if (verbose & n==1) { setTxtProgressBar(pb, i) }
    }
    
    ss <- splinefun(dec, dens)
    dec <- seq(min.dec, max.dec, .001)
    dens <- ss(dec)
    
    # normalise
    dens <- dens/MESS::auc(dec,dens)
    
    aux <- list(
      az = data.frame(x=azs, y=az.pdf(azs,k)),
      dec = data.frame(x=dec, y=dens) )
    
    out[[k]] <- aux
    if (verbose & n>1) { setTxtProgressBar(pb, k) }
  }
  
  mtdta <- list(
    name = name,
    pdf = pdf,
    az = az,
    unc = unc)
  mtdta$param <- list(refraction = refraction, atm=atm, temp=temp, .cutoff=.cutoff, .prec_az=.prec_az, .prec_dec=.prec_dec)
  mtdta$horizon = hor
  if (n > 1) { dta <- out } else { dta <- list(aux) }
  out <- list(metadata = mtdta, data = dta)
  class(out) <- 'skyscapeR.pdf'
  
  if (verbose) { cat('\nDone.') }
  return(out)
}
coordtrans <- compiler::cmpfun(coordtrans_unc, options=list(enableJIT = 3))





spd <- function(pdf, normalise = F, xrange, .cutoff = 1e-5, .prec_dec = 0.01) {
  if (class(pdf)!='skyscapeR.pdf') { stop('Please provide a valid skyscapeR.pdf object.') }
  
  dd <- pdf$data; aux <- lapply(dd, "[[", 'dec')
  x <- lapply(aux, "[[", 'x')
  y <- lapply(aux, "[[", 'y')
  
  if (missing(xrange)) { xrange <- range(x)}
  
  decs <- seq(xrange[1], xrange[2], .prec_dec)
  spd <- rep(0, length(decs))
  for (i in 1:NROW(x)) {
    spd <- spd + approx(x[[i]], y[[i]], xout=decs, yleft=0, yright=0)$y
  }
  
  spd[spd < .cutoff] <- 0
  
  if (normalise) { spd <- spd/sum(spd, na.rm = T) }
  out <- c()
  out$metadata <- c()
  out$metadata$names <- pdf$metadata$name
  out$metadata$measurements <- pdf$metadata$az
  out$metadata$uncertainty <- pdf$metadata$unc
  out$metadata$xrange <- xrange
  out$metadata$.cutoff <- .cutoff
  out$metadata$.prec_dec <- .prec_dec
  out$data <- data.frame(x = decs, y = spd)
  
  class(out) <- 'skyscapeR.spd'
  return(out)
}



# Declination Significance Test -------------------------------------------
sigTestDec_unc <- function(decs, nsims=1000, conf=.95, tails=2, normalise=F, ncores=parallel::detectCores()-1, save.sim=F, verbose=T) {
  
  pdf <- decs$metadata$pdf
  az <- decs$metadata$az
  unc <- decs$metadata$unc
  param <- decs$metadata$param
  hor <- decs$metadata$horizon
  
  if (save.sim & ncores>1) { stop('Cannot save simulated values when running on more than one core.')}
  
  ## empirical SPD
  if (verbose) { cat('Creating Empirical SPD...') }
  empirical <- spd(decs, xrange=c(-90,90), normalise = normalise)
  if (verbose) { cat('Done.\n') }
  
  
  if (ncores > 1) {
    ## bootstrapping
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    parallel::clusterEvalQ(cl, source("src.R"))
    if (verbose) { cat(paste0('Running ', nsims,' simulations on ', ncores, ' processing cores. This may take a while...')) }
    
    res <- matrix(NA, nsims, length(empirical$data$y))
    res <- foreach (i = 1:nsims, .combine=rbind, .inorder = F) %dopar% {
      simAz <- sample(seq(0, 360, param$.prec_az), length(az), replace=T)
      simUnc <- sample(unc, replace=F)
      cot <- coordtrans(pdf, simAz, simUnc, hor, refraction=param$refraction, atm=param$atm, temp=param$temp, verbose=F, .cutoff=param$.cutoff, .prec_az=param$.prec_az, .prec_dec=param$.prec_dec)
      spd(cot, xrange=empirical$metadata$xrange, normalise = normalise)$data$y
      
    }
    parallel::stopCluster(cl)
    if (verbose) { cat('Done.\n') }
    
  } else {
    if (verbose) { cat(paste0('Running ', nsims,' simulations on a single processing core. This may take a while...')) }
    if (save.sim) { 
      sim.az <- matrix(NA, nsims, length(az))
      sim.unc <- matrix(NA, nsims, length(az))
    }
    res <- matrix(NA, nsims, length(empirical$data$y))
    if (verbose) { pb <- txtProgressBar(max=nsims, style=3) }
    for (i in 1:nsims) {
      simAz <- sample(seq(0, 360, param$.prec_az), length(az), replace=T)
      simUnc <- sample(unc, replace=F)
      cot <- coordtrans(pdf, simAz, simUnc, hor, refraction=param$refraction, atm=param$atm, temp=param$temp, verbose=F, .cutoff=param$.cutoff, .prec_az=param$.prec_az, .prec_dec=param$.prec_dec)
      res[i,] <- spd(cot, xrange=empirical$metadata$xrange, normalise = normalise)$data$y
      
      if (save.sim) { sim.az[i,] <- simAz; sim.unc[i,] <- simUnc }
      if (verbose) { setTxtProgressBar(pb, i) }
    }
    if (verbose) { cat('Done.\n') }
  }
  
  ## confidence envelope
  zScore.sim <- matrix(0, nrow=nsims, ncol=length(empirical$data$y))
  zMean <- colMeans(res)
  zStd <- apply(res, 2, sd)
  
  zScore.emp <- (empirical$data$y - zMean)/zStd
  zScore.sim <- apply(res, 2, scale)
  
  if (verbose) { cat(paste0('Performing a ',tails,'-tailed test at the ', conf*100, '% significance level.\n')) }
  if (tails==2) {
    lvl.up <- 1-(1-conf)/2; lvl.dn <- (1-conf)/2
  } else if (tails==1) {
    lvl.up <- conf; lvl.dn <- 0
  } else { stop() }
  
  upper <- apply(zScore.sim, 2, quantile, probs = lvl.up, na.rm = T)
  upCI <- zMean + upper*zStd
  upCI[is.na(upCI)] <- 0
  if (tails==2) {
    lower <- apply(zScore.sim, 2, quantile, probs = lvl.dn, na.rm = T)
    loCI <- zMean + lower*zStd
    loCI[is.na(loCI)] <- 0
  }
  
  ## global p-value
  above <- which(zScore.emp > upper); emp.stat <- sum(zScore.emp[above] - upper[above])
  if (tails==2) { below <- which(zScore.emp < lower); emp.stat <- emp.stat + sum(lower[below] - zScore.emp[below]) }
  
  sim.stat <- abs(apply(zScore.sim, 1, function(x,y) { a = x-y; i = which(a>0); return(sum(a[i])) }, y = upper))
  if (tails==2) { sim.stat <- sim.stat + abs(apply(zScore.sim, 1, function(x,y) { a = y-x; i = which(a>0); return(sum(a[i])) }, y = lower)) }
  global.p <- round( 1 - ( length(sim.stat[sim.stat < emp.stat]) + 1 ) / ( nsims + 1 ), 5)
  
  ## local p-value
  ind <- split(above, cumsum(c(1,diff(above) > 1))); ind <- ind[which(lengths(ind) > 1)]
  local <- data.frame(type=NA, startDec=0, endDec=0, p.value=0); j <- 0
  if (length(ind)>0) {
    for (j in 1:NROW(ind)) {
      emp.stat <- sum(zScore.emp[ind[[j]]] - upper[ind[[j]]])
      sim.stat <- abs(apply(zScore.sim[,ind[[j]]], 1, function(x,y) { a = x-y; i = which(a>0); return(sum(a[i])) }, y = upper[ind[[j]]]))
      
      local[j,]$type <- '+'
      local[j,]$startDec <- min(empirical$data$x[ind[[j]]])
      local[j,]$endDec <- max(empirical$data$x[ind[[j]]])
      local[j,]$p.value <- round( 1 - ( length(sim.stat[sim.stat < emp.stat]) + 1 ) / ( nsims + 1 ), 5)
    }
  }
  
  if (tails==2) {
    ind <- split(below, cumsum(c(1,diff(below) > 1))); ind <- ind[which(lengths(ind) > 1)]
    if (length(ind)>0) {
      for (k in 1:NROW(ind)) {
        emp.stat <- sum(lower[ind[[k]]] - zScore.emp[ind[[k]]])
        sim.stat <- abs(apply(zScore.sim[,ind[[k]]], 1, function(x,y) { a = y-x; i = which(a>0); return(sum(a[i])) }, y = lower[ind[[k]]]))
        
        local[j+k,]$type <- '-'
        local[j+k,]$startDec <- min(empirical$data$x[ind[[k]]])
        local[j+k,]$endDec <- max(empirical$data$x[ind[[k]]])
        local[j+k,]$p.value <- round( 1 - ( length(sim.stat[sim.stat < emp.stat]) + 1 ) / ( nsims + 1 ), 5)
      }
    }
  }
  
  ## cleanup
  rownames(local) <- c()
  aux <- apply(local[,c(2,3)], 1, diff)
  ind <- which(aux <= 5*empirical$metadata$.prec_dec + 1000*.Machine$double.eps)
  local <- local[-ind,]
  rownames(local) <- c()
  
  ## output
  out <- c()
  out$metadata$nsims <- nsims
  out$metadata$conf <- conf
  out$metadata$tails <- tails
  out$metadata$normalise <- normalise
  out$metadata$global.pval <- global.p
  out$metadata$local.pval <- local
  
  out$data$empirical <- empirical$data
  out$data$null.hyp <- list(x = empirical$data$x, CE.mean = zMean, CE.upper = upCI)
  if (tails==2) { out$data$null.hyp$CE.lower <- loCI }
  if (save.sim) { out$data$simulated$Az <- sim.az; out$data$simulated$Unc <- sim.unc;}
  class(out) <- 'skyscapeR.sigTest'
  
  return(out)
}
sigTestDec <- compiler::cmpfun(sigTestDec_unc, options=list(enableJIT = 3))




# Additions for Plotting --------------------------------------------------
hpdi <- function(x, mass=0.954){
  if (class(x) == 'skyscapeR.spd') { grd <- x$data }
  if (class(x) == 'skyscapeR.pdf' | class(x) == 'list') { grd <- x$dec }
  sorted <- sort(grd$y, decreasing=TRUE)
  heightIdx = min( which( cumsum( sorted) >= sum(grd$y, na.rm=T) * mass ) )
  height = sorted[heightIdx]
  indices = which( grd$y >= height )
  gaps <- which(diff(indices) > 1)
  starts <- indices[c(1, gaps + 1)]
  ends <- indices[c(gaps, length(indices))]
  result <- cbind(startDec = grd$x[ends], endDec = grd$x[starts]) 
  if (result[2] < result[1]) { aux <- result[1]; result[1] <- result[2]; result[2] <- aux }
  return(result)
}


pval2stars <- function(p.value) {
  if (class(p.value)=='character') { p.value <- as.numeric(substr(p.value, 3,nchar(p.value))) }
  out <- symnum(p.value, corr = FALSE, na = FALSE,
                cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                symbols = c("***", "**", "*", "+", "ns"), legend=F)
  return(as.character(out))
}



plot.skyscapeR.pdf <- function(pdf, index, hdr=0.954, show.az=T, xlim, ylim, col=MESS::col.alpha('blue',0.5)) {
  original <- pdf 
  n <- length(pdf$metadata$az)
  pdf <- pdf$data
  
  if (missing(index)) { index <- 1:n }
  
  for (i in index) {
    x <- pdf[[i]]
    
    hd <- hpdi(x, hdr)
    
    if (!show.az) {
      par(mar=c(4,4,1,1))
      if(!missing(xlim)) { xl <- xlim } else { xl <- range(hpdi(x, 1)) }
      if(!missing(ylim)) { yl <- ylim } else { yl <- range(x$dec$y) * 1.2 }
      plot.default(x$dec$x, x$dec$y, type='l', main='', xlab='Declination (º)', ylab='Density', lwd=2, col='black', xlim=xl, ylim=yl, xaxs='i', yaxs='i', axes=F); box()
      
      for (j in 1:NROW(hd)) {
        ind <- which(x$dec$x >= min(hd[j,]) & x$dec$x <= max(hd[j,]))
        xvals <- x$dec$x[ind]; xvals <- c(xvals, rev(xvals))
        yvals <- x$dec$y[ind]; yvals <- c(yvals, rep(0, length(yvals)))
        polygon(xvals, yvals, col=col, border=NA)
      }
      lines(x$dec$x, x$dec$y, lwd=2, col='black')
      
      axis(1, at=pretty(seq(par('usr')[1],par('usr')[2])))
      axis(1, at=0, labels = 0)
      scale <- mean(diff(pretty(seq(par('usr')[1],par('usr')[2]))))
      if (scale <= 1) { axis(1, at=seq(-90,90,0.1), lwd=0.2, labels=F) }
      if (scale <= 2 & scale > 1) { axis(1, at=seq(-90,90,0.5), lwd=0.2, labels=F) }
      if (scale <= 5 & scale > 2) { axis(1, at=seq(-90,90,1), lwd=0.5, labels=F) }
      if (scale <= 20 & scale > 5) { axis(1, at=seq(-90,90,5), lwd=0.5, labels=F) }
      if (scale > 10) { axis(1, at=seq(-90,90,10), lwd=0.5, labels=F) }
      axis(2)
      
    } else {
      xx <- x$az$x
      dens <- x$az$y
      xl <- range(x$dec$x)
      yl <- range(xx);  yl[1] <- max(yl[1],0); yl[2] <- min(yl[2],360)
      
      layout(t(matrix(c(2,3,1,4), nrow=2)), widths = c(1,3), heights=c(3,1))
      
      # Information Panel
      par(mar=c(0,0,0,0))
      plot(-999,-999, xlim=c(0,1), ylim=c(0,1), xaxs='i', yaxs='i',axes=F, xlab='', ylab='')
      
      if (exists('name', where=original$metadata)) { text(.3, .8, labels=original$metadata$name[i], pos=4, font=2, cex=1.2)  }
      prob <- paste0(hdr*100, '% probability')
      for (j in 1:NROW(hd)) {
        prob <- paste0(prob,'\n   ',round(hd[j,1],2),'º  ::  ',round(hd[j,2],2),'º')
      }
      text(.3, .6, labels=prob, pos=4, font=1, cex=0.8)
      
      # Azimuth Distribution
      par(mar=c(0,4,1,0))
      plot(dens, xx, type='l', ylab='Azimuth (º)', xlab='', lwd=2, col='black', xlim=c(0, 1.2*max(dens)), ylim=yl, xaxs='i', yaxs='i', axes=F); box()
      polygon(c(dens,rep(0,length(dens))), c(xx, rev(xx)), col='grey', bordern=NA)
      axis(2, at=pretty(seq(par('usr')[3],par('usr')[4])))
      axis(2, at=0, labels = 0)
      scale <- mean(diff(pretty(seq(par('usr')[3],par('usr')[4]))))
      if (scale <= 1) { axis(2, at=seq(-45,360+45,0.1), lwd=0.2, labels=F) }
      if (scale <= 2 & scale > 1) { axis(2, at=seq(-45,360+45,0.5), lwd=0.2, labels=F) }
      if (scale <= 5 & scale > 2) { axis(2, at=seq(-45,360+45,1), lwd=0.5, labels=F) }
      if (scale <= 20 & scale > 5) { axis(2, at=seq(-45,360+45,5), lwd=0.5, labels=F) }
      if (scale > 10) { axis(2, at=seq(-45,360+45,10), lwd=0.5, labels=F) }
      
      # Transformation Curve
      par(mar=c(0,0,1,1))
      plot(-999,999, axes=F, xlim=xl, ylim=yl, xaxs='i', yaxs='i'); box()
      
      ## process horizon profile
      azs <- x$az$x
      hh <- original$metadata$horizon[[i]]$data
      alt <- approx(hh$az, hh$alt, xout=azs)$y
      alt.unc <- approx(hh$az, hh$alt.unc, xout=azs)$y
      refraction <- original$metadata$param$refraction
      atm <- original$metadata$param$atm
      temp <- original$metadata$param$temp
      
      dec0 <- az2dec(azs, original$metadata$horizon[[i]], alt, refraction=refraction, atm=atm, temp=temp)
      dec1 <- az2dec(azs, original$metadata$horizon[[i]], alt+alt.unc, refraction=refraction, atm=atm, temp=temp)
      dec2 <- az2dec(azs, original$metadata$horizon[[i]], alt-alt.unc, refraction=refraction, atm=atm, temp=temp)
      dec3 <- az2dec(azs, original$metadata$horizon[[i]], alt+2*alt.unc, refraction=refraction, atm=atm, temp=temp)
      dec4 <- az2dec(azs, original$metadata$horizon[[i]], alt-2*alt.unc, refraction=refraction, atm=atm, temp=temp)
      
      yp <- c(azs, rev(azs))
      xp <- c(dec3, rev(dec4))
      polygon(xp, yp, border=NA, col=MESS::col.alpha('blue', 0.25))
      yp <- c(azs, rev(azs))
      xp <- c(dec1, rev(dec2))
      polygon(xp, yp, border=NA, col=MESS::col.alpha('blue', 0.5))
      lines(dec0, azs, lwd=0.5, col='black')
      
      
      # Declination Distribution
      par(mar=c(4,0,0,1))
      plot.default(x$dec$x, x$dec$y, type='l', main='', xlab='Declination (º)', lwd=1.2, col='black', xlim=xl, ylim=c(0, 1.2*max(x$dec$y)), xaxs='i', yaxs='i', axes=F); box()
      
      for (j in 1:NROW(hd)) {
        ind <- which(x$dec$x >= min(hd[j,]) & x$dec$x <= max(hd[j,]))
        xvals <- x$dec$x[ind]; xvals <- c(xvals, rev(xvals))
        yvals <- x$dec$y[ind]; yvals <- c(yvals, rep(0, length(yvals)))
        polygon(xvals, yvals, col=col, border=NA)
        
      }
      lines(x$dec$x, x$dec$y)
      
      axis(1, at=pretty(seq(par('usr')[1],par('usr')[2])))
      axis(1, at=0, labels = 0)
      scale <- mean(diff(pretty(seq(par('usr')[1],par('usr')[2]))))
      if (scale <= 1) { axis(1, at=seq(-45,360+45,0.1), lwd=0.2, labels=F) }
      if (scale <= 2 & scale > 1) { axis(1, at=seq(-90,90,0.5), lwd=0.2, labels=F) }
      if (scale <= 5 & scale > 2) { axis(1, at=seq(-90,90,1), lwd=0.5, labels=F) }
      if (scale <= 20 & scale > 5) { axis(1, at=seq(-90,90,5), lwd=0.5, labels=F) }
      if (scale > 10) { axis(1, at=seq(-90,90,10), lwd=0., labels=F) }
      
      box()
      # par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
    }
  }
}



plot.skyscapeR.spd <- function(x, xlim, ylim, title=NULL, col='blue', shading=T) {
  par(mar=c(5, 4, 1, 1) + 0.1)
  if (missing(xlim)) { xlim <- sort(x$data$x[c(min(which(x$data$y >= 1e-12)), max(which(x$data$y >= 1e-12)))]) }
  if (missing(ylim)) { ylim <- c(0, max(x$data$y)) }
  plot.default(x$data$x, x$data$y, type='l', main=title, xlab='Declination (º)', ylab='Density', lwd=2, col=col, xlim=xlim, ylim=ylim, xaxs='i', yaxs='i', axes=F); box()
  if (shading) {
    xp <- c(x$data$x, rev(x$data$x))
    yp <- c(x$data$y, rep(0, length(x$data$x)))
    polygon(xp, yp, col=MESS::col.alpha(col,0.5), border=NA)
  }
  axis(1, at=pretty(seq(par('usr')[1],par('usr')[2])))
  axis(1, at=0, labels = 0)
  scale <- mean(diff(pretty(seq(par('usr')[1],par('usr')[2]))))
  if (scale <= 1) { axis(1, at=seq(-45,360+45,0.1), lwd=0.2, labels=F) }
  if (scale <= 2 & scale > 1) { axis(1, at=seq(-90,90,0.5), lwd=0.2, labels=F) }
  if (scale <= 5 & scale > 2) { axis(1, at=seq(-90,90,1), lwd=0.5, labels=F) }
  if (scale <= 20 & scale > 5) { axis(1, at=seq(-90,90,5), lwd=0.5, labels=F) }
  if (scale > 10) { axis(1, at=seq(-90,90,10), lwd=0.5, labels=F) }
  axis(2)
}



plot.skyscapeR.sigTest <- function(sig, xlim, title=NULL, show.pval=T, show.local=F) {
  # empirical spd
  spd <- sig$data$empirical
  if (missing(xlim)) { xlim <- sort(spd$x[c(min(which(spd$y >= 1e-12)), max(which(spd$y >= 1e-12)))]) }
  if (show.local) { ylim <- c(-max(spd$y)*.05,max(spd$y)) } else { ylim <- c(0,max(spd$y))}
  plot.default(spd$x, spd$y, type='l', main=title, xlab='Declination (º)', ylab='Density', lwd=2, col='blue', xlim=xlim, ylim=ylim, xaxs='i', yaxs='i', axes=F); box()
  xp <- c(spd$x, rev(spd$x))
  yp <- c(spd$y, rep(0, length(spd$x)))
  polygon(xp, yp, col=MESS::col.alpha('blue',0.5), border=NA)
  axis(1, at=pretty(seq(par('usr')[1],par('usr')[2])))
  axis(1, at=0, labels = 0)
  scale <- mean(diff(pretty(seq(par('usr')[1],par('usr')[2]))))
  if (scale <= 1) { axis(1, at=seq(-45,360+45,0.1), lwd=0.2, labels=F) }
  if (scale <= 2 & scale > 1) { axis(1, at=seq(-90,90,0.5), lwd=0.2, labels=F) }
  if (scale <= 5 & scale > 2) { axis(1, at=seq(-90,90,1), lwd=0.5, labels=F) }
  if (scale <= 20 & scale > 5) { axis(1, at=seq(-90,90,5), lwd=0.5, labels=F) }
  if (scale > 10) { axis(1, at=seq(-90,90,10), lwd=0.5, labels=F) }
  axis(2)
  
  # sigtest
  lines(spd$x, sig$data$null.hyp$CE.mean, col='grey')
  xp <- c(spd$x, rev(spd$x))
  if (sig$metadata$tails==1) {
    yp <- c(sig$data$null.hyp$CE.upper, rep(0, length(spd$x))) 
  } else {
    yp <- c(sig$data$null.hyp$CE.upper, rev(sig$data$null.hyp$CE.lower)) 
  }
  polygon(xp, yp, col=MESS::col.alpha('grey', 0.5), border=NA)
  
  # global p-value
  if (show.pval) { 
    if (sig$metadata$global.pval == 0 ) {
      pval <- paste0("global p-value < ", round(1/(sig$metadata$nsims+1),4))
    } else {
      pval <- paste0('global p-value = ', sig$metadata$global.pval)
    }
    text(par('usr')[2], abs(diff(par('usr')[3:4]))*.90, pos=2, pval, cex=1.2, font=2) 
  }
  
  # regions of significance
  if (show.local) {
    abline(0,0, lwd=1)
    aux <- as.matrix(sig$metadata$local.pval[,2:4])
    for (i in 1:NROW(aux)) {
      if (sig$metadata$local.pval[i,1] == '+') { col <- pal[2] } else { col <- pal[1] }
      if (aux[i,1] < xlim[1] & aux[i,2] > xlim[1]) { aux[i,1] <- xlim[1] }
      if (aux[i,2] > xlim[2] & aux[i,1] < xlim[2]) { aux[i,2] <- xlim[2] }
      xp <- c(aux[i,1], aux[i,1], aux[i,2], aux[i,2])
      yp <- c(ylim[1], 0, 0, ylim[1])
      polygon(xp, yp, col=MESS::col.alpha(col, .4), border=NA)
      text(mean(aux[i,1:2]), ylim[1]/2, labels=pval2stars(aux[i,3]), cex=0.9)
    }
  }
}


print.skyscapeR.sigTest <- function (x, ...) {
  cat("\n*** Results of Significance Test ***\n\n")
  cat(paste0(x$metadata$tails,'-tailed test at ',x$metadata$conf*100,'% confidence, based on ', x$metadata$nsims, ' simulations.\n'))
  
  if (x$metadata$global.pval == 0) { 
    p.value <- paste0("< ", round(1/(x$metadata$nsims+1),3)) 
  } else { p.value <- x$metadata$global.pval }
  cat(paste0('global p-value: ', p.value, ' (',pval2stars(p.value),')\n'))
  cat('local p-values:\n')
  for (i in 1:NROW( x$metadata$local.pval)) {
    if (x$metadata$local.pval[i,]$p.value == 0) { 
      p.value <- paste0("< ", round(1/(x$metadata$nsims+1),3)) 
    } else { p.value <- x$metadata$local.pval[i,]$p.value }
    cat(paste0('      ',x$metadata$local.pval[i,]$type,'  dec range [',round(x$metadata$local.pval[i,]$startDec,2), ', ',round(x$metadata$local.pval[i,]$endDec,2),'] :: p-value: ', p.value, ' (',pval2stars(p.value),')\n'))
  }
}
