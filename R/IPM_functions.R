# factor.mult: multiplicative factor to transform to m2/ha. Only needed for the first step.

# Survival, growth, recruitment... use the logarithm of the diameter as variable.

# It calculates the "survival X growth" integrand (but does not compute the integral) of the IPM.

# From the discrete-to-continuous case (i.e. when starting from IFN, where DBH is given
# for individual trees) no quadrature needs to be solved. After that step, further calculations
# require us to solve the quadrature, which is carried out with the so called "alternative extended
# Simpson's rule# (see "Numerical Recipes" by Press et al. 1989), which is a O(1/N^4) formula.
# library(Rcpp)
# sourceCpp("/home/david/CREAF/articulo/IPM/IPM_functions_v21.cpp")


ipm.adult.survival.times.growth <- function(ntrees,x,y,param.survival,param.growth,t.diff,max.diam,h=NULL,nx=NULL,y.minus.x=NULL) {
   q <-  ntrees * ipm.prob.survival(x=x,param=param.survival) * 
                  ipm.pdf.tree.growth(x,y,param.growth,t.diff,max.diam,nx,y.minus.x)
  if (!is.null(h)) q <- quad.trapez(q=q,h=h,nx=nx)
  return(q)
}


quad.trapez <- function(q,h,nx,dim.int=2) {
  if (dim.int==1) int.q <- sum(q) - .5*(q[1]+q[nx])
  else int.q <- colSums(q) - .5*(q[1,]+q[nx,])
  return(h*int.q)
}

ipm.adult.survival.times.growth.2 <- function(ntrees,x,y,param.survival,param.growth,t.diff,max.diam,h=NULL,nx=NULL,y.minus.x=NULL) {
  q <-  ntrees * #ipm.prob.survival(x=x,param=param.survival) * 
    ipm.pdf.tree.growth(x,y,param.growth,t.diff,max.diam,nx,y.minus.x)
  if (!is.null(h)) q <- quad.trapez(q=q,h=h,nx=nx)
  return(q)
}

# New version of the alternative extended Simpson's rule.
quad.ae.simpson <- function(q,h,nx,dim.int=2) {
  if (dim.int==1) int.q <- sum(q) - ((q[1]+q[nx])*17/48 + (q[2]+q[nx-1])*59/48 + 
    (q[3]+q[nx-2])*43/48 + (q[4]+q[nx-3])*49/48)
  else int.q <- colSums(q) - ((q[1,]+q[nx,])*17/48 + (q[2,]+q[nx-1,])*59/48 + 
    (q[3,]+q[nx-2,])*43/48 + (q[4,]+q[nx-3,])*49/48)
  return(h*int.q)
}

quad.trapez <- function(q,h,nx,dim.int=2) {
  if (dim.int==1) int.q <- sum(q) - .5*(q[1]+q[nx])
    else int.q <- colSums(q) - .5*(q[1,]+q[nx,])
  return(h*int.q)
}


# Survival function. We assume a logistic function to describe the dependence of tree survival on
# log diameter and the other parameters.
# param[1]: coefficient of log-diameter variable.
# param[2]: other terms, related to local climatic and tree stand characteristics.
ipm.prob.survival <- function(x,param) return(1/(1+exp(-(param[1]*log(x)+param[2]))))

f.exp <- function(x,param,t.diff,max.diam) {
  q <- (max.diam-x)*(1-exp(-t.diff*(param[1]*log(x)+param[2])))
  nq <- which(q>0)
  z <- q*0  # Fastest.
  z[nq] <- log(q[nq])
  if (length(nq)>0) z[-nq] <- mean(z[nq])
  return(z)
}
  

ipm.pdf.tree.growth <- function(x,y,param,t.diff,max.diam,nx=NULL,y.minus.x=NULL) {
  f.gr <- fExpCpp(x,param,t.diff,max.diam)
  variance <- exp(param[3] + param[4]*log(x))
  if (is.null(nx)) z <- dlnorm(y-x,meanlog=f.gr,sdlog=sqrt(variance))
    else z <- dlnorm(y.minus.x,meanlog=f.gr,sdlog=sqrt(variance))
  return(z)
}

# Ingrowth function, which is composed of two terms: one is an exponential PDF,
# which describes how new trees' sizes are distributed; and the other accounts
# for the number of trees.
# param[1]: lambda exponent for the exponential distribution.
# param[2]: number of seedlings and saplings times their respective coefficients.
# param[3]: local and environmental variables times their respective coefficients.

ipm.ingrowth.exp <- function(y,param.ingrowth) {
  return(param.ingrowth[1]*exp(-y*param.ingrowth[1])*exp(param.ingrowth[2]))
}
ipm.ingrowth.identity <- function(y,param.ingrowth) {
  return(param.ingrowth[1]*exp(-y*param.ingrowth[1])*param.ingrowth[2])
}
ipm.ingrowth.square <- function(y,param.ingrowth) {
  return(param.ingrowth[1]*exp(-y*param.ingrowth[1])*param.ingrowth[2]^2)
}

ipm.ingrowth.log <- function(y,param.ingrowth) {
  return(param.ingrowth[1]*exp(-y*param.ingrowth[1])*exp(param.ingrowth[2]))
}


# Sapling number calculation. methodology similar to ingrowth, in that
# we compute a probability * a number depending on the parameters

# actually we don't need the number of saplings, it's a remnant of earlier versions of the function

ipm.saplings <- function(nsapl,param) {
#   if (nsapl==0) {
#     pz <- param[1]
#     z <- param[2]
#   } else {
#     pz <- param[3]
#     z <- param[4]
#   }
  # undo anscombe
  #if (z>0) z <- (z/2)^2 - 1/8 else z<-0 
  return(param[1]*param[2])
}

# Multiply coefficients and variables.
# z[1]: rain
# z[2]: temper
# z[3]: basal area (m2 per ha).
# z[4]: rain anomaly.
# z[5]: temperature anomaly.
# m: stem number
coef.times.var <- function(param,z,m,s=NULL,spba,type) {
  
  if(type=="ingrowth"){
    q <-  param$intercept             + param$rain*z[1] +
      param$temper*z[2]           + param$basal.area*z[3] +
      param$anom.rain*z[4]        + param$anom.temper*z[5] +
      #param$temper.square*z[2]^2  + param$basal.area.square*z[3]^2 +
      #param$stems*m               + param$basal.area.stems*z[3]*m +
      #param$stems.square*m^2      
      param$rain.anom.temper*z[1]*z[5] 
      #param$rain.anom.rain*z[1]*z[4] + 
      #param$temper.anom.rain*z[2]*z[4] + 
      #param$temper.anom.temper*z[2]*z[5] 

      if (!is.null(s)) {
        q <- q + param$saplings*s
      }
    q <- q * (10000/(pi*25))
  }else if(type=="survival"){
    q <-  param$intercept             + param$rain*z[1] +
          param$temper*z[2]           + param$basal.area*z[3] +
          param$anom.rain*z[4]        + param$anom.temper*z[5] +
          #param$temper.square*z[2]^2  + param$basal.area.square*z[3]^2 +
          #param$stems*m               + param$basal.area.stems*z[3]*m +
          #param$stems.square*m^2      
          param$rain.anom.temper*z[1]*z[5]  
          #param$rain.anom.rain*z[1]*z[4] + 
          #param$temper.anom.rain*z[2]*z[4] + 
          #param$temper.anom.temper*z[2]*z[5] 
    
  }else if(type=="recruitment"){
    q <-  param$intercept             + param$rain*z[1] +
          param$temper*z[2]           + param$basal.area*z[3] +
          param$anom.rain*z[4]        + param$anom.temper*z[5] +
          #param$temper.square*z[2]^2  + param$basal.area.square*z[3]^2 +
          #param$stems*m               + param$basal.area.stems*z[3]*m +
          #param$stems.square*m^2      + 
          param$rain.anom.temper*z[1]*z[5] +
          param$basal.area.sp*spba    
          #param$rain.anom.rain*z[1]*z[4] + 
          #param$temper.anom.temper*z[2]*z[5] + 
          #param$temper.anom.rain*z[2]*z[4] 
    #if (!is.null(s)) q <- q + param$saplings*s # --V17
  }else if(type=="growth"){
    q <-  param$intercept +
          param$rain*z[1] +
          param$temper*z[2]           + param$basal.area*z[3] +
          param$anom.rain*z[4]        + param$anom.temper*z[5] +
          #param$temper.square*z[2]^2  + param$basal.area.square*z[3]^2 +
#           param$stems*m               + param$basal.area.stems*z[3]*m +
#           param$stems.square*m^2      + 
          param$rain.anom.temper*z[1]*z[5]
          #param$rain.anom.rain*z[1]*z[4] + 
          #param$temper.anom.rain*z[2]*z[4] + 
          #param$temper.anom.temper*z[2]*z[5] 
  }else q <- NULL
  return(q)
}

# Cutoff for ingrowth. Basal areas larger than 5 or 15 (conifers or planifolia, respectively)
# do not allow smaller trees to grow larger than 7.5 cm DBH.
ingrowth.cutoff <- function(species,basal.area,conifers) {
  i <- any(species==conifers)
  cutoff <- FALSE
  if (is.na(i)) {                 # It's not a conifers.
    if (basal.area>15) cutoff <- TRUE     # Basal area is too large. No ingrowth.
  } else {                        # It's a conifers.
    if (basal.area>5) cutoff <- TRUE      # Basal area is too large. No ingrowth.
  }
  return(cutoff)
}

################ Load results from survival regression.

load.survival <- function(name.file,species.name) {
  
  # Number of species.
  num <- length(species.name)
  
  # Labels.
  names.coef <- c("(Intercept)",
                  "log(dbh2)",
                  "lluvia",
                  "temper",
                  "abasal",
                  "anom_lluvia",
                  "anom_temper",
                  #"I(temper*temper)",
                  #"I(abasal*abasal)",
                  "I(lluvia * anom_temper)") ### CAREFUL: "I(rain*anom.temper)" DOES NOT WORK
                  #"I(lluvia*anom_lluvia)",
                  #"I(temper*anom_temper)",
                  #"I(temper*anom_lluvia)")
  
  # Create list with coefficients for species.
  dumb <- array(0,dim=num)
  survival.coef <- list(intercept=dumb,
                        log.dbh=dumb,
                        rain=dumb,
                        temper=dumb,
                        basal.area=dumb,
                        anom.rain=dumb,
                        anom.temper=dumb,
                        #temper.square=dumb,
                        #basal.area.square=dumb,
                        rain.anom.temper=dumb)
                        #rain.anom.rain=dumb,
                        #temper.anom.temper=dumb,
                        #temper.anom.rain=dumb)

  # Load data from regressions.
  load(name.file)
  
  # Save results in list.
  
  for (i in 1:length(species.name)) {
#     j <- match(species.name[i],survival.results$ESPECIE)
    j <- match(species.name[i],resultados$ESPECIE)
  
    # coef.out <- coefficients(survival.results$LOG.LOGIT[[j]])
    coef.out <- coefficients(resultados$LOG.LOGIT[[j]])
    
    # index.names <- match(names.coef,names(coefficients(survival.results$LOG.LOGIT[[j]])))
    index.names <- match(names.coef,names(coefficients(resultados$LOG.LOGIT[[j]])))
    
    if (!is.na(index.names[1])) survival.coef$intercept[i] <- coef.out[index.names[1]]
    if (!is.na(index.names[2])) survival.coef$log.dbh[i] <- coef.out[index.names[2]]
    if (!is.na(index.names[3])) survival.coef$rain[i] <- coef.out[index.names[3]]
    if (!is.na(index.names[4])) survival.coef$temper[i] <- coef.out[index.names[4]]
    if (!is.na(index.names[5])) survival.coef$basal.area[i] <- coef.out[index.names[5]]
    #if (!is.na(index.names[6])) survival.coef$stems[i] <- coef.out[index.names[6]]
    if (!is.na(index.names[6])) survival.coef$anom.rain[i] <- coef.out[index.names[6]]
    if (!is.na(index.names[7])) survival.coef$anom.temper[i] <- coef.out[index.names[7]]
    #if (!is.na(index.names[8])) survival.coef$temper.square[i] <- coef.out[index.names[8]]
    #if (!is.na(index.names[10])) survival.coef$stems.square[i] <- coef.out[index.names[10]]
    #if (!is.na(index.names[9])) survival.coef$basal.area.square[i] <- coef.out[index.names[9]]
    #if (!is.na(index.names[12])) survival.coef$basal.area.stems[i] <- coef.out[index.names[12]]
    if (!is.na(index.names[8])) survival.coef$rain.anom.temper[i] <- coef.out[index.names[8]]
    #if (!is.na(index.names[9])) survival.coef$rain.anom.rain[i] <- coef.out[index.names[9]]
    #if (!is.na(index.names[10])) survival.coef$temper.anom.temper[i] <- coef.out[index.names[10]]
    #if (!is.na(index.names[11])) survival.coef$temper.anom.rain[i] <- coef.out[index.names[11]]
  }
  survival.results <- NULL ; rm(survival.results) ; gc()
  return(survival.coef)
}


################ Maximum diameter per species.
# We give the option that all maximum diameters are located at the nodes
# of the quadrature mesh. In this way, the number of nodes used to solve
# the quadrature varies with the species. Otherwise, the maximum diameter
# is returned.

x.per.species <- function(min.dbh=min.dbh,n.intervals=NULL) {
#  setwd("C:\\Roberto\\CONSOLIDER\\Modelo forestal IPM Espa?a\\Regresiones\\")
  load("data/regresiones/growth_v6")
  max.dbh <- results.growth$MAXDBH[]
  rm(results.growth)
  gc()
  x <- sapply(1:length(max.dbh), function(i) min.dbh + c(0:n.intervals)*(max.dbh[i]-min.dbh)/n.intervals)
  return(x)
}


################ Load results from growth regression.

load.growth <- function(name.file,species.name) {
  
  # Number of species.
  num <- length(species.name)
  
  # Create list with coefficients for species.
  dumb <- array(0,dim=num)
  growth.coef <- list(intercept=dumb,
                      log.dbh=dumb,
                      rain=dumb,
                      temper=dumb,
                      #temper.square=dumb,
                      basal.area=dumb,
                      #basal.area.square=dumb,
                      anom.rain=dumb,
                      anom.temper=dumb,
                      rain.anom.temper=dumb,
                      #rain.anom.rain=dumb,
                      #temper.anom.temper=dumb,
                      #temper.anom.rain=dumb,
                      intercept.variance=dumb,
                      slope.variance=dumb)
  
  # Load data from regressions.
  load(name.file)

  # Save results in list.
  
  for (i in (1:length(species.name))) {
    j <- match(species.name[i],results.growth$SPECIES)    # Matches species name.
#    coef.out <- coefficients(results.growth$REGRES[[j]])  # Extracts coefficients.
    coef.out <- results.growth$REGRES[[j]]$par
    var.out <- coefficients(results.growth$VARIANCE[[j]])
    growth.coef$intercept.variance[i] <- var.out[1]
    growth.coef$slope.variance[i] <- var.out[2]
    kk <- 0
    if (is.null(results.growth$PARAMS.ELIM[[j]])) {
      dumb <- coef.out
    } else {
      dumb <- rep(0,8)
      for (k in 1:8) {
        if (length(which(results.growth$PARAMS.ELIM[[j]]==k))==0) { # This variable was not eliminated.
          kk <- kk + 1
          dumb[k] <- coef.out[kk]
        }      
      }
    }
    growth.coef$intercept[i] <- dumb[1]
    growth.coef$log.dbh[i] <- dumb[2]
    growth.coef$rain[i] <- dumb[3]
    growth.coef$temper[i] <- dumb[4]
    #growth.coef$temper.square[i] <- dumb[5]
#     growth.coef$stems[i] <- dumb[6]
#     growth.coef$stems.square[i] <- 0.
    growth.coef$basal.area[i] <- dumb[5]
    #growth.coef$basal.area.square[i] <- dumb[7]
#     growth.coef$basal.area.stems[i] <- dumb[8]
    growth.coef$anom.rain[i] <- dumb[6]
    growth.coef$anom.temper[i] <- dumb[7]
    growth.coef$rain.anom.temper[i] <- dumb[8]
    #growth.coef$rain.anom.rain[i] <- dumb[9]
    #growth.coef$temper.anom.temper[i] <- dumb[10]
    #growth.coef$temper.anom.rain[i] <- dumb[11]
  }
  return(growth.coef)
}

################ Load results from ingrowth regression.

load.ingrowth <- function(name.file,species.name) {
  
# Number of species.
  num <- length(species.name)
  
# Labels.
  names.coef <- c("(Intercept)",
                  "saplings",
                  "rain",
                  "temper",
                  "basal.area",
                  "anom.rain",
                  "anom.temper",
                  #"I(temper*temper)",
                  #"I(basal.area*basal.area)",
                  "I(rain * anom.temper)") ### CAREFUL: "I(rain*anom.temper)" DOES NOT WORK
                  #"I(rain*anom.rain)",
                  #"I(temper*anom.temper)",
                  #"I(temper*anom.rain)")
  
# Create list with coefficients for species.
  dumb <- array(0,dim=num)
  ingrowth.coef <- list(lambda=dumb,
                        intercept=dumb,
                        saplings=dumb,
                        rain=dumb,
                        temper=dumb,
                        basal.area=dumb,
                        anom.rain=dumb,
                        anom.temper=dumb,
                        #temper.square=dumb,
                        #basal.area.square=dumb,
                        rain.anom.temper=dumb)
                        #rain.anom.rain=dumb,
                        #temper.anom.temper=dumb,
                        #temper.anom.rain=dumb)
  
  # Load data from regressions.
  load(name.file)

# Save results in list.
  # ingrowth.coef$lambda <- ingrowth.results$LAMBDA
  ingrowth.coef$lambda <- out$LAMBDA
  
  for (i in 1:length(species.name)) {
    # j <- match(species.name[i],ingrowth.results$SPECIES)
    j <- match(species.name[i],out$SPECIES)
    
    # coef.out <- coefficients(ingrowth.results$GLM.IDENTITY[[j]])
    coef.out <- coefficients(out$GLM.IDENTITY[[j]])
    
    # index.names <- match(names.coef,names(coefficients(ingrowth.results$GLM.IDENTITY[[j]])))
    index.names <- match(names.coef,names(coefficients(out$GLM.IDENTITY[[j]])))
    
    if (!is.na(index.names[1])) ingrowth.coef$intercept[i] <- coef.out[index.names[1]]
    if (!is.na(index.names[2])) ingrowth.coef$saplings[i] <- coef.out[index.names[2]]
    if (!is.na(index.names[3])) ingrowth.coef$rain[i] <- coef.out[index.names[3]]
    if (!is.na(index.names[4])) ingrowth.coef$temper[i] <- coef.out[index.names[4]]
    if (!is.na(index.names[5])) ingrowth.coef$basal.area[i] <- coef.out[index.names[5]]
    #if (!is.na(index.names[6])) ingrowth.coef$stems[i] <- coef.out[index.names[6]]
    if (!is.na(index.names[6])) ingrowth.coef$anom.rain[i] <- coef.out[index.names[6]]
    if (!is.na(index.names[7])) ingrowth.coef$anom.temper[i] <- coef.out[index.names[7]]
    #if (!is.na(index.names[8])) ingrowth.coef$temper.square[i] <- coef.out[index.names[8]]
    #if (!is.na(index.names[10])) ingrowth.coef$stems.square[i] <- coef.out[index.names[10]]
    #if (!is.na(index.names[9])) ingrowth.coef$basal.area.square[i] <- coef.out[index.names[9]]
    #if (!is.na(index.names[12])) ingrowth.coef$basal.area.stems[i] <- coef.out[index.names[12]]
    if (!is.na(index.names[8])) ingrowth.coef$rain.anom.temper[i] <- coef.out[index.names[8]]
    #if (!is.na(index.names[9])) ingrowth.coef$rain.anom.rain[i] <- coef.out[index.names[9]]
    #if (!is.na(index.names[10])) ingrowth.coef$temper.anom.temper[i] <- coef.out[index.names[10]]
    #if (!is.na(index.names[11])) ingrowth.coef$temper.anom.rain[i] <- coef.out[index.names[11]]
  }
  return(ingrowth.coef)
}

################ Load results from sapling regression.

load.saplings <- function(name.file,species.name) {
  
  # Load data from regressions.
  load(name.file)

# Number of species.
  num <- length(species.name)
  
# Labels.
  names.coef <- c("(Intercept)",
                  #"saplings",
                  "rain",
                  "temper",
                  "basal.area",
                  "basal.area.sp",
                  #"stems",
                  "anom.rain",
                  "anom.temper",
                  #"I(temper*temper)",
                  #"I(stems*stems)",
                  #"I(basal.area*basal.area)",
                  #"I(basal.area*stems)",
                  "I(rain * anom.temper)") ### CAREFUL: "I(rain*anom.temper)" DOES NOT WORK
                  #"I(rain*anom.rain)",
                  #"I(temper*anom.temper)",
                  #"I(temper*anom.rain)")
 
# Create list with coefficients for species. First come the cases where there are large trees but no saplings.
  dumb <- array(0,dim=num)
  sapl.coef <- list(binom=dumb,
                      intercept=dumb,
                      rain=dumb,
                      temper=dumb,
                      basal.area=dumb,
                      basal.area.sp=dumb,
                      #stems=dumb,
                      anom.rain=dumb,
                      anom.temper=dumb,
                      #temper.square=dumb,
                      #stems.square=dumb,
                      #basal.area.square=dumb,
                      #basal.area.stems=dumb,
                      rain.anom.temper=dumb)
                      #rain.anom.rain=dumb,
                      #temper.anom.temper=dumb,
                      #temper.anom.rain=dumb)

# Save results in list.
  # sapl.coef$binom <- recruitment.results$binom
  sapl.coef$binom <- out$binom
  
  for (i in 1:length(species.name)) {
    
    # j <- match(species.name[i],recruitment.results$SPECIES)
    j <- match(species.name[i],out$SPECIES)
    
    if (!is.na(j)) {
      # coef.out <- coefficients(recruitment.results$poisson[[j]])
      coef.out <- coefficients(out$poisson[[j]])
      
      # index.names <- match(names.coef,names(coefficients(recruitment.results$poisson[[j]])))
      index.names <- match(names.coef,names(coefficients(out$poisson[[j]])))
      
      if (!is.na(index.names[1])) sapl.coef$intercept[i] <- coef.out[index.names[1]]
      if (!is.na(index.names[2])) sapl.coef$rain[i] <- coef.out[index.names[2]]
      if (!is.na(index.names[3])) sapl.coef$temper[i] <- coef.out[index.names[3]]
      if (!is.na(index.names[4])) sapl.coef$basal.area[i] <- coef.out[index.names[4]]
      if (!is.na(index.names[5])) sapl.coef$basal.area.sp[i] <- coef.out[index.names[5]]
      #if (!is.na(index.names[6])) sapl.coef1$stems[i] <- coef.out[index.names[6]]
      if (!is.na(index.names[6])) sapl.coef$anom.rain[i] <- coef.out[index.names[6]]
      if (!is.na(index.names[7])) sapl.coef$anom.temper[i] <- coef.out[index.names[7]]
      #if (!is.na(index.names[9])) sapl.coef1$temper.square[i] <- coef.out[index.names[9]]
      #if (!is.na(index.names[8])) sapl.coef1$stems.square[i] <- coef.out[index.names[8]]
      #if (!is.na(index.names[11])) sapl.coef1$basal.area.square[i] <- coef.out[index.names[11]]
      #if (!is.na(index.names[10])) sapl.coef1$basal.area.stems[i] <- coef.out[index.names[10]]
      if (!is.na(index.names[8])) sapl.coef$rain.anom.temper[i] <- coef.out[index.names[8]]
      #if (!is.na(index.names[12])) sapl.coef1$rain.anom.rain[i] <- coef.out[index.names[12]]
      #if (!is.na(index.names[13])) sapl.coef1$temper.anom.temper[i] <- coef.out[index.names[13]]
      #if (!is.na(index.names[14])) sapl.coef1$temper.anom.rain[i] <- coef.out[index.names[14]]
    }
  }
  
return(sapl.coef)
  #return(list(sapl.coef1,sapl.coef2))
}

################ For the ingrowth regression

Evaluate_lambda_NewtonRaphson <- function(lambda,xtrunc,xmedia) {
  
  laexpon <- function(x) exp(-x*xtrunc)/(1-exp(-x*xtrunc))
  lafuncion <- function(x) 1/x-xtrunc*laexpon(x)-xmedia
  lafuncionderivada <- function(x) xtrunc^2*(laexpon(x)+laexpon(x)^2)-1/x^2
  lambda1 <- lambda
  lambda2 <- 1e32
  while (abs(lambda1-lambda2)>1e-7) {
    lambda2 <- lambda1
    lambda1 <- lambda2 - lafuncion(lambda2)/lafuncionderivada(lambda2)
  }
  return(lambda1)
}          




