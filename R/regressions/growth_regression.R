
# obtain growth as a function of covariates
# -------------------------------------------------------------------------

source("R/auxiliary_functions.R")

###########################################
# Read data

# Read in datasets. First the labels for the species, then the large dataset, next 
# the climatic data and, finally, the fire data.

sp.list <- read.csv2("data/lista_especies_v2.csv"
                     ,dec=".",sep=";",comment.char="")
newspecies <- unique(sp.list[,3])
newspecies <- newspecies[trim(newspecies)!=""]
tesauro <- data.frame(num.sp = c(1:length(newspecies)),
                      name.sp = newspecies)
# a = trees

trees <- read.csv2("data/PIES_MAYORES_IFN2_IFN3_v1.csv"
                   ,dec=".",sep=";",comment.char="")
trees$name.sp <- tesauro[match(trees$num.sp,tesauro[,1]),2]

# b = clima

clima <- read.table(file="data/clima/clima_referencia_v5.csv",header=T,sep=";",dec=".")
fire <- read.table(file="data/fire_damage_v2.csv",header=T,sep=";",dec=".")

# We select surviving trees in plots without management. We also eliminate surviving trees which show no measurable growth
# in 10 years: we assume that to be a mistake.

index <- which(trees$idestatus=="s" & (trees$dbh_ifn2 != trees$dbh_ifn3) & trees$gestion_ifn3=="N")
trees <- trees[index,]

# Calculates the basal area (in m2) per plot.

trees$abasal <- pi*(trees$dbh_ifn2/200)^2
aabasal <- tapply(trees$FACTOR_IFN2*trees$abasal,trees$ID,function(x) sum(x[!is.na(x)]))
aabasal <- list(IDPLOT=unlist(dimnames(aabasal)),AREA=unlist(aabasal[]))

# Added on 16-7-2012.

stems.all <- tapply(trees$FACTOR_IFN2,trees$ID,function(x) sum(x[!is.na(x)]))
stems.all <- list(IDPLOT=unlist(dimnames(stems.all)),NUM=unlist(stems.all[]))

# The analysis is done per species.

#results.growth <- list(SPECIES=newspecies,MAXDBH=array(0,dim=length(newspecies)),REGRES=array(list(),dim=length(newspecies)),
#                       VARIANCE=array(list(),dim=length(newspecies)),PARAMS.ELIM=array(list(),dim=length(newspecies)))

LM=array(list(),dim=21)
results.growth <- list(SPECIES=newspecies,MAXDBH=array(0,dim=length(newspecies)),REGRES=array(list(),dim=length(newspecies)),
                       VARIANCE=array(list(),dim=length(newspecies)),PARAMS.ELIM=array(list(),dim=length(newspecies)))

# Maximum diameters that can be reached, in cm. Info obtained mostly from English Wikipedia (WK),
# the "Arbres Monumentals. Tresors Naturals dels Ports" poster (AM), the IFN2 and IFN3 National Inventories (IFN),
# the http://www.monumentaltrees.com web (MT) or the "Boscos singulars de Catalunya" database (BS).
# I put here the largest of the five.

results.growth$MAXDBH[1] <- -1    # conifers.
results.growth$MAXDBH[2] <- -1    # Deciduous.
results.growth$MAXDBH[3] <- 300   # Fagus sylvatica.
results.growth$MAXDBH[4] <- 200   # Juniperus thurifera.
results.growth$MAXDBH[5] <- 146   # Pinus halepensis.
results.growth$MAXDBH[6] <- 148   # Pinus nigra.
results.growth$MAXDBH[7] <- 180   # Pinus pinaster.
results.growth$MAXDBH[8] <- 137   # Pinus pinea.
results.growth$MAXDBH[9] <- 149   # Pinus sylvestris.
results.growth$MAXDBH[10] <- 129  # Pinus uncinata.
results.growth$MAXDBH[11] <- 180  # Quercus faginea.
results.growth$MAXDBH[12] <- 152  # Quercus ilex.
results.growth$MAXDBH[13] <- 172  # Quercus pyrenaica.
results.growth$MAXDBH[14] <- 285  # Quercus robur/petraea.
results.growth$MAXDBH[15] <- 147  # Quercus suber
results.growth$MAXDBH[16] <- -1   # Sclerophyllous

# Now, for those tree groups for which a maximum diameter is not well defined (=-1),
# we look it up in the database. We use both managed and unmanaged plots in this step.
# Then we multiply this maximum diameter by an ad-hoc 1.1, guessing that the truly
# maximum diameter for every species will be slightly larger.

for (i in 1:16) {
  results.growth$SPECIES[i] <- newspecies[i]
  if (results.growth$MAXDBH[i]==-1) {
    iaaa <- which(trees$name.sp==newspecies[i])
    trees.temp <- trees[iaaa,]
    results.growth$MAXDBH[i] <- max(trees.temp$dbh_ifn2[!is.na(trees.temp$dbh_ifn2)],trees.temp$dbh_ifn3[!is.na(trees.temp$dbh_ifn3)])
  }
  results.growth$MAXDBH[i] <- results.growth$MAXDBH[i]*1.1
}

levels.Akaike <- 3

####################################################################################
# Functions to be used below.
####################################################################################

# Akaike index.
Akaike.corrected <- function(k,n,LogLik) {
  ak <- 2*k - 2*LogLik          # General Akaike index.
  akc <- ak + 2*k*(k+2)/(n-k-1) # Correction for small sample size.
  return(akc)
}

f.rate <- function(param) {    
  index <- 1:11
  param2 <- index*0
  if (length(kc)>0) index <- index[-kc]
  param2[index] <- param
  dumb <- param2[1]
  dumb <- dumb + param2[2]*log(d2)           + param2[3]*rain       + param2[4]*temper       + param2[5]*temper^2
  dumb <- dumb + param2[6]*stems             + param2[7]*basal.area
  dumb <- dumb + param2[8]*basal.area*stems + param2[9]*anom.rain + param2[10]*anom.temper + param2[11]*rain*anom.temper
  return(dumb)
}

# Sinusoidal function to be fitted.
f.exp <- function(param) diam.max + (d2-diam.max)*exp(-t.diff*f.rate(param))

loglik.logn <- function(param) return(-sum(log(dlnorm(d32,meanlog=f.exp(param)-d2,sdlog=sqrt(Variance)))))

loglik.rate <- function(param) {
  x <- t.diff*f.rate(param)
  x <- 1-exp(-x)
  return(sum((y-x)^2)) 
}

# First estimation of variance and initial estimates for parameters.

init.param <- function() {
  param <- rep(0,11-length(kc))
  out <- optim(param,loglik.rate,control=list(maxit=1000))
  return(out$par)
}

calc.param <- function(param) {
  out <- optim(param,loglik.logn,control=list(maxit=1000))
  param <- out$par
  v1 <- out$value
  out <- optim(param,loglik.logn,control=list(maxit=1000))
  param <- out$par
  v2 <- out$value
  out <- optim(param,loglik.logn,control=list(maxit=1000))
  param <- out$par
  v3 <- out$value
  for (j in 1:50)  {
    out <- optim(param,loglik.logn,control=list(maxit=1000))
    param <- out$par
    v4 <- out$value
    if (abs(1-v4/v3)<.01 & abs(1-v3/v2)<.01 & abs(1-v2/v1)<.01) break
    v1 <- v2
    v2 <- v3
    v3 <- v4
  }
  return(out)
}

# Computes a linear regression of the log of the variance.
calc.variance <- function(param) {
  resid <- (log(d3-d2)-(f.exp(param)-d2))^2
  log.d2 <- log(d2)
  log.d2[resid==0] <- mean(log.d2)
  resid[resid==0] <- mean(resid)
  log.resid <- log(resid)
  return(lm(log.resid ~ log.d2))
}

#####################################################################################################

# Main loop.

for (i in (1:4)) {
  
  print(paste(i," - ",newspecies[i],sep=""))
  
  iaaa <- which(trees$name.sp==newspecies[i])
  my.trees <- trees[iaaa,]
  
  basal.area <- aabasal$AREA[match(my.trees$ID,aabasal$IDPLOT)]
  
# Datasets used for the regression. We eliminate some NA data.
  
  t.diff <- 10
  d2 <- my.trees$dbh_ifn2
  d3 <- my.trees$dbh_ifn3
  diam.max <- results.growth$MAXDBH[i]  
  rain <- clima$ref_pl[match(my.trees$ID,clima$ID)]
  temper <- clima$ref_mt[match(my.trees$ID,clima$ID)]
  anom.rain <- clima$anom_pl[match(my.trees$ID,clima$ID)]
  anom.temper <- clima$anom_mt[match(my.trees$ID,clima$ID)]
  stems <- stems.all$NUM[match(my.trees$ID,stems.all$IDPLOT)]
  
  i.NA <- which(!is.na(rain) & !is.na(stems) & !is.na(basal.area))
  d2 <- d2[i.NA]
  d3 <- d3[i.NA]
  d32 <- d3-d2
  log.d32 <- log(d32)
  rain <- rain[i.NA]
  temper <- temper[i.NA]
  anom.rain <- anom.rain[i.NA]
  anom.temper <- anom.temper[i.NA]
  stems <- stems[i.NA]
  basal.area <- basal.area[i.NA]

# To store results.
  out.regress <- array(list(),dim=1)
  param.regress <- out.regress
  variance.regress <- out.regress
  AIC.regress <- array(0,dim=1)
  
  ntotal <- 1
  for (j in 1:levels.Akaike) ntotal <- ntotal + dim(combn(c(2:11),j))[2]
  
  index.out <- 0
  for (j in 0:levels.Akaike) {
    if (j==0) {
      k.comb <- NULL
      n.k <- 1
    } else {
      k.comb <- combn(c(2:11),j)
      n.k <- dim(k.comb)[2]

    }
    for (k in 1:n.k) {
      if (j>0) {
        out.regress <- c(out.regress,array(list(),dim=1))
        param.regress <- c(param.regress,array(list(),dim=1))
        variance.regress <- c(variance.regress,array(list(),dim=1))
        AIC.regress <- c(AIC.regress,array(0,dim=1))
      }
      kc <- as.vector(k.comb[,k])
      Variance <- var(log(d3-d2))
      y <- (d3-d2)/(diam.max-d2)
      param <- init.param()
      out <- calc.param(param)
      dummy <- calc.variance(out$par)
      Variance <- exp(unname(unlist(predict(dummy))))
      index.out <- index.out + 1
      print(paste(sep="",newspecies[i],". Loop ",index.out," of ",ntotal,". Akaike full fit = ",Akaike.corrected(11-j,NROW(d2),loglik.logn(out$par))))
      AIC.regress[index.out] <- Akaike.corrected(11-j,NROW(d2),loglik.logn(out$par))
      out.regress[[index.out]] <- out
      param.regress[[index.out]] <- kc
      variance.regress[[index.out]] <- dummy
    }
  }
  results.growth$REGRES[[i]] <- out.regress[[which.min(AIC.regress)]]
  results.growth$PARAMS.ELIM[[i]] <- param.regress[[which.min(AIC.regress)]]
  results.growth$VARIANCE[[i]] <- variance.regress[[which.min(AIC.regress)]]
}
save(results.growth,file=paste("data/regresiones/growth_v6"),compress=TRUE)


