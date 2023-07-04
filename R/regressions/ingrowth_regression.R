
# obtain ingrowth coefficients

# -------------------------------------------------------------------------

source("R/auxiliary_functions.R")
#source("IPM_functions_v5.R")
  
  ###########################################
  # Read data
  
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
  fire$damage <- pmax(fire$DANYSEVER,fire$DANYMITJA,fire$DANYPETIT)
  fire <- subset(fire,DANYSEVER > 0 | DANYMITJA > 0 | DANYPETIT > 0)
  
  # We select surviving trees in plots without management. We also eliminate surviving trees which show no measurable growth
  # in 10 years: we assume that to be a mistake.
  
  j <- which(trees$idestatus != "n" & trees$gestion_ifn3 == "N") 
  trees <- trees[j,]
  j <- match(trees$ID,fire$ID,nomatch=0)
  trees <- trees[-j,]
  trees$abasal <- trees$FACTOR_IFN2*pi*(trees$dbh_ifn2/200)^2
  
  trees.basal <- tapply(trees$abasal,trees$ID,sum)
  trees.basal <- list(IDPLOT=unlist(dimnames(trees.basal)),AREA=unlist(trees.basal[]))
  
  stems <- tapply(trees$FACTOR_IFN2,trees$ID,function(x) sum(x[!is.na(x)]))
  stems <- list(IDPLOT=unlist(dimnames(stems)),NUM=unlist(stems[]))

# Ingrowth.
  i <- which(trees$gestion_ifn3=="N" & trees$dbh_ifn2<12.5 & trees$idestatus=="i" & (trees$dbh_ifn2 != trees$dbh_ifn3))
  ingrowth <- trees[i,]

# saplings
  pm2 <- read.csv2("data/SAPLINGS_IFN2_v8.csv"
                   ,dec=".",sep=";",comment.char="")
  j <- match(pm2$ID,fire$ID,nomatch=0)
  pm2 <- pm2[-j,]
  pm2$name.sp <- tesauro[match(pm2$num.sp,tesauro[,1]),2]
  
# El an?lisis se hace por especie.
  
  out <- list(SPECIES=array("",dim=length(newspecies)),LAMBDA=array(0,dim=length(newspecies)),
    GLM.IDENTITY=array(list(),dim=length(newspecies)))

  for (i in (1:length(newspecies))) {

    print(paste(i," - ",newspecies[i],sep=""))

# Calculates the ingrowth distribution first.
    j <- which(ingrowth$name.sp==newspecies[i])
    sub.ingrowth <- ingrowth[j,]
    sub.ingrowth <- sub.ingrowth[which(complete.cases(sub.ingrowth)),]
    id.plots.ingrowth <- unique(sub.ingrowth$ID)     # Just in case...
    incorporated <- sapply(id.plots.ingrowth,function(x) sum(sub.ingrowth$ID==x))

# Usamos m?xima verosimilitud para ajustar una exponencial a los datos. Dado que la
# distribuci?n de nuevos pies mayores con di?metros entre 7.5 y 12.5 est? truncada a 12.5                               
# tenemos que usar los resultados del trabajo de Deemer y Votaw (1955). Con esta lambda
# caracterizaremos la distribuci?n de di?metros de los ?rboles que se incorporan a
# los adultos desde los pies menores.

    mean_DN <- mean(sub.ingrowth$dbh_ifn3)- 7.5
    x0 <- 5
    out$LAMBDA[i] <- Evaluate_lambda_NewtonRaphson(1,x0,mean_DN)
    
# We also need to know the number of trees, and their basal area, of species "newspecies[i]"
# which are present in the stand.
    j <- which(trees$name.sp==newspecies[i])
    subtrees <- trees[j,]

# List to be used in the regression.
    data.ingrowth <- list()
    data.ingrowth$incorporated <- incorporated
    j <- match(id.plots.ingrowth,clima$ID)
    data.ingrowth$rain <- clima$ref_pl[j]
    data.ingrowth$temper <- clima$ref_mt[j]
    data.ingrowth$anom.rain <- clima$anom_pl[j]
    data.ingrowth$anom.temper <- clima$anom_mt[j]
    #data.ingrowth$stems <- stems$NUM[match(id.plots.ingrowth,stems$IDPLOT)]
    data.ingrowth$basal.area <- trees.basal$AREA[match(id.plots.ingrowth,trees.basal$IDPLOT)]
    j <- match(id.plots.ingrowth,pm2$ID)
    data.ingrowth$seedlings <- pm2$seedlings[j]
    data.ingrowth$saplings <- pm2$saplings[j]#*10000/(pi*25)
    data.ingrowth$seedlings[is.na(data.ingrowth$seedlings)] <- 0
    data.ingrowth$saplings[is.na(data.ingrowth$saplings)] <- 0

    data.ingrowth <- as.data.frame(data.ingrowth)
    out$SPECIES[i] <- newspecies[i]
    data.ingrowth$sqrt.incorp <- sqrt(data.ingrowth$incorporated)
    
# First, a lm for providing starting points. Then, glm with log link function. See IPM_functions_vX, the function "ipm.saplings".
# there we undo the link for the parameters.
    
    xnam <- c(names(data.ingrowth[1,c(-1,-7,-9)]),
              #"I(temper*temper)",
              #"I(basal.area*basal.area)",
              #"I(temper*anom.temper)",
              #"I(temper*anom.rain)",
              "I(rain*anom.temper)")
              #"I(rain*anom.rain)")
    naic <- sum(sapply(0:6,function(j) choose(7,j)))
    
# All variables.
    fmla1 <- as.formula(paste("sqrt.incorp ~ ", paste(xnam, collapse= "+")))
    fmla2 <- as.formula(paste("incorporated ~ ", paste(xnam, collapse= "+")))

    r.out.lm <- lm(fmla1, data=data.ingrowth)
    start.coef <- coefficients(r.out.lm)
    r.out.glm <- glm(fmla2,family=poisson(link="sqrt"), data=data.ingrowth, model=TRUE, start=start.coef)
    start.coef <- coefficients(r.out.glm)
    r.out.glm <- glm(fmla2,family=poisson(link="identity"), data=data.ingrowth, model=TRUE, start=start.coef)

    out.old <- r.out.glm
    aic.keep.old <- r.out.glm$aic
    print(paste(newspecies[i]," ",i," -- ",1," of ",naic," -- AIC=",aic.keep.old,sep=""))
    p <- 2

    for (j in 1:6) {
      var.remove <- combn(7,j)
      for (k in 1:dim(var.remove)[2]) {
        xnam.aic <- xnam[-var.remove[,k]]
        fmla1 <- as.formula(paste("sqrt.incorp ~ ", paste(xnam.aic, collapse= "+")))
        fmla2 <- as.formula(paste("incorporated ~ ", paste(xnam.aic, collapse= "+")))
        r.out.lm <- try(lm(fmla1, data=data.ingrowth))
        start.coef <- coefficients(r.out.lm)
        r.out.glm <- try(glm(fmla2,family=poisson(link="sqrt"), data=data.ingrowth, model=TRUE, start=start.coef))
        if ("try-error" %in% attr(r.out.glm,"class")) {
          aic.keep.new <- 1e64
        } else {
          start.coef <- coefficients(r.out.glm)
          r.out.glm <- try(glm(fmla2,family=poisson(link="identity"), data=data.ingrowth, model=TRUE, start=start.coef))
          if ("try-error" %in% attr(r.out.glm,"class")) {
            aic.keep.new <- 1e64
          } else {
            aic.keep.new <- r.out.glm$aic
            if (aic.keep.new<=aic.keep.old) {
              out.old <- r.out.glm
              aic.keep.old <- aic.keep.new
            }
          }
        }
        if (round(p/100)==p/100) print(paste(newspecies[i]," ",i," -- ",p," of ",naic," -- AIC=",aic.keep.old,sep=""))
        p <- p+1
      }
    }
    out$GLM.IDENTITY[[i]] <- out.old
  }
  # save(out,file="data/regresiones/ingrowth_v10",compress=TRUE)                
  