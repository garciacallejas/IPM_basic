
# recruitment coefficients

# -------------------------------------------------------------------------

# IPM functions.
source("R/auxiliary_functions.R")
source("R/IPM_functions.R")
library(Rcpp)
sourceCpp("R/IPM_functions.cpp")

##########################
##########################

# Read in datasets. First the labels for the species, then the large dataset, next 
# the climatic data and, finally, the fire data.

sp.list <- read.csv2("data/lista_especies_v2.csv"
                     ,dec=".",sep=";",comment.char="")
newspecies <- unique(sp.list[,3])
newspecies <- newspecies[trim(newspecies)!=""]
tesauro <- data.frame(num.sp = c(1:length(newspecies)),
                      name.sp = newspecies)

saplings.IFN2 <- read.table(file = "data/SAPLINGS_IFN2_v8.csv",header=T,sep=";",dec=".")
saplings.IFN3 <- read.table(file = "data/SAPLINGS_IFN3_v10.csv",header=T,sep=";",dec=".")

adult.trees <- read.table(file = "data/PIES_MAYORES_IFN2_IFN3_v1.csv",header=T,sep=";",dec=".")

map <- read.table(file = "data/MAP_v6.csv",header=T,sep=";",dec=".")

adult.trees$name.sp <- tesauro[match(adult.trees$num.sp,tesauro[,1]),2]

##############################
##############################

clima <- read.table(file = "data/clima/clima_referencia_v5.csv",header=T,sep=";",dec=".")
fire <- read.table(file = "data/fire_damage_v2.csv",header=T,sep=";",dec=".")

fire$damage <- pmax(fire$DANYSEVER,fire$DANYMITJA,fire$DANYPETIT)
fire <- subset(fire,DANYSEVER > 0 | DANYMITJA > 0 | DANYPETIT > 0)

######################
######################
# filter data

# We select surviving trees in plots without management. We also eliminate surviving trees which show no measurable growth
# in 10 years: we assume that to be a mistake.

adult.trees <- subset(adult.trees, idestatus != "n" & gestion_ifn3 == "N") 

j <- match(adult.trees$ID,fire$ID,nomatch=0)
adult.trees <- adult.trees[-j,]
adult.trees$abasal <- adult.trees$FACTOR_IFN2*pi*(adult.trees$dbh_ifn2/200)^2

# saplings

saplings.IFN2 <- subset(saplings.IFN2, num.sp != -1)
saplings.IFN3 <- subset(saplings.IFN3, num.sp != -1)

index <- match(saplings.IFN2$ID,fire$ID,nomatch=0)
saplings.IFN2 <- saplings.IFN2[-j,]
saplings.IFN2$name.sp <- tesauro[match(saplings.IFN2$num.sp,tesauro[,1]),2]

index <- match(saplings.IFN3$ID,fire$ID,nomatch=0)
saplings.IFN3 <- saplings.IFN3[-j,]
saplings.IFN3$name.sp <- tesauro[match(saplings.IFN3$num.sp,tesauro[,1]),2]

names(saplings.IFN2)[3] <- "saplings.IFN2"
names(saplings.IFN2)[4] <- "seedlings.IFN2"
names(saplings.IFN3)[3] <- "saplings.IFN3"
names(saplings.IFN3)[4] <- "seedlings.IFN3"

saplings.both <- merge(x=saplings.IFN2,y=saplings.IFN3,by=c("ID","num.sp"))
saplings.both <- saplings.both[,-8]
names(saplings.both)[5] <- "name.sp"

included.ifn2 <- match(paste(as.character(saplings.both$ID),as.character(saplings.both$num.sp),sep=""),
                       paste(as.character(saplings.IFN2$ID),as.character(saplings.IFN2$num.sp),sep=""))

saplings.IFN2 <- saplings.IFN2[-included.ifn2,]

included.ifn3 <- match(paste(as.character(saplings.both$ID),as.character(saplings.both$num.sp),sep=""),
                       paste(as.character(saplings.IFN3$ID),as.character(saplings.IFN3$num.sp),sep=""))

saplings.IFN3 <- saplings.IFN3[-included.ifn3,]

saplings.IFN2$saplings.IFN3 <- 0
saplings.IFN2$seedlings.IFN3 <- 0

saplings.IFN3$saplings.IFN2 <- 0
saplings.IFN3$seedlings.IFN2 <- 0

saplings.IFN2 <- saplings.IFN2[,c("ID","num.sp","name.sp","seedlings.IFN2","saplings.IFN2","seedlings.IFN3","saplings.IFN3")]
saplings.IFN3 <- saplings.IFN3[,c("ID","num.sp","name.sp","seedlings.IFN2","saplings.IFN2","seedlings.IFN3","saplings.IFN3")]
saplings.both <- saplings.both[,c("ID","num.sp","name.sp","seedlings.IFN2","saplings.IFN2","seedlings.IFN3","saplings.IFN3")]

pm <- rbind(saplings.IFN2,saplings.IFN3,saplings.both)
pm <- arrange(pm,ID)

pm <- pm[,c("ID","num.sp","name.sp","saplings.IFN2","saplings.IFN3")]
pm <- subset(pm,num.sp != -1 & !(saplings.IFN2 == 0 & saplings.IFN3 == 0))

pm$temper <- clima$ref_mt[match(pm$ID,clima$ID)]
pm$rain <- clima$ref_pl[match(pm$ID,clima$ID)]
pm$anom.temper <- clima$anom_mt[match(pm$ID,clima$ID)]
pm$anom.rain <- clima$anom_pl[match(pm$ID,clima$ID)]

adult.trees <- arrange(adult.trees,ID)
pm <- arrange(pm,ID)
basal.area.r <- tapply(X=adult.trees$abasal,INDEX=adult.trees$ID,FUN=sum)
pm$basal.area <- basal.area.r[match(pm$ID,names(basal.area.r))]
pm <- pm[which(complete.cases(pm)),]

out <- list(SPECIES = array("",dim=length(newspecies)),
            binom = array(0,dim=length(newspecies)),
            poisson = array(list(),dim=length(newspecies)))

fun <- function(id){
  index <- which(adult.trees.my.sp$ID == id)
  return (ifelse(length(index)>0,
                 sum(adult.trees.my.sp$FACTOR_IFN2[index]*pi*(adult.trees.my.sp$dbh_ifn2[index]/200)^2),
                 0))
}

for (i in 1:length(newspecies)) {
  
  pm.my.sp <- subset(pm,num.sp == i)
  
  index <- which(adult.trees$num.sp == i)
  adult.trees.my.sp <- adult.trees[index,]
  
  pm.my.sp$basal.area.sp <- sapply(X=pm.my.sp$ID,FUN=fun)
  
  # now, a unified glm, independent of saplings in IFN2
  subdata <- pm.my.sp
  
  out$SPECIES[i] <- newspecies[i]
  out$binom[[i]] <- sum(subdata$saplings.IFN3>0)/length(subdata$saplings.IFN3)
  
  j <- which(subdata$saplings.IFN3>0)
  subdata <- subdata[j,]
  subdata$sqrt.sapl <- sqrt(subdata$saplings.IFN3)
  
  #xnam <- c("rain","temper","anom.rain","anom.temper","basal.area","I(rain*anom.temper)", "saplings.IFN2","basal.area.sp")
  xnam <- c("rain","temper","anom.rain","anom.temper","basal.area","I(rain*anom.temper)","basal.area.sp")
  
  naic <- sum(sapply(0:(length(xnam)-1),function(j) choose(length(xnam),j)))
  
  fmla1 <- as.formula(paste("sqrt.sapl ~ ", paste(xnam, collapse= "+")))
  fmla2 <- as.formula(paste("saplings.IFN3 ~ ", paste(xnam, collapse= "+")))
  
  poisson2 <- lm(fmla1,data=subdata)
  start.coef <- coefficients(poisson2)
  if(i != 12){ 
    poisson2 <- glm(fmla2,family=poisson(link="sqrt"),data=subdata,trace=FALSE,maxit=100000,start=start.coef)
    start.coef <- coefficients(poisson2)
    poisson2 <- glm(fmla2,family=poisson(link="identity"),data=subdata,trace=FALSE,maxit=100000,start=start.coef)
  }
  if(i == 12){ # for quercus ilex, standard start.coef fail
    start.coef[1] <- 15
    start.coef[2] <- -0.002
    start.coef[3] <- -0.5
    start.coef[4] <- 2.5
    start.coef[5] <- 2.5
    start.coef[6] <- -0.1
    start.coef[7] <- -0.004
    start.coef[8] <- 0.1
    
    poisson2 <- glm(fmla2,family=poisson(link="sqrt"),data=subdata,trace=FALSE,maxit=100000,start=start.coef)
  }
  out.old <- poisson2
  aic.keep.old <- poisson2$aic
  print(paste(newspecies[i],"... ",sep=""))
  
  for (j in 1:(length(xnam)-1)) {
    var.remove <- combn(length(xnam),j)
    for (k in 1:dim(var.remove)[2]) {
      xnam.aic <- xnam[-var.remove[,k]]
      
      fmla1 <- as.formula(paste("sqrt.sapl ~ ", paste(xnam.aic, collapse= "+")))
      fmla2 <- as.formula(paste("saplings.IFN3 ~ ", paste(xnam.aic, collapse= "+")))
      
      lm.poisson2 <- try(lm(fmla1, data=subdata))
      start.coef <- coefficients(lm.poisson2)
      poisson2 <- try(glm(fmla2,family=poisson(link="sqrt"), data=subdata, model=TRUE, start=start.coef))
      if (attr(poisson2,"class")[1]=="try-error") {
        aic.keep.new <- 1e64
      } else {
        start.coef <- coefficients(poisson2)
        poisson2 <- try(glm(fmla2,family=poisson(link="identity"), data=subdata, model=TRUE, start=start.coef))
        if (attr(poisson2,"class")[1]=="try-error") {
          aic.keep.new <- 1e64
        } else {
          aic.keep.new <- poisson2$aic
          if (aic.keep.new<=aic.keep.old) {
            out.old <- poisson2
            aic.keep.old <- aic.keep.new
          }
        }
      }
    }
    
    out$poisson[[i]] <- out.old
  }
  
}# for species

save(out,file="data/regresiones/recruitment_regression_v15",compress=TRUE)
