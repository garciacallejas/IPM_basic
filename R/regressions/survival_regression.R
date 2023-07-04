
# obtain survival coefficients

source("R/auxiliary_functions.R")
source("R/IPM_functions.R")
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

index <- which((trees$idestatus=="s" | trees$idestatus == "d" | trees$idestatus == "m") & trees$gestion_ifn3 == "N")
trees <- trees[index,]

# Calculates the basal area (in m2) per plot.

trees$abasal <- pi*(trees$dbh_ifn2/200)^2
aabasal <- tapply(trees$FACTOR_IFN2*trees$abasal,trees$ID,function(x) sum(x[!is.na(x)]))
aabasal <- list(IDPLOT=unlist(dimnames(aabasal)),AREA=unlist(aabasal[]))

# Added on 16-7-2012.

stems <- tapply(trees$FACTOR_IFN2,trees$ID,function(x) sum(x[!is.na(x)]))
stems <- list(IDPLOT=unlist(dimnames(stems)),NUM=unlist(stems[]))

# El an?lisis se hace por especie.

  resultados <- list(ESPECIE=array("",dim=length(newspecies)),
                     LOGIT=array(list(),dim=length(newspecies)),
                     LOG.LOGIT=array(list(),dim=length(newspecies)))

  for (i in 1:length(newspecies)) {
  
    print(paste(i," - ",newspecies[i],sep=""))
    resultados$ESPECIE[i] <- newspecies[i]
    
    iaaa <- which(trees$name.sp==newspecies[i])

    if (length(iaaa)>10) {
      my.trees <- trees[iaaa,] 
      
# Lista a utilizar en la regresi?n.
                                                         
      datos <- list()
      datos$vivo <- as.numeric(my.trees$idestatus=="s")
      datos$dbh2 <- my.trees$dbh_ifn2
      datos$lluvia <- clima$ref_pl[match(my.trees$ID,clima$ID)]
      datos$temper <- clima$ref_mt[match(my.trees$ID,clima$ID)]
      datos$anom_lluvia <- clima$anom_pl[match(my.trees$ID,clima$ID)]
      datos$anom_temper <- clima$anom_mt[match(my.trees$ID,clima$ID)]
      #datos$pies <- stems$NUM[match(my.trees$ID,stems$IDPLOT)]
      datos$abasal <- aabasal$AREA[match(my.trees$ID,aabasal$IDPLOT)]

      r.out <- glm(vivo ~ dbh2+
                     lluvia+
                     temper+
                     abasal+
                     anom_lluvia+
                     anom_temper+
                     #I(temper*temper)+
                     #I(abasal*abasal)+
                     I(lluvia*anom_temper),family=binomial(link="logit"), data = datos)
      
      r.out <- stepAIC(r.out,scope=list(lower= ~ 1,upper=~.),steps=1000,trace=FALSE)
      resultados$LOGIT[[i]] <- r.out

      r.out.log <- glm(vivo ~ log(dbh2)+ 
                         lluvia+
                         temper+
                         abasal+
                         anom_lluvia+
                         anom_temper+
                         #I(temper*temper)+
                         #I(abasal*abasal)+
                         I(lluvia*anom_temper),family=binomial(link="logit"), data = datos)
      
      r.out.log <- stepAIC(r.out.log,scope=list(lower= ~ 1,upper=~.),steps=1000,trace=FALSE)
      resultados$LOG.LOGIT[[i]] <- r.out.log    
    }
    
#     temp <- predict(r.out.log,type="response")
#     plot(r.out.log$model[[2]],temp)#plot survival percentage as a function of log(dbh)
#     myplot <- data.frame(vivo = numeric(1),
#                          dbh2 = numeric(1),
#                          lluvia = numeric(1),
#                          temper = numeric(1),
#                          anom_lluvia = numeric(1),
#                          anom_temper = numeric(1),
#                          abasal = numeric(1))
#     myplot$vivo <- 1
#     myplot$dbh2 <- 22.80714
#     myplot$lluvia <- 564.2
#     myplot$temper <- 11.92078
#     myplot$anom_lluvia <- -0.1
#     myplot$anom_temper <- 1.6
#     myplot$abasal <- 14.89532
#     temp2 <- predict(object=r.out.log,newdata=myplot,type="response")
  }
  save(resultados,file="data/regresiones/survival_v8",compress=TRUE)