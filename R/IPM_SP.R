rm(list=ls(all=TRUE))
gc()
assign("last.warning", NULL, envir = baseenv())

# IPM functions.
setwd(wd)
source("auxiliary_functions_v9.R")
source("IPM_functions_v20_old.R")
library(Rcpp)
sourceCpp("/home/david/CREAF/articulo/IPM/IPM_functions_v21.cpp")

library(ff)
options(fftempdir = "/home/david/CREAF/articulo/IPM/")

# for loading the whole environment, and the adult.trees matrix:
# 1st load the environment
# 2nd rm(adult.trees) - because we need to load it with ffload, otherwise it will not work
# 3rd ffload("trees_array_original_v34-CAT_R")

###########################################
# Read data

scenario <- "R"

# do we want just tree distribution and abundance or also DBH and demographic dynamics
# note that this can potentially consume a lot of memory

ALL.RESULTS <- F

# final year of simulation. 
# NOTE: if > 2090, it is not compatible with ALL.RESULTS
# NOTE: > 2090 only valid for the CONSTANT CLIMATE SCENARIO (R)

final.year <- 2090

##### KEEP TRACK OF THE SIMULATION NUMBER
##### UPDATE WITH EACH NEW SIMULATION

SIM_NUMBER <- "v39-SP"
sp.tracked <- 1:16

# what do we want to store?
store.original.dbh.dist <- T
store.original.abundances <- T
store.decadal.dbh.dist <- F
store.decadal.abundances <- T

#### colonization parameters
####
max.dist <- 1500

colonization.threshold <- 0.05

# number of saplings to colonize a given plot
# new.saplings <- 100

new.saplings <- c(297.9471, 
                  353.4785,
                  264.4421,
                  339.7690,
                  485.3332,
                  324.3320,
                  485.3627,
                  322.2076,
                  322.2275,
                  254.6479,
                  279.6451,
                  516.2391,
                  605.3193,
                  309.4185,
                  242.8038,
                  502.7275)

load.workspace <- F

if(!load.workspace){
  
  BA_threshold <- Readcsv("BASAL_AREA_THRESHOLD_v3")
  #BA_threshold$perc_95[which(is.na(BA_threshold$perc_95))] <- 15
  BA_threshold$perc_95 <- 15
  
  ##### do we want a log file?
  
  STORE_LOG <- F
  if(STORE_LOG) sink(paste("./resultados/log_",SIM_NUMBER,".txt",sep=""))
  
  ##### do we want colonization?
  
  COLONIZATION <- T
  
  ##### do we want basal area threshold?
  
  BASAL.AREA.THRESHOLD <- T
  
  map <- read.table("/home/david/CREAF/articulo/IPM/MAP_v6.csv",header=T,sep=";",dec=".")

  # bcn.map <- read.table(file="/home/david/CREAF/articulo/IPM/MAP_BCN_v1.csv",header=T,sep=";",dec=".")
  #map <- read.table(file="/home/david/CREAF/articulo/IPM/MAP_CAT_v1.csv",header=T,sep=",",dec=".")
  
  map <- subset(map, valid)
  # map <- subset(map, ID %in% bcn.map$ID)
  # shapefile for plotting
  spain <- readOGR(dsn=".",layer="inland.spain_UTM")
  #spain <- readOGR(dsn=".",layer="CAT_UTM_v2")
  # spain <- readOGR(dsn=".",layer="bcn_UTM_v1")
  spain.df <- fortify(spain)
  ###########################################
  
  ### if the files are already sorted by species, i.e. with a correct "num.sp" field,
  ### we do not need to load the files in any special manner
  
  trees.sorted <- TRUE
  
  if(trees.sorted){
    sp.list <- Readcsv("lista_especies_v2")
    newspecies <- levels(unique(sp.list[,3]))
    newspecies <- newspecies[trim(newspecies)!=""]
    tesauro <- data.frame(num.sp = c(1:length(newspecies)),name.sp = as.character(levels(sp.list[,3])))
    trees.orig <- Readcsv("PIES_MAYORES_IFN3_v10")
    
    # plots recorded in the map
    trees.orig <- trees.orig[which(trees.orig$ID %in% map$ID),]
    
    #   # invalid plots: those with presence of plantation, riverside species.
    #   # we do not model these species, and we exclude plots where they are present.
    #   
    #   invalid.id <- unique(append(map$ID[which(map$LAND_USE %in% c(1,2,5))],unique(trees.orig$ID[which(trees.orig$num.sp == -1)])))
    #   
    #   trees.orig <- trees.orig[which(!(trees.orig$ID %in% invalid.id)),]
    
    saplings.orig <- Readcsv("SAPLINGS_IFN3_v10")
    saplings.orig <- saplings.orig[which(saplings.orig$num.sp %in% tesauro[,1]),]
    
    saplings.orig <- saplings.orig[which(saplings.orig$ID %in% map$ID),]
    #   saplings.orig <- saplings.orig[which(!(saplings.orig$ID %in% invalid.id)),]
    
    # we only need the plots with presence of saplings. In the original files,
    # plots with seedlings but no saplings are stored.
    
    saplings.orig <- saplings.orig[which(saplings.orig$saplings > 0),]
    
    NUM_SP <- length(newspecies)
    
    ########### ONLY TREES OF SP 5
    #trees.orig <- trees.orig[trees.orig$num.sp == MY.SP,]
    #saplings.orig <- saplings.orig[saplings.orig$num.sp == MY.SP,]
    ###########
    
  }else{
    sp.list <- Readcsv("lista_especies_v2")
    newspecies <- levels(unique(sp.list[,3]))
    newspecies <- newspecies[trim(newspecies)!=""]
    tesauro <- data.frame(num.sp = c(1:length(newspecies)),name.sp = as.character(levels(sp.list[,3])))
    NUM_SP <- length(newspecies)
    trees.orig <- load.trees.IFN(3,sp.list)
    saplings.orig <- load.saplings.IFN(3,sp.list)
  }
  
  NUM_PLOTS <- nrow(map)
  
  ###################### some constants
  
  year <- 2000
  
  # IPM Parameters.
  n.time.intervals <- 20
  t.diff <- 10
  initialize.num <- vector("list",NUM_SP)
  
  min.DBH <- 7.5  # in cm.
  n.intervals.mesh <- 1500
  if (n.intervals.mesh<7) stop("Mesh size must be larger than 8")
  
  #x(cm) is the interval of dbh for every sp.
  x <- x.per.species(min.dbh=min.DBH,n.intervals=n.intervals.mesh)
  y <- x
  x2 <- (x/200)^2
  max.diam <- sapply(1:NUM_SP, function(i) max(x[,i]))
  h <- x[2,]-x[1,]
  nx <- n.intervals.mesh+1
  
  # This matrix is used inside the growth function to save CPU time
  # in calculating the quadrature.
  #### A FASTER WAY?!?!?!?!?!?
  y.minus.x <- array(0,dim=c(dim(x)[1],dim(y)[1],NUM_SP))
  for (k in 1:NUM_SP) {
    for (i in 1:dim(x)[1]) for (j in 1:dim(y)[1]) y.minus.x[i,j,k] <- y[j,k] - x[i,k]
  }
  y.minus.x[y.minus.x<0] <- 0
  
  MAXDBH <- numeric(16)
  
  MAXDBH[1] <- 155.54    # conifers.
  MAXDBH[2] <- 290.73    # Deciduous.
  MAXDBH[3] <- 330   # Fagus sylvatica.
  MAXDBH[4] <- 220   # Juniperus thurifera.
  MAXDBH[5] <- 160.6   # Pinus halepensis.
  MAXDBH[6] <- 162.8   # Pinus nigra.
  MAXDBH[7] <- 198   # Pinus pinaster.
  MAXDBH[8] <- 150.7   # Pinus pinea.
  MAXDBH[9] <- 163.9   # Pinus sylvestris.
  MAXDBH[10] <- 141.9  # Pinus uncinata.
  MAXDBH[11] <- 198  # Quercus faginea.
  MAXDBH[12] <- 167.2  # Quercus ilex.
  MAXDBH[13] <- 189.2  # Quercus pyrenaica.
  MAXDBH[14] <- 313.5  # Quercus robur/petraea.
  MAXDBH[15] <- 161.7  # Quercus suber
  MAXDBH[16] <- 114.84   # Sclerophyllous
  
  ###################### climate
  
  clima <- read.table(file="/home/david/CREAF/articulo/IPM/clima/CLIMA_2000.csv",header=T,sep=";",dec=".")
  
  clima.ref <- read.table(file="/home/david/CREAF/articulo/IPM/clima/clima_referencia_v5.csv",header=T,sep=";",dec=".")
  
  # double check...
  
  clima <- clima[which(clima$ID %in% map$ID),]
  clima.ref <- clima.ref[which(clima.ref$ID %in% map$ID),]
  
  my.id <- match(map$ID,clima$ID)
  my.id.ref <- match(map$ID,clima.ref$ID)
  
  ###################### 
  
  # trees: ff 3d array
  # saplings: matrix
  # basal area: matrix
  
  adult.trees <- list()
  for(i.sp in 1:length(newspecies)){
    adult.trees[[newspecies[i.sp]]] <- ff(0,dim = c(NUM_PLOTS,nx))
  }
  ba <- matrix(0,nrow = NUM_PLOTS,ncol = NUM_SP)
  saplings <- matrix(0,nrow = NUM_PLOTS,ncol = NUM_SP)
  
  ################## ALL RESULTS
  ##################
  
  if(ALL.RESULTS){
    
    results <- list()
    for(i in 1:NUM_SP){
      results[[i]] <- data.frame(ID = integer(nrow(map)),
                                 adults_2000 = integer(nrow(map)),
                                 deaths_2000 = integer(nrow(map)),
                                 new_adults_2000 = integer(nrow(map)),
                                 saplings_2000 = integer(nrow(map)),
                                 basal_area_2000 = numeric(nrow(map)),
                                 #
                                 adults_2010 = integer(nrow(map)),
                                 deaths_2010 = integer(nrow(map)),
                                 new_adults_2010 = integer(nrow(map)),
                                 saplings_2010 = integer(nrow(map)),
                                 basal_area_2010 = numeric(nrow(map)),
                                 #
                                 adults_2020 = integer(nrow(map)),
                                 deaths_2020 = integer(nrow(map)),
                                 new_adults_2020 = integer(nrow(map)),
                                 saplings_2020 = integer(nrow(map)),
                                 basal_area_2020 = numeric(nrow(map)),
                                 #
                                 adults_2030 = integer(nrow(map)),
                                 deaths_2030 = integer(nrow(map)),
                                 new_adults_2030 = integer(nrow(map)),
                                 saplings_2030 = integer(nrow(map)),
                                 basal_area_2030 = numeric(nrow(map)),
                                 #
                                 adults_2040 = integer(nrow(map)),
                                 deaths_2040 = integer(nrow(map)),
                                 new_adults_2040 = integer(nrow(map)),
                                 saplings_2040 = integer(nrow(map)),
                                 basal_area_2040 = numeric(nrow(map)),
                                 #
                                 adults_2050 = integer(nrow(map)),
                                 deaths_2050 = integer(nrow(map)),
                                 new_adults_2050 = integer(nrow(map)),
                                 saplings_2050 = integer(nrow(map)),
                                 basal_area_2050 = numeric(nrow(map)),
                                 #
                                 adults_2060 = integer(nrow(map)),
                                 deaths_2060 = integer(nrow(map)),
                                 new_adults_2060 = integer(nrow(map)),
                                 saplings_2060 = integer(nrow(map)),
                                 basal_area_2060 = numeric(nrow(map)),
                                 #
                                 adults_2070 = integer(nrow(map)),
                                 deaths_2070 = integer(nrow(map)),
                                 new_adults_2070 = integer(nrow(map)),
                                 saplings_2070 = integer(nrow(map)),
                                 basal_area_2070 = numeric(nrow(map)),
                                 #
                                 adults_2080 = integer(nrow(map)),
                                 deaths_2080 = integer(nrow(map)),
                                 new_adults_2080 = integer(nrow(map)),
                                 saplings_2080 = integer(nrow(map)),
                                 basal_area_2080 = numeric(nrow(map)),
                                 #  
                                 adults_2090 = integer(nrow(map)),
                                 deaths_2090 = integer(nrow(map)),
                                 new_adults_2090 = integer(nrow(map)),
                                 saplings_2090 = integer(nrow(map)),
                                 basal_area_2090 = numeric(nrow(map)),
                                 #
                                 temperature = numeric(nrow(map)),
                                 precipitation = numeric(nrow(map)))
      
      results[[i]]$ID <- map$ID
    }
  }# if ALL.RESULTS
  ###################### IPM functions coefficients
  
  survival.coef <- load.survival("survival_v8",newspecies)
  growth.coef <- load.growth("growth_v6",newspecies)
  ingrowth.coef <- load.ingrowth("ingrowth_v10",newspecies)
  saplings.coef <- load.saplings("recruitment_regression_v15",newspecies)
  
  param.survival1 <- survival.coef$log.dbh
  param.growth1 <- growth.coef$log.dbh
  param.growth3 <- growth.coef$intercept.variance
  param.growth4 <- growth.coef$slope.variance
  param.sapl1 <- saplings.coef$binom
  param.ingrowth1 <- ingrowth.coef$lambda
  
  load("regresiones/colonization_v3")
  conifers <- c(1,4,5,6,7,8,9,10)
  deciduous <- c(2,3)
  quercus <- c(11,12,13,14,15,16)
  
  ###################### other data structures
  
  #load("distancias_bcn_v1")
  load("distancias_hasta_2236_v2")
  
  ###############################################################
  ###############################################################
  
  ### first, fill up the trees and saplings distributions
  # for that, we infer a pdf from the discrete original data
  
  print(paste(date()," - ",SIM_NUMBER," - ",scenario," - converting IFN data to continuous pdf...",sep=""))
  
  aux.data <- data.frame(ID = sort(rep(map$ID,NUM_SP)),
                         num.sp = rep(1:NUM_SP,NUM_PLOTS),
                         num.plot = sort(rep(1:NUM_PLOTS,NUM_SP)), # just an auxiliary index
                         predicted.basal.area = 0)
  
  trees <- trees.orig %>% group_by(ID,num.sp) %>% summarize(adult.trees = sum(factor))
  aux.data <- merge(aux.data,trees,by.x = c("ID","num.sp"),by.y = c("ID","num.sp"), all.x = T)
  aux.data$adult.trees[is.na(aux.data$adult.trees)] <- 0
  
  aux.data <- merge(aux.data,saplings.orig,by.x = c("ID","num.sp"),by.y = c("ID","num.sp"), all.x = T)
  aux.data$saplings[is.na(aux.data$saplings)] <- 0
  aux.data <- aux.data[,c("ID","num.sp","num.plot","adult.trees","saplings","predicted.basal.area")]
  
  dbh.my.sp <- list()
  for(i.sp in 1:NUM_SP){
    dbh.my.sp[[i.sp]] <- 7.5 + c(0:n.intervals.mesh)*(MAXDBH[i.sp]-7.5)/n.intervals.mesh
  }
  
  for(i in 1:nrow(aux.data)){
    if(aux.data$adult.trees[i] > 0){
      sp <- aux.data$num.sp[i]
      aux.index <- which(trees.orig$num.sp == aux.data$num.sp[i] & trees.orig$ID == aux.data$ID[i])
      orig.data <- rep(trees.orig$dbh[aux.index],round(trees.orig$factor[aux.index]))
      adult.trees[[sp]][aux.data$num.plot[i],] <- EstimateKernel(orig.data = orig.data,
                                                                 sp = sp,
                                                                 nx = nx,
                                                                 h = h,
                                                                 range.x = c(min(dbh.my.sp[[sp]]),max(dbh.my.sp[[sp]])),
                                                                 bandwidth.range = c(0.03,0.5))
      aux.data$predicted.basal.area[i] <- quadTrapezCpp_1(adult.trees[[sp]][aux.data$num.plot[i],]*x2[,sp],h[sp],nx)*pi
      
    }else{adult.trees[[aux.data$num.sp[i]]][aux.data$num.plot[i],] <- 0}
    if (round(i/10000)*10000==i) print(paste(date(),"- kernel smoothing: running loop",i,"of",nrow(aux.data),sep=" "))
    
  }# for i
  
  plot.predicted.ba <- aux.data %>% group_by(ID) %>% summarize(ba = sum(predicted.basal.area))
  
  for(i in 1:nrow(aux.data)){
    if(aux.data$predicted.basal.area[i]>0 & aux.data$saplings[i]>0){
      
      my.plot <- aux.data$num.plot[i]
      my.sp <- aux.data$num.sp[i]
      
      q <- c(clima$precip[my.plot],clima$temp[my.plot],plot.predicted.ba$ba[my.plot],clima$anom_pl[my.plot],clima$anom_temp[my.plot])
      
      param.sapl2 <- CoefTimesVarCpp(param=saplings.coef,
                                     z=q,
                                     # s=aux.data$saplings[aux.data$num.plot == my.plot],
                                     spba=aux.data$predicted.basal.area[aux.data$num.plot == my.plot],
                                     type="recruitment")
      
      saplings[my.plot,my.sp] <- ifelse(IPMSaplingsCpp(aux.data$saplings[i],c(param.sapl1[my.sp],param.sapl2[my.sp]))>100,
                                        100,
                                        IPMSaplingsCpp(aux.data$saplings[i],c(param.sapl1[my.sp],param.sapl2[my.sp])))
      if(saplings[my.plot,my.sp] < 0){
        saplings[my.plot,my.sp] <- 0
      }
    }# if basal area and saplings
    #     if (round(i/1000)*1000==i) print(paste(date(),"- saplings: running loop",i,"of",nrow(aux.data),sep=" "))
  }# for plots
  
  # put back ba values
  for(i.plot in 1:NUM_PLOTS){
    ba[i.plot,] <- aux.data$predicted.basal.area[aux.data$num.plot == i.plot]
  }
  
  # delete auxiliary structures
  
  rm(plot.predicted.ba)
  rm(aux.data)
  rm(trees)
  print(paste(date()," - ",SIM_NUMBER," - ",scenario," - converting IFN data to continuous pdf... done!",sep=""))
  
  #########################################################
  #########################################################
  
  if(store.original.dbh.dist){
    save(list=ls(all=T),file=paste(scenario,"_",SIM_NUMBER,".RData",sep=""))
    for(i.sp in 1:length(newspecies)){
      temp <- adult.trees[[i.sp]]
      ffsave(temp,file = paste("./resultados/trees_array_original_sp_",i.sp,"_",SIM_NUMBER,"_",scenario,sep=""))
    }
    rm(temp)  
  }
  ### store original results
  
  if(store.original.abundances){
    SaveAbundance(map = map,
                  trees = adult.trees,
                  saplings = saplings,
                  tesauro = tesauro,
                  year = "IFN3",
                  scenario = scenario,
                  SIM_NUMBER = SIM_NUMBER,
                  plot.distribution = T,
                  plot.abundance = T,
                  plot.richness = F,
                  plot.ba = F,
                  spain.df,
                  path="./resultados/")
  }
  
}else{
  load(file=paste(scenario,"_",SIM_NUMBER,".RData",sep=""))
  adult.trees <- list()
  for(i.sp in 1:16){
    ffload(file = paste("./resultados/trees_array_original_sp_",i.sp,"_",SIM_NUMBER,"_",scenario,sep=""),overwrite = T)
    adult.trees[[i.sp]] <- temp
    rm(temp)
  }
}

#########################################################
#########################################################
# IPM loop

for(iyear in seq(2000,final.year,10)){
  
  print(paste(date(),"- starting IPM - year",iyear,sep=" "))
  
  ### Read climate data if necessary
  if(iyear != 2000){
    if(scenario != "R"){
      if(iyear != 2010){
        if(scenario == "A2"){
          clima <- read.table(file=paste("/home/david/Data/Climate/Spain/averaged_decadal/MEAN_A2_",iyear,".csv",sep=""),header=T,sep=";",dec=".")
        }else{
          clima <- read.table(file=paste("/home/david/Data/Climate/Spain/averaged_decadal/MEAN_B2_",iyear,".csv",sep=""),header=T,sep=";",dec=".")
        }
      }else{
        clima <- read.table(file="/home/david/Data/Climate/Spain/averaged_decadal/MEAN_2010.csv",header=T,sep=";",dec=".")
      }
    }
  }
  my.id <- match(map$ID,clima$ID)
  
  if(ALL.RESULTS){
    # indexes for storing results
    ###################################################################################
    ###################################################################################
    my.new.adults <- which(names(results[[1]]) == paste("new_adults_",iyear,sep=""))
    my.adults <- which(names(results[[1]]) == paste("adults_",iyear,sep=""))
    my.deaths <- which(names(results[[1]]) == paste("deaths_",iyear,sep=""))
    my.saplings <- which(names(results[[1]]) == paste("saplings_",iyear,sep=""))
    my.ba <- which(names(results[[1]]) == paste("basal_area_",iyear,sep=""))
    
    previous.adults <- which(names(results[[1]]) == paste("adults_",(iyear - 10),sep=""))
    ###################################################################################
    ###################################################################################
  }
  
  # IPM loop        
  for (i in 1:NUM_PLOTS) {
    
    if (round(i/1000)*1000==i) print(paste(date(),"- IPM: running loop",i,"of",NUM_PLOTS," - year",iyear,sep=" "))
    
    sapl <- saplings[i,]
    
    # are there adult trees or saplings?
    if (sum(ba[i,])>0 | sum(sapl)>0){
      
      temper <- clima$temp[my.id[i]]
      rain <- clima$precip[my.id[i]]
      anom.temper <- clima$anom_temp[my.id[i]]  
      anom.rain <- clima$anom_pl[my.id[i]]
      
      q <- c(rain,temper,sum(ba[i,]),anom.rain,anom.temper)
      
      ########################
      ########################
      
      # if adult trees, apply IPM
      if(sum(ba[i,])>0){
        
        sum.ntrees <- numeric(length(newspecies))
        
        for(i.sp in 1:length(newspecies)){
          sum.ntrees[i.sp] <- sum(adult.trees[[i.sp]][i,])#apply(X = adult.trees[[i.sp]][i,],MARGIN = 1,FUN = sum) # this is not the number of trees, it's just for knowing if there are trees of that sp
        }
        
        param.survival2 <- CoefTimesVarCpp(param=survival.coef,z=q,type="survival")
        param.growth2 <- CoefTimesVarCpp(param=growth.coef,z=q,type="growth")
        param.sapl2 <- CoefTimesVarCpp(param=saplings.coef,z=q,spba=ba[i,],type="recruitment")
        
        # IPM, continuous case.
        order <- sample(x = 1:NUM_SP,size = NUM_SP,replace = F)
        for (j in order) {
          if (sum.ntrees[j]>0) {    # Are there trees of jth species?
            dummy <- IPMadultSurvivalTimesGrowthCpp(ntrees=adult.trees[[j]][i,],#trees[[i]]$num[[j]],
                                                     x=x[,j],
                                                     y=y[,j],
                                                     param_survival=c(param.survival1[j],param.survival2[j]),
                                                     param_growth=c(param.growth1[j],param.growth2[j],param.growth3[j],param.growth4[j]),
                                                     t_diff=t.diff,
                                                     max_diam=max.diam[j],
                                                     h=h[j],
                                                     nx=nx,
                                                     y_minus_x=y.minus.x[,,j])          
            adult.trees[[j]][i,] <- dummy
            p.sapl <- c(param.sapl1[j],param.sapl2[j])
            
            # careful with the saplings prediction. If not controlled, it can both explode and drop below 0
            saplings[i,j] <- ifelse(IPMSaplingsCpp(sapl[j],p.sapl)>100,100,IPMSaplingsCpp(sapl[j],p.sapl))
            if(saplings[i,j] < 0){
              saplings[i,j] <- 0
            }
            if(ALL.RESULTS){
              ###################################################################################
              ###################################################################################
              
              # adults and saplings are computed the same way in each step, because
              # of the data conversion to pdfs
              
              results[[j]][i,my.adults] <- quadTrapezCpp_1(dummy,h[j],nx)
              results[[j]][i,my.saplings] <- ifelse(is.null(saplings[i,j]),
                                                    0,
                                                    saplings[i,j]*(10000/(pi*25)))
              
              # deaths, on the other hand, are different for the first step
              if(iyear == 2000){
                results[[j]][i,my.deaths] <- ifelse(sum(trees.orig$factor[trees.orig$num.sp == j & trees.orig$ID == map$ID[i]]) - results[[j]]$adults_2000[i] > 0,
                                                    sum(trees.orig$factor[trees.orig$num.sp == j & trees.orig$ID == map$ID[i]]) - results[[j]]$adults_2000[i],
                                                    0)
              }else{
                results[[j]][i,my.deaths] <- ifelse(results[[j]][i,previous.adults] - results[[j]][i,my.adults] > 0,
                                                    results[[j]][i,previous.adults] - results[[j]][i,my.adults],
                                                    0)
              }
              ###################################################################################
              ###################################################################################
            }
          }#sum(ntrees[[j]])
        }# for j in NUM_SP
      }# if(sum(unlist(trees[[i]]$ba))>0)
      
      ########################
      ########################
      
      # if saplings on the previous step, apply ingrowth IPM
      if(sum(sapl)>0){
        
        param.ingrowth2 <-  CoefTimesVarCpp(param=ingrowth.coef,z=q,s=sapl,type="ingrowth") 
        
        order <- sample(x = 1:NUM_SP,size = NUM_SP,replace = F)
        for(j in order){
          
        if(sapl[j] > 0 & (!BASAL.AREA.THRESHOLD | (sum(ba[i,]) < BA_threshold$perc_95[j]))){

          dummy <- IPMIngrowthIdentityCpp(y[,j],c(param.ingrowth1[j],param.ingrowth2[j]))
          if(ALL.RESULTS){
            ###################################################################################
            results[[j]][i,my.new.adults] <- quadTrapezCpp_1(dummy,h[j],nx)
            ###################################################################################
          }
          adult.trees[[j]][i,] <- adult.trees[[j]][i,] + dummy
          if(ALL.RESULTS){
            ###################################################################################
            results[[j]][i,my.adults] <- quadTrapezCpp_1(adult.trees[[j]][i,],h[j],nx)
            ###################################################################################
          }
          }# if basal area threshold is met for j-th species and there are saplings
        } # for all j-th species
      } # if there are saplings of any sp
      
      ########################
      ########################
      
      # if any of the two conditions is met (saplings or adults), update basal area of the plot
      order <- sample(x = 1:NUM_SP,size = NUM_SP,replace = F)
      for (j in order) {
        dummy <- adult.trees[[j]][i,]
        num.trees <- quadTrapezCpp_1(dummy,h[j],nx)
        if(num.trees < 0.1){
          adult.trees[[j]][i,] <- 0
          ba[i,j] <- 0
        }else{
          ba[i,j] <- quadTrapezCpp_1(dummy*x2[,j],h[j],nx)*pi
        }
        if(ALL.RESULTS){
          ###################################################################################
          results[[j]][i,my.ba] <- ba[i,j] #ifelse(is.null(trees[[i]]$ba[[j]]),0,trees[[i]]$ba[[j]])
          ###################################################################################
        }
      }# for j in NUM_SP
    }# if saplings or adults
    
  }# for i in plots

  #########################################################
  #########################################################
  
  if(COLONIZATION){
    
    # keep a dataframe with the basal area for potential colonization
    temp.suitable <- data.frame(ID = map$ID,
                                saplings = integer(nrow(map)),
                                plot_basal_area = integer(nrow(map)),
                                #neigh = integer(nrow(map)),
                                neigh_ba = integer(nrow(map)),
                                suitable = logical(nrow(map)))
    temp.suitable$plot_basal_area <- rowSums(ba) #sapply(X=trees,FUN=function(x) sum(unlist(x$ba)))
    # default 
    temp.suitable$suitable <- T
    order <- sample(x = 1:NUM_SP,size = NUM_SP,replace = F)
    for(i in order){
      
      # for every species:
      # - select suitable plots
      # - load regression coefficients
      # - predict new number of saplings
      # - append new saplings to the saplings file
      
      # suitable plots are these without presence,
      # closer than 2236 to a plot with presence,
      # with suitable land use,
      # and total basal area lower than the specific threshold
      
      print(paste(date()," - ",SIM_NUMBER," - ",scenario," - year ",iyear," - colonization for sp ",i,"...",sep=""))
      suitable <- GetSuitablePlots(map=map,
                                   temp.suitable=temp.suitable,
                                   adult.trees=adult.trees,
                                   ba = ba,
                                   #trees.index=trees.index,
                                   distancias=distancias,
                                   sp=tesauro[i,1],
                                   max.dist = max.dist,
                                   verbose=T)   
      
      #apply colonization
      if(!is.null(nrow(suitable))){
        if(nrow(suitable)>0){
          if(i %in% conifers) my.model <- colonization.glm[[1]]
          else if(i %in% quercus) my.model <- colonization.glm[[2]]
          else if(i %in% deciduous) my.model <- colonization.glm[[3]]
          
          suitable$colonized <- predict(my.model,newdata = suitable,type = "response")
          suitable$colonized <- ifelse(suitable$colonized > colonization.threshold,1,0)
          
          k <- match(suitable$ID,map$ID)
          #add new saplings to saplings list
          for(j in 1:nrow(suitable)){
            # saplings[k[j],tesauro[i,1]] <- suitable$saplings[j] + saplings[k[j],tesauro[i,1]]
            saplings[k[j],tesauro[i,1]] <- new.saplings[tesauro[i,1]] + saplings[k[j],tesauro[i,1]]
          }
        }
      }# if there are suitable plots
    }#for sp.
  }
  #########################################################
  #########################################################
  print(paste(date()," - ",SIM_NUMBER," - ",scenario," - year ",iyear," - saving results...",sep=""))
  
  ### store results
  if(store.decadal.abundances){
    SaveAbundance(map = map,
                  trees = adult.trees,
                  saplings = saplings,
                  tesauro = tesauro,
                  year = iyear,
                  scenario = scenario,
                  SIM_NUMBER = SIM_NUMBER,
                  plot.distribution = T,
                  plot.abundance = T,
                  plot.richness = F,
                  plot.ba = F,
                  spain.df,
                  path="./resultados/")
  }
  if(store.decadal.dbh.dist){
    save(list=ls(all=T),file=paste(SIM_NUMBER,"_",scenario,"_",iyear,".RData",sep=""))
    for(i.sp in 1:length(newspecies)){
      temp <- adult.trees[[i.sp]]
      ffsave(temp,file = paste("./resultados/trees_array_year_",iyear,"_sp_",i.sp,"_",SIM_NUMBER,"_",scenario,sep=""))
    }
    rm(temp)  
  }
}#for iyear

###############################################################
###############################################################

if(ALL.RESULTS){
  
  setwd("./resultados/images/")
  
  temperature <- clima$temp[match(results[[1]]$ID,clima$ID)]
  precipitation <- clima$precip[match(results[[1]]$ID,clima$ID)]
  
  for(i in sp.tracked){
    
    results[[i]]$temperature <- temperature
    results[[i]]$precipitation <- precipitation
    
    results.my.sp <- results[[i]]
    
    index <- which(results.my.sp$adults_2000 != 0 | results.my.sp$adults_2090 != 0)
    
    if(length(index)>0){
      results.my.sp <- results.my.sp[index,]
      
      NUM_PLOTS <- nrow(results.my.sp)
      
      results2 <- data.frame(year = integer(nrow(results.my.sp) * 10),
                             adults = numeric(nrow(results.my.sp) * 10),
                             new_adults = numeric(nrow(results.my.sp) * 10),
                             deaths = numeric(nrow(results.my.sp) * 10),
                             saplings = numeric(nrow(results.my.sp) * 10),
                             basal.area = numeric(nrow(results.my.sp) * 10),
                             ID = rep(results.my.sp$ID,10))
      
      for(k in 1:nrow(results.my.sp)){
        temp <- which(results2$ID == results.my.sp$ID[k])
        for(j in 1:length(temp)){
          results2$year[temp[j]] <- 1990 + (j*10)
          ###
          adults.column <- which(names(results.my.sp) == paste("adults_",results2$year[temp[j]],sep=""))
          results2$adults[temp[j]] <- results.my.sp[k,adults.column]
          ###
          new.adults.column <- which(names(results.my.sp) == paste("new_adults_",results2$year[temp[j]],sep=""))
          results2$new_adults[temp[j]] <- results.my.sp[k,new.adults.column]
          ###
          deaths.column <- which(names(results.my.sp) == paste("deaths_",results2$year[temp[j]],sep=""))
          results2$deaths[temp[j]] <- results.my.sp[k,deaths.column]
          ###
          saplings.column <- which(names(results.my.sp) == paste("saplings_",results2$year[temp[j]],sep=""))
          results2$saplings[temp[j]] <- results.my.sp[k,saplings.column]
          ###
          ba.column <- which(names(results.my.sp) == paste("basal_area_",results2$year[temp[j]],sep=""))
          results2$basal.area[temp[j]] <- results.my.sp[k,ba.column]
        }
      }
      
      results2$deaths[results2$deaths < 0] <- 0
      write.table(x=results2,file=paste("sp",i,"_dynamics_",SIM_NUMBER,scenario,".csv",sep=""),append=F,sep=";",dec=".")
      #Writecsv(results2,paste("sp",i,"_dynamics_",SIM_NUMBER,scenario,sep=""))
    }
  }
}
# close log file
if(STORE_LOG) sink()
