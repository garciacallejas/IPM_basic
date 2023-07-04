
# obtain colonization coefficients

# -------------------------------------------------------------------------

# IPM functions
source("R/auxiliary_functions.R")
source("R/IPM_functions.R")


###########################################
###########################################

# 1 - build the table for the parameterization

# 1.1 select plots colonized in IFN3

map <- read.table("data/MAP_v6.csv",header=T,sep=";",dec=".")
map <- subset(map,valid)

IFN2.trees <- read.table(file = "data/PIES_MAYORES_IFN2_v8.csv",header = T,sep = ";",dec = ".")
IFN2.trees <- subset(IFN2.trees, ID %in% map$ID & num.sp != -1 & idestatus != "d")
IFN3.trees <- read.table(file = "data/PIES_MAYORES_IFN3_v10.csv",header = T,sep = ";",dec = ".")
IFN3.trees <- subset(IFN3.trees, ID %in% map$ID & num.sp != -1 & idestatus != "d")

IFN2.saplings <- read.table(file = "data/SAPLINGS_IFN2_v8.csv",header = T,sep = ";",dec = ".")
IFN2.saplings <- subset(IFN2.saplings, ID %in% map$ID & num.sp != -1)
IFN3.saplings <- read.table(file = "data/SAPLINGS_IFN3_v10.csv",header = T,sep = ";",dec = ".")
IFN3.saplings <- subset(IFN3.saplings, ID %in% map$ID & num.sp != -1)

load("data/distancias_hasta_2236_v2")

###########################
# my.sp <- 5
max.dist <- 1500
BA_threshold <- data.frame(perc_95 = rep(15,16))
############################

# if group-wise
# create column "group" in 
# - IFN2.trees
# - IFN2.saplings
# - IFN3.saplings
# collapse saplings data by this column (with group_by + summarize)

conifers <- c(1,4,5,6,7,8,9,10)
deciduous <- c(2,3)
quercus <- c(11,12,13,14,15,16)

IFN2.trees$group <- 0
IFN2.trees$group[IFN2.trees$num.sp %in% conifers] <- 1
IFN2.trees$group[IFN2.trees$num.sp %in% quercus] <- 2
IFN2.trees$group[IFN2.trees$num.sp %in% deciduous] <- 3

IFN3.trees$group <- 0
IFN3.trees$group[IFN3.trees$num.sp %in% conifers] <- 1
IFN3.trees$group[IFN3.trees$num.sp %in% quercus] <- 2
IFN3.trees$group[IFN3.trees$num.sp %in% deciduous] <- 3

IFN2.saplings$group <- 0
IFN2.saplings$group[IFN2.saplings$num.sp %in% conifers] <- 1
IFN2.saplings$group[IFN2.saplings$num.sp %in% quercus] <- 2
IFN2.saplings$group[IFN2.saplings$num.sp %in% deciduous] <- 3

IFN2.saplings.grouped <- IFN2.saplings %>% group_by(ID,group) %>% summarize(saplings = sum(saplings))

IFN3.saplings$group <- 0
IFN3.saplings$group[IFN3.saplings$num.sp %in% conifers] <- 1
IFN3.saplings$group[IFN3.saplings$num.sp %in% quercus] <- 2
IFN3.saplings$group[IFN3.saplings$num.sp %in% deciduous] <- 3

IFN3.saplings.grouped <- IFN3.saplings %>% group_by(ID,group) %>% summarize(saplings = sum(saplings))

# get the parameters for the logistic regression group-wise
for(i.group in 1:3){

# presence in IFN2 of saplings or adults
IFN2.presence.ids <- unique(append(IFN2.trees$ID[IFN2.trees$group == i.group],IFN2.saplings.grouped$ID[IFN2.saplings.grouped$group == i.group & IFN2.saplings.grouped$saplings > 0]))
# presence of adults
IFN2.trees.ids <- unique(IFN2.trees$ID[IFN2.trees$group == i.group])

# potentially colonizable plots are those that:
# DO NOT have IFN3 adults of the species (group) -- this would indicate incorrect IFN2 sampling
# do not show IFN2 presence
suitable.ids <- map$ID[!(map$ID %in% unique(IFN3.trees$ID[IFN3.trees$group == i.group])) & !(map$ID %in% IFN2.presence.ids)]

colonization.data <- data.frame(ID = suitable.ids, IFN3.saplings = 0, IFN2.basal.area = 0, IFN2.neigh.ba = 0)

# add the saplings on the potentially colonizable plots
colonized.plots <- IFN3.saplings.grouped[IFN3.saplings.grouped$group == i.group & 
                                   IFN3.saplings.grouped$saplings > 0 &
                                   IFN3.saplings.grouped$ID %in% suitable.ids,]

colonization.data$IFN3.saplings[colonization.data$ID %in% colonized.plots$ID] <- colonized.plots$saplings

# add the total basal area of every plot

temp.ba <- IFN2.trees %>% group_by(group,ID) %>% summarize(basal.area = sum(factor)) 

IFN2.basal.area <- data.frame(ID = map$ID, basal.area = 0, group.ba = 0)
IFN2.basal.area$group.ba[IFN2.basal.area$ID %in% temp.ba$ID[temp.ba$group == i.group]] <- temp.ba$basal.area[temp.ba$group == i.group]

temp.ba.2 <- temp.ba %>% group_by(ID) %>% summarize(basal.area = sum(basal.area))
IFN2.basal.area$basal.area[IFN2.basal.area$ID %in% temp.ba.2$ID] <- temp.ba.2$basal.area

rm(temp.ba, temp.ba.2)

colonization.data$IFN2.basal.area <- IFN2.basal.area$basal.area[IFN2.basal.area$ID %in% colonization.data$ID]

### Trim the suitable plots, select only those within max.distance ###

# 1 - plots with IFN2.presence are IFN2.trees.ids

# 2 - basal area lower than threshold
# in this version, all thresholds are the same
colonization.data <- colonization.data[which(colonization.data$IFN2.basal.area < BA_threshold$perc_95[1]),]

#which have nearby presence of the species
index <- match(colonization.data$ID,distancias[1,])

if(length(index)>1){
  
  if(max.dist < 1414){
    close <- lapply(index,function(x) c(distancias[,x]$d1000))
  }else if(max.dist < 2000){
    close <- lapply(index,function(x) c(distancias[,x]$d1000, distancias[,x]$d1414))
  }else if(max.dist < 2236){
    close <- lapply(index,function(x) c(distancias[,x]$d1000, distancias[,x]$d1414, distancias[,x]$d2000))
  }else{
    close <- lapply(index,function(x) c(distancias[,x]$d1000, distancias[,x]$d1414, distancias[,x]$d2000, distancias[,x]$d2236))
  }
  
  if(class(close) == "list"){
    colonization.data$suitable <- sapply(close,function(x) sum(match(x,IFN2.trees.ids,nomatch=0))>0)
    
    index <- which(colonization.data$suitable == TRUE)
    
    colonization.data <- colonization.data[index,]

  } else {
    colonization.data <- NULL
  }
}# if length(index)>1

# sum up the the specific basal area of the adjacent plots
for(i.plot in 1:nrow(colonization.data)){
  
  close.ids <- UTMDistance2(map[map$ID == colonization.data$ID[i.plot],],map$X_UTM,map$Y_UTM)
  close.ids <- which(close.ids > 0 & close.ids < max.dist)
  
  if(length(close.ids)>0){
    colonization.data$IFN2.neigh.ba[i.plot] <- sum(IFN2.basal.area$group.ba[IFN2.basal.area$ID %in% map$ID[close.ids]])
  }
}

write.table(x = colonization.data,file = paste("data/regresiones/colonization_data_group",i.group,".csv",sep=""), sep=";",dec=".",append = F)
}

########################################
########################################

# validation

# colonization.data <- read.table(file = "data/regresiones/colonization_data_group1.csv",header = T,sep = ";",dec = ".")
# 
# n.folds <- 5
# 
# # training and test sets - 70-30, of both presences and absences of saplings
# training.percentage <- 0.7
# 
# colonization.data$IFN3.saplings.presence <- ifelse(colonization.data$IFN3.saplings > 0,1,0)
# 
# training.presence <- sample(x = colonization.data$ID[colonization.data$IFN3.saplings.presence > 0],size = sum(colonization.data$IFN3.saplings.presence > 0)*training.percentage, replace = F)
# # training.absence <- sample(x = colonization.data$ID[colonization.data$IFN3.saplings.presence == 0],size = length(training.presence), replace = F)
# training.absence <- sample(x = colonization.data$ID[colonization.data$IFN3.saplings.presence == 0],size = sum(colonization.data$IFN3.saplings.presence == 0)*training.percentage, replace = F)
# 
# training.set <- rbind(colonization.data[colonization.data$ID %in% training.presence,],colonization.data[colonization.data$ID %in% training.absence,])
# test.set <- colonization.data[!(colonization.data$ID %in% training.set$ID),] #sample(x = colonization.data$ID[colonization.data$IFN3.saplings.presence > 0],size = sum(colonization.data$IFN3.saplings.presence > 0)*(1 - training.percentage), replace = F)
# 
# # build the logistic regression
# 
# colonization.model <- glm(IFN3.saplings.presence ~ IFN2.basal.area + IFN2.neigh.ba,family=binomial(link='logit'),data=training.set)
# # summary(colonization.model)
# # anova(colonization.model, test="Chisq")
# 
# # try it on the test set
# 
# fitted.results <- predict(colonization.model,newdata=subset(test.set,select=c(3,4)),type='response')
# test.set$predicted.presence <- ifelse(fitted.results > 0.05,1,0)
# 
# # misClasificError <- mean(test.set$predicted.presence != test.set$IFN3.saplings.presence)
# # print(paste('Accuracy',1-misClasificError))
# 
# # store contingency table
# 
# a <- sum(test.set$IFN3.saplings.presence == 1 & test.set$predicted.presence == 1)
# b <- sum(test.set$IFN3.saplings.presence == 0 & test.set$predicted.presence == 1)
# c <- sum(test.set$IFN3.saplings.presence == 1 & test.set$predicted.presence == 0)
# d <- sum(test.set$IFN3.saplings.presence == 0 & test.set$predicted.presence == 0)
# 
# contingency.table <- matrix(c(a, b, c, d), ncol=2, byrow=T)
# # contingency.table
# phi.coef <- phi(t = contingency.table,digits = 4)

########################################
########################################

colonization.glm = list() 

mydata1 <- read.table(file = "data/regresiones/colonization_data_group1.csv",header = T,sep = ";",dec = ".")
mydata2 <- read.table(file = "data/regresiones/colonization_data_group2.csv",header = T,sep = ";",dec = ".")
mydata3 <- read.table(file = "data/regresiones/colonization_data_group3.csv",header = T,sep = ";",dec = ".")

mydata1$IFN3.saplings.presence <- ifelse(mydata1$IFN3.saplings > 0,1,0)
mydata1$plot_basal_area <- mydata1$IFN2.basal.area
mydata1$neigh_ba <- mydata1$IFN2.neigh.ba
mydata1 <- mydata1[,c("IFN3.saplings.presence","plot_basal_area","neigh_ba")]

mydata2$IFN3.saplings.presence <- ifelse(mydata2$IFN3.saplings > 0,1,0)
mydata2$plot_basal_area <- mydata2$IFN2.basal.area
mydata2$neigh_ba <- mydata2$IFN2.neigh.ba
mydata2 <- mydata2[,c("IFN3.saplings.presence","plot_basal_area","neigh_ba")]

mydata3$IFN3.saplings.presence <- ifelse(mydata3$IFN3.saplings > 0,1,0)
mydata3$plot_basal_area <- mydata3$IFN2.basal.area
mydata3$neigh_ba <- mydata3$IFN2.neigh.ba
mydata3 <- mydata3[,c("IFN3.saplings.presence","plot_basal_area","neigh_ba")]

fit1 <- glm(IFN3.saplings.presence ~ plot_basal_area + neigh_ba,family=binomial(link='logit'),data=mydata1)
fit2 <- glm(IFN3.saplings.presence ~ plot_basal_area + neigh_ba,family=binomial(link='logit'),data=mydata2)
fit3 <- glm(IFN3.saplings.presence ~ plot_basal_area + neigh_ba,family=binomial(link='logit'),data=mydata3)

# summary(fit1)
# summary(fit2)
# summary(fit3)

colonization.glm[[1]] <- fit1
colonization.glm[[2]] <- fit2
colonization.glm[[3]] <- fit3

save(colonization.glm,file="data/regresiones/colonization_v3")

########################################
########################################

# mean number of saplings appearing in IFN3 per species

# read IFN2.trees, IFN2.saplings, IFN3.trees, IFN3.saplings

# for each species, get plots with IFN3 saplings and without IFN2 presence

# mean.new.saplings <- numeric(16)
# 
# for(i.sp in 1:16){
#   IFN3.saplings.presence <- IFN3.saplings[!(IFN3.saplings$ID %in% IFN3.trees$ID[IFN3.trees$num.sp == i.sp]) &
#                                          !(IFN3.saplings$ID %in% IFN2.trees$ID[IFN2.trees$num.sp == i.sp]) &
#                                          !(IFN3.saplings$ID %in% IFN2.saplings$ID[IFN2.saplings$num.sp == i.sp & IFN2.saplings$saplings > 0]) & 
#                                          IFN3.saplings$num.sp == i.sp & 
#                                          IFN3.saplings$saplings > 0,]
#   mean.new.saplings[i.sp] <- mean(IFN3.saplings.presence$saplings) #* 10000/(pi*25)
# }
























