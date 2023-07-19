rm(list=ls())

set.seed(1)

rainfall <- read.csv("weatherAUS_withcoords.csv",header=T)

library(igraph) 

### rescale coordinates linearly to be within 0 and 1. 
coords <- norm_coords(cbind(rainfall$Latitude,rainfall$Longitude),
                      xmin=0,xmax=1,ymin=0,ymax=1)
head(coords)
class(coords)

### transform the date into a numeric value

dates <- as.Date(rainfall$Date)
head(dates)

min(as.numeric(dates))

# convert to numeric
head(as.numeric(dates))-13817
tail(as.numeric(dates))-13817

t <- as.numeric(dates)-13817
sum((t==1))
dim(coords)

temp <- cbind(rainfall$Date,t)
temp[temp[,1]=="2017-01-01",]

### Transform the data

#plot(rainfall$Rainfall[rainfall$Location==rainfall$Location[1]],type="l")
#plot(diff(rainfall$Rainfall[rainfall$Location==rainfall$Location[1]]),type="l")

cities <- unique(rainfall$Location)
rain_diff <- c()
for(c in 1:49)
{
  rain_temp <- rainfall$Rainfall[rainfall$Location==cities[c]]
  rain_diff <- c(rain_diff,rain_temp[1],diff(rain_temp))
}

dim(rainfall)
length(rain_diff)

### Divide into training and testing dataset

data_tot_temp <- data.frame(xcoords=coords[,1],ycoords=coords[,2],t=t,data=rain_diff)
nrow(data_tot_temp)
which(rainfall$Date=="2017-01-01")
nrow(data_tot_temp[which(rainfall$Date=="2017-01-01"),])
data_tot_temp[which(rainfall$Date=="2017-01-01"),]
data_tot_temp[data_tot_temp$t==3350,]

temp <- rowSums(is.na(data_tot_temp))
data_tot <- data_tot_temp[temp==0,]

head(data_tot_temp[data_tot_temp$t==3350,])
head(data_tot[data_tot$t==3350,])
nrow(data_tot_temp[data_tot_temp$t==3350,])
nrow(data_tot[data_tot$t==3350,])

rm(data_tot_temp)
rm(temp)

data_train <- data_tot[data_tot$t<3350,]
data_test <- data_tot[data_tot$t>=3350,]

nrow(data_train)
nrow(data_test)

ind <- sample(1:nrow(data_train),size=50000,rep=FALSE)
ind <- sort(ind)

data_train <- data_train[ind,]

nrow(data_train)
nrow(data_test)

head(data_train)
tail(data_train)

coords_train <- data_train[,c(1,2)]
class(coords_train)

coords_test <- data_test[,c(1,2)]
class(coords_test)

# X is matrix of regressors
X_train <- as.matrix(rep(1,nrow(data_train)),ncol=1)
X_test <- as.matrix(rep(1,nrow(data_test)),ncol=1)

source("spatiotempo_SB.R")

burn = 10000
runs = 20000

obj_stSB <- STSB(y=data_train$data,x=X_train,z=coords_train,
                 t=data_train$t,DOF=1,
                 mx.sige=1,mx.sigs=1,n.terms=500,pred=T,x_new=X_test,z_new=coords_test,
                 t_new=data_test$t,runs=runs,burn=burn,display=100,
                 kernel="gneiting",a=1)

y_pred_stSB <- apply(obj_stSB$pred[burn:runs,],2,mean)

MSE_stSB <- mean((data_test$data - y_pred_stSB)^2)

save.image("rainfall_stSB.RData")

load("rainfall_stSB.RData")

MSE_stSB 

y_pred_stSB_q025 <- apply(obj_stSB$pred[burn:runs,],2,quantile,probs=0.025)
y_pred_stSB_q975 <- apply(obj_stSB$pred[burn:runs,],2,quantile,probs=0.975)
y_pred_stSB_q500 <- apply(obj_stSB$pred[burn:runs,],2,quantile,probs=0.50)

#plot(data_test$data[data_test$xcoords==data_test$xcoords[1] & data_test$ycoords==data_test$ycoords[1]],col="red")
#lines(y_pred_stSB_q025[data_test$xcoords==data_test$xcoords[1] & data_test$ycoords==data_test$ycoords[1]],col="blue",lty=2)
#lines(y_pred_stSB_q500[data_test$xcoords==data_test$xcoords[1] & data_test$ycoords==data_test$ycoords[1]],col="blue",lty=1)
#lines(y_pred_stSB_q975[data_test$xcoords==data_test$xcoords[1] & data_test$ycoords==data_test$ycoords[1]],col="blue",lty=2)

#table(obj_stSB$betagn)/1001
