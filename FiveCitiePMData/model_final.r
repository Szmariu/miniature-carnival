setwd("C:\\Users\\Michalina\\Desktop\\ekonometria\\FiveCitiePMData")
dane <- read.csv(file.choose(), dec = ".")

library(dplyr)
library(xts)
library(lmtest)
library(tseries)
library(urca)
library(fUnitRoots)
library(systemfit)
library(plm)
library(uroot)
library(forecast)
library(psych)  
library("Hmisc")
library(lubridate)
library(rlist)
library(uroot)

#######SARIMA AND LIEAR TIME SERIES MODELS

#DATA CLEANING and Manipulation

dane$date <- NA
dane$date <- paste(dane$year, dane$month, dane$day, dane$hour, sep =",")
x <- "%Y, %m, %d, %H"
dane$date<- as.POSIXct(strptime(dane$date, format = x)) 
dane <-dane %>% select(-c('PM_Jingan', 'PM_Xuhui', 'No', 'cbwd'))
dane <- na.omit(dane)

write.csv(dane, file = "Shanghai2.csv",row.names=FALSE)

data2 <- read.csv('Shanghai2.csv')

data2 <- xts(data2[,c('PM_US.Post',"DEWP", "HUMI", 'PRES', 'TEMP', 'Iws', 'precipitation', 'Iprec', 
                      'season')], order.by = as.POSIXct(data2[,"date"]), frequency = 24)

write.csv(as.data.frame(data2),"data2.csv")

##PLOTTING to glimpse the data

dane.zoo <- as.zoo(data2)
plot(dane.zoo)

plot(data2$PM_US.Post, type='l', main= "PM2.5")
plot(data2$DEW, type = 'l', main = "DEW")
plot(data2$HUMI, type = 'l', main = "Humidity")
plot(data2$PRES, type = 'l', main = "Pressure")
plot(data2$TEMP, type = 'l', main = "Temperature")
plot(data2$Iws, type = 'l', main = "Iws")
plot(data2$precipitation, type = 'l', main = "precipitation")
plot(data2$Iprec, type = 'l', main = "Iprec")
plot(data2$season, type = 'l', main = "season")

#Correlations
round(cor(data2, method = 'pearson'),2)
rcorr(data2, type = 'pearson')

#Creating Time Series
m.season <- msts(data2$PM_US.Post,seasonal.period=c(7,30,365.25),start = c(2012,1,1), end = c(2015,1,1))
plot(m.season)

m.season2 <- msts(data2$PM_US.Post,seasonal.period=c(7,30,365.25),start = c(2012,1,1), end = c(2015,3,1))
plot(m.season2)

#Seasonality test
hegy.out1 <- hegy.test(m.season, deterministic = c(1,0,1),lag.method = c("AIC"), maxlag = 31)
print(hegy.out1)

#Automatic differentiations needed
nsdiffs(m.season)

#ACF and PACF inspection
acf(m.season)
pacf(m.season)

#Modeling first ARIMA
ms.arima <- auto.arima(m.season, trace = T, seasonal = TRUE, test = 'kpss', ic = 'bic')
summary(ms.arima)

#Decomposition Using TBATS
tbats <- tbats(m.season)
plot(tbats)

#Prediction Using TBATS
predict.season <- predict(tbats,h=60)
plot(predict.season, main = "TBATS Forecast", include=60)
print(predict.season)

#Prediction Using Linear Models
time.series.data <- ts(data2$PM_US.Post, start = c(2012,1,1), frequency = 12)
mytslm <- tslm(time.series.data ~ trend + season)
mytslm2 <- tslm(m.season ~ trend + season)
print(mytslm2)

#Inspecting ARIMA residuals
residuals.arima <- auto.arima(mytslm2$residuals)
residuals.Arima.Forecast <- forecast(residuals.arima, h=60)
residualsF <- as.numeric(residuals.Arima.Forecast$mean)
plot(residuals.Arima.Forecast)

#Inspecting regression residuals
regression.Forecast <- forecast(mytslm2,h=60)
plot(regression.Forecast)
regressionF <- as.numeric(regression.Forecast$mean)

#Combined Forecast
forecastR <- regressionF+residualsF
print(forecastR)
plot(forecastR)
lines(forecastR)

#Cross Validation
tbats.predictions.accuracy.list = list()
tbats.predictions.model.stats = list()
tbats.plots = list()
plots.recorder = list()

for (i in 1:1){ 
  nTest <- 60*i  
  nTrain <- length(m.season) - nTest 
  train <- window(m.season,start=decimal_date(as.Date("2012-01-01")),end=c(decimal_date(as.Date("2012-01-01")),nTrain))
  test <- window(m.season, start=c(decimal_date(as.Date("2012-01-01")),nTrain+1), end=c(decimal_date(as.Date("2012-01-01")),nTrain+60))
  
  sample.tbats <- tbats(train)
  sample.predict <- predict(sample.tbats,h=60)
  sample.predict.graphs <- plot(sample.predict, main = "TBATS Forecast", include=60)
  plot(sample.predict)
  tbats.plots.info <- c(sample.predict.graphs)
  plots.recorder <- recordPlot()
  #plot.new()
  
  cat("----------------------------------
      
      Data Partition",i,"
      
      Training Set includes",nTrain," time periods. Observations 1 to", nTrain, "
      Test Set includes 14 time periods. Observations", nTrain+1, "to", nTrain+14,"
      
      ")
  
  print(accuracy(sample.predict,test))
  tbats.predictions.accuracy.list <- c(accuracy(sample.predict,test))
  
  cat("
      
      ")
  
  print(sample.predict$model)
  tbats.predictions.model.stats <- c(paste0("AIC interation no.",i,": ",sample.predict$model$AIC))

}

print(tbats.predictions.accuracy.list)
print(tbats.predictions.model.stats)
print(tbats.plots)
print(plots.recorder)

#Performing predictions
lmARIMA.predictions.accuracy.list = list()
ARIMA.predictions.model.stats = list()

for (i in 1:1){
  nTest <- 60*i  
  nTrain <- length(m.season)- nTest 
  train <- window(m.season,start=decimal_date(as.Date("2012-01-01")),end=c(decimal_date(as.Date("2012-01-01")),nTrain))
  test <- window(m.season, start=c(decimal_date(as.Date("2012-01-01")),nTrain+1), end=c(decimal_date(as.Date("2012-01-01")),nTrain+60))
  
  trainlm <- tslm(train ~ trend + season)
  trainlmf <- forecast(trainlm,h=60)
  
  residauto <- auto.arima(trainlm$residuals)
  residf <- forecast(residauto,h=60)
  
  y <- as.numeric(trainlmf$mean)
  x <- as.numeric(residf$mean)
  sp <- x+y
  
  plot(sp)
  lines(sp)
  
  cat("----------------------------------
      
      Data Partition",i,"
      
      Training Set includes",nTrain," time periods. Observations 1 to", nTrain, "
      Test Set includes 14 time periods. Observations", nTrain+1, "to", nTrain+14,"
      
      ")
  #print(accuracy(sp,test))
  #lmARIMA.predictions.accuracy.list <- c(print(accuracy(sp,test)))
  #print(residauto)
  #ARIMA.predictions.model.stats <- c(residauto)
  
  cat("
      
      ")
}

print(lmARIMA.predictions.accuracy.list)
print(ARIMA.predictions.model.stats)

#########YOUR SPACE FOR ARDL MODEL









#########################DO NOT INCLUDE###################################################

#Creating Time Series
data3 <- ts(dane$PM_US.Post, frequency=365, start = c(2016,1,1,1), end = c(2018,1,1,1))

#Automatic differentiations needed
nsdiffs(data3)

#fitting best ARIMA
acf(data3)
pacf(data3)

#Modeling first ARIMA
arima <- auto.arima(data3, trace = T, seasonal = TRUE, test = 'kpss', ic = 'bic')

#Decomposing time series with respect to seasonality
decomp <- stl(data3, s.window='periodic')
decomposition <- decompose(data3)
deseasonal_cnt <- seasadj(decomp)

#Plotting results
plot(decomp)
plot(decomposition)
plot(deseasonal_cnt)

#Fitting arima for deseasoned time series
acf(deseasonal_cnt)
pacf(deseasonal_cnt)

arima.decomposed <- auto.arima(deseasonal_cnt, trace = T, seasonal = FALSE, test = 'kpss', ic = 'bic')

#ADF Test for stationarity in deasesoned data
adf.test(deseasonal_cnt, alternative='stationary')

#display Time series for modeled ARIMAs
tsdisplay(residuals(arima), lag.max=32, main='1,1,0 Model Residuals')
tsdisplay(residuals(arima.decomposed), lag.max=32, main='1,1,0 Model Residuals')

#ACF and PACF by ARIMA Residuals
acf(arima$residuals, lag.max = 32)
pacf(arima$residuals, lag.max = 32)

#Tests
jotest=ca.jo(data2, type="trace", K=2, ecdet="none", spec="longrun")
summary(jotest)

source('function_testdf.R')
   
  testdf(variable = diff(arima), # vector tested
         max.augmentations = 5,  # maximum number of augmentations added
         max.order=5)           # maximum order of residual lags for BG test

###################################################################################

