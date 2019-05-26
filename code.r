# Libraries 
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
library(dynlm)
library(akima)
library(ardl)


################## DATA CLEANING and Manipulation ######################
dane <- read.csv("data/Shanghai.csv", dec = ".")
dane$date <- NA
dane$date <- paste(dane$year, dane$month, dane$day, dane$hour, sep =",")
x <- "%Y, %m, %d, %H"
dane$date<- as.POSIXct(strptime(dane$date, format = x)) 
dane <-dane %>% select(-c('PM_Jingan', 'PM_Xuhui', 'No', 'cbwd'))
dane <- na.omit(dane)
write.csv(dane, file = "data/Shanghai2.csv",row.names=FALSE)

data2 <- read.csv('data/Shanghai2.csv')
data2 <- xts(data2[,c('PM_US.Post',"DEWP", "HUMI", 'PRES', 'TEMP', 'Iws', 'precipitation', 'Iprec', 
                      'season')], order.by = as.POSIXct(data2[,"date"]), frequency = 24)
write.csv(as.data.frame(data2),"data/data2.csv")

data3 <- ts(data2[,c('PM_US.Post', "HUMI", 'PRES', 'TEMP', 'Iws', 'precipitation', 'Iprec' )], frequency = 24)

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


#### Convert to daily
periodicity(data2)
dailyData <- apply.daily(data2, mean) # Możemy tego używać w miejsce data2
plot(as.zoo(dailyData))
periodicity(dailyData)

#################### SARIMA - KORNEL ######################
#Correlations
round(cor(data2, method = 'pearson'),2)
rcorr(data2, type = 'pearson')
 
#Creating Time Series
m.season <- msts(data2$PM_US.Post,seasonal.period=c(7,30,365.25),start = c(2012,1,1,1), end = c(2018,1,1,1))

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
predict.season <- predict(tbats,h=365)
plot(predict.season, main = "TBATS Forecast", include=365)
print(predict.season)

#Prediction Using Linear Models
time.series.data <- ts(data2$PM_US.Post, start = c(2016,1,1,1), frequency = 7)
mytslm <- tslm(time.series.data ~ trend + season)
print(mytslm)

#Inspecting ARIMA residuals
residuals.arima <- auto.arima(mytslm$residuals)
residuals.Arima.Forecast <- forecast(residuals.arima, h=365)
residualsF <- as.numeric(residuals.Arima.Forecast$mean)

#Inspecting regression residuals
regression.Forecast <- forecast(mytslm,h=365)
regressionF <- as.numeric(regression.Forecast$mean)

#Combined Forecast
forecastR <- regressionF+residualsF
print(forecastR)

#Cross Validation
tbats.predictions.accuracy.list = list()
tbats.predictions.model.stats = list()
tbats.plots = list()
plots.recorder = list()

for (i in 1:10){ 
  nTest <- 14*i  
  nTrain <- length(m.season) - nTest 
  train <- window(m.season,start=decimal_date(as.Date("2012-01-01")),end=c(decimal_date(as.Date("2012-01-01")),nTrain))
  test <- window(m.season, start=c(decimal_date(as.Date("2012-01-01")),nTrain+1), end=c(decimal_date(as.Date("2012-01-01")),nTrain+14))
  
  sample.tbats <- tbats(train)
  sample.predict <- predict(sample.tbats,h=14)
  sample.predict.graphs <- plot(sample.predict, main = "TBATS Forecast", include=14)
  tbats.plots.info <- c(sample.predict.graphs)
  plots.recorder <- recordPlot()
  plot.new()
  
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

for (i in 1:10){
  nTest <- 14*i  
  nTrain <- length(time.series.data)- nTest 
  train <- window(time.series.data,start=decimal_date(as.Date("2012-01-01")),end=c(decimal_date(as.Date("2012-01-01")),nTrain))
  test <- window(time.series.data, start=c(decimal_date(as.Date("2012-01-01")),nTrain+1), end=c(decimal_date(as.Date("2012-01-01")),nTrain+14))
  
  trainlm <- tslm(train ~ trend + season)
  trainlmf <- forecast(trainlm,h=14)
  
  residauto <- auto.arima(trainlm$residuals)
  residf <- forecast(residauto,h=14)
  
  y <- as.numeric(trainlmf$mean)
  x <- as.numeric(residf$mean)
  sp <- x+y
  
  cat("----------------------------------
      
      Data Partition",i,"
      
      Training Set includes",nTrain," time periods. Observations 1 to", nTrain, "
      Test Set includes 14 time periods. Observations", nTrain+1, "to", nTrain+14,"
      
      ")
  print(accuracy(sp,test))
  lmARIMA.predictions.accuracy.list <- c(print(accuracy(sp,test)))
  print(residauto)
  ARIMA.predictions.model.stats <- c(residauto)
  
  cat("
      
      ")
}

print(lmARIMA.predictions.accuracy.list)
print(ARIMA.predictions.model.stats)









######### ARDL - sesonally adjusted monthly data from jDemetra (MICHALINA)

PM <- read.csv('data/adjusted/pm.txt', sep = "\t", dec = ',')
tsPM <- ts(PM[,'Seasonally.adjusted'], start = c(2011,12,1), end = c(2015,12, 1), frequency = 12)

DEWP <- read.csv('data/adjusted/DEWP.txt', sep = "\t", dec = ',')
tsDEWP <- ts(DEWP[,'Seasonally.adjusted'], start = c(2011,12,1), end = c(2015,12, 1), frequency = 12)

HUMI <- read.csv('data/adjusted/HUMI.txt', sep = "\t", dec = ',')
tsHUMI <- ts(HUMI[,'Seasonally.adjusted'], start = c(2011,12,1), end = c(2015,12, 1), frequency = 12)

lprec <- read.csv('data/adjusted/lprec.txt', sep = "\t", dec = ',')
tslprec <- ts(lprec[,'Seasonally.adjusted'], start = c(2011,12,1), end = c(2015,12, 1), frequency = 12)

lws <- read.csv('data/adjusted/lws.txt', sep = "\t", dec = ',')
tslws <- ts(lws[,'Seasonally.adjusted'], start = c(2011,12,1), end = c(2015,12, 1), frequency = 12)

precipitation <- read.csv('data/adjusted/precipitation.txt', sep = "\t", dec = ',')
tsprecipitation <- ts(precipitation[,'Seasonally.adjusted'], start = c(2011,12,1), end = c(2015,12, 1), frequency = 12)

#PRES <- read.csv('PRES.txt', sep = "\t", dec = ',')
#tsPRES <- ts(PRES[,'Seasonally.adjusted'], start = c(2011,12,1), end = c(2015,12, 1), frequency = 12)

#season <- read.csv('PRES.txt', sep = "\t", dec = ',')
#daneseason <- xts(season[,c('Seasonally.adjusted')], order.by = as.Date(season[,'X']), frequency = 12)

TEMP <- read.csv('data/adjusted/TEMP.txt', sep = "\t", dec = ',')
tsTEMP <- ts(TEMP[,'Seasonally.adjusted'], start = c(2011,12,1), end = c(2015,12, 1), frequency = 12)


timeseries <- ts.union(tsPM, tsDEWP, tsHUMI, tslprec, tslws, tsprecipitation, tsTEMP)



##plots
plot(timeseries)
dane.zoo <- as.zoo(timeseries)


##correlation
round(cor(timeseries, method = 'pearson'),2)
round(cor(timeseries, method = 'spearman'),2)

##white noise, h0 - values are independent
Box.test(timeseries[,'tsPM'])
Box.test(timeseries[,'tsDEWP'])
Box.test(timeseries[,'tsHUMI'])
Box.test(timeseries[,'tsTEMP'])
Box.test(timeseries[,'tslws'])
Box.test(timeseries[,'tsprecipitation'])
Box.test(timeseries[,'tslprec'])


#############################
###stationarity
#without differences
adfs1 <- ur.df(timeseries[,c("tsPM")], type = "none", lags = 3)
resids_ <- adfs1@testreg$residuals
bgtest(resids_~1, order = 1)
bgtest(resids_~1, order = 2)
bgtest(resids_~1, order = 3)
bgtest(resids_~1, order = 4)
bgtest(resids_~1, order = 5)
summary(adfs1) #non-stationary

adfs11 <- ur.df(timeseries[,c("tsPM")], type = "drift", lags = 0)
resids_ <- adfs11@testreg$residuals
bgtest(resids_~1, order = 1)
bgtest(resids_~1, order = 2)
bgtest(resids_~1, order = 3)
bgtest(resids_~1, order = 4)
bgtest(resids_~1, order = 5)
summary(adfs11) #stationary
################################

## stationarity on differences
adfs1 <- ur.df(diff(timeseries[,c("tsPM")]), type = "none", lags = 2)
resids_ <- adfs1@testreg$residuals
bgtest(resids_~1, order = 1)
bgtest(resids_~1, order = 2)
bgtest(resids_~1, order = 3)
bgtest(resids_~1, order = 4)
bgtest(resids_~1, order = 5)
summary(adfs1) #stationary

adfs2 <- ur.df(diff(timeseries[,c("tsHUMI")]), type = "none", lags = 1)
resids_ <- adfs2@testreg$residuals
bgtest(resids_~1, order = 1)
bgtest(resids_~1, order = 2)
bgtest(resids_~1, order = 3)
bgtest(resids_~1, order = 4)
bgtest(resids_~1, order = 5)
summary(adfs2) #stationary

adfs3 <- ur.df(diff(timeseries[,c("tsDEWP")]), type = "none", lags = 1)
resids_ <- adfs3@testreg$residuals
bgtest(resids_~1, order = 1)
bgtest(resids_~1, order = 2)
bgtest(resids_~1, order = 3)
bgtest(resids_~1, order = 4)
bgtest(resids_~1, order = 5)
summary(adfs3) #stationary

adfs4 <- ur.df(diff(timeseries[,c("tsTEMP")]), type = "none", lags = 0)
resids_ <- adfs4@testreg$residuals
bgtest(resids_~1, order = 1)
bgtest(resids_~1, order = 2)
bgtest(resids_~1, order = 3)
bgtest(resids_~1, order = 4)
bgtest(resids_~1, order = 5)
summary(adfs4) #stationary

adfs5 <- ur.df(diff(timeseries[,c("tslprec")]), type = "none", lags = 2)
resids_ <- adfs5@testreg$residuals
bgtest(resids_~1, order = 1)
bgtest(resids_~1, order = 2)
bgtest(resids_~1, order = 3)
bgtest(resids_~1, order = 4)
bgtest(resids_~1, order = 5)
summary(adfs5) #stationary

adfs6 <- ur.df(diff(timeseries[,c("tslws")]), type = "none", lags = 2)
resids_ <- adfs6@testreg$residuals
bgtest(resids_~1, order = 1)
bgtest(resids_~1, order = 2)
bgtest(resids_~1, order = 3)
bgtest(resids_~1, order = 4)
bgtest(resids_~1, order = 5)
summary(adfs6) #stationary

adfs7 <- ur.df(diff(timeseries[,c("tsprecipitation")]), type = "none", lags = 2)
resids_ <- adfs7@testreg$residuals
bgtest(resids_~1, order = 1)
bgtest(resids_~1, order = 2)
bgtest(resids_~1, order = 3)
bgtest(resids_~1, order = 4)
bgtest(resids_~1, order = 5)
summary(adfs7) #stationary

###lags
timeseries <- as.zoo(timeseries)
timeseries$diff_PM<-diff(timeseries[,c("tsPM")])
timeseries$diff_HUMI<-diff(timeseries[,c("tsHUMI")])
timeseries$diff_TEMP<-diff(timeseries[,c("tsTEMP")])
timeseries$diff_lprec<-diff(timeseries[,c("tslprec")])
timeseries$diff_lws<-diff(timeseries[,c("tslws")])
timeseries$diff_DEWP<-diff(timeseries[,c("tsDEWP")])
timeseries$diff_precipitation<-diff(timeseries[,c("tsprecipitation")])
#ardl_auto <- auto.ardl(diff_PM ~ diff_HUMI  + diff_TEMP + diff_lprec + diff_lws + diff_DEWP + diff_precipitation, 
#                      data=timeseries[-1,], xmax = c(10,10,10,10,10,10))
#summary(ardl_auto)

ardl_auto2 <- auto.ardl(diff_PM ~ diff_HUMI, 
                        data=timeseries[-1,], xmax = 10, ymax = 10)

summary(ardl_auto2) ##diff_HUMI too many lags, diff_PM 1,2,3 (0.01) & 4,5 (0.1)


ardl_auto3 <- auto.ardl(diff_PM ~ diff_TEMP, 
                        data=timeseries[-1,], xmax = 10, ymax = 10)

summary(ardl_auto3) ##diff_TEMP 4 (0.1), diff_PM 1,2,3 (0.01) 


ardl_auto4 <- auto.ardl(diff_PM ~ diff_DEWP, 
                        data=timeseries[-1,], xmax = 10, ymax = 10)

summary(ardl_auto4) ##diff_DEWP 8 (0.1), diff_PM 2,3 (0.05) & 1 (0.1) 


ardl_auto5 <- auto.ardl(diff_PM ~ diff_lprec, 
                        data=timeseries[-1,], xmax = 10, ymax = 10)

summary(ardl_auto5) ##diff_lprec 2 (0.05), diff_PM 1 (0.01) & 2,3 (0.05)

ardl_auto6 <- auto.ardl(diff_PM ~ diff_lws, 
                        data=timeseries[-1,], xmax = 10, ymax = 10)

summary(ardl_auto6) ##diff_lws 10 (0.1), diff_PM 2,3 (0.01) & 1,4 (0.05) & 10 (0.1)

ardl_auto7 <- auto.ardl(diff_PM ~ diff_precipitation, 
                        data=timeseries[-1,], xmax = 10, ymax = 10)

summary(ardl_auto7) ##diff_precipitation 2 (0.05) & 3 (0.1), diff_PM 1 (0.001) & 2,3 (0.01) & 4 (0.05) & 5 (0.1)

#without lags
ARDL <- dynlm(d(tsPM) ~ L(d(tsPM)) + d(tsDEWP) + L(d(tsDEWP)) + d(tsHUMI) + L(d(tsHUMI)) + d(tslws) + L(d(tslws)) + d(tsprecipitation) + 
                L(d(tsprecipitation)) + d(tslprec) + L(d(tslprec)) + d(tsTEMP) + L(d(tsTEMP)) , data = timeseries,start = c(2011, 12))
summary(ARDL)
AIC(ARDL)
BIC(ARDL)

#without tsHUMI
ARDL2 <- dynlm(d(tsPM) ~ L(d(tsPM)) + d(tsDEWP) + L(d(tsDEWP)) + d(tslws) + L(d(tslws)) + d(tsprecipitation) + 
                 L(d(tsprecipitation)) + d(tslprec) + L(d(tslprec)) + d(tsTEMP) + L(d(tsTEMP)) , data = timeseries,start = c(2011, 12))
summary(ARDL2)
AIC(ARDL2)
BIC(ARDL2)

anova(ARDL, ARDL2)

ARDL3 <- dynlm(d(tsPM) ~ L(d(tsPM), c(1:5)) + d(tsDEWP) + L(d(tsDEWP), 8) + d(tslws) + L(d(tslws),10) + d(tsprecipitation) + 
                 L(d(tsprecipitation),c(2,3)) + d(tslprec) + L(d(tslprec), 2) + d(tsTEMP) + L(d(tsTEMP), 4) , data = timeseries,start = c(2011, 12))
summary(ARDL3)
AIC(ARDL3)
BIC(ARDL3)

anova(ARDL2, ARDL3)

ARDL4 <- dynlm(d(tsPM) ~ L(d(tsPM), c(1:5)) + d(tsDEWP) + L(d(tsDEWP), 8) + d(tslws) + L(d(tslws),10) + d(tsprecipitation) + 
                 L(d(tsprecipitation),2) + d(tslprec) + L(d(tslprec), 2) + d(tsTEMP) + L(d(tsTEMP), 4) , data = timeseries,start = c(2011, 12))
summary(ARDL4)
AIC(ARDL4)
BIC(ARDL4)

anova(ARDL3, ARDL4)

ARDL5 <- dynlm(d(tsPM) ~ L(d(tsPM), c(1:3)) + d(tsDEWP) + L(d(tsDEWP), 8) + d(tslws) + L(d(tslws),10) + d(tsprecipitation) + 
                 L(d(tsprecipitation),2) + d(tslprec) + L(d(tslprec), 2) + d(tsTEMP) + L(d(tsTEMP), 4) , data = timeseries,start = c(2011, 12))
summary(ARDL5)
AIC(ARDL5)
BIC(ARDL5)

anova(ARDL4, ARDL5)

acf(ARDL4$residuals, type='correlation')
#has this model sorted out autocorrelation issue? - Yes!!
bgtest(residuals(ARDL4)~1, order = 1)
bgtest(residuals(ARDL4)~1, order = 2)
bgtest(residuals(ARDL4)~1, order = 3)
bgtest(residuals(ARDL4)~1, order = 4)
bgtest(residuals(ARDL4)~1, order = 5)

acf(ARDL5$residuals, type='correlation')
#has this model sorted out autocorrelation issue? - Yes!!
bgtest(residuals(ARDL5)~1, order = 1)
bgtest(residuals(ARDL5)~1, order = 2)
bgtest(residuals(ARDL5)~1, order = 3)
bgtest(residuals(ARDL5)~1, order = 4)
bgtest(residuals(ARDL5)~1, order = 5)


ARDL6 <- dynlm(d(tsPM)~L(d(tsPM), c(1:3)) + d(tsDEWP) + L(d(tsDEWP), 8) + L(d(tslws), 10)
               + L(d(tsTEMP), 4), data = timeseries, start = c(2011,12))
summary(ARDL6)
AIC(ARDL6)
BIC(ARDL6)
acf(ARDL6$residuals, type='correlation')
bgtest(residuals(ARDL6)~1, order = 1)
bgtest(residuals(ARDL6)~1, order = 2)
bgtest(residuals(ARDL6)~1, order = 3)
bgtest(residuals(ARDL6)~1, order = 4)
bgtest(residuals(ARDL6)~1, order = 5)
##no autocorrelation 

jbTest(as.matrix(residuals(ARDL6))) # residuals normally distributed 


bptest(ARDL6,data=timeseries, studentize=FALSE)
bptest(ARDL6,data=timeseries)  ##homoscedasticity 

resettest(ARDL6)  ##model well fitted/ correctly specified

vif(ARDL6, data = timeseries)

#forecast(ARDL , x , h = 1 , interval = FALSE, level = 0.95 , nSim = 500)
#predict



############## ARDL - MICHAŁ ##############

dane.zoo <- as.zoo(data2)
plot(dane.zoo)

m.season <- msts(data2$PM_US.Post,seasonal.period=c(7,30,365.25),start = c(2012,1,1,1), end = c(2018,1,1,1))
#m.season <- msts(data2$PM_US.Post,seasonal.period=c(24 , 24 * 7, 24 * 30,24 * 365.25))

#tbats <- tbats(m.season, use.box.cox = TRUE, use.trend = TRUE, use.damped.trend = TRUE, seasonal.periods = c(24 , 24 * 7, 24 * 30, 24 * 365.25), use.parallel = TRUE, num.cores = NULL)
#plot(tbats)

#tbats2 <- tbats.components(tbats)
#View(tbats2)
#plot(tbats2)


source("function_testdf.R")

testdf(variable = data2$PM_US.Post, # vector tested
       max.augmentations = 3, max.order=5)  # maximum number of augmentations added

testdf(variable = diff(USA[,c("consumption")]), # vector tested
       max.augmentations = 3,  # maximum number of augmentations added
       max.order=5)           # maximum order of residual lags for BG test

# Na różnicach
ARDL <- dynlm( d(PM_US.Post) ~ 
                 #d(PM_US.Post, 2) +
                 L(d(PM_US.Post), c(1, 24, 24 * 7, 24 * 30, 24 * 365)) +
                 d(DEWP) + 
                 d(DEWP, 24) + 
                 d(DEWP, 24 * 7) + 
                 d(DEWP, 24 * 30) + 
                 d(DEWP, 24 * 365) + 
                 L(d(DEWP), c(1, 24, 24 * 7, 24 * 30, 24 * 365)) +
                 d(HUMI) + 
                 d(HUMI, 24) + 
                 d(HUMI, 24 * 7) + 
                 d(HUMI, 24 * 30) + 
                 d(HUMI, 24 * 365) + 
                 L(d(HUMI), c(1, 24, 24 * 7, 24 * 30, 24 * 365)) +
                 d(PRES) + 
                 d(PRES, 24) + 
                 d(PRES, 24 * 7) + 
                 d(PRES, 24 * 30) + 
                 d(PRES, 24 * 365) + 
                 L(d(PRES), c(1, 24, 24 * 7, 24 * 30, 24 * 365)) +
                 d(TEMP) + 
                 d(TEMP, 24) + 
                 d(TEMP, 24 * 7) + 
                 d(TEMP, 24 * 30) + 
                 d(TEMP, 24 * 365) + 
                 L(d(TEMP), c(1, 24, 24 * 7, 24 * 30, 24 * 365)) +
                 d(Iws) + 
                 d(Iws, 24) + 
                 d(Iws, 24 * 7) + 
                 d(Iws, 24 * 30) + 
                 d(Iws, 24 * 365) + 
                 L(d(Iws), c(1, 24, 24 * 7, 24 * 30, 24 * 365)) +
                 d(precipitation) + 
                 d(precipitation, 24) + 
                 d(precipitation, 24 * 7) + 
                 d(precipitation, 24 * 30) + 
                 d(precipitation, 24 * 365) + 
                 L(d(precipitation), c(1, 24, 24 * 7, 24 * 30, 24 * 365)) +
                 d(Iprec) +
                 d(Iprec, 24) + 
                 d(Iprec, 24 * 7) + 
                 d(Iprec, 24 * 30) + 
                 d(Iprec, 24 * 365) + 
                 L(d(Iprec), c(1, 24, 24 * 7, 24 * 30, 24 * 365)), data = data2)


summary(ARDL)
AIC(ARDL)
BIC(ARDL)

# Na wartościach
ARDL2 <- dynlm( PM_US.Post~ 
                  L(PM_US.Post, c(1, 24, 24 * 7, 24 * 30, 24 * 365)) +
                  DEWP+ 
                  L(DEWP, c(1, 24, 24 * 7, 24 * 30, 24 * 365)) +
                  HUMI+ 
                  L(HUMI, c(1, 24, 24 * 7, 24 * 30, 24 * 365))+
                  PRES+ 
                  L(PRES, c(1, 24, 24 * 7, 24 * 30, 24 * 365))+
                  TEMP+ 
                  L(TEMP, c(1, 24, 24 * 7, 24 * 30, 24 * 365))+
                  Iws+ 
                  L(Iws, c(1, 24, 24 * 7, 24 * 30, 24 * 365))+
                  precipitation+ 
                  L(precipitation, c(1, 24, 24 * 7, 24 * 30, 24 * 365))+
                  Iprec +
                  L(Iprec, c(1, 24, 24 * 7, 24 * 30, 24 * 365)), data = data2)


summary(ARDL2)
AIC(ARDL2)
BIC(ARDL2)















#########################DO NOT INCLUDE - nie wiem do czego to?? ###################################################

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

