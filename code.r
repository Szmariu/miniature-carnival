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
library(ggfortify)
library(magrittr)
library(multipanelfigure)
library(formattable)



################## Data Cleaning ######################
### Load the data
dane <- read.csv("data/Shanghai.csv", dec = ".")

# Create the date variable
dane$date <- NA
dane$date <- paste(dane$year, dane$month, dane$day, dane$hour, sep =",")
x <- "%Y, %m, %d, %H"
dane$date<- as.POSIXct(strptime(dane$date, format = x)) 

# Remove unwanted variables
dane <-dane %>% select(-c('PM_Jingan', 'PM_Xuhui', 'No', 'cbwd'))

# Remove NAs
dane <- na.omit(dane)

# Convert to XTS 
data2 <- xts(data[,c('PM_US.Post',"DEWP", "HUMI", 'PRES', 'TEMP', 'Iws', 'precipitation', 'Iprec', 'season')],
             order.by = as.POSIXct(data2[,"date"]), 
             frequency = 24)

# Write the data to a CSV
data2 %>% as.data.frame() %>% write.csv("data/data2.csv") 





################## Data Manipulation ######################

#### Generate daily data using mean values from each day. Median was also considered.
periodicity(data2)
dailyData <- data2 %>% apply.daily(mean)
dailyData %>% as.zoo() %>% plot()
periodicity(dailyData) # Make sure the daya is daily

# Also generate weekly data
weeklyData <- data2 %>% apply.weekly(mean) 
weeklyData %>% as.zoo() %>% plot()
periodicity(weeklyData)

## Creating Time Series for ARIMA
m.season <- msts(data2$PM_US.Post,seasonal.period=c(7,30,365.25),start = c(2012,1,1,1), end = c(2018,1,1,1))






################## Basic EDA ######################
# All variables
data2 %>% as.zoo() %>% plot()

# Our main variable - PM levels
autoplot(data2$PM_US.Post, ts.colour = 'coral4', xlab = 'Year', ylab = 'Values')

# Other significant variables
x <- autoplot(data2$DEWP, main = "DEW point (C)", ts.colour = 'darkolivegreen3', xlab = 'Year', ylab = 'Values')
y <- autoplot(data2$TEMP, main = "Temperature (C)", ts.colour = 'indianred4', xlab = 'Year', ylab = 'Values')
z <- autoplot(data2$HUMI, main = "Humidity (%)", ts.colour = 'grey34', xlab = 'Year', ylab = 'Values')
t <- autoplot(data2$PRES, main = "Pressure (hPa)", ts.colour = 'goldenrod3', xlab = 'Year', ylab = 'Values')
figure <- multi_panel_figure(columns = 2, rows = 2, panel_label_type = "none")
figure %<>%
  fill_panel(x, column = 1, row = 1) %>%
  fill_panel(y, column = 2, row = 1) %>%
  fill_panel(z, column = 1, row = 2) %>%
  fill_panel(t, column = 2, row = 2)
print(figure)

# The rest ot the variables
autoplot(data2[,c('Iws', 'Iprec', 'precipitation')], xlab = 'Year', ylab = 'Values', ts.colour = 'black')






################## Sesonality removal  ######################

#### Decomposition the hourly PM2.5 Using TBATS
tbats <- tbats(m.season)

# Plot the results
plot(tbats)
tbats %>% 
  tbats.components() %>% 
  autoplot() 


#### Let's try to decompose the DAILY PM variable to remove the sesonality
tbatsDaily <-  dailyData$PM_US.Post %>%
  msts(seasonal.period=c(7,30,365.25)) %>% # Define the sesonal periods 
  tbats() 

# Results look promising
plot(tbatsDaily) 

# Hoever, when plotting on the same graph a flaw with the level is visible
# It's variance when compared to the original vairable is minimal - it might as well be a straing line.
# This means that there is no point to use this decomposed variable in modeling.
# It has to be noted, that this may be a result of the implementation of tbats. 
tbatsDaily %>% 
  tbats.components() %>% 
  autoplot() 





################## First tests ######################

#### Perform Unit Root Test
data2$PM_US.Post %>% 
  ur.pp(type = c("Z-tau"), model = c("trend")) %>%
  summary()
# It looks like the first differences are stationary


#### Fow many differences is needed to make PM2.5 stationary?
# 0 
nsdiffs(m.season)


#### Correlation of all variables
# There is a significant negative correlation between TEMP and PRES
data2 %>% cor(method = 'pearson') %>% round(2)


#### ACF and PACF inspection
acf(m.season)
pacf(m.season)


#### Breusch-Godfrey test for serial correlation of the PM2.5 variable
# Function from classes
source("function_testdf.R")

# For the PM2.5 variable
a <- testdf(variable = data2$PM_US.Post, max.augmentations = 15, max.order = 15) %>% round(2)
View(a)

# For the differenced in the PM2.5 variable
b <- testdf(variable = diff(data2$PM_US.Post), max.augmentations = 15, max.order = 15) %>% round(2)
View(b)




#################### ARIMA ######################
#Modeling first ARIMA
ms.arima <- auto.arima(m.season, trace = T, seasonal = TRUE, test = 'kpss', ic = 'bic')
summary(ms.arima)

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




#################### ARDL on sesonally adjusted monthly data from jDemetra ####################

# Loading the data exported from jDemetra
PM <- read.csv('pm.txt', sep = "\t", dec = ',')
tsPM <- ts(PM[,'Seasonally.adjusted'], start = c(2011,12,1), end = c(2015,12, 1), frequency = 12)

DEWP <- read.csv('DEWP.txt', sep = "\t", dec = ',')
tsDEWP <- ts(DEWP[,'Seasonally.adjusted'], start = c(2011,12,1), end = c(2015,12, 1), frequency = 12)

HUMI <- read.csv('HUMI.txt', sep = "\t", dec = ',')
tsHUMI <- ts(HUMI[,'Seasonally.adjusted'], start = c(2011,12,1), end = c(2015,12, 1), frequency = 12)

lprec <- read.csv('lprec.txt', sep = "\t", dec = ',')
tslprec <- ts(lprec[,'Seasonally.adjusted'], start = c(2011,12,1), end = c(2015,12, 1), frequency = 12)

lws <- read.csv('lws.txt', sep = "\t", dec = ',')
tslws <- ts(lws[,'Seasonally.adjusted'], start = c(2011,12,1), end = c(2015,12, 1), frequency = 12)

PRES <- read.csv('PRES.txt', sep = "\t", dec = ',')
tsPRES <- ts(PRES[,'Seasonally.adjusted'], start = c(2011,12,1), end = c(2015,12, 1), frequency = 12)

TEMP <- read.csv('TEMP.txt', sep = "\t", dec = ',')
tsTEMP <- ts(TEMP[,'Seasonally.adjusted'], start = c(2011,12,1), end = c(2015,12, 1), frequency = 12)

precipitation <- read.csv('precipitation.txt', sep = "\t", dec = ',')
tsprecipitation <- ts(precipitation[,'Seasonally.adjusted'], start = c(2011,12,1), end = c(2015,12, 1), frequency = 12)

timeseries <- ts.union(tsPM, tsDEWP, tsHUMI, tslprec, tslws, tsPRES, tsTEMP, tsprecipitation)

##correlation
round(cor(timeseries, method = 'pearson'),2)
round(cor(timeseries, method = 'spearman'),2)

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

adfs1 <- ur.df(timeseries[,c("tsPM")], type = "drift", lags = 0)
resids_ <- adfs1@testreg$residuals
bgtest(resids_~1, order = 1)
bgtest(resids_~1, order = 2)
bgtest(resids_~1, order = 3)
bgtest(resids_~1, order = 4)
bgtest(resids_~1, order = 5)
summary(adfs1) #stationary
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

adfs7 <- ur.df(diff(timeseries[,c("tsPRES")]), type = "none", lags = 1)
resids_ <- adfs7@testreg$residuals
bgtest(resids_~1, order = 1)
bgtest(resids_~1, order = 2)
bgtest(resids_~1, order = 3)
bgtest(resids_~1, order = 4)
bgtest(resids_~1, order = 5)
summary(adfs7) #stationary

adfs8 <- ur.df(diff(timeseries[,c("tsprecipitation")]), type = "none", lags = 2)
resids_ <- adfs8@testreg$residuals
bgtest(resids_~1, order = 1)
bgtest(resids_~1, order = 2)
bgtest(resids_~1, order = 3)
bgtest(resids_~1, order = 4)
bgtest(resids_~1, order = 5)
summary(adfs8) #stationary

timeseries <- as.zoo(timeseries)

#1 lag
ARDL <- dynlm(d(tsPM) ~ L(d(tsPM)) + d(tsDEWP) + L(d(tsDEWP)) + d(tsHUMI) + L(d(tsHUMI)) + d(tslws) + L(d(tslws)) + d(tsprecipitation) + 
                L(d(tsprecipitation)) + d(tslprec) + L(d(tslprec)) + d(tsTEMP) + L(d(tsTEMP)) + d(tsPRES) + L(d(tsPRES)), data = timeseries,start = c(2011, 12))
summary(ARDL)
AIC(ARDL)
BIC(ARDL)

#1:3 lags
ARDL2 <- dynlm(d(tsPM) ~ L(d(tsPM), c(1:3)) + d(tsDEWP) + L(d(tsDEWP), c(1:3)) + d(tsHUMI) + L(d(tsHUMI), c(1:3)) + d(tslws) + L(d(tslws), c(1:3)) + d(tsprecipitation) + 
                 L(d(tsprecipitation), c(1:3)) + d(tslprec) + L(d(tslprec), c(1:3)) + d(tsTEMP) + L(d(tsTEMP), c(1:3)) + d(tsPRES) + L(d(tsPRES), c(1:3)), data = timeseries,start = c(2011, 12))
summary(ARDL2)
AIC(ARDL2)
BIC(ARDL2)

#1:2 lags
ARDL3 <- dynlm(d(tsPM) ~ L(d(tsPM), c(1:2)) + d(tsDEWP) + L(d(tsDEWP), c(1:2)) + d(tsHUMI) + L(d(tsHUMI), c(1:2)) + d(tslws) + L(d(tslws), c(1:2)) + d(tsprecipitation) + 
                 L(d(tsprecipitation), c(1:2)) + d(tslprec) + L(d(tslprec), c(1:2)) + d(tsTEMP) + L(d(tsTEMP), c(1:2)) + d(tsPRES) + L(d(tsPRES), c(1:2)), data = timeseries,start = c(2011, 12))
summary(ARDL3)
AIC(ARDL3)
BIC(ARDL3)

##lowest AIC & BIC in model with 1:3 lags ARDL2 (AIC = 317.07, BIC = 376.69)

#1:3 lags in PM & 1:2 lags in DEWP
ARDL4 <- dynlm(d(tsPM) ~ L(d(tsPM), c(1,3)) + d(tsDEWP) + L(d(tsDEWP), c(1:2)) + d(tsHUMI) + L(d(tsHUMI), c(1:3)) + d(tslws) + L(d(tslws), c(1:3)) + d(tsprecipitation) + 
                 L(d(tsprecipitation), c(1:3)) + d(tslprec) + L(d(tslprec), c(1:3)) + d(tsTEMP) + L(d(tsTEMP), c(1:3)) + d(tsPRES) + L(d(tsPRES), c(1:3)), data = timeseries,start = c(2011, 12))
summary(ARDL4)
AIC(ARDL4)
BIC(ARDL4)

#1:3 lags in PM & 1 lags in DEWP
ARDL5 <- dynlm(d(tsPM) ~ L(d(tsPM), c(1:3)) + d(tsDEWP) + L(d(tsDEWP)) + d(tsHUMI) + L(d(tsHUMI), c(1:3)) + d(tslws) + L(d(tslws), c(1:3)) + d(tsprecipitation) + 
                 L(d(tsprecipitation), c(1:3)) + d(tslprec) + L(d(tslprec), c(1:3)) + d(tsTEMP) + L(d(tsTEMP), c(1:3)) + d(tsPRES) + L(d(tsPRES), c(1:3)), data = timeseries,start = c(2011, 12))
summary(ARDL5)
AIC(ARDL5)
BIC(ARDL5)

## better ARDL 5 (AIC = 313.8541 & BIC = 369.8607)

#1:2 lags in PM & 1 lag in DEWP
ARDL6 <- dynlm(d(tsPM) ~ L(d(tsPM), c(1:2)) + d(tsDEWP) + L(d(tsDEWP)) + d(tsHUMI) + L(d(tsHUMI), c(1:3)) + d(tslws) + L(d(tslws), c(1:3)) + d(tsprecipitation) + 
                 L(d(tsprecipitation), c(1:3)) + d(tslprec) + L(d(tslprec), c(1:3)) + d(tsTEMP) + L(d(tsTEMP), c(1:3)) + d(tsPRES) + L(d(tsPRES), c(1:3)), data = timeseries,start = c(2011, 12))
summary(ARDL6)
AIC(ARDL6)
BIC(ARDL6)

#1 lags in PM & 1 lag in DEWP
ARDL7 <- dynlm(d(tsPM) ~ L(d(tsPM)) + d(tsDEWP) + L(d(tsDEWP)) + d(tsHUMI) + L(d(tsHUMI), c(1:3)) + d(tslws) + L(d(tslws), c(1:3)) + d(tsprecipitation) + 
                 L(d(tsprecipitation), c(1:3)) + d(tslprec) + L(d(tslprec), c(1:3)) + d(tsTEMP) + L(d(tsTEMP), c(1:3)) + d(tsPRES) + L(d(tsPRES), c(1:3)), data = timeseries,start = c(2011, 12))
summary(ARDL7)
AIC(ARDL7)
BIC(ARDL7)

###best ARDL7 (AIC = 312.6063, BIC = 364.9995)

#1 lags in PM & 1 lag in DEWP & 1:2 lags in HUMI
ARDL8 <- dynlm(d(tsPM) ~ L(d(tsPM)) + d(tsDEWP) + L(d(tsDEWP)) + d(tsHUMI) + L(d(tsHUMI), c(1:2)) + d(tslws) + L(d(tslws), c(1:3)) + d(tsprecipitation) + 
                 L(d(tsprecipitation), c(1:3)) + d(tslprec) + L(d(tslprec), c(1:3)) + d(tsTEMP) + L(d(tsTEMP), c(1:3)) + d(tsPRES) + L(d(tsPRES), c(1:3)), data = timeseries,start = c(2011, 12))
summary(ARDL8)
AIC(ARDL8)
BIC(ARDL8)

#1 lags in PM & 1 lag in DEWP & 1 lag in HUMI
ARDL9 <- dynlm(d(tsPM) ~ L(d(tsPM)) + d(tsDEWP) + L(d(tsDEWP)) + d(tsHUMI) + L(d(tsHUMI)) + d(tslws) + L(d(tslws), c(1:3)) + d(tsprecipitation) + 
                 L(d(tsprecipitation), c(1:3)) + d(tslprec) + L(d(tslprec), c(1:3)) + d(tsTEMP) + L(d(tsTEMP), c(1:3)) + d(tsPRES) + L(d(tsPRES), c(1:3)), data = timeseries,start = c(2011, 12))
summary(ARDL9)
AIC(ARDL9)
BIC(ARDL9)

#1 lags in PM & 1 lag in DEWP & 1:3 lags in HUMI & 1:2 lags in lws 
ARDL10 <- dynlm(d(tsPM) ~ L(d(tsPM)) + d(tsDEWP) + L(d(tsDEWP)) + d(tsHUMI) + L(d(tsHUMI), c(1:3)) + d(tslws) + L(d(tslws), c(1:2)) + d(tsprecipitation) + 
                 L(d(tsprecipitation), c(1:3)) + d(tslprec) + L(d(tslprec), c(1:3)) + d(tsTEMP) + L(d(tsTEMP), c(1:3)) + d(tsPRES) + L(d(tsPRES), c(1:3)), data = timeseries,start = c(2011, 12))
summary(ARDL10)
AIC(ARDL10)
BIC(ARDL10)

#1 lags in PM & 1 lag in DEWP & 1:3 lags in HUMI & 1 lag in lws 
ARDL11 <- dynlm(d(tsPM) ~ L(d(tsPM)) + d(tsDEWP) + L(d(tsDEWP)) + d(tsHUMI) + L(d(tsHUMI), c(1:3)) + d(tslws) + L(d(tslws)) + d(tsprecipitation) + 
                  L(d(tsprecipitation), c(1:3)) + d(tslprec) + L(d(tslprec), c(1:3)) + d(tsTEMP) + L(d(tsTEMP), c(1:3)) + d(tsPRES) + L(d(tsPRES), c(1:3)), data = timeseries,start = c(2011, 12))
summary(ARDL11)
AIC(ARDL11)
BIC(ARDL11)

###best ARDL10 (AIC = 311.0076, BIC = 361.5941)

#1 lags in PM & 1 lag in DEWP & 1:3 lags in HUMI & 1 lag in lws & 1:2 lags in precipitation 
ARDL12 <- dynlm(d(tsPM) ~ L(d(tsPM)) + d(tsDEWP) + L(d(tsDEWP)) + d(tsHUMI) + L(d(tsHUMI), c(1:3)) + d(tslws) + L(d(tslws)) + d(tsprecipitation) + 
                  L(d(tsprecipitation), c(1:2)) + d(tslprec) + L(d(tslprec), c(1:3)) + d(tsTEMP) + L(d(tsTEMP), c(1:3)) + d(tsPRES) + L(d(tsPRES), c(1:3)), data = timeseries,start = c(2011, 12))
summary(ARDL12)
AIC(ARDL12)
BIC(ARDL12)

#1 lags in PM & 1 lag in DEWP & 1:3 lags in HUMI & 1 lag in lws & 1 lags in precipitation 
ARDL13 <- dynlm(d(tsPM) ~ L(d(tsPM)) + d(tsDEWP) + L(d(tsDEWP)) + d(tsHUMI) + L(d(tsHUMI), c(1:3)) + d(tslws) + L(d(tslws)) + d(tsprecipitation) + 
                  L(d(tsprecipitation)) + d(tslprec) + L(d(tslprec), c(1:3)) + d(tsTEMP) + L(d(tsTEMP), c(1:3)) + d(tsPRES) + L(d(tsPRES), c(1:3)), data = timeseries,start = c(2011, 12))
summary(ARDL13)
AIC(ARDL13)
BIC(ARDL13)

###best ARDL10 (AIC = 311.0076, BIC = 361.5941)

#1 lags in PM & 1 lag in DEWP & 1:3 lags in HUMI & 1:2 lags in lws & 1:3 lags in precipitation & 1:2 lags in lprec
ARDL14 <- dynlm(d(tsPM) ~ L(d(tsPM)) + d(tsDEWP) + L(d(tsDEWP)) + d(tsHUMI) + L(d(tsHUMI), c(1:3)) + d(tslws) + L(d(tslws), c(1:2)) + d(tsprecipitation) + 
                  L(d(tsprecipitation), c(1:3)) + d(tslprec) + L(d(tslprec), c(1:2)) + d(tsTEMP) + L(d(tsTEMP), c(1:3)) + d(tsPRES) + L(d(tsPRES), c(1:3)), data = timeseries,start = c(2011, 12))
summary(ARDL14)
AIC(ARDL14)
BIC(ARDL14)

#1 lags in PM & 1 lag in DEWP & 1:3 lags in HUMI & 1:2 lags in lws & 1:3 lags in precipitation & 1 lag in lprec
ARDL15 <- dynlm(d(tsPM) ~ L(d(tsPM)) + d(tsDEWP) + L(d(tsDEWP)) + d(tsHUMI) + L(d(tsHUMI), c(1:3)) + d(tslws) + L(d(tslws), c(1:2)) + d(tsprecipitation) + 
                  L(d(tsprecipitation), c(1:3)) + d(tslprec) + L(d(tslprec)) + d(tsTEMP) + L(d(tsTEMP), c(1:3)) + d(tsPRES) + L(d(tsPRES), c(1:3)), data = timeseries,start = c(2011, 12))
summary(ARDL15)
AIC(ARDL15)
BIC(ARDL15)

###best ARDL15 (AIC = 308.5407, BIC = 355.5139)

#1 lags in PM & 1 lag in DEWP & 1:3 lags in HUMI & 1:2 lags in lws & 1:3 lags in precipitation & 1 lag in lprec & 1:2 lags in TEMP
ARDL16 <- dynlm(d(tsPM) ~ L(d(tsPM)) + d(tsDEWP) + L(d(tsDEWP)) + d(tsHUMI) + L(d(tsHUMI), c(1:3)) + d(tslws) + L(d(tslws), c(1:2)) + d(tsprecipitation) + 
                  L(d(tsprecipitation), c(1:3)) + d(tslprec) + L(d(tslprec)) + d(tsTEMP) + L(d(tsTEMP), c(1:2)) + d(tsPRES) + L(d(tsPRES), c(1:3)), data = timeseries,start = c(2011, 12))
summary(ARDL16)
AIC(ARDL16)
BIC(ARDL16)

#1 lags in PM & 1 lag in DEWP & 1:3 lags in HUMI & 1:2 lags in lws & 1:3 lags in precipitation & 1 lag in lprec & 1 lag in TEMP
ARDL17 <- dynlm(d(tsPM) ~ L(d(tsPM)) + d(tsDEWP) + L(d(tsDEWP)) + d(tsHUMI) + L(d(tsHUMI), c(1:3)) + d(tslws) + L(d(tslws), c(1:2)) + d(tsprecipitation) + 
                  L(d(tsprecipitation), c(1:3)) + d(tslprec) + L(d(tslprec)) + d(tsTEMP) + L(d(tsTEMP)) + d(tsPRES) + L(d(tsPRES), c(1:3)), data = timeseries,start = c(2011, 12))
summary(ARDL17)
AIC(ARDL17)
BIC(ARDL17)

## best ARDL17 (AIC = 307.2816 & BIC = 350.6415)

#1 lags in PM & 1 lag in DEWP & 1:3 lags in HUMI & 1:2 lags in lws & 1:3 lags in precipitation & 1 lag in lprec & 1 lag in TEMP & 1:2 lags in PRES
ARDL18 <- dynlm(d(tsPM) ~ L(d(tsPM)) + d(tsDEWP) + L(d(tsDEWP)) + d(tsHUMI) + L(d(tsHUMI), c(1:3)) + d(tslws) + L(d(tslws), c(1:2)) + d(tsprecipitation) + 
                  L(d(tsprecipitation), c(1:3)) + d(tslprec) + L(d(tslprec)) + d(tsTEMP) + L(d(tsTEMP)) + d(tsPRES) + L(d(tsPRES), c(1:2)), data = timeseries,start = c(2011, 12))
summary(ARDL18)
AIC(ARDL18)
BIC(ARDL18)

#1 lags in PM & 1 lag in DEWP & 1:3 lags in HUMI & 1:2 lags in lws & 1:3 lags in precipitation & 1 lag in lprec & 1 lag in TEMP & 1 lag in PRES
ARDL19 <- dynlm(d(tsPM) ~ L(d(tsPM)) + d(tsDEWP) + L(d(tsDEWP)) + d(tsHUMI) + L(d(tsHUMI), c(1:3)) + d(tslws) + L(d(tslws), c(1:2)) + d(tsprecipitation) + 
                  L(d(tsprecipitation), c(1:3)) + d(tslprec) + L(d(tslprec)) + d(tsTEMP) + L(d(tsTEMP)) + d(tsPRES) + L(d(tsPRES)), data = timeseries,start = c(2011, 12))
summary(ARDL19)
AIC(ARDL19)
BIC(ARDL19)

##best ARDL18 (AIC = 305.763 & BIC = 347.3162)

##ARDL with no lags
ARDL20 <- dynlm(d(tsPM) ~ d(tsDEWP) + d(tsHUMI) + d(tslws) +d(tsprecipitation) + 
                  d(tslprec) + d(tsTEMP) + d(tsPRES), data = timeseries,start = c(2011, 12))
summary(ARDL20)
AIC(ARDL20)
BIC(ARDL20)

##ARDL with no lags in PM  & 1 lag for the rest
ARDL21 <- dynlm(d(tsPM) ~ d(tsDEWP) + L(d(tsDEWP)) + d(tsHUMI) + L(d(tsHUMI)) + d(tslws) + L(d(tslws)) + d(tsprecipitation) + 
                L(d(tsprecipitation)) + d(tslprec) + L(d(tslprec)) + d(tsTEMP) + L(d(tsTEMP)) + d(tsPRES) + L(d(tsPRES)), data = timeseries,start = c(2011, 12))
summary(ARDL21)
AIC(ARDL21)
BIC(ARDL21)

##ARDL with no lags in PM & DEWP & 1 lag for the rest
ARDL22 <- dynlm(d(tsPM) ~ d(tsDEWP) + d(tsHUMI) + L(d(tsHUMI)) + d(tslws) + L(d(tslws)) + d(tsprecipitation) + 
                  L(d(tsprecipitation)) + d(tslprec) + L(d(tslprec)) + d(tsTEMP) + L(d(tsTEMP)) + d(tsPRES) + L(d(tsPRES)), data = timeseries,start = c(2011, 12))
summary(ARDL22)
AIC(ARDL22)
BIC(ARDL22)

##ARDL with no lags in PM & DEWP & HUMI & 1 lag for the rest
ARDL23 <- dynlm(d(tsPM) ~ d(tsDEWP) + d(tsHUMI) + d(tslws) + L(d(tslws)) + d(tsprecipitation) + 
                  L(d(tsprecipitation)) + d(tslprec) + L(d(tslprec)) + d(tsTEMP) + L(d(tsTEMP)) + d(tsPRES) + L(d(tsPRES)), data = timeseries,start = c(2011, 12))
summary(ARDL23)
AIC(ARDL23)
BIC(ARDL23)

##ARDL with no lags in PM & DEWP & HUMI & lws & 1 lag for the rest
ARDL24 <- dynlm(d(tsPM) ~ d(tsDEWP) + d(tsHUMI) + d(tslws) + d(tsprecipitation) + 
                  L(d(tsprecipitation)) + d(tslprec) + L(d(tslprec)) + d(tsTEMP) + L(d(tsTEMP)) + d(tsPRES) + L(d(tsPRES)), data = timeseries,start = c(2011, 12))
summary(ARDL24)
AIC(ARDL24)
BIC(ARDL24)

##ARDL with no lags in PM & DEWP & HUMI & lws & pecipitation & 1 lag for the rest
ARDL25 <- dynlm(d(tsPM) ~ d(tsDEWP) + d(tsHUMI) + d(tslws) + d(tsprecipitation) + 
                  d(tslprec) + L(d(tslprec)) + d(tsTEMP) + L(d(tsTEMP)) + d(tsPRES) + L(d(tsPRES)), data = timeseries,start = c(2011, 12))
summary(ARDL25)
AIC(ARDL25)
BIC(ARDL25)

##ARDL with no lags in PM & DEWP & HUMI & lws & pecipitation & lprec & 1 lag for the rest
ARDL26 <- dynlm(d(tsPM) ~ d(tsDEWP) + d(tsHUMI) + d(tslws) + d(tsprecipitation) + 
                  d(tslprec) + d(tsTEMP) + L(d(tsTEMP)) + d(tsPRES) + L(d(tsPRES)), data = timeseries,start = c(2011, 12))
summary(ARDL26)
AIC(ARDL26)
BIC(ARDL26)

##ARDL with no lags in PM & DEWP & HUMI & lws & pecipitation & lprec & TEMP & 1 lag for the rest
ARDL27 <- dynlm(d(tsPM) ~ d(tsDEWP) + d(tsHUMI) + d(tslws) + d(tsprecipitation) + 
                  d(tslprec) + d(tsTEMP) + d(tsPRES) + L(d(tsPRES)), data = timeseries,start = c(2011, 12))
summary(ARDL27)
AIC(ARDL27)
BIC(ARDL27)

##best ARDL18 (AIC = 305.763 & BIC = 347.3162)
##model statystyczny istotny p-Value < 0.05

#1 lags in PM & 1 lag in DEWP & 1:3 lags in HUMI & 1:2 lags in lws & 1:3 lags in precipitation & 1 lag in lprec & 1 lag in TEMP & 1:2 lags in PRES
ARDL18 <- dynlm(d(tsPM) ~ L(d(tsPM)) + d(tsDEWP) + L(d(tsDEWP)) + d(tsHUMI) + L(d(tsHUMI), c(1:3)) + d(tslws) + L(d(tslws), c(1:2)) + d(tsprecipitation) + 
                  L(d(tsprecipitation), c(1:3)) + d(tslprec) + L(d(tslprec)) + d(tsTEMP) + L(d(tsTEMP)) + d(tsPRES) + L(d(tsPRES), c(1:2)), data = timeseries,start = c(2011, 12))
summary(ARDL18)
AIC(ARDL18)
BIC(ARDL18)

#no lags in PM & 1 lag in DEWP & 1:3 lags in HUMI & 1:2 lags in lws & 1:3 lags in precipitation & 1 lag in lprec & 1 lag in TEMP & 1:2 lags in PRES
ARDL28 <- dynlm(d(tsPM) ~ d(tsDEWP) + L(d(tsDEWP)) + d(tsHUMI) + L(d(tsHUMI), c(1:3)) + d(tslws) + L(d(tslws), c(1:2)) + d(tsprecipitation) + 
                  L(d(tsprecipitation), c(1:3)) + d(tslprec) + L(d(tslprec)) + d(tsTEMP) + L(d(tsTEMP)) + d(tsPRES) + L(d(tsPRES), c(1:2)), data = timeseries,start = c(2011, 12))
summary(ARDL28)
AIC(ARDL28)
BIC(ARDL28)


#no lags in PM & DEWP & 1:3 lags in HUMI & 1:2 lags in lws & 1:3 lags in precipitation & 1 lag in lprec & 1 lag in TEMP & 1:2 lags in PRES
ARDL29 <- dynlm(d(tsPM) ~ d(tsDEWP) + d(tsHUMI) + L(d(tsHUMI), c(1:3)) + d(tslws) + L(d(tslws), c(1:2)) + d(tsprecipitation) + 
                  L(d(tsprecipitation), c(1:3)) + d(tslprec) + L(d(tslprec)) + d(tsTEMP) + L(d(tsTEMP)) + d(tsPRES) + L(d(tsPRES), c(1:2)), data = timeseries,start = c(2011, 12))
summary(ARDL29)
AIC(ARDL29)
BIC(ARDL29)

#no lags in PM & DEWP & HUMI & lws & 1:3 lags in precipitation & 1 lag in lprec & 1 lag in TEMP & 1:2 lags in PRES
ARDL30 <- dynlm(d(tsPM) ~ d(tsDEWP) + d(tsHUMI) + d(tslws) + d(tsprecipitation) + 
                  L(d(tsprecipitation), c(1:3)) + d(tslprec) + L(d(tslprec)) + d(tsTEMP) + L(d(tsTEMP)) + d(tsPRES) + L(d(tsPRES), c(1:2)), data = timeseries,start = c(2011, 12))
summary(ARDL30)
AIC(ARDL30)
BIC(ARDL30)

#no lags in PM & DEWP & HUMI & lws & precipitation & 1 lag in lprec & 1 lag in TEMP & 1:2 lags in PRES
ARDL31 <- dynlm(d(tsPM) ~ d(tsDEWP) + d(tsHUMI) + d(tslws) + d(tsprecipitation) + 
                  d(tslprec) + L(d(tslprec)) + d(tsTEMP) + L(d(tsTEMP)) + d(tsPRES) + L(d(tsPRES), c(1:2)), data = timeseries,start = c(2011, 12))
summary(ARDL31)
AIC(ARDL31)
BIC(ARDL31)

#no lags in PM & DEWP & HUMI & lws & precipitation & lprec & 1 lag in TEMP & 1:2 lags in PRES
ARDL32 <- dynlm(d(tsPM) ~ d(tsDEWP) + d(tsHUMI) + d(tslws) + d(tsprecipitation) + 
                  d(tslprec) + d(tsTEMP) + L(d(tsTEMP)) + d(tsPRES) + L(d(tsPRES), c(1:2)), data = timeseries,start = c(2011, 12))
summary(ARDL32)
AIC(ARDL32)
BIC(ARDL32)

#no lags in PM & DEWP & HUMI & lws & precipitation & lprec & TEMP & 1:2 lags in PRES
ARDL33 <- dynlm(d(tsPM) ~ d(tsDEWP) + d(tsHUMI) + d(tslws) + d(tsprecipitation) + 
                  d(tslprec) + d(tsTEMP) +d(tsPRES) + L(d(tsPRES), c(1:2)), data = timeseries,start = c(2011, 12))
summary(ARDL33)
AIC(ARDL33)
BIC(ARDL33)

##best ARDL18 (AIC = 305.763 & BIC = 347.3162)
##model statystyczny istotny p-Value < 0.05
###AIC i BIC im mniejszy tym lepszy

autoplot(acf(ARDL18$residuals, type='correlation', plot = FALSE))
autoplot(ARDL18)

jbTest(as.matrix(residuals(ARDL18))) # residuals normally distributed 


bptest(ARDL18,data=timeseries, studentize=FALSE)
bptest(ARDL18,data=timeseries)  ##homoscedasticity 

resids_ <- ARDL18$residuals
bgtest(resids_~1, order = 1)
bgtest(resids_~1, order = 2)
bgtest(resids_~1, order = 3)
bgtest(resids_~1, order = 4)
bgtest(resids_~1, order = 5)


#resettest(fitted(ARDL18), data = timeseries)  ##model well fitted/ correctly specified
#resettest(dynlm(d(tsPM)~L(d(tsPM), c(1:3)) + d(tsDEWP) + L(d(tsDEWP), 8) + L(d(tslws), 10)
    #            + L(d(tsTEMP), 4), data = timeseries, start = c(2011,12)))

#vif(ARDL18, data = timeseries)

ARDL18v1 <- dynlm(d(tsPM) ~ L(d(tsPM)) + L(d(tsDEWP)) + L(d(tsHUMI), 3) + d(tslws) + L(d(tsprecipitation), 2) + L(d(tsPRES), c(1:2)), data = timeseries,start = c(2011, 12))

ARDL18v2 <- dynlm(d(tsPM) ~ L(d(tsPM)) + L(d(tsHUMI), 3) + d(tslws) + L(d(tsprecipitation), 2) + L(d(tsPRES), c(1:2)), data = timeseries,start = c(2011, 12))

#jbTest(as.matrix(residuals(ARDL18v2))) # residuals normally distributed 


#bptest(ARDL18v2,data=timeseries, studentize=FALSE)
#bptest(ARDL18v2,data=timeseries)  ##homoscedasticity 

#resids_ <- ARDL18v2$residuals
#bgtest(resids_~1, order = 1)
#bgtest(resids_~1, order = 2)
#bgtest(resids_~1, order = 3)
#bgtest(resids_~1, order = 4)
#bgtest(resids_~1, order = 5)



timeseries.end <- floor(0.7*length(timeseries)/8) 
timeseries.train <- timeseries[(1:timeseries.end),] 
timeseries.test <- timeseries[(timeseries.end+1):(length(timeseries)/8),]

timeseries.test <- as.data.frame(timeseries.test)
timeseries <- as.zoo(timeseries)

ARDL18v3 <- dynlm(d(tsPM) ~ L(d(tsPM)) + L(d(tsHUMI), 3) + d(tslws) + L(d(tsprecipitation), 2) + L(d(tsPRES)) + L(d(tsPRES), 2), data = timeseries,start = c(2011, 12))

pred <- predict(ARDL18v3, newdata = timeseries.test)
t0 <- tail(diff(timeseries$tsPM)*(-1),14)

pred
t0

mean_t0 <- mean(t0)

SSE <- sum((t0 - pred) ^ 2)
SST <- sum((t0 - mean_t0) ^ 2)
r2 <- 1 - SSE/SST

r_squared <- function(vals, preds) {
  1 - (sum((vals*(-1) - preds)^2) / sum((vals - mean(preds))^2))
}
r_squared(t0, pred)



############## ARDL on hourly, daily and weekly data ##############
#### The metodology is as follows
# We start with hourly, daily and weekly data
# And for each one we try to estimate both the value of PM2.5 as well as the differences
# For start we use every variable, lagged in all suspected sesonalities
# Finally, we choose one model in each type of data that performed the best
# And we optimize it by removing insignificant variables
# In the end, we run tests on these models

# Hourly - Differences
dynlm(d(PM_US.Post) ~
        L(d(PM_US.Post), c(1, 24 * 7, 24 * 30, 24 * 365)) +
        d(DEWP) + 
        d(DEWP, 24) + 
        d(DEWP, 24 * 7) + 
        d(DEWP, 24 * 30) + 
        d(DEWP, 24 * 365) + 
        L(d(DEWP), c(1, 24 * 7, 24 * 30, 24 * 365)) +
        d(HUMI) + 
        d(HUMI, 24) + 
        d(HUMI, 24 * 7) + 
        d(HUMI, 24 * 30) + 
        d(HUMI, 24 * 365) + 
        L(d(HUMI), c(1, 24 * 7, 24 * 30, 24 * 365)) +
        d(PRES) + 
        d(PRES, 24) + 
        d(PRES, 24 * 7) + 
        d(PRES, 24 * 30) + 
        d(PRES, 24 * 365) + 
        L(d(PRES), c(1, 24 * 7, 24 * 30, 24 * 365)) +
        d(TEMP) + 
        d(TEMP, 24) + 
        d(TEMP, 24 * 7) + 
        d(TEMP, 24 * 30) + 
        d(TEMP, 24 * 365) + 
        L(d(TEMP), c(1, 24 * 7, 24 * 30, 24 * 365)) +
        d(Iws) + 
        d(Iws, 24) + 
        d(Iws, 24 * 7) + 
        d(Iws, 24 * 30) + 
        d(Iws, 24 * 365) + 
        L(d(Iws), c(1, 24 * 7, 24 * 30, 24 * 365)) +
        d(precipitation) + 
        d(precipitation, 24) + 
        d(precipitation, 24 * 7) + 
        d(precipitation, 24 * 30) + 
        d(precipitation, 24 * 365) + 
        L(d(precipitation), c(1, 24 * 7, 24 * 30, 24 * 365)) +
        d(Iprec) +
        d(Iprec, 24) + 
        d(Iprec, 24 * 7) + 
        d(Iprec, 24 * 30) + 
        d(Iprec, 24 * 365) + 
        L(d(Iprec), c(1, 24 * 7, 24 * 30, 24 * 365)), data = data2) %>% summary() # 0.02

# Hourly - Values
dynlm( PM_US.Post ~ 
          L(PM_US.Post, c(1, 24 * 7, 24 * 30, 24 * 365)) +
          L(DEWP, c(0, 1, 24 * 7, 24 * 30, 24 * 365)) +
          L(HUMI, c(0, 1, 24 * 7, 24 * 30, 24 * 365))+
          L(PRES, c(0, 1, 24 * 7, 24 * 30, 24 * 365))+
          L(TEMP, c(0, 1, 24 * 7, 24 * 30, 24 * 365))+
          L(Iws, c(0, 1, 24 * 7, 24 * 30, 24 * 365))+
          L(precipitation, c(0, 1, 24 * 7, 24 * 30, 24 * 365))+
          L(Iprec, c(0, 1, 24 * 7, 24 * 30, 24 * 365)), data = data2) %>% summary() # 0.88

# Daily - Differences
dynlm( d(PM_US.Post) ~ 
         d(PM_US.Post, 2) +
         L(d(PM_US.Post), c(1, 7, 30, 365)) +
         d(DEWP) + 
         d(DEWP, 7) + 
         d(DEWP, 30) + 
         d(DEWP, 365) + 
         L(d(DEWP), c(1, 7, 30, 365)) +
         d(HUMI) + 
         d(HUMI, 7) + 
         d(HUMI, 30) + 
         d(HUMI, 365) + 
         L(d(HUMI), c(1, 7, 30, 365)) +
         d(PRES) + 
         d(PRES, 7) + 
         d(PRES, 30) + 
         d(PRES, 365) + 
         L(d(PRES), c(1, 7, 30, 365)) +
         d(TEMP) + 
         d(TEMP, 7) + 
         d(TEMP, 30) + 
         d(TEMP, 365) + 
         L(d(TEMP), c(1, 7, 30, 365)) +
         d(Iws) + 
         d(Iws, 7) + 
         d(Iws, 30) + 
         d(Iws, 365) + 
         L(d(Iws), c(1, 7, 30, 365)) +
         d(precipitation) +  
         d(precipitation, 7) + 
         d(precipitation, 30) + 
         d(precipitation, 365) + 
         L(d(precipitation), c(1, 7, 30, 365)) +
         d(Iprec) + 
         d(Iprec, 7) + 
         d(Iprec, 30) + 
         d(Iprec, 365) + 
         L(d(Iprec), c(1, 7, 30, 365)), data = dailyData) %>% summary() # 0.15

# Daily - Values
dynlm( PM_US.Post~ 
          L(PM_US.Post, c(1, 7, 30, 365)) +
          L(DEWP, c(0, 1, 7, 30, 365)) +
          L(HUMI, c(0, 1, 7, 30, 365))+
          L(PRES, c(0, 1, 7, 30, 365))+
          L(TEMP, c(0, 1, 7, 30, 365))+
          L(Iws, c(0, 1, 7, 30, 365))+
          L(precipitation, c(0, 1, 7, 30, 365))+
          L(Iprec, c(0, 1, 7, 30, 365)), data = dailyData) %>% summary() # 0.45


# Daily - Values only from lag(1)
dynlm( PM_US.Post ~ L(PM_US.Post), data = dailyData) %>% summary() # 0.35

# Weekly - Differences
dynlm(d(PM_US.Post) ~
        d(PM_US.Post, 2) +
        L(d(PM_US.Post), c(1, 4, 4 * 12)) +
        d(DEWP, 4) + 
        d(DEWP, 4 * 12) + 
        L(d(DEWP), c(0, 1, 4, 4 * 12)) +
        d(HUMI, 4) + 
        d(HUMI, 4 * 12) + 
        L(d(HUMI), c(0, 1, 4, 4 * 12)) +
        d(PRES, 4) + 
        d(PRES, 4 * 12) +  
        L(d(PRES), c(0, 1, 4, 4 * 12)) +
        d(TEMP, 4) + 
        d(TEMP, 4 * 12) +  
        L(d(TEMP), c(0, 1, 4, 4 * 12)) +
        d(Iws, 4) + 
        d(Iws, 4 * 12) +  
        L(d(Iws), c(0, 1, 4, 4 * 12)) +
        d(precipitation, 4) + 
        d(precipitation, 4 * 12) +  
        L(d(precipitation), c(0, 1, 4, 4 * 12)) +
        d(Iprec, 4) + 
        d(Iprec, 4 * 12) + 
        L(d(Iprec), c(0, 1, 4, 4 * 12)), data = weeklyData) %>% summary() # 0.5219

# Weekly - Values
dynlm( PM_US.Post~ 
         L(PM_US.Post, c(1, 4, 4 * 12)) +
         L(DEWP, c(0, 1, 4, 4 * 12)) +
         L(HUMI, c(0, 1, 4, 4 * 12))+
         L(PRES, c(0, 1, 4, 4 * 12))+
         L(TEMP, c(0, 1, 4, 4 * 12))+
         L(Iws, c(0, 1, 4, 4 * 12))+
         L(precipitation, c(0, 1, 4, 4 * 12))+
         L(Iprec, c(0, 1, 4, 4 * 12)), data = weeklyData) %>% summary() # 0.5228



#### Final models
## Hourly - Values (Optimized)
hourlyARDL <- dynlm( PM_US.Post ~ 
         L(PM_US.Post) +
         L(DEWP, c(0, 1, 24 * 365)) +
         L(HUMI, c(0, 1))+
         L(PRES, c(24 * 30))+
         L(TEMP, c(0, 1, 24 * 7))+
         Iws, data = data2) 
hourlyARDL %>% summary() # 0.89

## Daily - Values (Optimized)
dailyARDL <- dynlm( PM_US.Post ~ 
         L(PM_US.Post, c(1, 10)) +
         L(HUMI) +
         L(TEMP) +
         Iws+ 
         L(Iws)+
         precipitation, data = dailyData) 
dailyARDL %>% summary() # 0.44

## Weekly - Values (Optimized)
weeklyARDL <- dynlm( PM_US.Post~ 
         L(PM_US.Post) +
         DEWP +
         L(DEWP, 4) +
         PRES +
         L(TEMP)+
         Iws +
         precipitation + 
         L(precipitation, 4), data = weeklyData) 
weeklyARDL %>% summary() # 0.54



#### Plots of the final models
# Hourly
par(mfrow=c(2,2))
plot(hourlyARDL)
# Residuals vs Fitted - mostly linear, some outliers
# Normal Q-Q - follows strictly the normal dist, exept for the most extreme values that are way off
# Scale-Location - values lie mostly in the left-hand side corner
# Residuals vs Leverage - there is one observation outside the Cook's distance

# Daily
plot(dailyARDL)
# Residuals vs Fitted - no extra patterns, the residuals form a pack in the middle
# Normal Q-Q - some deviation from the dist in the extreme positive side
# Scale-Location - The points are not spread out evenly, but sid mostly in the middle
# Residuals vs Leverage - there is one observation outside the Cook's distance (different than the one in hourly)

# Weekly
plot(weeklyARDL)
# Residuals vs Fitted - mostly linear
# Normal Q-Q - very nice distribution, except for some outliers
# Scale-Location - no patterns observed
# Residuals vs Leverage - One significant outlier
par(mfrow=c(1,1)) # Don't forget to reset the plots grid 



#### Comparing predictions with real data - on graphs
## Hourly
# The values are ver simmilar, only small differences can be seen
ggplot() +
  geom_line(data = as.data.frame(hourlyARDL$fitted.values),
            aes(x = c(1:length(hourlyARDL$fitted.values)), 
                y = hourlyARDL$fitted.values),
            color = "red", 
            alpha=0.8) +
  geom_line(data = data2[1:length(hourlyARDL$fitted.values),],
            aes(x = c(1:length(hourlyARDL$fitted.values)),
                y = PM_US.Post), 
            color = "green", 
            alpha=0.5)

## Daily
# The model predicted the general trend well.
# It wasn't able to predict the significant outliers, and it generated it's own that were not present in the data (mostly lower that average)
ggplot() +
  geom_line(data = as.data.frame(dailyARDL$fitted.values),
            aes(x = c(1:length(dailyARDL$fitted.values)),
                y = dailyARDL$fitted.values), 
            color = "blue", 
            alpha=0.8) +
  geom_line(data = dailyData[1:length(dailyARDL$fitted.values),],
            aes(x = c(1:length(dailyARDL$fitted.values)),
                y = PM_US.Post), 
            color = "green", 
            alpha=0.5)

## Weekly
# The model predicted the trend very well
# Like previously, it wasn't able to predict the 3 outliers. It did not generate it's own
ggplot() +
  geom_line(data = as.data.frame(weeklyARDL$fitted.values), 
            aes(x = c(1:length(weeklyARDL$fitted.values)),
                y = weeklyARDL$fitted.values), 
            color = "red",
            alpha=0.8) +
  geom_line(data = weeklyData[1:length(weeklyARDL$fitted.values),], 
            aes(x = c(1:length(weeklyARDL$fitted.values)), 
                y = PM_US.Post), 
            color = "green", 
            alpha=0.5)


#### Let's plot the residuals
## Hourly
# Some serious outliers, otherwise the residuals seem to form a stationary process
hourlyARDL %>% 
  residuals() %>% 
  autoplot() 

## Daily
# Around 0, some sesonality can be observed
dailyARDL %>% 
  residuals() %>% 
  autoplot()

## Weekly
# Some mild outliers, otherwise stationary
weeklyARDL %>% 
  residuals() %>% 
  autoplot()


#### Breusch-Godfrey test for serial correlation of the residuals
## Hourly 
# Correlation
bgtest(residuals(hourlyARDL)~1, order = 1)
bgtest(residuals(hourlyARDL)~1, order = 2)
bgtest(residuals(hourlyARDL)~1, order = 3)
bgtest(residuals(hourlyARDL)~1, order = 4)
bgtest(residuals(hourlyARDL)~1, order = 5)

## Daily 
# No Correlation
bgtest(residuals(dailyARDL)~1, order = 1)
bgtest(residuals(dailyARDL)~1, order = 2)
bgtest(residuals(dailyARDL)~1, order = 3)
bgtest(residuals(dailyARDL)~1, order = 4)
bgtest(residuals(dailyARDL)~1, order = 5)

## Weekly
# No Correlation
bgtest(residuals(weeklyARDL)~1, order = 1)
bgtest(residuals(weeklyARDL)~1, order = 2)
bgtest(residuals(weeklyARDL)~1, order = 3)
bgtest(residuals(weeklyARDL)~1, order = 4)
bgtest(residuals(weeklyARDL)~1, order = 5)


#### Jarque - Bera Normality Test of the residuals
# Residuals are not normally distributed for all models
hourlyARDL %>% 
  residuals() %>% 
  as.matrix() %>% 
  jbTest()

dailyARDL %>% 
  residuals() %>% 
  as.matrix() %>%
  jbTest()

weeklyARDL %>% 
  residuals() %>%
  as.matrix() %>%
  jbTest()


#### Breusch-Pagan test for homoskedacity
hourlyARDL %>% bptest(data=data2) # heteroscedasticity
dailyARDL  %>% bptest(data=dailyData) # heteroscedasticity
weeklyARDL %>% bptest(data=weeklyData) # near the deciosion boundary, may be homoskedastic


#### Final 






