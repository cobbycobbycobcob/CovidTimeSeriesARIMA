###############################################################################################################################################################
#### COVID-19
# Time series analysis of COVID-19 infections, deaths, and vaccinations
###############################################################################################################################################################

#### Naming Convention
# pre = pre processing
# ts = time series

###############################################################################################################################################################

#### Time Series Analysis
### Models a dependent variable based on one explanatory variable: time
# Y = dependent
# X = time
### Components:
# Trend
# Seasonality
# Irregular fluctuations
# ARIMA 

###############################################################################################################################################################

### Import Libraries
library(AER)
library(tseries)
library(dynlm)
library(forecast)
library(readxl)
library(stargazer)
library(scales)
library(quantmod)
library(urca)
library(tidyverse) # Data management and manipulation
library(fpp2) # Forecasting
library(prophet) # Time series analyses
library(xts) # Extensible time series, i.e. time series matrices
library(ggplot2) # Visuals
library(ggfortify) # Interactive visualizations
library(plotly) # Interactive visualizations
library(leaflet) # Maps
library(lubridate)

### Import Data
data_pre = read.csv('https://covid.ourworldindata.org/data/owid-covid-data.csv')

### Select Columns
data_pre <- select(data_pre, c("date", "location", "continent", "population",
                         "new_cases", "new_deaths", "new_tests",
                         "new_vaccinations"))

#### Convert Character Columns to Factors
data_pre <- as.data.frame(unclass(data_pre),
                          stringsAsFactors = TRUE)

#### Remove NA and extraneous data rows
data_pre <- data_pre[!(is.na(data_pre$continent) | data_pre$continent == ""), ]

## Create day, month, year columns
# for time series and seasonality analyses
data_pre$day <- as.numeric(format(data_pre$date, '%d'))
data_pre$month <- as.numeric(format(data_pre$date, '%m'))
data_pre$year <- as.numeric(format(data_pre$date, '%Y'))

### Separate By Country
data_pre_can <- subset(data_pre, data_pre['location']  == 'Canada')
### Remove Country and Location columns
data_pre_can <- subset(data_pre_can, select = -c(location, continent))

##############################################################################################################################################################
### Explore the Data
### Plot each variable
data_pre_can$date <- as.Date(data_pre_can$date) ### GGPlot date_trans works with Date class objects

# New Cases
newcases_pl <- ggplot(data_pre_can, aes(x=date, y=new_cases)) +
  geom_line() + scale_x_date(date_breaks = "3 months", date_labels = "%b") +
  xlab("Months") + ylab("New Cases") + ggtitle("New Covid Cases in Canada")+geom_area(fill="lightblue", color="black") + 
  geom_smooth(method = lm, col = "red", se = FALSE)
newcases_pl

# New Deaths
newdeaths_pl <- ggplot(data_pre_can, aes(x=date, y=new_deaths)) +
  geom_line() + scale_x_date(date_breaks = "3 months", date_labels = "%b") +
  xlab("Months") + ylab("New Deaths") + ggtitle("New Covid Deaths in Canada")+geom_area(fill="lightblue", color="black") + 
  geom_smooth(method = lm, col = "red", se = FALSE)
newdeaths_pl

# New Vaccinations
newvax_pl <- ggplot(data_pre_can, aes(x=date, y=new_vaccinations)) +
  geom_line() + 
  scale_x_date(date_breaks = "3 months", date_labels = "%b") +
  xlab("Months") + ylab("New Vaccinations") + ggtitle("New Covid Vaccinations in Canada")+geom_area(fill="lightblue", color="black") + 
  geom_smooth(method = lm, col = "red", se = FALSE)
newvax_plot

# New Tests
newtests_pl <- ggplot(data_pre_can, aes(x=date, y=new_tests)) +
  geom_line() +
  scale_x_date(date_breaks = '3 months', date_labels = '%b') +
  xlab('Months') + ylab('New Tests') + ggtitle('New Covid Tests in Canada') + geom_area(fill='lightblue', color = 'black') +
  geom_smooth(method = lm, col = 'red', se = FALSE)
newtests_pl

# None of the variables are stationary, will need to be transformed before analysis
# The trends of Covid cases, tests, and vaccinations are increasing, while Covid deaths are decreasing

##############################################################################################################################################################

#### Convert Data to Time Series Object
### NOTE: R cannot handle more than 350 periods in a time series
## Therefore, this analysis will be from the last 8 months

## Set the time frame
data_timeframed <- data_pre_can %>%
  filter(data_pre_can$date > data_pre_can[length(data_pre_can$date), 'date'] %m-% months(8), 
         data_pre_can$date < data_pre_can[length(data_pre_can$date), 'date'])

## Create time serires as zoo objects
newcases_zoo <- zoo(data_timeframed[ , 'new_cases'], 
                order.by = data_timeframed$date)

newdeaths_zoo <- zoo(data_timeframed[ , 'new_deaths'], 
                    order.by = data_timeframed$date)

newtests_zoo <- zoo(data_timeframed[ , 'new_tests'], 
                    order.by = data_timeframed$date)

newvax_zoo <- zoo(data_timeframed[ , 'new_vaccinations'], 
                    order.by = data_timeframed$date)

###
## Decompose time series into to component waveforms
can_newcases_decomp <- decompose(newcases_zoo)

can_newdeaths_decomp <- decompose(data_ts_can[ , 'new_deaths'])

can_newvax_decomp <- decompose(data_ts_can[ , 'new_vaccinations'])

can_newtests_decomp <- decompose(data_ts_can[ , 'new_tests'])

### Earlier, trends were observed in the data
### Trends in the data violate the assumption for an ARIMA analysis
### Stationary data can be acheived by applying difference transformation to the data
### With trends differenced out, the data is ready for ARIMA analysis

newcases_stn <- diff(newcases_zoo)
newdeaths_stn <- diff(newdeaths_zoo)
newvax_stn <- diff(newvax_zoo)
newtests_stn <- diff(newtests_zoo)

### Evaluate if data is stationary

## Visually
# Blue is original
# Red is stable
plot.zoo(cbind(newcases_stn, newcases_zoo), 
         plot.type = "single", 
         col = c("red", "blue"))

plot.zoo(cbind(newdeaths_stn, newdeaths_zoo),
                             plot.type = 'single',
                             col = c('red', 'blue'))

plot.zoo(cbind(newtests_stn, newtests_zoo),
                             plot.type = 'single',
                             col = c('red', 'blue'))

plot.zoo(cbind(newvax_stn, newvax_zoo),
                             plot.type = 'single',
                             col = c('red', 'blue'))

## Statistically
# Augmented Dickey-Fuller Test.
# H_0 = The null hypothesis for this test is that there is some unit root (i.e. trend).
# H_A = The alternative hypothesis is that the time series is stationary (or trend-stationary).

adf.test(as.matrix(newcases_stn))
adf.test(as.matrix(newdeaths_stn))
adf.test(as.matrix(newvax_stn))
adf.test(as.matrix(newvax_stn))

### Evaluate Autocorrelation within the data
## Given a trend, data will be correlated with itself
## By removing the implicit trends within data, remove the self-correlation
acf(newcases_stn)
acf(newdeaths_stn)
acf(newtests_stn)
acf(newvax_stn)

### Construct ARIMA Models
## New Cases
newcases_arima <- auto.arima(newcases_stn)
## New Deaths
newdeaths_arima <- auto.arima(newdeaths_stn)
## New Tests
newtests_arima <- auto.arima(newtests_stn)
## New Vax
newvax_arima <- auto.arima(newvax_stn)

### Evaluate Models
checkresiduals(newcases_arima)
checkresiduals(newdeaths_arima)
checkresiduals(newtests_arima)
checkresiduals(newvax_arima)

### Generate & Plot Forecasts
newcases_fcst <- forecast(newcases_arima, h = 21)
newdeaths_fcst <- forecast(newdeaths_arima, h = 21)
newtests_fcst <- forecast(newtests_arima, h = 21)
newvax_fcst <- forecast(newvax_arima, h = 21)

plot(newcases_fcst)
plot(newdeaths_fcst)
plot(newvax_fcst)
plot(newtests_fcst)

### Point Estimates
# Cumulative amounts, i.e. sum()
newcases_3wk <- round(sum(newcases_fcst$upper[,2]),0)
newdeaths_3wk <- round(sum(newdeaths_fcst$upper[,2]),0)
newtests_3wk <- round(sum(newtests_fcst$upper[,2]),0)
newvax_3wk <- round(sum(newvax_fcst$upper[,2]),0)

### Collate Results
predictions_3weeks <- data.frame(
  NewCases = newcases_3wk,
  NewDeaths = newdeaths_3wk,
  NewVaccinations = newvax_3wk,
  NewTests = newtests_3wk
)
  
predictions_3weeks 

### Preliminary Evaluation
## The forecasts seem relatively flat
# This may be due to the long seasonality of the daily data, i.e. greater than 350 periods
# To correct for this, we can include a seasonal dummy variable extracted from the data via fourier transformations
# For the fourier transformation parameter k is 4, for 4 cycles, corresponding to the 4 yearly seasons in Canada
# This has been advised by Rob Hyndman

newcases_arimafourier <- auto.arima(newcases_zoo, seasonal=FALSE, xreg=fourier(newcases_zoo, 4))

newcases_7mo <- round(sum(forecast1$upper[,2]),0)

arima_newdeaths <- auto.arima(newdeaths_stn)
forecast2 <- forecast(arima_newdeaths, h = 49)
plot(forecast2)

