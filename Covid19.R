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

# New Deaths
newdeaths_pl <- ggplot(data_pre_can, aes(x=date, y=new_deaths)) +
  geom_line() + scale_x_date(date_breaks = "3 months", date_labels = "%b") +
  xlab("Months") + ylab("New Deaths") + ggtitle("New Covid Deaths in Canada")+geom_area(fill="lightblue", color="black") + 
  geom_smooth(method = lm, col = "red", se = FALSE)

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

### Preliminary Data Review

# None of the variables are stationary, will need to be transformed before analysis
# The trends of Covid cases, tests, and vaccinations are increasing, while Covid deaths are decreasing

##############################################################################################################################################################

#### Convert Data to Time Series Object
### NOTE: R cannot handle more than 350 periods in a timeseries
## Therefore, this analysis will be from the last 90 days (3 months)
### Method 1: base ts()
## Start takes c(year, period)
# Here, period is a day

# Convert Date to POSIXct class (storing time information)
data_pre$date <- as.POSIXct(data_pre$date)
data_ts_can <- ts(data_pre_can, 
                  start = c(2020, strftime(data_pre_can[1, 'date'], format = "%j")), 
                  end = c(2022, strftime(data_pre_can[length(data_pre_can$date), 'date'], format = "%j")), 
                  frequency = 365)

data_ts_can[is.na(data_ts_can)] <- 0

### Method 2: zoo()
newcases_analysis <- data_pre_can[ , c('date', 'new_cases')]
newcases_zoo <- zoo(newcases_analysis$new_cases, 
                seq(from = as.Date(data_pre_can[length(data_pre_can$date), 'date']) - 90, 
                to = as.Date(data_pre_can[length(data_pre_can$date), 'date']),
                by = 1))
plot(newcases_zoo)

newdeaths_analysis <- data_pre_can[ , c('date', 'new_deaths')]
newdeaths_analysis[is.na(newdeaths_analysis)] <- 0
newdeaths_zoo <- zoo(newdeaths_analysis$new_deaths, 
                    seq(from = as.Date(data_pre_can[length(data_pre_can$date), 'date']) - 90, 
                        to = as.Date(data_pre_can[length(data_pre_can$date), 'date']),
                        by = 1))

newvax_analysis <- data_pre_can[ , c('date', 'new_vaccinations')]
newvax_analysis[is.na(newvax_analysis)] <- 0
newvax_zoo <- zoo(newvax_analysis$new_vaccinations,
                  seq(from = as.Date(data_pre_can[length(data_pre_can$date), 'date']) - 90, 
                      to = as.Date(data_pre_can[length(data_pre_can$date), 'date']),
                      by = 1))

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

### Evaluate if data is stationary

## Visually
plot(newcases_stn)
plot(newdeaths_stn)
plot(newvax_stn)

## Statistically
# Augmented Dickey-Fuller Test.
# H_0 = The null hypothesis for this test is that there is some unit root (i.e. trend).
# H_A = The alternative hypothesis is that the time series is stationary (or trend-stationary).

adf.test(as.matrix(newcases_stn))
adf.test(as.matrix(newdeaths_stn))
adf.test(as.matrix(newvax_stn))

### Evaluate Autocorrelation within the data
## Given a trend, data will be correlated with itself
## By removing the implicit trends within data, remove the self-correlation
acf(newcases_stn)
pacf(newcases_stn)
acf(newdeaths_stn)

newcases_arima <- auto.arima(newcases_stn)
checkresiduals(newcases_arima)
autoplot(forecast(newcases_arima))
newcases_fcst <- forecast(newcases_arima, h = 7) # 7 week horizon
plot(newcases_fcst)

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

## use auto.arima to choose ARIMA terms
fit <- auto.arima(data_ts_can[ , 'new_cases'])

checkresiduals(fit)

### The model is a poor fit
### Use a Box-Cox transformation to correct for inconsistent variance
plot.ts(BoxCox(data_ts_can[ , 'new_cases'], 
       lambda = BoxCox.lambda(data_ts_can[ , 'new_cases'])))

## forecast for next 60 time points
fore <- forecast(fit, h = 365)
## plot it
plot(fore)

### Plot the Seasonal Data
ts.stl <- stl(TS,"periodic")  # decompose the TS
ts.sa <- seasadj(ts.stl)  # de-seasonalize
plot(AirPassengers, type="l")  # original series
plot(ts.sa, type="l")  # seasonal adjusted

### Create Regressors
regs <- cbind(NewDeaths = data_ts_can[ , 'new_deaths'],
             NewTests = data_ts_can[ , 'new_tests'],
             NewVax = data_ts_can[ , 'new_vaccinations'])

### Make the ARIMA Model
## Will allegorically determine the best model

arimafit <- auto.arima(data_ts_can[ , 'new_cases'],
                       xreg = regs)




arimaforecast <- forecast(arimafit, xreg = regs, h = 10)

plot(arimaforecast)

data_xts <- xts(data_pre_can, order.by = data_pre_can$date)

View(data_ts)
#### Explore data
covid_deaths_plt <- ggplot(data_ts, aes(date, total_deaths, color = location)) +
  geom_point() +
  facet_grid(facet = vars(continent), color ~ location) +
  theme(legend.position = "none")

covid_deaths_intplt <- ggplotly(covid_deaths_plt)
covid_deaths_intplt
covid_deaths_plt

### Declare Time Series 1
CovidDeaths_can <- ts(data2_canada$total_deaths, start = c(2020,2), frequency = 365)
CovidDeaths_can

#### Plot the Data: Total Deaths in Canada Over Time

autoplot(CovidDeaths_can) +
  ggtitle('COVID-19 Deaths in Canada (Feb 2020 - Feb 2022)') +
  ylab('Covid-19 Deaths')

#### NOTE:
#### As Total Deaths entails a positive trend, investigate transformations for stationary data

#### Transformations

#### Differences
#### Take the first difference of the data
#### Differencing a time series means, to subtract each data point in the series from its successor.
#### The first difference is the change in data over the unit of frequency
#### In this case, the increase in total deaths each day
CovidDeaths_can_diff <- diff(CovidDeaths_can)

#### Plot the Data: Daily Change in Covid-19 Deaths in Canada
autoplot(CovidDeaths_can_diff) +
  ggtitle('Daily Change COVID-19 Deaths in Canada (Feb 2020 - Feb 2022)') +
  ylab('Covid-19 Deaths')

#### The series appears trend stationary, and can be used to investigate seasonality
#### Use Augmented Dickey-Fuller Test (adf test). A p-Value of less than 0.05 in adf.test() indicates that it is stationary.
adf.test(na.omit(CovidDeaths_can_stn))
## ADF p-value = 0.6557
## The data is not stationary. Investigate further.

#### Make the data stationary; Another method:
ndiffs(CovidDeaths_can)  # number for seasonal differencing needed
# Number of differences = 2
CovidDeaths_can_stn <- diff(CovidDeaths_can, differences = 2)
plot(CovidDeaths_can_stn, type="l", main="Covid Deaths Stationary (d=2)")

#### Seasonality
ggseasonplot(diff_Y_can)
### There appear to be seasonal patterns to the data
### Such that there is a yearly rise in the daily change of total deaths in Nov-Dec 
### which briefly dips in the 2-3 quarter of the year

### Plot the daily change across years
ggsubseriesplot(diff_Y_can)

### Declare Time Series 2
covid_newdeaths_can <- xts(x = data2_canada$new_deaths, order.by = data2_canada$date, frequency = 365)

plot(covid_newdeaths_can)

