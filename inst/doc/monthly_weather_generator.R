## ----libraries, warning=FALSE, message=FALSE-----------------------------
library(forecast)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(ggplot2)
theme_set(theme_bw())
library(weatherGen)
data(climate)

# set RNG seed for reproducibility
set.seed(2345)

## ----db grid-------------------------------------------------------------
DB_PATH <- '/Users/jeff/Projects/UMass/climate-data/maurer/db/maurer_mon.db'
db <- src_sqlite(DB_PATH, create = FALSE)

df.grid <- tbl(db, "grid") %>%
  collect

## ----map east, warning=FALSE---------------------------------------------
library(ggmap)
map <- get_map(location=c(lon=mean(range(df.grid$LON)), lat=mean(range(df.grid$LAT))),
               zoom=4, maptype="satellite", color='bw')
ggmap(map, darken=c(0.25, "white"), extent="device") +
  geom_point(aes(x=LON, y=LAT), data=df.grid, color='red', size=1)

## ----map new england, warning=FALSE--------------------------------------
map <- get_map(location=c(lon=-72, lat=42),
               zoom=7, maptype="satellite", color='bw')
ggmap(map, darken=c(0.25, "white"), extent="device") +
  geom_point(aes(x=LON, y=LAT, color=SELECT), 
             data=mutate(df.grid, SELECT=(LAT>=41 & LAT<=43 & LON>=-74 & LON<=-70)), size=2) +
  geom_vline(aes(xintercept=LON), color='green', linetype=2, 
             data=data.frame(LON=seq(-74, -70))) +
  geom_hline(aes(yintercept=LAT), color='green', linetype=2, 
             data=data.frame(LAT=seq(41, 43))) +
  scale_color_manual(values=c("TRUE"="green", "FALSE"="red"), guide=FALSE)

## ----climate load--------------------------------------------------------
clim.mon <- tbl(db, "data") %>%
  filter(LAT>=41, LAT<=43, LON>=-74, LON<=-70) %>%
  collect %>%
  select(-REGION, -WIND) %>%
  mutate(TAVG=(TMIN+TMAX)/2) %>%
  gather(VAR, VALUE, PRCP, TMIN, TMAX, TAVG)
summary(clim.mon)

## ----climate yr----------------------------------------------------------
clim.yr <- spread(clim.mon, VAR, VALUE) %>%
  mutate(MONTHYEAR=ymd(paste(YEAR,MONTH,1,sep='-')),
         N_DAY=days_in_month(MONTHYEAR)) %>%
  group_by(LAT, LON, YEAR) %>%
  summarise(PRCP=sum(PRCP),
            TMAX=sum(TMAX*N_DAY)/sum(N_DAY),
            TMIN=sum(TMIN*N_DAY)/sum(N_DAY),
            TAVG=sum(TAVG*N_DAY)/sum(N_DAY)) %>%
  gather(VAR, VALUE, PRCP:TAVG)
summary(clim.yr)

## ----climate reg mon-----------------------------------------------------
clim.reg.mon <- clim.mon %>%
  group_by(YEAR, MONTH, VAR) %>%
  summarise(MEAN=mean(VALUE),
            SD=sd(VALUE))
summary(clim.reg.mon)

## ----plot clim reg mon---------------------------------------------------
mutate(clim.reg.mon, DATE=ymd(paste(YEAR,MONTH,1,sep='-'))) %>%
  ggplot(aes(DATE, MEAN)) +
  geom_ribbon(aes(ymin=MEAN-SD, ymax=MEAN+SD), fill='grey50') +
  geom_line(color='blue') +
  facet_wrap(~VAR, ncol=1, scales='free_y') +
  labs(x="Month/Year", y="Mean +/- StDev")

## ----climate reg yr------------------------------------------------------
clim.reg.yr <-  clim.yr %>%
  group_by(VAR, YEAR) %>%
  summarise(MEAN=mean(VALUE),
            SD=sd(VALUE))
summary(clim.reg.yr)

## ----plot clim reg yr----------------------------------------------------
ggplot(clim.reg.yr, aes(YEAR, MEAN)) +
  geom_ribbon(aes(ymin=MEAN-SD, ymax=MEAN+SD), fill='grey80') +
  geom_line(color='blue') +
  facet_wrap(~VAR, ncol=1, scales='free_y') +
  labs(x="Year", y="Mean +/- StDev")

## ----clim reg wide-------------------------------------------------------
clim.reg.mon <- select(clim.reg.mon, -SD) %>%
  spread(VAR, MEAN)
head(clim.reg.mon)

clim.reg.yr <- select(clim.reg.yr, -SD) %>%
  spread(VAR, MEAN)
head(clim.reg.yr)

## ----plot wave, message=FALSE, warning=FALSE-----------------------------
wt.reg <- wavelet_analysis(clim.reg.yr$PRCP, sig.level=0.90, noise.type='white')

par(mfrow=c(1,2))
plot(wt.reg, plot.cb=TRUE, plot.phase=FALSE)
plot(wt.reg$gws, wt.reg$period, type="b",
     xlab="Global Wavelet Spectrum", ylab="",
     log="y", ylim=rev(range(wt.reg$period)), 
     xlim=range(c(wt.reg$gws, wt.reg$gws.sig$signif)))
lines(wt.reg$gws.sig$signif, wt.reg$period, lty=2, col='red')

## ----arima model---------------------------------------------------------
models <- list(PRCP=auto.arima(clim.reg.yr$PRCP,
                               max.p=2,max.q=2,max.P=0,max.Q=0,
                               stationary=TRUE))
models[["PRCP"]]

## ----plot arima forecast-------------------------------------------------
par(mfrow=c(1,1))
plot(forecast(models[['PRCP']], h=20), main="ARIMA Forecast of Regional Annual Precip",
     xlab="Timestep (yr)", ylab='Annual Precip (mm/yr)')

## ----mc ar, warning=FALSE------------------------------------------------
sim.prcp <- weatherGen::mc_arimas(models=models, n=nrow(clim.reg.yr), n.iter=50)
str(sim.prcp)

## ----plot gws stats, fig.height=8----------------------------------------
par(mfrow=c(2,2))
x1 <- 1
x2 <- min(length(wt.reg$period), nrow(sim.prcp$gws))
ymin <- min(wt.reg$period[x1:x2], wt.reg$gws[x1:x2],
            sim.prcp$gws.stat[x1:x2, 'MEAN'],
            sim.prcp$gws.stat[x1:x2, 'Q025'])
ymax <- max(wt.reg$period[x1:x2], wt.reg$gws[x1:x2],
            sim.prcp$gws.stat[x1:x2, 'MEAN'],
            sim.prcp$gws.stat[x1:x2, 'Q975'])
plot(wt.reg$period[x1:x2], wt.reg$gws[x1:x2],
     type="n", ylim=c(ymin,ymax), xlab="Period (years)", ylab="",
     main="", log="y")
mtext(side=2, expression(paste("Power (",mm^2,")")), line=2.5)
polygon(c(wt.reg$period[x1:x2],
          rev(wt.reg$period[x1:x2])),
        c(sim.prcp$gws.stat[x1:x2, 'Q025'],
          rev(sim.prcp$gws.stat[x1:x2, 'Q975'])),
        col="grey")
lines(wt.reg$period[x1:x2], wt.reg$gws[x1:x2])
lines(wt.reg$period[x1:x2], sim.prcp$gws.stat[x1:x2, 'MEAN'], lty=2, col="blue")
lines(wt.reg$period[x1:x2], wt.reg$gws.sig$signif[x1:x2], col="red", lty=3)

boxplot(colMeans(sim.prcp$x), main="Mean")
points(mean(clim.reg.yr$PRCP),col="red",pch=19)
points(mean(colMeans(sim.prcp$x)), col="blue", pch=19)

boxplot(apply(sim.prcp$x, 2, sd), main="Standard Deviation")
points(sd(clim.reg.yr$PRCP), col="red", pch=19)
points(mean(apply(sim.prcp$x, 2, sd)), col="blue", pch=19)

boxplot(apply(sim.prcp$x, 2, skewness), main="Skew")
points(skewness(clim.reg.yr$PRCP), col="red", pch=19)
points(mean(apply(sim.prcp$x, 2, skewness)),col="blue",pch=19)

## ----plot ar sims--------------------------------------------------------
as.data.frame(sim.prcp$x) %>%
  mutate(TIMESTEP=row_number()) %>% 
  gather(TRIAL, VALUE, -TIMESTEP) %>% 
  ggplot(aes(TIMESTEP, VALUE)) + 
  geom_line(aes(group=TRIAL), alpha=0.3) + 
  stat_summary(fun.y=mean, geom='line', color='red', size=1) +
  geom_line(aes(x=TIMESTEP, y=PRCP), data=mutate(clim.reg.yr, TIMESTEP=row_number()), color='blue', size=1) +
  labs(x="Year", y="Annual Precipitation (mm/yr)", 
       title="WARM Annual Precipitation Simulations\nRed Line=Mean of Simulations, Blue Line=Observd Regional Mean")

## ----select loc----------------------------------------------------------
clim.locs <- select(clim.yr, LAT, LON) %>%
  unique()
loc <- clim.locs[sample(nrow(clim.locs), size=1),] %>% as.numeric
names(loc) <- c('LAT', 'LON')
loc

## ----loc datasets--------------------------------------------------------
clim.loc.yr <- filter(clim.yr, LAT==loc[['LAT']], LON==loc[['LON']]) %>%
  spread(VAR, VALUE) %>%
  arrange(YEAR)
clim.loc.mon <- filter(clim.mon, LAT==loc[['LAT']], LON==loc[['LON']]) %>%
  spread(VAR, VALUE) %>%
  arrange(YEAR, MONTH)

## ----plot clim loc reg, fig.height=8-------------------------------------
clim.loc.mon %>%
  select(-LAT, -LON) %>%
  mutate(SOURCE="Location") %>%
  rbind(clim.reg.mon %>% mutate(SOURCE="Region")) %>%
  mutate(SOURCE=factor(SOURCE)) %>%
  gather(VAR, VALUE, PRCP:TAVG) %>%
  mutate(DATE=ymd(paste(YEAR, MONTH, 1, sep='-'))) %>%
  ggplot(aes(DATE, VALUE, color=SOURCE)) +
  geom_line() +
  labs(x="Month/Year", y="") +
  scale_color_manual('', values=c('steelblue', 'orangered')) +
  facet_wrap(~VAR, ncol=1, scales='free_y')

## ----distances-----------------------------------------------------------
j <- 1 # year number
samp <- 1 # simulation trial number

# compute distances between simulated WARM timeseries and observed regional precip
distances <- cbind(clim.reg.yr, DISTANCE=sqrt((sim.prcp$x[j,samp] - clim.reg.yr$PRCP)^2))

ggplot(distances, aes(YEAR, PRCP)) + 
  geom_line() +
  geom_point(aes(color=DISTANCE), size=3) +
  geom_hline(yint=sim.prcp$x[j,samp], linetype=2, color='red') +
  geom_text(x=max(distances$YEAR), y=sim.prcp$x[j,samp], label="Sim", hjust=1, vjust=-0.5) +
  scale_color_gradient(low='green', high='red') +
  labs(x="Year", y="Annual Precipitation (mm/yr)", 
       title="Regional and Current Simulation Annual Precipitation")

## ----plot top distances--------------------------------------------------
distances <- mutate(distances,
                    IN_TOP_8=rank(DISTANCE) %in% seq(1,8))
distances %>%
  ggplot(aes(YEAR, PRCP)) + 
  geom_line() +
  geom_point(aes(color=DISTANCE), size=3) +
  geom_hline(yint=sim.prcp$x[j,samp], linetype=2, color='red') +
  geom_point(mapping=aes(alpha=IN_TOP_8), color='black', shape=1, size=4) +
  geom_text(aes(x=YEAR, y=PRCP, label=YEAR), data=filter(distances, IN_TOP_8), vjust=-1) +
  geom_text(x=max(distances$YEAR), y=sim.prcp$x[j,samp], label="Sim", hjust=1, vjust=-0.5) +
  scale_alpha_manual(values=c("TRUE"=1, "FALSE"=0), guide=FALSE) +
  scale_color_gradient(low='green', high='red') +
  labs(x="Year", y="Annual Precipitation (mm)", title="Regional Annual Precipitation with 8 Nearest Neighbors Selected")

## ----top distances-------------------------------------------------------
distances.top <- filter(distances, IN_TOP_8) %>%
  mutate(WEIGHT=(1/DISTANCE)^2,
         WEIGHT=WEIGHT/sum(WEIGHT))
distances.top %>%
  ggplot(aes(PRCP, WEIGHT)) +
  geom_point() +
  geom_vline(xint=sim.prcp$x[j,samp], linetype=2, color='red') +
  geom_text(aes(label=YEAR), hjust=-0.1) +
  geom_text(x=sim.prcp$x[j,samp], y=0.4, label="Sim Precip", hjust=-0.1) +
  labs(x="Annual Precipitation (mm)", y="Sample Weight")

## ----sample years--------------------------------------------------------
sampled.years <- sample(distances.top$YEAR,
                        size=1000,
                        prob=distances.top$WEIGHT,
                        replace=TRUE)
sampled.years.tbl <- prop.table(table(sampled.years))
merge(distances.top, 
      data.frame(YEAR=names(sampled.years.tbl),
                 FREQ=as.vector(sampled.years.tbl)),
      all.x=TRUE) %>%
  select(YEAR, WEIGHT, FREQ) %>%
  ggplot(aes(WEIGHT, FREQ)) +
  geom_point(size=3) +
  geom_abline(linetype=2) +
  labs(x="Sample Weight", y="Sampled Frequency")

## ----sample year---------------------------------------------------------
sampled.year <- sample(distances.top$YEAR, size=1, prob=distances.top$WEIGHT, replace=TRUE)
sampled.year

## ----sampled monthly climate---------------------------------------------
filter(clim.loc.mon, YEAR==sampled.year)

## ----gen set up----------------------------------------------------------
num_sim <- ncol(sim.prcp$x)
num_year <- nrow(sim.prcp$x) # number of composited years
num_months <- num_year*12

climate_variables <- select(clim.reg.yr, -YEAR) %>%
  names()

sampled.year <- array(NA, c(num_year, num_sim))
gen.loc.mon.arr <- array(NA, c(num_months, length(climate_variables), num_sim))

## ----gen loop------------------------------------------------------------
# loop through simulations (trials)
for (samp in 1:num_sim) {
    # loop through years
    for (j in 1:num_year) {
        # compute distances between simulated WARM timeseries and observed regional precip
      distances <- mutate(clim.reg.yr, 
                          DISTANCE=sqrt((sim.prcp$x[j, samp] - clim.reg.yr$PRCP)^2),
                          IN_TOP_8=rank(DISTANCE) %in% seq(1,8))
      distances.top <- filter(distances, IN_TOP_8) %>%
                       mutate(WEIGHT=(1/DISTANCE)^2,
                              WEIGHT=WEIGHT/sum(WEIGHT))
      sampled.year[j,samp] <- sample(distances.top$YEAR,
                                     size=1,
                                     prob=distances.top$WEIGHT,
                                     replace=FALSE)
      row.idx <- ((1+12*(j-1)):(12+12*(j-1)))
      gen.loc.mon.arr[row.idx, , samp] <- data.matrix(filter(clim.loc.mon, 
                                                             YEAR==sampled.year[j, samp]) %>%
                                                      select(PRCP, TMAX, TMIN, TAVG))
    }
}

## ----gen loc mon---------------------------------------------------------
gen.loc.mon <- sapply(climate_variables, function(v) {
  i.var <- which(climate_variables==v)
  df <- as.data.frame(gen.loc.mon.arr[,i.var,]) %>%
    mutate(TIMESTEP=row_number(),
           SIM_YEAR=rep(seq(1, num_year), each=12),
           SIM_MONTH=rep(seq(1, 12), times=num_year),
           VAR=v) %>%
    gather(TRIAL, VALUE, -TIMESTEP, -SIM_YEAR, -SIM_MONTH, -VAR) %>%
    mutate(TRIAL=as.numeric(TRIAL))
  return(list(df))
})
str(gen.loc.mon)

## ----gen loc yr----------------------------------------------------------
gen.loc.yr <- lapply(gen.loc.mon, function(df) {
  if ('PRCP' %in% unique(df$VAR)) {
    df.yr <- group_by(df, VAR, TRIAL, SIM_YEAR) %>%
      summarise(VALUE=sum(VALUE))
  } else {
    df.yr <- group_by(df, VAR, TRIAL, SIM_YEAR) %>%
      summarise(VALUE=mean(VALUE))
  }
  df.yr <- df.yr %>%
    ungroup() %>%
    as.data.frame()
  return(df.yr)
})
str(gen.loc.yr)

## ----plot gen loc yr-----------------------------------------------------
do.call(rbind, gen.loc.yr) %>%
  mutate(VAR=factor(VAR)) %>%
  ggplot() +
    geom_line(aes(SIM_YEAR, VALUE, group=TRIAL), alpha=0.5) +
    geom_line(aes(SIM_YEAR, VALUE), 
              data=mutate(clim.reg.yr, SIM_YEAR=row_number()) %>%
                gather(VAR, VALUE, PRCP:TAVG), 
              color='red') +
    facet_wrap(~VAR, scales='free_y') +
    labs(x="Simulation Year", y="Value", title="All Generated Datasets with Observed Location Dataset ")

## ----plot gen loc yr one-------------------------------------------------
do.call(rbind, gen.loc.yr) %>%
  mutate(VAR=factor(VAR)) %>%
  filter(TRIAL==1) %>%
  ggplot() +
    geom_line(aes(SIM_YEAR, VALUE, group=TRIAL), alpha=0.5) +
    geom_line(aes(SIM_YEAR, VALUE), 
              data=mutate(clim.reg.yr, SIM_YEAR=row_number()) %>%
                gather(VAR, VALUE, PRCP:TAVG), 
              color='red') +
    facet_wrap(~VAR, scales='free_y') +
    labs(x="Simulation Year", y="Value", title="Single Generated Dataset with Observed Location Dataset ")

## ----plot gen loc yr box mean, fig.width=4-------------------------------
do.call(rbind, gen.loc.yr) %>%
  mutate(VAR=factor(VAR)) %>%
  group_by(VAR, TRIAL) %>%
  summarise(MEAN=mean(VALUE)) %>%
  ggplot(aes(x=1, y=MEAN)) +
    geom_boxplot() +
    stat_summary(fun.y=mean, geom="point", color="red", size=3) +
    geom_point(aes(x=1, y=MEAN), 
               data=gather(clim.loc.yr, VAR, VALUE, PRCP:TAVG) %>% 
                      group_by(VAR) %>% 
                      summarise(MEAN=mean(VALUE)),
               color='blue', size=3) +
    facet_wrap(~VAR, scales='free_y') +
    labs(x="", y="Mean", title="Mean of Annual Precipitation for Generated (Red) and Observed (Blue)") +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

## ----plot gen loc yr box sd, fig.width=4---------------------------------
do.call(rbind, gen.loc.yr) %>%
  mutate(VAR=factor(VAR)) %>%
  group_by(VAR, TRIAL) %>%
  summarise(SD=sd(VALUE)) %>%
  ggplot(aes(x=1, y=SD)) +
    geom_boxplot() +
    stat_summary(fun.y=mean, geom="point", color="red", size=3) +
    geom_point(aes(x=1, y=SD), 
               data=gather(clim.loc.yr, VAR, VALUE, PRCP:TAVG) %>% 
                      group_by(VAR) %>% 
                      summarise(SD=sd(VALUE)),
               color='blue', size=3) +
    facet_wrap(~VAR, scales='free_y') +
    labs(x="", y="Standard Deviation", title="StDev of Annual Precipitation for Generated (Red) and Observed (Blue)") +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

## ----plot gen loc yr box skew, fig.width=4-------------------------------
do.call(rbind, gen.loc.yr) %>%
  mutate(VAR=factor(VAR)) %>%
  group_by(VAR, TRIAL) %>%
  summarise(SKEW=skewness(VALUE)) %>%
  ggplot(aes(x=1, y=SKEW)) +
    geom_boxplot() +
    stat_summary(fun.y=mean, geom="point", color="red", size=3) +
    geom_point(aes(x=1, y=SKEW), 
               data=gather(clim.loc.yr, VAR, VALUE, PRCP:TAVG) %>% 
                      group_by(VAR) %>% 
                      summarise(SKEW=skewness(VALUE)),
               color='blue', size=3) +
    facet_wrap(~VAR, scales='free_y') +
    labs(x="", y="Skewness", title="Skewness of Annual Precipitation for Generated (Red) and Observed (Blue)") +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

## ----wgen mon arma, warning=FALSE----------------------------------------
wgen_mon <- weatherGen::gen_month_arma(x.reg.yr=clim.reg.yr, x.loc.mon=clim.loc.mon,
                                       n.iter=50, n.year=nrow(clim.reg.yr))
str(wgen_mon, max.level=1)

## ----wgen mon example----------------------------------------------------
head(wgen_mon$sim.mon[[1]])

## ----plot wgen year compare----------------------------------------------
wgen_mon$sim.mon[[1]] %>%
  group_by(YEAR_SIM) %>%
  summarise(PRCP=sum(PRCP),
            YEAR_OBS=mean(YEAR_OBS),
            PRCP_OBS=mean(PRCP_OBS_YEAR)) %>%
  gather(VAR, VALUE, PRCP, PRCP_OBS) %>%
  ggplot(aes(YEAR_SIM, VALUE, color=VAR)) +
  geom_line() +
  labs(x="Year", y="Annual Precipitation (mm)")

## ----plot wgen mon-------------------------------------------------------
wgen_mon$sim.mon[[1]] %>%
  select(YEAR_SIM, MONTH, PRCP:TAVG) %>%
  gather(VAR, VALUE, PRCP:TAVG) %>%
  mutate(DATE=ymd(paste(YEAR_SIM, MONTH, 1, sep='-'))) %>%
  ggplot(aes(DATE, VALUE)) +
  geom_line() +
  facet_wrap(~VAR, scales='free_y', ncol=1) +
  labs(x="Month/Year", y="")

## ----compute trends------------------------------------------------------
trend.factors <- list(PRCP=1.5, # multiplicative trend factor
                      TEMP=2)   # additive trend factor

# select one trial at random
sim.gen <- wgen_mon$sim.mon[[sample(length(wgen_mon$sim.mon), 1)]]

# create linear trend lines
trend.lines <- list(PRCP=seq(1, trend.factors[['PRCP']], length.out=nrow(sim.gen)),
                    TEMP=seq(0, trend.factors[['TEMP']], length.out=nrow(sim.gen)))

sim.trend <- mutate(sim.gen,
                    PRCP.TREND=PRCP * trend.lines[['PRCP']],
                    TMAX.TREND=TMAX + trend.lines[['TEMP']],
                    TMIN.TREND=TMIN + trend.lines[['TEMP']],
                    TAVG.TREND=TAVG + trend.lines[['TEMP']])

## ----clim gen trend, fig.height=8----------------------------------------
sim.trend %>%
  mutate(PRCP.GEN=PRCP, TMAX.GEN=TMAX, TMIN.GEN=TMIN, TAVG.GEN=TAVG) %>%
  gather(VAR.GROUP, VALUE, PRCP.GEN:TAVG.GEN, PRCP.TREND:TAVG.TREND) %>%
  separate(VAR.GROUP, c("VAR", "GROUP")) %>%
  mutate(DATE=ymd(paste(YEAR_SIM, MONTH, 1, sep='-'))) %>%
  ggplot(aes(DATE, VALUE, color=GROUP)) +
  geom_line() +
  facet_wrap(~VAR, ncol=1, scales='free_y') +
  labs(x="Month/Year", y="") +
  scale_color_manual('', values=c(GEN='grey20', TREND='chartreuse3'),
                     labels=c(GEN='w/o trend', TREND='w/ trend'))

## ------------------------------------------------------------------------
trend.factors <- list(PRCP=seq(0.75, 1.25, by=0.05),
                      TEMP=seq(0, 5, by=0.5))

sim.trends <- trends_mon(sim.mon.list=wgen_mon$sim.mon, 
                         temp.factors=trend.factors[['TEMP']],
                         prcp.factors=trend.factors[['PRCP']])
str(sim.trends[[1]])

## ----trend df------------------------------------------------------------
df.trends <- lapply(sim.trends, function(x) {
  df <- x$df
  df$PRCP_FACTOR <- x$trend.factors[['PRCP']]
  df$TEMP_FACTOR <- x$trend.factors[['TEMP']]
  return(df)
}) %>%
  do.call(rbind, .)
summary(df.trends)

## ----plot prcp trends----------------------------------------------------
df.trends %>%
  filter(TEMP_FACTOR==0) %>%
  group_by(PRCP_FACTOR) %>%
  mutate(TIMESTEP=row_number()) %>%
  ggplot(aes(TIMESTEP, PRCP.TREND, color=PRCP_FACTOR)) +
  geom_line(alpha=0.2) +
  scale_color_gradient('Precip Factor', low='green', high='red') +
  labs(x="Simulation Month", y="Monthly Precipitation (mm)")

## ----plot tavg trends----------------------------------------------------
df.trends %>%
  filter(TEMP_FACTOR==0) %>%
  group_by(TEMP_FACTOR, PRCP_FACTOR, YEAR_SIM) %>%
  summarise(PRCP=sum(PRCP.TREND)) %>%
  ggplot(aes(YEAR_SIM, PRCP, color=PRCP_FACTOR, group=PRCP_FACTOR)) +
  geom_line(alpha=0.7) +
  scale_color_gradient('Precip Factor', low='green', high='red') +
  labs(x="Simulation Year", y="Annual Precipitation (mm)")

## ----plot temp trends----------------------------------------------------
df.trends %>% 
  filter(PRCP_FACTOR==1) %>%
  mutate(DATE=ymd(paste(2001, MONTH, 1, sep='-'))) %>%
  group_by(TEMP_FACTOR, PRCP_FACTOR, YEAR_SIM) %>%
  summarise(TAVG=sum(TAVG.TREND*days_in_month(DATE))/sum(days_in_month(DATE))) %>%
  ggplot(aes(YEAR_SIM, TAVG, color=TEMP_FACTOR, group=TEMP_FACTOR)) +
  geom_line(alpha=0.7) +
  scale_color_gradient('Precip Factor', low='green', high='red') +
  labs(x="Simulation Year", y="Annual Mean Temperature (degC)")

## ----session-------------------------------------------------------------
sessionInfo()

