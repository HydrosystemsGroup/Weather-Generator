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

## ----get file list-------------------------------------------------------
DATA_DIR <- '/Users/jeff/Projects/UMass/virtue/data/mon/east'
files <- list.files(DATA_DIR) %>%
  str_split('_') %>% 
  do.call(rbind, .) %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  mutate(FILE=paste(V1, V2, V3, sep="_")) %>%
  select(FILE, LAT=V2, LON=V3) %>%
  mutate(LAT=as.numeric(LAT), LON=as.numeric(LON))

## ----map east, warning=FALSE---------------------------------------------
library(ggmap)
map <- get_map(location=c(lon=mean(range(files$LON)), lat=mean(range(files$LAT))),
               zoom=4, maptype="satellite", color='bw')
ggmap(map, darken=c(0.25, "white"), extent="device") +
  geom_point(aes(x=LON, y=LAT), data=files, color='red', size=1)

## ----map new england, warning=FALSE--------------------------------------
map <- get_map(location=c(lon=-72, lat=42),
               zoom=7, maptype="satellite", color='bw')
ggmap(map, darken=c(0.25, "white"), extent="device") +
  geom_point(aes(x=LON, y=LAT, color=SELECT), 
             data=mutate(files, SELECT=(LAT>=41 & LAT<=43 & LON>=-74 & LON<=-70)), size=2) +
  geom_vline(aes(xintercept=LON), color='green', linetype=2, 
             data=data.frame(LON=seq(-74, -70))) +
  geom_hline(aes(yintercept=LAT), color='green', linetype=2, 
             data=data.frame(LAT=seq(41, 43))) +
  scale_color_manual(values=c("TRUE"="green", "FALSE"="red"), guide=FALSE)

## ----climate load--------------------------------------------------------
clim.raw <- filter(files, LAT>=41, LAT<=43, LON>=-74, LON<=-70) %>%
  apply(1, function(f) {
    x <- read.table(file.path(DATA_DIR, f[['FILE']]))
    names(x) <- c("YEAR", "MONTH", "PRCP","TMAX","TMIN","TAVG")
    x$LAT <- f['LAT']
    x$LON <- f['LON']
    return(x)
  }) %>%
  do.call(rbind, .)

clim.mon <- gather(clim.raw, VAR, VALUE, PRCP:TAVG) %>%
  mutate(VALUE=as.numeric(VALUE))

clim.yr <- group_by(clim.raw, YEAR, LAT, LON) %>%
  summarise(PRCP=sum(PRCP),
            TMAX=mean(TMAX),
            TMIN=mean(TMIN),
            TAVG=mean(TAVG)) %>%
  gather(VAR, VALUE, PRCP:TAVG)

## ----climate reg mon-----------------------------------------------------
clim.reg.mon <- clim.mon %>%
  group_by(YEAR, MONTH, VAR) %>%
  summarise(MEAN=mean(VALUE),
            SD=sd(VALUE))

mutate(clim.reg.mon, DATE=ymd(paste(YEAR,MONTH,1,sep='-'))) %>%
  ggplot(aes(DATE, MEAN)) +
  geom_ribbon(aes(ymin=MEAN-SD, ymax=MEAN+SD), fill='grey80') +
  geom_line(color='red') +
  facet_wrap(~VAR, ncol=1, scales='free_y') +
  labs(x="Month/Year", y="Mean +/- StDev")

## ----climate reg yr------------------------------------------------------
clim.reg.yr <-  clim.yr %>%
  group_by(VAR, YEAR) %>%
  summarise(MEAN=mean(VALUE),
            SD=sd(VALUE))

ggplot(clim.reg.yr, aes(YEAR, MEAN)) +
  geom_ribbon(aes(ymin=MEAN-SD, ymax=MEAN+SD), fill='grey80') +
  geom_line(color='red') +
  facet_wrap(~VAR, ncol=1, scales='free_y') +
  labs(x="Year", y="Mean +/- StDev")

## ----clim wide-----------------------------------------------------------
clim.reg.mon <- select(clim.reg.mon, -SD) %>%
  spread(VAR, MEAN)
head(clim.reg.mon)

clim.reg.yr <- select(clim.reg.yr, -SD) %>%
  spread(VAR, MEAN)
head(clim.reg.yr)

## ----plot wave, message=FALSE, warning=FALSE-----------------------------
wave <- wavelet_analysis(clim.reg.yr$PRCP, sig.level=0.90, noise.type='white')

par(mfrow=c(1,2))
plot(wave, plot.cb=TRUE, plot.phase=FALSE)
plot(wave$gws, wave$period, type="b",
     xlab="Global Wavelet Spectrum", ylab="Fourier Period (Years)",
     log="y", ylim=rev(range(wave$period)), 
     xlim=range(c(wave$gws, wave$gws.sig$signif)))
lines(wave$gws.sig$signif, wave$period, lty=2, col='red')

## ----arima model---------------------------------------------------------
models <- list(PRCP=auto.arima(clim.reg.yr$PRCP,
                               max.p=2,max.q=2,max.P=0,max.Q=0,
                               stationary=TRUE))
models[["PRCP"]]

## ----plot arima forecast-------------------------------------------------
par(mfrow=c(1,1))
plot(forecast(models[['PRCP']], h=20), main="ARIME Forecast of Regional Annual Precip",
     xlab="Timestep (yr)", ylab='Annual Precip (mm/yr)')

## ----mc ar, warning=FALSE------------------------------------------------
sim.prcp <- weatherGen::mc_arimas(models=models, n=nrow(clim.reg.yr), n.iter=50)
str(sim.prcp)

## ----plot gws stats, fig.height=8----------------------------------------
par(mfrow=c(2,2))
x1 <- 1
x2 <- min(length(wave$period), nrow(sim.prcp$gws))
ymin <- min(wave$period[x1:x2], wave$gws[x1:x2],
            sim.prcp$gws.stat[x1:x2, 'MEAN'],
            sim.prcp$gws.stat[x1:x2, 'Q025'])
ymax <- max(wave$period[x1:x2], wave$gws[x1:x2],
            sim.prcp$gws.stat[x1:x2, 'MEAN'],
            sim.prcp$gws.stat[x1:x2, 'Q975'])
plot(wave$period[x1:x2], wave$gws[x1:x2],
     type="n", ylim=c(ymin,ymax), xlab="Period (years)", ylab="",
     main="", log="y")
mtext(side=2, expression(paste("Power (",mm^2,")")), line=2.5)
polygon(c(wave$period[x1:x2],
          rev(wave$period[x1:x2])),
        c(sim.prcp$gws.stat[x1:x2, 'Q025'],
          rev(sim.prcp$gws.stat[x1:x2, 'Q975'])),
        col="grey")
lines(wave$period[x1:x2], wave$gws[x1:x2])
lines(wave$period[x1:x2], sim.prcp$gws.stat[x1:x2, 'MEAN'], lty=2, col="blue")
lines(wave$period[x1:x2], wave$gws.sig$signif[x1:x2], col="red", lty=3)

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

## ------------------------------------------------------------------------
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

## ----session-------------------------------------------------------------
sessionInfo()

