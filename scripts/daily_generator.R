library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)
theme_set(theme_bw())
library(gridExtra)
library(weathergen)
library(forecast)

# configuration ----

num_year_sim <- 40
water_yr_change <- TRUE
month_list <- c(10:12, 1:9)
START_WYEAR_SIM <- 2001
k1 <- 1
wet_spell_change <- rep(1, 12)
dry_spell_change <- rep(1, 12)

# functions ----
KNN_ANNUAL <- function(sim_annual_prcp,ANNUAL_PRCP,WATER_YEAR_A,k1,y) {
  # returns array of years (length of 100) by weighted sampling of kk nearest neighbors
  kk <- max(round(sqrt(length(ANNUAL_PRCP)),0),round(length(ANNUAL_PRCP)*0.5,0))
  var_order <- 1:length(ANNUAL_PRCP)
  distance <- sqrt((sim_annual_prcp - ANNUAL_PRCP)^2)
  ordered_distances <- matrix(cbind(var_order,distance)[order(distance),],ncol=2)
  K_Distances <- matrix(ordered_distances[1:kk,],ncol=2)
  PROBS <- (1/row(K_Distances)[,1]) / sum((1/row(K_Distances)[,1]))
  set.seed(k1*y)
  selection <- sample(row(K_Distances)[,1],size=100,prob=PROBS,replace=TRUE)
  FINAL_YEARS <- WATER_YEAR_A[K_Distances[selection,1]]
}

GET_PI <- function(m) {
  P <- rbind(m, rep(1, 3))
  P[1,1] <- P[1,1] - 1
  P[2,2] <- P[2,2] - 1
  P[3,3] <- P[3,3] - 1
  B <- c(0,0,0,1)
  PI <- solve(P[c(1,2,4),],B[c(1,2,4)])
  return(PI)
}

adjust_spell <- function(p, dry_change=1, wet_change=1) {
  pnew <- p

  # adjust dry spell
  p01 <- (p['1', '0'] + p['2', '0'])/dry_change - p['2', '0']
  p00 <- p['0', '0'] + (p['1', '0'] - p01)
  pnew['1','0'] <- p01
  pnew['0','0'] <- p00

  #adjust wet spell
  p10 <- (p['0', '1'] + p['2', '1'])/wet_change - p['2', '1']
  p11 <- p['1', '1'] + (p['0', '1'] - p10)
  pnew['0','1'] <- p10
  pnew['1','1'] <- p11

  pnew
}

next_state <- function(prev_state, p) {
  # given previous state and transition probability matrix p
  # simulate next state
  pp1 <- p['0', as.character(prev_state)]
  pp2 <- p['1', as.character(prev_state)] + pp1
  rn <- runif(1)
  if (rn < pp1) {
    next_state <- 0
  } else if (rn >= pp1 & rn < pp2) {
    next_state <-1
  } else {
    next_state <- 2
  }
  next_state
}

wday_range <- function(j, n=7) {
  if (n %% 2 == 1) {
    l <- (n-1)/2
    u <- (n-1)/2
  } else {
    l <- n/2
    u <- n/2 - 1
  }
  rng <- seq(j-l, j+u)
  rng <- ifelse(rng<=0, 365+rng, rng)
  rng <- ifelse(rng>365, rng-365, rng)
  rng
}

KNN <- function(cur_sim_PRCP,cur_sim_TEMP,PRCP_TODAY,TEMP_TODAY,PRCP_TOMORROW,TEMP_TOMORROW,TMAX_TOMORROW,TMIN_TOMORROW,DATE_TOMORROW,k,sd_monthly_PRCP,sd_monthly_TEMP,mean_monthly_PRCP,mean_monthly_TEMP,k1,count) {
  w_PRCP <- 100/sd_monthly_PRCP
  w_TEMP <- 10/sd_monthly_TEMP
  var_order <- 1:length(PRCP_TODAY)
  distance <- sqrt(w_PRCP*((cur_sim_PRCP-mean_monthly_PRCP) - (PRCP_TODAY-mean_monthly_PRCP))^2 + w_TEMP*((cur_sim_TEMP-mean_monthly_TEMP) - (TEMP_TODAY-mean_monthly_TEMP))^2)
  ordered_distances <- matrix(cbind(var_order,distance)[order(distance),],ncol=2)
  K_Distances <- matrix(ordered_distances[1:k,],ncol=2)
  PROBS <- (1/row(K_Distances)[,1]) / sum((1/row(K_Distances)[,1]))
  set.seed(k1*count)
  selection <- sample(row(K_Distances)[,1],size=1,prob=PROBS,replace=TRUE)
  FINAL_PRCP <- PRCP_TOMORROW[K_Distances[selection,1]]
  FINAL_TEMP <- TEMP_TOMORROW[K_Distances[selection,1]]
  FINAL_TMAX <- TMAX_TOMORROW[K_Distances[selection,1]]
  FINAL_TMIN <- TMIN_TOMORROW[K_Distances[selection,1]]
  FINAL_DATE <- DATE_TOMORROW[K_Distances[selection,1]]
  return(c(FINAL_PRCP,FINAL_TEMP,FINAL_TMAX,FINAL_TMIN,FINAL_DATE))
}

# load data ----
data(climate)
clim.site.da <- climate
clim.da <- group_by(climate, DATE) %>%
  summarise(PRCP=mean(PRCP),
            TMIN=mean(TMIN),
            TMAX=mean(TMAX)) %>%
  filter(wyear(DATE)>=1983) %>% # skip WY 1982 which is missing 4 sites
  mutate(TEMP=(TMIN+TMAX)/2)
clim.mon <- group_by(clim.da, DATE=floor_date(DATE, 'month')) %>%
  summarise(PRCP=sum(PRCP),
            TMAX=mean(TMAX),
            TMIN=mean(TMIN),
            TEMP=mean(TEMP))
clim.wyr <- group_by(clim.da, WYEAR=wyear(DATE)) %>%
  summarise(PRCP=sum(PRCP),
            TMAX=mean(TMAX),
            TMIN=mean(TMIN),
            TEMP=mean(TEMP))

# annual sim ----
wgen_wave <- wavelet_analysis(clim.wyr$PRCP, sig.level=0.90, noise.type='white')
signif.periods.idx <- which(wgen_wave$gws > wgen_wave$gws.sig$signif)
signif.periods <- wgen_wave$period[signif.periods.idx]
wgen_wave.comp <- wavelet_components(x=clim.wyr$PRCP,
                                     wt=wgen_wave,
                                     n.periods=2,
                                     sig.periods=signif.periods.idx,
                                     n.comp.periods=c(1,3))
wgen_models <- arima_models(components=wgen_wave.comp)
wgen_sim.warm <- simulate_arimas(models=wgen_models, n=num_year_sim)
sim.wyr <- data.frame(WYEAR=seq(START_WYEAR_SIM, length.out=num_year_sim),
                      PRCP=wgen_sim.warm$sum)

clim.da <- mutate(clim.da,
                  WYEAR=wyear(DATE),
                  YEAR=year(DATE),
                  MONTH=month(DATE),
                  DAY=day(DATE)) %>%
  select(DATE, WYEAR, YEAR, MONTH, DAY, PRCP, TMIN, TMAX, TEMP)
clim.mon <- mutate(clim.mon,
                   WYEAR=wyear(DATE),
                   YEAR=year(DATE),
                   MONTH=month(DATE),
                   DAY=day(DATE)) %>%
  select(DATE, WYEAR, YEAR, MONTH, DAY, PRCP, TMIN, TMAX, TEMP)

# daily sim ----
sim.da <- data.frame(DATE=seq(as.Date(paste(min(sim.wyr$WYEAR)-1, '-10-01', sep='')),
                              as.Date(paste(max(sim.wyr$WYEAR), '-09-30', sep='')),
                              by='day')) %>%
  mutate(MONTH=month(DATE),
         WYEAR=wyear(DATE),
         JDAY=yday(DATE),
         WDAY=JDAY-(273+leap_year(DATE)), # water year day (1=Oct 1)
         WDAY=ifelse(WDAY<=0, 365+WDAY+leap_year(DATE), WDAY),
         STATE=NA,
         TRANSITION=NA,
         PRCP=NA,
         TEMP=NA,
         TMIN=NA,
         TMAX=NA,
         DATE_SAMPLED=NA)
sim.da[1, 'STATE'] <- 0
sim.da[1, 'PRCP'] <- 0
sim.da[1, 'TEMP'] <- 25
sim.da[1, 'TMAX'] <- 30
sim.da[1, 'TMIN'] <- 20

i <- 1
for (i_wyr in seq_along(sim.wyr$WYEAR)) {
  # KNN of Annual Precip
  CUR_YEARS <- KNN_ANNUAL(sim.wyr[i_wyr, 'PRCP'], clim.wyr$PRCP, clim.wyr$WYEAR, k1, i_wyr)

  # Create current daily population
  conditional_selection <- NULL
  for (yy in 1:length(CUR_YEARS)) {
    conditional_selection <- c(conditional_selection, which(clim.da$WYEAR==CUR_YEARS[yy]))
  }

  current <- data.frame(DATE=clim.da$DATE[conditional_selection]) %>%
    left_join(clim.da)

  # Get thresholds
  thresh1 <- 0.3
  thresh2_quantile <- 0.8

  thresholds <- mutate(current, MONTH=ordered(MONTH, levels=month_list)) %>%
    filter(PRCP>thresh1) %>%
    group_by(MONTH) %>%
    summarise(THRESH_2=quantile(PRCP, thresh2_quantile),
              THRESH_MAX=max(PRCP)) %>%
    mutate(THRESH_1=thresh1) %>%
    arrange(MONTH) %>%
    select(MONTH, THRESH_1, THRESH_2, THRESH_MAX)

  # Merge current with thresholds
  current <- current %>%
    left_join(select(thresholds, MONTH, THRESH_1, THRESH_2) %>%
                mutate(MONTH=as.numeric(as.character(MONTH)))) %>%
    mutate(JDAY=yday(DATE),
           WDAY=JDAY-(273+leap_year(DATE)), # water year day (1=Oct 1)
           WDAY=ifelse(WDAY<=0, 365+WDAY+leap_year(DATE), WDAY)) %>%
    select(DATE, MONTH, JDAY, WDAY, PRCP, TEMP, TMAX, TMIN, THRESH_1, THRESH_2) %>%
    mutate(STATE=ifelse(PRCP<=THRESH_1, 0,
                        ifelse(PRCP<=THRESH_2, 1, 2)),
           STATE_NEXT=lead(STATE),
           TRANSITION=ifelse(is.na(STATE_NEXT), NA, paste(STATE, STATE_NEXT, sep='')))

  mc_p <- lapply(month_list, function(m) {
    x <- filter(current, MONTH==m)
    with(x, prop.table(table(STATE_NEXT, STATE), 2))
  })
  names(mc_p) <- as.character(month_list)

  current.mon <- select(current, MONTH, PRCP, TEMP) %>%
    gather(VAR, VALUE, PRCP, TEMP) %>%
    plyr::dlply(c("VAR"), function(x) {
      plyr::dlply(x, c("MONTH"), function(x) {
        c(MEAN=mean(x$VALUE), SD=sd(x$VALUE))
      })
    })

  for (j in 1:(nrow(filter(sim.da, WYEAR==sim.wyr$WYEAR[i_wyr])))) {
    i <- i+1
    if (i > nrow(sim.da)) {
      break
    }
    #     print(paste(i_wyr, i, sim.da$WYEAR[i], sim.da$WDAY[i]))

    # Simulate state transition
    sim.da$STATE[i] <- next_state(sim.da$STATE[i-1], mc_p[[as.character(sim.da$MONTH[i-1])]])
    sim.da$TRANSITION[i-1] <- paste0(sim.da$STATE[i-1], sim.da$STATE[i])

    #     print(paste(sim.da$DATE[i], j, sim.da$WDAY[i]))

    current.sel <- mutate(current, SELECTED=(TRANSITION==sim.da$TRANSITION[i-1] & WDAY %in% wday_range(sim.da$WDAY[i]))) %>%
      filter(!is.na(SELECTED))

    if (sum(current.sel$SELECTED)==0) {
      current.sel <- mutate(current, SELECTED=(TRANSITION==sim.da$TRANSITION[i-1] & WDAY %in% wday_range(sim.da$WDAY[i], n=61))) %>%
        filter(!is.na(SELECTED))
    }
    if (sum(current.sel$SELECTED)==0) {
      current.sel <- mutate(current, SELECTED=(TRANSITION==sim.da$TRANSITION[i-1])) %>%
        filter(!is.na(SELECTED))
    }

    PRCP_TODAY <- (filter(current.sel, SELECTED))$PRCP
    TEMP_TODAY <- (filter(current.sel, SELECTED))$TEMP
    PRCP_TOMORROW <- (mutate(current.sel, PRCP=lead(PRCP)) %>% filter(SELECTED))$PRCP
    TEMP_TOMORROW <- (mutate(current.sel, TEMP=lead(TEMP)) %>% filter(SELECTED))$TEMP
    TMAX_TOMORROW <- (mutate(current.sel, TMAX=lead(TMAX)) %>% filter(SELECTED))$TMAX
    TMIN_TOMORROW <- (mutate(current.sel, TMIN=lead(TMIN)) %>% filter(SELECTED))$TMIN
    DATE_TOMORROW <- (mutate(current.sel, DATE=lead(DATE)) %>% filter(SELECTED))$DATE

    cur_sim_PRCP <- sim.da$PRCP[i-1]
    cur_sim_TEMP <- sim.da$TEMP[i-1]

    k <- round(sqrt(length(PRCP_TODAY)))
    RESULT <- KNN(cur_sim_PRCP,cur_sim_TEMP,
                  PRCP_TODAY,TEMP_TODAY,
                  PRCP_TOMORROW,TEMP_TOMORROW,TMAX_TOMORROW,TMIN_TOMORROW,DATE_TOMORROW,
                  k,
                  current.mon[["PRCP"]][[sim.da[i, 'MONTH']]][['SD']],
                  current.mon[["TEMP"]][[sim.da[i, 'MONTH']]][['SD']],
                  current.mon[["PRCP"]][[sim.da[i, 'MONTH']]][['MEAN']],
                  current.mon[["TEMP"]][[sim.da[i, 'MONTH']]][['MEAN']],
                  k1,i)

    sim.da$PRCP[i] <- RESULT[1]
    sim.da$TEMP[i] <- RESULT[2]
    sim.da$TMAX[i] <- RESULT[3]
    sim.da$TMIN[i] <- RESULT[4]
    sim.da$DATE_SAMPLED[i] <- RESULT[5]
  }
}
sim.da$DATE_SAMPLED <- as.Date(sim.da$DATE_SAMPLED)

select(sim.da, DATE, WYEAR, WDAY, PRCP, TEMP, TMIN, TMAX) %>%
  write.csv(file='sim.csv', row.names=FALSE)
