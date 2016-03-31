## ----libraries, message=FALSE--------------------------------------------
library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(zoo)
library(moments)
library(weathergen)
theme_set(theme_bw())

# set random seed for reproducibility
set.seed(1)

## ----load-climate-cities-------------------------------------------------
data(climate_cities)
str(climate_cities)

## ----load-obs------------------------------------------------------------
# extract dataset for boston
obs_day <- climate_cities[['boston']]

# subset complete water years using complete_water_years() function
obs_day <- obs_day[complete_water_years(obs_day$DATE),]

obs_day <- mutate(obs_day,
                  WYEAR=wyear(DATE, start_month=10), # extract water year
                  MONTH=month(DATE),                 # extract month
                  TEMP=(TMIN+TMAX)/2) %>%            # compute mean temperature
  select(WYEAR, MONTH, DATE, PRCP, TEMP, TMIN, TMAX, WIND)

# compute annual timeseries by water year
obs_wyr <- group_by(obs_day, WYEAR) %>%
  summarise(PRCP=sum(PRCP),
            TEMP=mean(TEMP),
            TMIN=mean(TMIN),
            TMAX=mean(TMAX),
            WIND=mean(WIND))

# save the daily and annual timeseries to a list
obs <- list(day=obs_day, wyr=obs_wyr)

# summary of the daily timeseries
summary(obs[['day']])

## ----plot-obs-da, fig.height=8, fig.width=10-----------------------------
gather(obs[['day']], VAR, VALUE, PRCP:WIND) %>%
  ggplot(aes(DATE, VALUE)) +
  geom_line() +
  facet_wrap(~VAR, scales='free_y', ncol=1) +
  labs(x='', y='', title='Historical Daily Weather Data for Boston, MA')

## ----plot-obs-wyr, fig.height=8, fig.width=10----------------------------
gather(obs[['wyr']], VAR, VALUE, PRCP:WIND) %>%
  ggplot(aes(WYEAR, VALUE)) +
    geom_line() +
    facet_wrap(~VAR, scales='free_y', ncol=1) +
    labs(x='', y='', title='Historical Annual Weather Data for Boston, MA')

## ----create-zoo----------------------------------------------------------
zoo_day <- zoo(x = obs[['day']][, c('PRCP', 'TEMP', 'TMIN', 'TMAX', 'WIND')], 
             order.by = obs[['day']][['DATE']])
class(zoo_day)
summary(zoo_day)

## ----wgen-daily, warning=FALSE-------------------------------------------
sim <- wgen_daily(zoo_day, 
                  n_year=2, 
                  start_month=10, 
                  start_water_year=2000, 
                  include_leap_days=FALSE,
                  n_knn_annual=100,
                  dry_wet_threshold=0.3, 
                  wet_extreme_quantile_threshold=0.8,
                  adjust_annual_precip=TRUE,
                  annual_precip_adjust_limits=c(0.9, 1.1),
                  dry_spell_changes=1,
                  wet_spell_changes=1,
                  prcp_mean_changes=1,
                  prcp_cv_changes=1,
                  temp_mean_changes=0)

## ----str-sim-------------------------------------------------------------
str(sim, max=1)

## ----sim-obs-------------------------------------------------------------
head(sim[['obs']])

## ----sim-state_thresholds------------------------------------------------
sim[['state_thresholds']]

## ----sim-transition_matrices---------------------------------------------
sim[['transition_matrices']][[6]]

## ----sim-state_equilibrium-----------------------------------------------
sim[['state_equilibria']][[6]]

## ----sim-change_factors--------------------------------------------------
str(sim[['change_factors']])

## ----sim-out-------------------------------------------------------------
str(sim[['out']])

## ----sim-out-plot, fig.width=6, fig.height=8-----------------------------
select(sim[['out']], DATE, PRCP, TEMP, TMIN, TMAX, WIND) %>%
  gather(VAR, VALUE, -DATE) %>%
  ggplot(aes(DATE, VALUE)) +
  geom_line() +
  facet_wrap(~VAR, ncol=1, scales='free_y') +
  labs(x='Simulation Date', y='Simulated Value', title='Simulated Timeseries') +
  theme_bw()

## ----sim-annual----------------------------------------------------------
zoo_precip_wyr <- zoo(obs[['wyr']][['PRCP']], order.by = obs[['wyr']][['WYEAR']])
sim_annual_prcp <- sim_annual_arima(x = zoo_precip_wyr, start_year = 2000, n_year = 40)

## ----sim-annual-str------------------------------------------------------
str(sim_annual_prcp, max.level=1)

## ----sim-annual-prcp-model-----------------------------------------------
summary(sim_annual_prcp[['model']])

## ----sim-annual-prcp-out-------------------------------------------------
sim_annual_prcp[['out']]

## ----plot-sim-annual, fig.width=6, fig.height=3--------------------------
data.frame(WYEAR=time(sim_annual_prcp[['out']]),
           PRCP=zoo::coredata(sim_annual_prcp[['out']])) %>%
  ggplot(aes(WYEAR, PRCP)) +
  geom_line() +
  geom_hline(yintercept=mean(obs[['wyr']][['PRCP']]), color='red') +
  geom_text(aes(label=TEXT), data=data.frame(WYEAR=2000, PRCP=mean(obs[['wyr']][['PRCP']]), TEXT="Historical Mean"),
            hjust=0, vjust=-1, color='red') +
  labs(x="Simulation Water Year", y="Annual Precip (mm/yr)", title="Simulated Annual Precipitation")

## ----annual-mc-----------------------------------------------------------
# batch simulation to compare mean/sd/skew
sim_wyr_batch <- lapply(1:50, function(i) { 
    sim_annual_arima(x = obs[['wyr']][['PRCP']], start_year = 2000, n_year = 20)[['out']]
  }) 
names(sim_wyr_batch) <- paste0(seq(1, length(sim_wyr_batch)))
sim_wyr_batch <- do.call(merge, sim_wyr_batch) %>% 
  apply(MARGIN=2, FUN=function(x) { c(mean=mean(x), sd=sd(x), skew=skewness(x))}) %>% 
  t
sim_wyr_stats <- rbind(t(colMeans(sim_wyr_batch)),
                       c(mean=mean(obs[['wyr']][['PRCP']]), 
                         sd=sd(obs[['wyr']][['PRCP']]), 
                         skew=skewness(obs[['wyr']][['PRCP']])))
sim_wyr_stats <- as.data.frame(sim_wyr_stats)
sim_wyr_stats$dataset <- c('obs', 'sim')
sim_wyr_stats <- gather(sim_wyr_stats, stat, value, mean, sd, skew)

## ----plot-annual-mc, fig.width=6, fig.height=4---------------------------
as.data.frame(sim_wyr_batch) %>% 
  mutate(trial=row_number()) %>% 
  gather(stat, value, mean:skew) %>% 
  ggplot(aes(x=stat, value)) + 
  geom_boxplot() + 
  geom_point(aes(x=stat, y=value, color=dataset), data=sim_wyr_stats, size=3) +
  scale_color_manual('', values=c('sim'='deepskyblue', 'obs'='orangered'),
                     labels=c('sim'='Mean of Simulated Trials', 'obs'='Mean of Historical')) +
  labs(x='', y='') +
  facet_wrap(~stat, scales='free')

## ----knn_annual----------------------------------------------------------
n_knn_annual <- 100

# generate vector of sampled years
sampled_years <- knn_annual(prcp=coredata(sim_annual_prcp[['out']])[1],
                        obs_prcp=zoo_precip_wyr,
                        n=n_knn_annual)
sampled_years

## ----plot-knn_annual, fig.width=8, fig.height=4--------------------------
ggplot() +
  geom_point(aes(WYEAR, PRCP, color=SAMPLED), data=obs[['wyr']] %>% mutate(SAMPLED=WYEAR %in% sampled_years)) +
  geom_hline(yintercept=sim_annual_prcp[['out']][1], linetype=2, color='red') +
  geom_text(aes(x=x, y=y, label=label),
            data=data.frame(x=min(obs[['wyr']][['WYEAR']]),
                            y=sim_annual_prcp[['out']][1],
                            label='Target Simulated Annual Precip'),
            vjust=-0.5, hjust=0) +
  scale_color_manual('Sampled in KNN', values=c('TRUE'='orangered', 'FALSE'='grey30')) +
  labs(x='Water Year', y='Annual Precipitation (mm/yr)')

## ----plot-knn_annual-days, fig.width=8, fig.height=4---------------------
# loop through population years and extract historical daily values
sampled_days <- lapply(sampled_years, function(yr) {
  obs[['day']][which(obs[['day']][['WYEAR']]==yr), ]
}) %>%
  rbind_all() %>%
  mutate(SAMPLE_INDEX=row_number())

# plot combined daily timeseries based on the KNN sampled years
ggplot(sampled_days, aes(SAMPLE_INDEX, PRCP, color=WYEAR)) +
  geom_point(size=1) +
  labs(x='Index Day', y='Daily Precip (mm/day)', title='Synthetic Daily Timeseries from KNN Sampled Years')

## ----mc-state-threshold--------------------------------------------------
thresh <- mc_state_threshold(sampled_days[['PRCP']], sampled_days[['MONTH']],
                             dry_wet_threshold=0.3, wet_extreme_quantile_threshold=0.8)

thresh_df <- do.call(rbind, thresh) %>%
  as.data.frame() %>%
  mutate(month=row_number()) %>%
  select(month, dry_wet, wet_extreme)
thresh_df

## ----plot-mc-state-threshold---------------------------------------------
gather(thresh_df, var, value, -month) %>%
  ggplot(aes(factor(month), value, fill=var)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_manual('Threshold', values=c('dry_wet'='orangered', 'wet_extreme'='deepskyblue'),
                    labels=c('dry_wet'='Dry-Wet', 'wet_extreme'='Wet-Extreme')) +
  labs(x='Month', y='Threshold Amount (mm/day)', title='Markov State Thresholds by Month')

## ----assign-states-------------------------------------------------------
sampled_days$STATE <- mc_assign_states(sampled_days$PRCP, sampled_days$MONTH, c('d', 'w', 'e'), thresh)
head(sampled_days)

## ----plot-states, fig.width=6, fig.height=4------------------------------
filter(sampled_days, SAMPLE_INDEX <= 365*5) %>%
ggplot(aes(SAMPLE_INDEX, PRCP, color=STATE)) +
  geom_point(size=2) +
  scale_color_manual('Markov State',
                     values=c('d'='orangered', 'w'='chartreuse3', 'e'='deepskyblue'),
                     labels=c('d'='Dry', 'w'='Wet', 'e'='Extreme')) +
  ylim(0, 50) +
  labs(x='Sample Date Index', y='Daily Precipitation (mm/day)', title='Markov State of Synthetic Daily Timeseries')

## ----fit-transitions-----------------------------------------------------
transition_matrices <- mc_fit(states=sampled_days[['STATE']], months=sampled_days[['MONTH']])

## ----transition-june-----------------------------------------------------
transition_matrices[[6]]

## ----equil-june----------------------------------------------------------
mc_state_equilibrium(transition_matrices[[6]])

## ----sim-days------------------------------------------------------------
sampled_days <- dplyr::mutate(sampled_days,
                              WDAY=waterday(DATE, start_month=10),
                              STATE_PREV=lag(STATE),
                              PRCP_PREV=lag(PRCP),
                              TEMP_PREV=lag(TEMP),
                              TMAX_PREV=lag(TMAX),
                              TMIN_PREV=lag(TMIN),
                              WIND_PREV=lag(WIND))

sim_days <- sim_mc_knn_day(x=sampled_days, n_year=1, states=c('d', 'w', 'e'), transitions=transition_matrices,
                           start_month=10, start_water_year=2000, include_leap_days=FALSE)

str(sim_days, max=1)

## ----plot-sim-days, fig.width=8, fig.height=8----------------------------
select(sim_days, DATE, PRCP:WIND) %>%
  gather(VAR, VALUE, -DATE) %>%
  ggplot(aes(DATE, VALUE)) +
  geom_line() +
  facet_wrap(~VAR, scales='free_y', ncol=1) +
  labs(x='Simulate Date', y='Value', title='Daily Simulation for One Year')

## ----adjust-prcp---------------------------------------------------------
orig_prcp <- sim_days[['PRCP']]
adj_prcp <- adjust_daily_to_annual(values_day=sim_days[['PRCP']], 
                                   years_day=wyear(sim_days[['DATE']]), 
                                   values_yr=coredata(sim_annual_prcp[['out']])[1], 
                                   years_yr=time(sim_annual_prcp[['out']])[1], 
                                   min_ratio=0.5, max_ratio=1.5)
c(Original=sum(orig_prcp), Adjusted=sum(adj_prcp$adjusted), Annual=coredata(sim_annual_prcp[['out']])[1])

## ----plot-adj-prcp, fig.width=8, fig.height=4----------------------------
data.frame(DATE=sim_days[['DATE']],
           ORIGINAL=orig_prcp,
           ADJUSTED=adj_prcp$adjusted) %>%
  gather(VAR, VALUE, -DATE) %>%
  ggplot(aes(DATE, VALUE, color=VAR, linetype=VAR)) +
  geom_line() +
  scale_color_manual('', values=c('ORIGINAL'='black', 'ADJUSTED'='orangered')) +
  scale_linetype_manual('', values=c('ORIGINAL'=1, 'ADJUSTED'=2)) +
  labs(x='Simulation Date', y='Daily Precipitation (mm/day)')

## ----adj-sim, warning=FALSE----------------------------------------------
sim <- wgen_daily(zoo_day, 
                  n_year=5, 
                  start_month=10, 
                  start_water_year=2000, 
                  include_leap_days=FALSE,
                  n_knn_annual=100,
                  dry_wet_threshold=0.3, 
                  wet_extreme_quantile_threshold=0.8,
                  adjust_annual_precip=TRUE,
                  annual_precip_adjust_limits=c(0.5, 1.5),
                  dry_spell_changes=1,
                  wet_spell_changes=1,
                  prcp_mean_changes=1,
                  prcp_cv_changes=1,
                  temp_mean_changes=0)

## ----adj-prcp, warning=FALSE---------------------------------------------
prcp.unadj <- sim[['out']][, c('DATE','MONTH','PRCP')]
prcp.adj.mean <- adjust_daily_gamma(prcp.unadj[['PRCP']], prcp.unadj[['MONTH']],
                                    mean_change=1.2, cv_change=1)
prcp.adj.cv <- adjust_daily_gamma(prcp.unadj[['PRCP']], prcp.unadj[['MONTH']],
                                  mean_change=1, cv_change=1.2)
prcp.adj.both <- adjust_daily_gamma(prcp.unadj[['PRCP']], prcp.unadj[['MONTH']],
                                    mean_change=1.2, cv_change=1.2)
prcp.df <- data.frame(DATE=prcp.unadj[['DATE']], 
                      UNADJESTED=prcp.unadj[['PRCP']],
                      ADJUST_MEAN=prcp.adj.mean[['adjusted']], 
                      ADJUST_CV=prcp.adj.cv[['adjusted']], 
                      ADJUST_BOTH=prcp.adj.both[['adjusted']]) %>%
  gather(GROUP, VALUE, -DATE)

## ----plot-prcp-adj, fig.width=8, fig.height=6----------------------------
ggplot(prcp.df, aes(DATE, VALUE)) +
  geom_line() +
  facet_wrap(~GROUP) +
  labs(x='', y='Daily Precipitation (mm/day)', title='Unadjusted and Adjusted Daily Precipitation')

## ----plot-prcp-adj-wyr, fig.width=6, fig.height=4------------------------
mutate(prcp.df, WYEAR=wyear(DATE, start_month=10)) %>%
  group_by(WYEAR, GROUP) %>%
  summarise(VALUE=sum(VALUE)) %>%
  ggplot(aes(factor(WYEAR), VALUE, fill=GROUP)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_discrete('') +
  labs(x='Water Year', y='Annual Precipitation (mm/yr)', title='Annual Precipitation With Adjustments')

## ----plot-prcp-adj-cdf, fig.width=6, fig.height=4------------------------
prcp.df %>%
  filter(VALUE>0) %>%
  group_by(GROUP) %>%
  arrange(VALUE) %>%
  mutate(ROW=row_number(),
         PROB=(ROW-0.5)/n()) %>%
  ggplot(aes(PROB, VALUE, color=GROUP)) +
  geom_line() +
  ylim(0, 60) +
  scale_color_discrete('') +
  labs(y='Daily Precipitation (mm/day)', x='Non-Exceedence Probability', title='Distributions of Adjusted Precipitation Timeseries')

## ----temp-adj------------------------------------------------------------
temp.unadj <- sim[['out']][, c('DATE','MONTH','TEMP')]
temp.adj <- adjust_daily_additive(temp.unadj[['TEMP']], temp.unadj[['MONTH']],
                                  mean_change=5)
temp.df <- data.frame(DATE=temp.unadj[['DATE']], 
                      UNADJESTED=temp.unadj[['TEMP']],
                      ADJUSTED=temp.adj[['adjusted']]) %>%
  gather(GROUP, VALUE, -DATE)

## ----plot-temp-adj, fig.width=6, fig.height=4----------------------------
ggplot(temp.df, aes(DATE, VALUE, color=GROUP)) +
  geom_line() +
  scale_color_discrete('') +
  labs(x='', y='Daily Mean Temperature (degC)', title='Unadjusted and Adjusted Daily Mean Temperature')

## ----plot-temp-adj-wyr, fig.width=6, fig.height=4------------------------
mutate(temp.df, WYEAR=wyear(DATE, start_month = 10)) %>%
  group_by(WYEAR, GROUP) %>%
  summarise(VALUE=mean(VALUE)) %>%
  ggplot(aes(factor(WYEAR), VALUE, fill=GROUP)) +
  geom_bar(stat='identity', position='dodge') +
  scale_fill_discrete('') +
  labs(x='', y='Annual Mean Temperature (degC)', title='Unadjusted and Adjusted Annual Mean Temperature')

