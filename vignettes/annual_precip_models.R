## ----load data, warning=FALSE, message=FALSE-----------------------------
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
theme_set(theme_bw())
library(weatherGen)
data(climate)

# compute daily averages across sites
clim.da <- group_by(climate, DATE) %>%
  summarise(N=n(),
            PRCP=mean(PRCP),
            TMIN=mean(TMIN),
            TMAX=mean(TMAX))

# aggregate by water year
clim.wyr <- clim.da %>%
  mutate(WYEAR=ifelse(month(DATE)>=10, year(DATE)+1, year(DATE))) %>%
  group_by(WYEAR) %>%
  summarise(N=n(),
            PRCP=sum(PRCP),
            TMIN=mean(TMIN),
            TMAX=mean(TMAX))

## ----plot annual data, fig.cap="Annual Timeseries"-----------------------
ggplot(clim.wyr, aes(WYEAR, PRCP)) +
  geom_line() +
  labs(x="Water Year", y="Annual Precipitation (mm)")

