library(dplyr)
library(RPostgreSQL)
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, dbname = "cst")

locations <- list(boston=c(42.358056, -71.063611),
                  providence=c(41.823611, -71.422222),
                  worcester=c(42.266667, -71.8))

climate_cities <- lapply(locations, function(loc) {
  df <- dbGetQuery(con, statement=paste(
    "SELECT date, prcp, tmax, tmin, wind",
    "FROM location_data(",
    loc[1], ",", loc[2],
    ")",
    "ORDER BY date"))
  names(df) <- toupper(names(df))
  df
})

lapply(climate_cities, head)

save(climate_cities, file='../../data/climate_cities.rda')
