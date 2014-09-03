# creates an rda file containing regional and local monthly timeseries
# from the Maurer et al (2002) dataset

# get and check command line arguments
argv <- commandArgs(TRUE)

if (length(argv) < 2) {
  stop('Missing directory arguments, expecting:\n$ RScript make_rdata.R <input_db> <output_rda>')
}

DB_PATH <- argv[1]
RDA_PATH <- argv[2]

if (!file.exists(DB_PATH)) {
  stop(paste0("Input database does not exist: ", DB_PATH))
}

# load libraries
library(dplyr)
library(lubridate)

# connect to sqlite database
cat(paste0(rep('-', 80), collapse=''), '\n')
cat('Connecting to database:',DB_PATH,'\n')
db <- src_sqlite(DB_PATH, create = FALSE)

# extract grid points
grid <- tbl(db, "grid") %>%
  collect %>%
  as.data.frame

# function to compute regional average
compute.regional <- function(name, lat_limits) {
  reg.df <- tbl(db, "data") %>%
    filter(LAT>=lat_limits[1], LAT<lat_limits[2]) %>%
    group_by(YEAR, MONTH) %>%
    summarise(PRCP=mean(PRCP),
              TMAX=mean(TMAX),
              TMIN=mean(TMIN),
              WIND=mean(WIND)) %>%
    collect %>%
    mutate(DATE=ymd(paste(YEAR, MONTH, 1, sep='-'))) %>%
    as.data.frame
  reg.grid <- grid %>%
    filter(LAT>=lat_limits[1], LAT<=lat_limits[2])
  return(list(name=name,
              lat_limits=name,
              data=reg.df,
              grid=reg.grid))
}

# compute regional average monthly timeseries
cat('Computing regional average monthly timeseries for lower northeast\n')
reg.ne.lower <- compute.regional(name='Northeast-Lower', lat_limits=c(40, 43))

cat('Computing regional average monthly timeseries for upper northeast\n')
reg.ne.upper <- compute.regional(name='Northeast-Upper', lat_limits=c(43, 49))

cat('Computing regional average monthly timeseries for all of northeast\n')
reg.ne <- compute.regional(name='Northeast', lat_limits=c(40, 49))

cat('Extracting monthly timeseries for 100 random locations in northeast\n')
# choose 100 random locations
set.seed(1234)
locs <- select(reg.ne[['grid']], LAT, LON) %>%
  unique %>%
  sample_n(size=100)

# extract monthly timeseries for random locations
# stored as unnamed list of data frames
local <- apply(locs, 1, function(loc) {
  df.loc <- tbl(db, "data") %>%
    filter(LAT==loc[['LAT']], LON==loc[['LON']]) %>%
    collect %>%
    mutate(DATE=ymd(paste(YEAR, MONTH, 1, sep='-'))) %>%
    as.data.frame
  return(df.loc)
}) %>%
  unname

local <- tbl(db, "data") %>%
  filter(LAT>=40) %>%
  collect %>%
  mutate(DATE=ymd(paste(YEAR, MONTH, 1, sep='-'))) %>%
  as.data.frame

cat('Saving regional and local timeseries to rda file:', RDA_PATH, '\n')
maurer <- list(regional=list('ne-lower'=reg.ne.lower,
                             'ne-upper'=reg.ne.upper,
                             'ne'=reg.ne),
               local=local,
               grid=grid)

save(maurer, file=RDA_PATH)
