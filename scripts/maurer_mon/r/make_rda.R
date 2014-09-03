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

# set lat limits
LAT_LIMITS <- c(40, 90)

# load libraries
library(dplyr)
library(lubridate)

# connect to sqlite database
cat(paste0(rep('-', 80), collapse=''), '\n')
cat('Connecting to database:',DB_PATH,'\n')
db <- src_sqlite(DB_PATH, create = FALSE)

# extract grid points
grid <- tbl(db, "grid") %>%
  filter(LAT>=LAT_LIMITS[1], LAT<=LAT_LIMITS[2]) %>%
  collect %>%
  as.data.frame

# compute regional average monthly timeseries
cat('Computing regional average monthly timeseries\n')
regional <- tbl(db, "data") %>%
  filter(LAT>=LAT_LIMITS[1], LAT<=LAT_LIMITS[2]) %>%
  group_by(YEAR, MONTH) %>%
  summarise(PRCP=mean(PRCP),
            TMAX=mean(TMAX),
            TMIN=mean(TMIN),
            WIND=mean(WIND)) %>%
  collect %>%
  mutate(DATE=ymd(paste(YEAR, MONTH, 1, sep='-'))) %>%
  as.data.frame


cat('Extracting monthly timeseries for 100 random locations\n')
# choose 100 random locations
set.seed(1234)
locs <- select(grid, LAT, LON) %>%
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

cat('Saving regional and local timeseries to rda file:', RDA_PATH, '\n')
maurer <- list(regional=regional,
               local=local,
               grid=grid)

save(maurer, file=RDA_PATH)
