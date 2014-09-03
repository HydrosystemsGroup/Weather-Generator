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

# compute regional average monthly timeseries
cat('Computing regional average monthly timeseries\n')
regional <- tbl(db, "data") %>%
  group_by(YEAR, MONTH) %>%
  summarise(PRCP=mean(PRCP),
            TMAX=mean(TMAX),
            TMIN=mean(TMIN),
            WIND=mean(WIND)) %>%
  collect %>%
  mutate(DATE=ymd(paste(YEAR, MONTH, 1, sep='-')))


cat('Extracting monthly timeseries for 10 random locations\n')
# choose 10 random locations
set.seed(1234)
locs <- tbl(db, "data") %>%
  select(LAT, LON) %>%
  collect %>%
  unique %>%
  sample_n(size=10)

# extract monthly timeseries for random locations
# stored as unnamed list of data frames
local <- apply(locs, 1, function(loc) {
  df.loc <- tbl(db, "data") %>%
    filter(LAT==loc[['LAT']], LON==loc[['LON']]) %>%
    collect %>%
    mutate(DATE=ymd(paste(YEAR, MONTH, 1, sep='-')))
  return(df.loc)
})

cat('Saving regional and local timeseries to rda file:', RDA_PATH, '\n')
maurer <- list(regional=regional,
               local=local)

save(maurer, file=RDA_PATH)
