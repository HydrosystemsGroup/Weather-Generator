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

cat('Extracting grid points for the northeast\n')
grid <- tbl(db, "grid") %>%
  filter(LAT>=40) %>%
  collect %>%
  as.data.frame

cat('Extracting datasets for all points for the northeast\n')
df <- tbl(db, "data") %>%
  filter(LAT>=40) %>%
  collect %>%
  mutate(DATE=ymd(paste(YEAR, MONTH, 1, sep='-'))) %>%
  as.data.frame

cat('Saving datasets to rda file:', RDA_PATH, '\n')
maurer <- list(data=df,
               grid=grid)

save(maurer, file=RDA_PATH)
