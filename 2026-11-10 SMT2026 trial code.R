### trying stuff with SMT2026 data ###

## will need to make new cmat with just peles !!
## code for creating cmat is in script "make-capture-history.R"

# call in data sets (original excel sheets in SMT 2026 Grid Data)
cap <- read.csv("capture-events.csv")
day <- read.csv("session-days.csv")
ind <- read.csv("capture-individuals.csv")
site.dat <- read.csv("SMT 2026 Grid Data(smt26-site).csv")

## call in RMark
library(RMark)

## creating time interval
n_primary <- 5
n_secondary <- 4

build_time_intervals <- function(n_primary, n_secondary) {
  block     <- c(rep(0, n_secondary - 1), 1)
  intervals <- rep(block, n_primary)
  intervals[-length(intervals)]
}

time_intervals <- build_time_intervals(n_primary, n_secondary)
time_intervals

# 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0

## organizing code into proper formatting
# uses cmat to create string of capture histories and places strings together in a data frame
ch_strings <- apply(cmat, 1, paste, collapse = "")
ch_df <- data.frame(ch = ch_strings, stringsAsFactors = FALSE)

# need to add urban level to ind so can add to ch_df
cap$site <- ind$site[match(cap$tag, ind$tag)]
ind$site.type <- site.dat$type[match(ind$site, site.dat$site)]

# changes column name of ch_strings to just ch 
# adding tag, site, and urban level as individual covariates
ch_df <- data.frame(
  ch             = ch_strings,
  tag   = ind$tag,
  site = ind$site,
  urban.level =ind$site.type,
  stringsAsFactors = FALSE
)

## trouble shooting

# make sure ch lengths are all 20
nchar(ch_df$ch)

## try and run model
rd_processed <- process.data(
  ch_df,
  model          = "Robust",
  time.intervals = time_intervals
)

rd_ddl <- make.design.data(rd_processed)

rd_fit <- mark(
  rd_processed,
  rd_ddl,
  model.parameters = list(
    S                = list(formula = ~ site.type),
    p                = list(formula = ~ 1),
    GammaPrime       = list(formula = ~ 1),
    GammaDoublePrime = list(formula = ~ 1)
  ),
  output = FALSE,
  silent = TRUE
)