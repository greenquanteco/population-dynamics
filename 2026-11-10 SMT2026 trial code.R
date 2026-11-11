### trying stuff with SMT2026 data ###

## will need to make new cmat with just peles !!
## code for creating cmat is in script "make-capture-history.R"

# call in data sets (original excel sheets in SMT 2026 Grid Data)
cap <- read.csv("capture-events.csv")
day <- read.csv("session-days.csv")
ind <- read.csv("capture-individuals.csv")
site.dat <- read.csv("site-metadata.csv")
cmat <- read.csv("capture-history-matrix.csv")

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

# need to create an individual-specific time interval matrix 
# because hunt has different time interval lengths

##### need help on this :(

## organizing code into proper formatting

# need to delete tag column for creating character string
# but want to also keep original
cmat.str <- cmat
cmat.str$X <- NULL

# delete neg tags from ind
ind <- ind[ind$tag > 0,]

# remove tags 1055 and 1082 because 
#they are marked as two different species
# in original data sheets
ind <- ind[!ind$tag %in% c(1055, 1082),]

# adding site to cap
cap$site <- ind$site[match(cap$tag, ind$tag)]

# remove a possible duplicate, based on
# coin flip which entry from that day to keep
cap <- cap[-which(cap$date=="2026-07-22" & cap$tag==1065)[1],]

# remove any rows where site column is NA (issue for if individual does not have
# tag number to link to site)
cap <- cap[!is.na(cap$site),]

# remove entry for tag 1033 because it was recorded as different site (noon) than its actual
# site location (field)
cap <- cap[-which(cap$tag==1033 & cap$date=="2026-07-28"),]

# need to add urban level to ind so can add to ch_df
ind$site.type <- site.dat$type[match(ind$site, site.dat$site)]

# sort ind by tag number so it matches order in cmat
ind <- ind[order(ind$tag),]

# uses cmat to create string of capture histories and places strings together in a data frame
ch_strings <- apply(cmat.str, 1, paste, collapse = "")
ch_df <- data.frame(ch = ch_strings, stringsAsFactors = FALSE)

# need to remove NAs from Hunt tags
# do this by keeping only last 16 characters in ch
# for all rows where site is hunt
ch_df$ch <- ifelse(
  ch_df$site == "hunt",
  substr(ch_df$ch, nchar(ch_df$ch) - 15, nchar(ch_df$ch)),
  ch_df$ch
)

# changes column name of ch_strings to just ch 
# adding tag, site, and urban level as individual covariates
ch_df <- data.frame(
  ch             = ch_strings,
  tag   = cmat$X,
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

