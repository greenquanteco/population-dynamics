### Trying with Bayesian method ###
# this code seems like it is running for single site
# have to figure out to include urbanization

## loading in necessary packages
# have to make sure using the proper version of jags
# Sys.setenv(JAGS_HOME="C:/Program Files/JAGS/JAGS-4.3.2")
library(RMark)
library(tidyverse)
library(rjags)
library(R2jags)

## call in data sets (original excel sheets in SMT 2026 Grid Data)
cap <- read.csv("capture-events.csv")
day <- read.csv("session-days.csv")
ind <- read.csv("capture-individuals.csv")
site.dat <- read.csv("site-metadata.csv")
cmat <- read.csv("capture-history-matrix.csv")

### Starting with frequentist approach to compare ###

## creating time interval
# changing primary to 4 to test code, should be 5
n_primary <- 4
n_secondary <- 4

build_time_intervals <- function(n_primary, n_secondary) {
  block     <- c(rep(0, n_secondary - 1), 1)
  intervals <- rep(block, n_primary)
  intervals[-length(intervals)]
}

time_intervals <- build_time_intervals(n_primary, n_secondary)
time_intervals

## editing cmat to create data frame in proper format

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


# changes column name of ch_strings to just ch 
# adding tag, site, and urban level as individual covariates
ch_df <- data.frame(
  ch             = ch_strings,
  tag   = cmat$X,
  site = ind$site,
  urban.level =ind$site.type,
  stringsAsFactors = FALSE
)

## now have ch_df in proper format but need to remove issues from hunt
# making copy of ch_df 
ch_df.2 <- ch_df

# need to remove NAs from Hunt tags
# do this by keeping only last 16 characters in ch
# for all rows where site is hunt
#ch_df$ch <- ifelse(
# ch_df$site =="hunt",
# substr(ch_df$ch, nchar(ch_df$ch) - 15, nchar(ch_df$ch)),
#ch_df$ch
#)

# removing first occasion for all sites to test code
ch_df.2$ch <- substr(ch_df.2$ch, nchar(ch_df.2$ch) - 15, nchar(ch_df.2$ch))

# need to remove any individuals with no captures after first session was removed
is_all_zeros <- sapply(ch_df.2, grepl, pattern = "^0+$")
ch_df.2 <- ch_df.2[!apply(is_all_zeros, 1, any), , drop = FALSE]

## trouble shooting

# make sure ch lengths are all 20
nchar(ch_df.2$ch)

## fitting a model with constant survival, constant detection and random temporary emigration
p.dot <- list(formula=~1, share = TRUE)
S.dot <- list(formula=~1)
GammaDoublePrime.dot <- list(formula=~1, 
                             share = TRUE)
model.1 <- mark(data = ch_df.2, 
                model = "Robust",
                time.intervals = time_intervals,
                model.parameters = list(S = S.dot,
                                        GammaDoublePrime = GammaDoublePrime.dot,
                                        p = p.dot),
                threads=2,
                output = FALSE)

# print parameter estimates
model.1$results$real

# print population size estimates
model.1$results$derived$`N Population Size`

## if we want to change to Markovian temporary emigration
p.dot <- list(formula=~1, share = TRUE)
S.dot <- list(formula=~1)
GammaDoublePrime.dot <- list(formula=~1)
GammaPrime.dot <- list(formula=~1)
model.2 <- mark(data = ch_df.2, 
                model = "Robust",
                time.intervals = time_intervals,
                model.parameters = list(S = S.dot,
                                        GammaPrime = GammaPrime.dot,
                                        GammaDoublePrime = GammaDoublePrime.dot,
                                        p = p.dot),
                threads=2,
                output = FALSE)

# print parameter estimates
model.2$results$real

# print population size estimates
model.2$results$derived$`N Population Size`

### Now for the Bayesian approach!! ### 

## this data should be in form of cmat (capture matrix)
# need to move tag number column to side label
rownames(cmat) <- cmat$X
cmat$X <- NULL

## how will this work with NAs from Hunt? **
## look at code for duplicating encounter histories. Is that necessary? **

## need to compute some starting quantities
n.ind <- nrow(cmat) # number of individuals
n.secondary <- c(4, 4, 4, 4, 4) # number of secondary occasions per primary occasion
n.primary <- length(n.secondary) # number of primary occasions
index <- list(1:4,
              5:8,
              9:12,
              13:16,
              17:20) # the secondary occasions

## calculate the number of individuals caught in each primary occasion, 
## which we will need to get an estimate of population size
caught <- rep(NA, n.primary) # creates vector to catch results
for (i in 1:n.primary){ # iterate through each primary session
  tmp <- cmat[,index[[i]]] # iterate through each secondary within each primary
  caught[i] <- nrow(tmp[rowSums(tmp)!=0,]) # sum number of ind captures and put into "caught" vector
}
caught
## need to check these values are correct ^^

## format the data as an array with dimensions the number of individuals
## times the number of primary occasions times the number of secondary occasions
obs <- array(NA, dim = c(n.ind, n.primary, max(n.secondary)))
dim(obs)

# have to convert cmat from data frame to matrix
cmat.matrix <- as.matrix(cmat)
for (i in 1:n.primary){
  obs[,i,1:n.secondary[i]] <- cmat.matrix[,index[[i]]]
 }
dim(obs)


## format the data as required in the Bayesian implementation of the robust design
ch <- matrix(NA, n.ind, n.primary)
for (i in 1:n.ind){
  for (t in 1:n.primary){
    ifelse(any(obs[i,t,1:n.secondary[t]] == 1), ch[i,t] <- 1, ch[i,t] <- 2)
  }
}

## summarize detection by primary and secondary occasions
test <- matrix(NA, n.ind, n.primary)
for (i in 1:nrow(test)){
  for (j in 1:ncol(test)){
    test[i,j] <- sum(obs[i,j,], na.rm = TRUE)
  }
}

seen <- array(NA, c(n.ind, n.primary, max(n.secondary)))
missed <- array(NA, c(n.ind, n.primary, max(n.secondary)))

for (i in 1:nrow(test)){
  for (t in 1:ncol(test)){
    for (j in 1:n.secondary[t]){
      if(test[i,t] > 1 & obs[i,t,j] == 1){seen[i,t,j] <- 1}
      if(test[i,t] >= 1 & obs[i,t,j] == 0){missed[i,t,j] <- 1}
    }
  }
}

yes <- matrix(NA, n.primary, max(n.secondary))
no <- matrix(NA, n.primary, max(n.secondary))

for (i in 1:nrow(yes)){
  for (j in 1:ncol(yes)){
    yes[i,j] <- sum(seen[,i,j], na.rm = TRUE)
    no[i,j] <- sum(missed[,i,j], na.rm = TRUE)
  }
}

total <- yes + no
total

## get first occasions of capture
get.first <- function(x)min(which (x != 2))
first <- apply(ch,1,get.first); first[first == "Inf"] <- NA

## cut individuals released in last primary occasion (?)
ch <- subset(ch, first != n.primary)
first <- subset(first, first != n.primary)

## define initial values for latent states
z.init <- matrix(NA, nrow(ch), ncol(ch))
for (i in 1:nrow(ch)){
  if(first[i] < ncol(z.init)){
    z.init[i,(first[i] + 1):ncol(z.init)] <- 1
  }
}

## using random immigration




