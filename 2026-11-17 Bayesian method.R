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

 #######################################################################################
### Now for the Bayesian approach!! ### 

## this data should be in form of cmat (capture matrix)
# need to move tag number column to side label
rownames(cmat) <- cmat$X
cmat$X <- NULL

## how will this work with NAs from Hunt? **
## look at code for duplicating encounter histories. Is that necessary? **

## need to compute some starting quantities

# separates out timing of primary vs secondary occasions
n.ind <- nrow(cmat) # number of individuals
n.secondary <- c(4, 4, 4, 4, 4) # number of secondary occasions per primary occasion
n.primary <- length(n.secondary) # number of primary occasions
index <- list(1:4,
              5:8,
              9:12,
              13:16,
              17:20) # the secondary occasions

## calculate the number of individuals caught in each primary occasion, 
## which we will need to get an estimate of population size ##
caught <- rep(NA, n.primary) # creates vector to catch results
for (i in 1:n.primary){ # iterate through each primary session
  tmp <- cmat[,index[[i]]] # iterate through each secondary within each primary
  caught[i] <- nrow(tmp[rowSums(tmp)!=0,]) # sum number of ind captures and put into "caught" vector
}
caught # prints vector

# format the data as an array with dimensions the number of individuals
# times the number of primary occasions times the number of secondary occasions
obs <- array(NA, dim = c(n.ind, n.primary, max(n.secondary)))
dim(obs)

# have to convert cmat from data frame to matrix
cmat.matrix <- as.matrix(cmat)

# input values from cmat into obs array format
# this is creating a 3D array which holds the capture histories for all individuals
for (i in 1:n.primary){
  obs[,i,1:n.secondary[i]] <- cmat.matrix[,index[[i]]]
 }
dim(obs)
# dim = rows, columns, layers
# so there are 219 rows (for each individual), 
# 5 columns (for each primary session), and 4 layers (for each secondary session)
# moving forward: i = ind, t = primary, j = secondary for "obs:

## format the data as required in the Bayesian implementation of the robust design
ch <- matrix(NA, n.ind, n.primary) # creates an empty matrix with rows = num inds and columns = num primary sessions

# using for loop to fill in empty matrix with 1s and 2s
# if in the obs array, any values equaling 1 means the animal was observed at 
# least once during that primary ("any" returns TRUE)
# any values equaling 2 means that animal was not observed during that primary ("any" returns FALSE)
# for Hunt primary 1, all are returning NAs
for (i in 1:n.ind){ # iterating through every individual
  for (t in 1:n.primary){ # iterating through each primary session (and each secondary within)
    ifelse(any(obs[i,t,1:n.secondary[t]] == 1), ch[i,t] <- 1, ch[i,t] <- 2) # checking if the ind was observed
  }
}

## summarize detection by primary and secondary occasions

# creates another empty matrix with rows = num inds and columns = num primary sessions
test <- matrix(NA, n.ind, n.primary)

# filling in test matrix with actual number of times each individual was caught
# within each primary session
# issue! Counting Hunt NAs as 0s for detection. Will this cause issues even if NAS for ch (primary sessions)
for (i in 1:nrow(test)){
  for (j in 1:ncol(test)){
    test[i,j] <- sum(obs[i,j,], na.rm = TRUE) # sums number of detections within each primary session
  }
}

# created two empty arrays with same dimensions as obs (219x5x4) to catch data
seen <- array(NA, c(n.ind, n.primary, max(n.secondary)))
missed <- array(NA, c(n.ind, n.primary, max(n.secondary)))

# iterating through "test" matrix. If individual was detected more than once at a given primary session
# and was seen at the secondary within that primary, then 1 is assigned to "seen" array
# If individual was detected one or more times at a given primary session and was NOT seen 
# at a secondary within that primary, then 1 is assigned to "missed" array
for (i in 1:nrow(test)){
  for (t in 1:ncol(test)){
    for (j in 1:n.secondary[t]){
      if(test[i,t] > 1 & obs[i,t,j] == 1){seen[i,t,j] <- 1}
      if(test[i,t] >= 1 & obs[i,t,j] == 0){missed[i,t,j] <- 1}
    }
  }
}

# creating two empty matrices which represent 
# rows = number of primary sessions and columns = number of secondary sessions
yes <- matrix(NA, n.primary, max(n.secondary))
no <- matrix(NA, n.primary, max(n.secondary))

# yes sums number of 1s within "seen" array
# no sums number of 1s within "missed" array
for (i in 1:nrow(yes)){
  for (j in 1:ncol(yes)){
    yes[i,j] <- sum(seen[,i,j], na.rm = TRUE)
    no[i,j] <- sum(missed[,i,j], na.rm = TRUE)
  }
}

# totaling yes and no. This is important because can assume if animal was seen within later secondary sessions
# but was not seen in other secondary sessions, it was likely there but just not captured
# so totaling missed and seen gives a more actual estimate of actual population size
# at each given occasion
total <- yes + no
total

## now we have an estimate of pop size at each occasion! ##

## Next step: get first occasions of capture ##

# creates a function which takes input "x" then iterates through values to find those 
# NOT equal to 2 (which in ch means not detected) then "min" finds the smallest
# position which is TRUE (meaning the first time the animal is observed)
get.first <- function(x)min(which (x != 2))

# applies "get.first" function to ch to determine what was the first primary session
# where the individual was seen
# ";" separates multiple commands on a single line
first <- apply(ch,1,get.first); first[first == "Inf"] <- NA

## cut individuals released in last primary occasion (?) ##
# this cut 17 individuals when subsetting ch. Condition to keep is value is in first
# is not equal to total number number of primary sampling occasions.
# meaning, any animal first captured in the last session is removed
ch <- subset(ch, first != n.primary)

# does same thing but for "first"
first <- subset(first, first != n.primary)

## define initial values for latent states (true biological condition of animal)
# such as dead or alive or in a specific location
z.init <- matrix(NA, nrow(ch), ncol(ch)) # creating empty matrix with row = ind and col = primary sessions

# code is modifying matrix "z.init"
# iterates through all ind and if first observation was before final sampling period 
# (determined by columns in z.init)
# then replace NA with 1 for every sampling period after first observation
for (i in 1:nrow(ch)){ # iterating through all individuals
  if(first[i] < ncol(z.init)){ # checking whether the first observation occurred before the final sampling period
    z.init[i,(first[i] + 1):ncol(z.init)] <- 1 #replace NA with 1 for every sampling period after first observation
  }
}

## now have initial states for all individuals!! ##

### using random immigration and running Bayesian model ###
model <- function() {
  
  # priors
  phi ~ dunif(0,1)     # survival; dunif calculates PDF of a cont uniform distribution
  gamma ~ dunif(0,1)   # gamma; temporary emigration probability
  gam <- 1 - gamma     # MARK parameterization
  mean.p ~ dunif(0,1)  # detection probability
  
  # secondary occasions p's
  for (t in 1:n.years){
    for (j in 1:max(n.sec[1:n.years])){
      p[t,j] <- mean.p
    }
  }
  
  for (t in 1:n.years){
    for (j in 1:n.sec[t]){
      yes[t,j] ~ dbin(p[t,j], total[t,j])
    }
  }   
  
  # primary occasions p's or pooled detection probability
  for (t in 1:n.years){
    pstar[t] <- 1 - prod(1 - p[t,1:n.sec[t]])
  }
  
  # likelihood
  
  for (i in 1:n.ind){
    z[i,first[i]] <- ch[i,first[i]]
    for (t in (first[i]+1):n.years){
      mu1[i,t] <- z[i,t-1] * phi
      mu2[i,t] <- z[i,t] * (gamma) * pstar[t]
      z[i,t] ~ dbern(mu1[i,t])  # state equations
      ch[i,t] ~ dbern(mu2[i,t]) # observation equations
    } 
  }
}

ch[ch == 2] <- 0 # Bernoulli likelihood
dat <- list(first = first, 
            ch = ch, 
            n.sec = n.secondary, 
            n.years = ncol(ch), 
            n.ind = nrow(ch),
            yes = yes, 
            total = total)

# initial values, parameters to monitor, MCMC settings
# (number of iterations for burnin and post-inference probably need to be increased)
inits <- function(){list(z = z.init)}  
pars <- c('pstar','mean.p','phi','gam')
n.chains <- 1
n.iter <- 1000
n.burnin <- 500

# fitting model to data
res_random <- jags(data = dat, 
                   inits = inits, 
                   parameters.to.save = pars,
                   model.file = model, 
                   n.chains = n.chains,
                   n.iter = n.iter, 
                   n.burnin = n.burnin)

# posterior density distribution of the parameters
jagsfit.mcmc <- as.mcmc(res_random)
densityplot(jagsfit.mcmc)

# display results
summary(jagsfit.mcmc)

# posterior means for population size
Nmcmc <- matrix(NA, nrow(res_random$BUGSoutput$sims.list$pstar), n.primary)
for (i in 1:n.primary){
  Nmcmc[,i] <- caught[i] / res_random$BUGSoutput$sims.list$pstar[,i]
}
apply(Nmcmc,2,mean)
