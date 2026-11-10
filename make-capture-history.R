# call in data sets (original excel sheets in SMT 2026 Grid Data)
cap <- read.csv("capture-events.csv")
day <- read.csv("session-days.csv")
ind <- read.csv("capture-individuals.csv")

# here is a silly comment
# here is another

# only need first four columns because excel decided to add extra
day <- day[,1:4]
head(ind)
head(cap)
head(day)

# add site column to match dates when at those sites
cap$site <- ind$site[match(cap$tag, ind$tag)]

# check for problem tags
sort(unique(ind$tag))
sort(unique(cap$tag))

# remove any rows where site column is NA (issue for if individual does not have
# tag number to link to site)
cap <- cap[!is.na(cap$site),]

# added empty columns to "cap" data frame
cap$primary <- NA
cap$secondary <- NA

# adds columns to data set "day": now have columns to denote when each primary session began ("first")
# and ended ("last") as well as the days 1-4 (secondary sessions) within each primary session
# "as.Date" converts dates in csv file to be identified as dates in R
day$day1 <- as.Date(day$first)
day$day2 <- day$day1+1
day$day3 <- day$day2+1
day$day4 <- day$day3+1

# remove entry for tag 1033 because it was recorded as different site (noon) than its actual
# site location (field)
cap <- cap[-which(cap$tag==1033 & cap$date=="2026-07-28"),]

# load reshape2 package to use function "melt"
library(reshape2)

# need to tell function which columns to focus on
day.to.melt <- day[,c(1,2,5:8)]

# melt function creates a "molten" data frame
# converts from wide format to long format 
# (where multiple measurement columns are stacked into key-value pairs)
meltday <- melt(day.to.melt,
                id=c("site", "primary"))

# reorder so all sites are grouped together and go down in date order 
# (ex shows every date of beach then every date of field, etc)
meltday <- meltday[order(meltday$site, meltday$primary, meltday$variable),]

# added a column "secondary" that corresponds to which day within each primary session:
# does this by creating a column with values that correspond to the 4th character of the
# "variable" column
meltday$secondary <- as.numeric(substr(meltday$variable,4,4))

# data checks for "cap"
# look for NAs & make sure every date in "cap"
# matches with date(called value) in "meltday"
summary(cap)
length(which(is.na(cap$date)))
length(which(is.na(cap$site)))
cap$date <- as.Date(cap$date)

# must be TRUE
all(cap$date %in% meltday$value)

# filled in blank "primary" and "secondary" columns which were made in earlier code
# for every row within "cap"
# isites looks at the site present at each row
# idays looks at "meltday" to determine secondary and primary values 
#for each site
# irow creates row for every instance where site in metlday 
# matches date in cap
# last portion of code links everything together
for(i in 1:nrow(cap)){
  isite <- cap$site[i]
  idays <- meltday[meltday$site==isite,]
  irow <- which(idays$value == cap$date[i])
  cap$primary[i] <- idays$primary[irow]
  cap$secondary[i] <- idays$secondary[irow]
}

table(cap$primary, cap$secondary)
ftable(secondary~primary+site, data=cap)

# more data cleaning:

# removed "negative tags"
ind <- ind[ind$tag > 0,]

# remove tags 1055 and 1082 because 
#they are marked as two different species
# in original data sheets
ind <- ind[!ind$tag %in% c(1055, 1082),]

# remove a possible duplicate, based on
# coin flip which entry from that day to keep
nrow(cap)
cap <- cap[-which(cap$date=="2026-07-22" & cap$tag==1065)[1],]

# check row number in "cap"
nrow(cap)

# table of NAs (?)
table(cap$tag, useNA="ifany")

# check is there are any NAs
length(which(is.na(cap$tag)))

# create object with same number of rows as ind
nind <- nrow(ind)

# create empty matrix to late store values in for capture history
cmat <- matrix(0, nrow=nind, ncol=20)
nsecondaries <- 4

# create capture history and put it into cmat matrix
# icap connects tag numbers from cap and ind
# icols is an equation we created to number columns (session)
for(i in 1:nind){
  icap <- cap[cap$tag==ind$tag[i],]
  icols <- nsecondaries*(icap$primary-1)+icap$secondary
  cmat[i,icols] <- 1
}#i

# rename rows and columns to make it more intuitive
rownames(cmat) <- paste("tag", ind$tag, sep="_")
colnames(cmat) <- paste(rep(1:5,each=4), rep(1:4, times=5), sep="_")

# blank out session 1 for hunt because
# didn't start trapping until the week
# of every other grid's session 2

# find hunt tags
htags <- ind$tag[which(ind$site=="hunt")]
hrows <- paste("tag", htags, sep="_")
all(hrows %in% rownames(cmat))

rm.cols <- grep("1_", colnames(cmat))
cmat[hrows,rm.cols] <- NA
cmat[!rownames(cmat) %in% hrows,1:5]

# fix rownames so they will sort correctly
tagnums <- as.numeric(sapply(strsplit(rownames(cmat), "_"), "[", 2))
newrownames <- paste("tag", 
                     formatC(tagnums, width=4, flag="0"),
                     sep="_")
rownames(cmat) <- newrownames
cmat <- cmat[sort(rownames(cmat)),]
cmat[,1:4]

# create a csv of capture history and saves into directory
write.csv(cmat, "capture-history-matrix.csv", row.names=rownames(cmat))
