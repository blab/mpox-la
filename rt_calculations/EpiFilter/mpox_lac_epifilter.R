######################################################################
## Apply EpiFilter to mpox in LAC
######################################################################

# Notes and assumptions
# this is based on the COVID in New zealand vignette found on the EpiFilter github
# - estimate reproduction numbers with EpiFilter 

# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()

# Set working directory to source
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
# Folder path for results
folres = paste0("./results/mpox/")

# Main functions to run EpiFilter
files.sources = list.files(path = "./main")
for (i in 1:length(files.sources)) {
  source(paste0(c("./main/", files.sources[i]), collapse = ''))
}

# Load data from WHO for New Zealand
alldata = read.csv('../../mpox_la/multitree_coalescent/data/monkeypox_data.csv')

#alldata <- alldata %>% filter(date > "2023-04-13")

# Incidence and dates
Iday = alldata$cases
dates  = alldata$date
# Time series lengths
nday = length(dates); tday = 1:nday

# Approxumate serial interval distribution from https://wwwnc.cdc.gov/eid/article/30/6/23-1095_article
wdist = dgamma(tday, shape = 2.9, scale = 2.9)

# Total infectiousness
Lday = rep(0, nday) 
for(i in 2:nday){
  # Total infectiousness
  Lday[i] = sum(Iday[seq(i-1, 1, -1)]*wdist[1:(i-1)])    
}

######################################################################
## EpiFilter: provides formally smoothed and exact estimates
# Method based on Bayesian recursive filtering and smoothing
######################################################################

# Setup grid and noise parameters
Rmin = 0.01; Rmax = 10; eta = 0.1

# Uniform prior over grid of size m
m = 200; pR0 = (1/m)*rep(1, m)
# Delimited grid defining space of R
Rgrid = seq(Rmin, Rmax, length.out = m)

# Filtered (causal) estimates as list [Rmed, Rhatci, Rmean, pR, pRup, pstate]
Rfilt = epiFilter(Rgrid, m, eta, pR0, nday, Lday[tday], Iday[tday], 0.025)
# Causal predictions from filtered estimates [pred predci]
Ifilt = recursPredict(Rgrid, Rfilt[[4]], Lday[tday], Rfilt[[3]], 0.025)

# Smoothed estimates as list of [Rmed, Rhatci, Rmean, qR]
Rsmooth = epiSmoother(Rgrid, m, Rfilt[[4]], Rfilt[[5]], nday, Rfilt[[6]], 0.025)
# Smoothed predictions from filtered estimates [pred predci]
Ismooth = recursPredict(Rgrid, Rsmooth[[4]], Lday[tday], Rsmooth[[3]], 0.025)

 # Plot estimates and predictions from filtering
plotEpiFilter(Rfilt[[3]][2:nday], Rfilt[[2]][, 2:nday], Ifilt[[1]], Ifilt[[2]],
              'EpiFilter', Iday[2:nday], folres, eta)

# Plot estimates and predictions from smoothing
plotEpiFilter(Rsmooth[[3]][2:nday], Rsmooth[[2]][, 2:nday], Ismooth[[1]], Ismooth[[2]],
              'EpiSmooth', Iday[2:nday], folres, eta)


rhat <-Rsmooth[[3]][2:nday]
Rhatci <-  Rsmooth[[2]][, 2:nday]

write.csv(rhat, "epifilter_r_estimates.csv", row.names=F)
write.csv(Rhatci, "epifilter_r_ci_estimates.csv", row.names=F)


