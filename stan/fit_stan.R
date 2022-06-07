# AIMS:


################
# DEPENDENCIES #
################

library(data.table)
library(rstan)
library(ggplot2)

################
#    PATHS     #
################

usr <- Sys.info()[['user']]
if(usr == 'andrea')
{
        indir <- '~/git/covid_mobility/'

}else{

}

data.dir <- file.path(indir, 'data')
path.who.deaths = file.path(data.dir, 'WHO_coviddeaths_bycountry_processed.csv')
path.mobility.data = file.path(data.dir, 'apple_mobility_bycountry_processed.csv')

################
#    HELPERS   #
################

.truncate_before_n_consecutive_reports <- function(DT, n)
{
        cat('Starting epidemic at first', n ,' days of consecutive reported deaths...\n')
        DT[, DUMMY := dplyr::if_else(New_deaths > 0, TRUE, FALSE)]
        tmp <- rle(DT$DUMMY)
        idx <- which(tmp$values == TRUE & tmp$lengths >= 7)[1]
        tmp$lengths <- cumsum(tmp$lengths)
        idx <- ifelse(idx > 0, tmp$lengths[ idx-1 ] + 1, 1)
        DT[, DUMMY := 1:.N]

        DT <- DT[DUMMY >= idx]
        DT[, DUMMY := NULL]
        DT
}

make_dataset <- function( country_code, T=100)
{

        .ma <- function(x, n = 5) stats::filter(x, rep(1 / n, n), sides = 2)

        # 
        cols <- c('Country_code', 'date', 'New_deaths')
        dcov_tmp <- dcovid[ Country_code == country_code , ..cols ]

        cols <- c('Country_code', 'date', 'mobility', 'transportation_type')
        dmob_tmp <- dmobility[ Country_code == country_code , ..cols ]

        cat('Interpolation to fix NAs in mobility... \n')
        dmob_tmp[, mobility := zoo::na.approx(mobility) , by='transportation_type']

        cat('Taking averages of the 3 transportation type...\n')
        dmob_tmp[ , mobility := mean(mobility), by='date' ]
        dmob_tmp[, mobility := as.numeric(.ma(mobility))]

        ddata <- merge(dmob_tmp, dcov_tmp, by=c('Country_code', 'date'))
        ddata <- ddata[!is.na(mobility)]
        ddata <- .truncate_before_n_consecutive_reports(ddata, 7)
        
        stan_data <- list(
                N = ddata[, .N],
                deaths = ddata[, New_deaths] ,
                mob = ddata[, mobility],
                T = T
        )
        
        cnd <- lapply(stan_data, length)
        stopifnot(cnd$deaths == cnd$mob)

        stan_data
}

################
#     MAIN     #
################

dcovid <- fread(path.who.deaths)
dmobility <- fread(path.mobility.data)


stan.model.path <- file.path(indir, 'stan', 'mobility.stan')
model <- stan_model(stan.model.path)

stan_data <- make_dataset('IT')

stan_data

fit <- sampling(mod, data=stan_data)



