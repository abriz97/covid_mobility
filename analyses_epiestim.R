# AIMS:


################
# DEPENDENCIES #
################

library(data.table)
library(ggplot2)
library(EpiEstim)

################
#    PATHS     #
################

usr <- Sys.info()[['user']]
if(usr == 'andrea')
{
        indir <- '~/git/covid_mobility'
}else{

}

stopifnot(dir.exists(indir))
outdir <- file.path(indir, 'results')
data.dir <- file.path(indir, 'data')

# TODO: change data to gov + phi data (proportion of alpha deaths)
path.who.deaths=file.path(data.dir, 'WHO_coviddeaths_bycountry_processed.csv')
path.mobility.data=file.path(data.dir, 'apple_mobility_bycountry_processed.csv')

################
#    HELPERS   #
################

# Make (empty) list to fill in with serial interval specifications for each variant.
variants <- c('wildtype','alpha', 'beta', 'gamma', 'delta', 'omicron')
.mk.empty.config <- function(i) list(mean_si=NA_real_, std_si=NA_real_)
serial_intervals <- lapply(variants, .mk.empty.config)
names(serial_intervals) <- variants

serial_intervals[['alpha']]$mean_si <- 5.5
serial_intervals[['alpha']]$std_si <- 3.3

# From the nouvellet paper
serial_intervals[['wildtype']]$mean_si <- 6.48
serial_intervals[['wildtype']]$std_si <- 3.83

if(0)
{
        # Also used as sensitivity
        serial_intervals[['wildtype']]$mean_si <- 4.8
        serial_intervals[['']]$std_si <- 2.7
}




################
#     MAIN     #
################

dcovid <- fread(path.who.deaths)
dmob <- fread(path.mobility.data)


dcovid


# In order to fit EpiEstim need:
# - incidence (easy to get)
# - depending on the method:
# - * si with config
data(Flu2009)
Flu2009

make_config(serial_intervals$beta)

estimate_Rt_by_country <- function(DT=dcovid, country_code='IT', truncate_n=7)
{
        # country_code = 'ES'
        # Truncate til first deaths
        .truncate_before_n_consecutive_reports <- function(n)
        {
                cat('Starting epidemic at first', n ,' days of consecutive reported deaths\n')
                dcov_country[, DUMMY := dplyr::if_else(New_deaths > 0, TRUE, FALSE)]
                tmp <- rle(dcov_country$DUMMY)
                idx <- which(tmp$values == TRUE & tmp$lengths >= 7)[1]
                tmp$lengths <- cumsum(tmp$lengths)
                idx <- ifelse(idx > 0, tmp$lengths[ idx-1 ] + 1, 1)
                dcov_country[, DUMMY := 1:.N]

                dcov_country <- dcov_country[DUMMY >= idx]
                dcov_country[, DUMMY := NULL]
                dcov_country
        }

        # DT <- copy(dcovid)
        # country_code = 'IT'
        # truncate_n = 7
        stopifnot( country_code %in% DT$Country_code )

        cols <- c('date', 'New_deaths')
        dcov_country <- DT[Country_code == country_code, ..cols]
        head(dcov_country,50)

        nrow(dcov_country)
        dcov_country <- .truncate_before_n_consecutive_reports(truncate_n)
        print(head(dcov_country))

        setnames(dcov_country, c('date', 'New_deaths'), c('dates', 'I'))
        dcov_country <- as.data.frame(dcov_country)
        dcov_country$dates <- as.Date(dcov_country$dates)

        # Want weekly sliding windows:
        config_data = serial_intervals$wildtype
        config_data$t_start <- seq(2,(nrow(dcov_country)-6))
        config_data$t_end <- config_data$t_start + 6

        res <- EpiEstim:::estimate_R(as.data.frame(dcov_country), 
                             method='parametric_si',
                             config=make_config(config_data))
        plot(res)
}

estimate_Rt_by_country(country_code='IT')
estimate_Rt_by_country(country_code='ES')
estimate_Rt_by_country(country_code='FR')
estimate_Rt_by_country(country_code='GE')

