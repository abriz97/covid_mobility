#AIMS:
# prepare the datasets to run stan on
# need to load data on:
# mobility, UK government data, variant frequency

# Run as 

################
# DEPENDENCIES #
################

library(rstan)
library(data.table)
library(ggplot2)
library(bayesplot)

################
#  ARGUMENTS   #
################

option_list <- list(
  optparse::make_option(
    c("c", "--country"),
    type = "character",
    default = 'England',
    help = "Country for which to run analysis[must be one of England, Scotland, NothernIreland, Wales]",
    dest = 'country'
  ),
  optparse::make_option(
    c("--ratio_alpha_mortality"),
    type = "numeric",
    default = 1.6,
    help = "Ratio between Fatality ratio of Alpha and Fatality Ratio of Wildtype[default=1.6]",
    dest = "ratio_alpha_mortality"
  ),
  optparse::make_option(
    "--max_date",
    type = 'character',
    default = '2021-04-01',
    metavar = '"YYYY-MM-DD"',
    help = 'Last date for analysis',
    dest = 'max_date'
  ),
  optparse::make_option(
    c("-f", "--figures"),
    action = "store_true",
    default = TRUE,
    help = "Save Figures",
    dest = "make.figures"
  )
)

args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
args$max_date <- as.Date(args$max_date, format='%Y-%m-%d')

################
#    PATHS     #
################

indir <- getwd()
indir <- gsub('covid_mobility/.*$', 'covid_mobility', indir)

data.dir <- file.path(indir, 'data')
results.dir <- file.path(indir, 'results')
# path.who.deaths = file.path(data.dir, 'WHO_coviddeaths_bycountry_processed.csv')
path.mobility.data.uk <- file.path(data.dir, 'apple_mobility_bycountry_processed_UK.csv')
path.gov.data <- file.path(data.dir, 'UKgovernment_deaths_cases.csv')
# make.figures <- 1
path.alpha.proportions <- file.path(data.dir,'alpha_proportion.csv')



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

convolve.past <- function(vec1, vec2)
{
        L1 <- length(vec1)
        L2 <- length(vec2)
        out <- rep(NA, length(vec1))

        vec1_pad <- c(rep(0, L2 - 1 ), vec1)
        vec2_rev <- rev(vec2)
        
        for(i in seq_along(vec1))
        {
                vec2_pad <- c(rep(0, i-1), vec2_rev, rep(0, L1 - i ) )
                out[i] <- as.numeric(vec1_pad %*% vec2_pad)
        }

        return(c(NA, out[-length(out)]))
}

make_dataset <- function(country=args$country, max_date=args$max_date, ratio_alpha_mortality = args$ratio_alpha_mortality)
{

        # ratio_alpha_mortality = 1.6
        # country= 'England'


        .make_discrete_gamma_pmf <- function(mean, std)
        {
                # r uses scale rate parametrisation, where
                # mean = shape / rate
                # std = sqrt(shape) / rate
                s = mean^2/std^2
                r = mean/std^2

                if(0) # test
                {
                        tmp <- rgamma(1000, shape=s, rate=r)
                        mean(tmp); sqrt(var(tmp))
                }
                probs <- sapply(0:100, pgamma, shape=s, rate=r)
                probs <- c(probs[-1], NA) - probs
                idx <- which(1 - cumsum(probs) <= 10e-4)[1]
                probs <- probs[1:idx]
                probs/sum(probs)
        }

        .ma <- function(x, n = 5) as.numeric(stats::filter(x, rep(1 / n, n), sides = 2))

        # 
        cols <- c('Country', 'date', 'deaths', 'cumdeaths')
        dcov_tmp <- dcovid[ Country == country , ..cols ]

        cols <- c('Country', 'date', 'mobility', 'transportation_type')
        dmob_tmp <- dmobility[ Country == country, ..cols ]

        cat('Interpolation to fix NAs in mobility... \n')
        dmob_tmp[, mobility := zoo::na.approx(mobility) , by='transportation_type']

        cat('Taking averages of the 3 transportation type...\n')
        dmob_tmp <- dmob_tmp[ , list(mobility = mean(mobility)), by=c('Country', 'date') ]
        dmob_tmp[, mobility := .ma(mobility)]

        # merge
        ddata <- merge(dmob_tmp, dcov_tmp, by=c('Country', 'date'))
        ddata <- ddata[!is.na(mobility)]

        cat('Selecting study period...\n')
        cat('Start epidemic from day of 30th cum death as in Flaxman et al...: ')
        ddata <- ddata[cumdeaths >= 30]
        cat(ddata[, as.character(min(date))], '\n')
        cat('Truncating epidemic on:',  as.character(max_date), '\n')
        ddata <- ddata[date <= max_date]

        if(0)
        {
                ddata <- .truncate_before_n_consecutive_reports(ddata, 7)
        }

        # Select serial interval for wildtype and alpha
        serial_interval <- list(wildtype=list(),alpha=list())
        serial_interval$wildtype <- list(mean = 6.48, std=3.83)
        serial_interval$alpha <-  list(mean = 6.48, std=3.83)

        case2death_interval <- list(wildtype=list(),alpha=list())
        case2death_interval$wildtype <- list(mean= 18.8, std=8.46)
        case2death_interval$alpha <- list(mean= 18.8, std=8.46)

        stan_data <- list(
                N = ddata[, .N],
                deaths = ddata[, deaths] ,
                mob = ddata[, mobility/100]
        )

        # add distribution of interval from case to death
        stan_data$h_wildtype <- .make_discrete_gamma_pmf(mean=case2death_interval$wildtype$mean,
                                               std=case2death_interval$wildtype$std)
        stan_data$L_h_wildtype <- length(stan_data$h_wildtype)

        stan_data$h_alpha <- .make_discrete_gamma_pmf(mean=case2death_interval$alpha$mean,
                                               std=case2death_interval$alpha$std)
        stan_data$L_h_alpha <- length(stan_data$h_alpha)

        # add serial interval distribution
        stan_data$w_wildtype <- .make_discrete_gamma_pmf(mean=serial_interval$wildtype$mean,
                                               std=case2death_interval$wildtype$std)
        stan_data$L_w_wildtype <- length(stan_data$w_wildtype)

        stan_data$w_alpha <- .make_discrete_gamma_pmf(mean=case2death_interval$alpha$mean,
                                               std=case2death_interval$alpha$std)
        stan_data$L_w_alpha <- length(stan_data$w_alpha)

        # find estimated proportion of deaths per variant 
        dalpha_tmp <- dalpha[Country==country,,]
        dalpha_tmp <- merge(dalpha_tmp, 
                            dcovid[, .(Country, date,cases)],
                            by=c('Country', 'date'))

        dalpha_tmp[, cases := .ma(cases, n=7)]
        dalpha_tmp[!is.na(cases)]
        ggplot(dalpha_tmp, aes(x=date, y=cases)) + geom_point()

        dalpha_tmp <- dalpha_tmp[!is.na(cases)]
        dalpha_tmp[, alpha_conv := convolve.past(proportion*cases, stan_data$h_alpha )]
        dalpha_tmp[, wildtype_conv := convolve.past((1-proportion)*cases, stan_data$h_wildtype )]

        dalpha_tmp[, tmp := ratio_alpha_mortality * alpha_conv / wildtype_conv ]
        dalpha_tmp[, phi := 1/(1+tmp^(-1))]
        dalpha_tmp[, tmp := NULL]

        ggplot(dalpha_tmp, aes(x=date, y=phi)) + geom_line()

        if(args$make.figures) # Show deaths attribution by variant
        {

                tmp <- merge(ddata, dalpha_tmp, by=c('Country','date'))
                tmp <- tmp[, list(
                                  proportion=proportion,
                                  deaths_alpha = deaths * phi,
                                  deaths_wildtype = deaths * (1 - phi)
                                  ), by=c('Country', 'date')]

                tmp <- melt(tmp, id.vars=c('Country' , 'date', 'proportion'),
                            value.name='deaths',
                            variable.name='variant')
                tmp[, variant:=gsub('deaths_', '', variant)]
                setkey(tmp, date)
                tmp

                p <- ggplot(tmp, aes(x=date)) +
                        geom_bar(position='stack', stat='identity', aes(fill=variant, y=deaths)) + 
                        geom_line(aes(x=date, y=proportion*max(ddata$deaths))) +
                        scale_x_date(breaks='2 weeks' , date_labels = "%d-%b-%y") +
                        scale_y_continuous('COVID-19 deaths reported by UK government',
                                           sec.axis = sec_axis(~. /max(ddata$deaths), 
                                                               name='Variant frequencies among reported cases',
                                                               labels=scales::percent) ) + 
                        theme_bw() + 
                        theme(legend.position='bottom',
                              axis.text.x=element_text(angle=45, vjust=0.5, hjust=0.5)) + 
                        labs(fill='attributed variant', x='',
                             title=paste0('Estimated deaths by variant: ', country))
                p

                filename=file.path(results.dir, 'deaths_attributed_to_variant.png')
                ggsave(filename, p, width=10, height=8)
        }

        # merge with ddata to subset to correct time window
        dalpha_tmp <- merge(ddata[, .(date)], dalpha_tmp, all.x=TRUE)

        stan_data$phi <- dalpha_tmp$phi

        stan_data$dates <- unique(dalpha_tmp$date)

        # Some checks
        cnd <- lapply(stan_data, length)
        stopifnot( uniqueN(c(cnd$deaths, cnd$mob, cnd$phi)) == 1)
        stopifnot( cnd$w_wildtype == stan_data$L_w_wildtype)
        stopifnot( cnd$w_alpha == stan_data$L_w_alpha)
        stopifnot( cnd$h_wildtype == stan_data$L_h_wildtype)
        stopifnot( cnd$h_alpha == stan_data$L_h_alpha)
 
        return(stan_data)
}

################
#     MAIN     #
################

## GET DATA
#__________

# dcovid <- fread(path.who.deaths)
dcovid <- fread(path.gov.data)
dmobility <- fread(path.mobility.data.uk)
setnames(dcovid, 'country', 'Country')
setnames(dmobility, 'country', 'Country')

# while waiting for alpha prevalence, assume it increases linearly from 0 to 1 linearly.
dalpha <- fread(path.alpha.proportions)
setnames(dalpha, 'country', 'Country')
dalpha[, proportion := proportion/100]
dalpha <- merge(dalpha, dcovid[, .(Country, date)], all.y=T)


cat('Assume alpha can t go down after maximum...\n')
dalpha[!is.na(proportion), {
        z <- max(proportion, na.rm=TRUE);
        idx <- which.max(proportion);
        list(proportion = c(proportion[1:idx], rep(z, .N-idx)), date)
        }, by='Country'] -> dalpha
cat('Assume it is 1 by end of study period:', as.character(args$max_date) ,'...\n\n')
dalpha <- merge(dalpha, dcovid[, .(Country, date)], all.y=T, by=c('Country', 'date'))
dalpha[ date >= args$max_date, proportion := 1 ]
dalpha[ is.na(proportion) & date < '2021-01-01', proportion := 0 ]
dalpha[ , proportion := zoo::na.approx(proportion) ]

p <- ggplot(dalpha,  aes(x=date, y=proportion, color=Country)) +
        geom_line() +
        theme_bw() + 
        scale_y_continuous(labels=scales::percent_format(scale=100))+
        labs(x='', y='Proportion',
             title='Estimated proportion of reported cases with alpha variant')

filename <- file.path(results.dir,
                      paste0(args$country ,'_alpha_by_country.png'))
ggsave(filename, p, width=6, height=5)

## FIT STAN
#__________

stan_data <- make_dataset(country=args$country)
options(mc.cores = parallel::detectCores())

stan.model.path <- file.path(indir, 'stan', 'mobility.stan')
model <- stan_model(stan.model.path)
fit <- sampling(model,
                data=stan_data, 
                chains = 4,
                iter=1000
)
print(fit, pars = 'RD_alpha[1]')
print(fit, pars = 'RD_wildtype[1]')
# print(fit, pars = 'mu[1]')

if(0)
{
        # test whether convolutions are working:
        # seems ok except for the fact that some 0 are NA...
        nms_w <- grep('^D_wildtype',names(fit), value=T)
        nms_a <- grep('^D_alpha',names(fit), value=T)

        tmp <- as.array(fit)
        tmp_w <- tmp[,, dimnames(tmp)$parameters %in% nms_w ]
        tmp_a <- tmp[,, dimnames(tmp)$parameters %in% nms_a ]

        tmp_w
        tmp_a

        tmp_w <- unname(tmp_w)[2, 2, ]
        tmp_a <- unname(tmp_a)[2, 2, ]

        stan_data$phi
        tail(tmp_w)
        tail(tmp_a)
}

if(0)
{
        # save fit
        tmp <- format(Sys.Date(), '%y%m%d')
        filename = file.path(results.dir,
                             paste0('chains_', args$country, '_', tmp, '.rds'))
        saveRDS(fit, filename)
}else{
        tmp <- format(Sys.Date(), '%y%m%d')
        filename = file.path(results.dir,
                             paste0('chains_', args$country, '_', tmp, '.rda')
        )
        save.image(filename)
}

# Get parameter names stan
summary(fit)
prms <- grep('^R0|^beta|delta', names(fit), value=TRUE)


mcmc_pairs(fit, regex_pars = '^R0')
mcmc_pairs(fit, regex_pars = '^beta')
mcmc_pairs(fit, regex_pars = '^R0')

# First plots
plot(fit, pars = prms)
plot(fit, pars = prms,
     show_density = TRUE, ci_level = 0.5, fill_color = "purple")
plot(fit, pars = prms[1],
     plotfun = "hist", include = FALSE)

plot(fit, pars = grep('R0', prms, value=TRUE),  
     plotfun = "trace" , inc_warmup = TRUE)
plot(fit, pars = grep('beta', prms, value=TRUE),  
     plotfun = "trace" , inc_warmup = TRUE)

plot(fit, pars = prms,
     plotfun = "rhat") + ggtitle("Example of adding title to plot")
