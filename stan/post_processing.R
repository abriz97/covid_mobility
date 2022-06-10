# AIMS:
# Load latest model and do plots

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
    c("-f", "--figures"),
    action = "store_true",
    default = TRUE,
    help = "Save Figures",
    dest = "make.figures"
  )
)

args <- optparse::parse_args(optparse::OptionParser(option_list))

if(is.null(args$country))
{
        args$country <- 'England'
        args$make.figures <- TRUE
}

################
#    PATHS     #
################

indir <- getwd()
indir <- gsub('covid_mobility/.*$', 'covid_mobility', indir)

data.dir <- file.path(indir, 'data')
results.dir <- file.path(indir, 'results')
path.mobility.data.uk <- file.path(data.dir, 'apple_mobility_bycountry_processed_UK.csv')
path.gov.data <- file.path(data.dir, 'UKgovernment_deaths_cases.csv')

################
#    HELPERS   #
################


################
#     MAIN     #
################

# find latest chain for a given country
if(0)
{
        path.chains <- grep(args$country, list.files(results.dir, full.names=TRUE), value=TRUE)
        path.chains <- grep('.rds$', path.chains, value=TRUE)
        tmp <- as.numeric(gsub('[A-z]|\\.|\\/', '', path.chains))
        path.chains <- path.chains[which.max(tmp)]
        fit <- readRDS(path.chains)
}else{
        path.images <- grep(args$country, list.files(results.dir, full.names=TRUE), value=TRUE)
        path.images <- grep('.rda$', path.images, value=TRUE)
        path.images
        tmp <- as.numeric(gsub('[A-z]|\\.|\\/', '', path.images))
        path.images <- path.images[which.max(tmp)]
        load(path.images)
}

pdraws <- as.array(fit)
# pdraws <- fit$draws(format='array', inc_warmup=FALSE)
prms <- dimnames(pdraws)[["parameters"]]
str(pdraws)

select.chains <- seq_along(dimnames(pdraws)[['chains']])

mcmc_pairs(x = pdraws,  regex_pars = 'R0|beta',
           off_diag_args = list(size = 0.75))

mcmc_trace(pdraws, regex_pars = 'R0') 

prms
pos <- list()
{
	# bring samples into long array format 
        gsub('\\[.*?\\]', '', prms) |> unique()

        tmp <- pdraws[, select.chains, grepl('beta_alpha', prms)]
        pos$beta_alpha <- as.vector(tmp)
        tmp <- pdraws[, select.chains, grepl('beta_wildtype', prms)]
        pos$beta_wildtype <- as.vector(tmp)

        tmp <- pdraws[, select.chains, grepl('R0_alpha', prms)]
        pos$R0_alpha <- as.vector(tmp)
        tmp <- pdraws[, select.chains, grepl('R0_wildtype', prms)]
        pos$R0_wildtype <- as.vector(tmp)

        tmp <- pdraws[, select.chains, grepl('delta', prms)]
        pos$delta <- as.vector(tmp)

        tmp <- pdraws[, select.chains, grepl('^R_wildtype\\[', prms)]
        tmp <- unname(apply(tmp[,select.chains,], 3, rbind))
        tmp <- setDT(reshape2::melt(tmp))
        names(tmp) <- c('iter', 'day', 'value')
        pos$R_wildtype <- tmp

        tmp <- pdraws[, select.chains, grepl('^R_alpha\\[', prms)]
        tmp <- unname(apply(tmp[,select.chains,], 3, rbind))
        tmp <- setDT(reshape2::melt(tmp))
        names(tmp) <- c('iter', 'day', 'value')
        pos$R_alpha <- tmp

        tmp <- pdraws[, select.chains, grepl('^RD_wildtype\\[', prms)]
        tmp <- unname(apply(tmp[,select.chains,], 3, rbind))
        tmp <- setDT(reshape2::melt(tmp))
        names(tmp) <- c('iter', 'day', 'value')
        pos$RD_wildtype <- tmp

        tmp <- pdraws[, select.chains, grepl('^RD_alpha\\[', prms)]
        tmp <- unname(apply(tmp[,select.chains,], 3, rbind))
        tmp <- setDT(reshape2::melt(tmp))
        names(tmp) <- c('iter', 'day', 'value')
        pos$RD_alpha <- tmp
}

# translate days to date:
dict_daydate <- data.table(date=stan_data$dates)
dict_daydate[, day := seq(.N)]
dict_daydate




# Make custom posterior plots
#____________________________

make.quantiles <- function(DT, prm = NA_character_)
{
        by_cols <- intersect(names(DT), 'day')
        DT[, {
                z <- quantile(value, probs=c(.025,.25,.5,.75,.975));
                list(CL = z[1], IL = z[2], M = z[3], IU=z[4], CU=z[5])
             }, by=by_cols] -> tmp

        if( !is.na(prm) ) tmp[, prm := prm]

        if( 'day' %in% names(tmp)) tmp <- merge(dict_daydate, tmp,  all.y=TRUE)
        tmp
}

# Plot posterior RD's 
pa <- rbind(
            make.quantiles(pos$RD_alpha, prm='RD_alpha'), 
            make.quantiles(pos$RD_wildtype, prm='RD_wildtype')
)

ggplot(pa, aes(x=date, color=prm, fill=prm)) + 
        geom_ribbon(aes(ymin=IL, ymax=IU), alpha=.5) + 
        geom_line(aes(y=M)) 


# Plot posterior R's 
pa <- rbind(
            make.quantiles(pos$RD_alpha, prm='R_alpha'), 
            make.quantiles(pos$RD_wildtype, prm='R_wildtype')
)

# TODO: maybe add mobility...
ggplot(pa, aes(x=date, color=prm, fill=prm)) + 
        # geom_ribbon(aes(ymin=IL, ymax=IU), alpha=.5) + 
        geom_line(aes(y=M))

pos$RD_alpha[value > 100]
pos$R_alpha[value > 100]

# There definitely are problems, possibly with the convolutions...

?convolve
convolve(c(1,2,3), c(1,2), type='open')
convolve()

tmp <- data.table(
        deaths_wildtype = stan_data$deaths * (1-stan_data$phi),
        deaths_alpha = stan_data$deaths * stan_data$phi
)
tmp[, day := 1:.N]
tmp <- merge(dict_daydate, tmp)


convolve(stan_data$deaths, stan_data$)


pos$beta_alpha
