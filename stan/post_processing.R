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

# calculate convolutions for mean of NB:
with(stan_data,
     {
             tmp <- convolve(deaths * (1-phi),  w_wildtype, type='open');
             data.table(
                        date = dates,
                        wildtype = head(c(0, tmp), length(deaths))
             )
     }
) -> conv_wildtype

with(stan_data,
     {
             tmp <- convolve(deaths * phi,  w_alpha, type='open');
             data.table(
                        date = dates,
                        alpha = head(c(0, tmp), length(deaths))
             )
     }
) -> conv_alpha
convs_phiDw <- merge(conv_wildtype, conv_alpha)
rm(conv_wildtype, conv_alpha)


# Make custom posterior plots
#____________________________


# Plot posterior RD's 
pa <- rbind(
            make.quantiles(pos$RD_alpha, prm='RD_alpha'), 
            make.quantiles(pos$RD_wildtype, prm='RD_wildtype')
)

g <- ggplot(pa, aes(x=date, color=prm, fill=prm)) + 
        geom_ribbon(aes(ymin=IL, ymax=IU), alpha=.5) + 
        geom_line(aes(y=M)) +
        scale_x_date(breaks='1 month' , date_labels = "%b-%y") +
        theme_bw() +
        theme(legend.position = 'bottom') + 
        labs(x='', y='variant-specific Deaths reproduction number')

# mean of NB against deaths:
setkey(pa, prm, day, date)
cols <- c('CL', 'IL', 'M', 'IU', 'CU')

tmp1 <- pa[ prm == 'RD_alpha', lapply(.SD, `*`, convs_phiDw$alpha ) , .SDcols=cols]
tmp1[, `:=` (day = 1:.N, variant='alpha')]
tmp2 <- pa[ prm == 'RD_wildtype', lapply(.SD, `*`, convs_phiDw$wildtype ) , .SDcols=cols]
tmp2[, `:=` (day = 1:.N, variant='wildtype')]
pas <- rbind(tmp2, tmp1)
pas <- merge(dict_daydate, pas, all.y=T)
rm(tmp1, tmp2)

# get deaths data
with(stan_data,
     rbind(
        data.table(date=dates, deaths=deaths*(1-phi),variant='wildtype'),
        data.table(date=dates, deaths=deaths*phi, variant='alpha')
     )
) -> tmp
tmp

# Joint
p0 <- ggplot(data=pas, aes(x=date, color=variant, fill=variant)) + 
        geom_col(data=tmp, aes(y=deaths, fill=variant)) + 
        geom_line(aes(y=M)) + 
        geom_ribbon(aes(ymin=IL, ymax=IU), alpha=.5) +
        facet_grid(~'Overall') +
        scale_x_date(breaks='1 month' , date_labels = "%b-%y") +
        theme_bw() + 
        theme(legend.position='bottom') +
        labs(x='', y='COVID-19 attributable deaths')

# one var
p1 <- ggplot(data=pas, aes(x=date, color=variant, fill=variant)) + 
        geom_col(data=tmp, aes(y=deaths), color='grey80') + 
        geom_line(aes(y=M)) + 
        geom_ribbon(aes(ymin=IL, ymax=IU), alpha=.5) + 
        facet_grid(~variant) +
        theme_bw() + 
        scale_x_date(breaks='1 month' , date_labels = "%b-%y") +
        theme(legend.position='bottom') + 
        labs(x='', y='COVID-19 attributable deaths')

require(ggpubr)
ggarrange(p0, p1,
          nrow=2,
          legend='bottom',
          common.legend=TRUE)


# Plot posterior R's 
pa <- rbind(
            make.quantiles(pos$R_alpha, prm='R_alpha'), 
            make.quantiles(pos$R_wildtype, prm='R_wildtype')
)

# TODO: maybe add mobility...
ggplot(pa, aes(x=date, color=prm, fill=prm)) + 
        geom_ribbon(aes(ymin=IL, ymax=IU), alpha=.5) + 
        geom_line(aes(y=M)) + 
        scale_x_date(breaks='1 month' , date_labels = "%b-%y") +
        theme_bw() +
        theme(legend.position = 'bottom') + 
        labs(x='', y='variant-specific reproduction number')
