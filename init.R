################
# DEPENDENCIES #
################

library(data.table)
library(ggplot2)
require(ggpubr)

################
#    PATHS     #
################

usr <- Sys.info()[['user']]
if(usr == 'andrea')
{
        indir <- '~/git/covid_mobility'
}else{
        indir <- '~/git/covid_mobility'
}

stopifnot(dir.exists(indir))
outdir <- file.path(indir, 'results')
data.dir <- file.path(indir, 'data')

path.apple.mobility <- file.path(data.dir, "applemobilitytrends-2022-04-03.csv")
path.WHO.covid <- file.path(data.dir,"WHO-COVID-19-global-data.csv")

selected.countries <- c("England", "NorthernIreland", "Scotland", "Wales")

################
#    HELPERS   #
################

empty2na <- function(x)
{
        x <- gsub(' ', '', x)
        x[x==''] <- NA
        x
}

make_dcountries <- function(DT)
{
        cols <- c('Country', 'Country_code')
        tmp <- unique(DT[, ..cols])
        tmp[, lapply(.SD, uniqueN) , .SDcols=cols]
        tmp[, Country2 := gsub(' ', '', Country)]
        tmp[Country_code == 'RU', Country2 := 'Russia' ]
        tmp[Country_code == 'CZ', Country2 := 'CzechRepublic' ]
        tmp[Country_code == 'GB', Country2 := 'UnitedKingdom' ]
        tmp[Country_code == 'US', Country2 := 'UnitedStates' ]
        tmp
}

# could have some subsetting functions here

################
#     MAIN     #
################

dmob <- fread(path.apple.mobility, header=TRUE)
dcovid <- fread(path.WHO.covid)

# reshape dmob
idx <- grepl('^20', names(dmob))
cols0 <- names(dmob)[idx]
cols1 <- names(dmob)[!idx]
dmob <- melt(dmob, id.vars=cols1, measure.vars=cols0)
setnames(dmob, c('variable', 'value'), c('date', 'mobility'))
dmob[, date:=as.Date(date, format='%Y-%m-%d')]
dmob[, transportation_type := as.factor(transportation_type)]
# dmob[, mobility:=as.numeric(mobility)]
dmob[, alternative_name := NULL ]
cols <- sapply(dmob, is.character)
cols <- names(cols)[cols]
dmob <- dmob[, (cols):=lapply(.SD, empty2na) , .SDcols=cols]

# SUBSET TO UK + IRELAND
# Check: are there mobility measures specific to countries within the UK
dmob_sub <- dmob[ region %in% selected.countries ]

# dcovid seems better presented a
setnames(dcovid, 'Date_reported', 'date')
dcountry <- make_dcountries(dcovid)

# however there are some issues, ie negative daily cases
cat('Setting reports with negative numbers of cases to 0...\n')
dcovid[New_deaths < 0, New_deaths := 0]

dcovid[Country %in% selected.countries]
dcovid[, unique(Country)]



dmob[, table(geo_type)]
#           city country/region         county     sub-region 
#         641480         124236        2142056         901320 
dmob[geo_type == 'sub-region', ]
dmob[geo_type == 'country/region']


# To match the 2 datasets, maybe just subset to country
#_____
dmob <- dmob[geo_type == 'country/region']
cols <- c('geo_type', 'sub-region', 'country')
set(dmob, NULL, cols, NULL)
setnames(dmob, 'region', 'Country')

# add country code to both datasets
dmob <- merge(dmob, dcountry[, .(Country, Country_code)])
dcovid <- merge(dcovid, dcountry[, .(Country2, Country_code)])

# study for which countries we have both data
# we have 237 countries in the WHO death dataset while mobility only for 50 countries.
dcovid[, uniqueN(Country)]
dmob[, uniqueN(Country)]
tmp1 <- dmob[, unique(Country_code)]
tmp2 <- dcovid[, unique(Country_code)]
tmp0 <- intersect(tmp1, tmp2)
setdiff(tmp1, tmp0)
setdiff(tmp2, tmp0)
# namibia is not in dmob
dcovid[is.na(Country_code), unique(Country2)]



plot_mobility_deaths <- function(state)
{
        # TODO:
        # Improve x axis
        # secondary axis needes percentage labels.

        # need to change
        stopifnot(state %in% dcountry$Country_code)

        tmp_mob <- dmob[Country_code == state]
        tmp_cov <- dcovid[Country_code == state]
        tmp_cov[, date:=as.Date(date, format='%Y-%m-%d')]

        tmp_cov[, TMP:=log(New_deaths)]
        tmp_cov[is.na(TMP)]

        daxis <- rbind(
                tmp_cov[, range(log(New_deaths + 0.05))],
                tmp_mob[, range(mobility, na.rm=TRUE)]
        ) |> as.data.table()
        names(daxis) <- c('min', 'max' )
        daxis[, max:=as.numeric(max)]
        daxis[, min:=as.numeric(min)]
        daxis[, delta := max-0]
        daxis[, delta2 := c(delta[1]/delta[2], delta[2]/delta[1])]

        tmp_mob[, mob_scaled := mobility / daxis$delta2[2]]
        tmp_mob[, date:=as.Date(date, format='%Y-%m-%d')]

        span <- 1/21
        ggplot(data=tmp_cov, aes(x=date) ) + 
                geom_line(data=tmp_cov, aes(y=log(New_deaths + 0.05)), alpha=.2) +
                geom_line(data=tmp_mob, aes(x=date, y=mob_scaled, color=transportation_type), alpha=.2) +
                scale_y_continuous(sec.axis = sec_axis( ~ . * daxis[, delta2[2]]),
                                   labels = scales::percent) +
                geom_smooth(data=tmp_cov, aes(y=log(New_deaths + 0.05)), color='black', span=span) +
                geom_smooth(data=tmp_mob, aes(x=date, y=mob_scaled, color=transportation_type, fill=transportation_type), span=span) +
                scale_y_continuous(
                                   sec.axis = sec_axis( ~ . * daxis[, delta2[2]],
                                                       name='mobility indices')) +
                theme_bw() +
                labs(x='Date', y='New deaths', title=paste0('Mobility indices and new deaths: ',state)) +
                theme(legend.position='bottom')
}

p0 <- plot_mobility_deaths('IT')
p1 <- plot_mobility_deaths('ES')
p2 <- plot_mobility_deaths('DE')
p3 <- plot_mobility_deaths('FR')

p <- ggarrange(p0, p1, p2, p3,
               ncol=2, nrow=2,
               common.legend=TRUE,
               legend='bottom')

filename = file.path(outdir, 'logdeaths_by_mobility_ITESDEFR.png')
ggsave(p, filename=filename, w=9, h=8)

filename=file.path(data.dir, 'WHO_coviddeaths_bycountry_processed.csv')
fwrite(dcovid, filename)
filename=file.path(data.dir, 'apple_mobility_bycountry_processed.csv')
fwrite(dmob, filename)

dmob_sub[, `:=` (geo_type = NULL, `sub-region` = NULL, country=NULL)]
setnames(dmob_sub, 'region', 'country')
filename=file.path(data.dir, 'apple_mobility_bycountry_processed_UK.csv')
fwrite(dmob_sub, filename)
