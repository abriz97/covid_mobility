################
# DEPENDENCIES #
################

library(data.table)
library(ggplot2)

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
data.dir <- file.path(indir, 'data')

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

list.files(data.dir)

path.apple.mobility <- file.path(data.dir, "applemobilitytrends-2022-04-03.csv")
path.WHO.covid <- file.path(data.dir,"WHO-COVID-19-global-data.csv")

dmob <- fread(path.apple.mobility, header=TRUE)
dcovid <- fread(path.WHO.covid)


# reshape dmob
idx <- grepl('^20', names(dmob))
cols0 <- names(dmob)[idx]
cols1 <- names(dmob)[!idx]
dmob <- melt(dmob, id.vars=cols1, measure.vars=cols0)
setnames(dmob, c('variable', 'value'), c('date', 'mobility'))
dmob[, date:=as.Date(date, format='%Y-%m-%d')]
dmob[, mobility:=as.numeric(mobility)]
dmob <- dmob[, lapply(.SD, empty2na) , ]
dmob[, alternative_name := NULL ]

# for some reason need to re-specify type
# Do I now?
dmob[, date:=as.Date(date, format='%Y-%m-%d')]
dmob[, mobility:=as.numeric(mobility)]
dmob[, transportation_type := as.factor(transportation_type)]


# dcovid seems better presented 
setnames(dcovid, 'Date_reported', 'date')
dcountry <- make_dcountries(dcovid)

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

idx <- unique(dmob$Country) %in% dcountry$Country2
unique(dmob$Country)[!idx]

dmob <- merge(dmob, dcountry[, .(Country, Country_code)])
dcovid <- merge(dcovid, dcountry[, .(Country2, Country_code)])

# namibia is not in dmob
dcovid[is.na(Country_code), unique(Country2)]

dcovid
dmob

state='IT'

plot_mobility_cases <- function(state)
{
        # need to change
        stopifnot(state %in% dcountry$Country_code)

        str(dmob)
        tmp_mob <- dmob[Country_code == state]
        tmp_cov <- dcovid[Country_code == state]
        tmp_cov[, date:=as.Date(date, format='%Y-%m-%d')]

        daxis <- rbind(
                tmp_cov[, range(log(New_cases))],
                tmp_mob[, range(mobility, na.rm=TRUE)]
        ) |> as.data.table()
        names(daxis) <- c('min', 'max' )
        daxis[, max:=as.numeric(max)]
        daxis[, min:=as.numeric(min)]
        daxis[, delta := max-0]
        daxis[, delta2 := c(delta[1]/delta[2], delta[2]/delta[1])]
        # force(daxis)

        str(tmp_mob)
        tmp_mob[, mobility]
        daxis$delta2[2]
        tmp_mob
        tmp_mob[, mob_scaled := mobility / daxis$delta2[2]]

        str(tmp_mob)
        str(tmp_cov)
        tmp_mob[, date:=as.Date(date, format='%Y-%m-%d')]

        span <- 1/21
        ggplot(data=tmp_cov, aes(x=date) ) + 
                geom_line(data=tmp_cov, aes(y=log(New_cases)), alpha=.2) +
                geom_line(data=tmp_mob, aes(x=date, y=mob_scaled, color=transportation_type), alpha=.2) +
                scale_y_continuous(sec.axis = sec_axis( ~ . * daxis[, delta2[2]]),
                                   labels = scales::percent) +
                geom_smooth(data=tmp_cov, aes(y=log(New_cases)), color='black', span=span) +
                geom_smooth(data=tmp_mob, aes(x=date, y=mob_scaled, color=transportation_type, fill=transportation_type), span=span) +
                scale_y_continuous(
                                   sec.axis = sec_axis( ~ . * daxis[, delta2[2]],
                                                       name='mobility indices')) +
                theme_bw() +
                labs(x='Date', y='New cases', title=paste0('Mobility indices and new cases: ',state)) +
                theme(legend.position='bottom')
}

p0 <- plot_mobility_cases('IT')
p1 <- plot_mobility_cases('ES')
p2 <- plot_mobility_cases('DE')
p3 <- plot_mobility_cases('FR')

ggpubr::ggarrange(p0, p1, p2, p3, ncol=2, nrow=2)

