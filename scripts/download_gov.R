# AIMS:
# Download WHO COVID-19 confirmed cases and deaths from their github repo
# from: https://coronavirus.data.gov.uk/details/cases?areaType=nation&areaName=Wales
# and similar.
# cases by date reported and 
# Daily deaths with COVID-19 on the death certificate by date of death

################
# DEPENDENCIES #
################

library(data.table)

################
#    PATHS     #
################

usr <- Sys.info()[['user']]
if(usr == 'andrea')
{
        indir <- '~/git/covid_mobility/'
}else{
        indir <- '~/git/covid_mobility/'
}

data.dir <- file.path(indir, 'data')
selected.countries <- c("England", "NorthernIreland", "Scotland", "Wales")

################
#     MAIN     #
################

# Merge the downloaded cases in a unique dataset.
tmp <- file.path(data.dir, paste0('gov_cases_', tolower(selected.countries), '.csv'))
stopifnot(all(file.exists(tmp)))
dcases <- lapply(tmp, fread)
dcases <- do.call('rbind', dcases)

# same for deaths
tmp <- file.path(data.dir, paste0('gov_deaths_', tolower(selected.countries), '.csv'))
stopifnot(all(file.exists(tmp)))
ddeaths <- lapply(tmp, fread)
ddeaths <- do.call('rbind', ddeaths)

# merge, rename cols
ddata <- merge(ddeaths, dcases)
ddata[, `:=` (areaType = NULL, areaCode = NULL)]
names(ddata)
setnames(ddata,
        c("areaName","date","newDailyNsoDeathsByDeathDate","cumDailyNsoDeathsByDeathDate","newCasesByPublishDate","cumCasesByPublishDate"), 
        c('country', 'date', 'deaths', 'cumdeaths', 'cases', 'cumcases')
)

# wont need after May 2021
ddata <- ddata[date < '2021-05-01']

# save
filename = file.path(data.dir, 'UKgovernment_deaths_cases.csv')
fwrite(ddata, filename)
