# divide the deaths by region into wildtype and alpha
# uses total deaths and alpha case prevalence as inputs

library(dplyr)

my_path = "~/Desktop/covid_mobility"
alpha_prop = read.csv(file.path(my_path, "/data/alpha_proportion.csv"))
uk_deaths = read.csv(file.path(my_path, "/data/UKgovernment_deaths_cases.csv"))
alpha_mult = 1.6 # assume that alpha is 1.6 times as deadly as wildtype

alpha_prop['death_prop'] = 1/(1 + (100/(alpha_prop['proportion'] * alpha_mult)) - 1/alpha_mult)
output = left_join(alpha_prop, uk_deaths, by = c('country', 'date'))

output['alpha_deaths'] = round(output['death_prop'] * output['deaths'])
output['wildtype_deaths'] = output['deaths'] - output['alpha_deaths']
output['alpha_cases'] = round(output['proportion']/100 * output['cases'])
output['wildtype_cases'] = output['cases'] - output['alpha_cases']

output = subset(output, select = -c(3:8))
write.csv(output, file.path(my_path, "/data/alpha_breakdown.csv"), row.names=FALSE)
