#### translation to covid cases

library(EpiEstim)
library(ggplot2)
library(incidence)
library(data.table)

path_res = '/Users/annamenacher/Downloads/results_epi/'

data_org = read.csv('https://raw.githubusercontent.com/abriz97/covid_mobility/main/data/UKgovernment_deaths_cases.csv')
data_phi = read.csv('https://raw.githubusercontent.com/abriz97/covid_mobility/main/data/ratio_alpha_mortality.csv')
data_variant = read.csv('https://raw.githubusercontent.com/abriz97/covid_mobility/main/data/alpha_breakdown.csv')

country = 'England'

data_org = data_org[which(data_org$country == country, arr.ind = T),]
data_phi = data_phi[which(data_phi$Country == country, arr.ind = T),]
data_variant = data_variant[which(data_variant$country == country, arr.ind = T),]

###################
country = 'England'
variant = 'alpha'
##################

ind_phi = which(data_phi$date == '2020-10-07', arr.ind = T)
ind_phi2 = which(data_phi$date == '2021-04-01', arr.ind = T)

ind_org = which(data_org$date == '2020-10-07', arr.ind = T)
ind_org2 = which(data_org$date == '2021-04-01', arr.ind = T)


# dates from 2020-10-07 to 2021-04-01
countries = data_phi$Country[ind_phi:ind_phi2]
dates = data_phi[ind_phi:ind_phi2,2]
deaths = round(data_phi[ind_phi:ind_phi2,7]*data_org[ind_org:ind_org2,3],0)
cases = round(data_phi[ind_phi:ind_phi2,3]*data_phi[ind_phi:ind_phi2,4],0)

data_alpha = data.frame('date' = dates, 'cases' = cases, 'deaths' = deaths)

plot(as.incidence(data_alpha$cases, dates = as.Date(data_alpha$date)))


# use infection interval with mean 6.48 and sd 3.83
data = list('incidence' = data.frame(dates = as.Date(data_alpha$date[27:nrow(data_alpha)]), I = data_alpha$deaths[27:nrow(data_alpha)]))

res_parametric_si <- estimate_R(data$incidence, 
                                method="parametric_si",
                                config = make_config(list(
                                  mean_si = 6.48, 
                                  std_si = 3.83))
)

head(res_parametric_si$R)
pdf(paste0(path_res, 'plots/', variant, '_', country, '.pdf'), width=10, height = 10)
plot(res_parametric_si, legend = FALSE)
dev.off()

col1 = countries
col2 = dates
col3 = rep(NA, length(col1))
col3[(length(col1)-nrow(res_parametric_si$R) + 1):length(col1)] = res_parametric_si$R$`Mean(R)`
col4 = rep(variant, length(col1))

results_england_alpha = data.frame('Country' = col1, 'date' = col2, 'R' = col3, 'variant' = col4)

save(results_england_alpha, res_parametric_si, data_alpha, file = paste0(path_res, "RData/alpha_england.RData"))
write.csv(results_england_alpha, paste0(path_res, 'csv/alpha_england.csv'))

###################
country = 'England'
variant = 'wildtype'
##################

ind_phi = which(data_phi$date == '2020-03-19', arr.ind = T)
ind_phi2 = which(data_phi$date == '2021-04-01', arr.ind = T)
ind_phi3 = which(data_phi$date == '2020-10-08', arr.ind = T)


ind_org = which(data_org$date == '2020-03-19', arr.ind = T)
ind_org2 = which(data_org$date == '2021-04-01', arr.ind = T)
ind_org3 = which(data_org$date == '2020-10-07', arr.ind = T)
ind_org4 = which(data_org$date == '2020-10-08', arr.ind = T)

# dates from 2020-03-19 to 2021-04-01
countries = data_phi$Country[ind_phi:ind_phi2]
dates = data_phi[ind_phi:ind_phi2,2]
# dates from 2020-03-08 to 2020-10-07 AND 2020-10-08 to 2021-04-01
deaths = c(data_org[ind_org:ind_org3,3], round(data_org[ind_org4:ind_org2,3]*(1-data_phi[ind_phi3:ind_phi2,7]),0))
cases = c(data_org[ind_org:ind_org3,5], round(data_org[ind_org4:ind_org2,5]*(1-data_phi[ind_phi3:ind_phi2,3]),0))

data_wildtype = data.frame('date' = dates, 'cases' = cases, 'deaths' = deaths)

plot(as.incidence(data_wildtype$cases, dates = as.Date(data_wildtype$date)))


# use infection interval with mean 6.48 and sd 3.83
data = list('incidence' = data.frame(dates = as.Date(data_wildtype$date[1:nrow(data_wildtype)]), I = data_wildtype$deaths[1:nrow(data_wildtype)]))

res_parametric_si <- estimate_R(data$incidence, 
                                method="parametric_si",
                                config = make_config(list(
                                  mean_si = 6.48, 
                                  std_si = 3.83))
)

head(res_parametric_si$R)
pdf(paste0(path_res, 'plots/', variant, '_', country, '.pdf'), width=10, height = 10)
plot(res_parametric_si, legend = FALSE)
dev.off()

col1 = countries
col2 = dates
col3 = rep(NA, length(col1))
col3[(length(col1)-nrow(res_parametric_si$R) + 1):length(col1)] = res_parametric_si$R$`Mean(R)`
col4 = rep(variant, length(col1))

results_england_wildtype = data.frame('Country' = col1, 'date' = col2, 'R' = col3, 'variant' = col4)

save(results_england_wildtype, res_parametric_si, data_wildtype, file = paste0(path_res, "RData/wildtype_england.RData"))
write.csv(results_england_wildtype, paste0(path_res, 'csv/wildtype_england.csv'))

###################
country = 'England'
variant = 'joint'
##################


ind_phi = which(data_phi$date == '2020-03-19', arr.ind = T)
ind_phi2 = which(data_phi$date == '2021-04-01', arr.ind = T)


ind_org = which(data_org$date == '2020-03-19', arr.ind = T)
ind_org2 = which(data_org$date == '2021-04-01', arr.ind = T)

# dates from 2020-03-19 to 2021-04-01
countries = data_phi$Country[ind_phi:ind_phi2]
dates = data_phi[ind_phi:ind_phi2,2]
deaths = data_org[ind_org:ind_org2,3]
cases = data_org[ind_org:ind_org2,5]

data_joint = data.frame('date' = dates, 'cases' = cases, 'deaths' = deaths)

plot(as.incidence(data_joint$cases, dates = as.Date(data_joint$date)))


# use infection interval with mean 6.48 and sd 3.83
data = list('incidence' = data.frame(dates = as.Date(data_joint$date[1:nrow(data_joint)]), I = data_joint$deaths[1:nrow(data_joint)]))

res_parametric_si <- estimate_R(data$incidence, 
                                method="parametric_si",
                                config = make_config(list(
                                  mean_si = 6.48, 
                                  std_si = 3.83))
)

head(res_parametric_si$R)
pdf(paste0(path_res, 'plots/', variant, '_', country, '.pdf'), width=10, height = 10)
plot(res_parametric_si, legend = FALSE)
dev.off()

col1 = countries
col2 = dates
col3 = rep(NA, length(col1))
col3[(length(col1)-nrow(res_parametric_si$R) + 1):length(col1)] = res_parametric_si$R$`Mean(R)`
col4 = rep(variant, length(col1))

results_england_joint = data.frame('Country' = col1, 'date' = col2, 'R' = col3, 'variant' = col4)

save(results_england_joint, res_parametric_si, data_joint, file = paste0(path_res, "RData/joint_england.RData"))
write.csv(results_england_joint, paste0(path_res, 'csv/joint_england.csv'))








#############
### Wales ###
#############

data_org = read.csv('https://raw.githubusercontent.com/abriz97/covid_mobility/main/data/UKgovernment_deaths_cases.csv')
data_phi = read.csv('https://raw.githubusercontent.com/abriz97/covid_mobility/main/data/ratio_alpha_mortality.csv')
data_variant = read.csv('https://raw.githubusercontent.com/abriz97/covid_mobility/main/data/alpha_breakdown.csv')

country = 'Wales'

data_org = data_org[which(data_org$country == country, arr.ind = T),]
data_phi = data_phi[which(data_phi$Country == country, arr.ind = T),]
data_variant = data_variant[which(data_variant$country == country, arr.ind = T),]

###################
country = 'Wales'
variant = 'alpha'
##################

ind_phi = which(data_phi$date == '2020-10-07', arr.ind = T)
ind_phi2 = which(data_phi$date == '2021-04-01', arr.ind = T)

ind_org = which(data_org$date == '2020-10-07', arr.ind = T)
ind_org2 = which(data_org$date == '2021-04-01', arr.ind = T)

# dates from 2020-10-07 to 2021-04-01
countries = data_phi$Country[ind_phi:ind_phi2]
dates = data_phi[ind_phi:ind_phi2,2]
deaths = round(data_phi[ind_phi:ind_phi2,7]*data_org[ind_org:ind_org2,3],0)
cases = round(data_phi[ind_phi:ind_phi2,3]*data_phi[ind_phi:ind_phi2,4],0)

data_alpha = data.frame('date' = dates, 'cases' = cases, 'deaths' = deaths)

plot(as.incidence(data_alpha$cases, dates = as.Date(data_alpha$date)))


# use infection interval with mean 6.48 and sd 3.83
data = list('incidence' = data.frame(dates = as.Date(data_alpha$date[36:nrow(data_alpha)]), I = data_alpha$deaths[36:nrow(data_alpha)]))

res_parametric_si <- estimate_R(data$incidence, 
                                method="parametric_si",
                                config = make_config(list(
                                  mean_si = 6.48, 
                                  std_si = 3.83))
)

head(res_parametric_si$R)
pdf(paste0(path_res, 'plots/', variant, '_', country, '.pdf'), width=10, height = 10)
plot(res_parametric_si, legend = FALSE)
dev.off()

col1 = countries
col2 = dates
col3 = rep(NA, length(col1))
col3[(length(col1)-nrow(res_parametric_si$R) + 1):length(col1)] = res_parametric_si$R$`Mean(R)`
col4 = rep(variant, length(col1))

results_wales_alpha = data.frame('Country' = col1, 'date' = col2, 'R' = col3, 'variant' = col4)

save(results_wales_alpha, res_parametric_si, data_alpha, file = paste0(path_res, "RData/alpha_wales.RData"))
write.csv(results_wales_alpha, paste0(path_res, 'csv/alpha_wales.csv'))

###################
country = 'Wales'
variant = 'wildtype'
##################

ind_phi = which(data_phi$date == '2020-03-19', arr.ind = T)
ind_phi2 = which(data_phi$date == '2021-04-01', arr.ind = T)
ind_phi3 = which(data_phi$date == '2020-10-08', arr.ind = T)


ind_org = which(data_org$date == '2020-03-19', arr.ind = T)
ind_org2 = which(data_org$date == '2021-04-01', arr.ind = T)
ind_org3 = which(data_org$date == '2020-10-07', arr.ind = T)
ind_org4 = which(data_org$date == '2020-10-08', arr.ind = T)

# dates from 2020-03-19 to 2021-04-01
countries = data_phi$Country[ind_phi:ind_phi2]
dates = data_phi[ind_phi:ind_phi2,2]
# dates from 2020-03-08 to 2020-10-07 AND 2020-10-08 to 2021-04-01
deaths = c(data_org[ind_org:ind_org3,3], round(data_org[ind_org4:ind_org2,3]*(1-data_phi[ind_phi3:ind_phi2,7]),0))
cases = c(data_org[ind_org:ind_org3,5], round(data_org[ind_org4:ind_org2,5]*(1-data_phi[ind_phi3:ind_phi2,3]),0))

data_wildtype = data.frame('date' = dates, 'cases' = cases, 'deaths' = deaths)

plot(as.incidence(data_wildtype$cases, dates = as.Date(data_wildtype$date)))


# use infection interval with mean 6.48 and sd 3.83
data = list('incidence' = data.frame(dates = as.Date(data_wildtype$date[1:nrow(data_wildtype)]), I = data_wildtype$deaths[1:nrow(data_wildtype)]))

res_parametric_si <- estimate_R(data$incidence, 
                                method="parametric_si",
                                config = make_config(list(
                                  mean_si = 6.48, 
                                  std_si = 3.83))
)

head(res_parametric_si$R)
pdf(paste0(path_res, 'plots/', variant, '_', country, '.pdf'), width=10, height = 10)
plot(res_parametric_si, legend = FALSE)
dev.off()

col1 = countries
col2 = dates
col3 = rep(NA, length(col1))
col3[(length(col1)-nrow(res_parametric_si$R) + 1):length(col1)] = res_parametric_si$R$`Mean(R)`
col4 = rep(variant, length(col1))

results_wales_wildtype = data.frame('Country' = col1, 'date' = col2, 'R' = col3, 'variant' = col4)

save(results_wales_wildtype, res_parametric_si, data_wildtype, file = paste0(path_res, "RData/wildtype_wales.RData"))
write.csv(results_wales_wildtype, paste0(path_res, 'csv/wildtype_wales.csv'))

###################
country = 'Wales'
variant = 'joint'
##################

ind_phi = which(data_phi$date == '2020-03-19', arr.ind = T)
ind_phi2 = which(data_phi$date == '2021-04-01', arr.ind = T)


ind_org = which(data_org$date == '2020-03-19', arr.ind = T)
ind_org2 = which(data_org$date == '2021-04-01', arr.ind = T)

# dates from 2020-03-19 to 2021-04-01
countries = data_phi$Country[ind_phi:ind_phi2]
dates = data_phi[ind_phi:ind_phi2,2]
deaths = data_org[ind_org:ind_org2,3]
cases = data_org[ind_org:ind_org2,5]

data_joint = data.frame('date' = dates, 'cases' = cases, 'deaths' = deaths)

plot(as.incidence(data_joint$cases, dates = as.Date(data_joint$date)))


# use infection interval with mean 6.48 and sd 3.83
data = list('incidence' = data.frame(dates = as.Date(data_joint$date[1:nrow(data_joint)]), I = data_joint$deaths[1:nrow(data_joint)]))

res_parametric_si <- estimate_R(data$incidence, 
                                method="parametric_si",
                                config = make_config(list(
                                  mean_si = 6.48, 
                                  std_si = 3.83))
)

head(res_parametric_si$R)
pdf(paste0(path_res, 'plots/', variant, '_', country, '.pdf'), width=10, height = 10)
plot(res_parametric_si, legend = FALSE)
dev.off()

col1 = countries
col2 = dates
col3 = rep(NA, length(col1))
col3[(length(col1)-nrow(res_parametric_si$R) + 1):length(col1)] = res_parametric_si$R$`Mean(R)`
col4 = rep(variant, length(col1))

results_wales_joint = data.frame('Country' = col1, 'date' = col2, 'R' = col3, 'variant' = col4)

save(results_wales_joint, res_parametric_si, data_joint, file = paste0(path_res, "RData/joint_wales.RData"))
write.csv(results_wales_joint, paste0(path_res, 'csv/joint_wales.csv'))








################
### Scotland ###
################

data_org = read.csv('https://raw.githubusercontent.com/abriz97/covid_mobility/main/data/UKgovernment_deaths_cases.csv')
data_phi = read.csv('https://raw.githubusercontent.com/abriz97/covid_mobility/main/data/ratio_alpha_mortality.csv')
data_variant = read.csv('https://raw.githubusercontent.com/abriz97/covid_mobility/main/data/alpha_breakdown.csv')

country = 'Scotland'

data_org = data_org[which(data_org$country == country, arr.ind = T),]
data_phi = data_phi[which(data_phi$Country == country, arr.ind = T),]
data_variant = data_variant[which(data_variant$country == country, arr.ind = T),]

###################
country = 'Scotland'
variant = 'alpha'
##################

ind_phi = which(data_phi$date == '2020-10-07', arr.ind = T)
ind_phi2 = which(data_phi$date == '2021-04-01', arr.ind = T)

ind_org = which(data_org$date == '2020-10-07', arr.ind = T)
ind_org2 = which(data_org$date == '2021-04-01', arr.ind = T)

# dates from 2020-10-07 to 2021-04-01
countries = data_phi$Country[ind_phi:ind_phi2]
dates = data_phi[ind_phi:ind_phi2,2]
deaths = round(data_phi[ind_phi:ind_phi2,7]*data_org[ind_org:ind_org2,3],0)
cases = round(data_phi[ind_phi:ind_phi2,3]*data_phi[ind_phi:ind_phi2,4],0)

data_alpha = data.frame('date' = dates, 'cases' = cases, 'deaths' = deaths)

plot(as.incidence(data_alpha$cases, dates = as.Date(data_alpha$date)))


# use infection interval with mean 6.48 and sd 3.83
data = list('incidence' = data.frame(dates = as.Date(data_alpha$date[45:nrow(data_alpha)]), I = data_alpha$deaths[45:nrow(data_alpha)]))

res_parametric_si <- estimate_R(data$incidence, 
                                method="parametric_si",
                                config = make_config(list(
                                  mean_si = 6.48, 
                                  std_si = 3.83))
)

head(res_parametric_si$R)
pdf(paste0(path_res, 'plots/', variant, '_', country, '.pdf'), width=10, height = 10)
plot(res_parametric_si, legend = FALSE)
dev.off()

col1 = countries
col2 = dates
col3 = rep(NA, length(col1))
col3[(length(col1)-nrow(res_parametric_si$R) + 1):length(col1)] = res_parametric_si$R$`Mean(R)`
col4 = rep(variant, length(col1))

results_scotland_alpha = data.frame('Country' = col1, 'date' = col2, 'R' = col3, 'variant' = col4)

save(results_scotland_alpha, res_parametric_si, data_alpha, file = paste0(path_res, "RData/alpha_scotland.RData"))
write.csv(results_scotland_alpha, paste0(path_res, 'csv/alpha_scotland.csv'))

###################
country = 'Scotland'
variant = 'wildtype'
##################

ind_phi = which(data_phi$date == '2020-03-19', arr.ind = T)
ind_phi2 = which(data_phi$date == '2021-04-01', arr.ind = T)
ind_phi3 = which(data_phi$date == '2020-10-08', arr.ind = T)


ind_org = which(data_org$date == '2020-03-19', arr.ind = T)
ind_org2 = which(data_org$date == '2021-04-01', arr.ind = T)
ind_org3 = which(data_org$date == '2020-10-07', arr.ind = T)
ind_org4 = which(data_org$date == '2020-10-08', arr.ind = T)

# dates from 2020-03-19 to 2021-04-01
countries = data_phi$Country[ind_phi:ind_phi2]
dates = data_phi[ind_phi:ind_phi2,2]
# dates from 2020-03-19 to 2020-10-07 AND 2020-10-08 to 2021-04-01
deaths = c(data_org[ind_org:ind_org3,3], round(data_org[ind_org4:ind_org2,3]*(1-data_phi[ind_phi3:ind_phi2,7]),0))
cases = c(data_org[ind_org:ind_org3,5], round(data_org[ind_org4:ind_org2,5]*(1-data_phi[ind_phi3:ind_phi2,3]),0))

data_wildtype = data.frame('date' = dates, 'cases' = cases, 'deaths' = deaths)

plot(as.incidence(data_wildtype$cases, dates = as.Date(data_wildtype$date)))


# use infection interval with mean 6.48 and sd 3.83
data = list('incidence' = data.frame(dates = as.Date(data_wildtype$date[1:nrow(data_wildtype)]), I = data_wildtype$deaths[1:nrow(data_wildtype)]))

res_parametric_si <- estimate_R(data$incidence, 
                                method="parametric_si",
                                config = make_config(list(
                                  mean_si = 6.48, 
                                  std_si = 3.83))
)

head(res_parametric_si$R)
pdf(paste0(path_res, 'plots/', variant, '_', country, '.pdf'), width=10, height = 10)
plot(res_parametric_si, legend = FALSE)
dev.off()

col1 = countries
col2 = dates
col3 = rep(NA, length(col1))
col3[(length(col1)-nrow(res_parametric_si$R) + 1):length(col1)] = res_parametric_si$R$`Mean(R)`
col4 = rep(variant, length(col1))

results_scotland_wildtype = data.frame('Country' = col1, 'date' = col2, 'R' = col3, 'variant' = col4)

save(results_scotland_wildtype, res_parametric_si, data_wildtype, file = paste0(path_res, "RData/wildtype_scotland.RData"))
write.csv(results_scotland_wildtype, paste0(path_res, 'csv/wildtype_scotland.csv'))

###################
country = 'Scotland'
variant = 'joint'
##################

ind_phi = which(data_phi$date == '2020-03-19', arr.ind = T)
ind_phi2 = which(data_phi$date == '2021-04-01', arr.ind = T)


ind_org = which(data_org$date == '2020-03-19', arr.ind = T)
ind_org2 = which(data_org$date == '2021-04-01', arr.ind = T)

# dates from 2020-03-19 to 2021-04-01
countries = data_phi$Country[ind_phi:ind_phi2]
dates = data_phi[ind_phi:ind_phi2,2]
deaths = data_org[ind_org:ind_org2,3]
cases = data_org[ind_org:ind_org2,5]

data_joint = data.frame('date' = dates, 'cases' = cases, 'deaths' = deaths)

plot(as.incidence(data_joint$cases, dates = as.Date(data_joint$date)))


# use infection interval with mean 6.48 and sd 3.83
data = list('incidence' = data.frame(dates = as.Date(data_joint$date[1:nrow(data_joint)]), I = data_joint$deaths[1:nrow(data_joint)]))

res_parametric_si <- estimate_R(data$incidence, 
                                method="parametric_si",
                                config = make_config(list(
                                  mean_si = 6.48, 
                                  std_si = 3.83))
)

head(res_parametric_si$R)
pdf(paste0(path_res, 'plots/', variant, '_', country, '.pdf'), width=10, height = 10)
plot(res_parametric_si, legend = FALSE)
dev.off()

col1 = countries
col2 = dates
col3 = rep(NA, length(col1))
col3[(length(col1)-nrow(res_parametric_si$R) + 1):length(col1)] = res_parametric_si$R$`Mean(R)`
col4 = rep(variant, length(col1))

results_scotland_joint = data.frame('Country' = col1, 'date' = col2, 'R' = col3, 'variant' = col4)

save(results_scotland_joint, res_parametric_si, data_joint, file = paste0(path_res, "RData/joint_scotland.RData"))
write.csv(results_scotland_joint, paste0(path_res, 'csv/joint_scotland.csv'))



