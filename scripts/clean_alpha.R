# raw data copy-pasted into CSV from
# https://en.wikipedia.org/wiki/SARS-CoV-2_Alpha_variant#Spread_in_Europe

my_path = "~/Desktop/covid_mobility"
raw_data = read.csv(file.path(my_path, "/data/alpha_wikipedia.csv"), row.names = 1)

first_val <- function(x){
  return(as.numeric(strsplit(x, '%')[[1]][1]))
}
data = data.frame(apply(raw_data, c(1, 2), first_val))
data$'Week.10'[2] = 100.0
data['Test'] = 0 
for(col in 2:6){
  uk_ratio = data[1, col] / data[1, 7]
  for(row in 3:5){
    data[row, col] = round(data[row, 7] * uk_ratio, 1)
  }
}
data = t(data[2:5, ])

output = data.frame(matrix(ncol = 3, nrow = 0))
date = as.Date('2020-10-05')
for(week in 1:dim(data)[1]){
  for(country in c('England', 'Northern Ireland', 'Scotland', 'Wales')){
    output = rbind(output, c(country, as.character(date), data[week, country]))
    if(week < dim(data)[1]){  
      for(daystep in 1:6){
        output = rbind(output, c(country, as.character(date+daystep), NA))
      }
    }
  }
  date = date + 7
}
colnames(output) = c("country", "date", "proportion")
output = output[order(output$country), ]
output = output %>% mutate(proportion = na.approx(proportion))
output$proportion = round(output$proportion, 2)
write.csv(output, "~/Desktop/alpha_proportion.csv", row.names=FALSE)
