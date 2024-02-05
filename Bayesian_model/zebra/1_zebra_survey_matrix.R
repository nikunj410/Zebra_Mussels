library(readr)
start_year = 1980
end_year = 2019
dt = 1/3
Time = as.integer((end_year-start_year+1)/dt)

presence_dat = read_csv("Bayesian_model/processing_raw_data/data_files/presence_data/zebra_presence.csv",show_col_types = FALSE)
Adj_hex = read_csv("Bayesian_model/processing_raw_data/data_files/adj/Adj_hex.csv",show_col_types = FALSE)
Hex_id = as.integer(colnames(Adj_hex))

survey = matrix(0, nrow = Time, ncol = length(Hex_id))
colnames(survey) = colnames(Adj_hex)
for (i in 1:Time){ 
  index = which(presence_dat$date>= (start_year +(i-1)*dt) & presence_dat$date<(start_year+ i*dt))
  if (length(index)!=0) survey[i,match(unique(presence_dat$hex_id[index]),Hex_id)]=1
}
print(sum(survey))
write_csv(as.data.frame(survey),"Bayesian_model/zebra/data_files/zebra_survey_matrix.csv")
