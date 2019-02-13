co2 = read.csv('CO2.csv') 
co2_vec = co2[,2:13]
vec = as.vector(t(co2_vec))

#data cleaning
idx = which(vec < 0)
vec[idx] = NA   
nrow(co2) 
date = seq(1958,(2008+11/12),by=1/12)
x = date[-idx]
y = vec[-idx]

