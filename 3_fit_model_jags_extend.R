jags_fit <- extend.jags(jags_fit, sample=1000, add.monitor=c("int_jk", "int_k", "int_t", "B_g"))

jags_fit <- extend.jags(jags_fit, sample=1000, add.monitor=c("B_g"))

jags_fit <- extend.jags(jags_fit, sample=400, summarise=FALSE)

# Sample including latent diameters
jags_fit <- extend.jags(jags_fit, sample=200, summarise=FALSE, add.monitor=c("dbh_latent"))
