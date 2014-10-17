library(dplyr)
library(ggplot2)
library(lubridate)
library(MuMIn)
library(lme4)
library(parallel)

load("model_data_long.RData")

sitecode_key <- read.csv('H:/Data/TEAM/Sitecode_Key/sitecode_key.csv')
sitecode_key <- select(sitecode_key, sitecode, sitetype, continent)
model_data_long <- merge(model_data_long, sitecode_key, by.x="site_ID", by.y="sitecode")

model_data_long$growth_rgr <- with(model_data_long,
                                   (growth_rgr - mean(growth_rgr)) / sd(growth_rgr))
model_data_long$diameter_end <- with(model_data_long,
                                       (diameter_end - mean(diameter_end)) / sd(diameter_end))
model_data_long$diameter_start <- with(model_data_long,
                                       (diameter_start - mean(diameter_start)) / sd(diameter_start))
# WD is already standardized
# hist(model_data_long$growth_rgr)
# hist(model_data_long$WD)
# hist(model_data_long$diameter_start)
# hist(model_data_long$diameter_end)

img_height <- 4
img_width <- 3
img_dpi <- 300

model_data_long <- filter(model_data_long, n_days > 200)
#model_data_long <- filter(model_data_long, n_days < 550)
# model_data_long <- filter(model_data_long, growth_rgr < 10)
# model_data_long <- filter(model_data_long, growth_rgr > -10)

max(model_data_long$growth_rgr)
min(model_data_long$growth_rgr)

# Exclude NAK, CSN, and YAN - too little data at these sites
#model_data_long <- filter(model_data_long, !(sitecode %in% c("NAK", "CSN", "YAN")))

# Exclude sites with less than 3 years of data:
#model_data_long <- filter(model_data_long, !(sitecode %in% c("BCI", "COU", "KRP", "PSH", 
#"YAS")))

var_order <-c("(Intercept)",
              "Density(0.4,0.6]",
              "Density(0.6,0.8]",
              "Density(> 0.8)",
              "DBH(20,40]",
              "DBH(> 40)",
              "SPI",
              "SPI:Density(0.4,0.6]",
              "SPI:Density(0.6,0.8]",
              "SPI:Density(> 0.8)",
              "DBH(20,40]:SPI",
              "DBH(> 40):SPI",
              "Asia",
              "Latin America")

# Plot whole dataaset dbh_class versus WD_class versus SPI
# mean_growth <- group_by(model_data_long, site_ID, SamplingPeriodEnd, plot_ID) %>%
#     summarise(mean_growth_rgr=mean(growth_rgr),
#               mean_spi_6=mean(spi_6),
#               mean_spi_12=mean(spi_12),
#               mean_spi_24=mean(spi_24))
# mean_growth <- data.frame(mean_growth)
# ggplot(mean_growth, aes(mean_spi_6, mean_growth_rgr, colour=site_ID)) +
#     geom_point()
# ggplot(mean_growth, aes(mean_spi_12, mean_growth_rgr, colour=site_ID)) +
#     geom_point()
# ggplot(mean_growth, aes(mean_spi_24, mean_growth_rgr, colour=site_ID)) +
#     geom_point()

# Ensure that models are not fit to different datasets when performing model 
# selection
options(na.action = "na.fail")

model_data_long$mcwd_6 <- log(1 + model_data_long$mcwd_6)
model_data_long$mcwd_12 <- log(1 + model_data_long$mcwd_12)
model_data_long$mcwd_24 <- log(1 + model_data_long$mcwd_24)

model_data_long$mcwd_6 <- with(model_data_long, (mcwd_6 - mean(mcwd_6)) / sd(mcwd_6))
model_data_long$mcwd_12 <- with(model_data_long, (mcwd_12 - mean(mcwd_12)) / sd(mcwd_12))
model_data_long$mcwd_24 <- with(model_data_long, (mcwd_24 - mean(mcwd_24)) / sd(mcwd_24))

model_data_long$SPI_12_class <- "normal"
model_data_long$SPI_12_class[model_data_long$spi_12 < -1] <- "drought"
model_data_long$SPI_12_class[model_data_long$spi_12 > 1] <- "wet"
table(model_data_long$SPI_12_class)

###############################################################################
# SPI 24 with plot-level random effect
m1 <- lmer(growth_rgr ~ (1|tree_ID), data=model_data_long)
m2 <- lmer(growth_rgr ~ (1|genus_ID), data=model_data_long)
m3 <- lmer(growth_rgr ~ (1|plot_ID) + (1|tree_ID), data=model_data_long)
m4 <- lmer(growth_rgr ~ (1|site_ID) + (1|plot_ID) + (1|tree_ID), data=model_data_long)
m5 <- lmer(growth_rgr ~ (1|genus_ID) + (1|site_ID) + (1|plot_ID) + (1|tree_ID), data=model_data_long)

m6 <- lm(growth_rgr ~ continent +
         spi_24 +
         WD + I(WD^2)+ WD*spi_24 +
         diameter_start + I(diameter_start^2), data=model_data_long)

m7 <- lm(growth_rgr ~ continent +
         mcwd_12 +
         WD + I(WD^2)+ WD*mcwd_12 +
         diameter_start + I(diameter_start^2), data=model_data_long)

m7 <- lmer(growth_rgr ~ continent +
           mcwd_12 + I(mcwd_12^2) +
           WD + I(WD^2)+ WD*mcwd_12 +
           diameter_start + I(diameter_start^2) +
           (mcwd_12 + I(mcwd_12)^2 + diameter_start + I(diameter_start^2)|genus_ID) +
           (1|site_ID) + (1|plot_ID) + (1|tree_ID), data=model_data_long)

m8 <- lmer(growth_rgr ~ continent +
           SPI_24_class +
           WD + I(WD^2)+ WD*spi_24 +
           diameter_start + I(diameter_start^2) +
           (1|genus_ID) + (1|site_ID) + (1|plot_ID) + (1|tree_ID), data=model_data_long)

m8 <- lmer(growth_rgr ~ continent +
           mcwd_12 +
           WD + I(WD^2)+ WD*mcwd_12 +
           diameter_start + I(diameter_start^2) +
           (mcwd_12 + diameter_start + I(diameter_start^2)|genus_ID) +
           (1|site_ID) + (1|plot_ID) + (1|tree_ID), data=model_data_long)


m2 <- lmer(growth_rgr ~ spi_24 + WD_class + dbh_class + (1|sitecode) + (1|plot_ID) + (1|SamplingUnitName), data=model_data_long)
m3 <- lmer(growth_rgr ~ spi_24 + dbh_class + (1|sitecode) + (1|plot_ID) + (1|SamplingUnitName), data=model_data_long)
m4 <- lmer(growth_rgr ~ spi_24 + WD_class + (1|sitecode) + (1|plot_ID) + (1|SamplingUnitName), data=model_data_long)
m5 <- lmer(growth_rgr ~ dbh_class + WD_class + (1|sitecode) + (1|plot_ID) + (1|SamplingUnitName), data=model_data_long)
m6 <- lmer(growth_rgr ~ spi_24 + (1|sitecode) + (1|plot_ID) + (1|SamplingUnitName), data=model_data_long)
m7 <- lmer(growth_rgr ~ WD_class + (1|sitecode) + (1|plot_ID) + (1|SamplingUnitName), data=model_data_long)
m8 <- lmer(growth_rgr ~ dbh_class + (1|sitecode) + (1|plot_ID) + (1|SamplingUnitName), data=model_data_long)
m9 <- lmer(growth_rgr ~ continent + (1|sitecode) + (1|plot_ID) + (1|SamplingUnitName), data=model_data_long)
mlist <- list(m1, m2, m3, m4, m5, m6, m7, m8, m9)
msel <- model.sel(mlist)
msel

models_95 <- get.models(msel, cumsum(weight) <= .95)
models_95_avg <- model.avg(models_95, beta=TRUE, fit=TRUE)
summary(models_95_avg)

models_coefs <- data.frame(var_name=models_95_avg$term.names,
                                  est=summary(models_95_avg)$coefmat[, 1])
models_coefs_conf <- data.frame(confint(models_95_avg))
names(models_coefs_conf) <- c("lower", "upper")
models_coefs <- cbind(models_coefs, models_coefs_conf)
row.names(models_coefs) <- NULL
models_coefs <- models_coefs[models_coefs$var_name != "(Intercept)", ]
models_coefs$var_name <- gsub("_class", "", models_coefs$var_name)
models_coefs$var_name <- gsub("continentAS", "Asia", models_coefs$var_name)
models_coefs$var_name <- gsub("continentLA", "Latin America", models_coefs$var_name)
models_coefs$var_name <- gsub("spi", "SPI", models_coefs$var_name)
models_coefs$var_name <- gsub("dbh", "DBH", models_coefs$var_name)
models_coefs$var_name <- gsub("WD", "Density", models_coefs$var_name)
models_coefs$var_name <- gsub("_24", "", models_coefs$var_name)
models_coefs$var_name <- ordered(models_coefs$var_name, levels=rev(var_order))

ggplot(models_coefs, aes(est, var_name)) +
    theme_bw(base_size=10) +
    geom_point() + 
    geom_errorbarh(aes(xmin=lower, xmax=upper), height=0.5) +
    geom_vline(aes(xintercept=0), color="grey") +
    xlim(c(-.1, .1)) +
    xlab("") + ylab("")
ggsave("growth_model_spi_24_coefs_plot_re.png", width=img_width, 
       height=img_height, dpi=img_dpi)

importance(models_95)

###############################################################################
# Dredge

global_m <- lmer(growth_rgr_std ~ spi_24 + I(spi_24^2) + I(spi_24^3) + 
                 continent + WD_class*spi_24 + dbh_class*spi_24 + (1|sitecode) 
                 + (1|plot_ID) + (1|SamplingUnitName), data=model_data_long, REML=FALSE)
# global_m <- lmer(growth_rgr_std ~ spi_24 + I(spi_24^2) + continent + 
#                  WD_class*spi_24 + dbh_class*spi_24 + (1|sitecode) + 
#                  (1|plot_ID) + (1|SamplingUnitName), data=model_data_long, REML=FALSE)
# global_m <- lmer(growth_rgr_std ~ spi_24 + continent + WD_class*spi_24 + 
#                  dbh_class*spi_24 + (1|sitecode) + (1|plot_ID) + 
#                  (1|SamplingUnitName), data=model_data_long, REML=FALSE)

# empty <- lmer(growth_rgr ~ 1 + (1|sitecode) + (1|plot_ID) + 
#               (1|SamplingUnitName), data=model_data_long, REML=FALSE)

#dredged <- dredge(global_m)

cl <- makeCluster(4)
clusterExport(cl, "model_data_long")
clusterEvalQ(cl, library(lme4))
dredged <- pdredge(global_m, cluster=cl, extra=c("R^2", "adjR^2"))
dredged

# dredge_models_delta_lt_6 <- get.models(dredged, delta <= 6)
# dredge_models_65_avg <- model.avg(dredge_models_delta_lt_6, beta=TRUE, fit=TRUE)
#
# dredge_models_delta_lt_9 <- get.models(dredged, delta <= 9)
# dredge_models_95_avg <- model.avg(dredge_models_delta_lt_9, beta=TRUE, fit=TRUE)

#dredge_models_95 <- get.models(dredged, cumsum(weight) <= .95)
#dredge_models_95 <- get.models(dredged, 1:6)
dredge_models_95 <- get.models(dredged, 1:4)
dredge_models_95_avg <- model.avg(dredge_models_95, beta=TRUE)
summary(dredge_models_95_avg)

dredge_models_coefs <- data.frame(var_name=dredge_models_95_avg$term.names,
                                  est=summary(dredge_models_95_avg)$coefmat[, 1])
dredge_models_coefs_conf <- data.frame(confint(dredge_models_95_avg))
names(dredge_models_coefs_conf) <- c("lower", "upper")
dredge_models_coefs <- cbind(dredge_models_coefs, dredge_models_coefs_conf)
row.names(dredge_models_coefs) <- NULL
dredge_models_coefs <- dredge_models_coefs[dredge_models_coefs$var_name != "(Intercept)", ]
dredge_models_coefs$var_name <- gsub("_class", "", dredge_models_coefs$var_name)
dredge_models_coefs$var_name <- gsub("continentAS", "Asia", dredge_models_coefs$var_name)
dredge_models_coefs$var_name <- gsub("continentLA", "Latin America", dredge_models_coefs$var_name)
dredge_models_coefs$var_name <- gsub("spi", "SPI", dredge_models_coefs$var_name)
dredge_models_coefs$var_name <- gsub("dbh", "DBH", dredge_models_coefs$var_name)
dredge_models_coefs$var_name <- gsub("WD", "Density", dredge_models_coefs$var_name)
dredge_models_coefs$var_name <- gsub("_24", "", dredge_models_coefs$var_name)
dredge_models_coefs$var_name <- gsub("_24", "", dredge_models_coefs$var_name)
dredge_models_coefs$var_name <- gsub("\\(40,500\\]", "(> 40)", dredge_models_coefs$var_name)
dredge_models_coefs$var_name <- gsub("\\(0.8,2\\]", "(> 0.8)", dredge_models_coefs$var_name)
dredge_models_coefs$var_name <- ordered(dredge_models_coefs$var_name, 
                                        levels=rev(var_order))

ggplot(dredge_models_coefs, aes(est, var_name)) +
    theme_bw(base_size=10) +
    geom_point() + 
    geom_errorbarh(aes(xmin=lower, xmax=upper), height=0.5) +
    geom_vline(aes(xintercept=0), color="grey") +
    #xlim(c(-.15, .15)) +
    theme(panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())

ggsave("growth_model_spi_24_dredge_coefs.png", width=img_width, 
       height=img_height, dpi=img_dpi)

save.image(file="1_model_growth_v_climate_lmer.RData")
