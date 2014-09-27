library(dplyr)
library(ggplot2)
library(lubridate)
library(MuMIn)
library(lme4)
library(parallel)

load("growth_filtered_detrended.RData")

img_height <- 4
img_width <- 3
img_dpi <- 300

growth <- filter(growth, ctfsaccept)
growth <- filter(growth, growth_ann_index < 10)
growth <- filter(growth, growth_ann_index > -10)
growth <- filter(growth, n_days > 200)
growth <- filter(growth, n_days < 550)

# Exclude NAK, CSN, and YAN - too little data at these sites
#growth <- filter(growth, !(sitecode %in% c("NAK", "CSN", "YAN")))

growth$SamplingPeriodID <- with(growth, factor(factor(sitecode):factor(SamplingPeriodNumber)))

# Exclude sites with less than 3 years of data:
#growth <- filter(growth, !(sitecode %in% c("BCI", "COU", "KRP", "PSH", 
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

###############################################################################
# Merge climate data
load('../CHIRPS/vg_plot_spis.RData')
spis <- tbl_df(spis)

spis$plot_ID <- as.character(spis$plot_ID)
growth$plot_ID <- gsub("-", "", growth$OnehaPlotNumber)
table(spis$plot_ID)
table(growth$plot_ID[!(growth$plot_ID %in% spis$plot_ID)])
# Note that some PSH plot IDs are missing - check this with Jimmy. The spatial 
# data has the PSH plot IDs ranging from 4-9.

# Calculate closest SPI end point to growth period end point
growth$period_end_month <- round_date(growth$period_end, "month")
spis$period_end_month <- spis$date

spi_6 <- filter(spis, spi_period == 6) %>% select(plot_ID, period_end_month, spi)
names(spi_6)[names(spi_6) == "spi"] <- "spi_6"
spi_12 <- filter(spis, spi_period == 12) %>% select(plot_ID, period_end_month, spi)
names(spi_12)[names(spi_12) == "spi"] <- "spi_12"
spi_24 <- filter(spis, spi_period == 24) %>% select(plot_ID, period_end_month, spi)
names(spi_24)[names(spi_24) == "spi"] <- "spi_24"
growth <- merge(growth, spi_6)
growth <- merge(growth, spi_12)
growth <- merge(growth, spi_24)

# Plot whole dataaset dbh_class versus WD_class versus SPI
group_by(growth, sitecode, SamplingPeriodEnd, plot_ID, WD_class, dbh_class) %>%
    summarise(mean_growth_ann_index=mean(growth_ann_index), mean_spi_12=mean(spi_12)) %>%
    ggplot(aes(mean_spi_12, mean_growth_ann_index, colour=sitecode)) +
    geom_point() + facet_grid(WD_class~dbh_class)

# table(growth$growth_ann_index > 10)
# table(growth$growth_ann_index > 20)
# table(growth$growth_ann_index < -10)
# max(growth$growth_ann_index)
# min(growth$growth_ann_index)

growth <- filter(growth, growth_ann_index < 10)
growth <- filter(growth, growth_ann_index > -10)
growth <- filter(growth, n_days > 200)
growth <- filter(growth, n_days < 550)

growth$growth_ann_index_std <- (growth$growth_ann_index - mean(growth$growth_ann_index)) /  sd(growth$growth_ann_index)

# Ensure that models are not fit to different datasets when performing model 
# selection
options(na.action = "na.fail")

###############################################################################
# SPI 24 with plot-level random effect
m1 <- lmer(growth_ann_index ~ spi_24 + continent + WD_class*spi_24 + dbh_class*spi_24 + (1|sitecode) + (1|plot_ID) + (1|SamplingUnitName), data=growth)
m2 <- lmer(growth_ann_index ~ spi_24 + WD_class + dbh_class + (1|sitecode) + (1|plot_ID) + (1|SamplingUnitName), data=growth)
m3 <- lmer(growth_ann_index ~ spi_24 + dbh_class + (1|sitecode) + (1|plot_ID) + (1|SamplingUnitName), data=growth)
m4 <- lmer(growth_ann_index ~ spi_24 + WD_class + (1|sitecode) + (1|plot_ID) + (1|SamplingUnitName), data=growth)
m5 <- lmer(growth_ann_index ~ dbh_class + WD_class + (1|sitecode) + (1|plot_ID) + (1|SamplingUnitName), data=growth)
m6 <- lmer(growth_ann_index ~ spi_24 + (1|sitecode) + (1|plot_ID) + (1|SamplingUnitName), data=growth)
m7 <- lmer(growth_ann_index ~ WD_class + (1|sitecode) + (1|plot_ID) + (1|SamplingUnitName), data=growth)
m8 <- lmer(growth_ann_index ~ dbh_class + (1|sitecode) + (1|plot_ID) + (1|SamplingUnitName), data=growth)
m9 <- lmer(growth_ann_index ~ continent + (1|sitecode) + (1|plot_ID) + (1|SamplingUnitName), data=growth)
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
# Model mean changes at plot-level
growth_plot_mean <- group_by(growth, sitecode, plot_ID, SamplingPeriodNumber) %>%
    summarize(spi_6=mean(spi_6),
              spi_12=mean(spi_12),
              spi_24=mean(spi_24),
              continent=continent[1],
              growth_ann_index_mean=mean(growth_ann_index),
              growth_ann_mean=mean(growth_ann),
              WD_mean=mean(WD),
              diameter_start_mean=mean(diameter_start),
              growth_ann_index_mean=mean(growth_ann_index))

ggplot(growth_plot_mean) +
    geom_point(aes(spi_24, growth_ann_mean, color=sitecode)) +
    facet_wrap(~sitecode) +
    geom_smooth(aes(spi_24, growth_ann_mean))

ggplot(growth_plot_mean) +
    geom_point(aes(spi_24, growth_ann_index_mean, color=sitecode)) +
    geom_smooth(aes(spi_24, growth_ann_index_mean))

summary(lmer(growth_ann_index_mean ~ spi_6  + continent + (1|plot_ID), data=growth_plot_mean))
summary(lmer(growth_ann_index_mean ~ spi_12  + continent + (1|plot_ID), data=growth_plot_mean))
summary(lmer(growth_ann_index_mean ~ spi_24  + (1|plot_ID), data=growth_plot_mean))

mean_m1 <- lmer(growth_ann_index_mean ~ spi_24 + continent + WD_mean*spi_24 + diameter_start_mean*spi_24 + (1|sitecode), data=growth_plot_mean)
mean_m2 <- lmer(growth_ann_index_mean ~ spi_24 + WD_mean + diameter_start_mean + (1|sitecode), data=growth_plot_mean)
mean_m3 <- lmer(growth_ann_index_mean ~ spi_24 + diameter_start_mean + (1|sitecode), data=growth_plot_mean)
mean_m4 <- lmer(growth_ann_index_mean ~ spi_24 + WD_mean + (1|sitecode), data=growth_plot_mean)
mean_m5 <- lmer(growth_ann_index_mean ~ diameter_start_mean + WD_mean + (1|sitecode), data=growth_plot_mean)
mean_m6 <- lmer(growth_ann_index_mean ~ spi_24 + (1|sitecode), data=growth_plot_mean)
mean_m7 <- lmer(growth_ann_index_mean ~ WD_mean + (1|sitecode), data=growth_plot_mean)
mean_m8 <- lmer(growth_ann_index_mean ~ diameter_start_mean + (1|sitecode), data=growth_plot_mean)
mean_m9 <- lmer(growth_ann_index_mean ~ continent + (1|sitecode), data=growth_plot_mean)
mean_mlist <- list(mean_m1, mean_m2, mean_m3, 
                          mean_m4, mean_m5, mean_m6, 
                          mean_m7, mean_m8, mean_m9)
mean_msel <- model.sel(mean_mlist)
mean_msel

mean_models_95 <- get.models(mean_msel, cumsum(weight) <= .95)
mean_models_95_avg <- model.avg(mean_models_95, beta=TRUE, fit=TRUE)
summary(mean_models_95_avg)

mean_models_coefs <- data.frame(var_name=mean_models_95_avg$term.names,
                                  est=summary(mean_models_95_avg)$coefmat[, 1])
mean_models_coefs_conf <- data.frame(confint(mean_models_95_avg))
names(mean_models_coefs_conf) <- c("lower", "upper")
mean_models_coefs <- cbind(mean_models_coefs, mean_models_coefs_conf)
row.names(mean_models_coefs) <- NULL
mean_models_coefs <- mean_models_coefs[mean_models_coefs$var_name != "(Intercept)", ]

ggplot(mean_models_coefs, aes(est, var_name)) +
    theme_bw(base_size=10) +
    geom_point() + 
    geom_errorbarh(aes(xmin=lower, xmax=upper), height=0.5) +
    geom_vline(aes(xintercept=0), color="grey") +
    #xlim(c(-.15, .15)) +
    xlab("") + ylab("")
ggsave("growth_model_spi_24_mean_coefs.png", width=img_width, 
       height=img_height, dpi=img_dpi)

###############################################################################
# Dredge

global_m <- lmer(growth_ann_index_std ~ spi_24 + I(spi_24^2) + I(spi_24^3) + 
                 continent + WD_class*spi_24 + dbh_class*spi_24 + (1|sitecode) 
                 + (1|plot_ID) + (1|SamplingUnitName), data=growth, REML=FALSE)
# global_m <- lmer(growth_ann_index_std ~ spi_24 + I(spi_24^2) + continent + 
#                  WD_class*spi_24 + dbh_class*spi_24 + (1|sitecode) + 
#                  (1|plot_ID) + (1|SamplingUnitName), data=growth, REML=FALSE)
# global_m <- lmer(growth_ann_index_std ~ spi_24 + continent + WD_class*spi_24 + 
#                  dbh_class*spi_24 + (1|sitecode) + (1|plot_ID) + 
#                  (1|SamplingUnitName), data=growth, REML=FALSE)

# empty <- lmer(growth_ann_index ~ 1 + (1|sitecode) + (1|plot_ID) + 
#               (1|SamplingUnitName), data=growth, REML=FALSE)

#dredged <- dredge(global_m)

cl <- makeCluster(4)
clusterExport(cl, "growth")
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
