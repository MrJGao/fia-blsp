# ******************************************************************************
# Train a random forest model to investigate the relative effects of climate,
# FIA site characteristics, and other explanatory variables on driving the
# spring and autumn phenology trends.
# 
# Author: Xiaojie Gao
# Date: 2024-12-16
# ******************************************************************************
rm(list = ls())
source("src/base.R")
library(data.table)
library(magrittr)
library(randomForest)



site_trend_dt <- fread("pipe/fia_dom_blsp_trend_ols.csv")[Type == "deci"]

temp_trend <- fread("pipe/fia_domspec_daymet_temp_trends_ols.csv")[, vari := "temp"]
temp_trend <- dcast(temp_trend, ID ~ season + vari, value.var = c("slp", "pval"))

prcp_trend <- fread("pipe/fia_domspec_daymet_prcp_trends_ols.csv")[, vari := "prcp"]
prcp_trend <- dcast(prcp_trend, ID ~ season + vari, value.var = c("slp", "pval"))

srad_trend <- fread("pipe/fia_domspec_daymet_srad_trends_ols.csv")[, vari := "srad"]
srad_trend <- dcast(srad_trend, ID ~ season + vari, value.var = c("slp", "pval"))

vp_trend <- fread("pipe/fia_domspec_daymet_vp_trends_ols.csv")[, vari := "vp"]
vp_trend <- dcast(vp_trend, ID ~ season + vari, value.var = c("slp", "pval"))

# Merge w/ climate trends
trend_dt <- Reduce(
    f = function(x, y) {
        merge(x, y, by = "ID")
    },
    x = list(site_trend_dt, temp_trend, prcp_trend, srad_trend, vp_trend)
)

# Merge w/ FIA attributes
fia_attr_dt <- fread("data/fia_attr_dt_r9.csv")

# Compute species richness
fia_spec_cols <- names(fia_attr_dt) %>%
    grep("spBAm2ha*.*_live", ., value = TRUE)

spec_rich <- lapply(fia_spec_cols, function(spec_col) {
    spec_ba <- fia_attr_dt[, get(spec_col)]
    spec_ba <- ifelse(is.na(spec_ba), 0, 1)
    return(spec_ba)
})
spec_rich <- do.call(cbind, spec_rich)
spec_rich <- apply(spec_rich, 1, sum)

fia_attr_dt[, richness := spec_rich]

fia_attr_dt <- fia_attr_dt[, .(
    ID = paste0("X", X),
    MEASYEAR,
    nstem = plotGroupnStemsha_live,
    ba = plotGroupBAm2ha_live,
    richness
)]

# Get mean record
fia_attr_mn_dt <- fia_attr_dt[, 
    .SD[, .(
        nstem = median(as.numeric(nstem)), 
        ba = mean(ba), 
        richness = median(richness))
    ], 
    by = ID
]


trend_dt <- merge(trend_dt, fia_attr_mn_dt, by = "ID")



# # Classify slopes into `sig. positive`, `insig. positive`, `sig. negative`,
# # `insig. negative`, and `no trend`
# ClassifyTrendClass <- function(slp, pval) {
#     trendclass <- rep("no_trend", length(slp))
#     trendclass[slp > 0] <- "pos"
#     trendclass[slp > 0 & pval < 0.05] <- "sig_pos"
#     trendclass[slp < 0] <- "neg"
#     trendclass[slp < 0 & pval < 0.05] <- "sig_neg"
    
#     return(trendclass)
# }
# spr_nona_dt <- trend_dt[, .(ID, spec, nstem, richness)]
# spr_nona_dt$spr_slp <- ClassifyTrendClass(trend_dt$spr_slp, trend_dt$spr_pval)
# spr_nona_dt$slp_spr_temp <- ClassifyTrendClass(trend_dt$slp_spr_temp, trend_dt$pval_spr_temp)
# spr_nona_dt$slp_spr_prcp <- ClassifyTrendClass(trend_dt$slp_spr_prcp, trend_dt$pval_spr_prcp)
# spr_nona_dt$slp_spr_srad <- ClassifyTrendClass(trend_dt$slp_spr_srad, trend_dt$pval_spr_srad)
# spr_nona_dt <- na.omit(spr_nona_dt)

# spr_mod <- randomForest(
#     as.factor(spr_slp) ~ slp_spr_temp + slp_spr_prcp + slp_spr_srad +
#         nstem + richness +
#         spec,
#     data = spr_nona_dt,
#     importance = TRUE
# )
# varImpPlot(spr_mod)

# summary(lm(spr_slp ~ spec + slp_spr_prcp, data = spr_nona_dt))


# ~ Fit model ####
# ~ ----------------------------------------------------------------------------

# ~ Spring ----------------------------------------

spr_nona_dt <- na.omit(
    trend_dt[, .(
        spr_slp, slp_spr_temp, slp_spr_prcp, slp_spr_srad, nstem, richness, spec
    )]
)
# # Tune random forest
# ntree_rg <- seq(100, 2000, by = 100)
# tune_dt <- NULL
# pb <- txtProgressBar(min = 0, max = 100, style = 3)
# for (i in seq_along(ntree_rg)) {
#     the_ntree <- ntree_rg[i]
#     tunerf <- tuneRF(
#         spr_nona_dt[, .(
#             spr_slp, slp_spr_temp, slp_spr_prcp, slp_spr_srad, nstem, richness, spec
#         )],
#         spr_nona_dt[, spr_slp],
#         mtryStart = 6,
#         improve = 0.00001,
#         stepFactor = 1.5,
#         ntreeTry = the_ntree,
#         trace = FALSE, plot = FALSE
#     )

#     tune_dt <- rbind(tune_dt, data.table(
#         ntree = the_ntree,
#         mtry = tunerf[which.min(tunerf[[2]]), 1],
#         ooberr = tunerf[which.min(tunerf[[2]]), 2]
#     ))
    
#     setTxtProgressBar(pb, i * 100 / length(ntree_rg)) # update progress
# }
# close(pb)

# tune_dt[which.min(ooberr), ]

spr_res_dt <- NULL
for (i in 1:100) {
    spr_mod <- randomForest(
        spr_slp ~ slp_spr_temp + slp_spr_prcp + slp_spr_srad + spec,
        data = spr_nona_dt,
        mtry = 1,
        ntree = 500,
        importance = TRUE
    )
    spr_imp <- importance(spr_mod)
    spr_res_dt <- cbind(spr_res_dt, spr_imp[, 1])
}
# Summarize feature importance
spr_imp <- apply(spr_res_dt, 1, mean)
names(spr_imp) <- c(
    "Spr_Temp", "Spr_Prcp", "Spr_Srad", "Species"
)
spr_imp <- sort(spr_imp, decreasing = TRUE)

# print(spr_mod)
# varImpPlot(spr_mod)

# spr_mod2 <- randomForest(
#     spr_slp ~ 
#         slp_spr_temp + spec + slp_spr_prcp,
#     data = spr_nona_dt,
#     mtry = 1, 
#     ntree = 500,
#     importance = TRUE
# )
# print(spr_mod2)

# The correlation between the predictors
# GGally::ggpairs(
#     spr_nona_dt[, .(
#         slp_spr_temp, slp_spr_prcp, slp_spr_srad,
#         nstem, richness, DEM
#     )]
# )

# For spring, it's obvious that `slp_spr_temp` had the largest effect, then,
# species, then, `slp_spr_prcp`. The sub model can achieve 95% of the accuracy
# of the full model. And, the predictors are not highly correlated.




# vis_base$LinXY(spr_nona_dt$spr_slp, predict(spr_mod))



# ~ Fall ---------------------------------------- 
atm_nona_dt <- na.omit(trend_dt)

# # Tune random forest
# tune_dt <- NULL
# for (the_ntree in seq(100, 2000, by = 100)) {
#     tunerf <- tuneRF(
#         atm_nona_dt[, .(
#             slp_smr_temp, slp_atm_srad, slp_smr_prcp,
#             nstem, richness, spec, spr_slp
#         )],
#         atm_nona_dt[, atm_slp],
#         mtryStart = 7,
#         improve = 0.0001,
#         stepFactor = 1.5,
#         ntreeTry = the_ntree
#         trace = FALSE, plot = FALSE
#     )
    
#     tune_dt <- rbind(tune_dt, data.table(
#         ntree = the_ntree,
#         mtry = tunerf[which.min(tunerf[[2]]), 1],
#         ooberr = tunerf[which.min(tunerf[[2]]), 2]
#     ))
# }

# tune_dt[which.min(ooberr), ]


atm_res_dt <- NULL
for (i in 1:100) {
    atm_mod <- randomForest(
        atm_slp ~
            slp_smr_prcp + slp_smr_temp + slp_smr_srad +
            slp_atm_temp + slp_atm_prcp + slp_atm_srad +
            spec +
            nstem +
            spr_slp,
        data = atm_nona_dt,
        mtry = 1,
        ntree = 500,
        importance = TRUE
    )
    atm_imp <- importance(atm_mod)
    atm_res_dt <- cbind(atm_res_dt, atm_imp[, 1])
}
# Summarize feature importance
atm_imp <- apply(atm_res_dt, 1, mean)
names(atm_imp) <- c(
    "Smr_Prcp", "Smr_Temp", "Smr_Srad",
    "Atm_Temp", "Atm_Prcp", "Atm_Srad",
    "Species", "Density", "SOS"
)
atm_imp <- sort(atm_imp, decreasing = TRUE)
# barplot(atm_imp)


# vis_base$LinXY(atm_nona_dt$atm_slp, predict(atm_mod))








# plot(predict(atm_mod), atm_nona_dt$atm_slp)
# atm_mod$importance

# plot(atm_nona_dt[, .(nstem, atm_slp)])

# lapply(unique(atm_nona_dt$spec), function(sp) {
#     spec_dt <- atm_nona_dt[spec == sp]
#     ppcor::pcor.test(
#         spec_dt$ba, spec_dt$atm_slp,
#         z = spec_dt[, .(slp_atm_temp, slp_atm_prcp)],
#         method = "spearman"
#     )

# })


# pdf("zzz.pdf", width = 20, height = 20)
# pairs(atm_nona_dt[, .(
#     spr_slp, atm_slp, gsl_slp, lon, 
#     slp_ann_temp, slp_atm_temp, slp_spr_temp, 
#     slp_ann_prcp, slp_atm_prcp, slp_spr_prcp,
#     slp_ann_srad, slp_atm_srad, slp_spr_srad,
#     slp_ann_vp, slp_atm_vp, slp_spr_vp,
#     nstem, ba, richness
# )])
# dev.off()

