# ******************************************************************************
# Model BLSP anomalies using a Bayesian hierarchical model
# 
# Author: Xiaojie Gao
# Date: 2025-05-25
# ******************************************************************************
rm(list = ls())
source("src/base.R")
source("src/03_pheno_variation_explain/hlp_fit_bayes.R")
library(data.table)
library(lubridate)
library(magrittr)
library(rjags)



CalSeasonClimateAnom <- function(mths, dm_dt, season_label) {
    season_dm <- dm_dt[mth %in% mths, .(
        season_Tavg = mean((`tmax (deg c)` + `tmin (deg c)`) / 2),
        season_Prcp = mean(`prcp (mm/day)`),
        season_Srad = mean(`srad (W/m^2)`)
    ), by = c("ID", "year")]

    # Anomalies
    season_dm_ano <- lapply(unique(season_dm$ID), function(id) {
        iddt <- season_dm[ID == id, ]
        ano_dt <- iddt[, .(
            ID, year,
            season_Tavg_ano = scale(season_Tavg, scale = FALSE),
            season_Tavg_scl = scale(season_Tavg, scale = TRUE),
            season_Tavg_mn = mean(season_Tavg),
            season_Prcp_ano = scale(season_Prcp, scale = FALSE),
            season_Prcp_scl = scale(season_Prcp, scale = TRUE),
            season_Prcp_mn = mean(season_Prcp),
            season_Srad_ano = scale(season_Srad, scale = FALSE),
            season_Srad_scl = scale(season_Srad, scale = TRUE),
            season_Srad_mn = mean(season_Srad)
        )]

        return(ano_dt)
    })
    season_dm_ano <- rbindlist(season_dm_ano)

    # Rename columns
    newcolnames <- gsub(".V1", "", colnames(season_dm_ano))
    newcolnames <- gsub("season", season_label, newcolnames)
    colnames(season_dm_ano) <- newcolnames

    return(season_dm_ano)
}


PrepareAtmData <- function(pipedir, blsp_ano, dm_dt, fia_attr_mn_dt) {
    atm_dt_file <- file.path(pipedir, "atm_dt4bayes_new.csv")
    if (file.exists(atm_dt_file)) {
        atm_dt <- fread(atm_dt_file)
        return(atm_dt)
    }

    spr_dm_ano <- CalSeasonClimateAnom(c(3, 4, 5), dm_dt, "spr")
    smr_dm_ano <- CalSeasonClimateAnom(c(6, 7, 8), dm_dt, "smr")
    atm_dm_ano <- CalSeasonClimateAnom(c(9, 10, 11), dm_dt, "atm")


    atm_dt <- merge(
        blsp_ano, spr_dm_ano, 
        by.x = c("ID", "Year"), by.y = c("ID", "year")
    )
    atm_dt <- merge(
        atm_dt, smr_dm_ano, 
        by.x = c("ID", "Year"), by.y = c("ID", "year")
    )
    atm_dt <- merge(
        atm_dt, atm_dm_ano, 
        by.x = c("ID", "Year"), by.y = c("ID", "year")
    )
    atm_dt <- merge(atm_dt, fia_attr_mn_dt, by = "ID")

    atm_dt <- na.omit(atm_dt)

    # out:
    fwrite(atm_dt, atm_dt_file)

    return(atm_dt)
}



pipedir <- "pipe/species_pheno"
dir.create(pipedir, showWarnings = FALSE, recursive = TRUE)


data_list <- LoadAnoData()
blsp_ano <- data_list$blsp_ano
dm_dt <- data_list$dm_dt
fia_attr_mn_dt <- data_list$fia_attr_mn_dt

atm_dt <- PrepareAtmData(pipedir, blsp_ano, dm_dt, fia_attr_mn_dt)
atm_dt <- atm_dt[Year >= 1985,]


# DEBUG ========================================================================
# Sample some sites per species to save some model fitting time
# atm_dt <- lapply(unique(atm_dt$spec), function(sp) {
#     sp_dt <- atm_dt[spec == sp, ]
#     sub_ids <- sample(unique(sp_dt$ID), 50)
#     sub_dt <- sp_dt[ID %in% sub_ids,]
    
#     return(sub_dt)
# }) %>%
#     rbindlist()
# ============== delete above ==================================================



# Data in JAGS format
Y <- atm_dt$EOS_scl
n <- length(Y)
Yr <- atm_dt$Year - 1985 + 1

SOS <- atm_dt$SOS_scl

spr_Tavg <- atm_dt$spr_Tavg_scl
spr_Prcp <- atm_dt$spr_Prcp_scl
spr_Srad <- atm_dt$spr_Srad_scl

smr_Tavg <- atm_dt$smr_Tavg_scl
smr_Prcp <- atm_dt$smr_Prcp_scl
smr_Srad <- atm_dt$smr_Srad_scl
smr_mn_Prcp <- atm_dt$smr_Prcp_mn

atm_Tavg <- atm_dt$atm_Tavg_scl
atm_Prcp <- atm_dt$atm_Prcp_scl
atm_Srad <- atm_dt$atm_Srad_scl

site <- as.numeric(as.factor(atm_dt$ID))
nsite <- length(unique(site))

nstem <- as.vector(scale(unique(atm_dt[, .(ID, nstem)])$nstem))
ba <- as.vector(scale(unique(atm_dt[, .(ID, ba)])$ba))
richness <- as.vector(scale(unique(atm_dt[, .(ID, richness)])$richness))

sp <- as.numeric(as.factor(unique(atm_dt[, .(ID, spec)])$spec))
nsp <- length(unique(sp))



# ~ Model 1 ####
# ~ ----------------------------------------------------------------------------

# We can compare two variants of the model: one uses spring climates; the other
# one uses SOS as the predictor. 

# 1 -----------------
X_mat <- as.matrix(cbind(
    rep(1, n),
    spr_Tavg, spr_Prcp, spr_Srad,
    atm_Tavg, atm_Prcp, atm_Srad,
    smr_Tavg, smr_Prcp, smr_Srad
))
atm_model_1 <- FitBHM(
    atm_mod_per_site,
    data = list(
        X = X_mat,
        Y = Y, n = n, site = site,
        J = nsite, sp = sp, P = ncol(X_mat),
        K = nsp
    ),
    inits = list(
        B = matrix(
            rep(runif(nsite, 0, 5), ncol(X_mat)),
            nrow = ncol(X_mat)
        )
    ),
    var_names = c(
        "sigma_y",
        "B", "sigma_b",
        "gamma", "sigma_g",
        "mu",
        "f", "like"
    )
)
atm_model_1$DIC
# Mean deviance:  252473 
# penalty 4357 
# Penalized deviance: 256830 

atm_model_1$WAIC
# 256899.3

pred_1 <- apply(
    atm_model_1$samp[, grep("^f", colnames(atm_model_1$samp))], 
    2, 
    median
)
base$FitLm(pred_1, Y)
# R2 = 0.23


gamma <- apply(
    atm_model_1$samp[, grep("^gamma", colnames(atm_model_1$samp))],
    2,
    quantile, c(0.025, 0.5, 0.975)
)
colnames(gamma) <- rep(colnames(X_mat), 10)

mu <- apply(
    atm_model_1$samp[, grep("^mu", colnames(atm_model_1$samp))],
    2,
    quantile, c(0.025, 0.5, 0.975)
)
colnames(mu) <- colnames(X_mat)

print(gamma)
print(mu)


# 1a -----------------
X_mat <- as.matrix(cbind(
    rep(1, n),
    SOS,
    atm_Tavg, atm_Prcp, atm_Srad,
    smr_Tavg, smr_Prcp, smr_Srad
))
atm_model_1a <- FitBHM(
    atm_mod_per_site,
    data = list(
        X = X_mat,
        Y = Y, n = n, site = site,
        J = nsite, sp = sp, P = ncol(X_mat),
        K = nsp
    ),
    inits = list(
        B = matrix(
            rep(runif(nsite, 0, 5), ncol(X_mat)),
            nrow = ncol(X_mat)
        )
    ),
    var_names = c(
        "sigma_y",
        "B", "sigma_b",
        "gamma", "sigma_g",
        "mu",
        "f", "like"
    )
)
atm_model_1a$DIC
# Mean deviance:  257355 
# penalty 3858 
# Penalized deviance: 261212 

atm_model_1a$WAIC
# 261393.8


pred_1a <- apply(
    atm_model_1a$samp[, grep("^f", colnames(atm_model_1a$samp))],
    2,
    median
)
vis$GetLmFit(pred_1a, Y)
# R2 = 0.19


gamma <- apply(
    atm_model_1a$samp[, grep("^gamma", colnames(atm_model_1a$samp))],
    2,
    quantile, c(0.025, 0.5, 0.975)
)
mu <- apply(
    atm_model_1a$samp[, grep("^mu", colnames(atm_model_1a$samp))],
    2,
    quantile, c(0.025, 0.5, 0.975)
)
colnames(mu) <- colnames(X_mat)

saveRDS(atm_model_1a, file.path(pipedir, "atm_model_1a.Rds"))



# 1b -----------------
# Another model could test whether there are interactions between variables.
X_mat <- as.matrix(cbind(
    rep(1, n),
    SOS, 
    atm_Tavg, atm_Prcp, atm_Srad,
    smr_Tavg, smr_Prcp, smr_Srad,
    SOS * smr_Prcp
))
atm_model_1b <- FitBHM(
    atm_mod_per_site,
    data = list(
        X = X_mat,
        Y = Y, n = n, site = site,
        J = nsite, sp = sp, P = ncol(X_mat),
        K = nsp
    ),
    inits = list(
        B = matrix(
            rep(runif(nsite, 0, 5), ncol(X_mat)),
            nrow = ncol(X_mat)
        )
    ),
    var_names = c(
        "sigma_y",
        "B", "sigma_b",
        "gamma", "sigma_g",
        "mu",
        "f", "like"
    )
)
atm_model_1b$DIC
# Mean deviance:  256917
# penalty 4220
# Penalized deviance: 261138

atm_model_1b$WAIC
# 261323.1

pred_1b <- apply(
    atm_model_1b$samp[, grep("^f", colnames(atm_model_1b$samp))],
    2,
    median
)
mu <- apply(
    atm_model_1b$samp[, grep("^mu", colnames(atm_model_1b$samp))],
    2,
    quantile, c(0.025, 0.5, 0.975)
)
colnames(mu) <- colnames(X_mat)

base$FitLm(pred_1b, Y)
# R2 = 0.2


# 1c -----------------
# Is there an interaction between smr_Prcp and mn_Prcp?
X_mat <- as.matrix(cbind(
    rep(1, n),
    SOS,
    atm_Tavg, atm_Prcp, atm_Srad,
    smr_Tavg, smr_Prcp, smr_Srad
))
atm_model_1c <- FitBHM(
    atm_mod_per_site,
    data = list(
        X = X_mat,
        Y = Y, n = n, site = site,
        J = nsite, sp = sp, P = ncol(X_mat),
        K = nsp
    ),
    inits = list(
        B = matrix(
            rep(runif(nsite, 0, 5), ncol(X_mat)),
            nrow = ncol(X_mat)
        )
    ),
    var_names = c(
        "sigma_y",
        "B", "sigma_b",
        "gamma", "sigma_g",
        "mu",
        "f", "like"
    )
)
atm_model_1c$DIC


atm_model_1c$WAIC


pred_1c <- apply(
    atm_model_1c$samp[, grep("^f", colnames(atm_model_1c$samp))],
    2,
    median
)
base$FitLm(pred_1c, Y)
# R2 = 0.26

gamma <- apply(
    atm_model_1c$samp[, grep("^gamma", colnames(atm_model_1c$samp))],
    2,
    quantile, c(0.025, 0.5, 0.975)
)
colnames(gamma) <- rep(colnames(X_mat), 10)

sigma_g <- apply(
    atm_model_1c$samp[, grep("^sigma_g", colnames(atm_model_1c$samp))],
    2,
    quantile, c(0.025, 0.5, 0.975)
)


mu <- apply(
    atm_model_1c$samp[, grep("^mu", colnames(atm_model_1c$samp))],
    2,
    quantile, c(0.025, 0.5, 0.975)
)
colnames(mu) <- colnames(X_mat)


mu[2, ] - 1 / sqrt(sigma_g[2, ])

plot(1:ncol(mu), mu[2, ], ylim = c(-0.5, 0.5))
abline(h = 0, lty = 2)
segments(
    1:ncol(mu), mu[1, ],
    1:ncol(mu), mu[3, ]
)
segments(
    1:ncol(mu) + 0.1, mu[2, ] + 1 / sqrt(sigma_g[2, ]),
    1:ncol(mu) + 0.1, mu[2, ] - 1 / sqrt(sigma_g[2, ]),
    col = "blue"
)


# ~ Model 2 ####
# ~ ----------------------------------------------------------------------------






# ~ Autumn model 1 ####
# ~ ----------------------------------------------------------------------------
atm_model1 <- FitBHM(
    atmBHM1,
    data = list(
        Y = Y, n = n, site = site,
        atm_Tavg = atm_Tavg, atm_Prcp = atm_Prcp, atm_Srad = atm_Srad,
        smr_Tavg = smr_Tavg, smr_Prcp = smr_Prcp, smr_Srad = smr_Srad,
        SOS = SOS,
        J = nsite, sp = sp,
        K = nsp
    ),
    inits = list(
        beta1 = runif(nsite, 0, 5),
        beta2 = runif(nsite, 0, 5),
        beta3 = runif(nsite, 0, 5),
        beta4 = runif(nsite, 0, 5),
        beta5 = runif(nsite, 0, 5),
        beta6 = runif(nsite, 0, 5),
        beta7 = runif(nsite, 0, 5),
        beta8 = runif(nsite, 0, 5)
    ),
    var_names = c(
        "beta1", "sigma_b1",
        "beta2", "sigma_b2",
        "beta3", "sigma_b3",
        "beta4", "sigma_b4",
        "beta5", "sigma_b5",
        "beta6", "sigma_b6",
        "beta7", "sigma_b7",
        "beta8", "sigma_b8",
        "gamma1", "sigma_g1",
        "gamma2", "sigma_g2",
        "gamma3", "sigma_g3",
        "gamma4", "sigma_g4",
        "gamma5", "sigma_g5",
        "gamma6", "sigma_g6",
        "gamma7", "sigma_g7",
        "gamma8", "sigma_g8",
        "mu1", "mu2", "mu3", "mu4", "mu5", "mu6", "mu7", "mu8",
        "f", "like"
    )
)



pdf(file.path(pipedir, "Autumn_BHM_fit.pdf"))

pred <- apply(atm_model1$samp[, grep("^f", colnames(atm_model1$samp))], 2, median)
vis_base$smoothScatter_J2(Y, pred, main = "Autumn model 1")



# ~ Autumn model 3 ####
# ~ ----------------------------------------------------------------------------
atm_model3 <- FitBHM(
    atmBHM3,
    data = list(
        Y = Y, n = n, site = site,
        atm_Tavg = atm_Tavg, atm_Prcp = atm_Prcp, atm_Srad = atm_Srad,
        smr_Tavg = smr_Tavg, smr_Prcp = smr_Prcp, smr_Srad = smr_Srad,
        SOS = SOS,
        J = nsite, sp = sp, nstem = nstem, smr_mn_Prcp = smr_mn_Prcp,
        K = nsp
    ),
    inits = list(
        beta1 = runif(nsite, 0, 5),
        beta2 = runif(nsite, 0, 5),
        beta3 = runif(nsite, 0, 5),
        beta4 = runif(nsite, 0, 5),
        beta5 = runif(nsite, 0, 5),
        beta6 = runif(nsite, 0, 5)
    ),
    var_names = c(
        "beta1", "sigma_b1",
        "beta2", "sigma_b2",
        "beta3", "sigma_b3",
        "beta4", "sigma_b4",
        "beta5", "sigma_b5",
        "beta6", "sigma_b6",
        "gamma1", "sigma_g1",
        "gamma2", "sigma_g2",
        "gamma31", "sigma_g31",
        "gamma32", "sigma_g32",
        "gamma33", "sigma_g33",
        "gamma4", "sigma_g4",
        "gamma5", "sigma_g5",
        "gamma6", "sigma_g6",
        "gamma7", "sigma_g7",
        "gamma8", "sigma_g8",
        "mu1", "mu2", "mu31", "mu32", "mu33", "mu4", "mu5", "mu6",
        "f", "like"
    )
)
tend <- Sys.time()
ttake <- tend - tstart
message(round(ttake, 2), " ", units(ttake))

pred <- apply(atm_model3$samp[, grep("^f", colnames(atm_model3$samp))], 2, median)
vis_base$smoothScatter_J2(Y, pred, main = "Autumn model 3")

dev.off()



# ~ Autumn model 4 ####
# ~ ----------------------------------------------------------------------------
# Include site-level predictors
atmBHM4 <- "model {
    # observation-level
    for (i in 1:n) {
        Y[i] ~ dnorm(f[i], sigma_y)
        f[i] <- beta1[site[i]] +
            beta2[site[i]] * atm_Tavg[i] +
            beta3[site[i]] * atm_Prcp[i] +
            beta4[site[i]] * atm_Srad[i] +
            beta5[site[i]] * smr_Tavg[i] +
            beta6[site[i]] * smr_Prcp[i] +
            beta7[site[i]] * smr_Srad[i] +
            beta8[site[i]] * SOS[i]

    }
    # site-level
    for (j in 1:J) {
        beta1[j] ~ dnorm(gamma1[sp[j]], sigma_b1)
        beta2[j] ~ dnorm(gamma2[sp[j]], sigma_b2)
        beta3[j] ~ dnorm(gamma3[sp[j]], sigma_b3)
        beta4[j] ~ dnorm(gamma4[sp[j]], sigma_b4)
        beta5[j] ~ dnorm(gamma5[sp[j]], sigma_b5)

        beta6[j] ~ dnorm(b6[j], sigma_b6)
        b6[j] <-
            gamma61[sp[j]] +
            gamma62[sp[j]] * nstem[j]

        beta7[j] ~ dnorm(gamma7[sp[j]], sigma_b7)
        beta8[j] ~ dnorm(gamma8[sp[j]], sigma_b8)
    }
    # species-level
    for (k in 1:K) {
        gamma1[k] ~ dnorm(mu1, sigma_g1)
        gamma2[k] ~ dnorm(mu2, sigma_g2)
        gamma3[k] ~ dnorm(mu3, sigma_g3)
        gamma4[k] ~ dnorm(mu4, sigma_g4)
        gamma5[k] ~ dnorm(mu5, sigma_g5)

        gamma61[k] ~ dnorm(mu61, sigma_g61)
        gamma62[k] ~ dnorm(mu62, sigma_g62)

        gamma7[k] ~ dnorm(mu7, sigma_g7)
        gamma8[k] ~ dnorm(mu8, sigma_g8)
    }

    # priors
    sigma_y ~ dgamma(0.1, 0.1)
    sigma_b1 ~ dgamma(0.1, 0.1)
    sigma_b2 ~ dgamma(0.1, 0.1)
    sigma_b3 ~ dgamma(0.1, 0.1)
    sigma_b4 ~ dgamma(0.1, 0.1)
    sigma_b5 ~ dgamma(0.1, 0.1)
    sigma_b6 ~ dgamma(0.1, 0.1)
    sigma_b7 ~ dgamma(0.1, 0.1)
    sigma_b8 ~ dgamma(0.1, 0.1)

    sigma_g1 ~ dgamma(0.1, 0.1)
    sigma_g2 ~ dgamma(0.1, 0.1)
    sigma_g3 ~ dgamma(0.1, 0.1)
    sigma_g4 ~ dgamma(0.1, 0.1)
    sigma_g5 ~ dgamma(0.1, 0.1)

    sigma_g61 ~ dgamma(0.1, 0.1)
    sigma_g62 ~ dgamma(0.1, 0.1)

    sigma_g7 ~ dgamma(0.1, 0.1)
    sigma_g8 ~ dgamma(0.1, 0.1)

    mu1 ~ dnorm(0, 0.0001)
    mu2 ~ dnorm(0, 0.0001)
    mu3 ~ dnorm(0, 0.0001)
    mu4 ~ dnorm(0, 0.0001)
    mu5 ~ dnorm(0, 0.0001)

    mu61 ~ dnorm(0, 0.0001)
    mu62 ~ dnorm(0, 0.0001)

    mu7 ~ dnorm(0, 0.0001)
    mu8 ~ dnorm(0, 0.0001)

    # WAIC calculations
    for (i in 1:n) {
        like[i] <- dnorm(Y[i], f[i], sigma_y)
    }
}"

tstart <- Sys.time()
atm_model4 <- FitBHM(
    atmBHM4,
    data = list(
        Y = Y, n = n, site = site,
        atm_Tavg = atm_Tavg, atm_Prcp = atm_Prcp, atm_Srad = atm_Srad,
        smr_Tavg = smr_Tavg, smr_Prcp = smr_Prcp, smr_Srad = smr_Srad,
        SOS = SOS,
        J = nsite, sp = sp, nstem = nstem, richness = richness, ba = ba,
        K = nsp
    ),
    inits = list(
        beta1 = runif(nsite, 0, 5),
        beta2 = runif(nsite, 0, 5),
        beta3 = runif(nsite, 0, 5),
        beta4 = runif(nsite, 0, 5),
        beta5 = runif(nsite, 0, 5),
        beta6 = runif(nsite, 0, 5),
        beta7 = runif(nsite, 0, 5),
        beta8 = runif(nsite, 0, 5)
    ),
    var_names = c(
        "beta1", "sigma_b1",
        "beta2", "sigma_b2",
        "beta3", "sigma_b3",
        "beta4", "sigma_b4",
        "beta5", "sigma_b5",
        "beta6", "sigma_b6",
        "gamma1", "sigma_g1",
        "gamma2", "sigma_g2",
        "gamma3", "sigma_g3",
        "gamma4", "sigma_g4",
        "gamma5", "sigma_g5",
        "gamma61", "sigma_g61",
        "gamma62", "sigma_g62",
        "gamma7", "sigma_g7",
        "gamma8", "sigma_g8",
        "mu1", "mu2", "mu3", "mu4", "mu5", "mu61", "mu62",
        "mu7", "mu8",
        "f", "like"
    )
)
tend <- Sys.time()
ttake <- tend - tstart
message(round(ttake, 2), " ", units(ttake))

pred <- apply(atm_model4$samp[, grep("^f", colnames(atm_model4$samp))], 2, median)
vis_base$smoothScatter_J2(Y, pred, main = "Autumn model 4")

atm_model1$DIC
atm_model1$WAIC
atm_model3$DIC
atm_model3$WAIC
atm_model4$DIC
atm_model4$WAIC



# ~ Atumn model 5 ####
# ~ ----------------------------------------------------------------------------
# Test SOS-Smr_Prec interaction
atmBHM5 <- "model {
    # observation-level
    for (i in 1:n) {
        Y[i] ~ dnorm(f[i], sigma_y[site[i]])
        f[i] <- beta1[site[i]] +
            beta2[site[i]] * atm_Tavg[i] +
            beta3[site[i]] * atm_Prcp[i] +
            beta4[site[i]] * atm_Srad[i] +
            beta5[site[i]] * smr_Tavg[i] +
            beta6[site[i]] * smr_Prcp[i] +
            beta7[site[i]] * smr_Srad[i] +
            beta8[site[i]] * SOS[i] +
            beta9[site[i]] * SOS[i] * atm_Tavg[i]
    }
    # site-level
    for (j in 1:J) {
        beta1[j] ~ dnorm(gamma1[sp[j]], sigma_b1)
        beta2[j] ~ dnorm(gamma2[sp[j]], sigma_b2)
        beta3[j] ~ dnorm(gamma3[sp[j]], sigma_b3)
        beta4[j] ~ dnorm(gamma4[sp[j]], sigma_b4)
        beta5[j] ~ dnorm(gamma5[sp[j]], sigma_b5)
        beta6[j] ~ dnorm(gamma6[sp[j]], sigma_b6)
        beta7[j] ~ dnorm(gamma7[sp[j]], sigma_b7)
        beta8[j] ~ dnorm(gamma8[sp[j]], sigma_b8)
        beta9[j] ~ dnorm(gamma9[sp[j]], sigma_b9)

        sigma_y[j] ~ dgamma(0.1, 0.1)
    }
    # species-level
    for (k in 1:K) {
        gamma1[k] ~ dnorm(mu1, sigma_g1)
        gamma2[k] ~ dnorm(mu2, sigma_g2)
        gamma3[k] ~ dnorm(mu3, sigma_g3)
        gamma4[k] ~ dnorm(mu4, sigma_g4)
        gamma5[k] ~ dnorm(mu5, sigma_g5)
        gamma6[k] ~ dnorm(mu6, sigma_g6)
        gamma7[k] ~ dnorm(mu7, sigma_g7)
        gamma8[k] ~ dnorm(mu8, sigma_g8)
        gamma9[k] ~ dnorm(mu9, sigma_g9)
    }

    # priors
    
    sigma_b1 ~ dgamma(0.1, 0.1)
    sigma_b2 ~ dgamma(0.1, 0.1)
    sigma_b3 ~ dgamma(0.1, 0.1)
    sigma_b4 ~ dgamma(0.1, 0.1)
    sigma_b5 ~ dgamma(0.1, 0.1)
    sigma_b6 ~ dgamma(0.1, 0.1)
    sigma_b7 ~ dgamma(0.1, 0.1)
    sigma_b8 ~ dgamma(0.1, 0.1)
    sigma_b9 ~ dgamma(0.1, 0.1)

    sigma_g1 ~ dgamma(0.1, 0.1)
    sigma_g2 ~ dgamma(0.1, 0.1)
    sigma_g3 ~ dgamma(0.1, 0.1)
    sigma_g4 ~ dgamma(0.1, 0.1)
    sigma_g5 ~ dgamma(0.1, 0.1)
    sigma_g6 ~ dgamma(0.1, 0.1)
    sigma_g7 ~ dgamma(0.1, 0.1)
    sigma_g8 ~ dgamma(0.1, 0.1)
    sigma_g9 ~ dgamma(0.1, 0.1)

    mu1 ~ dnorm(0, 0.0001)
    mu2 ~ dnorm(0, 0.0001)
    mu3 ~ dnorm(0, 0.0001)
    mu4 ~ dnorm(0, 0.0001)
    mu5 ~ dnorm(0, 0.0001)
    mu6 ~ dnorm(0, 0.0001)
    mu7 ~ dnorm(0, 0.0001)
    mu8 ~ dnorm(0, 0.0001)
    mu9 ~ dnorm(0, 0.0001)

    # WAIC calculations
    for (i in 1:n) {
        like[i] <- dnorm(Y[i], f[i], sigma_y[site[i]])
    }
}"


# ! Takes 1.2 hours to run
tstart <- Sys.time()
atm_model5 <- FitBHM(
    atmBHM5,
    data = list(
        Y = Y, n = n, site = site,
        atm_Tavg = atm_Tavg, atm_Prcp = atm_Prcp, atm_Srad = atm_Srad,
        smr_Tavg = smr_Tavg, smr_Prcp = smr_Prcp, smr_Srad = smr_Srad,
        SOS = SOS,
        J = nsite, sp = sp, nstem = nstem, richness = richness, ba = ba,
        K = nsp
    ),
    inits = list(
        beta1 = runif(nsite, 0, 5),
        beta2 = runif(nsite, 0, 5),
        beta3 = runif(nsite, 0, 5),
        beta4 = runif(nsite, 0, 5),
        beta5 = runif(nsite, 0, 5),
        beta6 = runif(nsite, 0, 5),
        beta7 = runif(nsite, 0, 5),
        beta8 = runif(nsite, 0, 5),
        beta9 = runif(nsite, 0, 5)
    ),
    var_names = c(
        "beta1", "sigma_b1",
        "beta2", "sigma_b2",
        "beta3", "sigma_b3",
        "beta4", "sigma_b4",
        "beta5", "sigma_b5",
        "beta6", "sigma_b6",
        "beta7", "sigma_b7",
        "beta8", "sigma_b8",
        "beta9", "sigma_b9",
        "gamma1", "sigma_g1",
        "gamma2", "sigma_g2",
        "gamma3", "sigma_g3",
        "gamma4", "sigma_g4",
        "gamma5", "sigma_g5",
        "gamma6", "sigma_g6",
        "gamma7", "sigma_g7",
        "gamma8", "sigma_g8",
        "gamma9", "sigma_g9",
        "mu1", "mu2", "mu3", "mu4", "mu5", "mu6",
        "mu7", "mu8", "mu9",
        "f", "like"
    )
)
tend <- Sys.time()
ttake <- tend - tstart
saveRDS(atm_model5, "zzz_atm_model5.Rds")
message(round(ttake, 2), " ", units(ttake))


pred <- apply(atm_model5$samp[, grep("^f", colnames(atm_model5$samp))], 2, median)
vis_base$smoothScatter_J2(pred, Y, main = "Autumn model 5")

atm_model1$DIC
atm_model1$WAIC
atm_model3$DIC
atm_model3$WAIC
atm_model4$DIC
atm_model4$WAIC
atm_model5$DIC
atm_model5$WAIC


gamma <- apply(
    atm_model5$samp[, grep(paste0("^gamma", 9), colnames(atm_model5$samp))],
    2,
    quantile, c(0.025, 0.5, 0.975)
)




# ~ Model 6 ####
# ~ ----------------------------------------------------------------------------
atmBHM6 <- "model {
    # observation-level
    for (i in 1:n) {
        Y[i] ~ dnorm(f[i], sigma_y)
        f[i] <- inprod(X[i,], B[, site[i]])
    }
    # site-level
    for (p in 1:P) {
        for (j in 1:J) {
            B[p, j] ~ dnorm(gamma[p, sp[j]], sigma_b[p])
        }
    }
    # species-level
    for (p in 1:P) {
        for (k in 1:K) {
            gamma[p, k] ~ dnorm(mu[p], sigma_g[p])
        }
    }

    # priors
    sigma_y ~ dgamma(0.1, 0.1)
    for (p in 1:P) {
        sigma_b[p] ~ dgamma(0.1, 0.1)
        mu[p] ~ dnorm(0, 0.0001)
        sigma_g[p] ~ dgamma(0.1, 0.1)
    }

    # WAIC calculations
    for (i in 1:n) {
        like[i] <- dnorm(Y[i], f[i], sigma_y)
    }
}"

X_mat <- as.matrix(cbind(
    rep(1, n),
    atm_Tavg, atm_Prcp, atm_Srad,
    smr_Tavg, smr_Prcp, smr_Srad,
    SOS, SOS * atm_Tavg
))

tstart <- Sys.time()
atm_model6 <- FitBHM(
    atmBHM6,
    data = list(
        Y = Y, n = n, site = site,
        X = X_mat,
        J = nsite, sp = sp, P = ncol(X_mat),
        K = nsp
    ),
    inits = list(
        B = matrix(
            rep(runif(nsite, 0, 5), ncol(X_mat)), 
            nrow = ncol(X_mat)
        )
    ),
    var_names = c(
        "sigma_y",
        "B", "sigma_b",
        "gamma", "sigma_g",
        "mu",
        "f", "like"
    )
)
tend <- Sys.time()
ttake <- tend - tstart
message(round(ttake, 2), " ", units(ttake))

pred <- apply(atm_model6$samp[, grep("^f", colnames(atm_model6$samp))], 2, median)
vis_base$smoothScatter_J2(pred, Y, main = "Autumn model 6")


