# ******************************************************************************
# Final models for spring and autumn phenology
# 
# Author: Xiaojie Gao
# Date: 2025-01-04
# ******************************************************************************
rm(list = ls())
source("src/base.R")
library(data.table)
library(lubridate)
library(magrittr)
library(rjags)



# This function does not compute DIC and WAIC to save some computing time
FitBHM2 <- function(
    model_str, data, inits, var_names,
    nchain = 2, ifplot = FALSE, quiet = FALSE, progress.bar = "text") {
    message("Init model...")
    mod <- jags.model(
        textConnection(model_str),
        data = data, inits = inits, n.chains = nchain,
        quiet = quiet
    )
    message("Update model...")
    update(mod, 5000, progress.bar = progress.bar)
    message("Sample parameters...")
    samp <- coda.samples(
        mod,
        variable.names = var_names,
        n.iter = 5000, thin = 5,
        progress.bar = progress.bar
    )
    samp <- Reduce(rbind, samp)

    return(list(
        samp = samp
    ))
}



pipedir <- "pipe/species_pheno"
dir.create(pipedir, showWarnings = FALSE, recursive = TRUE)



# ~ Spring ####
# ~ ----------------------------------------------------------------------------
spr_dt <- fread(file.path(pipedir, "spr_dt4bayes.csv"))


# Data in JAGS format
Y <- spr_dt$SOS_scl.V1
n <- length(Y)

Tavg <- spr_dt$spr_Tavg_scl.V1
Prcp <- spr_dt$spr_Prcp_scl.V1
Srad <- spr_dt$spr_Srad_scl.V1

site <- as.numeric(as.factor(spr_dt$ID))
nsite <- length(unique(site))

nstem <- as.vector(scale(unique(spr_dt[, .(ID, nstem)])$nstem))
ba <- as.vector(scale(unique(spr_dt[, .(ID, ba)])$ba))
richness <- as.vector(scale(unique(spr_dt[, .(ID, richness)])$richness))

sp <- as.numeric(as.factor(unique(spr_dt[, .(ID, spec)])$spec))
nsp <- length(unique(sp))


# This model include `Tavg`, `Prcp`, and `Srad`
sprBHM3 <- "model {
    # observation-level
    for (i in 1:n) {
        Y[i] ~ dnorm(f[i], sigma_y)
        f[i] <- beta1[site[i]] +
            beta2[site[i]] * Tavg[i] +
            beta3[site[i]] * Prcp[i] +
            beta4[site[i]] * Srad[i]
    }
    # site-level
    for (j in 1:J) {
        beta1[j] ~ dnorm(gamma1[sp[j]], sigma_b1)
        beta2[j] ~ dnorm(gamma2[sp[j]], sigma_b2)
        beta3[j] ~ dnorm(gamma3[sp[j]], sigma_b3)
        beta4[j] ~ dnorm(gamma4[sp[j]], sigma_b4)
    }
    # species-level
    for (k in 1:K) {
        gamma1[k] ~ dnorm(mu1, sigma_g1)
        gamma2[k] ~ dnorm(mu2, sigma_g2)
        gamma3[k] ~ dnorm(mu3, sigma_g3)
        gamma4[k] ~ dnorm(mu4, sigma_g4)
    }

    # priors
    sigma_y ~ dgamma(0.1, 0.1)
    sigma_b1 ~ dgamma(0.1, 0.1)
    sigma_b2 ~ dgamma(0.1, 0.1)
    sigma_b3 ~ dgamma(0.1, 0.1)
    sigma_b4 ~ dgamma(0.1, 0.1)

    sigma_g1 ~ dgamma(0.1, 0.1)
    sigma_g2 ~ dgamma(0.1, 0.1)
    sigma_g3 ~ dgamma(0.1, 0.1)
    sigma_g4 ~ dgamma(0.1, 0.1)

    mu1 ~ dnorm(0, 0.0001)
    mu2 ~ dnorm(0, 0.0001)
    mu3 ~ dnorm(0, 0.0001)
    mu4 ~ dnorm(0, 0.0001)

    # WAIC calculations
    for (i in 1:n) {
        like[i] <- dnorm(Y[i], f[i], sigma_y)
    }
}"

tstart <- Sys.time()
spr_model3 <- FitBHM2(
    sprBHM3,
    data = list(
        Y = Y, n = n, site = site,
        Tavg = Tavg, Prcp = Prcp, Srad = Srad,
        J = nsite, sp = sp,
        K = nsp
    ),
    inits = list(
        beta1 = runif(nsite, 0, 1000),
        beta2 = runif(nsite, 0, 1000),
        beta3 = runif(nsite, 0, 1000)
    ),
    var_names = c(
        "beta1", "sigma_b1",
        "beta2", "sigma_b2",
        "beta3", "sigma_b3",
        "beta4", "sigma_b4",
        "gamma1", "sigma_g1",
        "gamma2", "sigma_g2",
        "gamma3", "sigma_g3",
        "gamma4", "sigma_g4",
        "mu1", "mu2", "mu3", "mu4",
        "f", "like"
    )
)
tend <- Sys.time()
ttake <- tend - tstart
message(round(ttake, 2), " ", units(ttake))

pred <- apply(spr_model3$samp[, grep("^f", colnames(spr_model3$samp))], 2, median)
vis_base$smoothScatter_J2(Y, pred, main = "Spring model 3")

saveRDS(spr_model3, file.path(pipedir, "spr_model3.Rds"))



# ~ Autumn model ####
# ~ ----------------------------------------------------------------------------
atm_dt <- fread(file.path(pipedir, "atm_dt4bayes.csv"))

# Data in JAGS format
Y <- atm_dt$EOS_scl.V1
n <- length(Y)

SOS <- atm_dt$SOS_scl.V1

smr_Tavg <- atm_dt$smr_Tavg_scl.V1
smr_Prcp <- atm_dt$smr_Prcp_scl.V1
smr_Srad <- atm_dt$smr_Srad_scl.V1

smr_mn_Prcp <- atm_dt$smr_Prcp_mn

atm_Tavg <- atm_dt$atm_Tavg_scl.V1
atm_Prcp <- atm_dt$atm_Prcp_scl.V1
atm_Srad <- atm_dt$atm_Srad_scl.V1

site <- as.numeric(as.factor(atm_dt$ID))
nsite <- length(unique(site))

nstem <- as.vector(scale(unique(atm_dt[, .(ID, nstem)])$nstem))
ba <- as.vector(scale(unique(atm_dt[, .(ID, ba)])$ba))
richness <- as.vector(scale(unique(atm_dt[, .(ID, richness)])$richness))

sp <- as.numeric(as.factor(unique(atm_dt[, .(ID, spec)])$spec))
nsp <- length(unique(sp))


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
atm_model4 <- FitBHM2(
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
        "beta7", "sigma_b7",
        "beta8", "sigma_b8",
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

saveRDS(atm_model4, file.path(pipedir, "atm_model4_new.Rds"))
# saveRDS(atm_model4, file.path(pipedir, "atm_model4.Rds"))



# ~ Autumn model 2 ####
# ~ ----------------------------------------------------------------------------
atm_dt <- fread(file.path(pipedir, "atm_dt4bayes.csv"))

# Data in JAGS format
Y <- atm_dt$EOS_scl.V1
n <- length(Y)

SOS <- atm_dt$SOS_scl.V1

m8_Tavg <- atm_dt$m8_Tavg_scl.V1
m8_Prcp <- atm_dt$m8_Prcp_scl.V1
m8_Srad <- atm_dt$m8_Srad_scl.V1

m9_Tavg <- atm_dt$m9_Tavg_scl.V1
m9_Prcp <- atm_dt$m9_Prcp_scl.V1
m9_Srad <- atm_dt$m9_Srad_scl.V1

m10_Tavg <- atm_dt$m10_Tavg_scl.V1
m10_Prcp <- atm_dt$m10_Prcp_scl.V1
m10_Srad <- atm_dt$m10_Srad_scl.V1

m8_mn_Prcp <- atm_dt$m8_Prcp_mn

site <- as.numeric(as.factor(atm_dt$ID))
nsite <- length(unique(site))

nstem <- as.vector(scale(unique(atm_dt[, .(ID, nstem)])$nstem))
ba <- as.vector(scale(unique(atm_dt[, .(ID, ba)])$ba))
richness <- as.vector(scale(unique(atm_dt[, .(ID, richness)])$richness))

sp <- as.numeric(as.factor(unique(atm_dt[, .(ID, spec)])$spec))
nsp <- length(unique(sp))


# Include site-level predictors
atmBHM4 <- "model {
    # observation-level
    for (i in 1:n) {
        Y[i] ~ dnorm(f[i], sigma_y)
        f[i] <- beta1[site[i]] +
            beta2[site[i]] * m8_Tavg[i] +
            beta3[site[i]] * m8_Prcp[i] +
            beta4[site[i]] * m8_Srad[i] +
            beta5[site[i]] * m9_Tavg[i] +
            beta6[site[i]] * m9_Prcp[i] +
            beta7[site[i]] * m9_Srad[i] +
            beta8[site[i]] * m10_Tavg[i] +
            beta9[site[i]] * m10_Prcp[i] +
            beta10[site[i]] * m10_Srad[i] +
            beta11[site[i]] * SOS[i]
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
        beta9[j] ~ dnorm(gamma9[sp[j]], sigma_b9)
        beta10[j] ~ dnorm(gamma10[sp[j]], sigma_b10)
        beta11[j] ~ dnorm(gamma11[sp[j]], sigma_b11)
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
        gamma9[k] ~ dnorm(mu9, sigma_g9)
        gamma10[k] ~ dnorm(mu10, sigma_g10)
        gamma11[k] ~ dnorm(mu11, sigma_g11)
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
    sigma_b9 ~ dgamma(0.1, 0.1)
    sigma_b10 ~ dgamma(0.1, 0.1)
    sigma_b11 ~ dgamma(0.1, 0.1)

    sigma_g1 ~ dgamma(0.1, 0.1)
    sigma_g2 ~ dgamma(0.1, 0.1)
    sigma_g3 ~ dgamma(0.1, 0.1)
    sigma_g4 ~ dgamma(0.1, 0.1)
    sigma_g5 ~ dgamma(0.1, 0.1)

    sigma_g61 ~ dgamma(0.1, 0.1)
    sigma_g62 ~ dgamma(0.1, 0.1)

    sigma_g7 ~ dgamma(0.1, 0.1)
    sigma_g8 ~ dgamma(0.1, 0.1)
    sigma_g9 ~ dgamma(0.1, 0.1)
    sigma_g10 ~ dgamma(0.1, 0.1)
    sigma_g11 ~ dgamma(0.1, 0.1)

    mu1 ~ dnorm(0, 0.0001)
    mu2 ~ dnorm(0, 0.0001)
    mu3 ~ dnorm(0, 0.0001)
    mu4 ~ dnorm(0, 0.0001)
    mu5 ~ dnorm(0, 0.0001)

    mu61 ~ dnorm(0, 0.0001)
    mu62 ~ dnorm(0, 0.0001)

    mu7 ~ dnorm(0, 0.0001)
    mu8 ~ dnorm(0, 0.0001)
    mu9 ~ dnorm(0, 0.0001)
    mu10 ~ dnorm(0, 0.0001)
    mu11 ~ dnorm(0, 0.0001)

    # WAIC calculations
    for (i in 1:n) {
        like[i] <- dnorm(Y[i], f[i], sigma_y)
    }
}"

tstart <- Sys.time()
atm_model4 <- FitBHM2(
    atmBHM4,
    data = list(
        Y = Y, n = n, site = site,
        m8_Tavg = m8_Tavg, m8_Prcp = m8_Prcp, m8_Srad = m8_Srad,
        m9_Tavg = m9_Tavg, m9_Prcp = m9_Prcp, m9_Srad = m9_Srad,
        m10_Tavg = m10_Tavg, m10_Prcp = m10_Prcp, m10_Srad = m10_Srad,
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
        beta9 = runif(nsite, 0, 5),
        beta10 = runif(nsite, 0, 5)
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
        "beta10", "sigma_b10",
        "gamma1", "sigma_g1",
        "gamma2", "sigma_g2",
        "gamma3", "sigma_g3",
        "gamma4", "sigma_g4",
        "gamma5", "sigma_g5",
        "gamma61", "sigma_g61",
        "gamma62", "sigma_g62",
        "gamma7", "sigma_g7",
        "gamma8", "sigma_g8",
        "gamma9", "sigma_g9",
        "gamma10", "sigma_g10",
        "mu1", "mu2", "mu3", "mu4", "mu5", "mu61", "mu62",
        "mu7", "mu8", "mu9", "mu10",
        "f", "like"
    )
)
tend <- Sys.time()
ttake <- tend - tstart
message(round(ttake, 2), " ", units(ttake))

pred <- apply(atm_model4$samp[, grep("^f", colnames(atm_model4$samp))], 2, median)
vis_base$smoothScatter_J2(Y, pred, main = "Autumn model 4")

saveRDS(atm_model4, file.path(pipedir, "atm_model44.Rds"))




