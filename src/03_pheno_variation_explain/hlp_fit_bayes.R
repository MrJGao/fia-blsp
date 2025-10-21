# ******************************************************************************
# Functions for fit and extracting coefs from Bayesian models
#
# Author: Xiaojie Gao
# Date: 2025-05-25
# ******************************************************************************
library(data.table)
library(lubridate)
library(magrittr)


LoadAnoData <- function() {
    # ~ Load BLSP and calculate anomalies ------------------------------------

    blsp_fit_dt <- fread("data/fia_blsp_fit.csv")[, .(
        ID, spec, x, y, Year,
        MidGreenup,
        MidGreenup_unc = MidGreenup_upr - MidGreenup_lwr,
        MidGreendown,
        MidGreendown_unc = MidGreendown_upr - MidGreendown_lwr
    )]
    # Remove records with high uncertainty
    blsp_fit_dt[MidGreenup_unc > 30, MidGreenup := NA]
    blsp_fit_dt[MidGreendown_unc > 30, MidGreendown := NA]

    # Calculate GSL
    blsp_fit_dt[, GSL := MidGreendown - MidGreenup]

    blsp_fit_dt <- na.omit(blsp_fit_dt)

    # Calculate site annomalies
    blsp_ano <- lapply(unique(blsp_fit_dt$ID), function(id) {
        iddt <- blsp_fit_dt[ID == id, ]
        ano_dt <- iddt[, .(
            ID, spec, x, y, Year,
            SOS_ano = scale(MidGreenup, scale = FALSE),
            SOS_scl = scale(MidGreenup, scale = TRUE),
            SOS_mn = mean(MidGreenup),
            EOS_ano = scale(MidGreendown, scale = FALSE),
            EOS_scl = scale(MidGreendown, scale = TRUE),
            EOS_mn = mean(MidGreendown),
            GSL_ano = scale(GSL, scale = FALSE),
            GSL_scl = scale(GSL, scale = TRUE),
            GSL_mn = mean(GSL)
        )]

        return(ano_dt)
    })
    blsp_ano <- rbindlist(blsp_ano)
    colnames(blsp_ano) <- gsub(".V1", "", colnames(blsp_ano))



    # ~ Load climate data and calculate anomalies ----------------------------

    dm_dt <- fread("data/fia_domspec_daymet.csv")
    dm_dt[, ID := paste0("X", site)]
    dm_dt[, date := as_date(paste0(year, "-01-01")) + yday - 1]
    dm_dt[, mth := month(date)]

    dm_dt <- dm_dt[ID %in% unique(blsp_ano$ID)]


    # ~ Load FIA attributes --------------------------------------------------

    # Merge w/ FIA attributes
    fia_attr_dt <- fread(
        file.path(base$imax_fia_proj_dir, "data/fia_attr_dt_r9.csv")
    )

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
            richness = median(richness)
        )],
        by = ID
    ]

    return(list(
        blsp_ano = blsp_ano, dm_dt = dm_dt, fia_attr_mn_dt = fia_attr_mn_dt
    ))
}



FitBHM <- function(model_str, data, inits, var_names,
    nchain = 3, ifplot = FALSE, quiet = FALSE, progress.bar = "text"
) {
    tstart <- Sys.time()

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
    # Check convergence
    if (ifplot) {
        plot(samp)
    }
    # message("Compute effective size...")
    # eff <- round(effectiveSize(samp), 1)
    eff <- NULL
    # Compute DIC
    message("Compute DIC...")
    DIC <- dic.samples(mod, n.iter = 5000, n.thin = 5, progress.bar = "none")
    # Compute WAIC
    message("Compute WAIC...")
    samp <- Reduce(rbind, samp)
    like <- samp[, grep("^like", colnames(samp))]
    fbar <- colMeans(like)
    Pw <- sum(apply(log(like), 2, var))
    WAIC <- -2 * sum(log(fbar)) + 2 * Pw

    tend <- Sys.time()
    ttake <- tend - tstart
    message("Time taken:", round(ttake, 2), " ", units(ttake))

    return(list(
        eff = eff, DIC = DIC, WAIC = WAIC, Pw = Pw, samp = samp
    ))
}



# ~ Autumm models ####
# ~ ----------------------------------------------------------------------------
# This model only has interannual variables, no site effects
atm_mod_per_site <- "model {
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


# This is essentially the above model but explicitly write out paramters
atmBHM1 <- "model {
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
        beta6[j] ~ dnorm(gamma6[sp[j]], sigma_b6)
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
        gamma6[k] ~ dnorm(mu6, sigma_g6)
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
    sigma_g6 ~ dgamma(0.1, 0.1)
    sigma_g7 ~ dgamma(0.1, 0.1)
    sigma_g8 ~ dgamma(0.1, 0.1)

    mu1 ~ dnorm(0, 0.0001)
    mu2 ~ dnorm(0, 0.0001)
    mu3 ~ dnorm(0, 0.0001)
    mu4 ~ dnorm(0, 0.0001)
    mu5 ~ dnorm(0, 0.0001)
    mu6 ~ dnorm(0, 0.0001)
    mu7 ~ dnorm(0, 0.0001)
    mu8 ~ dnorm(0, 0.0001)

    # WAIC calculations
    for (i in 1:n) {
        like[i] <- dnorm(Y[i], f[i], sigma_y)
    }
}"


# Include site-level predictors
atmBHM3 <- "model {
    # observation-level
    for (i in 1:n) {
        Y[i] ~ dnorm(f[i], sigma_y)
        f[i] <- beta1[site[i]] +
            beta2[site[i]] * atm_Tavg[i] +
            beta3[site[i]] * smr_Prcp[i] +
            beta4[site[i]] * atm_Srad[i] +
            beta5[site[i]] * smr_Tavg[i] +
            beta6[site[i]] * atm_Prcp[i] +
            beta7[site[i]] * smr_Srad[i] +
            beta8[site[i]] * SOS[i]

    }
    # site-level
    for (j in 1:J) {
        beta1[j] ~ dnorm(gamma1[sp[j]], sigma_b1)
        beta2[j] ~ dnorm(gamma2[sp[j]], sigma_b2)

        beta3[j] ~ dnorm(b3[j], sigma_b3)
        b3[j] <- gamma31[sp[j]] +
            gamma32[sp[j]] * nstem[j] +
            gamma33[sp[j]] * smr_mn_Prcp[j]

        beta4[j] ~ dnorm(gamma4[sp[j]], sigma_b4)
        beta5[j] ~ dnorm(gamma5[sp[j]], sigma_b5)
        beta6[j] ~ dnorm(gamma6[sp[j]], sigma_b6)
        beta7[j] ~ dnorm(gamma7[sp[j]], sigma_b7)
        beta8[j] ~ dnorm(gamma8[sp[j]], sigma_b8)
    }
    # species-level
    for (k in 1:K) {
        gamma1[k] ~ dnorm(mu1, sigma_g1)
        gamma2[k] ~ dnorm(mu2, sigma_g2)

        gamma31[k] ~ dnorm(mu31, sigma_g31)
        gamma32[k] ~ dnorm(mu32, sigma_g32)
        gamma33[k] ~ dnorm(mu33, sigma_g33)

        gamma4[k] ~ dnorm(mu4, sigma_g4)
        gamma5[k] ~ dnorm(mu5, sigma_g5)
        gamma6[k] ~ dnorm(mu6, sigma_g6)
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

    sigma_g31 ~ dgamma(0.1, 0.1)
    sigma_g32 ~ dgamma(0.1, 0.1)
    sigma_g33 ~ dgamma(0.1, 0.1)

    sigma_g4 ~ dgamma(0.1, 0.1)
    sigma_g5 ~ dgamma(0.1, 0.1)
    sigma_g6 ~ dgamma(0.1, 0.1)
    sigma_g7 ~ dgamma(0.1, 0.1)
    sigma_g8 ~ dgamma(0.1, 0.1)

    mu1 ~ dnorm(0, 0.0001)
    mu2 ~ dnorm(0, 0.0001)

    mu31 ~ dnorm(0, 0.0001)
    mu32 ~ dnorm(0, 0.0001)
    mu33 ~ dnorm(0, 0.0001)

    mu4 ~ dnorm(0, 0.0001)
    mu5 ~ dnorm(0, 0.0001)
    mu6 ~ dnorm(0, 0.0001)
    mu7 ~ dnorm(0, 0.0001)
    mu8 ~ dnorm(0, 0.0001)

    # WAIC calculations
    for (i in 1:n) {
        like[i] <- dnorm(Y[i], f[i], sigma_y)
    }
}"




