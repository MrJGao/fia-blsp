# ******************************************************************************
# Model BLSP anomalies by biotic and abiotic factors using a Bayesian
# hierarchical model.
# 
# Author: Xiaojie Gao
# Date: 2024-12-19
# ******************************************************************************
rm(list = ls())
source("src/base.R")
source("src/species_phenology/03_pheno_variation_explain/hlp_fit_bayes.R")
library(data.table)
library(lubridate)
library(magrittr)
library(rjags)





pipedir <- "pipe/species_pheno"
dir.create(pipedir, showWarnings = FALSE, recursive = TRUE)



# ~ Spring ####
# ~ ----------------------------------------------------------------------------
spr_dm <- dm_dt[mth %in% c(3, 4, 5), .(
    spr_Tavg = mean((`tmax (deg c)` + `tmin (deg c)`) / 2),
    spr_Prcp = mean(`prcp (mm/day)`), 
    spr_Srad = mean(`srad (W/m^2)`)
), by = c("ID", "year")]

# Anomalies
spr_dm_ano <- lapply(unique(spr_dm$ID), function(id) {
    iddt <- spr_dm[ID == id, ]
    ano_dt <- iddt[, .(
        ID, year,
        spr_Tavg_ano = scale(spr_Tavg, scale = FALSE),
        spr_Tavg_scl = scale(spr_Tavg, scale = TRUE),
        spr_Tavg_mn = mean(spr_Tavg),
        spr_Prcp_ano = scale(spr_Prcp, scale = FALSE),
        spr_Prcp_scl = scale(spr_Prcp, scale = TRUE),
        spr_Prcp_mn = mean(spr_Prcp),
        spr_Srad_ano = scale(spr_Srad, scale = FALSE),
        spr_Srad_scl = scale(spr_Srad, scale = TRUE),
        spr_Srad_mn = mean(spr_Srad)
    )]
    
    return(ano_dt)
})
spr_dm_ano <- rbindlist(spr_dm_ano)


spr_dt <- merge(blsp_ano, spr_dm_ano, by.x = c("ID", "Year"), by.y = c("ID", "year"))
spr_dt <- merge(spr_dt, fia_attr_mn_dt, by = "ID")
spr_dt <- na.omit(spr_dt)

# out:
fwrite(spr_dt, file.path(pipedir, "spr_dt4bayes.csv"))


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



# ~ Spring model 1 ####
# ~ ----------------------------------------------------------------------------
# This model only include `Tavg` as the sole predictor
sprBHM1 <- "model {
    # observation-level
    for (i in 1:n) {
        Y[i] ~ dnorm(f[i], sigma_y)
        f[i] <- beta1[site[i]] +
            beta2[site[i]] * Tavg[i]
    }
    # site-level
    for (j in 1:J) {
        beta1[j] ~ dnorm(gamma1[sp[j]], sigma_b1)
        beta2[j] ~ dnorm(gamma2[sp[j]], sigma_b2)
    }
    # species-level
    for (k in 1:K) {
        gamma1[k] ~ dnorm(mu1, sigma_g1)
        gamma2[k] ~ dnorm(mu2, sigma_g2)
    }

    # priors
    sigma_y ~ dgamma(0.1, 0.1)
    sigma_b1 ~ dgamma(0.1, 0.1)
    sigma_b2 ~ dgamma(0.1, 0.1)
    sigma_g1 ~ dgamma(0.1, 0.1)
    sigma_g2 ~ dgamma(0.1, 0.1)
    
    mu1 ~ dnorm(0, 0.0001)
    mu2 ~ dnorm(0, 0.0001)

    # WAIC calculations
    for (i in 1:n) {
        like[i] <- dnorm(Y[i], f[i], sigma_y)
    }
}"

tstart <- Sys.time()
spr_model1 <- FitBHM(
    sprBHM1,
    data = list(
        Y = Y, n = n, site = site, Tavg = Tavg, 
        J = nsite, sp = sp, 
        K = nsp
    ), 
    inits = list(
        beta1 = runif(nsite, 0, 1),
        beta2 = runif(nsite, 0, 5)
    ), 
    var_names = c(
        "beta1", "beta2", "sigma_b1", "sigma_b2",
        "gamma1", "gamma2", "sigma_g1", "sigma_g2",
        "mu1", "mu2",
        "f", "like"
    ),
    nchain = 3
)
tend <- Sys.time()
ttake <- tend - tstart
message(round(ttake, 2), " ", units(ttake))

pipedir <- "pipe/species_pheno"
dir.create(pipedir, showWarnings = FALSE, recursive = TRUE)

pdf(file.path(pipedir, "spring_BHM_fit.pdf"))
pred <- apply(spr_model1$samp[, grep("^f", colnames(spr_model1$samp))], 2, median)
vis_base$smoothScatter_J2(Y, pred, main = "Spring model 1")



# ~ Spring model 2 ####
# ~ ----------------------------------------------------------------------------
# This model include both `Tavg` and `Prcp`
sprBHM2 <- "model {
    # observation-level
    for (i in 1:n) {
        Y[i] ~ dnorm(f[i], sigma_y)
        f[i] <- beta1[site[i]] +
            beta2[site[i]] * Tavg[i] +
            beta3[site[i]] * Prcp[i]
    }
    # site-level
    for (j in 1:J) {
        beta1[j] ~ dnorm(gamma1[sp[j]], sigma_b1)
        beta2[j] ~ dnorm(gamma2[sp[j]], sigma_b2)
        beta3[j] ~ dnorm(gamma3[sp[j]], sigma_b3)
    }
    # species-level
    for (k in 1:K) {
        gamma1[k] ~ dnorm(mu1, sigma_g1)
        gamma2[k] ~ dnorm(mu2, sigma_g2)
        gamma3[k] ~ dnorm(mu3, sigma_g3)
    }

    # priors
    sigma_y ~ dgamma(0.1, 0.1)
    sigma_b1 ~ dgamma(0.1, 0.1)
    sigma_b2 ~ dgamma(0.1, 0.1)
    sigma_b3 ~ dgamma(0.1, 0.1)
    
    sigma_g1 ~ dgamma(0.1, 0.1)
    sigma_g2 ~ dgamma(0.1, 0.1)
    sigma_g3 ~ dgamma(0.1, 0.1)

    mu1 ~ dnorm(0, 0.0001)
    mu2 ~ dnorm(0, 0.0001)
    mu3 ~ dnorm(0, 0.0001)

    # WAIC calculations
    for (i in 1:n) {
        like[i] <- dnorm(Y[i], f[i], sigma_y)
    }
}"

tstart <- Sys.time()
spr_model2 <- FitBHM(
    sprBHM2,
    data = list(
        Y = Y, n = n, site = site, 
        Tavg = Tavg, Prcp = Prcp,
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
        "gamma1", "sigma_g1",
        "gamma2", "sigma_g2",
        "gamma3", "sigma_g3",
        "mu1", "mu2", "mu3",
        "f", "like"
    ),
    nchain = 3
)
tend <- Sys.time()
ttake <- tend - tstart
message(round(ttake, 2), " ", units(ttake))

pred <- apply(spr_model2$samp[, grep("^f", colnames(spr_model2$samp))], 2, median)
vis_base$smoothScatter_J2(Y, pred, main = "Spring model 2")



# ~ Spring model 3 ####
# ~ ----------------------------------------------------------------------------
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
spr_model3 <- FitBHM(
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
        "gamma1", "sigma_g1",
        "gamma2", "sigma_g2",
        "gamma3", "sigma_g3",
        "mu1", "mu2", "mu3",
        "f", "like"
    ),
    nchain = 3
)
tend <- Sys.time()
ttake <- tend - tstart
message(round(ttake, 2), " ", units(ttake))

pred <- apply(
    spr_model3$samp[, grep("^f", colnames(spr_model3$samp))], 
    2, 
    median
)
vis_base$smoothScatter_J2(Y, pred, main = "Spring model 3")

dev.off()

# DIC comparison
spr_model1$DIC
spr_model2$DIC
spr_model3$DIC
# WAIC comparison
spr_model1$WAIC
spr_model2$WAIC
spr_model3$WAIC

# ^^ Based on the comparison of the previous 3 models:
#   Both DIC and WAIC chose the model w/ all three predictors. R2 is about 0.37
#   with no site-level predictors.



# ~ Spring model 4 ####
# ~ ----------------------------------------------------------------------------
# This model includes site-level predictor `nstem` to model temperature
# sensitivity 
sprBHM4 <- "model {
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
        beta2[j] ~ dnorm(b2[j], sigma_b2)
        b2[j] <- gamma21[sp[j]] + gamma22[sp[j]] * nstem[j]

        beta3[j] ~ dnorm(gamma3[sp[j]], sigma_b3)
        beta4[j] ~ dnorm(gamma4[sp[j]], sigma_b4)
    }
    # species-level
    for (k in 1:K) {
        gamma1[k] ~ dnorm(mu1, sigma_g1)
        gamma21[k] ~ dnorm(mu21, sigma_g21)
        gamma22[k] ~ dnorm(mu22, sigma_g22)
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
    sigma_g21 ~ dgamma(0.1, 0.1)
    sigma_g22 ~ dgamma(0.1, 0.1)
    sigma_g3 ~ dgamma(0.1, 0.1)
    sigma_g4 ~ dgamma(0.1, 0.1)

    mu1 ~ dnorm(0, 0.0001)
    mu21 ~ dnorm(0, 0.0001)
    mu22 ~ dnorm(0, 0.0001)
    mu3 ~ dnorm(0, 0.0001)
    mu4 ~ dnorm(0, 0.0001)

    # WAIC calculations
    for (i in 1:n) {
        like[i] <- dnorm(Y[i], f[i], sigma_y)
    }
}"

tstart <- Sys.time()
spr_model4 <- FitBHM(
    sprBHM4,
    data = list(
        Y = Y, n = n, site = site,
        Tavg = Tavg, Prcp = Prcp, Srad = Srad, 
        J = nsite, sp = sp,
        nstem = nstem,
        K = nsp
    ),
    inits = list(
        beta1 = runif(nsite, 0, 10),
        beta2 = runif(nsite, 0, 10),
        beta3 = runif(nsite, 0, 10),
        gamma21 = runif(nsp, 0, 10),
        gamma22 = runif(nsp, 0, 10)
    ),
    var_names = c(
        "beta1", "sigma_b1",
        "beta2", "sigma_b2",
        "beta3", "sigma_b3",
        "gamma1", "sigma_g1",
        "gamma21", "sigma_g21",
        "gamma22", "sigma_g22",
        "gamma3", "sigma_g3",
        "mu1", "mu21", "mu22", "mu3",
        "f", "like"
    ),
    nchain = 3
)
tend <- Sys.time()
ttake <- tend - tstart
message(round(ttake, 2), " ", units(ttake))

pred <- apply(spr_model4$samp[, grep("^f", colnames(spr_model4$samp))], 2, median)
vis_base$smoothScatter_J2(Y, pred, main = "Spring model 4")



# ~ Spring model 5 ####
# ~ ----------------------------------------------------------------------------
# Starting from here, I fit the full model and then reduce the insignificant
# parameters b/c fitting these models really takes a long time.
sprBHM5 <- "model {
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
        beta1[j] ~ dnorm(b1[j], sigma_b1)
        b1[j] <- gamma11[sp[j]] +
            gamma12[sp[j]] * nstem[j] +
            gamma13[sp[j]] * ba[j] +
            gamma14[sp[j]] * richness[j]

        beta2[j] ~ dnorm(b2[j], sigma_b2)
        b2[j] <- gamma21[sp[j]] + 
            gamma22[sp[j]] * nstem[j] +
            gamma23[sp[j]] * ba[j] +
            gamma24[sp[j]] * richness[j]

        beta3[j] ~ dnorm(b3[j], sigma_b3)
        b3[j] <- gamma31[sp[j]] + 
            gamma32[sp[j]] * nstem[j] +
            gamma33[sp[j]] * ba[j] +
            gamma34[sp[j]] * richness[j]

        beta4[j] ~ dnorm(b4[j], sigma_b4)
        b4[j] <- gamma41[sp[j]] +
            gamma42[sp[j]] * nstem[j] +
            gamma43[sp[j]] * ba[j] +
            gamma44[sp[j]] * richness[j]
    }
    # species-level
    for (k in 1:K) {
        gamma11[k] ~ dnorm(mu11, sigma_g11)
        gamma12[k] ~ dnorm(mu12, sigma_g12)
        gamma13[k] ~ dnorm(mu13, sigma_g13)
        gamma14[k] ~ dnorm(mu14, sigma_g14)
        
        gamma21[k] ~ dnorm(mu21, sigma_g21)
        gamma22[k] ~ dnorm(mu22, sigma_g22)
        gamma23[k] ~ dnorm(mu23, sigma_g23)
        gamma24[k] ~ dnorm(mu24, sigma_g24)
        
        gamma31[k] ~ dnorm(mu31, sigma_g31)
        gamma32[k] ~ dnorm(mu32, sigma_g32)
        gamma33[k] ~ dnorm(mu33, sigma_g33)
        gamma34[k] ~ dnorm(mu34, sigma_g34)
        
        gamma41[k] ~ dnorm(mu41, sigma_g41)
        gamma42[k] ~ dnorm(mu42, sigma_g42)
        gamma43[k] ~ dnorm(mu43, sigma_g43)
        gamma44[k] ~ dnorm(mu44, sigma_g44)
    }

    # priors
    sigma_y ~ dgamma(0.1, 0.1)
    sigma_b1 ~ dgamma(0.1, 0.1)
    sigma_b2 ~ dgamma(0.1, 0.1)
    sigma_b3 ~ dgamma(0.1, 0.1)
    sigma_b4 ~ dgamma(0.1, 0.1)

    sigma_g11 ~ dgamma(0.1, 0.1)
    sigma_g12 ~ dgamma(0.1, 0.1)
    sigma_g13 ~ dgamma(0.1, 0.1)
    sigma_g14 ~ dgamma(0.1, 0.1)
    
    sigma_g21 ~ dgamma(0.1, 0.1)
    sigma_g22 ~ dgamma(0.1, 0.1)
    sigma_g23 ~ dgamma(0.1, 0.1)
    sigma_g24 ~ dgamma(0.1, 0.1)
    
    sigma_g31 ~ dgamma(0.1, 0.1)
    sigma_g32 ~ dgamma(0.1, 0.1)
    sigma_g33 ~ dgamma(0.1, 0.1)
    sigma_g34 ~ dgamma(0.1, 0.1)
    
    sigma_g41 ~ dgamma(0.1, 0.1)
    sigma_g42 ~ dgamma(0.1, 0.1)
    sigma_g43 ~ dgamma(0.1, 0.1)
    sigma_g44 ~ dgamma(0.1, 0.1)

    mu11 ~ dnorm(0, 0.0001)
    mu12 ~ dnorm(0, 0.0001)
    mu13 ~ dnorm(0, 0.0001)
    mu14 ~ dnorm(0, 0.0001)
    
    mu21 ~ dnorm(0, 0.0001)
    mu22 ~ dnorm(0, 0.0001)
    mu23 ~ dnorm(0, 0.0001)
    mu24 ~ dnorm(0, 0.0001)
    
    mu31 ~ dnorm(0, 0.0001)
    mu32 ~ dnorm(0, 0.0001)
    mu33 ~ dnorm(0, 0.0001)
    mu34 ~ dnorm(0, 0.0001)
    
    mu41 ~ dnorm(0, 0.0001)
    mu42 ~ dnorm(0, 0.0001)
    mu43 ~ dnorm(0, 0.0001)
    mu44 ~ dnorm(0, 0.0001)
    
    # WAIC calculations
    for (i in 1:n) {
        like[i] <- dnorm(Y[i], f[i], sigma_y)
    }
}"

tstart <- Sys.time()
spr_model5 <- FitBHM(
    sprBHM5,
    data = list(
        Y = Y, n = n, site = site,
        Tavg = Tavg, Prcp = Prcp, Srad = Srad, 
        J = nsite, sp = sp,
        nstem = nstem, ba = ba, richness = richness,
        K = nsp
    ),
    inits = list(
        beta1 = runif(nsite, 0, 10),
        beta2 = runif(nsite, 0, 10),
        beta3 = runif(nsite, 0, 10),
        gamma21 = runif(nsp, 0, 10),
        gamma22 = runif(nsp, 0, 10)
    ),
    var_names = c(
        "beta1", "sigma_b1",
        "beta2", "sigma_b2",
        "beta3", "sigma_b3",
        "beta4", "sigma_b4",
        
        "gamma11", "sigma_g11",
        "gamma12", "sigma_g12",
        "gamma13", "sigma_g13",
        "gamma14", "sigma_g14",
        
        "gamma21", "sigma_g21",
        "gamma22", "sigma_g22",
        "gamma23", "sigma_g23",
        "gamma24", "sigma_g24",
        
        "gamma31", "sigma_g31",
        "gamma32", "sigma_g32",
        "gamma33", "sigma_g33",
        "gamma34", "sigma_g34",
        
        "gamma41", "sigma_g41",
        "gamma42", "sigma_g42",
        "gamma43", "sigma_g43",
        "gamma44", "sigma_g44",
        
        "mu11", "mu12", "mu13", "mu14",
        "mu21", "mu22", "mu23", "mu24",
        "mu31", "mu32", "mu33", "mu34",
        "mu41", "mu42", "mu43", "mu44",
        
        "f", "like"
    ),
    nchain = 3
)
tend <- Sys.time()
ttake <- tend - tstart
message(round(ttake, 2), " ", units(ttake))

pred <- apply(spr_model5$samp[, grep("^f", colnames(spr_model5$samp))], 2, median)
vis_base$smoothScatter_J2(Y, pred, main = "Spring model 5")

spr_model4$DIC
spr_model5$DIC

# ^^ So, for the spring model, group level predictors do not improve the fit. 
