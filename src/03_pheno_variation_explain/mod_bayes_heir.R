# ******************************************************************************
# Build a Bayesian heirarchical model to quantify species-specific environmental
# effects on phenology.
# 
# Author: Xiaojie Gao
# Date: 2024-04-06
# ******************************************************************************
library(data.table)
library(rjags)
library(magrittr)



PlotSpecFit <- function(i) {
    spname <- levels(sp_fac)[i]

    spec_dt <- dt[spec == spname, ]
    b <- c(beta1[i, med], beta2[i, med], beta3[i, med], beta4[i, med], beta5[i, med])
    x <- as.matrix(spec_dt[, .(1, lat, (tmax + tmin) / 2, prcp, swe)])
    pred <- x %*% b

    # RMSE
    rmse <- round(sqrt(mean((pred - spec_dt$avg50PCGI)^2)), 2)
    rmse_str <- bquote(RMSE == .(rmse))

    plot(pred, spec_dt$avg50PCGI,
        xlim = range(pred, spec_dt$avg50PCGI), xlab = "Pred.",
        ylim = range(pred, spec_dt$avg50PCGI), ylab = "Obs.",
        pch = 16,
        main = spname
    )
    mtext(side = 3, spec_table[spec == spname, gsub("\r\n", " ", CommonName)])
    abline(1, 1, lty = 2)

    legend("bottomright", legend = rmse_str, bty = "n", y.intersp = 3)
}


PlotSpecEff <- function(b, e, ...) {
    plot(NA,
        xlim = c(0, nsp),
        ylim = range(b),
        xlab = "", ylab = "Effects",
        xaxt = "n",
        bty = "L",
        ...
    )
    abline(h = 0, lty = 2)
    axis(side = 1, at = 1:nsp, labels = levels(sp_fac), las = 2, gap.axis = 0.01)

    # Beta
    points(1:nsp, b[, med], pch = 16, cex = 1.5)
    segments(
        x0 = 1:nsp, y0 = b[, lwr],
        x1 = 1:nsp, y1 = b[, upr]
    )
    # Eta
    points(0, e[, med], pch = 16, cex = 2, col = "blue")
    segments(
        x0 = 0, y0 = e[, lwr],
        x1 = 0, y1 = e[, upr],
        col = "blue"
    )
}



mod_jags <- "model {
    # likelihood
    for (i in 1:n) {
        Y[i] ~ dnorm(f[i], tau_y)
        f[i] <- beta1[sp[i]] + 
            beta2[sp[i]] * lat[i] +
            beta3[sp[i]] * tavg[i] + 
            beta4[sp[i]] * prcp[i] + 
            beta5[sp[i]] * swe[i] +
            beta6[sp[i]] * dem[i]
    }

    # priors
    # nsp - number of species
    for (i in 1:nsp) {
        beta1[i] ~ dnorm(eta[1], tau[1])
        beta2[i] ~ dnorm(eta[2], tau[2])
        beta3[i] ~ dnorm(eta[3], tau[3])
        beta4[i] ~ dnorm(eta[4], tau[4])
        beta5[i] ~ dnorm(eta[5], tau[5])
        beta6[i] ~ dnorm(eta[6], tau[6])
    }

    tau_y ~ dgamma(0.1, 0.1)

    for (i in 1:6) {
        eta[i] ~ dnorm(0, 0.0001)
        tau[i] ~ dgamma(0.1, 0.1)
    }

    sigma_y <- 1 / sqrt(tau_y)
    sigma_beta <- 1 / sqrt(tau)
}"



# ~ Format data for JAGS ----------------------------------------
dt <- fread("data/sp_60_mslsp_dm.csv")
dt <- dt[Type == "deci"]

# Make a avg through 2016-2019
dt <- dt[, lapply(.SD, mean), 
    by = c("X", "spec", "Genus", "Species", "Type")
]
spec_table <- fread("data/fia_domspec_truexy.csv") %>%
    .[, .(spec, CommonName, Genus)] %>%
    unique()

sp_fac <- as.factor(dt$spec)
nsp <- length(levels(sp_fac))
gn_fac <- as.factor(dt$Genus)
ngn <- length(gn_fac)

dt[, .N, by = list(spec, Genus)] %>%
    setorder(Genus) %>%
    print()

data <- list(
    Y = dt$avg50PCGI,
    sp = as.numeric(sp_fac), nsp = nsp,
    lat = dt$lat, tavg = (dt$tmax + dt$tmin) / 2, 
    prcp = dt$prcp, swe = dt$swe,
    dem = dt$DEM,
    n = nrow(dt)
)

inits <- list(
    beta1 = runif(nsp, 0, 30), beta2 = runif(nsp, 0, 20),
    beta3 = runif(nsp, 0, 30), beta4 = runif(nsp, 0, 30),
    beta5 = runif(nsp, 0, 30)
)

model <- jags.model(textConnection(mod_jags),
    data = data, inits = inits,
    n.chains = 2,
    quiet = TRUE
)
update(model, 20000, progress.bar = "none")

iter <- 1
repeat({
    samp <- coda.samples(model,
        variable.names = c(
            "beta1", "beta2", "beta3", "beta4", "beta5", "beta6",
            "eta",
            "sigma_y", "sigma_beta"
        ),
        n.iter = 50000, thin = 10,
        progress.bar = "none"
    )
    gelman_diag <- gelman.diag(samp)
    mpsrf <- gelman_diag$mpsrf
    iter <- iter + 1
    if (mpsrf < 1.1) {
        break
    }

    if (iter > 6) {
        break
        warning("May not converged!")
    }
})


# ~ Summarize model results ----------------------------------------
SummarySample <- function(samp) {
    if (is.null(dim(samp))) {
        samp_dt <- quantile(samp, c(0.025, 0.5, 0.975))
    } else {
        samp_dt <- apply(samp, 2, quantile, c(0.025, 0.5, 0.975))
    }
    samp_dt <- data.table(t(samp_dt))
    colnames(samp_dt) <- c("lwr", "med", "upr")
    return(samp_dt)
}

samp_df <- as.data.frame(rbind(samp[[1]], samp[[2]]))
samp_names <- colnames(samp_df)

beta1 <- SummarySample(samp_df[, grep("^beta1", samp_names)])
beta2 <- SummarySample(samp_df[, grep("^beta2", samp_names)])
beta3 <- SummarySample(samp_df[, grep("^beta3", samp_names)])
beta4 <- SummarySample(samp_df[, grep("^beta4", samp_names)])
beta5 <- SummarySample(samp_df[, grep("^beta5", samp_names)])
beta6 <- SummarySample(samp_df[, grep("^beta6", samp_names)])

sigma_y <- SummarySample(samp_df[, grep("^sigma_y", samp_names)])
sigma_beta <- SummarySample(samp_df[, grep("^sigma_beta", samp_names)])

eta <- SummarySample(samp_df[, grep("^eta", samp_names)])



pdf(paste0("pipe/domspec_bayes_avg.pdf"), width = 10, height = 5)
par(mfrow = c(2, 5), mgp = c(1.5, 0.5, 0), mar = c(3, 3, 3, 2))
for (i in 1:nsp) {
    PlotSpecFit(i)
}

par(mfrow = c(1, 1))
PlotSpecEff(beta2, eta[2, ], main = "Lat")
PlotSpecEff(beta3, eta[3, ], main = "Avg temperature")
PlotSpecEff(beta4, eta[4, ], main = "Precipitation")
PlotSpecEff(beta5, eta[5, ], main = "Snow Water Equivalent")
PlotSpecEff(beta6, eta[6, ], main = "DEM")
dev.off()




# ~ For each year ####
# ~ ----------------------------------------------------------------------------

for (yr in 2016:2019) {
    dt <- fread("data/sp_60_mslsp_dm.csv")
    spec_table <- fread("data/fia_domspec_truexy.csv") %>%
        .[, .(spec, CommonName, Genus)] %>%
        unique()

    dt <- dt[Year == yr & Type == "deci"]

    sp_fac <- as.factor(dt$spec)
    nsp <- length(levels(sp_fac))
    gn_fac <- as.factor(dt$Genus)
    ngn <- length(gn_fac)

    dt[, .N, by = list(spec, Genus)] %>%
        setorder(Genus) %>%
        print()

    data <- list(
        Y = dt$avg50PCGI,
        sp = as.numeric(sp_fac), nsp = nsp,
        lat = dt$lat, tavg = (dt$tmax + dt$tmin) / 2, 
        prcp = dt$prcp, swe = dt$swe,
        n = nrow(dt)
    )

    inits <- list(
        beta1 = runif(nsp, 0, 30), beta2 = runif(nsp, 0, 20),
        beta3 = runif(nsp, 0, 30), beta4 = runif(nsp, 0, 30),
        beta5 = runif(nsp, 0, 30)
    )

    model <- jags.model(textConnection(mod_jags),
        data = data, inits = inits,
        n.chains = 2,
        quiet = TRUE
    )
    update(model, 20000, progress.bar = "none")

    iter <- 1
    repeat({
        samp <- coda.samples(model,
            variable.names = c(
                "beta1", "beta2", "beta3", "beta4", "beta5",
                "eta",
                "sigma_y", "sigma_beta"
            ),
            n.iter = 50000, thin = 10,
            progress.bar = "none"
        )
        gelman_diag <- gelman.diag(samp)
        mpsrf <- gelman_diag$mpsrf
        iter <- iter + 1
        if (mpsrf < 1.1) {
            break
        }

        if (iter > 6) {
            break
            warning("May not converged!")
        }
    })


    # ~ Summarize model results ----------------------------------------
    SummarySample <- function(samp) {
        if (is.null(dim(samp))) {
            samp_dt <- quantile(samp, c(0.025, 0.5, 0.975))
        } else {
            samp_dt <- apply(samp, 2, quantile, c(0.025, 0.5, 0.975))
        }
        samp_dt <- data.table(t(samp_dt))
        colnames(samp_dt) <- c("lwr", "med", "upr")
        return(samp_dt)
    }

    samp_df <- as.data.frame(rbind(samp[[1]], samp[[2]]))
    samp_names <- colnames(samp_df)

    beta1 <- SummarySample(samp_df[, grep("^beta1", samp_names)])
    beta2 <- SummarySample(samp_df[, grep("^beta2", samp_names)])
    beta3 <- SummarySample(samp_df[, grep("^beta3", samp_names)])
    beta4 <- SummarySample(samp_df[, grep("^beta4", samp_names)])
    beta5 <- SummarySample(samp_df[, grep("^beta5", samp_names)])

    sigma_y <- SummarySample(samp_df[, grep("^sigma_y", samp_names)])
    sigma_beta <- SummarySample(samp_df[, grep("^sigma_beta", samp_names)])

    eta <- SummarySample(samp_df[, grep("^eta", samp_names)])




    pdf(paste0("pipe/domspec_bayes_", yr, ".pdf"), width = 10, height = 5)
    par(mfrow = c(2, 5), mgp = c(1.5, 0.5, 0), mar = c(3, 3, 3, 2))
    for (i in 1:nsp) {
        PlotSpecFit(i)
    }

    par(mfrow = c(1, 1))
    PlotSpecEff(beta2, eta[2, ], main = "Latitude")
    PlotSpecEff(beta3, eta[3, ], main = "Avg temperature")
    PlotSpecEff(beta4, eta[4, ], main = "Precipitation")
    PlotSpecEff(beta5, eta[5, ], main = "Snow Water Equivalent")
    dev.off()
}

