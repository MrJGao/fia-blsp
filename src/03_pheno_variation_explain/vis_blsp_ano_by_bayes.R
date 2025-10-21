# ******************************************************************************
# Visualize coefficients of the Bayesian hierarchical models for spring and
# autumn phenology, respectively
# 
# Author: Xiaojie Gao
# Date: 2025-01-04
# ******************************************************************************
rm(list = ls())
source("src/base.R")
library(data.table)
library(magrittr)



pipedir <- "pipe/species_pheno"
outdir <- "out/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)



# ~ Load spring model ####
# ~ ----------------------------------------------------------------------------
spr_mod_file <- file.path(pipedir, "spr_model3.Rds")
if (!file.exists(spr_mod_file)) {
    source("src/mod_blsp_ano_by_bayes2.R")
}
spr_mod <- readRDS(spr_mod_file)
spr_dt <- fread(file.path(pipedir, "spr_dt4bayes.csv"))



# ~ Load autumn model ####
# ~ ----------------------------------------------------------------------------
# atm_mod_file <- file.path(pipedir, "atm_model4.Rds")
atm_mod_file <- file.path(pipedir, "atm_model_1a.Rds")
if (!file.exists(atm_mod_file)) {
    source("src/mod_blsp_ano_by_bayes2.R")
}

atm_mod <- readRDS(atm_mod_file)
atm_dt <- fread(file.path(pipedir, "atm_dt4bayes.csv"))



# ~ Goodness-of-fit ####
# ~ ----------------------------------------------------------------------------
{
    png(
        file.path(outdir, "bayes_goodness_of_fit.png"),
        width = 2000, height = 1000, res = 300
    )
    par(
        mfrow = c(1, 2), mar = c(3, 3, 2, 1), mgp = c(1.5, 0.5, 0),
        bty = "L", cex.axis = 1.2, cex.lab = 1.2
    )

    # Plot goodness-of-fit
    spr_pred <- apply(spr_mod$samp[, grep("^f", colnames(spr_mod$samp))], 2, median)
    smoothScatter_J2(
        spr_pred, spr_dt$SOS_scl.V1,
        colorful = TRUE,
        main = "Spring model fit",
        xlab = "Pred.", ylab = "Obs."
    )

    atm_pred <- apply(atm_mod$samp[, grep("^f", colnames(atm_mod$samp))], 2, median)
    smoothScatter_J2(
        atm_pred, atm_dt$EOS_scl.V1,
        colorful = TRUE,
        main = "Autumn model fit",
        xlab = "Pred.", ylab = "Obs."
    )

    dev.off()
}


# ~ Spring ####
# ~ ----------------------------------------------------------------------------
mu2 <- quantile(
    spr_mod$samp[, grep("^mu2", colnames(spr_mod$samp))], 
    c(0.025, 0.5, 0.975)
)
mu3 <- quantile(
    spr_mod$samp[, grep("^mu3", colnames(spr_mod$samp))],
    c(0.025, 0.5, 0.975)
)
mu4 <- quantile(
    spr_mod$samp[, grep("^mu4", colnames(spr_mod$samp))],
    c(0.025, 0.5, 0.975)
)

gamma2 <- apply(
    spr_mod$samp[, grep("^gamma2", colnames(spr_mod$samp))],
    2,
    quantile, c(0.025, 0.5, 0.975)
)
gamma3 <- apply(
    spr_mod$samp[, grep("^gamma3", colnames(spr_mod$samp))],
    2,
    quantile, c(0.025, 0.5, 0.975)
)
gamma4 <- apply(
    spr_mod$samp[, grep("^gamma4", colnames(spr_mod$samp))],
    2,
    quantile, c(0.025, 0.5, 0.975)
)


# Get species names
latin_names <- fread("data/common10.csv")
latin_names[, latin := paste(GENUS, SPECIES_NAME)]
spr_dt <- merge(spr_dt, latin_names, by.x = "spec", by.y = "gesp")

spnames <- as.factor(unique(spr_dt[, .(ID, spec)])$spec) %>%
    levels()
spnames_latin <- latin_names[gesp %in% spnames, latin]

# reorder to reflect genus
# cols <- RColorBrewer::brewer.pal(10, "Paired")
cols <- c(
    "#A6CEE3", "#1F78B4", 
    "#B2DF8A", 
    "#33A02C", 
    "#FDBF6F",
    "#ff787b", "#E31A1C", "#960808", "#b70e6e", 
    "#a25feb"
)

{
    svglite::svglite(file.path(outdir, "spr_ano_BHM.svg"), width = 6, height = 4)
    par(mar = c(2, 3, 2, 1), mgp = c(1.5, 0.5, 0), cex.lab = 1.2)
    plot(
        NA,
        bty = "L", xlim = c(0.5, 3.5), ylim = c(-0.8, 0.2),
        xlab = "", ylab = "Normalized Coefficients",
        xaxt = "n"
    )
    axis(1, at = 1:3, labels = c("Temp", "Prcp", "Srad"))
    abline(h = 0, lty = 2)

    # Species-level
    xjitter <- jitter(rep(1.2, ncol(gamma2)), factor = 5)
    points(xjitter, gamma2[2, ], pch = 16, col = cols)
    segments(xjitter, gamma2[1, ], xjitter, gamma2[3, ], col = cols)

    xjitter <- jitter(rep(2.2, ncol(gamma2)), factor = 4)
    points(xjitter, gamma3[2, ], pch = 16, col = cols)
    segments(xjitter, gamma3[1, ], xjitter, gamma3[3, ], col = cols)

    xjitter <- jitter(rep(3.2, ncol(gamma2)), factor = 3)
    points(xjitter, gamma4[2, ], pch = 16, col = cols)
    segments(xjitter, gamma4[1, ], xjitter, gamma4[3, ], col = cols)

    # Overall
    points(1:3, c(mu2[2], mu3[2], mu4[2]), pch = 16, cex = 1.5)
    segments(
        1:3, c(mu2[1], mu3[1], mu4[1]), 1:3, c(mu2[3], mu3[3], mu4[3])
    )

    names <- lapply(spnames_latin, function(x) {
        bquote(italic(.(x)))
    })

    legend(
        grconvertX(0.45, "ndc"), grconvertY(0.5, "ndc"),
        xjust = 0,
        bty = "n", ncol = 2,
        legend = do.call(expression, names),
        pch = 16, col = cols,
        xpd = NA,
        cex = 0.9
    )

    dev.off()
}



# ~ Autumn ####
# ~ ----------------------------------------------------------------------------
# ! This section may only works for atm model 1a

PlotSpecEff <- function(mod, idx, at, fact) {
    gamma <- apply(
        mod$samp[, grep(paste0("^gamma\\[", idx, ",\\d+\\]"), colnames(mod$samp))],
        2,
        quantile, c(0.025, 0.5, 0.975)
    )
    xjitter <- jitter(rep(at + 0.2, ncol(gamma)), factor = fact)
    segments(xjitter, gamma[1, ], xjitter, gamma[3, ], col = cols)
    points(xjitter, gamma[2, ], pch = 16, col = cols)
}

PlotMnEff <- function(mod, idx, at) {
    mu <- quantile(
        mod$samp[, grep(paste0("^mu\\[", idx, "\\]$"), colnames(mod$samp))],
        c(0.025, 0.5, 0.975)
    )
    segments(at, mu[1], at, mu[3])
    points(at, mu[2], pch = 16, cex = 1.5)
}


sp <- as.numeric(as.factor(unique(atm_dt[, .(ID, spec)])$spec))
nsp <- length(unique(sp))

{
    svglite::svglite(file.path(outdir, "atm_ano_BHM.svg"), width = 6, height = 5)
    par(mar = c(6, 3, 2, 1), mgp = c(1.5, 0.5, 0), cex.lab = 1.2)
    plot(
        NA,
        bty = "L", xlim = c(1.5, 9.5), ylim = c(-0.3, 0.5),
        xlab = "", ylab = "Normalized Coefficients",
        xaxt = "n"
    )
    axis(
        1, at = 2:8, 
        labels = c(
            "SOS",
            "atm_Temp", "atm_Prcp", "atm_Srad",
            "smr_Temp", "smr_Prcp", "smr_Srad"
        ),
        las = 2
    )
    abline(h = 0, lty = 2)

    for (i in c(2:8)) {
        PlotSpecEff(atm_mod, idx = i, at = i, fact = 2.5)
        PlotMnEff(atm_mod, idx = i, at = i)
    }
    
    dev.off()
}
