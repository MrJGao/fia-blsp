# ******************************************************************************
# Visualize BLSP trends for the FIA dominant species sites.
# 
# Author: Xiaojie Gao
# Date: 2024-10-02
# ******************************************************************************
rm(list=ls())
source("src/base.R")
library(data.table)
library(magrittr)
library(ggridges)
library(ggplot2)



VisTrendSlp <- function(site_trend_dt, figname) {
    p1 <- ggplot() +
        geom_density_ridges(
            data = site_trend_dt[atm_pval < 0.05, ],
            aes(x = atm_slp, y = latin),
            fill = adjustcolor("orange", 0.5),
            color = NA,
            scale = 0.8
        ) +
        geom_density_ridges(
            data = site_trend_dt[spr_pval < 0.05, ],
            aes(x = spr_slp, y = latin),
            fill = adjustcolor("seagreen", 0.5),
            color = NA,
            scale = 0.8
        ) +
        geom_density_ridges(
            data = site_trend_dt[gsl_pval < 0.05, ],
            aes(x = gsl_slp, y = latin),
            fill = NA,
            scale = 0.8,
            color = "blue"
        ) +
        theme_ridges(center_axis_labels = TRUE, font_size = 20) +
        xlim(-1, 1) +
        ylab("") +
        xlab("Slope") +
        geom_vline(xintercept = 0, linetype = "dashed") +
        theme(
            legend.position = "right",
            axis.text.y = element_text(face = "italic")
        )

    ggsave(
        filename = figname, p1,
        width = 5, height = 8
    )
}


VisSigProportion <- function(site_trend_dt, figname) {
    total <- site_trend_dt[, .(TotN = .N),
        keyby = c("spec", "CommonName")
    ]

    sig_dt <- cbind(
        site_trend_dt[spr_pval < 0.05, .(Nspr = .N), keyby = "spec"],
        site_trend_dt[atm_pval < 0.05, .N, keyby = "spec"][, .(Natm = N)],
        site_trend_dt[gsl_pval < 0.05, .N, keyby = "spec"][, .(Ngsl = N)]
    )
    sig_dt <- merge(sig_dt, total, by = "spec")
    sig_dt[, ":="(
        spr_pct = Nspr / TotN * 100,
        atm_pct = Natm / TotN * 100,
        gsl_pct = Ngsl / TotN * 100)]


    AddBars <- function(at, width, height, border = NA, col = "grey") {
        for (i in seq(at)) {
            x1 <- 0
            x2 <- height[i]
            y1 <- at[i] - (width / 2)
            y2 <- at[i] + (width / 2)
            rect(x1, y1, x2, y2, border = border, col = col)
        }
    }


    png(figname,
        width = 1300, height = 1800, res = 300
    )
    par(mar = c(4, 5, 1, 1), cex.lab = 1.2, cex.axis = 1.2, mgp = c(2, 0.7, 0))
    plot(NA,
        xlim = c(0, 35), ylim = c(0.5, 10), bty = "L",
        yaxt = "n", ylab = "", xlab = "Sig. proportion (%)"
    )
    axis(side = 2, at = 1:10, labels = sig_dt$spec, las = 2)

    # Spring
    AddBars(
        at = c(1:10) + 0.2, width = 0.15,
        height = sig_dt[, spr_pct],
        col = "seagreen"
    )
    # Autumn
    AddBars(
        at = c(1:10), width = 0.15,
        height = sig_dt[, atm_pct],
        col = "orange"
    )
    # GSL
    AddBars(
        at = c(1:10) - 0.2, width = 0.15,
        height = sig_dt[, gsl_pct],
        col = "blue"
    )
    dev.off()
}



# ~ Load data ####
# ~ ----------------------------------------------------------------------------

# ~ OLS ----------------------------------------
# The input file was produced by `src/species_phenology/mod_fia_blsp_trend_ols.R`
site_trend_dt_ols <- fread("pipe/fia_dom_blsp_trend_ols.csv")
# Remove evergreen
site_trend_dt_ols <- site_trend_dt_ols[Type == "deci"]

# Get latin names
latin_names <- fread("data/common10.csv")
site_trend_dt_ols <- merge(site_trend_dt_ols, latin_names, by.x = "spec", by.y = "gesp")
site_trend_dt_ols[, latin := paste(GENUS, SPECIES_NAME)]


# ~ MK ----------------------------------------
# The input file was produced by `src/species_phenology/mod_fia_blsp_trend_mk.R`
site_trend_dt <- fread("pipe/fia_dom_blsp_trend.csv")
# Remove evergreen
site_trend_dt <- site_trend_dt[Type == "deci"]

# Get latin names
latin_names <- fread("data/common10.csv")
site_trend_dt <- merge(site_trend_dt, latin_names, by.x = "spec", by.y = "gesp")
site_trend_dt[, latin := paste(GENUS, SPECIES_NAME)]


# ~ Robust SE ----------------------------------------
# The input file was produced by `src/species_phenology/mod_fia_blsp_trend_robust.R`
site_trend_dt_ols_robust <- fread("pipe/fia_dom_blsp_trend_ols_robust.csv")
# Remove evergreen
site_trend_dt_ols_robust <- site_trend_dt_ols_robust[Type == "deci"]

# Get latin names
latin_names <- fread("data/common10.csv")
site_trend_dt_ols_robust <- merge(site_trend_dt_ols_robust, latin_names, by.x = "spec", by.y = "gesp")
site_trend_dt_ols_robust[, latin := paste(GENUS, SPECIES_NAME)]




figoutdir <- "out/"
dir.create(figoutdir, showWarnings = FALSE, recursive = TRUE)



# ~ Trend slopes ####
# ~ ----------------------------------------------------------------------------
# density(site_trend_dt[spec == "ACRU", spr_slp])
# hist(site_trend_dt[spec == "ACRU", spr_slp], breaks = 300)

# OLS
VisTrendSlp(
    site_trend_dt = site_trend_dt_ols,
    figname = file.path(figoutdir, "fia_blsp_trend_ridge_ols.png")
)
# Mann-kendall
VisTrendSlp(
    site_trend_dt = site_trend_dt,
    figname = file.path(figoutdir, "fia_blsp_trend_ridge.png")
)



# ~ Sig. proportions ####
# ~ ----------------------------------------------------------------------------

# OLS
VisSigProportion(
    site_trend_dt = site_trend_dt_ols,
    figname = file.path(figoutdir, "fia_blsp_trend_sigpct_ols.png")
)

# MK
VisSigProportion(
    site_trend_dt = site_trend_dt_ols,
    figname = file.path(figoutdir, "fia_blsp_trend_sigpct_mk.png")
)

# OLS robust SE
VisSigProportion(
    site_trend_dt = site_trend_dt_ols,
    figname = file.path(figoutdir, "fia_blsp_trend_sigpct_ols_robust.png")
)



# ~ Absolute correlations ####
# ~ ----------------------------------------------------------------------------
spr_cor_dt <- site_trend_dt[, .(corr = abs(cor(gsl_slp, spr_slp, use = "complete.obs"))), keyby = "spec"]
atm_cor_dt <- site_trend_dt[, .(corr = abs(cor(gsl_slp, atm_slp, use = "complete.obs"))), keyby = "spec"]

png(
    file.path(figoutdir, "fia_blsp_trend_corr_ols.png"), 
    width = 1200, height = 1800, res = 300
)

par(mar = c(4, 5, 1, 1), cex.lab = 1.2, cex.axis = 1.2, mgp = c(2, 0.7, 0))
plot(NA,
    xlim = c(0, 1), ylim = c(0.5, 10), bty = "L",
    yaxt = "n", ylab = "", xlab = "Correlation with GSL"
)
axis(side = 2, at = 1:10, labels = sig_dt$spec, las = 2)

# Spring
AddBars(
    at = c(1:10) + 0.1, width = 0.15,
    height = spr_cor_dt[, corr],
    col = "seagreen"
)
# Autumn
AddBars(
    at = c(1:10) - 0.1, width = 0.15,
    height = atm_cor_dt[, corr],
    col = "orange"
)

dev.off()


