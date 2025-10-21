# ******************************************************************************
# BLSP variation among and within species, sites, and years
# 
# Author: Xiaojie Gao
# Date: 2024-09-03
# ******************************************************************************
rm(list = ls())

source("src/base.R")
library(data.table)
library(terra)
library(magrittr)
library(ggplot2)



tmpfigdir <- "pipe/mod_blsp_variation"
dir.create(tmpfigdir, showWarnings = FALSE)


# Load BLSP
blsp_fit_dt <- fread("data/fia_blsp_fit.csv")[, .(
    ID, spec, x, y, Year,
    MidGreenup,
    MidGreenup_unc = MidGreenup_upr - MidGreenup_lwr,
    MidGreendown, MidGreendown_unc = MidGreendown_upr - MidGreendown_lwr
)]

# Remove records with high uncertainty
blsp_fit_dt[MidGreenup_unc > 30, MidGreenup := NA]
blsp_fit_dt[MidGreendown_unc > 30, MidGreendown := NA]

# Calculate GSL
blsp_fit_dt[, GSL := MidGreendown - MidGreenup]

blsp_fit_dt <- na.omit(blsp_fit_dt)



pdf(file.path(tmpfigdir, "1.pdf"), width = 10, height = 10)

specs <- unique(blsp_fit_dt$spec)
cols <- hcl.colors(10, "Spectral")

# SOS
plot(NA, xlim = range(blsp_fit_dt$Year), ylim = c(100, 200), 
    xlab = "Year", ylab = "SOS"
)
for (i in seq_along(specs)) {
    spec_dt <- blsp_fit_dt[spec == specs[i],]
    spec_var_dt <- spec_dt[, .(
        SOS = mean(MidGreenup), sd = sd(MidGreenup)
    ), by = "Year"]
    
    setorder(spec_var_dt, Year)
    lines(spec_var_dt[, .(Year, SOS)], pch = 16, type = "o", col = cols[i])
    # segments(spec_var_dt$Year, spec_var_dt[, SOS + sd], 
    #     spec_var_dt$Year, spec_var_dt[, SOS - sd], 
    #     col = cols[i]
    # )
}

# EOS
plot(NA,
    xlim = range(blsp_fit_dt$Year), ylim = c(220, 300),
    xlab = "Year", ylab = "EOS"
)
for (i in seq_along(specs)) {
    spec_dt <- blsp_fit_dt[spec == specs[i], ]
    spec_var_dt <- spec_dt[, .(
        EOS = mean(MidGreendown), sd = sd(MidGreendown)
    ), by = "Year"]

    setorder(spec_var_dt, Year)
    lines(spec_var_dt[, .(Year, EOS)], pch = 16, type = "o", col = cols[i])
    # segments(spec_var_dt$Year, spec_var_dt[, EOS + sd],
    #     spec_var_dt$Year, spec_var_dt[, EOS - sd],
    #     col = cols[i]
    # )
}

# GSL
plot(NA,
    xlim = range(blsp_fit_dt$Year), ylim = c(80, 200),
    xlab = "Year", ylab = "GSL"
)
for (i in seq_along(specs)) {
    spec_dt <- blsp_fit_dt[spec == specs[i], ]
    spec_var_dt <- spec_dt[, .(
        GSL = mean(GSL), sd = sd(GSL)
    ), by = "Year"]

    setorder(spec_var_dt, Year)
    lines(spec_var_dt[, .(Year, GSL)], pch = 16, type = "o", col = cols[i])
    # segments(spec_var_dt$Year, spec_var_dt[, GSL + sd],
    #     spec_var_dt$Year, spec_var_dt[, GSL - sd],
    #     col = cols[i]
    # )
}

dev.off()


# ~ Vary among species ----------------------------------------
# Mean phenology of each site, varies among species
mn_pheno <- blsp_fit_dt[, 
    .(
        SOS = mean(MidGreenup), sd_SOS = sd(MidGreenup),
        EOS = mean(MidGreendown), sd_EOS = sd(MidGreendown),
        GSL = mean(GSL), sd_GSL = sd(GSL)
    ), 
    by = c("ID", "spec")
]

pdf(file.path(tmpfigdir, "2.pdf"), width = 10, height = 10)

par(mfrow = c(3, 1))
boxplot(SOS ~ spec, data = mn_pheno)
boxplot(EOS ~ spec, data = mn_pheno)
boxplot(GSL ~ spec, data = mn_pheno)

par(mfrow = c(3, 1))
boxplot(sd_SOS ~ spec, data = mn_pheno)
boxplot(sd_EOS ~ spec, data = mn_pheno)
boxplot(sd_GSL ~ spec, data = mn_pheno)

dev.off()

# So, spatial variation is larger than interannual variation for all species.