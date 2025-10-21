# ******************************************************************************
# Species-specific phenology trends using HB ground data.
# 
# Author: Xiaojie Gao
# Date: 2024-07-09
# ******************************************************************************
rm(list = ls())
source("src/base.R")
library(data.table)
library(lubridate)
library(trend)
library(ggridges)
library(ggplot2)



DoTrend_backup <- function(dt, season = "spring") {
    specnames <- unique(dt$SPECIES)

    trend_dt <- lapply(specnames, function(spec) {
        spec_dt <- dt[SPECIES == spec]

        ind_trend_dt <- lapply(unique(spec_dt$SITE), function(ind) {
            ind_dt <- spec_dt[SITE == ind]
            
            # Fit a linear regression
            fit <- tryCatch(
                {
                    if (season == "spring") {
                        plot(ind_dt[, .(YEAR, spring)], 
                            main = ind_dt$SPECIES[1],
                            type = "o"
                        )
                        vis_base$LinXY(ind_dt$YEAR, ind_dt$spring, add = TRUE)
                        lm(spring ~ YEAR, data = ind_dt)
                    } else if (season == "fall") {
                        plot(ind_dt[, .(YEAR, fall)], 
                            main = ind_dt$SPECIES[1],
                            type = "o"
                        )
                        vis_base$LinXY(ind_dt$YEAR, ind_dt$fall, add = TRUE)
                        lm(fall ~ YEAR, data = ind_dt)
                    } else {
                        plot(ind_dt[, .(YEAR, gsl)],
                            main = ind_dt$SPECIES[1],
                            type = "o"
                        )
                        vis_base$LinXY(ind_dt$YEAR, ind_dt$gsl, add = TRUE)
                        lm(gsl ~ YEAR, data = ind_dt)
                    }
                },
                error = function(e) {
                    return(NULL)
                }
            )

            if (is.null(fit) || is.na(coef(fit)[2])) {
                return(NULL)
            }

            slp <- round(coef(fit)[2], 3)
            f <- summary(fit)$fstatistic
            p_val <- pf(f[1], f[2], f[3], lower.tail = FALSE)

            arow <- data.table(
                spec = spec, slp = slp, pval = p_val
            )

            return(arow)
        })
        ind_trend_dt <- do.call(rbind, ind_trend_dt)
    })
    trend_dt <- do.call(rbind, trend_dt)

    return(trend_dt)
}


# Trend analysis in Mann-Kendall
DoTrendMk <- function(dt, season = "spring") {
    specnames <- unique(dt$SPECIES)

    trend_dt <- lapply(specnames, function(spec) {
        spec_dt <- dt[SPECIES == spec]

        ind_trend_dt <- lapply(unique(spec_dt$SITE), function(ind) {
            ind_dt <- spec_dt[SITE == ind]
            setorder(ind_dt, YEAR)
            ind_dt <- na.omit(ind_dt)

            arow <- tryCatch(
                {
                    if (season == "spring") {
                        pval <- mk.test(ind_dt$spring)$p.value
                        slp <- sens.slope(ind_dt$spring)$estimates
                    } else if (season == "fall") {
                        pval <- mk.test(ind_dt$fall)$p.value
                        slp <- sens.slope(ind_dt$fall)$estimates
                    } else {
                        pval <- mk.test(ind_dt$gsl)$p.value
                        slp <- sens.slope(ind_dt$gsl)$estimates
                    }
                    
                    arow <- data.table(
                        spec = spec, site = ind, slp = slp, pval = pval
                    )
                    return(arow)
                },
                error = function(e) {
                    return(NULL)
                }
            )

            return(arow)
        })
        ind_trend_dt <- do.call(rbind, ind_trend_dt)
    })
    trend_dt <- do.call(rbind, trend_dt)

    return(trend_dt)
}


# Trend analysis in ordinary least square
DoTrendOLS <- function(dt, season = "spring") {
    specnames <- unique(dt$SPECIES)

    trend_dt <- lapply(specnames, function(spec) {
        spec_dt <- dt[SPECIES == spec]

        ind_trend_dt <- lapply(unique(spec_dt$SITE), function(ind) {
            ind_dt <- spec_dt[SITE == ind]
            setorder(ind_dt, YEAR)
            ind_dt <- na.omit(ind_dt)

            arow <- tryCatch(
                {
                    if (season == "spring") {
                        lmfit <- base$FitLm(ind_dt$YEAR, ind_dt$spring)
                    } else if (season == "fall") {
                        lmfit <- base$FitLm(ind_dt$YEAR, ind_dt$fall)
                    } else {
                        lmfit <- base$FitLm(ind_dt$YEAR, ind_dt$gsl)
                    }
                    pval <- lmfit$pval
                    slp <- lmfit$cf[2]

                    arow <- data.table(
                        spec = spec, site = ind, slp = slp, pval = pval
                    )
                    return(arow)
                },
                error = function(e) {
                    return(NULL)
                }
            )

            return(arow)
        })
        ind_trend_dt <- do.call(rbind, ind_trend_dt)
    })
    trend_dt <- do.call(rbind, trend_dt)

    return(trend_dt)
}


# Visualize trends in spring, fall, and growing season length
PlotTrends <- function(spr_trend, atm_trend, gsl_trend, figname) {
    png(
        figname, 
        width = 1200, height = 1600, res = 300
    )
    # par(
    #     bg = NA, fg = "white", col.axis = "white", col.lab = "white",
    #     col.main = "white", col.sub = "white"
    # )
    par(mfrow = c(3, 1), mar = c(2, 3, 1, 1), mgp = c(1.7, 0.7, 0))

    # Spring -----------------
    nonsig <- spr_trend$slp
    sig <- spr_trend[pval < 0.05, slp]
    hist(
        nonsig,
        border = NA, breaks = length(nonsig), xlim = c(-1, 1),
        bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = ""
    )
    if (length(sig) > 0) {
        hist(sig, border = NA, breaks = length(nonsig), col = "seagreen", add = TRUE)
    }
    axis(side = 2, at = seq(0, 10, by = 2), line = -1, cex.axis = 1.5)
    mtext("Frequency", side = 2, line = 1)
    text(
        grconvertX(0.8), grconvertY(0.8, "npc"),
        labels = paste0(
            round(length(sig) / length(nonsig) * 100), "% (", 
            length(sig), " / ", length(nonsig), ")"
        ),
        cex = 1.5, xpd = NA
    )
    # abline(v = 0, col = "red", lwd = 2)

    # Fall -----------------
    nonsig <- atm_trend$slp
    sig <- atm_trend[pval < 0.05, slp]
    hist(
        nonsig,
        border = NA, breaks = length(nonsig), xlim = c(-1, 1),
        bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = ""
    )
    hist(sig, border = NA, breaks = length(nonsig), col = "orange", add = TRUE)
    axis(side = 2, at = seq(0, 10, by = 2), line = -1, cex.axis = 1.5)
    mtext("Frequency", side = 2, line = 1)
    text(
        grconvertX(0.8), grconvertY(0.8, "npc"),
        labels = paste0(
            round(length(sig) / length(nonsig) * 100), "% (",
            length(sig), " / ", length(nonsig), ")"
        ),
        cex = 1.5, xpd = NA
    )
    # abline(v = 0, col = "red", lwd = 2)

    # GSL -----------------
    nonsig <- gsl_trend$slp
    sig <- gsl_trend[pval < 0.05, slp]
    par(mar = c(4, 3, 0, 1))
    hist(
        nonsig,
        border = NA, breaks = length(nonsig), xlim = c(-1, 1),
        bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = ""
    )
    hist(sig, border = NA, breaks = length(nonsig), col = "blue", add = TRUE)
    axis(side = 2, at = seq(0, 10, by = 2), line = -1, cex.axis = 1.5)
    mtext("Frequency", side = 2, line = 1)
    axis(side = 1, cex.axis = 1.5)
    mtext("Slope", side = 1, line = 2.5)
    text(
        grconvertX(0.8), grconvertY(0.8, "npc"),
        labels = paste0(
            round(length(sig) / length(nonsig) * 100), "% (",
            length(sig), " / ", length(nonsig), ")"
        ),
        cex = 1.5, xpd = NA
    )
    # abline(v = 0, col = "red", lwd = 2)

    dev.off()
}



# ~ Load data ####
# ~ ----------------------------------------------------------------------------
hb_pheno <- fread("data/raw/hb_pheno.csv")

# Linear intepolation
hb_pheno <- by(hb_pheno, list(hb_pheno$SITE, hb_pheno$SPECIES), function(sp) {
    # pdf("zzz.pdf", width = 12, height = 4)
    # plot(sp[, .(DATE, Phenology_Stage)])
    # lines(interpo$x, interpo$y, col = "blue")
    # dev.off()
    alldays <- seq(min(sp$DATE), max(sp$DATE), by = "day")
    interpo <- approx(sp$DATE, sp$Phenology_Stage, xout = alldays)

    interpo_dt <- data.table(DATE = interpo$x, Phenology_Stage = interpo$y)
    interpo_dt[, DAY := yday(as_date(DATE))]
    interpo_dt[, YEAR := year(DATE)]
    
    spr <- lapply(unique(interpo_dt$YEAR), function(yr) {
        yr_dt <- interpo_dt[YEAR == yr]
        spr <- yr_dt[DAY < 185][which.min(abs(Phenology_Stage - 3)), ]
        return(spr)
    })
    spr <- data.table(do.call(rbind, spr))[, .(
        YEAR, spring = DAY, spring_date = DATE
    )]

    atm <- lapply(unique(interpo_dt$YEAR), function(yr) {
        yr_dt <- interpo_dt[YEAR == yr]
        atm <- yr_dt[DAY > 185][which.min(abs(Phenology_Stage - 1)), ]
        return(atm)
    })
    atm <- data.table(do.call(rbind, atm))[, .(
        YEAR, fall = DAY, fall_date = DATE
    )]
    
    arow <- cbind(spr, atm)

    arow[, gsl := fall - spring]

    arow$SITE <- sp$SITE[1]
    arow$SPECIES <- sp$SPECIES[1]

    
    return(arow)
})
hb_pheno <- do.call(rbind, hb_pheno)


outdir <- "out/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)



# ~ OLS trends ####
# ~ ----------------------------------------------------------------------------
# OLS trend
spr_trend_ols <- DoTrendOLS(hb_pheno, season = "spring")
atm_trend_ols <- DoTrendOLS(hb_pheno, season = "fall")
gsl_trend_ols <- DoTrendOLS(hb_pheno, season = "gsl")

# fig: HB pheno trends - OLS
PlotTrends(
    spr_trend_ols, atm_trend_ols, gsl_trend_ols,
    figname = file.path(outdir, "hb_pheno_trends_ols.png")
)

# ~ Report some stats ----------------------------------------

# Spring -----------------
spr_trend_ols[pval < 0.05, ]

# Autumn -----------------
atm_trend_ols[pval < 0.05, ]
# Overall fall trends
atm_trend_ols[pval < 0.05, .(mean = mean(slp), sd = sd(slp))]
# Overall fall trends by species
atm_trend_ols[pval < 0.05, .(mean = mean(slp), sd = sd(slp)), keyby = spec]

# GSL -----------------
gsl_trend_ols[pval < 0.05, ]
# Overall fall trends
gsl_trend_ols[pval < 0.05, .(mean = mean(slp), sd = sd(slp))]
# Overall fall trends by species
gsl_trend_ols[pval < 0.05, .(mean = mean(slp), sd = sd(slp)), keyby = spec]


# ~ MK trends ####
# ~ ----------------------------------------------------------------------------
# Mann-Kendall trend
spr_trend_mk <- DoTrendMk(hb_pheno, season = "spring")
atm_trend_mk <- DoTrendMk(hb_pheno, season = "fall")
gsl_trend_mk <- DoTrendMk(hb_pheno, season = "gsl")

# fig: HB pheno trends - Mann-Kendall
PlotTrends(
    spr_trend_mk, atm_trend_mk, gsl_trend_mk,
    figname = file.path(outdir, "hb_pheno_trends_mk.png")
)

# ~ Report some stats ----------------------------------------
# Spring -----------------
spr_trend_mk[pval < 0.05, ]

# Autumn -----------------
atm_trend_mk[pval < 0.05, ]
# Overall fall trends
atm_trend_mk[pval < 0.05, .(mean = mean(slp), sd = sd(slp))]
# Overall fall trends by species
atm_trend_mk[pval < 0.05, .(mean = mean(slp), sd = sd(slp)), keyby = spec]

# GSL -----------------
gsl_trend_mk[pval < 0.05, ]
# Overall fall trends
gsl_trend_mk[pval < 0.05, .(mean = mean(slp), sd = sd(slp))]
# Overall fall trends by species
gsl_trend_mk[pval < 0.05, .(mean = mean(slp), sd = sd(slp)), keyby = spec]

