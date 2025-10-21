# ******************************************************************************
# Species-specific phenology trends using HF ground data.
# 
# Author: Xiaojie Gao
# Date: 2024-05-31
# ******************************************************************************
rm(list = ls())
source("src/base.R")
library(data.table)
library(trend)




DoTrendMK <- function(dt, season = "spring") {
    specnames <- unique(dt$species)

    trend_dt <- lapply(specnames, function(spec) {
        spec_dt <- dt[species == spec]

        ind_trend_dt <- lapply(unique(spec_dt$tree.id), function(ind) {
            ind_dt <- spec_dt[tree.id == ind]
            setorder(ind_dt, year)
            ind_dt <- na.omit(ind_dt)
            if (nrow(ind_dt) < 10) {
                return(NULL)
            }
            
            arow <- tryCatch(
                {
                    if (season == "spring") {
                        pval <- mk.test(ind_dt$l75.doy)$p.value
                        slp <- sens.slope(ind_dt$l75.doy)$estimates
                    } else if (season == "fall") {
                        pval <- mk.test(ind_dt$lf.doy)$p.value
                        slp <- sens.slope(ind_dt$lf.doy)$estimates
                    } else {
                        pval <- mk.test(ind_dt$gsl)$p.value
                        slp <- sens.slope(ind_dt$gsl)$estimates
                    }
                    
                    arow <- data.table(
                        spec = spec, slp = slp, pval = pval
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


DoTrendOLS <- function(dt, season = "spring") {
    specnames <- unique(dt$species)

    trend_dt <- lapply(specnames, function(spec) {
        spec_dt <- dt[species == spec]

        ind_trend_dt <- lapply(unique(spec_dt$tree.id), function(ind) {
            ind_dt <- spec_dt[tree.id == ind]
            setorder(ind_dt, year)
            ind_dt <- na.omit(ind_dt)
            if (nrow(ind_dt) < 10) {
                return(NULL)
            }

            arow <- tryCatch(
                {
                    if (season == "spring") {
                        lmfit <- base$FitLm(ind_dt$year, ind_dt$l75.doy)
                    } else if (season == "fall") {
                        lmfit <- base$FitLm(ind_dt$year, ind_dt$lf.doy)
                    } else {
                        lmfit <- base$FitLm(ind_dt$year, ind_dt$gsl)
                    }

                    pval <- lmfit$pval
                    slp <- lmfit$cf[2]

                    arow <- data.table(
                        spec = spec, tree.id = ind, slp = slp, pval = pval
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
        border = NA, breaks = length(nonsig), xlim = c(-5, 5),
        bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = ""
    )
    hist(sig, border = NA, breaks = length(nonsig), col = "seagreen", add = TRUE)
    axis(side = 2, at = c(0, 10, 20), line = -2, cex.axis = 1.5)
    mtext("Frequency", side = 2)
    text(
        grconvertX(3), grconvertY(0.8, "npc"),
        labels = paste0(
            round(length(sig) / length(nonsig) * 100), "% (",
            length(sig), " / ", length(nonsig), ")"
        ),
        cex = 1.5
    )

    # Fall -----------------
    nonsig <- atm_trend$slp
    sig <- atm_trend[pval < 0.05, slp]
    hist(
        nonsig,
        border = NA, breaks = length(nonsig), xlim = c(-5, 5),
        bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = ""
    )
    hist(sig, border = NA, breaks = length(nonsig), col = "orange", add = TRUE)
    axis(side = 2, at = seq(0, 25, by = 5), line = -2, cex.axis = 1.5)
    mtext("Frequency", side = 2)
    text(
        grconvertX(3), grconvertY(0.8, "npc"),
        labels = paste0(
            round(length(sig) / length(nonsig) * 100), "% (",
            length(sig), " / ", length(nonsig), ")"
        ),
        cex = 1.5
    )

    # GSL -----------------
    nonsig <- gsl_trend$slp
    sig <- gsl_trend[pval < 0.05, slp]
    par(mar = c(4, 3, 0, 1))
    hist(
        nonsig,
        border = NA, breaks = length(nonsig), xlim = c(-5, 5),
        bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = ""
    )
    hist(sig, border = NA, breaks = length(nonsig), col = "blue", add = TRUE)
    axis(side = 2, at = seq(0, 25, by = 5), line = -2, cex.axis = 1.5)
    mtext("Frequency", side = 2)
    axis(side = 1, cex.axis = 1.5)
    mtext("Slope", side = 1, line = 2.5)
    text(
        grconvertX(3), grconvertY(0.8, "npc"),
        labels = paste0(
            round(length(sig) / length(nonsig) * 100), "% (",
            length(sig), " / ", length(nonsig), ")"
        ),
        cex = 1.5
    )

    dev.off()
}



# ~ Load data ####
# ~ ----------------------------------------------------------------------------
spr_dt <- fread("https://harvardforest.fas.harvard.edu/data/p00/hf003/hf003-05-spring-mean-ind.csv")
atm_dt <- fread("https://harvardforest.fas.harvard.edu/data/p00/hf003/hf003-07-fall-mean-ind.csv")
com_dt <- merge(spr_dt, atm_dt, by = c("year", "tree.id", "species"))

# Select only trees w/ longer than 10 years of data
spr_tree_records <- spr_dt[, .(num_records = .N), by = tree.id]
spr_dt <- spr_dt[tree.id %in% spr_tree_records[num_records > 10, tree.id]]

atm_tree_records <- atm_dt[, .(num_records = .N), by = tree.id]
atm_dt <- atm_dt[tree.id %in% atm_tree_records[num_records > 10, tree.id]]

com_dt <- merge(spr_dt, atm_dt, by = c("year", "tree.id", "species"))
# Remove the trees that only have fall phenology
com_dt <- na.omit(com_dt[, .(year, tree.id, species, l75.doy, lf.doy)])
com_dt[, gsl := lf.doy - l75.doy]

# Select only deciduous tree species
deci_spec <- c(
    "ACPE", "ACRU", "ACSA", "AMSP", "BEAL", "BELE", "BEPA", "BEPO", "CADE",
    "FAGR", "FRAM", "NYSY", "POTR", "PRSE",
    "QUAL", "QURU", "QUVE"
)
com_dt <- com_dt[species %in% deci_spec,]

uniqueN(com_dt$species)
uniqueN(com_dt$tree.id)


figdir <- "out/"
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)



# ~ OLS trends ####
# ~ ----------------------------------------------------------------------------
spr_trend_ols <- DoTrendOLS(
    com_dt[, .(year, tree.id, species = species, l75.doy)], 
    season = "spring"
)
atm_trend_ols <- DoTrendOLS(
    com_dt[, .(year, tree.id, species = species, lf.doy)], 
    season = "fall"
)
gsl_trend_ols <- DoTrendOLS(
    com_dt[, .(year, tree.id, species = species, gsl)],
    season = "gsl"
)

# fig: HF pheno trends - OLS
PlotTrends(
    spr_trend_ols, atm_trend_ols, gsl_trend_ols,
    figname = file.path(figdir, "hf_pheno_trends_ols.png")
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
spr_trend_mk <- DoTrendMK(
    com_dt[, .(year, tree.id, species = species, l75.doy)],
    season = "spring"
)
atm_trend_mk <- DoTrendMK(
    com_dt[, .(year, tree.id, species = species, lf.doy)],
    season = "fall"
)
gsl_trend_mk <- DoTrendMK(com_dt[, .(year, tree.id, species = species, gsl)],
    season = "gsl"
)

# fig: HF pheno trends - OLS
PlotTrends(
    spr_trend_ols, atm_trend_ols, gsl_trend_ols,
    figname = file.path(figdir, "hf_pheno_trends_mk.png")
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




