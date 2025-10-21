# ******************************************************************************
# Trend analysis on MCD12Q2
# 
# Author: Xiaojie Gao
# Date: 2024-07-11
# ******************************************************************************
rm(list=ls())
source("src/base.R")
library(terra)
library(data.table)
library(magrittr)
library(lubridate)
library(ggridges)
library(ggplot2)
library(magrittr)
library(patchwork)
library(terra)
library(RColorBrewer)


# Calculate trends
SiteSpecTrend <- function(specname) {
    spec_dt <- mcd_dt[Category == specname]
    ids <- unique(spec_dt$ID)

    # Iterate all the fit, use a linear regression to get a slope
    trend_dt <- lapply(ids, function(id) {
        thefit <- spec_dt[ID == id,]
        yrs <- year(as_date(thefit$Date))
        thefit[, MidGreenup := yday(as_date(MidGreenup))]
        thefit[, MidGreendown := yday(as_date(MidGreendown))]
        
        thefit[!MidGreenup_QA %in% c("Best", "Good"), MidGreenup := NA]
        thefit[!MidGreendown_QA %in% c("Best", "Good"), MidGreendown := NA]


        # ~ Spring ----------------------------------------
        sos <- thefit[, MidGreenup]

        # Fit a linear regression
        fit <- tryCatch(
            {
                lm(sos ~ yrs)
            },
            error = function(e) {
                return(NULL)
            }
        )

        if (is.null(fit) || is.na(coef(fit)[2])) {
            return(NULL)
        }

        spr_slp <- round(coef(fit)[2], 3)
        f <- summary(fit)$fstatistic
        spr_pval <- pf(f[1], f[2], f[3], lower.tail = FALSE)


        # ~ Autumn ----------------------------------------
        eos <- thefit[, MidGreendown]

        # Fit a linear regression
        fit <- tryCatch(
            {
                lm(eos ~ yrs)
            },
            error = function(e) {
                return(NULL)
            }
        )

        if (is.null(fit) || is.na(coef(fit)[2])) {
            return(NULL)
        }

        atm_slp <- round(coef(fit)[2], 3)
        f <- summary(fit)$fstatistic
        atm_pval <- pf(f[1], f[2], f[3], lower.tail = FALSE)


        # ~ Growing season length ----------------------------------------
        gsl <- eos - sos
        # Fit a linear regression
        fit <- tryCatch(
            {
                lm(gsl ~ yrs)
            },
            error = function(e) {
                return(NULL)
            }
        )

        if (is.null(fit) || is.na(coef(fit)[2])) {
            return(NULL)
        }

        gsl_slp <- round(coef(fit)[2], 3)
        f <- summary(fit)$fstatistic
        gsl_pval <- pf(f[1], f[2], f[3], lower.tail = FALSE)


        # Get latitude
        lat <- spec_dt[ID == id, Latitude][1]
        lon <- spec_dt[ID == id, Longitude][1]

        arow <- data.table(
            spec = specname, ID = id, lat = lat, lon = lon,
            spr_slp = spr_slp, spr_pval = spr_pval,
            atm_slp = atm_slp, atm_pval = atm_pval,
            gsl_slp = gsl_slp, gsl_pval = gsl_pval
        )
        return(arow)
    })
    trend_dt <- do.call(rbind, trend_dt)

    return(trend_dt)
}


# Plot the trend as points on maps
PlotTrendPts <- function(site_trend_vect, specname, varname) {
    plot(r9_shp, col = "grey", lty = 2)
    mtext(paste(specname), line = 0, cex = 1.5)
    # mtext(paste0(
    #     gsub("\r\n", " ", spec_trend_vect$CommonName[1]),
    #     " (", spec_trend_vect$Type[1], ")"
    # ), line = -2, cex = 1.5)

    text(-80, 47, labels = varname)

    ptcol <- sapply(1:nrow(spec_trend_vect), function(i) {
        theslp <- spec_trend_vect[i, ][[paste0(varname, "_slp")]]
        if (theslp < 0) {
            thecol <- "blue"
        } else {
            thecol <- "red"
        }

        thepval <- spec_trend_vect[i, ][[paste0(varname, "_pval")]]
        if (!is.na(thepval) & thepval > 0.05) {
            thecol <- adjustcolor(thecol, 0.1)
        }

        return(thecol)
    })

    slp <- spec_trend_vect[[paste0(varname, "_slp"), drop = TRUE]]
    points(spec_trend_vect,
        pch = 21,
        cex = Scale2AnyRange(slp) * 2,
        bg = ptcol, col = NA
    )

    legend_pt_range <- c(0.1, 0.5, 0.9)
    legend_ptcex <- Scale2AnyRange(
        legend_pt_range,
        smin = min(slp),
        smax = max(slp)
    ) * 2
    legend(-74.5, 39,
        bty = "n", legend = legend_pt_range,
        pt.cex = legend_ptcex, pt.bg = "grey", pch = 21, xpd = NA
    )
    if (varname == "gsl") {
        legend(-71.5, 39,
            bty = "n", legend = c("Shortened", "Extended"),
            pch = 21, pt.bg = c("blue", "red"), xpd = NA
        )
    } else {
        legend(-71.5, 39,
            bty = "n", legend = c("Advanced", "Delayed"),
            pch = 21, pt.bg = c("blue", "red"), xpd = NA
        )
    }
}



mcd_dt <- fread("data/fia_dom_mcd12q2.csv")
specnames <- unique(mcd_dt$Category)

site_trend_dt <- NULL
for (specname in specnames) {
    site_trend_dt <- rbind(site_trend_dt, SiteSpecTrend(specname))
}

fwrite(site_trend_dt, "pipe/fia_dom_mcd12q2_trend.csv")

site_trend_dt <- fread("pipe/fia_dom_mcd12q2_trend.csv")


# ~ Map for each trend ----------------------------------------
r9_shp <- vect(base$region9_shpfile)

site_trend_vect <- vect(site_trend_dt,
    geom = c("lon", "lat"), crs = "epsg:4326"
)
cols <- hcl.colors(2, "dark 2")


# fig: Map mcd12q2 trends
pdf("out/spec_mcd12q2_trend_map.pdf", width = 15, height = 4)
par(mfrow = c(1, 3))
for (specname in sort(unique(site_trend_dt$spec))) {
    spec_trend_vect <- site_trend_vect[site_trend_vect$spec == specname]
    PlotTrendPts(spec_trend_vect, specname, "spr")
    PlotTrendPts(spec_trend_vect, specname, "atm")
    PlotTrendPts(spec_trend_vect, specname, "gsl")
}
dev.off()




# ~ Ridge plot for trends ----------------------------------------
# fig: Ridge plot for trends

# Spring
sig_trend_dt <- site_trend_dt[spr_pval < 0.05, ]
plt_spr <- ggplot(sig_trend_dt, aes(x = spr_slp, y = spec)) +
    geom_density_ridges(fill = "seagreen") +
    theme_ridges() +
    xlim(-1, 1) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme(legend.position = "none")

# Autumn
sig_trend_dt <- site_trend_dt[atm_pval < 0.05, ]
plt_atm <- ggplot(sig_trend_dt, aes(x = atm_slp, y = spec)) +
    geom_density_ridges(fill = "orange") +
    theme_ridges() +
    xlim(-1, 1) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme(legend.position = "none")

# GSL
sig_trend_dt <- site_trend_dt[gsl_pval < 0.05, ]
plt_gsl <- ggplot(sig_trend_dt, aes(x = gsl_slp, y = spec)) +
    geom_density_ridges(fill = "#e600ff") +
    theme_ridges() +
    xlim(-1, 1) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme(legend.position = "none")

plt <- plt_spr + plt_atm + plt_gsl


ggsave(filename = "out/spec_mcd12q2_trend.jpg", plt, width = 18, height = 6)



# ~ Summary ----------------------------------------

by(site_trend_dt, site_trend_dt$spec, function(sp_dt) {
    lm(gsl_slp ~ scale(spr_slp) + scale(atm_slp), 
        data = sp_dt[spr_pval < 0.05 | atm_pval < 0.05]
    ) %>%
        summary()
})

site_trend_dt[, cor(gsl_slp, spr_slp), keyby = "spec"]
site_trend_dt[, cor(gsl_slp, atm_slp), keyby = "spec"]


total <- site_trend_dt[, .(TotN = .N), keyby = "spec"]

spr_sig <- merge(
    site_trend_dt[spr_pval < 0.05, .N, keyby = "spec"],
    total,
    by = "spec"
)
spr_sig[, sig_pct := N / TotN * 100]

atm_sig <- merge(
    site_trend_dt[atm_pval < 0.05, .N, keyby = "spec"],
    total,
    by = "spec"
)
atm_sig[, sig_pct := N / TotN * 100]

gsl_sig <- merge(
    site_trend_dt[gsl_pval < 0.05, .N, keyby = "spec"],
    total,
    by = "spec"
)
gsl_sig[, sig_pct := N / TotN * 100]

spr_sig
atm_sig
gsl_sig


par(mfrow = c(2, 5))
null <- lapply(specnames, function(sp) {
    hist(site_trend_dt[spec == sp & spr_slp > 0, spr_slp], 
        xlim = c(-1, 1), ylim = c(0, 50),
        breaks = 50, border = NA, col = "red"
    )
    hist(site_trend_dt[spec == sp & spr_slp < 0, spr_slp], 
        breaks = 50, border = NA, add = TRUE, col = "seagreen"
    )
})

hist(site_trend_dt[spr_pval < 0.05 & spr_slp > 0, spr_slp], breaks = 100, 
    border = NA, col = "red", xlim = c(-2, 2), ylim = c(0, 100)
)
hist(site_trend_dt[spr_pval < 0.05 & spr_slp < 0, spr_slp], breaks = 200, 
    border = NA, col = "seagreen", add = TRUE
)

site_trend_dt[spr_pval < 0.05 & spr_slp > 0]
site_trend_dt[spr_pval < 0.05 & spr_slp < 0]
