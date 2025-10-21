# ******************************************************************************
# Visualize Daymet climate trends at FIA sites w/ dominant species.
# 
# Author: Xiaojie Gao
# Date: 2024-12-07
# ******************************************************************************
rm(list = ls())
library(data.table)
library(lubridate)
library(svglite)


PlotTrends <- function(trend_file, figname) {
    trend_dt <- fread(trend_file)

    # Make a figure
    svglite(figname, width = 8, height = 2)
    # png("pipe/temp_trends_dark.png", width = 2100, height = 700, res = 300)
    # par(
    #     bg = NA, fg = "white", col.axis = "white", col.lab = "white",
    #     col.main = "white", col.sub = "white"
    # )

    par(
        mfrow = c(1, 4), mgp = c(1.5, 0.5, 0), mar = c(3, 3, 2, 1),
        cex.lab = 1.2, cex.axis = 1.2
    )
    seasons <- c("spr", "smr", "atm", "ann")
    for (sea in seasons) {
        sea_trend <- trend_dt[season == sea]
        sea_col <- switch(sea, 
            "spr" = "seagreen",
            "smr" = "red",
            "atm" = "orange",
            "ann" = "blue"
        )
        sea_mth <- switch(sea, 
            "spr" = "March-May",
            "smr" = "Jun-Aug",
            "atm" = "Sep-Nov",
            "ann" = "Annual"
        )

        qt <- quantile(sea_trend$slp, c(0, 0.01, 0.99, 1))
        lim <- max(abs(qt[2:3]))
        # Summer Srad has large outliners
        nbreak <- ifelse(sea == "smr" && grepl("_srad_", trend_file), 100, 50)
        hist(
            sea_trend$slp,
            breaks = nbreak, border = "white",
            xlab = sea_mth, ylab = "Frequency", main = "",
            xaxt = "n",
            xlim = c(-lim, lim)
        )
        axis(side = 1, at = c(-round(lim, 2), 0, round(lim, 2)))
        abline(v = 0, lty = 2)
        hist(
            sea_trend[pval < 0.05]$slp,
            breaks = nbreak, col = sea_col, border = "white", border.lwd = 0.5,
            add = TRUE
        )
    }

    dev.off()
}



figdir <- "out/"
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

# Temperature
PlotTrends(
    trend_file = "pipe/fia_domspec_daymet_temp_trends_ols.csv",
    figname = file.path(figdir, "temp_trends.svg")
)
# Precipitation
PlotTrends(
    trend_file = "pipe/fia_domspec_daymet_prcp_trends_ols.csv",
    figname = file.path(figdir, "prcp_trends.svg")
)
# Srad
PlotTrends(
    trend_file = "pipe/fia_domspec_daymet_srad_trends_ols.csv",
    figname = file.path(figdir, "srad_trends.svg")
)

