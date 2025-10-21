# ******************************************************************************
# Visualize dominant species site locations and phenology from MSLSP.
# 
# Author: Xiaojie Gao
# Date: 2024-04-06
# ******************************************************************************
library(data.table)
library(terra)
library(tmap)
source("src/base.R")



PlotSpecPheno <- function(spname, vari) {
    # Use multi-year mean for phenology
    sp_dt <- dt[spec == spname, .(pheno = mean(get(vari))), by = "X"]
    site_vec <- vect(
        fia_domspec_sites[spec == spname, .(X, x, y, Type, CommonName)], 
        geom = c("x", "y"), crs = "EPSG:4326"
    ) %>%
        project(r9_shp) %>%
        merge(sp_dt, by = "X")
    
    # Number of sites
    npt <- nrow(site_vec)

    # Normalize pheno for drawing point sizes
    min_pheno <- min(site_vec$pheno)
    max_pheno <- max(site_vec$pheno)
    ptcex <- (site_vec$pheno - min_pheno) / 
        (max_pheno - min_pheno) * (3 - 0.5) + 0.5

    cols <- hcl.colors(2, "dark 2")
    bgcol <- ifelse(site_vec$Type[1] == "deci", cols[1], cols[2])

    par(fig = c(0, 1, 0, 1))
    plot(r9_shp, col = "grey", lty = 2)
    mtext(paste(spname, "-", vari), line = 0, cex = 1.5)
    mtext(
        paste0(
            gsub("\r\n", " ", site_vec$CommonName[1]), 
            " (", site_vec$Type[1], ")"
        ), 
        line = -2, cex = 1.5
    )
    points(site_vec, pch = 21, cex = ptcex, bg = bgcol, col = "white")

    # Draw legend
    qt <- quantile(site_vec$pheno, c(0.25, 0.75))
    qt_ptcex <- quantile(ptcex, c(0, 0.25, 0.75, 1))
    legend(grconvertX(0.5, "npc"), grconvertY(0.8, "npc"), bty = 'n', text.width = 2,
        legend = c(
            round(min_pheno), round(qt[1]), round(qt[2]), round(max_pheno)
        ),
        pch = rep(21, 4), col = rep("white", 4), pt.bg = rep(bgcol, 4), 
        pt.cex = qt_ptcex,
        ncol = 4, xpd = NA
    )

    # Add a boxplot here
    par(fig = c(0.8, 0.98, 0.05, 0.4), new = TRUE, bty = "n", 
        mar = c(4, 5, 1, 1)
    )
    boxplot(site_vec$pheno, 
        xaxt = "n", xlab = "", ylab = vari, 
        las = 1
    )
}




fia_domspec_sites <- fread("data/fia_domspec_truexy.csv")
sp_fac <- unique(fia_domspec_sites$spec)
sp_fac <- sort(sp_fac)

dt <- fread("data/sp_60_mslsp_dt.csv")


r9_shp <- vect(base$region9_shpfile)


# fig: pipe/domspec_mslsp_sites.pdf
pdf("pipe/domspec_mslsp_sites.pdf", width = 12, height = 8)

for (metric in colnames(dt)[12:21]) {
    for (spname in sp_fac) {
        PlotSpecPheno(spname, vari = metric)
    }
}

dev.off()



