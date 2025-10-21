# ******************************************************************************
# Make a map to show all deciduous species that are dominant species of FIA
# sites and their locations.
# 
# Author: Xiaojie Gao
# Date: 2024-06-17
# ******************************************************************************
rm(list=ls())

source("src/base.R")
library(data.table)
library(magrittr)
library(terra)
library(rnaturalearth)
library(RColorBrewer)



r9_shp <- vect(base$region9_shpfile)
r9_bbox <- as.polygons(ext(r9_shp), crs = crs(r9_shp))


blsp_fit_dt <- fread("data/fia_blsp_fit.csv")

# Get species names
latin_names <- fread("data/common10.csv")
latin_names[, latin := paste(GENUS, SPECIES_NAME)]

specnames <- unique(blsp_fit_dt$spec) %>%
    sort()
specnames_latin <- latin_names[gesp %in% specnames, latin]

blsp_fit_vect <- vect(blsp_fit_dt, geom = c("x", "y"), crs = "epsg:4326")


states <- ne_states("United States of America")
states <- states[!states$name %in% c("Alaska", "Hawaii"), ] %>%
    vect()



outdir <- "out/species_phenology/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)



# fig: out/fia_dom_spec_sites.png
png(file.path(outdir, "fia_dom_spec_sites.png"), width = 2000, height = 1300, res = 300)
par(mar = c(3, 3, 1, 1))

plot(r9_shp, col = "grey", lty = 2, xlab = "Longitude", ylab = "Latitude", box = FALSE)
# cols <- brewer.pal(length(specnames), "Paired")
cols <- c(
    "#A6CEE3", "#1F78B4",
    "#B2DF8A",
    "#33A02C",
    "#FDBF6F",
    "#ff787b", "#E31A1C", "#960808", "#b70e6e",
    "#a25feb"
)
for (i in seq_along(specnames)) {
    spec <- specnames[i]
    # points(
    #     blsp_fit_vect[blsp_fit_vect$spec == spec, ], 
    #     pch = 16, col = adjustcolor(cols[i], 0.1), cex = 0.5
    # )
    points(
        blsp_fit_vect[blsp_fit_vect$spec == spec, ], 
        pch = 21, bg = adjustcolor(cols[i], 0.1), col = "white", 
        cex = 0.6, lwd = 0.5
    )
}

latin_names <- lapply(specnames_latin, function(x) {
    bquote(italic(.(x)))
})

legend(
    grconvertX(0.2, "ndc"), grconvertY(0.98, "ndc"), bty = "n",
    legend = do.call(expression, latin_names), 
    pch = 21, pt.bg = cols, col = "white", 
    ncol = 4, cex = 0.65,
    xpd = NA
)

text(
    -83, 36.2, 
    paste(
        "Number of plots:", 
        format(length(unique(blsp_fit_dt$ID)), big.mark = ",")
    )
)


# Overview
par(fig = c(0.73, 0.92, 0.1, 0.4), new = TRUE)
plot(states, col = "grey", bty = "n", 
    axes = FALSE, mar = c(0, 0, 0, 0)
)
lines(r9_bbox, col = "blue")

dev.off()

