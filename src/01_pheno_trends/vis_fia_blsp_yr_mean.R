# ******************************************************************************
# Make maps of 40-year mean BLSP for all deciduous species
# 
# Author: Xiaojie Gao
# Date: 2024-06-20
# ******************************************************************************
library(terra)
library(data.table)
library(magrittr)
library(RColorBrewer)



Scale2AnyRange <- function(val, tmin = 0, tmax = 1, smin = NULL, smax = NULL) {
    smin <- ifelse(is.null(smin), min(val, na.rm = TRUE), smin)
    smax <- ifelse(is.null(smax), max(val, na.rm = TRUE), smax)

    sval <- (val - smin) / (smax - smin) * (tmax - tmin) + tmin
    attr(sval, "smin") <- smin
    attr(sval, "smax") <- smax

    return(sval)
}



r9 <- vect(base$region9_shpfile)

blsp <- fread("data/fia_blsp_fit.csv")

# Only do the 10 major decidous species
specnames <- c(
    "ACRU", "ACSA", "FRNI", "LALA", "POTR", "QUAL", "QUPR", "QURU",
    "QUVE", "TIAM"
)

blsp_yr_mean <- blsp[
    spec %in% specnames & 
    MidGreenup_upr - MidGreenup_lwr < 30 & 
    MidGreendown_upr - MidGreendown_lwr < 30, .(
        midgup_mean = mean(MidGreenup),
        midgd_mean = mean(MidGreendown),
        gsl_mean = mean(MidGreendown - MidGreenup),
        x, y, spec
), keyby = ID] %>%
    unique()

blsp_yr_mean <- vect(blsp_yr_mean, geom = c("x", "y"), crs = "epsg:4326")



# fig: 
cols <- brewer.pal(length(specnames), "Paired")
spec_cols <- data.table(cbind(specnames, cols))


{
png("out/blsp_yr_mean.png", width = 2500, height = 2400, res = 300)
par(mfrow = c(3, 2), mar = c(0, 0, 0, 0), oma = c(0, 3, 2, 0))

scale_factor <- 1.5

# ~ Spring ----------------------------------------
plot(r9, col = adjustcolor("grey", 0.2), axes = FALSE)
blsp_yr_mean$spr_cex <- Scale2AnyRange(blsp_yr_mean$midgup_mean) * scale_factor
for (spec in c("ACRU", "ACSA", "QURU", "QUPR", "QUVE")) {
    points(blsp_yr_mean[blsp_yr_mean$spec == spec, "midgup_mean"],
        cex = blsp_yr_mean[blsp_yr_mean$spec == spec]$spr_cex,
        pch = 21,
        bg = spec_cols[specnames == spec, cols], col = "white",
        alpha = 0.5
    )
}

plot(r9, col = adjustcolor("grey", 0.2), axes = FALSE)
for (spec in c("FRNI", "LALA", "QUAL", "POTR", "TIAM")) {
    points(blsp_yr_mean[blsp_yr_mean$spec == spec, "midgup_mean"],
        cex = blsp_yr_mean[blsp_yr_mean$spec == spec, ]$spr_cex,
        pch = 21,
        bg = spec_cols[specnames == spec, cols], col = "white",
        alpha = 0.5
    )
}


# ~ Autumn ----------------------------------------
plot(r9, col = adjustcolor("grey", 0.2), axes = FALSE)
blsp_yr_mean$atm_cex <- Scale2AnyRange(blsp_yr_mean$midgd_mean) * scale_factor
for (spec in c("ACRU", "ACSA", "QURU", "QUPR", "QUVE")) {
    points(blsp_yr_mean[blsp_yr_mean$spec == spec, "midgd_mean"],
        cex = blsp_yr_mean[blsp_yr_mean$spec == spec, ]$atm_cex,
        pch = 21,
        bg = spec_cols[specnames == spec, cols], col = "white",
        alpha = 0.5
    )
}

plot(r9, col = adjustcolor("grey", 0.2), axes = FALSE)
for (spec in c("FRNI", "LALA", "QUAL", "POTR", "TIAM")) {
    points(blsp_yr_mean[blsp_yr_mean$spec == spec, "midgd_mean"],
        cex = blsp_yr_mean[blsp_yr_mean$spec == spec, ]$atm_cex,
        pch = 21,
        bg = spec_cols[specnames == spec, cols], col = "white",
        alpha = 0.5
    )
}


# ~ GSL ----------------------------------------
plot(r9, col = adjustcolor("grey", 0.2), axes = FALSE)
blsp_yr_mean$gsl_cex <- Scale2AnyRange(blsp_yr_mean$gsl_mean) * scale_factor
for (spec in c("ACRU", "ACSA", "QURU", "QUPR", "QUVE")) {
    points(blsp_yr_mean[blsp_yr_mean$spec == spec, "gsl_mean"],
        cex = blsp_yr_mean[blsp_yr_mean$spec == spec, ]$gsl_cex,
        pch = 21,
        bg = spec_cols[specnames == spec, cols], col = "white",
        alpha = 0.5
    )
}

plot(r9, col = adjustcolor("grey", 0.2), axes = FALSE)
for (spec in c("FRNI", "LALA", "QUAL", "POTR", "TIAM")) {
    points(blsp_yr_mean[blsp_yr_mean$spec == spec, "gsl_mean"],
        cex = blsp_yr_mean[blsp_yr_mean$spec == spec, ]$gsl_cex,
        pch = 21,
        bg = spec_cols[specnames == spec, cols], col = "white",
        alpha = 0.5
    )
}


# ~ Annotations and legends ----------------------------------------

# Spring -----------------

text(grconvertX(0.05, "ndc"), grconvertY(0.82, "ndc"),
    pos = 1,
    labels = "Spring", xpd = NA, cex = 2, 
    srt = 90
)

spr_legend_labels <- c(120, 130, 140, 150, 160)
spr_scale <- Scale2AnyRange(blsp_yr_mean$midgup_mean)
spr_legend_cex <- Scale2AnyRange(spr_legend_labels, 
    smin = attr(spr_scale, "smin"),
    smax = attr(spr_scale, "smax")
) * scale_factor
legend(grconvertX(0.5, "ndc"), grconvertY(0.82, "ndc"),
    bty = "n", xpd = NA,
    legend = spr_legend_labels,
    pch = 16, pt.cex = spr_legend_cex,
    title = "DOY"
)


# Fall -----------------

text(grconvertX(0.05, "ndc"), grconvertY(0.5, "ndc"),
    pos = 1,
    labels = "Fall", xpd = NA, cex = 2,
    srt = 90
)
atm_legend_labels <- c(250, 260, 270, 280, 290)
atm_scale <- Scale2AnyRange(blsp_yr_mean$midgd_mean)
atm_legend_cex <- Scale2AnyRange(atm_legend_labels,
    smin = attr(atm_scale, "smin"),
    smax = attr(atm_scale, "smax")
) * scale_factor
legend(grconvertX(0.5, "ndc"), grconvertY(0.5, "ndc"),
    bty = "n", xpd = NA,
    legend = atm_legend_labels,
    pch = 16, pt.cex = atm_legend_cex,
    title = "DOY"
)

# GSL -----------------
gsl_legend_labels <- c(100, 120, 140, 160, 180)
gsl_scale <- Scale2AnyRange(blsp_yr_mean$gsl_mean)
gsl_legend_cex <- Scale2AnyRange(gsl_legend_labels,
    smin = attr(gsl_scale, "smin"),
    smax = attr(gsl_scale, "smax")
) * scale_factor
legend(grconvertX(0.5, "ndc"), grconvertY(0.18, "ndc"),
    bty = "n", xpd = NA,
    legend = gsl_legend_labels,
    pch = 16, pt.cex = gsl_legend_cex,
    title = "DOY"
)

text(grconvertX(0.05, "ndc"), grconvertY(0.18, "ndc"),
    pos = 1,
    labels = "Growing season length", xpd = NA, cex = 2,
    srt = 90
)


# Species legend
legend(grconvertX(0.5, "ndc"), grconvertY(0.99, "ndc"),
    bty = "n", xpd = NA, xjust = 0.5,
    legend = specnames, pch = 16, col = spec_cols$cols,
    ncol = 10, cex = 1.2
)



dev.off()

}
