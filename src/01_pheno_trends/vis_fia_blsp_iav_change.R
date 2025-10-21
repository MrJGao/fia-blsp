# ******************************************************************************
# Visualize site-specific trends of IAV for spring/fall/gsl
# 
# Author: Xiaojie Gao
# Date: 2025-05-21
# ******************************************************************************
rm(list = ls())
source("src/base.R")
library(data.table)
library(terra)
library(tmap)
library(rnaturalearth)



if (!exists("var_dt")) {
    source("src/01_pheno_trends/mod_fia_blsp_iav_change.R")
}


# Percent of ascending and descending
sos_iav_asc <- var_dt[sos_iav_pval < 0.05 & sos_iav_slp > 0, .N / nrow(var_dt) * 100]
sos_iav_desc <- var_dt[sos_iav_pval < 0.05 & sos_iav_slp < 0, .N / nrow(var_dt) * 100]

eos_iav_asc <- var_dt[eos_iav_pval < 0.05 & eos_iav_slp > 0, .N / nrow(var_dt) * 100]
eos_iav_desc <- var_dt[eos_iav_pval < 0.05 & eos_iav_slp < 0, .N / nrow(var_dt) * 100]

gsl_iav_asc <- var_dt[gsl_iav_pval < 0.05 & gsl_iav_slp > 0, .N / nrow(var_dt) * 100]
gsl_iav_desc <- var_dt[gsl_iav_pval < 0.05 & gsl_iav_slp < 0, .N / nrow(var_dt) * 100]


{ # fig: Var(IAV) trends
    svglite::svglite(
        "pipe/species_pheno/var_iav_trends.svg", 
        height = 2.5, width = 8, pointsize = 13
    )
    par(
        mfrow = c(1, 3), bty = "L", mgp = c(1.7, 0.5, 0), mar = c(3, 3, 1, 1),
        cex.lab = 1.2, tck = -0.03
    )

    # SOS
    hist(
        var_dt[abs(sos_iav_slp) < 5, sos_iav_slp],
        breaks = seq(-5, 5, 0.1),
        border = NA,
        xlim = c(-5, 5),
        xlab = "ΔIAV SOS",
        main = ""
    )
    hist(
        var_dt[sos_iav_pval < 0.05 & abs(sos_iav_slp) < 5, sos_iav_slp],
        breaks = seq(-5, 5, 0.1),
        include.lowest = FALSE, right = FALSE,
        border = NA, col = "seagreen",
        add = TRUE
    )
    box(bty = "L")
    abline(v = 0, lty = 2)
    legend(
        "topleft", bty = "n", 
        legend = c(
            paste0("INC.: ", round(sos_iav_asc), "%"), 
            paste0("DEC.: ", round(sos_iav_desc), "%")
        )
    )

    # EOS
    hist(
        var_dt[abs(eos_iav_slp) < 5, eos_iav_slp],
        breaks = seq(-5, 5, 0.1),
        border = NA,
        xlim = c(-5, 5),
        xlab = "ΔIAV EOS",
        main = ""
    )
    hist(
        var_dt[eos_iav_pval < 0.05 & abs(eos_iav_slp) < 5, eos_iav_slp],
        breaks = seq(-5, 5, 0.1),
        include.lowest = FALSE, right = FALSE,
        border = NA, col = "orange",
        add = TRUE
    )
    box(bty = "L")
    abline(v = 0, lty = 2)
    legend(
        "topleft", bty = "n", 
        legend = c(
            paste0("INC.: ", round(eos_iav_asc), "%"), 
            paste0("DEC.: ", round(eos_iav_desc), "%")
        )
    )
    
    
    # GSL
    hist(
        var_dt[abs(gsl_iav_slp) < 5, gsl_iav_slp],
        breaks = seq(-5, 5, 0.1),
        border = NA,
        xlim = c(-5, 5),
        xlab = "ΔIAV GSL",
        main = ""
    )
    hist(
        var_dt[gsl_iav_pval < 0.05 & abs(gsl_iav_slp) < 5, gsl_iav_slp],
        breaks = seq(-5, 5, 0.1),
        include.lowest = FALSE, right = FALSE,
        border = NA, col = "blue",
        add = TRUE
    )
    box(bty = "L")
    abline(v = 0, lty = 2)
    legend(
        "topleft", bty = "n", 
        legend = c(
            paste0("INC.: ", round(gsl_iav_asc), "%"), 
            paste0("DEC.: ", round(gsl_iav_desc), "%")
        )
    )

    dev.off()
}



# ~ Relative contributions of SOS and EOS to GSL ####
# ~ ----------------------------------------------------------------------------
contri_dt <- var_dt[
    gsl_iav_pval < 0.05,
    .(
        c_sos = sos_iav_slp / gsl_iav_slp * 100, 
        c_eos = eos_iav_slp / gsl_iav_slp * 100
    )
]
# tt <- var_dt[
#     ,
#     .(c_sos = sos_iav_slp / gsl_iav_slp, c_eos = eos_iav_slp / gsl_iav_slp)
# ]

contri_dt[, diff := c_sos - c_eos]

contri_dt[c_sos >= c_eos, .N / nrow(contri_dt) * 100]
contri_dt[c_sos < c_eos, .N / nrow(contri_dt) * 100]


sig_var_dt <- var_dt[gsl_iav_pval < 0.05, ]

pcor_sos <- ppcor::pcor.test(
    sig_var_dt$sos_iav_slp, 
    sig_var_dt$gsl_iav_slp, 
    sig_var_dt$eos_iav_slp
)

pcor_eos <- ppcor::pcor.test(
    sig_var_dt$eos_iav_slp, 
    sig_var_dt$gsl_iav_slp, 
    sig_var_dt$sos_iav_slp
)



{ #fig: Contribution of SOS and EOS
    svglite::svglite(
        "pipe/species_pheno/var_iav_trends_contri.svg",
        height = 5, pointsize = 17
    )
    
    divcols <- c("seagreen", "orange")
    
    hist(contri_dt[, c_sos - c_eos], breaks = 50, border = NA)
    abline(v = 0)

    par(
        mfrow = c(1, 2), mgp = c(1.5, 0.1, 0), mar = c(3, 3, 1, 1),
        tck = -0.01, cex.lab = 1.2
    )

    # Beta sensitivity
    hist(
        contri_dt[diff >= 0, diff],
        breaks = 50,
        col = divcols[1],
        xlim = c(-500, 500),
        ylim = c(0, 180),
        xlab = "",
        main = ""
    )
    mtext(
        expression(TC[SOS]^{Delta*IAV} ~ "-" ~ TC[EOS]^{Delta*IAV}), 1,
        line = 2, cex = 1.2
    )
    box(bty = "L")
    hist(
        contri_dt[diff < 0, diff],
        breaks = 50,
        col = divcols[2],
        add = TRUE
    )
    legend(
        grconvertX(0, "npc"), grconvertY(1, "npc"),
        fill = divcols,
        legend = c(
            expression(TC[SOS]^{Delta*IAV} ~ ">" ~ TC[EOS]^{Delta*IAV} ~ "(42%)"),
            expression(TC[SOS]^{Delta*IAV} ~ "<" ~ TC[EOS]^{Delta*IAV} ~ "(58%)")
        ),
        bty = "n",
        y.intersp = 1.2
    )


    # Partial correlation coefficients
    coord <- barplot(
        height = c(pcor_sos$estimate, pcor_eos$estimate), 
        width = 0.4, space = 1.2,
        col = divcols,
        ylim = c(0, 1), border = NA,
        xlim = c(0, 2),
        
        ylab = "Partial Correlation Coef."
    )
    text(
        coord, rep(-0.07, 2),
        labels = c(
            expression(Delta*IAV[SOS]~"~"~Delta*IAV[GSL]), 
            expression(Delta*IAV[EOS]~"~"~Delta*IAV[GSL])
        ),
        xpd = NA,
        srt = 20
    )
    text(
        coord, c(pcor_sos$estimate, pcor_eos$estimate),
        labels = c(
            paste(round(pcor_sos$estimate, 2), "**"),
            paste(round(pcor_eos$estimate, 2), "**")
        ),
        pos = 3
    )

    dev.off()
}




# ~ Species specific ####
# ~ ----------------------------------------------------------------------------
# Get species names
latin_names <- fread("data/common10.csv")
latin_names[, latin := paste(GENUS, SPECIES_NAME)]

specnames <- unique(var_dt$spec) %>%
    sort()
specnames_latin <- latin_names[gesp %in% specnames, latin]


{ # fig: Var(IAV) trends per species
    svglite::svglite(
        "pipe/var_iav_trends_per_spec.svg",
        height = 11, width = 10
    )
    par(
        mfrow = c(10, 3), bty = "L", mgp = c(1.7, 0.5, 0), mar = c(1, 3, 1, 1),
        cex.lab = 1.5, tck = -0.03, oma = c(4, 0, 0, 0)
    )

    for (sp in specnames) {
        sp_dt <- var_dt[spec == sp]

        # Percent of ascending and descending
        sos_iav_asc <- sp_dt[
            sos_iav_pval < 0.05 & sos_iav_slp > 0, 
            .N / nrow(sp_dt) * 100
        ]
        sos_iav_desc <- sp_dt[
            sos_iav_pval < 0.05 & sos_iav_slp < 0, 
            .N / nrow(sp_dt) * 100
        ]

        eos_iav_asc <- sp_dt[
            eos_iav_pval < 0.05 & eos_iav_slp > 0, 
            .N / nrow(sp_dt) * 100
        ]
        eos_iav_desc <- sp_dt[
            eos_iav_pval < 0.05 & eos_iav_slp < 0, 
            .N / nrow(sp_dt) * 100
        ]

        gsl_iav_asc <- sp_dt[
            gsl_iav_pval < 0.05 & gsl_iav_slp > 0, 
            .N / nrow(sp_dt) * 100
        ]
        gsl_iav_desc <- sp_dt[
            gsl_iav_pval < 0.05 & gsl_iav_slp < 0, 
            .N / nrow(sp_dt) * 100
        ]

        # SOS
        hist(
            sp_dt[abs(sos_iav_slp) < 5, sos_iav_slp],
            breaks = seq(-5, 5, 0.1),
            border = NA,
            xlim = c(-5, 5),
            xlab = "",
            main = ""
        )
        hist(
            sp_dt[sos_iav_pval < 0.05 & abs(sos_iav_slp) < 5, sos_iav_slp],
            breaks = seq(-5, 5, 0.1),
            include.lowest = FALSE, right = FALSE,
            border = NA, col = "seagreen",
            add = TRUE
        )
        box(bty = "L")
        abline(v = 0, lty = 2)
        legend(
            "topleft",
            bty = "n", legend = latin_names[gesp == sp, latin], text.font = 3,
            cex = 1.2
        )
        legend(
            "topright",
            bty = "n",
            legend = c(
                paste0("INC.: ", round(sos_iav_asc), "%"),
                paste0("DEC.: ", round(sos_iav_desc), "%")
            ),
            cex = 1.2
        )

        # EOS
        hist(
            sp_dt[abs(eos_iav_slp) < 5, eos_iav_slp],
            breaks = seq(-5, 5, 0.1),
            border = NA,
            xlim = c(-5, 5),
            xlab = "",
            main = ""
        )
        hist(
            sp_dt[eos_iav_pval < 0.05 & abs(eos_iav_slp) < 5, eos_iav_slp],
            breaks = seq(-5, 5, 0.1),
            include.lowest = FALSE, right = FALSE,
            border = NA, col = "orange",
            add = TRUE
        )
        box(bty = "L")
        abline(v = 0, lty = 2)
        legend(
            "topright",
            bty = "n",
            legend = c(
                paste0("INC.: ", round(eos_iav_asc), "%"),
                paste0("DEC.: ", round(eos_iav_desc), "%")
            ),
            cex = 1.2
        )

        # GSL
        hist(
            sp_dt[abs(gsl_iav_slp) < 5, gsl_iav_slp],
            breaks = seq(-5, 5, 0.1),
            border = NA,
            xlim = c(-5, 5),
            xlab = "Trend of IAV(GSL)",
            main = ""
        )
        hist(
            sp_dt[gsl_iav_pval < 0.05 & abs(gsl_iav_slp) < 5, gsl_iav_slp],
            breaks = seq(-5, 5, 0.1),
            include.lowest = FALSE, right = FALSE,
            border = NA, col = "blue",
            add = TRUE
        )
        box(bty = "L")
        abline(v = 0, lty = 2)
        legend(
            "topright",
            bty = "n",
            legend = c(
                paste0("INC.: ", round(gsl_iav_asc), "%"),
                paste0("DEC.: ", round(gsl_iav_desc), "%")
            ),
            cex = 1.2
        )
    }

    text(
        grconvertX(c(0.17, 0.51, 0.85), "ndc"),
        grconvertY(rep(0.025, 3), "ndc"),
        labels = c(
            "ΔIAV SOS", "ΔIAV EOS", "ΔIAV GSL"
        ),
        xpd = NA,
        cex = 1.5
    )

    dev.off()
}




# ~ Map ####
# ~ ----------------------------------------------------------------------------
r9_shp <- vect(base$region9_shpfile)

# Get lat and lon
fia_site_loc <- fread("data/fia_blsp_fit.csv")
fia_site_loc <- unique(fia_site_loc[, .(ID, x, y)])
var_dt <- merge(var_dt, fia_site_loc, by = "ID")
var_dt_vec <- vect(var_dt, geom = c("x", "y"), crs = "epsg:4326")


# For each species, make a map
plts <- list()
for (i in seq_along(sort(specnames))) {
    sp <- sort(specnames)[i]
    spec_vect <- var_dt_vec[var_dt_vec$spec == sp, ]
    
    p_sos <- tm_shape(r9_shp) +
        tm_polygons(col = "grey") +
        tm_shape(sf::st_as_sf(spec_vect)) +
        tm_dots(
            fill = "sos_iav_slp",
            fill.scale = tm_scale_continuous(
                n = 7, midpoint = 0, values = "pu_gn_div",
                limits = c(
                    quantile(spec_vect$sos_iav_slp, 0.025),
                    quantile(spec_vect$sos_iav_slp, 0.925)
                ),
                outliers.trunc = c(TRUE, TRUE)
            ),
            fill.legend = tm_legend(
                title = "",
                frame = FALSE,
                frame.lwd = NA,
                orientation = "landscape",
                position = tm_pos_out("center", "bottom", "center", "bottom"),
                margins = c(0, 0, 0, 0),
                title.size = 1,
                title.align = "center",
                text.size = 1,
                width = 12,
                height = 1.2,
                ticks.disable.na = TRUE,
                na.show = FALSE
            )
        ) +
        tm_credits(
            text = paste(sp, "SOS"), size = 1,
            position = tm_pos_in("center", "top", "center", "top")
        ) +
        tm_layout(
            frame = FALSE,
            component.autoscale = FALSE,
            legend.frame = FALSE
        )
    
    p_eos <- tm_shape(r9_shp) +
        tm_polygons(col = "grey") +
        tm_shape(sf::st_as_sf(spec_vect)) +
        tm_dots(
            fill = "eos_iav_slp",
            fill.scale = tm_scale_continuous(
                n = 7, midpoint = 0, values = "pu_gn_div",
                limits = c(
                    quantile(spec_vect$eos_iav_slp, 0.025),
                    quantile(spec_vect$eos_iav_slp, 0.925)
                ),
                outliers.trunc = c(TRUE, TRUE)
            ),
            fill.legend = tm_legend(
                title = "",
                frame = FALSE,
                frame.lwd = NA,
                orientation = "landscape",
                position = tm_pos_out("center", "bottom", "center", "bottom"),
                margins = c(0, 0, 0, 0),
                title.size = 1,
                title.align = "center",
                text.size = 1,
                width = 12,
                height = 1.2,
                ticks.disable.na = TRUE,
                na.show = FALSE
            )
        ) +
        tm_credits(
            text = paste(sp, "EOS"), size = 1,
            position = tm_pos_in("center", "top", "center", "top")
        ) +
        tm_layout(
            frame = FALSE,
            component.autoscale = FALSE,
            legend.frame = FALSE
        )
    
    p_gsl <- tm_shape(r9_shp) +
        tm_polygons(col = "grey") +
        tm_shape(sf::st_as_sf(spec_vect)) +
        tm_dots(
            fill = "gsl_iav_slp",
            fill.scale = tm_scale_continuous(
                n = 7, midpoint = 0, values = "pu_gn_div",
                limits = c(
                    quantile(spec_vect$gsl_iav_slp, 0.025),
                    quantile(spec_vect$gsl_iav_slp, 0.925)
                ),
                outliers.trunc = c(TRUE, TRUE)
            ),
            fill.legend = tm_legend(
                title = "",
                frame = FALSE,
                frame.lwd = NA,
                orientation = "landscape",
                position = tm_pos_out("center", "bottom", "center", "bottom"),
                margins = c(0, 0, 0, 0),
                title.size = 1,
                title.align = "center",
                text.size = 1,
                width = 12,
                height = 1.2,
                ticks.disable.na = TRUE,
                na.show = FALSE
            )
        ) +
        tm_credits(
            text = paste(sp, "GSL"), size = 1,
            position = tm_pos_in("center", "top", "center", "top")
        ) +
        tm_layout(
            frame = FALSE,
            component.autoscale = FALSE,
            legend.frame = FALSE
        )
    
    plts <- append(plts, list(p_sos))
    plts <- append(plts, list(p_eos))
    plts <- append(plts, list(p_gsl))
}

png(
    "pipe/var_iav_trend_map_spec.png",
    width = 2600, height = 700 * 10,
    res = 300
)
pp <- tmap_arrange(plts, ncol = 3, nrow = 10, outer.margins = 0.1)
print(pp)
dev.off()



{
    # fig: out/fia_dom_spec_sites.png
    png(
        "pipe/var_iav_trend_map.png", 
        width = 2000, height = 1300, res = 300
    )
    par(mar = c(3, 3, 1, 1))

    plot(
        r9_shp, col = "grey", 
        lty = 2, xlab = "Longitude", ylab = "Latitude", box = FALSE
    )
    sos_sig_idx <- which(var_dt_vec$sos_iav_pval < 0.05)
    points(
        var_dt_vec,
        pch = 16, col = adjustcolor(pt_colors, 0.5),
        cex = 0.5, lwd = 0.5
    )
    points(
        var_dt_vec[sos_sig_idx,],
        pch = 16, col = pt_colors[sos_sig_idx],
        cex = 0.5, lwd = 0.5
    )
    
    dev.off()
}






# ! How much proportion of IAV trends are significant? and
# ! Are they increasing or decreasing?

specN <- var_dt[, .(specN = .N), by = spec]
var_dt <- var_dt[specN, on = "spec"]

blsp_var_dt <- lapply(specnames, function(sp) {
    spec_dt <- var_dt[spec == sp]

    sos_iav_asc_pct <- spec_dt[
        sos_iav_pval < 0.05 & sos_iav_slp > 0,
        .N / nrow(spec_dt) * 100
    ]
    sos_iav_desc_pct <- spec_dt[
        sos_iav_pval < 0.05 & sos_iav_slp < 0,
        .N / nrow(spec_dt) * 100
    ]

    eos_iav_asc_pct <- spec_dt[
        eos_iav_pval < 0.05 & eos_iav_slp > 0,
        .N / nrow(spec_dt) * 100
    ]
    eos_iav_desc_pct <- spec_dt[
        eos_iav_pval < 0.05 & eos_iav_slp < 0,
        .N / nrow(spec_dt) * 100
    ]

    gsl_iav_asc_pct <- spec_dt[
        gsl_iav_pval < 0.05 & gsl_iav_slp > 0,
        .N / nrow(spec_dt) * 100
    ]
    gsl_iav_desc_pct <- spec_dt[
        gsl_iav_pval < 0.05 & gsl_iav_slp < 0,
        .N / nrow(spec_dt) * 100
    ]


    return(cbind(
        sp,
        sos_iav_asc_pct, sos_iav_desc_pct,
        eos_iav_asc_pct, eos_iav_desc_pct,
        gsl_iav_asc_pct, gsl_iav_desc_pct
    ))
}) %>%
    do.call(rbind, .) %>%
    as.data.table()


# out:
fwrite(blsp_var_dt, "pipe/species_pheno/blsp_var_dt.csv")


blsp_var_dt <- fread("pipe/species_pheno/blsp_var_dt.csv")

blsp_var_dt[, .(
    sos = as.numeric(sos_iav_asc_pct) + as.numeric(sos_iav_desc_pct),
    eos = as.numeric(eos_iav_asc_pct) + as.numeric(eos_iav_desc_pct)
)]
