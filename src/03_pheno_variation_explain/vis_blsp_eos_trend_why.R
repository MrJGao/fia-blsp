# ******************************************************************************
# Why SOS has no trend but EOS does given EOS is dominantely controlled by SOS?
# 
# Author: Xiaojie Gao
# Date: 2025-05-26
# ******************************************************************************
rm(list = ls())
source("src/base.R")
source("src/species_phenology/03_pheno_variation_explain/hlp_fit_bayes.R")
library(data.table)
library(lubridate)
library(magrittr)



pipedir <- "pipe/species_pheno"
dir.create(pipedir, showWarnings = FALSE, recursive = TRUE)


data_list <- LoadAnoData()
blsp_ano <- data_list$blsp_ano
dm_dt <- data_list$dm_dt
fia_attr_mn_dt <- data_list$fia_attr_mn_dt

atm_dt_file <- file.path(pipedir, "atm_dt4bayes.csv")
atm_dt <- fread(atm_dt_file)

# Also read in spring dataset to get spring temperatures
spr_dt <- fread(file.path(pipedir, "spr_dt4bayes.csv"))
atm_dt <- merge(
    atm_dt, 
    spr_dt[, .(ID, Year, spec, spr_Tavg_scl.V1, spr_Prcp_scl.V1, spr_Srad_scl.V1)], 
    by = c("ID", "Year", "spec")
)



# ~ Select sites where SOS has no trend but EOS does ####
# ~ ----------------------------------------------------------------------------
# The input file was produced by `src/species_phenology/mod_fia_blsp_trend.R`
site_trend_dt_ols <- fread("pipe/fia_dom_blsp_trend_ols.csv")[Type == "deci"]

tgt_sites <- site_trend_dt_ols[spr_pval > 0.05 & atm_pval < 0.05, ID]
eos_trend_sites <- site_trend_dt_ols[atm_pval < 0.05, .N]
#^^ 525 sites from xxx sites that EOS had a trend


#! EOS sig. trend all
eos_id <- site_trend_dt_ols[
    spr_pval > 0.05 & atm_pval < 0.05,
    ID
]
eos_dt <- atm_dt[ID %in% eos_id, ]


eos_mn_dt <- eos_dt[,
    lapply(.SD, mean),
    .SDcols = c(
        "EOS_scl.V1", "SOS_scl.V1", "smr_Tavg_scl.V1", "smr_Prcp_scl.V1",
        "atm_Tavg_scl.V1", "atm_Prcp_scl.V1", "atm_Srad_scl.V1",
        "spr_Tavg_scl.V1", "spr_Prcp_scl.V1", "spr_Srad_scl.V1"
    ),
    keyby = "Year"
]
# Rename
setnames(
    eos_mn_dt,
    colnames(eos_mn_dt[, -1]),
    gsub("_scl.V1", "", colnames(eos_mn_dt[, -1]))
)


#! EOS increasing; 319 sites
eos_asc_id <- site_trend_dt_ols[
    spr_pval > 0.05 & atm_pval < 0.05 & atm_slp > 0, 
    ID
]
eos_asc_dt <- atm_dt[ID %in% eos_asc_id, ]
eos_asc_mn_dt <- eos_asc_dt[,
    lapply(.SD, mean),
    .SDcols = c(
        "EOS_scl.V1", "SOS_scl.V1", "smr_Tavg_scl.V1", "smr_Prcp_scl.V1",
        "atm_Tavg_scl.V1", "atm_Prcp_scl.V1", "atm_Srad_scl.V1",
        "spr_Tavg_scl.V1", "spr_Prcp_scl.V1", "spr_Srad_scl.V1"
    ),
    keyby = "Year"
]
# Rename
setnames(
    eos_asc_mn_dt, 
    colnames(eos_asc_mn_dt[, -1]), 
    gsub("_scl.V1", "", colnames(eos_asc_mn_dt[, -1]))
)


#! EOS decreasing: 206 sites
eos_desc_id <- site_trend_dt_ols[
    spr_pval > 0.05 & atm_pval < 0.05 & atm_slp < 0,
    ID
]
eos_desc_dt <- atm_dt[ID %in% eos_desc_id, ]
eos_desc_mn_dt <- eos_desc_dt[,
    lapply(.SD, mean),
    .SDcols = c(
        "EOS_scl.V1", "SOS_scl.V1", "smr_Tavg_scl.V1", "smr_Prcp_scl.V1",
        "atm_Tavg_scl.V1", "atm_Prcp_scl.V1", "atm_Srad_scl.V1",
        "spr_Tavg_scl.V1", "spr_Prcp_scl.V1", "spr_Srad_scl.V1"
    ),
    keyby = "Year"
]
# Rename
setnames(
    eos_desc_mn_dt,
    colnames(eos_desc_mn_dt[, -1]),
    gsub("_scl.V1", "", colnames(eos_desc_mn_dt[, -1]))
)



# ~ Fit models ####
# ~ ----------------------------------------------------------------------------
# Fit a model using the averaged EOS variables, then predict increasing and
# decreasing mean records separately

mod <- lm(
    EOS ~ SOS + smr_Tavg + smr_Prcp + atm_Tavg + atm_Prcp + atm_Srad +
        smr_Prcp * atm_Tavg,
    data = eos_mn_dt
)
summary(mod)


# Predict EOS increasing
mod_asc <- update(mod, data = eos_asc_mn_dt)
summary(mod_asc)
mod_pred_eos_asc <- predict(mod_asc, eos_asc_mn_dt, interval = "confidence")
summary(lm(mod_pred_eos_asc[,1] ~ eos_asc_mn_dt$Year))


# Predict EOS decreasing
mod_desc <- update(mod, data = eos_desc_mn_dt)
summary(mod_desc)
mod_pred_eos_desc <- predict(mod_desc, eos_desc_mn_dt, interval = "confidence")
summary(lm(mod_pred_eos_desc[, 1] ~ eos_desc_mn_dt$Year))




CalPartialR2 <- function(themod) {
    anova_table <- car::Anova(themod, type = "II")
    ss_total <- sum(anova_table$"Sum Sq")

    partial_r2 <- anova_table$"Sum Sq" / ss_total
    names(partial_r2) <- rownames(anova_table)

    # Show partial R²s
    partial_r2_pct <- partial_r2 / sum(partial_r2)

    partial_r2_smry <- c(
        partial_r2_pct[1],
        climate = sum(partial_r2_pct[2:7]),
        partial_r2_pct[8]
    )
    return(partial_r2_smry)
}


(partial_r2_asc <- CalPartialR2(mod_asc))
(partial_r2_desc <- CalPartialR2(mod_desc))



PlotEff <- function(thedt, eff, partial_r2, ylab, textloc) {
    plot(
        thedt$Year, eff[, 1],
        type = "o", ylim = c(-2, 2), pch = 16,
        xlab = "Year", ylab = ylab
    )
    polygon(
        c(thedt[, Year], rev(thedt[, Year])),
        c(eff[, 2], rev(eff[, 3])),
        col = adjustcolor("grey", 0.5), border = NA
    )
    lmfit <- base$FitLm(thedt$Year, eff[, 1])
    abline(lmfit$fit, col = "blue")
    legend(
        textloc,
        legend = c(
            bquote(
                Partial ~ R^2~":"~ .(r2), 
                list(r2 = round(unname(partial_r2), 2))
            ),
            ifelse(lmfit$pval < 0.05,
                "Trend p-value < 0.05",
                paste("trend p-value: ", round(lmfit$pval, 3))
            )
        ),
        bty = "n"
    )
}



{ #fig: Make a figure
    svglite::svglite(
        file.path(pipedir, "eos_trend_interpret.svg"),
        width = 5, height = 5
    )

    par(mgp = c(1.5, 0.5, 0), bty = "L", mar = c(3, 3, 1, 1))
    layout(matrix(1:6, nrow = 3))

    
    # ~ EOS increassing ----------------------------------------
    # SOS effect
    partial_resid_sos <- predict(
        mod_asc,
        newdata = eos_asc_mn_dt[, .(
            Year, SOS,
            smr_Tavg = mean(smr_Tavg), smr_Prcp = mean(smr_Prcp),
            atm_Tavg = mean(atm_Tavg), atm_Prcp = mean(atm_Prcp),
            atm_Srad = mean(atm_Srad)
        )],
        interval = "confidence"
    )
    PlotEff(
        eos_asc_mn_dt, partial_resid_sos, 
        partial_r2_asc[1],
        ylab = "SOS effect", textloc = "bottomright"
    )

    # Climate effects
    partial_resid_clim <- predict(
        mod_asc,
        newdata = eos_asc_mn_dt[, .(
            Year,
            SOS = mean(SOS),
            smr_Tavg, smr_Prcp,
            atm_Tavg, atm_Prcp,
            atm_Srad
        )],
        interval = "confidence"
    )
    PlotEff(
        eos_asc_mn_dt, partial_resid_clim, 
        partial_r2_asc[2],
        ylab = "Climate Effects", textloc = "bottomright"
    )

    # Combined predictions
    PlotEff(
        eos_asc_mn_dt, mod_pred_eos_asc, 
        partial_r2_asc[1] + partial_r2_asc[2],
        ylab = "Combined Effects", textloc = "bottomright"
    )
    lines(eos_asc_mn_dt[, .(Year, EOS)], type = "l", col = "orange", lwd = 2)
    # legend(
    #     "topleft", bty = "n",
    #     legend = bquote(R^2==.(round(summary(mod_asc)$r.squared, 2)))
    # )


    # ~ EOS decreasing ----------------------------------------
    # SOS effect
    partial_resid_sos <- predict(
        mod_desc,
        newdata = eos_desc_mn_dt[, .(
            Year, SOS,
            smr_Tavg = mean(smr_Tavg), smr_Prcp = mean(smr_Prcp),
            atm_Tavg = mean(atm_Tavg), atm_Prcp = mean(atm_Prcp),
            atm_Srad = mean(atm_Srad)
        )],
        interval = "confidence"
    )
    PlotEff(
        eos_desc_mn_dt, partial_resid_sos, 
        partial_r2_desc[1],
        ylab = "SOS effect", textloc = "bottomleft"
    )


    # Climate effects
    partial_resid_clim <- predict(
        mod_desc,
        newdata = eos_desc_mn_dt[, .(
            Year,
            SOS = mean(SOS),
            smr_Tavg, smr_Prcp,
            atm_Tavg, atm_Prcp,
            atm_Srad
        )],
        interval = "confidence"
    )
    PlotEff(
        eos_desc_mn_dt, partial_resid_clim, 
        partial_r2_desc[2],
        ylab = "Climate Effects", textloc = "bottomleft"
    )


    # Combined predictions
    PlotEff(
        eos_desc_mn_dt, mod_pred_eos_desc, 
        partial_r2_desc[1] + partial_r2_desc[2],
        ylab = "Combined Effects", textloc = "bottomleft"
    )
    lines(eos_desc_mn_dt[, .(Year, EOS)], type = "l", col = "orange", lwd = 2)
    # legend(
    #     "topright", bty = "n",
    #     legend = bquote(R^2==.(round(summary(mod_desc)$r.squared, 2)))
    # )

    dev.off()
}



# ~ Replace SOS with spring temperature ####
# ~ ----------------------------------------------------------------------------
mod2 <- lm(
    EOS ~ spr_Tavg + smr_Tavg + smr_Prcp + atm_Tavg + atm_Prcp + atm_Srad +
        smr_Prcp * atm_Tavg,
    data = eos_mn_dt
)
summary(mod2)


# Predict EOS increasing
mod_asc2 <- update(mod2, data = eos_asc_mn_dt)
summary(mod_asc2)
mod_pred_eos_asc2 <- predict(mod_asc2, eos_asc_mn_dt, interval = "confidence")
summary(lm(mod_pred_eos_asc2[, 1] ~ eos_asc_mn_dt$Year))


# Predict EOS decreasing
mod_desc2 <- update(mod2, data = eos_desc_mn_dt)
summary(mod_desc2)
mod_pred_eos_desc2 <- predict(mod_desc2, eos_desc_mn_dt, interval = "confidence")
summary(lm(mod_pred_eos_desc2[, 1] ~ eos_desc_mn_dt$Year))



{ # fig: Make a figure
    svglite::svglite(
        file.path(pipedir, "eos_trend_interpret_sprTavg.svg"),
        width = 5, height = 5
    )

    par(mgp = c(1.5, 0.5, 0), bty = "L", mar = c(3, 3, 1, 1))
    layout(matrix(1:6, nrow = 3))


    # ~ EOS increassing ----------------------------------------
    # spr_Tavg effect
    partial_resid_sprTavg <- predict(
        mod_asc2,
        newdata = eos_asc_mn_dt[, .(
            Year, spr_Tavg,
            smr_Tavg = mean(smr_Tavg), smr_Prcp = mean(smr_Prcp),
            atm_Tavg = mean(atm_Tavg), atm_Prcp = mean(atm_Prcp),
            atm_Srad = mean(atm_Srad)
        )],
        interval = "confidence"
    )
    PlotEff(
        eos_asc_mn_dt, partial_resid_sprTavg, 
        ylab = "spr_Tavg effect", textloc = "bottomright"
    )

    # Climate effects
    partial_resid_clim <- predict(
        mod_asc2,
        newdata = eos_asc_mn_dt[, .(
            Year,
            spr_Tavg = mean(spr_Tavg),
            smr_Tavg, smr_Prcp,
            atm_Tavg, atm_Prcp,
            atm_Srad
        )],
        interval = "confidence"
    )
    PlotEff(
        eos_asc_mn_dt, partial_resid_clim, 
        ylab = "Climate Effects", textloc = "bottomright"
    )

    # Combined predictions
    PlotEff(
        eos_asc_mn_dt, mod_pred_eos_asc2, 
        ylab = "Combined Effects", textloc = "bottomright"
    )
    lines(eos_asc_mn_dt[, .(Year, EOS)], type = "l", col = "orange", lwd = 2)
    legend(
        "topleft", bty = "n",
        legend = bquote(R^2 == .(round(summary(mod_asc2)$r.squared, 2)))
    )


    # ~ EOS decreasing ----------------------------------------
    # spr_Tavg effect
    partial_resid_sprTavg <- predict(
        mod_desc2,
        newdata = eos_desc_mn_dt[, .(
            Year, spr_Tavg,
            smr_Tavg = mean(smr_Tavg), smr_Prcp = mean(smr_Prcp),
            atm_Tavg = mean(atm_Tavg), atm_Prcp = mean(atm_Prcp),
            atm_Srad = mean(atm_Srad)
        )],
        interval = "confidence"
    )
    PlotEff(
        eos_desc_mn_dt, partial_resid_sprTavg, 
        ylab = "spr_Tavg effect", textloc = "bottomleft"
    )


    # Climate effects
    partial_resid_clim <- predict(
        mod_desc2,
        newdata = eos_desc_mn_dt[, .(
            Year,
            spr_Tavg = mean(spr_Tavg),
            smr_Tavg, smr_Prcp,
            atm_Tavg, atm_Prcp,
            atm_Srad
        )],
        interval = "confidence"
    )
    PlotEff(
        eos_desc_mn_dt, partial_resid_clim, 
        ylab = "Climate Effects", textloc = "bottomleft"
    )


    # Combined predictions
    PlotEff(
        eos_desc_mn_dt, mod_pred_eos_desc2, 
        ylab = "Combined Effects", textloc = "bottomleft"
    )
    lines(eos_desc_mn_dt[, .(Year, EOS)], type = "l", col = "orange", lwd = 2)
    legend(
        "topright", bty = "n",
        legend = bquote(R^2 == .(round(summary(mod_desc2)$r.squared, 2)))
    )


    dev.off()
}

