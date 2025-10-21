# ******************************************************************************
# Analyze contributions of SOS and EOS on GSL in terms of trend and interannual
# variability. 
# 
# Author: Xiaojie Gao
# Date: 2025-02-23
# ******************************************************************************
rm(list = ls())
source("src/base.R")
library(data.table)
library(magrittr)
library(lubridate)
library(ggplot2)
library(ggridges)



DecomposeSignal <- function(pheno_dt) {
    # Fit a linear model to get trend
    model <- lm(pheno_dt$pheno ~ pheno_dt$yr)
    # Trend
    slp <- coef(model)[2]

    # Isolate IAV from trend and calculate variance
    trend_estimated <- predict(model)
    # var_trend <- var(trend_estimated)
    detrend <- pheno_dt$pheno - trend_estimated
    var_iav <- var(detrend)
    
    var_total <- var(pheno_dt$pheno)

    return(list(
        stat = data.table(
            slp, var_iav, var_total
        ),
        detrend = detrend
    ))
}


AnalyzeSiteSpecData <- function(specname) {
    # Read the species file
    spec_fit <- readRDS(file.path("pipe/blsp_fit", paste0(specname, ".Rds")))
    
    # Use the EVI2 dt to get the latitude for each site
    evi2_dt <- rbind(
        fread("data/fia_dom_landsat.csv"),
        fread("data/fia_dom_landsat_more.csv")
    )
    # B/c I updated the major species site table, here I need to filter out the
    # sites that were removed and add sites that were added!
    majorspec_dt <- fread("data/fia_domspec_truexy.csv")
    major_ids <- majorspec_dt[, paste0("X", X)]
    evi2_dt <- evi2_dt[ID %in% major_ids]

    yrs <- spec_fit[[1]]$phenos[, as.numeric(Year)]

    var_dt <- lapply(spec_fit, function(thefit) {
        phenos <- thefit$phenos
        if (nrow(phenos) < length(yrs) || ncol(phenos) < 22) {
            return(NULL)
        }
        setorder(phenos, Year)

        phenos[MidGreenup_upr - MidGreenup_lwr > 30, MidGreenup := NA]
        phenos[MidGreendown_upr - MidGreendown_lwr > 30, MidGreenup := NA]

        phenos <- na.omit(phenos[, .(MidGreenup, MidGreendown, yr = Year)])
        phenos[, yr := as.numeric(yr)]

        if (nrow(phenos) < 5) {
            return(NULL)
        }

        sos_de <- DecomposeSignal(phenos[, .(pheno = MidGreenup, yr)])
        eos_de <- DecomposeSignal(phenos[, .(pheno = MidGreendown, yr)])
        gsl_de <- DecomposeSignal(
            phenos[, .(pheno = MidGreendown - MidGreenup, yr)]
        )
        # Get latitude
        lat <- evi2_dt[ID == thefit$ID, Latitude][1]

        
        # Calculate trend cotribution
        sos_trend_pct <- -sos_de$stat$slp / gsl_de$stat$slp * 100
        eos_trend_pct <- eos_de$stat$slp / gsl_de$stat$slp * 100

        # Calculate IAV contribution
        rho <- cor(sos_de$detrend, eos_de$detrend)
        covar <- rho * sqrt(sos_de$stat$var_iav) * sqrt(eos_de$stat$var_iav)
        sos_iav_pct <- (sos_de$stat$var_iav - covar) / gsl_de$stat$var_iav * 100
        eos_iav_pct <- (eos_de$stat$var_iav - covar) / gsl_de$stat$var_iav * 100

        # Correlation contribution
        covar_pct <- covar / gsl_de$stat$var_iav * 100

        arow <- data.table(
            spec = specname, ID = thefit$ID, lat = lat,
            rho, covar_pct, 
            sos_slp = sos_de$stat$slp,
            sos_iav = sos_de$stat$var_iav,
            eos_slp = eos_de$stat$slp,
            eos_iav = eos_de$stat$var_iav,
            
            gsl_slp = gsl_de$stat$slp,
            gsl_iav = gsl_de$stat$var_iav,

            sos_trend_pct = sos_trend_pct, 
            eos_trend_pct = eos_trend_pct, 
            sos_iav_pct = sos_iav_pct,
            eos_iav_pct = eos_iav_pct
        )

        return(arow)
    }) %>%
        rbindlist()

    return(var_dt)
}



evi2_dt <- rbind(
    fread("data/fia_dom_landsat.csv"),
    fread("data/fia_dom_landsat_more.csv")
)
# B/c I updated the major species site table, here I need to filter out the
# sites that were removed and add sites that were added!
majorspec_dt <- fread("data/fia_domspec_truexy.csv")
major_ids <- majorspec_dt[, paste0("X", X)]
evi2_dt <- evi2_dt[ID %in% major_ids]

specnames <- majorspec_dt[Type == "deci", spec] %>%
    unique()


var_dt <- lapply(specnames, AnalyzeSiteSpecData)
var_dt <- rbindlist(var_dt)


# Remove non-sense values
var_dt <- var_dt[
    sos_iav_pct < 100 & eos_iav_pct < 100 & 
    (sos_iav_pct + eos_iav_pct - 100) < 0.00001
]

var_dt <- var_dt[
    (sos_trend_pct + eos_trend_pct - 100) < 0.001
]

var_dt[, ":=" (
    iav_diff = sos_iav_pct - eos_iav_pct,
    trend_diff = sos_trend_pct - eos_trend_pct
)]


# Get trend significance
site_trend_dt <- fread("pipe/fia_dom_blsp_trend_ols.csv")[Type == "deci"]
var_dt <- merge(var_dt, site_trend_dt[, .(ID, gsl_pval)], by = "ID")

# Get latin names
latin_names <- fread("data/common10.csv")
var_dt <- merge(var_dt, latin_names, by.x = "spec", by.y = "gesp")
var_dt[, latin := paste(GENUS, SPECIES_NAME)]

# Calculate some proportions
var_dt[iav_diff < 0, .N / nrow(var_dt)]
var_dt[iav_diff > 0, .N / nrow(var_dt)]
var_dt[gsl_pval < 0.05 & iav_diff < 0, .N / nrow(var_dt[gsl_pval < 0.05])]
var_dt[gsl_pval < 0.05 & iav_diff > 0, .N / nrow(var_dt[gsl_pval < 0.05])]

var_dt[trend_diff < 0, .N / nrow(var_dt)]
var_dt[trend_diff > 0, .N / nrow(var_dt)]
var_dt[gsl_pval < 0.05 & trend_diff < 0, .N / nrow(var_dt[gsl_pval < 0.05])]
var_dt[gsl_pval < 0.05 & trend_diff > 0, .N / nrow(var_dt[gsl_pval < 0.05])]

# boxplot(
#     var_dt[gsl_pval < 0.05, covar_pct], outline = FALSE, 
#     xlab = expression(
#         rho ~ "*" ~ sigma[SOS]^2 ~ "*" ~ sigma[EOS]^2 ~ 
#         "/" ~ sigma[GSL]^2 ~"*"~100
#     )
# )


# ~ Make SOS-EOS slope and IAV figures ####
# ~ ----------------------------------------------------------------------------

par(mfrow = c(1, 2))
plot(
    var_dt[, .(sos_slp, eos_slp)], 
    xlim = c(-2, 2), ylim = c(-2, 2),
    xlab = "SOS slope", ylab = "EOS slope"
)
points(
    var_dt[gsl_pval < 0.05, .(sos_slp, eos_slp)], 
    pch = 21, bg = adjustcolor("red", 0.5)
)
gsl_lt_0_pct <- round(
    var_dt[gsl_pval < 0.05 & gsl_slp < 0, .N / nrow(var_dt)] * 100
)
gsl_gt_0_pct <- round(
    var_dt[gsl_pval < 0.05 & gsl_slp < 0, .N / nrow(var_dt)] * 100
)
text(
    -1, 1,
    labels = paste0("Sig. longer GSL (", gsl_lt_0_pct, "%)"),
    col = "red"
)
text(
    1, -1,
    labels = paste0("Sig. shorter GSL (", gsl_gt_0_pct, "%)"),
    col = "red"
)
abline(h = 0)
abline(v = 0)
text(
    c(2, -2, -2, 2),
    c(2, 2, -2, -2),
    labels = c("QI", "QII", "QIII", "QIV")
)



plot(
    var_dt[, .(sqrt(sos_iav), sqrt(eos_iav))], 
    xlim = c(0, 30), ylim = c(0, 30),
    xlab = "SOS sd(IAV)", ylab = "EOS sd(IAV)"
)
points(
    var_dt[gsl_pval < 0.05 & sos_iav > eos_iav, .(sqrt(sos_iav), sqrt(eos_iav))],
    pch = 21, bg ="seagreen"
)
points(
    var_dt[gsl_pval < 0.05 & sos_iav < eos_iav, .(sqrt(sos_iav), sqrt(eos_iav))],
    pch = 21, bg ="orange"
)
var_dt[gsl_pval < 0.05 & sos_iav > eos_iav, .N / nrow(var_dt)]
var_dt[gsl_pval < 0.05 & sos_iav < eos_iav, .N / nrow(var_dt)]
abline(0, 1, lty = 2)



# ~ Make the figure ####
# ~ ----------------------------------------------------------------------------
outdir <- "out/"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


# divcols <- hcl.colors(2, "Berlin")
divcols <- c("seagreen", "orange")

{
    svglite::svglite(
        file.path(outdir, "gsl_decompose.svg"),
        width = 4, height = 6.5
    )
    par(
        mfrow = c(2, 1), mgp = c(1.5, 0.1, 0), mar = c(4, 3, 0.5, 0.2), 
        tck = -0.01, cex.lab = 1.2
    )

    hist(
        var_dt[iav_diff > 0 & gsl_pval < 0.05, iav_diff],
        breaks = 20,
        col = divcols[1],
        xlim = c(-100, 100),
        ylim = c(0, 40),
        xlab = "",
        main = ""
    )
    mtext(
        expression(C[SOS]^IAV~"-"~C[EOS]^IAV), 1, 
        line = 2, cex = 1.2
    )
    box(bty = "L")
    hist(
        var_dt[iav_diff < 0 & gsl_pval < 0.05, iav_diff],
        breaks = 20,
        col = divcols[2],
        add = TRUE
    )
    legend(
        grconvertX(0, "npc"), grconvertY(0.97, "npc"),
        fill = divcols,
        legend = c(
            expression(C[SOS]^IAV~">"~C[EOS]^IAV~"(55%)"), 
            expression(C[SOS]^IAV ~ "<" ~ C[EOS]^IAV ~ "(45%)")
        ),
        bty = "n", 
        y.intersp = 1.2
    )
    

    qt <- quantile(var_dt$trend_diff, c(0.1, 0.97))
    hist(
        var_dt[
            between(trend_diff, qt[1], qt[2]) & trend_diff > 0 & gsl_pval < 0.05,
            trend_diff
        ],
        breaks = 30,
        col = divcols[1],
        xlim = c(-300, 300),
        xlab = "",
        ylim = c(0, 65),
        main = ""
    )
    mtext(
        expression(C[SOS]^Trend~"-"~C[EOS]^Trend), 1,
        line = 2, cex = 1.2
    )
    box(bty = "L")
    hist(
        var_dt[
            between(trend_diff, qt[1], qt[2]) & trend_diff < 0 & gsl_pval < 0.05,
            trend_diff
        ],
        breaks = 20,
        col = divcols[2],
        add = TRUE
    )
    legend(
        grconvertX(0, "npc"), grconvertY(1, "npc"),
        fill = divcols,
        legend = c(
            expression(C[SOS]^Trend~">"~C[EOS]^Trend~"(32%)"), 
            expression(C[SOS]^Trend~"<"~C[EOS]^Trend~"(68%)")
        ),
        bty = "n",
        y.intersp = 1.2
    )
    
    dev.off()
}



# ~ SI figure ####
# ~ ----------------------------------------------------------------------------
g1 <- ggplot(var_dt, aes(x = iav_diff, y = latin)) +
    geom_density_ridges2(fill = "grey90") +
    theme_ridges(center_axis_labels = TRUE, font_size = 20) +
    xlim(-100, 100) +
    ylab("") +
    xlab(expression(C[SOS]^IAV~"-"~C[EOS]^IAV)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme(
        legend.position = "right",
        axis.text.y = element_text(face = "italic")
    )

ggsave(
    file.path(outdir, "gsl_decompose_iav_spec.png"), g1,
    width = 6, height = 8
)


g2 <- ggplot(
    var_dt[
        between(trend_diff, qt[1], qt[2]), .(trend_diff, latin)], 
    aes(x = trend_diff, y = latin)
) +
    geom_density_ridges2(fill = "grey90") +
    theme_ridges(center_axis_labels = TRUE, font_size = 20) +
    ylab("") +
    xlim(-400, 400) +
    xlab(expression(C[SOS]^Trend~"-"~C[EOS]^Trend)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme(
        legend.position = "right",
        axis.text.y = element_text(face = "italic")
    )

ggsave(
    file.path(outdir, "gsl_decompose_trend_spec.png"), g2,
    width = 6, height = 8
)
