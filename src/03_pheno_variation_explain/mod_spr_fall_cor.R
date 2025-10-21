# ******************************************************************************
# Check whether spring and fall phenology are correlated using partial
# correlation coefficients for controlling environmental factors.
# 
# Author: Xiaojie Gao
# Date: 2024-06-25
# ******************************************************************************
rm(list = ls())
source("src/base.R")
library(data.table)
library(magrittr)
library(ppcor)
library(parallel)
library(ggplot2)



# Prepare the site-specific climate data
PrepareClim <- function() {
    # Compute site-specific metrics
    # Spring variables
    spr_dm <- blsp_dm_dt[year > 1984 & month(date) %in% c(3, 4, 5), .(
        spr_tavg = mean((`tmax (deg c)` + `tmin (deg c)`) / 2),
        spr_prcp = sum(`prcp (mm/day)`),
        spr_vp = sum(`prcp (mm/day)`),
        spr_dayl = mean(`dayl (s)`),
        spr_srad = mean((`srad (W/m^2)` * `dayl (s)`) / 1e6)
    ), by = .(ID, year)]
    
    # Summer mean precipitation
    smr_dm <- blsp_dm_dt[year > 1984 & month(date) %in% c(6, 7, 8), .(
        smr_tavg = mean((`tmax (deg c)` + `tmin (deg c)`) / 2),
        smr_prcp = sum(`prcp (mm/day)`),
        smr_vp = sum(`prcp (mm/day)`),
        smr_dayl = mean(`dayl (s)`),
        smr_srad = mean((`srad (W/m^2)` * `dayl (s)`) / 1e6)
    ), by = .(ID, year)]

    # Summer-autumn mean precipitation
    stf_dm <- blsp_dm_dt[year > 1984 & month(date) %in% c(6:11), .(
        stf_tavg = mean((`tmax (deg c)` + `tmin (deg c)`) / 2),
        stf_prcp = sum(`prcp (mm/day)`),
        stf_vp = sum(`prcp (mm/day)`),
        stf_dayl = mean(`dayl (s)`),
        stf_srad = mean((`srad (W/m^2)` * `dayl (s)`) / 1e6)
    ), by = .(ID, year)]
    
    # Autumn mean temperature
    atm_dm <- blsp_dm_dt[year > 1984 & month(date) %in% c(9:11), .(
        atm_tavg = mean((`tmax (deg c)` + `tmin (deg c)`) / 2),
        atm_prcp = sum(`prcp (mm/day)`),
        atm_vp = sum(`prcp (mm/day)`),
        atm_dayl = mean(`dayl (s)`),
        atm_srad = mean((`srad (W/m^2)` * `dayl (s)`) / 1e6)
    ), by = .(ID, year)]
    
    cdd_dm <- blsp_dm_dt[year > 1984 & month(date) > 7, .(
        (`tmax (deg c)` + `tmin (deg c)`) / 2
    ), by = .(ID, year)]

    site_dm <- Reduce(
        function(x, y) {
            merge(x, y, by = c("ID", "year"))
        }, 
        list(spr_dm, smr_dm, stf_dm, atm_dm)
    )

    return(site_dm)
}

# Calculate cold-degree day
CalCDD <- function(st) {
    st_dm <- blsp_dm_dt[ID == st$ID[1] & year > 1984,]
    mn_eos <- mean(st$MidGreendown, na.rm = TRUE)

    cdd1 <- st_dm[month(date) > 7 & yday < mn_eos, .(
        cdd1 = sum(
            (25 - (`tmax (deg c)` + `tmin (deg c)`) / 2) * (`dayl (s)` / 3600 / 12)
        )
    ), by = .(ID, year)]
    
    cdd2 <- st_dm[month(date) > 7 & yday < mn_eos, .(
        cdd2 = sum(
            (25 - (`tmax (deg c)` + `tmin (deg c)`) / 2)^2 * (`dayl (s)` / 3600 / 12)
        )
    ), by = .(ID, year)]

    cdd_dt <- Reduce(
        function(x, y) {
            merge(x, y, by = c("ID", "year"))
        },
        list(cdd1, cdd2)
    )

    return(cdd_dt)
}



# For saving intermittent results
pipedir <- "pipe/species_pheno"
dir.create(pipedir, showWarnings = FALSE, recursive = TRUE)


blsp <- fread("data/fia_blsp_fit.csv")
# Remove large uncertainty observations
blsp[MidGreenup_upr - MidGreenup_lwr > 30, MidGreenup := NA]
blsp[MidGreendown_upr - MidGreendown_lwr > 30, MidGreendown := NA]



# ~ Site level ####
# ~ ----------------------------------------------------------------------------
blsp_dm_dt <- fread("data/blsp_dm_dt.csv")

# Calculate site-level means
site_dm <- PrepareClim()

# Compute spring-fall correlation. The following code may take a while to run

cl <- makeCluster(10, outfile = "")
calls <- clusterCall(cl, function() {
    suppressWarnings({
        source("src/base.R")
        library(data.table)
        library(magrittr)
        library(ppcor)
    })
})
clusterExport(cl, c("blsp_dm_dt", "site_dm", "blsp", "CalCDD"))

st_fit <- clusterApplyLB(cl, x = unique(blsp$ID), fun = function(id) {
    st <- blsp[ID == id, ]

    # Compute cold-degree days
    cdd_dt <- CalCDD(st)

    st <- merge(st, site_dm, by.x = c("ID", "Year"), by.y = c("ID", "year"))
    st <- merge(st, cdd_dt, by.x = c("ID", "Year"), by.y = c("ID", "year"))
    st <- na.omit(st)

    if (nrow(st) < 10 || var(st$MidGreenup) < 5) {
        return(NULL)
    }

    corr <- cor.test(st$MidGreenup, st$MidGreendown, method = "spearman")
    parcorr <- ppcor::pcor.test(
        st$MidGreenup, st$MidGreendown,
        st[, .(atm_tavg, cdd1)],
        method = "spearman"
    )

    return(data.table(
        ID = id,
        spec = st$spec[1], 
        corr = corr$estimate,
        corr_p = corr$p.value,
        parcorr = parcorr$estimate,
        parcorr_p = parcorr$p.value
    ))
})
st_fit <- rbindlist(st_fit)

stopCluster(cl)

# out: 
fwrite(st_fit, file.path(pipedir, "spr_fall_parcor.csv"))



# ~ Analyze the partial correlations ####
# ~ ----------------------------------------------------------------------------
st_fit <- fread(file.path(pipedir, "spr_fall_parcor.csv"))

# Get Latin name
latin_names <- fread("data/common10.csv")
st_fit <- merge(st_fit, latin_names, by.x = "spec", by.y = "gesp")
st_fit[, latin := paste(GENUS, SPECIES_NAME)]

# Make a figure that shows histograms for each species
{ #fig:
    pdf(file.path(pipedir, "spr_fall_cor.pdf"), width = 10, height = 4)
    par(mfrow = c(1, 2))
    
    for (sp in unique(st_fit$spec)) {
        hist(
            st_fit[spec == sp, parcorr], 
            breaks = 50, 
            border = "white", 
            xlab = "Partial Corr.",
            main = ""
        )
        hist(
            st_fit[spec == sp & parcorr_p < 0.05, parcorr], 
            breaks = 50, 
            border = "white", 
            col = "red", add = TRUE
        )

        hist(
            st_fit[spec == sp, corr], 
            breaks = 50, 
            border = "white", 
            xlab = "Corr.",
            main = ""
        )
        hist(
            st_fit[spec == sp & corr_p < 0.05, corr], 
            breaks = 50, 
            border = "white", 
            col = "red", 
            add = TRUE
        )

        text(
            grconvertX(0.5, "ndc"), grconvertY(0.9, "ndc"),
            labels = sp, 
            cex = 1.5,
            xpd = NA
        )

    }

    dev.off()
}
# ^^ It seems like spring and fall phenology are positively correlated whether
# using regular spearman correlation or the partial correlation methods.

# Make a figure to show the species-specific partial correlations for adding to
# the SI of the paper
{
    png(
        file.path("out/", "spr_fall_parcorr.png"),
        width = 1600, height = 2800, res = 300
    )
    par(mfrow = c(5, 2), mgp = c(1.5, 0.5, 0), mar = c(3, 3, 0.5, 1))

    for (sp in unique(st_fit$spec)) {
        breaks_vec <- quantile(st_fit[spec == sp, parcorr], c(0, 1))
        breaks_vec <- seq(breaks_vec[1], breaks_vec[2], len = 50)
        
        hist(
            st_fit[spec == sp, parcorr],
            breaks = breaks_vec,
            border = "white",
            xlab = "Partial Corr.",
            main = ""
        )
        hist(
            st_fit[spec == sp & parcorr_p < 0.05, parcorr],
            breaks = breaks_vec,
            border = "white",
            col = adjustcolor("red", 0.5), add = TRUE
        )

        num_sig <- length(st_fit[spec == sp & parcorr_p < 0.05, parcorr])
        num_tot <- length(st_fit[spec == sp, parcorr])
        sigpct <- num_sig / num_tot * 100
        
        legend(
            "topleft", bty = "n",
            legend = c(
                paste0("Sig.:", round(sigpct, 1), "%")
            ),
            fill = c(adjustcolor("red", 0.5)),
            title = st_fit[spec == sp, latin][1],
            cex = 1.2
        )
    }

    dev.off()
}



# ~ # What's the porportion of signifcance? ------------------------------------

lapply(sort(unique(st_fit$spec)), function(sp) {
    spec_dt <- st_fit[spec == sp]
    total <- nrow(spec_dt)
    pct <- nrow(spec_dt[parcorr_p < 0.05]) / total
    parcorr <- mean(spec_dt[parcorr_p < 0.05, parcorr])

    return(cbind(sp, total, pct, parcorr))
})
# ^^ The mean partial correlation are kind of conserved, ranging from 0.4-0.5.


# ~ Does FIA attributes explain the partial correlation? ####
# ~ ----------------------------------------------------------------------------
# Merge w/ FIA attributes
fia_attr_dt <- fread("data/fia_attr_dt.csv")
fia_attr_dt <- fia_attr_dt[nsubplotsForest == 4, ]

# Compute species richness
fia_spec_cols <- names(fia_attr_dt) %>%
    grep("spBAm2ha*.*_live", ., value = TRUE)

spec_rich <- lapply(fia_spec_cols, function(spec_col) {
    spec_ba <- fia_attr_dt[, get(spec_col)]
    spec_ba <- ifelse(is.na(spec_ba), 0, 1)
    return(spec_ba)
})
spec_rich <- do.call(cbind, spec_rich)
spec_rich <- apply(spec_rich, 1, sum)

fia_attr_dt[, richness := spec_rich]

st_fit <- merge(
    st_fit,
    fia_attr_dt[, .(
        ID = paste0("X", X),
        nstem = plotGroupnStemsha_live,
        ba = plotGroupBAm2ha_live,
        richness
    )],
    by = "ID"
)

par(mfrow = c(1, 3))
plot(st_fit[, .(parcorr, ba)])
plot(st_fit[, .(parcorr, nstem)])
plot(st_fit[, .(parcorr, richness)])

summary(
    lm(parcorr ~ nstem + ba + richness, data = st_fit)
)
# ^^ FIA attributes poorly explain the partical correlation values.