# ******************************************************************************
# Compare BLSP and MSLSP at species-dominant sites
# 
# Author: Xiaojie Gao
# Date: 2024-06-09
# ******************************************************************************
rm(list = ls())
source("src/base.R")
library(terra)
library(data.table)



mslsp_dt <- fread("data/sp_60_mslsp_dt.csv")

blsp_dir <- "pipe/blsp_fit"

evi2_dt <- fread("data/fia_dom_landsat.csv")
specnames <- unique(evi2_dt$Category)

# The following step may take some minutes
blsp_dt <- lapply(
    list.files(blsp_dir, ".Rds", full.names = TRUE),
    function(ff) {
        x <- readRDS(ff)
        phenos <- lapply(x, function(thefit) {
            if (ncol(thefit$phenos) < 22) {
                return(NULL)
            }
            pheno_dt <- thefit$phenos[, .(Year, 
                MidGreenup, MidGreenup_lwr, MidGreenup_upr,
                MidGreendown, MidGreendown_lwr, MidGreendown_upr
            )]
            pheno_dt$ID <- thefit$ID

            return(pheno_dt)
        })
        phenos <- do.call(rbind, phenos)

        return(phenos)
    }
)
blsp_dt <- do.call(rbind, blsp_dt)
blsp_dt[, Year := as.integer(Year)]

# Merge w/ MSLSP
mslsp_dt$ID <- paste0("X", mslsp_dt$X)

com_dt <- merge(mslsp_dt, blsp_dt, by.x = c("ID", "yr"), by.y = c("ID", "Year"))
com_dt <- com_dt[, .(ID, yr, spec, CommonName, Genus, Species, Type, 
    avg50PCGI, avg50PCGD, 
    MidGreenup, MidGreenup_lwr, MidGreenup_upr,
    MidGreendown, MidGreendown_lwr, MidGreendown_upr
)]

# Only do the 10 major decidous species
specnames <- c(
    "ACRU", "ACSA", "FRNI", "LALA", "POTR", "QUAL", "QUPR", "QURU",
    "QUVE", "TIAM"
)


com_dt <- com_dt[spec %in% specnames, ]


# ~ Overall comparison ####
# ~ ----------------------------------------------------------------------------
# Remove large uncertainty records
com_dt <- com_dt[MidGreenup_upr - MidGreenup_lwr < 14 & between(MidGreenup, 50, 250)]
com_dt <- com_dt[MidGreendown_upr - MidGreendown_lwr < 14 & MidGreendown > 100]
com_dt <- com_dt[Type == "deci"]


# Site-specific
rmse_dt <- by(com_dt, com_dt[, ID], function(iddt) {
    spr_rmse <- sqrt(mean((iddt$avg50PCGI - iddt$MidGreenup)^2))
    atm_rmse <- sqrt(mean((iddt$avg50PCGD - iddt$MidGreendown)^2))

    arow <- data.table(
        ID = iddt$ID[1],
        spec = iddt$spec[1],
        type = iddt$Type[1],
        commonName = iddt$CommonName[1],
        spr_rmse = spr_rmse,
        atm_rmse = atm_rmse
    )
    return(arow)
}) %>%
    do.call(rbind, .)

# A short summary of within-site RMSE values
rmse_dt[, .(
    mean_spr_rmse = mean(spr_rmse), 
    sd_spr_rmse = sd(spr_rmse), 
    mean_atm_rmse = mean(atm_rmse), 
    sd_atm_rmse = sd(atm_rmse))
]


figdir <- "out/species_phenology"
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

png(
    file.path(figdir, "fia_blsp_mslsp_corr.png"), 
    width = 2000, height = 1000, res = 300
)
par(mfrow = c(1, 2), mgp = c(1.5, 0.5, 0), mar = c(3, 3, 2, 0), bty = "l")
# png("out/fia_blsp_mslsp_corr_dark.png", width = 2000, height = 1000, res = 300)
# par(
#     bg = NA, fg = "white", col.axis = "white", col.lab = "white", 
#     col.main = "white", col.sub = "white"
# )
# pal <- colorRampPalette(
#     c(
#         rgb(0, 0, 0, 0), 
#         hcl.colors(100, "plasma")
#     ),
#     alpha = TRUE
# )
pal <- colorRampPalette(
    c(
        bg = "white", "#f9f9fd", "#a3a6e1",
        rev(brewer.pal(11, "Spectral"))
    ),
    alpha = TRUE
)

# Remove large uncertainty
vis_base$smoothScatter_J2(
    com_dt$avg50PCGI, com_dt$MidGreenup,
    pal = pal, nbin = 1000,
    main = "Spring",
    xlab = "MSLSP", ylab = "BLSP", colorful = TRUE
)

# Remove large uncertainty
vis_base$smoothScatter_J2(
    com_dt$avg50PCGD, com_dt$MidGreendown,
    pal = pal,
    main = "Autumn",
    xlab = "MSLSP", ylab = "BLSP", colorful = TRUE
)

par(new = TRUE, fig = c(0.05, 0.2, 0.45, 0.9), bty = "n", mar = c(0, 3, 0, 0))
boxplot(rmse_dt$spr_rmse, outline = FALSE, yaxt = "n", ylim = c(0, 10))
axis(2, tck = -0.1, at = c(0, 5, 10), cex.axis = 0.8, las = 1)

par(new = TRUE, fig = c(0.55, 0.7, 0.45, 0.9), bty = "n", mar = c(0, 3, 0, 0))
boxplot(rmse_dt$atm_rmse, outline = FALSE, yaxt = "n", ylim = c(0, 10))
axis(2, tck = -0.1, at = c(0, 5, 10), cex.axis = 0.8, las = 1)

dev.off()





# ~ Species-specific comparison ####
# ~ ----------------------------------------------------------------------------
pdf(file.path(figdir, "fia_blsp_mslsp_corr_spec.pdf"), width = 8, height = 4)
par(mfrow = c(1, 2), mgp = c(1.5, 0.5, 0), mar = c(3, 3, 3, 1))

# To summarize species-specific correlations
sum_dt <- NULL
for (specname in specnames) {
    spec_dt <- com_dt[spec == specname,]
    # Remove large uncertainty
    spr_dt <- spec_dt[MidGreenup_upr - MidGreenup_lwr < 14 &
        between(MidGreenup, 50, 250), ]
    vis_base$smoothScatter_J2(spr_dt$avg50PCGI, spr_dt$MidGreenup,
        main = paste0("Spring\n", 
            specname, " (", spec_dt$CommonName[1], ",", spec_dt$Type[1], ")"
        ),
        xlab = "MSLSP", ylab = "BLSP", colorful = FALSE
    )
    sum_spr <- data.table(
        spec = specname, type = spec_dt$Type[1], 
        commonName = spec_dt$CommonName[1],
        season = "spring", 
        corr = cor(spr_dt$avg50PCGI, spr_dt$MidGreenup)^2,
        rmse = sqrt(mean((spr_dt$avg50PCGI - spr_dt$MidGreenup)^2))
    )

    # Remove large uncertainty
    atm_dt <- spec_dt[MidGreendown_upr - MidGreendown_lwr < 14 & MidGreendown > 100, ]
    vis_base$smoothScatter_J2(atm_dt$avg50PCGD, atm_dt$MidGreendown,
        main = paste0("Fall\n", 
            specname, " (", spec_dt$CommonName[1], ",", spec_dt$Type[1], ")"
        ),
        xlab = "MSLSP", ylab = "BLSP", colorful = FALSE
    )
    sum_atm <- data.table(
        spec = specname, type = spec_dt$Type[1],
        commonName = spec_dt$CommonName[1],
        season = "fall",
        corr = cor(atm_dt$avg50PCGD, atm_dt$MidGreendown)^2,
        rmse = sqrt(mean((atm_dt$avg50PCGD - atm_dt$MidGreendown)^2))
    )

    sum_dt <- rbind(sum_dt, sum_spr, sum_atm)
}


# Make a summary figure
par(mfrow = c(1, 2), bty = "L", mgp = c(1.5, 0.5, 0), mar = c(3, 3, 1, 1))

# ~ Correlation ----------------------------------------
# Spring
boxplot(NA, xlim = c(0.5, 2.5), ylim = c(0, 1), ylab = expression(R^2))
axis(1, at = 1:2, labels = c("Spring", "Fall"))
boxplot(
    at = 1, sum_dt[type == "deci" & season == "spring", corr],
    add = TRUE, col = "seagreen", boxwex = 0.5
)
# # Fall
boxplot(
    at = 2, sum_dt[type == "deci" & season == "fall", corr],
    add = TRUE, col = "seagreen", boxwex = 0.5
)

legend("topright", bty = "n", legend = c("deci", "ever"), 
    fill = c("seagreen", "cyan")
)

# ~ RMSE ----------------------------------------
# Spring
boxplot(NA, xlim = c(0.5, 2.5), ylim = c(0, 10), ylab = "RMSE")
axis(1, at = 1:2, labels = c("Spring", "Fall"))
boxplot(
    at = 1, sum_dt[type == "deci" & season == "spring", rmse],
    add = TRUE, col = "seagreen", boxwex = 0.5
)
# Fall
boxplot(
    at = 2, sum_dt[type == "deci" & season == "fall", rmse],
    add = TRUE, col = "seagreen", boxwex = 0.5
)

dev.off()

