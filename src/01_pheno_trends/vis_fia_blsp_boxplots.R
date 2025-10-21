# ******************************************************************************
# FIA dominant species BLSP variations
# 
# Author: Xiaojie Gao
# Date: 2024-06-17
# ******************************************************************************
library(data.table)
library(magrittr)



blsp_fit_dt <- fread("data/fia_blsp_fit.csv")[, .(
    ID, spec, x, y, 
    MidGreenup, MidGreenup_unc = MidGreenup_upr - MidGreenup_lwr,
    MidGreendown, MidGreendown_unc = MidGreendown_upr - MidGreendown_lwr
)]

# Remove records with high uncertainty
blsp_fit_dt[MidGreenup_unc > 30, MidGreenup := NA]
blsp_fit_dt[MidGreendown_unc > 30, MidGreendown := NA]

# Calculate GSL
blsp_fit_dt[, GSL := MidGreendown - MidGreenup]


mn_blsp <- blsp_fit_dt[, .(
        mn_spr = mean(MidGreenup, na.rm = TRUE), 
        mn_atm = mean(MidGreendown, na.rm = TRUE),
        mn_gsl = mean(GSL, na.rm = TRUE)
    ), 
    by = c("ID", "spec")
]

specnames <- unique(mn_blsp$spec) %>%
    sort()



# fig: FIA blsp boxplots
png("out/fia_blsp_boxplots.png", width = 1300, height = 1300, res = 300)
par(mfrow = c(3, 1), mgp = c(2, 0.5, 0), bty = "L", mar = c(0, 4, 1, 1), 
    las = 1
)

# Spring
boxplot(NA, xlim = c(1, length(specnames)), ylim = c(100, 170), 
    ylab = "Spring"
)
for (i in seq_along(specnames)) {
    boxplot(mn_blsp[spec == specnames[i], mn_spr], at = i, add = TRUE)
}

# Autumn
boxplot(NA, xlim = c(1, length(specnames)), ylim = c(240, 300),
    ylab = "Fall"
)
for (i in seq_along(specnames)) {
    boxplot(mn_blsp[spec == specnames[i], mn_atm], at = i, add = TRUE)
}

par(mar = c(3, 4, 1, 1))
# GSL
boxplot(NA, xlim = c(1, length(specnames)), ylim = c(80, 200),
    ylab = "GSL"
)
for (i in seq_along(specnames)) {
    boxplot(mn_blsp[spec == specnames[i], mn_gsl], at = i, add = TRUE)
}
axis(side = 1, at = 1:length(specnames), labels = specnames, gap.axis = 0.1)

dev.off()

