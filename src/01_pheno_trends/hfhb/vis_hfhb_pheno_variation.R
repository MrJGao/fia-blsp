# ******************************************************************************
# Temporal variation of species-specific phenology at HF and HB.
# 
# The spatial coverage of trees at HF and HB are relatively small, especially
# for HF, so we can assume the trees experience similar climates. Then, their
# phenology variation is more species-specific response.
# 
# Author: Xiaojie Gao
# Date: 2025-05-30
# ******************************************************************************
rm(list = ls())
library(data.table)
library(magrittr)
library(RColorBrewer)


# ! It's hard to interpret the results here! 

# ~ HF ####
# ~ ----------------------------------------------------------------------------
spr_dt <- fread("https://harvardforest.fas.harvard.edu/data/p00/hf003/hf003-05-spring-mean-ind.csv")
atm_dt <- fread("https://harvardforest.fas.harvard.edu/data/p00/hf003/hf003-07-fall-mean-ind.csv")
com_dt <- merge(spr_dt, atm_dt, by = c("year", "tree.id", "species"))

# Select only trees w/ longer than 10 years of data
spr_tree_records <- spr_dt[, .(num_records = .N), by = tree.id]
spr_dt <- spr_dt[tree.id %in% spr_tree_records[num_records > 10, tree.id]]

atm_tree_records <- atm_dt[, .(num_records = .N), by = tree.id]
atm_dt <- atm_dt[tree.id %in% atm_tree_records[num_records > 10, tree.id]]

com_dt <- merge(spr_dt, atm_dt, by = c("year", "tree.id", "species"))
# Remove the trees that only have fall phenology
com_dt <- na.omit(com_dt[, .(year, tree.id, species, l75.doy, lf.doy)])
com_dt[, gsl := lf.doy - l75.doy]

# Select only deciduous tree species
deci_spec <- c(
    "ACPE", "ACRU", "ACSA", "AMSP", "BEAL", "BELE", "BEPA", "BEPO", "CADE",
    "FAGR", "FRAM", "NYSY", "POTR", "PRSE",
    "QUAL", "QURU", "QUVE"
)
com_dt <- com_dt[species %in% deci_spec, ]


uniqueN(com_dt$species)
uniqueN(com_dt$tree.id)


# Compute median values for the same tree
mn_dt <- com_dt[, .(
    mn_l75.doy = median(l75.doy, na.rm = TRUE),
    var_l75.doy = var(l75.doy, na.rm = TRUE),
    mn_lf.doy = median(lf.doy, na.rm = TRUE),
    var_lf.doy = var(lf.doy, na.rm = TRUE),
    mn_gsl = median(gsl, na.rm = TRUE),
    var_gsl = var(gsl, na.rm = TRUE)
), keyby = c("species", "tree.id")]


par(mfrow = c(2, 1))
boxplot(
    mn_l75.doy ~ species,
    data = mn_dt,
    xlab = "", ylab = "leaf-out",
    las = 2
)
boxplot(
    mn_lf.doy ~ species,
    data = mn_dt,
    xlab = "", ylab = "leaf-fall",
    las = 2
)



# Calculate anomalies for each year
yr_mn <- com_dt[, .(
    yr_mn_lo = mean(l75.doy, na.rm = TRUE),
    yr_mn_lf = mean(lf.doy, na.rm = TRUE)
), keyby = "year"]
com_dt <- merge(com_dt, yr_mn, by = "year")
com_dt[, leafout_ano := l75.doy - yr_mn_lo]
com_dt[, lf_ano := lf.doy - yr_mn_lf]


pdf("zzz.pdf")
par(mfrow = c(3, 1), mar = c(3, 3, 1, 1), mgp = c(1.5, 0.5, 0))
boxplot(
    leafout_ano ~ species,
    data = com_dt, 
    xlab = "", las = 2
)
abline(h = 0, lty = 2)
boxplot(
    lf_ano ~ species,
    data = com_dt, 
    xlab = "", las = 2
)
abline(h = 0, lty = 2)
null <- lapply(unique(com_dt$year), function(yr) {
    boxplot(
        leafout_ano ~ species,
        data = com_dt[year == yr, ],
        xlab = "", las = 2
    )
    abline(h = 0, lty = 2)
})
dev.off()
