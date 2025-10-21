# ******************************************************************************
# Visualize site-specific trends of IAV for spring/fall/gsl for HF & HB pheno
# 
# Author: Xiaojie Gao
# Date: 2025-06-20
# ******************************************************************************
library(data.table)




hfhb_var_dt <- fread("pipe/species_pheno/hfhb_var_dt.csv")

hf_var_dt <- hfhb_var_dt[site == "HF",]
hb_var_dt <- hfhb_var_dt[site == "HB",]


# Percent of ascending and descending
hf_sos_iav_asc <- hf_var_dt[sos_iav_pval < 0.05 & sos_iav_slp > 0, .N]
hf_sos_iav_desc <- hf_var_dt[sos_iav_pval < 0.05 & sos_iav_slp < 0, .N]

hf_eos_iav_asc <- hf_var_dt[eos_iav_pval < 0.05 & eos_iav_slp > 0, .N]
hf_eos_iav_desc <- hf_var_dt[eos_iav_pval < 0.05 & eos_iav_slp < 0, .N]

hf_gsl_iav_asc <- hf_var_dt[gsl_iav_pval < 0.05 & gsl_iav_slp > 0, .N]
hf_gsl_iav_desc <- hf_var_dt[gsl_iav_pval < 0.05 & gsl_iav_slp < 0, .N]

hb_sos_iav_asc <- hb_var_dt[sos_iav_pval < 0.05 & sos_iav_slp > 0, .N]
hb_sos_iav_desc <- hb_var_dt[sos_iav_pval < 0.05 & sos_iav_slp < 0, .N]

hb_eos_iav_asc <- hb_var_dt[eos_iav_pval < 0.05 & eos_iav_slp > 0, .N]
hb_eos_iav_desc <- hb_var_dt[eos_iav_pval < 0.05 & eos_iav_slp < 0, .N]

hb_gsl_iav_asc <- hb_var_dt[gsl_iav_pval < 0.05 & gsl_iav_slp > 0, .N]
hb_gsl_iav_desc <- hb_var_dt[gsl_iav_pval < 0.05 & gsl_iav_slp < 0, .N]


PlotHist <- function(st_dt, varname, varcol, asc, desc) {
    slp_varname <- paste0(varname, "_iav_slp")
    pval_varname <- paste0(varname, "_iav_pval")

    hist(
        st_dt[abs(get(slp_varname)) < 5, get(slp_varname)],
        breaks = seq(-5, 5, 0.1),
        border = NA,
        xlim = c(-5, 5),
        xlab = paste("ΔIAV", toupper(varname)),
        main = ""
    )
    hist(
        st_dt[
            get(pval_varname) < 0.05 & abs(get(slp_varname)) < 5, 
            get(slp_varname)
        ],
        breaks = seq(-5, 5, 0.1),
        include.lowest = FALSE, right = FALSE,
        border = NA, col = varcol,
        add = TRUE
    )
    box(bty = "L")
    abline(v = 0, lty = 2)
    
    legend(
        "topright",
        bty = "n",
        legend = c(
            paste0(
                "INC.: ", round(asc / nrow(st_dt) * 100), "%",
                "(", asc, " / ", nrow(st_dt), ")"
            ),
            paste0(
                "DEC.: ", round(desc / nrow(st_dt) * 100), "%",
                "(", desc, " / ", nrow(st_dt), ")"
            )
        )
    )
}


{ # fig: 
    svglite::svglite(
        "pipe/var_iav_trends_hfhb.svg",
        height = 4, width = 8
    )
    par(
        mfrow = c(2, 3), bty = "L", mgp = c(1.7, 0.5, 0), mar = c(3, 3, 1, 1),
        cex.lab = 1.2, tck = -0.03
    )

    PlotHist(hf_var_dt, "sos", "seagreen", hf_sos_iav_asc, hf_sos_iav_desc)
    PlotHist(hf_var_dt, "eos", "orange", hf_eos_iav_asc, hf_eos_iav_desc)
    PlotHist(hf_var_dt, "gsl", "blue", hf_gsl_iav_asc, hf_gsl_iav_desc)
    
    PlotHist(hb_var_dt, "sos", "seagreen", hb_sos_iav_asc, hb_sos_iav_desc)
    PlotHist(hb_var_dt, "eos", "orange", hb_eos_iav_asc, hb_eos_iav_desc)
    PlotHist(hb_var_dt, "gsl", "blue", hb_gsl_iav_asc, hb_gsl_iav_desc)

    dev.off()
}



# ~ Relative contributions of SOS and EOS to GSL ####
# ~ ----------------------------------------------------------------------------

# HF -----------------
hf_contri_dt <- hf_var_dt[
    gsl_iav_pval < 0.05,
    .(
        c_sos = sos_iav_slp / gsl_iav_slp * 100,
        c_eos = eos_iav_slp / gsl_iav_slp * 100
    )
]

hf_contri_dt[, diff := c_sos - c_eos]

hf_contri_dt[c_sos >= c_eos, .N / nrow(hf_contri_dt) * 100]
hf_contri_dt[c_sos < c_eos, .N / nrow(hf_contri_dt) * 100]


hf_pcor_sos <- ppcor::pcor.test(
    hf_var_dt$sos_iav_slp,
    hf_var_dt$gsl_iav_slp,
    hf_var_dt$eos_iav_slp
)

hf_pcor_eos <- ppcor::pcor.test(
    hf_var_dt$eos_iav_slp,
    hf_var_dt$gsl_iav_slp,
    hf_var_dt$sos_iav_slp
)


# HB -----------------
hb_contri_dt <- hb_var_dt[
    gsl_iav_pval < 0.05,
    .(
        c_sos = sos_iav_slp / gsl_iav_slp * 100,
        c_eos = eos_iav_slp / gsl_iav_slp * 100
    )
]

hb_contri_dt[, diff := c_sos - c_eos]

hb_contri_dt[c_sos >= c_eos, .N / nrow(hb_contri_dt) * 100]
hb_contri_dt[c_sos < c_eos, .N / nrow(hb_contri_dt) * 100]


hb_sig_var_dt <- hb_var_dt

hb_pcor_sos <- ppcor::pcor.test(
    hb_sig_var_dt$sos_iav_slp,
    hb_sig_var_dt$gsl_iav_slp,
    hb_sig_var_dt$eos_iav_slp
)

hb_pcor_eos <- ppcor::pcor.test(
    hb_sig_var_dt$eos_iav_slp,
    hb_sig_var_dt$gsl_iav_slp,
    hb_sig_var_dt$sos_iav_slp
)

