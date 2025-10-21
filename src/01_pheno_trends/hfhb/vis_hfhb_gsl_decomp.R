# ******************************************************************************
# Visualize the decomposition of GSL into spring and fall trends and IAV
# 
# Author: Xiaojie Gao
# Date: 2025-04-18
# ******************************************************************************
rm(list = ls())


SourcePartial <- function(filename, startLine = 1, endLine = -1) {
    lines <- scan(
        filename, 
        what = character(), sep = "\n", blank.lines.skip = FALSE,
        quiet = TRUE
    )
    st <- startLine
    en <- ifelse(endLine == -1, nrow(lines), endLine)
    tc <- textConnection(lines[(st):(en)])
    source(tc)
    close(tc)
}


# Need to manually run `vis_hf_pheno_trends.R` to get `gsl_trend_ols` and assign
# it to be `hf_gsl_trend_ols`
SourcePartial(
    filename = "src/01_pheno_trends/hfhb/vis_hf_pheno_trends.R",
    startLine = 8,
    endLine = 236
)
SourcePartial(
    filename = "src/01_pheno_trends/hfhb/vis_hf_gsl_decomp.R",
    startLine = 11,
    endLine = 132
)
hf_var_dt <- var_dt
hf_var_dt$site <- "HF"


# Need to manually run `vis_hb_pheno_trends.R` to get `gsl_trend_ols` and assign
# it to be `hb_gsl_trend_ols`
SourcePartial(
    filename = "src/01_pheno_trends/hfhb/vis_hb_pheno_trends.R",
    startLine = 8, 
    endLine = 311
)
SourcePartial(
    filename = "src/01_pheno_trends/hfhb/vis_hb_gsl_decomp.R",
    startLine = 11,
    endLine = 145
)
hb_var_dt <- var_dt
hb_var_dt$site <- "HB"



# ~ Calculate proportions ####
# ~ ----------------------------------------------------------------------------
hf_sos_iav_pct <- hf_var_dt[iav_diff > 0, .N / nrow(hf_var_dt)] * 100
hf_eos_iav_pct <- 100 - hf_sos_iav_pct

hf_eos_iav_sig_pct <- hf_var_dt[
    gsl_pval < 0.05 & iav_diff < 0, .N / nrow(hf_var_dt)
] * 100
hf_sos_iav_sig_pct <- hf_var_dt[
    gsl_pval < 0.05 & iav_diff > 0, .N / nrow(hf_var_dt)
] * 100

hf_eos_trend_pct <- hf_var_dt[trend_diff < 0, .N / nrow(hf_var_dt)] * 100
hf_sos_trend_pct <- 100 - hf_eos_trend_pct
hf_eos_trend_sig_pct <- hf_var_dt[
    gsl_pval < 0.05 & trend_diff < 0, .N / nrow(hf_var_dt)
] * 100
hf_sos_trend_sig_pct <- hf_var_dt[
    gsl_pval < 0.05 & trend_diff > 0, .N / nrow(hf_var_dt)
] * 100


hb_sos_iav_pct <- hb_var_dt[iav_diff > 0, .N / nrow(hb_var_dt)] * 100
hb_eos_iav_pct <- 100 - hb_sos_iav_pct

hb_eos_iav_sig_pct <- hb_var_dt[
    gsl_pval < 0.05 & iav_diff < 0, .N / nrow(hb_var_dt)
] * 100
hb_sos_iav_sig_pct <- hb_var_dt[
    gsl_pval < 0.05 & iav_diff > 0, .N / nrow(hb_var_dt)
] * 100

hb_eos_trend_pct <- hb_var_dt[trend_diff < 0, .N / nrow(hb_var_dt)] * 100
hb_sos_trend_pct <- 100 - hb_eos_trend_pct
hb_eos_trend_sig_pct <- hb_var_dt[
    gsl_pval < 0.05 & trend_diff < 0, .N / nrow(hb_var_dt)
] * 100
hb_sos_trend_sig_pct <- hb_var_dt[
    gsl_pval < 0.05 & trend_diff > 0, .N / nrow(hb_var_dt)
] * 100



{
    svglite::svglite(
        "out/hfhb_gsl_decomp.svg",
        width = 7, height = 7
    )
    par(mar = c(3, 3, 7, 1), mgp = c(1.5, 0.5, 0))

    xcoords <- barplot(
        matrix(
            c(
                hf_sos_iav_pct, hf_eos_iav_pct,
                hf_sos_trend_pct, hf_eos_trend_pct,
                hb_sos_iav_pct, hb_eos_iav_pct,
                hb_sos_trend_pct, hb_eos_trend_pct
            ),
            nrow = 2, byrow = FALSE
        ),
        col = c("seagreen", "orange"),
        ylab = "Percentage (%)"
    )
    abline(h = 50, lty = 2)
    axis(
        side = 1,
        at = xcoords,
        labels = c("IAV", "Trend", "IAV", "Trend"),
        lwd = 0, line = -0.5
    )
    axis(
        side = 1,
        at = c(mean(xcoords[1:2]), mean(xcoords[3:4])),
        labels = c("HF", "HB"),
        lwd = 0, line = 1,
        cex.axis = 1.5
    )

    # Add proportions with sig. GSL trends
    # HF - IAV
    rect(
        xcoords[1] - 0.5, hf_sos_iav_pct - hf_sos_iav_sig_pct,
        xcoords[1] + 0.5, hf_sos_iav_pct,
        density = 5, angle = 45
    )
    text(
        xcoords[1], (hf_sos_iav_pct - hf_sos_iav_sig_pct) + hf_sos_iav_sig_pct / 2,
        labels = paste0(round(hf_sos_iav_sig_pct, 2), "%"),
        col = "white"
    )
    rect(
        xcoords[1] - 0.5, hf_sos_iav_pct,
        xcoords[1] + 0.5, hf_sos_iav_pct + hf_eos_iav_sig_pct,
        density = 5, angle = 45
    )
    text(
        xcoords[1], hf_sos_iav_pct + hf_eos_iav_sig_pct / 2,
        labels = paste0(round(hf_eos_iav_sig_pct, 2), "%"),
        col = "white"
    )

    # HF - Trend
    rect(
        xcoords[2] - 0.5, hf_sos_trend_pct - hf_sos_trend_sig_pct,
        xcoords[2] + 0.5, hf_sos_trend_pct,
        density = 5, angle = 45
    )
    text(
        xcoords[2], (hf_sos_trend_pct - hf_sos_trend_sig_pct) + hf_sos_trend_sig_pct / 2,
        labels = paste0(round(hf_sos_trend_sig_pct, 2), "%"),
        col = "white"
    )
    rect(
        xcoords[2] - 0.5, hf_sos_trend_pct,
        xcoords[2] + 0.5, hf_sos_trend_pct + hf_eos_trend_sig_pct,
        density = 5, angle = 45
    )
    text(
        xcoords[2], hf_sos_trend_pct + hf_eos_trend_sig_pct / 2,
        labels = paste0(round(hf_eos_trend_sig_pct, 2), "%"),
        col = "white"
    )


    # HB - IAV
    rect(
        xcoords[3] - 0.5, hb_sos_iav_pct - hb_sos_iav_sig_pct,
        xcoords[3] + 0.5, hb_sos_iav_pct,
        density = 5, angle = 45
    )
    text(
        xcoords[3], (hb_sos_iav_pct - hb_sos_iav_sig_pct) + hb_sos_iav_sig_pct / 2,
        labels = paste0(round(hb_sos_iav_sig_pct, 2), "%"),
        col = "white"
    )
    rect(
        xcoords[3] - 0.5, hb_sos_iav_pct,
        xcoords[3] + 0.5, hb_sos_iav_pct + hb_eos_iav_sig_pct,
        density = 5, angle = 45
    )
    text(
        xcoords[3], hb_sos_iav_pct + hb_eos_iav_sig_pct / 2,
        labels = paste0(round(hb_eos_iav_sig_pct, 2), "%"),
        col = "white"
    )

    # HB - Trend
    rect(
        xcoords[4] - 0.5, hb_sos_trend_pct - hb_sos_trend_sig_pct,
        xcoords[4] + 0.5, hb_sos_trend_pct,
        density = 5, angle = 45
    )
    text(
        xcoords[4], (hb_sos_trend_pct - hb_sos_trend_sig_pct) + hb_sos_trend_sig_pct / 2,
        labels = paste0(round(hb_sos_trend_sig_pct, 2), "%"),
        col = "white"
    )
    rect(
        xcoords[4] - 0.5, hb_sos_trend_pct,
        xcoords[4] + 0.5, hb_sos_trend_pct + hb_eos_trend_sig_pct,
        density = 5, angle = 45
    )
    text(
        xcoords[4], hb_sos_trend_pct + hb_eos_trend_sig_pct / 2,
        labels = paste0(round(hb_eos_trend_sig_pct, 2), "%"),
        col = "white"
    )


    legend(
        grconvertX(0.5, "ndc"), grconvertY(1, "ndc"), 
        bty = "n",
        fill = c("seagreen", "orange"),
        legend = c(
            expression(
                C[leaf-out]^IAV ~ ">" ~ C[leaf-fall]^IAV ~ "or" ~ 
                C[leaf-out]^Trend ~ ">" ~ C[leaf-fall]^Trend
            ),
            expression(
                C[leaf-out]^IAV ~ "<" ~ C[leaf-fall]^IAV ~ "or" ~ 
                C[leaf-out]^Trend ~ "<" ~ C[leaf-fall]^Trend
            )
        ),
        xjust = 0.5,
        y.intersp = 1.5,
        cex = 1.5,
        xpd = NA
    )

    dev.off()

}



