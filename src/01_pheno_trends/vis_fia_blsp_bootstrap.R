# ******************************************************************************
# Visualize the bootstrap results
# 
# Author: Xiaojie Gao
# Date: 2026-06-24
# ******************************************************************************
library(data.table)



cachedir <- "pipe/blsp_trend_bootstrap"


# ~ Trends ####
# ~ ----------------------------------------------------------------------------
# For each iteration, summarize the proportions sig. trends
bootstrap_files <- list.files(cachedir, ".csv", full.names = TRUE)
bootstrap_sig_dt <- NULL
for (i in seq_along(bootstrap_files)) {
    site_trend_dt <- fread(bootstrap_files[i])

    total <- site_trend_dt[, .(TotN = .N),
        keyby = c("spec", "CommonName")
    ]

    sig_dt <- cbind(
        site_trend_dt[spr_pval < 0.05, .(Nspr = .N), keyby = "spec"],
        site_trend_dt[atm_pval < 0.05, .N, keyby = "spec"][, .(Natm = N)],
        site_trend_dt[gsl_pval < 0.05, .N, keyby = "spec"][, .(Ngsl = N)]
    )
    sig_dt <- merge(sig_dt, total, by = "spec")
    sig_dt[, ":="(
        spr_pct = Nspr / TotN * 100,
        atm_pct = Natm / TotN * 100,
        gsl_pct = Ngsl / TotN * 100)]

    sig_dt$rep_id <- i

    bootstrap_sig_dt <- rbind(bootstrap_sig_dt, sig_dt)
}


specnames <- unique(bootstrap_sig_dt$spec)

{ # fig: bootstrap trends
    png(
        "out/blsp_bootstrap_trends.png",
        width = 1300, height = 1800, res = 300
    )
    par(bty = "L", mar = c(4, 5, 1, 1), cex.lab = 1.2, cex.axis = 1.2, mgp = c(2, 0.7, 0))
    plot(NA,
        xlim = c(0, 35), ylim = c(0.5, 10), 
        yaxt = "n", ylab = "", xlab = "Sig. proportion (%)"
    )
    axis(side = 2, at = 1:10, labels = specnames, las = 2)

    for (i in seq_along(specnames)) {
        sp <- specnames[i]
        
        # Spring
        boxplot(
            bootstrap_sig_dt[spec == sp, spr_pct],
            at = i + 0.1, 
            boxwex = 0.2,
            horizontal = TRUE,
            col = "seagreen", 
            outline = FALSE,
            add = TRUE
        )
        # Autumn
        boxplot(
            bootstrap_sig_dt[spec == sp, atm_pct],
            at = i, 
            boxwex = 0.2,
            col = "orange",
            horizontal = TRUE,
            outline = FALSE,
            add = TRUE
        )
        # GSL
        boxplot(
            bootstrap_sig_dt[spec == sp, gsl_pct],
            at = i - 0.1, 
            boxwex = 0.2,
            col = "blue",
            horizontal = TRUE,
            outline = FALSE,
            add = TRUE
        )
    }
    dev.off()
}



# ~ IAV and Trend Contribution to GSL ####
# ~ ----------------------------------------------------------------------------
# For each iteration, summarize the var
bootstrap_var_dt <- NULL
for (i in seq_along(bootstrap_files)) {
    rep_dt <- fread(bootstrap_files[i])
    
    # Remove non-sense values
    rep_dt <- rep_dt[
        sos_iav_pct > 0 & sos_iav_pct < 100 &
        eos_iav_pct > 0 & eos_iav_pct < 100 &
        (sos_iav_pct + eos_iav_pct) > 100 &
        (sos_trend_pct + eos_trend_pct > 100)
    ]

    rep_dt[, ":="(
        iav_diff = sos_iav_pct - eos_iav_pct,
        trend_diff = sos_trend_pct - eos_trend_pct
    )]

    # Calculate some proportions
    sos_less_eos_iav_pct <- rep_dt[iav_diff < 0, .N / nrow(rep_dt)] * 100
    sos_greater_eos_iav_pct <- rep_dt[iav_diff > 0, .N / nrow(rep_dt)] * 100

    sos_less_eos_trend_pct <- rep_dt[trend_diff < 0, .N / nrow(rep_dt)] * 100
    sos_greater_eos_trend_pct <- rep_dt[trend_diff > 0, .N / nrow(rep_dt)] * 100

    bootstrap_var_dt <- rbind(bootstrap_var_dt, data.table(
        rep_id = i,
        sos_less_eos_iav_pct, sos_greater_eos_iav_pct,
        sos_less_eos_trend_pct, sos_greater_eos_trend_pct
    ))
}


mn_sos_greater_eos_iav <- round(mean(bootstrap_var_dt$sos_greater_eos_iav_pct))
sd_sos_greater_eos_iav <- round(sd(bootstrap_var_dt$sos_greater_eos_iav_pct))
mn_sos_less_eos_trend <- round(mean(bootstrap_var_dt$sos_less_eos_trend_pct))
sd_sos_less_eos_trend <- round(sd(bootstrap_var_dt$sos_less_eos_trend_pct))

message(
    "Prop. of SOS > EOS in IAV: ", 
    bootstrap_var_dt[sos_greater_eos_iav_pct > 50, .N / nrow(bootstrap_var_dt) * 100],
    "%"
)
message(
    "Prop. of SOS < EOS in Trend: ",
    bootstrap_var_dt[sos_less_eos_trend_pct < 50, .N / nrow(bootstrap_var_dt) * 100],
    "%"
)

{ # fig: Bootstrap contributions
    svglite::svglite(
        file.path("out", "bootstrap_gsl_decompose.svg"),
        width = 4, height = 6.5
    )
    divcols <- c("seagreen", "orange")
    par(
        mfrow = c(2, 1), mgp = c(1.7, 0.1, 0), mar = c(4, 3, 1, 1),
        tck = -0.01, cex.lab = 1.2
    )
    hist(
        bootstrap_var_dt$sos_greater_eos_iav_pct, 
        breaks = 100, xlim = c(0, 100), border = NA,
        xlab = expression(C[SOS]^IAV ~ ">" ~ C[EOS]^IAV), 
        main = ""
    )
    abline(v = 50, lwd = 2, lty = 2, col = "black")
    abline(v = mn_sos_greater_eos_iav, lwd = 2, lty = 1, col = "red")
    text(
        mn_sos_greater_eos_iav, grconvertY(0.8, "npc"), 
        paste0( 
            mn_sos_greater_eos_iav, "%", 
            "(\u00B1", sd_sos_greater_eos_iav, "%)"
        ),
        pos = 4, xpd = NA, 
        col = "red"
    )

    
    hist(
        bootstrap_var_dt$sos_less_eos_trend_pct, 
        breaks = 100, xlim = c(0, 100), ylim = c(0, 30), border = NA,
        xlab = expression(C[SOS]^Trend ~ "<" ~ C[EOS]^Trend), 
        main = ""
    )
    abline(v = 50, lwd = 2, lty = 2, col = "black")
    abline(v = mn_sos_less_eos_trend, lwd = 2, lty = 1, col = "red")
    text(
        mn_sos_less_eos_trend, grconvertY(0.8, "npc"),
        paste0(
            mn_sos_less_eos_trend, "%", 
            "(\u00B1", sd_sos_less_eos_trend, "%)"
        ),
        pos = 4, xpd = NA,
        col = "red"
    )
    dev.off()
}

