# ******************************************************************************
# Base variables and functions
# 
# Author: Xiaojie Gao
# Date: 2024-05-21
# ******************************************************************************


# Github Gist for my common functions
source("https://gist.githubusercontent.com/MrJGao/7bd6f978771746fd27a999c9a5c808c6/raw/58d6f9b2fd052bcd5c824aa36450c907b4e912ba/gao_fun.R")

base <- list()

# The original big root project directory
base$imax_fia_proj_dir <- "D:/Gao/projects/imax_fia"

# The Region 9 shapefile
base$region9_shpfile <- "Y:/FICE/Northeast20_Salata/SpatialData/Northern20.shp"


base$FitLm <- function(x, y) {
    if (all(is.na(x)) == TRUE || all(is.na(y)) == TRUE) {
        warning("Data is all NA!")
        return(NULL)
    }
    fit <- lm(y ~ x)
    if (any(is.na(coef(fit))) == TRUE) {
        warning("Coef NA!")
        return(NULL)
    }
    cf <- coef(fit)
    eq <- paste("y=", ifelse(sign(cf[2]) == 1, "", "-"), abs(cf[2]), "x",
        ifelse(sign(cf[1]) == 1, "+", "-"), abs(cf[1]),
        sep = ""
    )
    r_sqr <- round(summary(fit)$r.squared, 2)

    f <- summary(fit)$fstatistic
    p_val <- pf(f[1], f[2], f[3], lower.tail = FALSE)

    RSS <- c(crossprod(fit$residuals))
    MSE <- RSS / length(fit$residuals)
    rmse <- sqrt(MSE)

    return(list(
        fit = fit,
        cf = cf, eq = eq,
        r2 = r_sqr, pval = p_val,
        rmse = rmse
    ))
}

# Decompose time series into IAV and trend
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


# Calculate trend of IAV using a moving window
CalIAVTrend <- function(iav_signal, yrs) {
    # Remove outliers
    mn_iav <- mean(iav_signal)
    sd_iav <- sd(iav_signal)
    iav_signal[iav_signal > (mn_iav + 3 * sd_iav)] <- NA
    iav_signal[iav_signal < (mn_iav - 3 * sd_iav)] <- NA

    window_size <- 9
    side_size <- window_size / 2 - 0.5
    start_idx <- side_size + 1
    end_idx <- length(iav_signal) - side_size
    iav_var <- sapply(start_idx:end_idx, function(x) {
        var(iav_signal[(x - side_size):(x + side_size)])
    })
    iav_var <- c(rep(NA, side_size), iav_var, rep(NA, side_size))

    mod <- base$FitLm(yrs, iav_var)

    if (is.null(mod)) {
        return(data.table(
            slp = NA,
            pval = NA
        ))
    }

    return(data.table(
        slp = mod$cf[2],
        pval = mod$pval
    ))
}


smoothScatter_J2 <- function(
    x, y, colorful = FALSE,
    nbin = 300, nrpoints = 0, bgCol = rgb(1, 1, 1, 0),
    plotPoints = FALSE, pal = NULL,
    ...
) {
    # Use Type-II regression
    require(lmodel2)

    range <- range(x, y, na.rm = TRUE)

    if (plotPoints == TRUE) {
        plot(x, y,
            xlim = range, ylim = range,
            pch = 21, bg = "grey", mgp = c(1.5, 0.5, 0),
            ...
        )
    } else {
        if (isTRUE(colorful)) {
            require(RColorBrewer)
            if (is.null(pal)) {
                pal <- colorRampPalette(
                    c(
                        bgCol, "#f9f9fd", "#a3a6e1",
                        rev(brewer.pal(11, "Spectral"))
                    ),
                    alpha = TRUE
                )
            }
            smoothScatter(x, y,
                colramp = pal,
                xlim = range, ylim = range,
                nbin = nbin, nrpoints = nrpoints,
                ...
            )
        } else {
            smoothScatter(x, y,
                xlim = range, ylim = range,
                nbin = nbin, nrpoints = nrpoints,
                ...
            )
        }
    }
    abline(a = 0, b = 1, lty = 2)

    fit <- tryCatch(
        {
            suppressMessages(
                lmodel2(y ~ x, data.frame(x = x, y = y), "interval", "interval")
            )
        },
        error = function(x) {
            return(NULL)
        }
    )

    if (!is.null(fit) & any(is.na(coef(fit))) == FALSE) {
        fit_res <- data.table(fit$regression.results)
        inter_slp <- fit_res[Method == "RMA", .(Intercept, Slope)]
        abline(
            a = inter_slp$Intercept, b = inter_slp$Slope, lwd = 2,
            col = "grey"
        )

        # get stats
        r_sqr <- round(fit$rsquare, 2)
        r_sqr_str <- bquote(R^2 == .(r_sqr))

        rmse <- round(sqrt(mean((y - x)^2, na.rm = TRUE)), 2)
        rmse_str <- bquote(RMSE == .(rmse))

        if (is.na(r_sqr)) {
            rmse_str <- "RMSE=NaN"
        }

        # legend("bottomright",
        #     legend = bquote(n == .(length(x))),
        #     bty = "n", cex = 1.2, y.intersp = 5
        # )
        legend("bottomright",
            legend = r_sqr_str, bty = "n",
            cex = 1.2, y.intersp = 3
        )
        legend("bottomright",
            legend = rmse_str, bty = "n",
            cex = 1.2, y.intersp = 1
        )
    }
}