# ******************************************************************************
# Extract MSLSP30NA phenology for FIA points
# 
# Author: Xiaojie Gao
# Date: 2024-01-17
# ******************************************************************************
library(terra)
library(data.table)
library(parallel)


majorspec_dt <- fread("data/fia_domspec_truexy.csv")
majorspec_vect <- vect(majorspec_dt, geom = c("x", "y"), crs = "EPSG:4326")


hls_tiles <- vect("data/raw/HLS_tiles/sentinel2_tiles_north_america.shp")

msimgfiles <- extract(hls_tiles, majorspec_vect)
msimgfiles <- unique(msimgfiles$Name)
msimgfiles <- na.omit(msimgfiles)


ExtractImg <- function(imgname, majorspec_dt) {
    majorspec_vect <- vect(majorspec_dt, geom = c("x", "y"), crs = "EPSG:4326")

    msimg <- rast(imgname)
    yr <- substr(basename(imgname), 
        nchar(basename(imgname)) - 6, 
        nchar(basename(imgname)) - 3
    )

    sp_vec <- project(majorspec_vect, crs(msimg))
    cell_ids <- cells(msimg, sp_vec)

    res_dt <- apply(cell_ids, 1, function(cellid) {
        if (is.na(cellid[2])) {
            return(NULL)
        } 

        # Get 3x3 window
        rowcol <- rowColFromCell(msimg, cellid[2])
        upper_row <- rowcol[1] - 1
        botoom_row <- rowcol[1] + 1
        left_col <- rowcol[2] - 1
        right_col <- rowcol[2] + 1

        cells_to_use <- expand.grid(
            c(upper_row, rowcol[1], botoom_row),
            c(left_col, rowcol[2], right_col)
        )

        cells <- apply(cells_to_use, 1, function(x) {
            cellFromRowCol(msimg, x[1], x[2])
        })
        
        
        # Only use the first cycle
        dt <- as.data.table(msimg[cells])
        vals <- dt[gupQA %in% c(1, 2), .(
            id = cellid[1],
            yr = yr,

            avgOGI = mean(OGI, na.rm = TRUE),
            avg50PCGI = mean(`50PCGI`, na.rm = TRUE),
            avgOGMx = mean(OGMx, na.rm = TRUE),
            avgPeak = mean(Peak, na.rm = TRUE),
            avgOGD = mean(OGD, na.rm = TRUE),
            avg50PCGD = mean(`50PCGD`, na.rm = TRUE),
            avgOGMn = mean(OGMn, na.rm = TRUE),
            avgEVImax = mean(EVImax, na.rm = TRUE),
            avgEVIamp = mean(EVIamp, na.rm = TRUE),
            avgEVIarea = mean(EVIarea, na.rm = TRUE),
            
            medOGI = median(OGI, na.rm = TRUE),
            med50PCGI = median(`50PCGI`, na.rm = TRUE),
            medOGMx = median(OGMx, na.rm = TRUE),
            medPeak = median(Peak, na.rm = TRUE),
            medOGD = median(OGD, na.rm = TRUE),
            med50PCGD = median(`50PCGD`, na.rm = TRUE),
            medOGMn = median(OGMn, na.rm = TRUE),
            medEVImax = median(EVImax, na.rm = TRUE),
            medEVIamp = median(EVIamp, na.rm = TRUE),
            medEVIarea = median(EVIarea, na.rm = TRUE)
        )]
        
        return(vals)
    })
    res_dt <- do.call(rbind, res_dt)
    com_dt <- merge(cell_ids, res_dt, by.x = "ID", by.y = "id", all.x = TRUE)
    setDT(com_dt)

    final_dt <- na.omit(cbind(as.data.table(sp_vec), com_dt))
        

    return(final_dt)
}


# Extract each FIA site using a 3x3 window
cl <- makeCluster(30)
calls <- clusterCall(cl, function() {
    suppressWarnings({
        library(terra)
        library(data.table)
        library(parallel)
    })
})
clusterExport(cl, c("msimgfiles", "ExtractImg", "majorspec_dt"))

sp_pheno_dt <- clusterApplyLB(cl, 1:length(msimgfiles), function(i) {
    img_dt <- NULL
    for (yr in 2016:2019) {
        imgname <- paste0("MSLSP_", msimgfiles[i], "_", yr, ".nc")
        imgname <- file.path("data/raw/fia_mslsp30na", imgname)
        thedt <- ExtractImg(imgname, majorspec_dt)
    
        img_dt <- rbind(img_dt, thedt)
    }

    return(img_dt)
})
sp_pheno_dt <- do.call(rbind, sp_pheno_dt)
setDT(sp_pheno_dt)

stopCluster(cl)

# out: "data/sp_60_mslsp_dt.csv"
fwrite(sp_pheno_dt, "data/sp_60_mslsp_dt.csv")


