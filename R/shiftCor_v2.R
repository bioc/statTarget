#' @name shiftCor
#' @title shiftCor
#' @description  shiftCor provides the QC based signal correction for 
#' large scale metabolomics and targeted proteomics.
#' @param samPeno File path. The file with the meta information including the sample name,
#'  batches, class and order. 
#' @param samFile File path. The file with the expression information. 
#' @param Frule Modified n precent rule function. A variable will be kept if it has a non-zero value
#' for at least n precent of samples in any one group. 
#' (Default: 0.8)  
#' @param MLmethod The machine learning method for QC based signal correction, such as QC based random forest signal correction (QC-RFSC). QC-RLSC was deprecated .
#' @param ntree Number of trees to grow in random forest model.
#' @param imputeM The parameter for imputation method i.e., nearest neighbor 
#' averaging, 'KNN'; minimum values, 'min'; Half of minimum values, 'minHalf'; 
#' median values, 'median'.
#' @param coCV Define the cutoff value (0-100) of CV for controlling the number of features.
#' @param plot Defines if images of feature quality should be generated (TRUE) or not (FALSE). 
#' Defaults to FALSE.
#' @return the shiftCor files. See the details at https://stattarget.github.io
#' @examples 
#' datpath <- system.file('extdata',package = 'statTarget')
#' samPeno <- paste(datpath,'MTBLS79_sampleList.csv', sep='/')
#' samFile <- paste(datpath,'MTBLS79.csv', sep='/')
#' samPeno
#' samFile
#' shiftCor(samPeno,samFile, MLmethod = 'QCRFSC', imputeM = 'KNN',coCV = 30)
#' @keywords Quality Controls,Correction
#' @export 
shiftCor <- function(samPeno, samFile, Frule = 0.8, MLmethod = "QCRFSC", ntree = 500,imputeM = "KNN", coCV = 30, plot = FALSE) {
    cat("\n")
    
    
    
    cat("statTarget: Signal Correction Start... Time:", date(), "\n\n")
    cat("* Step 1: Data File Checking Start..., Time: ", date(), "\n")
    
    cat("\n", "Data Link", "\n")
    cat(" metaFile:", samPeno, "\n")
    cat(" profileFile:", samFile, "\n")
    
    # read the data
    samPeno <- read.csv(samPeno, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
    samPeno <- as.data.frame(samPeno)
    # check the title
    checkPeno <- match(c("sample", "batch",  "class",  "order" ),colnames(samPeno))
    if(sum(is.na(checkPeno)) > 0) stop("The names in column of metaFile should be `sample`, `batch`,  `class`,  `order`")
    samPeno <- plyr::arrange(samPeno, order)
    samFile <- read.csv(samFile, header = FALSE, check.names = FALSE, stringsAsFactors = FALSE)
    samFile <- t(samFile)
    colnames(samFile) <- samFile[1, ]
    samFile <- as.data.frame(samFile[-1, ])
    rownames(samFile) <- samFile[,1]

    
    ############## Checking the input file#############
    
    
    ## .............................................
    
    ## Data Matching
    
    ## .............................................
    
    cat("\n", " ", dim(samPeno)[1], " Meta Samples vs ", dim(samFile)[1], " Profile samples", sep = "")
    cat("\n", "The Meta samples list (*NA, missing data from the Profile File)")
    cat("\n\n")
    mcdat <- samFile[, 1][match(samPeno[, 1], samFile[, 1])]
    print(as.vector(mcdat))
    
    ## .............................................
    if (any(is.na(mcdat))) {
        stop("\n", "Missing data from the Profile File! Check your data please!!")
    }
    
    if (dim(samFile)[1] - dim(samPeno)[1] > 0) {
        cat("\n", "Warning: The sample size in Profile File is larger than Meta File! ")
        cat("\n")
    } else if (dim(samFile)[1] - dim(samPeno)[1] < 0) {
        stop("\n", "The sample size in Profile File should be no less than Meta File!", " Check your data please!!")
    }
    
    samFP <- samFile[samPeno$sample, ]
    
    if (sum(is.na(samPeno$class)) <= 0) {
        stop("\n", "There were no QC samples in your data! NA should be defined for QC samples in `class` ")
    }
    
    samPeno_stat <- samPeno
    rownames(samPeno_stat) <- samPeno[, 1]
    qc_seq <- rownames(samPeno_stat)
    qc_seq_tmp <- grep("QC|qc|Qc", qc_seq)
    
    
    
    sam_seq_tmp <- samPeno_stat[-c(qc_seq_tmp), ]
    
    if (sum(is.na(samPeno_stat$batch)) > 0) {
        stop("\n", "There were missing values in batch! Check your data please!\n")
    }
    if (sum(is.na(sam_seq_tmp$class)) > 0) {
        stop("\n", "There were missing values (Unclassified data) in sample class! 
         Check your data please!\n")
    }
    
    
    samPeno_stat$class[is.na(samPeno_stat$class)] <- "QC"
    data_stat = aggregate(samPeno_stat$class, by = list(Category = samPeno_stat$class), FUN = length)
    batch_stat = aggregate(samPeno_stat$class, by = list(Category = samPeno_stat$batch), FUN = length)
    colnames(data_stat) <- c("Class", "No.")
    colnames(batch_stat) <- c("Batch", "No.")
    cat("\n", "Meta-information:", "\n\n")
    print(data_stat)
    print(batch_stat)
    num_sam <- dim(samFile[, 2:ncol(samFile)])
    num_sam <- data.frame(num_sam)
    colnames(num_sam) <- c("no.")
    rownames(num_sam) <- c("QC and samples", "Metabolites")
    cat("\n", "Metabolic profile information:", "\n\n")
    print(num_sam)
    
    
    ################# data NA########
    samFP_temp <- as.matrix(samFP[, 2:ncol(samFP)])
    
    if (length(samFP_temp[samFP_temp < 0L]) > 0) {
        samFP_temp[samFP_temp < 0L] <- 0L
    }
    samFP_temp[samFP_temp == 0L] <- NA
    
    cat("\n")
    cat("* Step 2: Evaluation of Missing Value...", "\n")
    cat("\n", "The number of missing value before QC based signal correction: ", sum(is.na(samFP_temp)))
    imsamFP <- samFP_temp
    ############# Filter miss value###################
    
    classF <- as.factor(samPeno$class)
    classF = addNA(classF)
    imsamFPF = FilterMV(imsamFP, classF, rule = Frule)
    Frule_warning = paste("The number of filtered variables using the modified ", Frule * 100, "% rule :", 
        sep = " ")
    cat("\n", Frule_warning, " ", dim(imsamFP)[2] - dim(imsamFPF)[2], "\n")
    imsamFP = as.matrix(imsamFPF)
    ############## impute missing value#################
    cat("\n")
    cat("* Step 3: Imputation start...", "\n")
    
    
    if (imputeM == "KNN") {
        # require(impute)
        cat("\n", "The imputation method was set at 'KNN'")
        mvd <- impute::impute.knn(imsamFP, rowmax = 0.99, colmax = 0.99, maxp = 15000)
        inputedData <- mvd$data
    } else if (imputeM == "min") {
        cat("\n", "The imputation method was set at 'min'")
        

        inputedData = minValue(imsamFP, classF)
        # inputedData = inputedData[,-1]
        
    } else if (imputeM == "minHalf") {
        cat("\n", "The imputation method was set at 'minHalf'")
        
        inputedData = minHalfValue(imsamFP, classF)
        # inputedData = inputedData[,-1]
        
    } else if (imputeM == "median") {
        

        
        cat("\n", "The imputation method was set at 'median'")
        
        inputedData = medianvalue(imsamFP, classF)
        # inputedData = inputedData[,-1]
    }
    cat("\n", "The number of missing value after imputation: ", sum(is.na(inputedData)))
    
    
    if (sum(is.na(inputedData)) > 0) {
        inputedData = minHalfValue(inputedData, classF)
        cat("\n", "The number of missing value after the second imputation: ", sum(is.na(inputedData)))
    }
    
    
    cat("\n", "Imputation Finished!", "\n", "\n")
    
    cat("* Step 4: QC-based Signal Correction Start... Time: ", date(), "\n")
    
    dat <- as.matrix(t(inputedData))
    numX <- 1:dim(dat)[2]
    
    
    if (MLmethod == "QCRFSC") {
        # require(randomForest)
        cat("\n", "The Signal Correction method was set at QC-RFSC", "\n")
        

        loessDat <- REGfit(x = dat, y = numX, ntree = ntree)
    }
    
    #################################################### 

    ############### dataCheck
    loessDatmp <- apply(loessDat, 2, function(x) as.numeric(as.character(x)))
    loessDatT <- loessDatmp * 1000
    rownames(loessDatT) <- rownames(dat)
    
    datmp <- apply(dat, 2, function(x) as.numeric(as.character(x)))
    rownames(datmp) <- rownames(loessDatT)
    
    
    if (length(loessDatT[loessDatT < 0L]) > 0) {
        loessDatT[loessDatT < 0L] <- 0
    }
    loessDatT[loessDatT == 0L] <- NA
    if (sum(is.na(loessDatT) > 0)) {

        loessDatT <- minHalfValue(t(loessDatT), classF)
        loessDatT <- t(loessDatT)
        cat("\n", "The number of missing value after QC-based signal correction: ", sum(is.na(loessDatT)), 
            "\n")
    }
    dirout.uni = paste(getwd(), "/statTarget/", sep = "")
    dirsc.ID = getwd()
    dir.create(dirout.uni)
    dirout.w = paste(getwd(), "/statTarget/shiftCor", sep = "")
    dir.create(dirout.w)
    dirout.Bs = paste(getwd(), "/statTarget/shiftCor/Before_shiftCor", sep = "")
    dir.create(dirout.Bs)
    dirout.As = paste(getwd(), "/statTarget/shiftCor/After_shiftCor", sep = "")
    dir.create(dirout.As)
    ########### Output images of each peak quality ############
    if(plot){
      
      cat("\n", "High-resolution images output...", "\n")
      pbloplot <- txtProgressBar(min = 1, max = dim(datmp)[1], style = 3)
      for (i in 1:dim(datmp)[1]) {
        loplot(datmp, loessDatT, i)
        setTxtProgressBar(pbloplot, i)
      }
      close(pbloplot)
    }
    
    ############### Raw output###########
    
    raw_temp <- cbind(samPeno, inputedData)
    nam_qc <- rownames(raw_temp)
    
    ##### QC Cal################
    QC_temp_raw <- grep("QC|qc|Qc", nam_qc)
    QC_temp_raw <- raw_temp[c(QC_temp_raw), ]
    raw_temp_qc <- QC_temp_raw[, -c(3, 4)]
    rownames(raw_temp_qc) <- NULL
    RSD30_CV = paste("shift_QC_raw", ".csv", sep = "")
    write.csv(raw_temp_qc, paste(dirout.Bs, RSD30_CV, sep = "/"))
    cat("\n", "Calculation of CV distribution of raw peaks (QC)...\n\n")
    # raw_temp_qc = cbind(seq(1,dim(raw_temp_qc)[1],1),raw_temp_qc)
    Rsdist_QC_raw = RsdCal(raw_temp_qc, batch = TRUE, DistPattern = TRUE, output = FALSE)
    
    ##### Sample Cal################
    sam_temp_raw <- grep("QC|qc|Qc", nam_qc)
    sam_temp_raw <- raw_temp[-c(sam_temp_raw), ]
    raw_temp_sam <- sam_temp_raw[, -c(3, 4)]
    rownames(raw_temp_sam) <- NULL
    RSD30_CV = paste("shift_sam_raw", ".csv", sep = "")
    write.csv(raw_temp_sam, paste(dirout.Bs, RSD30_CV, sep = "/"))
    # Rsdist_sam_raw = RsdCal(raw_temp_sam,DistPattern = FALSE) raw_temp_sam =
    # cbind(seq(1,dim(raw_temp_sam)[1],1),raw_temp_sam)
    Rsdist_sam_raw = RsdCal(raw_temp_sam, batch = FALSE, DistPattern = FALSE, output = FALSE)
    
    
    ############### Loess output###########
    
    lo_temp <- cbind(samPeno, t(loessDatT))
    nam_qc <- rownames(lo_temp)
      
    ############ loess QC Cal#############
    QC_temp <- grep("QC|qc|Qc", nam_qc)
    QC_temp <- lo_temp[c(QC_temp), ]
    lo_temp_qc <- QC_temp[, -c(3, 4)]
    rownames(lo_temp_qc) <- NULL
    cat("\n", "Calculation of CV distribution of corrected peaks (QC)...\n\n")
    # Rsdist_QC_cor=RsdCal(lo_temp_qc,DistPattern = TRUE) lo_temp_qc =
    # cbind(seq(1,dim(lo_temp_qc)[1],1),lo_temp_qc)
    Rsdist_QC_cor = RsdCal(lo_temp_qc, batch = TRUE, DistPattern = TRUE, output = TRUE)
    
    # cutoff CV
    
    cfid <- which(Rsdist_QC_cor > coCV/100) + 2
    
    # filter lo_temp_qc
    if(length(cfid) == 0) {
      lo_temp_qc_filter <- lo_temp_qc
    } else {
    lo_temp_qc_filter <- lo_temp_qc[,-cfid]}
    
    #output
    RSD30_CV = paste("shift_QC_cor", ".csv", sep = "")
    
    lo_temp_qc_filter <- t(lo_temp_qc_filter)
    write.table(lo_temp_qc_filter, paste(dirout.As, RSD30_CV, sep = "/"), sep = ",",col.names = FALSE)
    
    
    ############ loess sam Cal#############
    QC_temp <- grep("QC|qc|Qc", nam_qc)
    QC_temp <- lo_temp[-c(QC_temp), ]
    lo_temp_sam <- QC_temp[, -c(2, 4)]
    rownames(lo_temp_sam) <- NULL
    
    # filter 
    if(length(cfid) == 0) {
      lo_temp_sam_filter <- lo_temp_sam
    } else {
      lo_temp_sam_filter <- lo_temp_sam[,-cfid]}
    #lo_temp_sam_filter <- lo_temp_sam[,-cfid]
    
    RSD30_CV = paste("shift_sample_cor", ".csv", sep = "")
    
    lo_temp_sam_filter <- t(lo_temp_sam_filter)
    write.table(lo_temp_sam_filter, paste(dirout.As, RSD30_CV, sep = "/"), sep = ",",col.names = FALSE)
    
    ############ loess sample Cal############# Rsdist_sam_cor = RsdCal(lo_temp_sam,DistPattern = FALSE)
    ############ lo_temp_sam = cbind(seq(1,dim(lo_temp_sam)[1],1),lo_temp_sam)
    Rsdist_sam_cor = RsdCal(lo_temp_sam, batch = FALSE, DistPattern = FALSE, output = FALSE)
    
    
    
    ################# 
    RSDdist(Rsdist_sam_raw, Rsdist_sam_cor, Rsdist_QC_raw, Rsdist_QC_cor)
    lo_temp$class[is.na(lo_temp$class)] <- "QC"
    lo_temp_all <- lo_temp[, -c(2, 4)]
    rownames(lo_temp_all) <- NULL
    
    # fitler
    if(length(cfid) == 0) {
      lo_temp_all_filter <- lo_temp_all
    } else {
      lo_temp_all_filter <- lo_temp_all[,-cfid]}
    #lo_temp_all_filter <- lo_temp_all[,-cfid]
    
    RSD30_CV = paste("shift_all_cor", ".csv", sep = "")
    
    lo_temp_all_filter <- t(lo_temp_all_filter)
    write.table(lo_temp_all_filter, paste(dirout.As, RSD30_CV, sep = "/"), sep = ",",col.names = FALSE)
    
    cat("\n* Step 5: ", paste("Removal of the features (CV% > ", coCV,"%)",sep=""),date(), "\n")
    
    dn <- length(colnames(lo_temp_qc)[cfid])
    cat("\n", "No. of removed features:", dn)
    cat("\n","Feature name:", colnames(lo_temp_qc)[cfid], "\n")
    
    cat("\n", "Output Link:", getwd(), "\n")
    cat("\n", "Correction Finished! Time: ", date(), "\n")
    cat("\n", "####################################", "\n")
    cat(" # Software Version: statTarget 2.0 + #", "\n")
    cat(" # Email: hemi.luan@gmail.com #", "\n")
    cat(" ####################################", "\n")
    
    ################## Loess Plot########################
    setwd(dirsc.ID)
    
    # parameter output
    scPam1 <- c("Frule", "MLmethod", "ntree", "imputeM","coCV")
    scPam2 <- c(Frule, MLmethod, ntree, imputeM, coCV)
    scpam <- data.frame(scPam1, scPam2)
    colnames(scpam) <- c("parameter", "value")
    par_sh = paste("statTarget/ParameterShiftCor", ".log", sep = "")
    write.table(scpam, paste(getwd(), par_sh, sep = "/"), row.names = FALSE)
    
    tmpfilesc = paste(getwd(), "/tmp", sep = "")
    unlink(tmpfilesc, recursive = TRUE)
}
