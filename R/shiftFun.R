FilterMV = function(m, degree,rule) {
  dx <- c()
  for (i in 1:ncol(m)) {
    freq <- as.vector(tapply(m[, i], degree, function(x) {
      sum(is.na(x))/length(x)
    }))
    if (sum(freq > 1 - rule) > 0) 
      dx <- c(dx, i)
  }
  if (length(dx) > 0) 
    m <- m[, -dx]
  return(m)
}
minValue <- function(x, group) {
  group = as.factor(as.numeric(group))
  for (i in 1:dim(x)[1]) {
    for (j in 1:dim(x)[2]) {
      if (is.na(x[i, j]) == TRUE) {
        x[i, j] <- tapply(as.numeric(x[, j]), group, min, na.rm = TRUE)[group[i]]
      }
    }
  }
  return(x)
}
minHalfValue <- function(x, group) {
  group = as.factor(as.numeric(group))
  for (i in 1:dim(x)[1]) {
    for (j in 1:dim(x)[2]) {
      if (is.na(x[i, j]) == TRUE) {
        x[i, j] <- tapply(as.numeric(x[, j]), group, min, na.rm = TRUE)[group[i]]/2
      }
    }
  }
  return(x)
}
medianvalue <- function(x, group) {
  group = as.factor(as.numeric(group))
  for (i in 1:dim(x)[1]) {
    for (j in 1:dim(x)[2]) {
      if (is.na(x[i, j]) == TRUE) {
        x[i, j] <- tapply(as.numeric(x[, j]), group, median, na.rm = TRUE)[group[i]]
      }
    }
  }
  return(x)
}


REGfit = function(x, y, ntree = ntree) {
  cn <- colnames(x)
  x <- as.matrix(x)
  #### Check########
  st_QC <- grep("QC|qc|Qc", cn[1])
  ed_QC <- grep("QC|qc|Qc", cn[length(cn)])
  if (length(st_QC) == 0) {
    stop("\nWrong: the first sample must be QC sample; please check ......")
  }
  if (length(ed_QC) == 0) {
    stop("\nWrong: the sample at the end of sequence must be QC sample; 
          please check ......")
  }
  qcid <- grep("QC|qc|Qc", cn)
  pb <- txtProgressBar(min = 1, max = dim(x)[1], style = 3)
  for (i in 1:dim(x)[1]) {
    temp <- randomForest(data.frame(qcid), as.numeric(x[i, qcid]), ntree = ntree)
    y <- data.frame(y)
    colnames(y) <- "qcid"
    rfP <- predict(temp, y)
    x[i, ] <- as.numeric(x[i, ])/rfP
    setTxtProgressBar(pb, i)
  }
  close(pb)
  loessDat = x
  return(loessDat)
}

loessFit = function(x, y, QCspan, degree) {
  cn <- colnames(x)
  #### Check########
  st_QC <- grep("QC|qc|Qc", cn[1])
  ed_QC <- grep("QC|qc|Qc", cn[length(cn)])
  if (length(st_QC) == 0) {
    stop("the first sample must be QC sample; please check ......")
  }
  if (length(ed_QC) == 0) {
    stop("the sample at the end of sequence must be QC sample; 
          please check ......")
  }
  qcid <- grep("QC|qc|Qc", cn)
  
  pb <- txtProgressBar(min = 1, max = dim(x)[1], style = 3)
  
  for (i in 1:dim(x)[1]) {
    loe <- stats::loess(x[i, qcid] ~ qcid, span = QCspan, degree = degree)
    yf <- stats::predict(loe, y)
    x[i, ] <- as.numeric(x[i, ])/yf
    setTxtProgressBar(pb, i)
  }
  close(pb)
  loessDat = x
  return(loessDat)
}

sploe <- function(sp) {
  loe2 <- get("loe1", envir = env)
  mod <- stats::update(loe2, span = sp)
  CVspan = loessGCV(mod)[["gcv"]]
}

minHalfValue <- function(x, group) {
  group = as.factor(as.numeric(group))
  for (i in 1:dim(x)[1]) {
    for (j in 1:dim(x)[2]) {
      if (is.na(x[i, j]) == TRUE) {
        x[i, j] <- tapply(as.numeric(x[, j]), group, min, na.rm = TRUE)[group[i]]/2
      }
    }
  }
  return(x)
}

autoFit <- function(xl, y) {
  cn <- colnames(xl)
  #### Check########
  st_QC <- grep("QC|qc|Qc", cn[1])
  ed_QC <- grep("QC|qc|Qc", cn[length(cn)])
  if (length(st_QC) == 0) {
    stop("the first sample must be QC sample; please check ......")
  }
  if (length(ed_QC) == 0) {
    stop("the sample at the end of sequence must be QC sample; 
          please check ......")
  }
  qcid <- grep("QC|qc|Qc", cn)
  pb <- txtProgressBar(min = 1, max = dim(xl)[1], style = 3)
  
  for (i in 1:dim(xl)[1]) {
    
    Sys.sleep(1e-06)
    
    loe1 <- loess(xl[i, qcid] ~ qcid)
    # loe2 <- loe1
    env <- environment()
    
    sp <- c(seq(0.5, 0.75, 0.01))
    CVspan = as.matrix(lapply(sp, sploe))
    CVspan[!is.finite(as.numeric(CVspan))] <- NA
    minG <- data.frame(sp, CVspan)
    minspan <- minG[which.min(minG[, 2]), 1]
    minspan
    # sp <- c(seq(0.05,0.75,0.01)) CVspan <-c() for(j in 1:length(sp)){ mod <- stats::update(loe1,
    # span = sp[j]) CVspan[j] = loessGCV(mod)[['gcv']] } minG <- as.matrix(data.frame(sp,CVspan))
    # minG[!is.finite(minG)] <- max(minG[,2],na.rm = TRUE) minspan <- minG[which.min(minG[,2]),1]
    # minspan
    loeN <- stats::update(loe1, span = minspan)
    yf <- predict(loeN, y)
    xl[i, ] <- as.numeric(xl[i, ])/yf
    setTxtProgressBar(pb, i)
  }
  close(pb)
  loessDat = xl
  return(loessDat)
}

