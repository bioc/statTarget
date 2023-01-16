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
