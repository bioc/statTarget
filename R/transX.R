#' @name transX
#' @title transX for statTarget inputs
#' @description transX is to generate statTarget inputs from Mass Spectrometry 
#' Data softwares, like XCMS.
#' @param data A transX objects. The output file from Mass Spectrometry Data 
#' softwares. 
#' @param type The output file formats from Mass Spectrometry Data software, 
#' including "XCMS" or "xcms",...; The softwares include XCMS (.tsv file 
#' generated from diffreport function), ...
#' @return A objects of transX
#' @examples
#' datpath <- system.file("extdata",package = "statTarget")
#' data <- paste(datpath,"xcmsOutput.tsv", sep="/")
#' transX(data,"xcms")
#' transCode(data,"xcms")
#' @keywords XCMS 
#' @keywords inputs 
#' @export 
transX <- function(data, type){

         dirout.uni = paste(getwd(), "/statTargetDirectory/", sep = "")
         if (!file.exists("statTargetDirectory")){
           dir.create(dirout.uni)
           #setwd(dirout.uni)
         }
         #dirout.w = paste(getwd(), "/statTargetDirectory/",type, sep="")
         #dir.create(dirout.w)
         transX <- transCode(data=data,type=type)
         write.csv(transX$PhenoFile,paste(dirout.uni,"/PhenoFile_",type,".csv",
                                          sep = ""), row.names = FALSE)
         write.csv(transX$ProfileFile,
                   paste(dirout.uni,"/ProfileFile_", type,".csv",
                         sep = ""),row.names = FALSE)
         write.csv(transX$StatFile,
                   paste(dirout.uni,"/StatFile_", type,".csv",
                         sep = ""),row.names = FALSE)
         #setwd(dirout.uni)
         cat("Note: The inputs have been generated!", 
             "Filling the missing info. please!\n")
}

#' @name transCode
#' @title transCode for statTarget inputs
#' @description transCode is to generate statTarget inputs from Mass Spectrometry 
#' Data softwares, like XCMS.
#' @param data A transCode objects. The output file from Mass Spectrometry Data 
#' softwares. 
#' @param type The output file formats from Mass Spectrometry Data software, 
#' including "XCMS" or "xcms",...; The softwares include XCMS (.tsv file 
#' generated from diffreport function), ...
#' @return A list of inputs 
#' @examples
#' datpath <- system.file("extdata",package = "statTarget")
#' data <- paste(datpath,"xcmsOutput.tsv", sep="/")
#' transCode(data,"xcms")
#' @keywords XCMS 
#' @keywords inputs 
#' @export 
transCode <- function(data, type) {
  #write.table(datR,"xcmsOutput_true.tsv",sep = "\t")
  if(type == "XCMS" | type =="xcms"){
    if(!grepl(".tsv",data) == TRUE){
      stop("Read-only .tsv file from diffreport in XCMS software")
    } 
    datR <- utils::read.delim(data)
    #ProfileFile
    ProfileFile <- cbind(datR$name,datR[,15:ncol(datR)])
    colnames(ProfileFile)[1] <- c("name")
    sampleT <- c(colnames(ProfileFile[,2:ncol(ProfileFile)]))
    #PhenoFile
    PhenoFile <- as.data.frame(matrix(data=NA, nrow=length(sampleT),ncol = 4))
    colnames(PhenoFile) <- c("sample","batch", "class","order")
    PhenoFile$sample <- sampleT
    #statFile
    statF <- t(ProfileFile)
    colnames(statF) <- statF[1,]
    statF <- data.frame(statF)
    statF <- cbind(PhenoFile$sample,PhenoFile$class,statF[-1,])
    colnames(statF)[1:2] <- c("name","group")
    }
  dataOutput <- list(PhenoFile = PhenoFile, ProfileFile = ProfileFile,
                   StatFile = statF)
  return(dataOutput)
}






