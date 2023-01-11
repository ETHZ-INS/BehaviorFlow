#' Reads labels from a folder containing one or multiple CSVs with label data. rows correspond to frames of the recording and column(s) to label(s). for each sample a separate file is needed. headers of the columns correspond to label names used downstream. 
#' 
#' @param path a character with the path to the folder containing the CSVs
#' @param sep character. the separator used in CSV files (defaults to german system ";"). for english system use "," instead.
#' @return An object of type USData
#' @examples
#' US <- LoadFromCSVs("C:/CSVfolder")
LoadFromCSVs <- function(path, sep = ";"){
  out <- list()
  out$files <- list()
  fnames <- list.files(path)
  
  #add data
  for(i in fnames){
    dat <- read.table(paste(path,i,sep = "/"), sep = sep, header = T)
    out$files[[i]] <- list()
    out$files[[i]]$data <- list()
    out$files[[i]]$raw_data <-list()
    out$files[[i]]$n_frames <- nrow(dat)
    for(j in colnames(dat)){
      out$files[[i]]$data[[j]] <- as.character(dat[,j])
      out$files[[i]]$raw_data[[j]] <- as.character(dat[,j])
    }
  }
  
  #add other variables
  out$label_names <- unique(unlist(lapply(out$files, FUN = function(x){names(x$data)})))
  out$meta <- data.frame(filename = names(out$files))
  out$smoothing <- 0
  rownames(out$meta) <- out$meta$filename
  out$file_names <- names(out$files)
  out$has.updated.analysis <- list(Transmat = FALSE, Onset = FALSE, Report = FALSE, Confmat = FALSE)
  return(out)
}

#' Will export all label data as individual .cvs files (on a per sample basis) into an output folder
#' 
#' @param us an object of type USData
#' @param path character. the path to the output folder where files will be saved. if nonexistent the folder will be created
#' @param raw.data a boolean. indicates if raw data (TRUE) or processed data (FALSE) (i.e smoothed data) should be exported. defaults to TRUE
#' @param sep character. the delimiter to be used in the csv file. defaults to ";"
#' @return files saved to the indicated folder
#' @examples
#' ExportAsCSVs(US, "C:/CSVfolder")
#' ExportAsCSVs(US, "C:/CSVfolder", raw.data = FALSE, sep = ",")
ExportAsCSVs <- function(us, path, raw.data = TRUE, sep = ";"){
  #check if directory exists, otherwise create it
  if(!dir.exists(path)){
    print(paste("directory does not exist. creating directory:", path))
    print(paste("current working directory:", getwd()))
    dir.create(path)
  }
  
  print(paste("saving data to directory: ", path, "/", sep = ""))
  for(i in us$file_names){
    if(raw.data){
      dat <- do.call(cbind,us$files[[i]]$raw_data)
    }else{
      dat <- do.call(cbind,us$files[[i]]$data)
    }
    f <- ifelse(grepl(".csv",i), i, paste(i,".csv",sep = ""))
    print(paste("saving file:",i))
    print(paste("saving to:", paste(path,f, sep = "/")))
    write.table(dat,file = paste(path,f, sep = "/"),sep = sep,col.names = TRUE, row.names = FALSE)
  }
}

#' Reads labels from a DLC Tracking object and returns an object of type USData
#' 
#' @param ts list of objects of type TrackingData (see DLCAnalyzer)
#' @return An object of type USData
#' @examples
#' US <- LoadFromTrackingObject(ts)
LoadFromTrackingObject <- function(ts){
  out <- list()
  out$files <- list()
  
  #add data
  for(i in names(ts)){
    if(!is.null(ts[[i]]$Report)){
      if(is.null(out$Report)){
        out$Report <- list()
        out$Report$raw <- data.frame()
      }
      out$Report$raw <- plyr::rbind.fill(out$Report$raw, data.frame(ts[[i]]$Report))
    }
    
    out$files[[i]] <- list()
    out$files[[i]]$data <- ts[[i]]$labels
    out$files[[i]]$raw_data <- ts[[i]]$labels
    out$files[[i]]$n_frames <- max(sapply(ts[[i]]$labels, FUN = function(x){length(x)}))
  }
  if(!is.null(out$Report)){
    rownames(out$Report$raw) <- names(out$files)
  }
  
  out$label_names <- unique(unlist(lapply(out$files, FUN = function(x){names(x$data)})))
  out$meta <- data.frame(filename = names(ts))
  out$smoothing <- 0
  rownames(out$meta) <- out$meta$filename
  out$file_names <- names(out$files)
  out$has.updated.analysis <- list(Transmat = FALSE, Onset = FALSE, Report = FALSE, Confmat = FALSE)
  return(out)
}

#' Smooths labels of an object of Type USData
#' 
#' @param us an object of type USData
#' @param integration_period the integration period (in frames) over which the data should be smoothed
#' @return An object of type USData
#' @examples
#' US <- SmoothLabels_US(us = US,integration_period = 5)
SmoothLabels_US <- function(us, integration_period){
  for(i in names(us$files)){
    for(j in names(us$files[[i]]$raw_data)){
      print(paste("smoothing",i,j,sep = " "))
      us$files[[i]]$data[[j]] <- SmoothLabel(us$files[[i]]$raw_data[[j]], integration_period)
    }
  }
  us$smoothing <- integration_period
  return(us)
}

#' Adds labels from a folder containing one or multiple CSVs with label data to an existing USData object. rows correspond to frames of the recording and column(s) to label(s). for each sample a separate file is needed that corresponds to the filename present in the USData object. headers of the columns correspond to label names used downstream. 
#' 
#' @param path a character with the path to the folder containing the CSVs
#' @param sep character. the separator used in CSV files (defaults to german system ";"). for english system use "," instead.
#' @return An object of type USData
#' @examples
#' US <- LoadFromCSVs("C:/CSVfolder_Newlabels")
AddFromCSVs <- function(us, path, sep = ";"){
  files <- list.files(path)
  backup <- us
  
  lab_names <- us$label_names
  for(i in us$file_names){
    if(length(grep(i, x = files)) == 1){
      print(paste("found file:",i, "in folder, adding labels" , sep = " "))
      dat <- read.table(paste(path, grep(i, x = files,value = TRUE), sep = "/"), sep = sep, header = T)
      lab_names <- append(lab_names,colnames(dat))
      for(j in colnames(dat)){
        us$files[[i]]$data[[j]] <- as.character(dat[,j])
        us$files[[i]]$raw_data[[j]] <- as.character(dat[,j])
      }
    }else{
      warning(paste("file:",i, "not found in folder, no data will be added!\n" , sep = " "))
    }
  }
  us$label_names <- unique(lab_names)
  
  # Add NA vectors for labels not present in certain files
  print("checking for data completeness")
  for(i in names(us$files)){
    for(j in us$label_names){
      if(is.null(us$files[[i]]$data[[j]])){
        print(paste("For file",i, "label",j ,"is missing. adding NA vector", sep = " "))
        us$files[[i]]$data[[j]] <- rep(NA,us$files[[i]]$n_frames)
        us$files[[i]]$raw_data[[j]] <- rep(NA,us$files[[i]]$n_frames)
      }
    }
  }
  
  #data integrity check. if failed return original object and print integrity report
  print("checking for data integrity")
  if(USDataCheck(us)){
    print("data integrity ok")
    return(us)
  }else{
    warning("data integrity compromised, returing original dataframe\n")
    return(backup)
  }
}


#' Adds data from a list of objects of Type TrackingData to an object of Type USData. there are 2 modes that are supported. 
#' if the list contains new files (not present in USData), these will be added to the USdata. 
#' For files already present in USData, labels that are present in the ts but not the USData will be added to the USData
#' 
#' @param us an object of type USData
#' @param ts list of objects of type TrackingData
#' @return An object of type USData
#' @examples
#' US <- AddFromTrackingObject(us = US,ts = TS)
AddFromTrackingObject <- function(us, ts){
  require(plyr)
  for(i in names(ts)){
    if(i %in% names(us$files)){
      for(j in names(ts[[i]]$labels)){
        if(is.null(us$files[[i]]$raw_data[[j]])){
          print(paste("adding label:",j,"to file:", i, sep = " "))
          us$files[[i]]$raw_data[[j]] <- ts[[i]]$labels[[j]]
          us$files[[i]]$data[[j]] <- ts[[i]]$labels[[j]]
        }
      }
    }else{
        print(paste("adding data for file:", i, sep = " "))
      us$files[[i]]$data <- ts[[i]]$labels
      us$files[[i]]$raw_data <- ts[[i]]$labels
      us$meta <- plyr::rbind.fill(us$meta,data.frame(filename = i))
    }
  }
  us$label_names <- unique(unlist(lapply(us$files, FUN = function(x){names(x$data)})))
  us$file_names <- names(us$files)
  rownames(us$meta) <- us$meta$filename
  us$has.updated.analysis <- list(Transmat = FALSE, Onset = FALSE, Report = FALSE, Confmat = FALSE)
  return(us)
}

  
#' sums a numeric vector x over a window (integer)
#' 
#' @param x a numeric or logical vector
#' @param window an integer that determines the size of the window
#' @return a numeric vector
#' @examples
#' x <- periodsum(x, window = 5)
periodsum <- function(x, window){
  res <- rep(0, length(x))
  for(i in 1:length(x)){
    res[i] <- sum(x[max(0,i-window):min(length(x), i + window)], na.rm = T)
  }
  return(res)
}

#' smooths a character vector over a integration period. will determine which entry is most abundant over a window of -integration_period to +integration_period and set the element the window is centered on to this
#' 
#' @param x a character or numeric vector
#' @param integration_period an integer that determines the size of the integration period to be used (+- integration_period centerd on the current element will be smoothed)
#' @return a character or numeric vector
#' @examples
#' x <- SmoothLabel(x, integration_period = 5)
SmoothLabel <- function(x, integration_period){
  types <-unique(x)
  mat <- NULL
  for (i in types){
    mat <- cbind(mat,periodsum(x == i, integration_period))
  }
  c <- apply(mat,1,FUN = which.max)
  return(types[c])
}


#' Adds BSOID data from a folder to an existing object of type USData. Will first determine if there exists a file with a comparable name for each sample in USData in the folder and then add it to the USData. 
#' if no file is found a NA vector will be added to this sample.
#' 
#' @param us an object of type USData
#' @param path a character. the path to the folder containing all the BSOID results
#' @param lab_name a character. name that will be used for the label in USData. Defaults to "BSOID"
#' @param cut_start an integer. defines if the first cut_start (integer) frames should be removed from BSOID
#' @return An object of type USData
#' @examples
#' US <- AddFromBsoid(us = US,path = "C:/PATHTOBSOIDFOLDER", lab_name = "BSOID", cut_start = 30)
AddFromBsoid <- function(us, path, lab_name = "BSOID", cut_start = 0){
  files <- list.files(path)
  
  for(i in us$file_names){
    us$files[[i]]$raw_data[[lab_name]] <- rep(NA,us$files[[i]]$n_frames)
    us$files[[i]]$data[[lab_name]] <- rep(NA,us$files[[i]]$n_frames)
    range_orig <- 1:us$files[[i]]$n_frames
    
    if(length(grep(i, x = files)) == 1){
      print(paste("found file:",i, "in folder, adding" , sep = " "))
      dat <- read.table(paste(path, grep(i, x = files,value = TRUE), sep = "/"), sep = ",", header = T)
      for(j in 1:nrow(dat)){
        dat$B.SOiD.labels <- as.character(dat$B.SOiD.labels)
        range <- (dat[j,]$Start.time..frames. + 1 - cut_start) : (dat[j,]$Start.time..frames. - cut_start + dat[j,]$Run.lengths) 
        range <- range[range %in% range_orig]
        if(length(range)!= 0){
        us$files[[i]]$raw_data[[lab_name]][range] <- dat[j,"B.SOiD.labels"]
        us$files[[i]]$data[[lab_name]][range] <- dat[j,"B.SOiD.labels"]
        }
      }
    }else{
      warning(paste("file:",i, "not found in folder, adding NA vector\n" , sep = " "))
    }
  }
  us$label_names <- append(us$label_names, lab_name)
  return(us)
}



#' Calculates a number of metrics for an object of type USData. these include a label usage report (time and occurences) as well all onset frames for each occurence
#' 
#' @param us an object of type USData
#' @return An object of type USData
#' @examples
#' US <- CalculateMetrics(US)
CalculateMetrics <- function(us){
  require(plyr)
  for(i in names(us$files)){
    d <- us$files[[i]]
    d$Report <- list()
      for(j in names(d$data)){
        d$Report[[j]] <- list()
        c <- na.omit(d$data[[j]])
        for(k in unique(c)){
        c2 <- c == k
        d$Report[[j]][[paste(k,"nframes", sep = ".")]] <- sum(c2)
        d$Report[[j]][[paste(k,"count", sep = ".")]] <- ceiling(sum((c2[2:length(c2)]!= c2[1:length(c2)-1])) / 2)
      }
      }
    us$files[[i]] <- d
  }
  if(is.null(us$Report)){
  us$Report <- list()
  }
  for(j in us$label_names){
    us$Report[[j]] <- data.frame()
  }
  for(i in names(us$files)){
    for(j in names(us$files[[i]]$Report)){
    add <- data.frame(file = i,us$files[[i]]$Report[[j]])
    names(add) <- c("file",names(us$files[[i]]$Report[[j]]))
    us$Report[[j]] <- plyr::rbind.fill(us$Report[[j]], data.frame(add))
    }
  }
  for(j in us$label_names){
    rownames(us$Report[[j]]) <- us$Report[[j]]$file
    us$Report[[j]]$file <- NULL
  }
  us <- CalculateOnSet(us)
  us$has.updated.analysis$Report <- TRUE
  return(us)
}


#' Adds TransitionMatrix data to an object of type USData
#' 
#' @param us an object of type USData
#' @param labels a list of characters. names of labels for which a transitionmatrix should be calculated. defaults to NULL which will do all labels
#' @return An object of type USData
#' @examples
#' US <- AddTransitionMatrixData(US,labels = c("kmeans.25","BSOID.27"))
AddTransitionMatrixData <- function(us, labels = NULL){
  if(!USDataCheck(us)){
    warning("Data check failed! see warning message for more information, retruning original data\n")
    return(us)
  }
  
  if(is.null(labels)){
    labels <- us$label_names
    print("no label specified, performing transition matrix calculations for all")
  }
  
  for(label in labels){
    lab <- NULL
    for(i in us$files){
      lab <- append(lab,i$data[[label]])
    }
      lev <- unique(na.omit(lab))
  
    #create ordered and complete transition matrix for each files
    for(i in names(us$files)){
      lab <- us$files[[i]]$data[[label]]

      transvector <- na.omit(lab[lab[2:length(lab)] != lab[1:(length(lab)-1)]])
      transmatrix <- table(transvector[1:(length(transvector)-1)],transvector[2:length(transvector)])
      for(j in lev){
        if(!(j %in% rownames(transmatrix))){
          transmatrix <- rbind(transmatrix, rep(0,ncol(transmatrix)))
          rownames(transmatrix)[nrow(transmatrix)] <- j
        }
      }
      for(j in lev){
        if(!(j %in% colnames(transmatrix))){
          transmatrix <- cbind(transmatrix, rep(0,nrow(transmatrix)))
          colnames(transmatrix)[ncol(transmatrix)] <- j
        }
      }
        transmatrix <- transmatrix[lev,lev]
        
        if(is.null(us$files[[i]]$transmatrix)){
          us$files[[i]]$transmatrix <- list()
          us$files[[i]]$transmatrix_norm <- list()
        }
        us$files[[i]]$transmatrix[[label]] <- transmatrix
        us$files[[i]]$transmatrix_norm[[label]] <- t(apply(transmatrix, MARGIN = 1, FUN = function(x){if(sum(x) > 0){x / sum(x)}else{x}}))
              }
  }
  us$has.updated.analysis$Transmat <- TRUE
  us$transitionmatrix_names <- labels
  return(us)
}

#' Drops one or multiple labels from an object of type USdata
#' 
#' @param us an object of type USData
#' @param labels a list of characters or single character. names of labels that should be dropped from the object
#' @return An object of type USData
#' @examples
#' US <- DropLabels(US,labels = c("kmeans.50"))
DropLabels <- function(us, labels){
  keep = us$label_names[!(us$label_names %in% labels)]
  
  for(i in us$file_names){
    for(j in names(us$files[[i]])){
      if(is.list(us$files[[i]][[j]])){
        us$files[[i]][[j]] <- us$files[[i]][[j]][keep]
      }
    }
  }
  if(us$has.updated.analysis$Report){
    us$Report <- us$Report[!(names(us$Report) %in% labels)]
  }
  if(us$has.updated.analysis$Transmat){
    us$transitionmatrix_names <- us$transitionmatrix_names[!(us$transitionmatrix_names %in% labels)]
  }
  us$label_names <- keep
  return(us)
}

#' Fuses two objects of type USData. will only keep labels for which it has an entry in each sample in the Fuse USdata object. can only fuse along the sample axis (no new labels for existing samples can be added!)
#' 
#' @param us1 an object of type USData
#' @param us2 an object of type USData
#' @return An object of type USData
#' @examples
#' US <- FuseUSData(US,US2)
FuseUSData <- function(us1, us2){
  require(plyr)
  if(sum(us1$file_names %in% us2$file_names) != 0){
    print("warning: some samples are present in both USdata objects, can not fuse. returning original us1 without fusing\n")
    return(us1)
  }
  
  #build a boolean matrix that indicates if data for a sample - label pair is present in the fused data. row = samples, col = labels
  fn <- unique(c(us1$file_names,us2$file_names))
  ln <- unique(c(us1$label_names,us2$label_names))
  mat <- matrix(data = FALSE, nrow = length(fn), ncol = length(ln),dimnames = list(fn,ln))
  mat[us1$file_names,us1$label_names] <- TRUE
  mat[us2$file_names,us2$label_names] <- TRUE
  
  #test which of the labels are present across all samples after fusing
  keep <- apply(mat, MARGIN = 2, FUN = function(x){as.logical(min(x))})
  
  #dropping labels that are not present in all files
  if(sum(keep) != length(keep)){
    print(paste("detected labels with no entries in some files, dropping labels:", names(keep)[!keep], sep = " "))
    us1 <- DropLabels(us1,names(keep)[!keep])
    us2 <- DropLabels(us2,names(keep)[!keep])
  }
  keep <- names(keep)[keep]
  
  #adding files
  for(i in names(us2$files)){
    us1$files[[i]] <- us2$files[[i]]
  }
  
  us1$label_names <- keep
  us1$file_names <- names(us1$files)
  
  if(us1$has.updated.analysis$Report){
    print("re-runing combined report")
    us1 <- CalculateMetrics(us1)
  }
  
  #combine raw report if present in both us1 and us2
  if(!is.null(us1$Report$raw) & !is.null(us2$Report$raw)){
    us1$Report$raw <- plyr::rbind.fill(us1$Report$raw, us2$Report$raw)
  }
  
  #combine metadata of us1 and us2
  rname <- append(rownames(us1$meta),rownames(us2$meta))
  us1$meta <- plyr::rbind.fill(us1$meta, us2$meta)
  rownames(us1$meta) <- rname
  
  #remove confusion matrix if present
  if(!is.null(us1$ConfusionMatrix)){
    print("deleting confusion matrix")
    us1$ConfusionMatrix <- NULL
    us1$ConfusionMatrix_norm <- NULL
    us1$has.updated.analysis$Confmat <- FALSE
  }
  
  #re-run transitionmatrix calculation in case us1 or us2 had transitionmatrix data (otherwise row and colnames of transmatrix may be inconsistent)
  if(us1$has.updated.analysis$Transmat | us2$has.updated.analysis$Transmat){
    print("re-runing transitionmatrix calculations")
    us1 <- AddTransitionMatrixData(us1)
  }
  
  #print warning if us1 and us2 had unequal smoothing parameters
  if(us1$smoothing != us2$smoothing){
    warning(paste("datasets have unequal smoothing. us1:", us1$smoothing, "us2:", us2$smoothing,"! highly keep recommended to rerun smoothing and other down stream analyses!\n", sep = " "))
  }
  
  return(us1)
}


#' Adds stabilized transitions to an object of type USData. (either for all data, or for a selected subset) 
#' 
#' @param us an object of type USData
#' @param group a character. determines which column (named group) from us$meta should be used for stabilization
#' @param control a character. determines the value that should be used as control group for us$meta$group
#' @param labels a list of characters. names of labels for which a stabilized transitions should be calculated. defaults to NULL which will do all labels
#' @param subset a list of characters. names of files that define the subset for which the stabilization should be performed
#' @return An object of type USData
#' @examples
#' US <- CalculateStabilizedTransitions(US, group = "Treatment", control = "Control")
#' US <- CalculateStabilizedTransitions(US, group = "Condition", control = "Control",labels = c("kmeans.25","BSOID.27"), subset = c("File1.csv","File2.csv","File3.csv"))
CalculateStabilizedTransitions <- function(us, group, control, labels = NULL, subset = NULL){
  if(!USDataCheck(us)){
    warning("Data check failed! see warning message for more information, retruning original data\n")
    return(us)
  }
  
  if(is.null(labels)){
    labels <- us$label_names
    print("no label specified, performing transition matrix stabilization for all")
  }
  if(is.null(subset)){
    subset <- names(us$files)
  }
  if(!(group %in% colnames(us$meta))){
    stop("invalid grouping variable. not contained in metadata")
  }
  
  grp <- us$meta[subset, group]
  
  if(!(control %in% grp)){
    stop("invalid value for the control group")
  }
  
  for(label in labels){
    mats <- list()
    for(i in names(us$files[subset])){
      mats[[paste(i)]] <- us$files[[i]]$transmatrix[[label]]
    }
    cmat <- mats[grp == control]
    cmat <- Reduce('+', cmat) / sum(grp == control)
    mats <- lapply(mats,FUN = function(x){x - cmat})
    for(i in names(us$files[subset])){
      if(is.null(us$files[[i]]$transmatrix_stabilized)){
        us$files[[i]]$transmatrix_stabilized <- list()
      }
      us$files[[i]]$transmatrix_stabilized[[label]] <- mats[[i]]
    }
  }
  return(us)
}

#' Adds confusion matrix (or many confusion matrices) to an object of type USData
#' 
#' @param us an object of type USData
#' @param files a character or list of characters. determines which files should be included into the confusion matrix calculation. defaults to NULL = all files
#' @param from a character or list of characters. determines from which labels confusion should be checked. defaults to NULL = all labels
#' @param to a character or list of characters. determines to which labels confusion should be checked. defaults to NULL = all labels
#' @return An object of type USData
#' @examples
#' US <- AddConfusionMatrix(US)
#' US <- AddConfusionMatrix(US, files = c("File1.csv","File2.csv","File3.csv"), from = "kmeans.25", to = "kmeans.50")
#' US <- AddConfusionMatrix(US, from = c("kmeans.25","kmeans.50"), to = c("kmeans.25","kmeans.50"))
AddConfusionMatrix <- function(us, files = NULL, from = NULL, to = NULL){
  if(!USDataCheck(us)){
    warning("Data check failed! see warning message for more information, retruning original data\n")
    return(us)
  }
  
  if(is.null(files)){
    files = names(us$files)
    print("no files specified, runing confusion matrix for all files")
  }
  if(is.null(from)){
    from <- us$label_names
  }
  if(is.null(to)){
    to <- us$label_names
  }
  
  for(f in from){
    for(t in to[!to %in% f]){
        name <- paste(f, t, sep = "-")
        tr <- NULL
        comp <- NULL
        for(i in us$files[files]){
          tr <- append(tr,i$data[[f]])
          comp <- append(comp,i$data[[t]])
        }
  
      Results <- table(from = tr, to = comp)
    if(is.null(us$ConfusionMatrix)){
      us$ConfusionMatrix <- list()
      us$ConfusionMatrix_norm <- list()
      }
      us$ConfusionMatrix[[name]] <- Results
      us$ConfusionMatrix_norm[[name]] <-  t(apply(Results,MARGIN = 1, FUN = function(x){x / sum(x)}))
    }
  }
  us$has.updated.analysis$Confmat <- TRUE
  return(us)
}

#' Plots a confusion matrix
#' 
#' @param mat a confusion matrix as produced by the function AddConfusionMatrix()
#' @param pointsize a numeric value. defines the size of the points
#' @param hclust boolean. decides if hierarchical clustering on rows and columns should be used (TRUE) or not (FALSE). defauts to TRUE
#' @param sort boolean. decides if rows and columns should be sorted by cluster number value (small to large). defauts to TRUE
#' @return A plot
#' @examples
#' PlotConfusionMatrix(US$ConfusionMatrix$`rear.classifier-kmeans.25`)
#' PlotConfusionMatrix(US$ConfusionMatrix$`rear.classifier-kmeans.25`, hclust = FALSE, sort = TRUE)
#' PlotConfusionMatrix(US$ConfusionMatrix$`rear.classifier-kmeans.25`, pointsize = 5)
PlotConfusionMatrix <- function(mat, pointsize = 10, hclust = T, sort = F){
  require(ggplot2)
  plotdat <- NULL
  for(i in rownames(mat)){
    for(j in colnames(mat)){
      plotdat <- rbind(plotdat, data.frame(from = i, to = j, value = mat[i,j]))
    }
  }
  plotdat$from <-factor(plotdat$from)
  plotdat$to <- factor(plotdat$to)
  
  if(hclust){
    r <- hclust(dist(mat))
    c <- hclust(dist(t(mat)))
    plotdat$from <- factor(plotdat$from, levels = levels(plotdat$from)[r$order])
    plotdat$to <- factor(plotdat$to, levels = levels(plotdat$to)[c$order])
  }
  if(sort){
    r <- order(as.integer(levels(factor(plotdat$from))))
    c <- order(as.integer(levels(factor(plotdat$to))))
    plotdat$from <- factor(plotdat$from, levels = levels(plotdat$from)[r])
    plotdat$to <- factor(plotdat$to, levels = levels(plotdat$to)[c])
  }
  
  ggplot(plotdat, aes(from,to, fill = value)) + 
    scale_fill_gradient(low = "white", high = "red") +
    geom_point(shape = 22, size = pointsize * sqrt(plotdat$value), color = "black") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}


#' Splits an Object of Type USData to create a new USdata of a reduced dataset
#' 
#' @param us us an object of type USData
#' @param include a character or list of characters. names of the files to be included in the new subset
#' @param omit a character or list of characters. names of the files to be excluded in the new subset
#' @param select boolean vector. a vector of boolean values (TRUE or FALSE) that is exactly as long as the number of files of the original USData object. TRUE indexes corresponds to the us$file_names entries to keep
#' @return An object of type USData
#' @examples
#' US <- SplitUSData(US, include = c("File1.csv","File2.csv","File3.csv"))
#' US <- SplitUSData(US, omit = c("File1.csv","File2.csv","File3.csv"))
#' US <- SplitUSData(US, select = c(TRUE,TRUE,FALSE,FALSE,TRUE,TRUE))
#' US <- SplitUSData(US, select = US$meta$group == "Control")
SplitUSData <- function(us, include = NULL, omit = NULL, select = NULL){
  keep <- names(us$files)
  if(!is.null(select)){
    if(length(select) != length(us$files) | !is.logical(select)){
      warning("for select a logical vector has to be supplied that is exactly as long as the number of files in the data object. returning original data\n")
      return(us)
    }
    select[is.na(select)] <- FALSE
    keep <- names(us$files)[select]
  }
  
  if(!is.null(include)){
    if(sum(include %in% names(us$files)) == 0){
      warning("ensure that include is a vector of names to keep, non of the filenames in include are present in the data. returning original data\n")
      return(us)
    }
    keep <- include[include %in% names(us$files)]
  }
  if(!is.null(omit)){
    keep <- keep[!(keep %in% omit)]
  }
  
  us$files <- us$files[keep]
  us$meta <- us$meta[keep,,drop = FALSE]
  us$file_names <- keep
  if(!is.null(us$Report)){
    for(j in names(us$Report)){
      us$Report[[j]] <- us$Report[[j]][keep,]
    }
  }
  
  return(us)
}


#' Cuts each file in an object of type USdata.
#' 
#' @param us us an object of type USData
#' @param maxlength a integer. maximum length a file may be (in frames). will cut the end until length =< maxlength. defaults to NULL = no defined maximum length
#' @param start integer. number of frames that should be cut at the start of each file. defaults to NULL = no cut at start
#' @param end integer. number of frames that should be cut at the end of each file. defaults to NULL = no cut at end
#' @return An object of type USData
#' @examples
#' US <- CutUSdata(US, maxlength = 10000)
#' US <- CutUSdata(US, start = 500)
#' US <- CutUSdata(US, end = 500)
#' US <- CutUSdata(US, start = 500, end = 500, maxlength = 2000)
CutUSdata <- function(us, maxlength = NULL, start = NULL, end = NULL){
  for(i in us$file_names){
    dat <- us$files[[i]]
    keep <- 1:dat$n_frames
    if(!is.null(start)){
      keep <- start:dat$n_frames
    }
    if(!is.null(end)){
      keep <- keep[1]:(dat$n_frames - end)
    }
    if(!is.null(maxlength)){
      keep <- keep[1]:keep[min(maxlength,length(keep))]
    }
    for(j in names(dat$data)){
      dat$data[[j]] <- dat$data[[j]][keep]
      dat$raw_data[[j]] <- dat$raw_data[[j]][keep]
    }
    dat$n_frames <- length(keep)
    us$files[[i]] <- dat[1:3]
  }
  
  if(us$has.updated.analysis$Transmat){
    warning("deleted all Transitionmatrices due to changes in data. re-run analysis!\n")
    us$has.updated.analysis$Transmat <- FALSE
    us$transitionmatrix_names <- NULL
    us$has.updated.analysis$Transmat <- FALSE
  }
  if(us$has.updated.analysis$Confmat){
    warning("deleted all Confusionmatrices due to changes in data. re-run analysis!\n")
    us$has.updated.analysis$Confmat <- FALSE
    us$ConfusionMatrix <- NULL
    us$ConfusionMatrix_norm <- NULL
    us$has.updated.analysis$Confmat <- FALSE
  }
  if(us$has.updated.analysis$Report){
    warning("deleted all Metrics (Onsetdata and per file Reports) due to changes in data. re-run analysis!\n")
    us$Report[us$label_names] <- NULL
    us$has.updated.analysis$Report <- FALSE
    us$has.updated.analysis$Onset <- FALSE
  }
  
  return(us)
}




#' Runs an USdata integrity check and returns TRUE (passed) or FALSE (failed). will also produce warning messages to highlight where the problem is located
#' 
#' @param us an object of type USData
#' @return A boolean value
#' @examples
#' USDataCheck(US)
USDataCheck <- function(us){
  Report <- USDataReport(us)
  
  iseq <- function(x){
    y <- x[1]
    sanity <- TRUE
    for(i in x){
      if(i != y){
        sanity <- FALSE
      }
      y = i
    }
    return(sanity)
  }
  
  if(sum(is.na(Report))){
    warning("Data check failed: there are files with missing data\n")
    message(Report)
    return(FALSE)
  }
  eq_num <- apply(Report,1,FUN = iseq)
  if(sum(!eq_num)){
    warning("Data check failed: data associated with the following files are of different lengths!\n")
    for(i in rownames(Report[!eq_num,])){
      message(paste(i, ":"))
      for(j in colnames(Report[!eq_num,])){
        message(paste("N =", Report[i,j], "for label:",j, sep = " "))
      }
    }
    return(FALSE)
  }
  return(TRUE)
}


#' Creates a USData Report as table that contains files in rows and labels in columns. entries will indicate the length of the data contained in each sample-label pair
#' 
#' @param us an object of type USData
#' @return a data frame
#' @examples
#' USDataReport(US)
USDataReport <- function(us){
  require(plyr)
  Report <- NULL
  for(i in names(us$files)){
    Report <- plyr::rbind.fill(Report, data.frame(file = i, t(unlist(lapply(us$files[[i]]$data, FUN = length)))))
  }
  rownames(Report) <- Report$file
  Report$file <- NULL
  return(Report)
}

Onset <- function(x){
  which((x[2:(length(x))] - x[1:(length(x) -1)]) == 1) + 1
}
Offset <- function(x){
  which((x[2:(length(x))] - x[1:(length(x) -1)]) == -1)
}


#' Calculates onset frames for each label in each sample and adds this data to an object of type USData
#' 
#' @param us an object of type USData
#' @return an object of type USData
#' @examples
#' US <- CalculateOnSet(US)
CalculateOnSet <- function(us){
  for(i in names(us$files)){
    us$files[[i]]$OnSetData <- list()
    for(j in names(us$files[[i]]$data)){
      us$files[[i]]$OnSetData[[j]] <- list()
      lev <- na.omit(unique(us$files[[i]]$data[[j]]))
      for(k in lev){
        us$files[[i]]$OnSetData[[j]][[k]] <- Onset(us$files[[i]]$data[[j]] == k)
      }
    }
  }
  us$has.updated.analysis$Onset <- TRUE
  return(us)
}

#' Creates Behavior Train Plot for selected labels and values. these plots will select random behavior trains of the selected cluster of a label class and plot the corresponding clusters they map to in all other label classes
#' 
#' @param us an object of type USData
#' @param lab a character. selects the label class to be used
#' @param val a character. selects the label value to be used
#' @param len an integer. selects the maximum length to be plotted
#' @param n an integer. selects the number of random behavior trains to be plotted
#' @param max_clust an integer. selects the maximum number of clusters to be colored for each other label
#' @param file a character or list of characters. name(s) of files to be used for behavior train plots
#' @return a plot
#' @examples
#' US <- BehaviorTrainPlot(US, lab = "kmeans.25", value = "11")
#' US <- BehaviorTrainPlot(US, lab = "kmeans.25", value = "11", len = 50, n = 50, max_clust = 7)
#' US <- BehaviorTrainPlot(US, lab = "kmeans.25", value = "11", files = c("File1.csv","File2.csv","File3.csv"))
BehaviorTrainPlot <- function(us,lab,val, len = 100, n = 100, max_clust = 5, file = NULL){
  require(ggplot2)
  require(cowplot)
  if(!us$has.updated.analysis$Onset){
    print("no updated behavior Onset data present. Generating Onset data  first")
    us <- CalculateOnSet(us)
  }
  
  Plotdat <- data.frame()
  
  if(is.null(file)){
    file <- names(us$files)
  }
  
  #creates dataframe with that defines where cluster examples can be found (file and frame) and how long they are
  examples <- data.frame()
  for(fn in file){
    f <- us$files[[fn]]
    if(!is.null(f$OnSetData[[lab]][[val]])){
      on <- f$OnSetData[[lab]][[val]]
      on <- on[(on + len) < f$n_frames]
      examples <- rbind(examples, data.frame(file = fn, onset = on))
    }
  }
  #shuffle dataframe to randomize order
  examples <- examples[sample(1:nrow(examples)),]
  
  #selects n examples and creates the plotdat
  idx = 0
  for(i in 1:min(n,nrow(examples))){
    idx = idx + 1
    traindat <- data.frame()
    for(j in names(us$files[[examples[i,"file"]]]$data)){
      traindat <- rbind(traindat,data.frame(idx = rep(idx,len), idx2 = c(1:len), value = us$files[[examples[i,"file"]]]$data[[j]][examples[i,"onset"]:(examples[i,"onset"] + len - 1)], type = j))
    }
    Plotdat <- rbind(Plotdat,traindat)
  }
  
  Plotdat$keep <- paste(Plotdat$idx,Plotdat$idx2)
  keep <- Plotdat[Plotdat$type == lab & Plotdat$value == val,"keep"]
  Plotdat <- Plotdat[Plotdat$keep %in% keep,]
  
  ps <- list()
  ps_d <- list()
  for(i in unique(Plotdat$type)){
    pd <- Plotdat[Plotdat$type == i,]
    keep <- names(sort(table(pd$value), decreasing = T))[1:min(max_clust, length(table(pd$value)))]
    pd$value[!(pd$value %in% keep)] <- NA
    ps[[i]] <- ggplot(pd ,aes(idx2,idx, color = value)) + 
      geom_point(shape = 15) + 
      theme_bw() + 
      ggtitle(i) + 
      xlab("frame") +
      ylab("example") + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    ps_d[[i]] <- ggplot(pd ,aes(idx2, color = value)) + 
      geom_density() + 
      theme_bw() +
      xlab("frame") +
      ylab("density") + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
  }
  p1 <- cowplot::plot_grid(plotlist = append(ps[lab],ps[names(ps) != lab]), nrow = 1)
  p2 <- cowplot::plot_grid(plotlist = append(ps_d[lab],ps_d[names(ps) != lab]), nrow = 1)
  cowplot::plot_grid(p1,p2, ncol = 1, rel_heights = c(1,0.3))
}


#' Creates Behavior Flow plots for a selected label class of an object of type USData
#' 
#' @param us an object of type USData
#' @param file a character or list of characters. name(s) of files to be used for behavior flow plot
#' @param lab a character. selects the label class for which behavior flow should be plotted
#' @return a plot
#' @examples
#' US <- PlotBehaviorFlow(US, lab = "kmeans.25")
#' US <- PlotBehaviorFlow(US, file = c("File1.csv","File2.csv","File3.csv") ,lab = "kmeans.25")
PlotBehaviorFlow <- function(us,file = NULL, lab){
  require(circlize)
  if(!us$has.updated.analysis$Transmat){
    print("no updated Tranisition data present. Generating Transition matrix data first")
    us <- AddTransitionMatrixData(us)
  }
  
  
  if(is.null(file)){
    file <- names(us$files)
  }
  
  #build summed transition matrix
  mat <- 0
  for(i in file){
    mat <- mat + us$files[[i]]$transmatrix[[lab]]
  }
  len <- nrow(mat)
  #colors with good contrast. note, they will start to repeat after the 25 entries supplied
  contrastcols <- c(
    "dodgerblue2", "#E31A1C", # red
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "black", "gold1",
    "skyblue2", "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown")
  cols <- contrastcols[(1:len %% length(contrastcols) + 1)]
  names(cols) <- rownames(mat)
  
  #plot chord
  circlize::chordDiagramFromMatrix(mat,grid.col = cols)
}

#' Creates Behavior Flow delta plots for a selected label class of an object of type USData
#' 
#' @param us an object of type USData
#' @param grouping a boolean vector that has exactly as long as the number of samples in USdata. indicates if a sample is part of the reference (FALSE) or target (TRUE) set
#' @param lab a character. selects the label class for which behavior flow delta should be plotted
#' @param method type of behavior flow to be plotted. "absolute" shows absolute value of up and downregulated clusters. "up" plots only clusters that are increased in the target set compared to the reference set. "down" plots only clusters that are decreased in the target set compared to the reference set. defaults to "absolute"
#' @return a plot
#' @examples
#' US <- PlotBehaviorFlow_Delta(US, grouping = c(TRUE,TRUE,TRUE,FALSE,FALSE,FALSE),lab = "kmeans.25")
#' US <- PlotBehaviorFlow_Delta(US, grouping = US$meta$group == "Test" ,lab = "kmeans.25")
#' US <- PlotBehaviorFlow_Delta(US, grouping = US$meta$group == "Test" ,lab = "kmeans.25", method = "up")
#' US <- PlotBehaviorFlow_Delta(US, grouping = US$meta$group == "Test" ,lab = "kmeans.25", method = "down")
PlotBehaviorFlow_Delta <- function(us,grouping, lab, method = "absolute"){
  require(circlize)
  if(!us$has.updated.analysis$Transmat){
    print("no updated Tranisition data present. Generating Transition matrix data first")
    us <- AddTransitionMatrixData(us)
  }
  if(length(us$file_names) != length(grouping) | !is.logical(grouping)){
    warning("Improper grouping vector. needs to be a logical (True/False) vector with one entry per sample\n")
    return(NULL)
  }

    
    
  #define groups
  g1 <- us$file_names[grouping]
  g2 <- us$file_names[!grouping]

  #build avearged transition matrix for group1 and group2
  mat1 <- 0
  mat2 <- 0
  for(i in g1){
    mat1 <- mat1 + us$files[[i]]$transmatrix[[lab]]
  }
  for(i in g2){
    mat2 <- mat2 + us$files[[i]]$transmatrix[[lab]]
  }
  mat1 <- mat1 / length(g1)
  mat2 <- mat2 / length(g2)
  
  len <- nrow(mat1)
  #colors with good contrast. note, they will start to repeat after the 25 entries supplied
  contrastcols <- c(
    "dodgerblue2", "#E31A1C", # red
    "green4",
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "black", "gold1",
    "skyblue2", "#FB9A99", # lt pink
    "palegreen2",
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "gray70", "khaki2",
    "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
    "darkturquoise", "green1", "yellow4", "yellow3",
    "darkorange4", "brown")
  cols <- contrastcols[(1:len %% length(contrastcols) + 1)]
  names(cols) <- rownames(mat1)
  
  #select method
  if(method == "up"){
    mat <- mat1 - mat2
    mat[mat < 0] <- 0
  }else if(method == "down"){
    mat <- mat2 - mat1
    mat[mat < 0] <- 0
  }else{
    mat <- abs(mat1-mat2)
  }
  
  #plot chord
  circlize::chordDiagramFromMatrix(mat,grid.col = cols)
}

#' Adds metadata from a data frame to an object of type USdata
#' 
#' @param us an object of type USData
#' @param s2c a data frame. sample names as rows and metadata variables as columns. needs to contain one row for each sample in the USdata object and the rowname needs to correspond to the sample name used in USData (US$file_names)
#' @return a plot
#' @examples
#' metadata <- data.frame(group = c("Control","Control","Test","Test"), sex = c("male","female","male","female"))
#' rownames(metadata) <- c("File1.csv","File2.csv","File3.csv","File4.csv)
#' US <- AddMetaData(US, s2c = metadata)
AddMetaData <- function(us,s2c){
  if(sum(!(rownames(us$meta) %in% rownames(s2c)))){
    warning("Can not add ,metadata. Incomplete list, following files are missing:\n")
    print((rownames(us$meta)[!(rownames(us$meta) %in% rownames(s2c))]))
    return(us)
  }
  assign <- match(rownames(us$meta), rownames(s2c))
  us$meta <- cbind(us$meta, s2c[assign,])
  return(us)
}


#' Runs a two group analysis on an object of type USData across raw data, cluster usage, transition usage and transition matrix differences. the USdata object needs to have only samples of the two groups to compare (use function SplitUSData() first if there are more samples in USData that are irrelevant for this specific analysis)
#' 
#' @param us an object of type USData
#' @param group a vector of characters. this vector needs to have an entry for each sample within the USdata and exactly two groups
#' @param name a character. the name under which the analysis will be saved. defaults to NULL which will automatically generate a name
#' @param n_bootstraps an integer. number of bootstraps that will be used to calculate transition matrix distance statistics. larger number will enable a better modeling of the null distribution at the cost of increased runtime. defaults to 1000
#' @param lab character or list of characters. the label group(s) for which the two group analysis will be performed. defaults to NULL which will do an independent analysis for each label group
#' @return an object of type USData
#' @examples
#' US <- TwoGroupAnalysis(US, group = c("control","control","control","test","test","test"))
#' US <- TwoGroupAnalysis(US, group = US$meta$group)
#' US <- TwoGroupAnalysis(US, group = US$meta$group, name = "Analysis_group_10000bootstraps", n_bootstraps = 10000, lab = "kmeans.25"))
#' print(US$Results$`Analysis_group_10000bootstraps`$Statistics)
#' print(US$Results$`Analysis_group_10000bootstraps`$TransitionStats)
TwoGroupAnalysis <- function(us,group, name = NULL, n_bootstraps = 1000, lab = NULL){
  require(imputeTS)
  require(prcma)
  if(length(us$files) != length(group)){
    warning("grouping vectorn needs to have same length as number of files. can not perform analysis\n")
    return(us)
  }
  if(length(unique(group)) != 2){
    warning("there need to be exactly two groups in the grouping vector. can not perform analysis\n")
    return(us)
  }
  if(is.null(name)){
    name <- paste(unique(group), collapse = "-vs-")
  }
  if(is.null(lab)){
    print("no label specified, performing analysis for all")
    lab <- us$label_names
  }
  
  if(is.null(us$Results)){
    us$Results <- list()
  }
  us$Results[[name]] <- list()
  
  # t.test based statistics on Report data (i.e raw data and cluster usage data)
  us$Results[[name]]$Statistics <- list()
  for(r in names(us$Report)){
    Report <- imputeTS::na_replace(us$Report[[r]],fill = 0)
    res <- NULL
    for(i in names(Report)){
      a <- Report[group == unique(group)[1],i]
      b <- Report[group == unique(group)[2],i]
      nline <- data.frame(i, t.test(a,b)$p.value, mean(a), mean(b), mean(a) / mean(b))
      names(nline) <- c("name","p",unique(group)[1], unique(group)[2], "FC")
      res <- rbind(res,nline)
    }
    res$FDR <- p.adjust(res$p, method = "BY")
    us$Results[[name]]$Statistics[[r]] <- res[order(res$p),]
  }
  
  us$Results[[name]]$TransitionStats <- list()
  for(r in lab){
    transmat <- list()
  for(i in names(us$files)){
    transmat[[i]] <- us$files[[i]]$transmatrix[[r]]
  }
  
    # bootstrapping based analysis of inter group distance of mean transition matrix values
  A <- transmat[group == unique(group)[1]]
  B  <- transmat[group == unique(group)[2]]
  
  a <- do.call(cbind, A)
  a <- array(a, dim=c(dim(A[[1]]), length(A)))
  a <- apply(a, c(1, 2), mean, na.rm = TRUE)
  
  b <- do.call(cbind, B)
  b <- array(b, dim=c(dim(B[[1]]), length(B)))
  b <- apply(b, c(1, 2), mean, na.rm = TRUE)
  distance <- sum(abs(a-b))
  
  bootstraps <- NULL
  for(i in 1:n_bootstraps){
    grp <- sample(group)
    A <- transmat[grp == unique(group)[1]]
    B  <- transmat[grp == unique(group)[2]]
    a <- do.call(cbind, A)
    a <- array(a, dim=c(dim(A[[1]]), length(A)))
    a <- apply(a, c(1, 2), mean, na.rm = TRUE)
    b <- do.call(cbind, B)
    b <- array(b, dim=c(dim(B[[1]]), length(B)))
    b <- apply(b, c(1, 2), mean, na.rm = TRUE)
    bootstraps <- append(bootstraps,sum(abs(a-b)))
  }
  
  # t.test based statistics on Transition data (i.e do specific transitions incerase or decrease in test vs control?)
  lev <- rownames(do.call(cbind, A))
  transitions <- NULL
  for(i in lev){
    for(j in lev){
      A <- lapply(transmat[group == unique(group)[1]],FUN = function(x){x[i,j]})
      B <- lapply(transmat[group == unique(group)[2]],FUN = function(x){x[i,j]})
      tryCatch({
      test <- t.test(unlist(A),unlist(B))
      transitions <- rbind(transitions, data.frame(from = i, to = j, p.value = test$p.value, A = mean(unlist(A)), B = mean(unlist(B)) ,FC = mean(unlist(A)) / mean(unlist(B))))
      }, error = function(e){
      transitions <- rbind(transitions, data.frame(from = i, to = j, p.value = NA, A = mean(unlist(A)), B = mean(unlist(B)) ,FC = NA))
      })
      }
  }
  names(transitions)[4:5] <- c(paste(unique(group)[1]),paste(unique(group)[2]))
  transitions[is.nan(transitions$p.value),c("p.value","FC")] <- c(NA,NA)
  transitions$FDR <- p.adjust(transitions$p.value, method = "BY")
  transitions <- transitions[order(transitions$p.value),]
  
  # Add Results to USData
  us$Results[[name]]$TransitionStats[[r]]$distance <- distance
  us$Results[[name]]$TransitionStats[[r]]$bootstraps <- bootstraps
  us$Results[[name]]$TransitionStats[[r]]$percentile <- sum(bootstraps < distance) / length(c(bootstraps,distance)) * 100
  us$Results[[name]]$TransitionStats[[r]]$sigma <- (distance - mean(bootstraps)) / sd(bootstraps)
  us$Results[[name]]$TransitionStats[[r]]$p.value <- (1 - pracma::erf(us$Results[[name]]$TransitionStats[[r]]$sigma / sqrt(2))) / 2
  us$Results[[name]]$TransitionStats[[r]]$transitions <- transitions
  }
  
  return(us)
}


#' Runs a PCA analysis on a Report group followed by a linear model analysis on PCs for group effect
#' 
#' @param us an object of type USData
#' @param group a vector of characters. this vector needs to have an entry for each sample within the USdata and exactly two groups
#' @param cumulative_proportion numeric. decides to what cumulative proportion explain variability top PCs are considere. defaults to 0.9 (=90%)
#' @param normalize boolean should data be normalized before PCA? defaults to TRUE
#' @param lab character. determines which US$Report element should be used for PCA. defaults to "raw"
#' @return a statistical report
#' @examples
#' US <- TwoGroupAnalysis_PCA_Report(US, group = c("control","control","control","test","test","test"))
#' US <- TwoGroupAnalysis_PCA_Report(US, group = US$meta$Condition, cumulative_proportion = 0.8, normalize = FALSE, lab = "kmeans.25")

TwoGroupAnalysis_PCA_Report <- function(us,group, cumulative_proportion = 0.9, normalize = TRUE, lab = "raw"){
  if(length(us$files) != length(group)){
    warning("grouping vectorn needs to have same length as number of files. can not perform analysis")
    return(us)
  }
  if(length(unique(group)) != 2){
    warning("there need to be exactly two groups in the grouping vector. can not perform analysis")
    return(us)
  }
  if(!(lab %in% names(us$Report))){
    warning("invalid label to run two group analysis on PCA results on. needs to be contained in us$Report. can not perform analysis")
    return(us)
  }
  
  cdat <- us$Report[[lab]]
  cdat <- na_replace(cdat,fill = 0)
  
  #normalize using Zscore
  if(normalize){
    cdat <- apply(cdat,MARGIN = 2, FUN = function(x){(x - mean(x)) / sd(x)})
  }
  #run PCA
  pcares <- prcomp(cdat)
  smr <- summary(pcares)
  #determine which number of components explain the set cumulative proportion
  n_comp <- min(which(smr$importance[3,] > cumulative_proportion))
  
  #set group to 0 and 1
  grp <- as.numeric(factor(group)) -1 
  
  #build dataframe for linear model
  lmdat <- data.frame(group = grp, pcares$x[,1:n_comp])
  
  mod <- lm(lmdat,formula = group ~ .)
  mod.null <- lm(lmdat, formula = group ~ 1)
  stats <- anova(mod,mod.null)
  
  out <- list(mod = mod,mod.null = mod.null,stats = stats)
  return(out)
}


#' Plots a 2D embedding based on transition matrices of all samples contained in USdata. to plot a subset use the function SplitUSData() first. embedding can be done on multiple combined label groups
#' 
#' @param us an object of type USData
#' @param reports character or list of characters. name(s) of the Report data to be added to the 2DEmbedding data. defaults to NULL (= do not use any report data)
#' @param transitionmatrices character or vector of characters. name(s) of  transitionmatrix data to be added to the 2DEmbedding data. defaults to NULL (= do not use any report data)
#' @param transitionmatrices_stabilized character or vector of characters. name(s) of stabilized transitionmatrix data to be added to the 2DEmbedding data. defaults to NULL (= do not use any stabilized transitionmatrix data)
#' @param normalize a boolean. selects if combined 2Dembedding data should be Z-score normalized prior to embedding, defaults to TRUE
#' @param method a character. can be either "umap" (default) or "tsen". selects the method to be used for 2D embedding
#' @param colorby a character. name of the column in US$meta that should define by what the plot should be colored. if column is a character vector a discrete coloring will be applied, if a numeric value a continuous coloring will be applied. defaults to NULL = no coloring / grouping and all black
#' @param colors a character vector. contains R color names (i.e "red") or hex color code (i.e #3eb489). discrete groups will be colored following this scheme. the length needs to correspond to the number of discrete groups!. defaults to NULL = automatic colors
#' @param plot.density a boolean. selects if densities for discrete groups should be added to the plot. defaults to FALSE
#' @param plot.sem a boolean. selects if sem (standard error of the mean) for discrete groups based on X1 and X2 coordinates should be added to the plot. defaults to FALSE
#' @param seed a integer. sets the seed that will be used for 2D embedding algorithms. defaults to 123
#' @return a 2d embedding plot
#' @examples
#' Plot2DEmbedding(US,transitionmatrices = "kmeans.25")
#' Plot2DEmbedding(US,report = "raw")
#' Plot2DEmbedding(US,transitionmatrices_stabilized = "kmeans.25")
#' Plot2DEmbedding(US,transitionmatrices = "kmeans.25", normalize = FALSE, method = "tsne", seed = 404)
#' 
#' #we can use one or multiple reports data and/or one or multiple transitionmatrices data for the same embedding if we want
#' Plot2DEmbedding(US, reports = "raw", transitionmatrices = c("kmeans.25","kmeans.50))
#' 
#' #in case we have a variable "Animal" in US$meta (US$meta$Animal) we can use it to color points by animals
#' Plot2DEmbedding(US, transitionmatrices = "kmeans.25", colorby = "Animal")
#' 
#' #we can specify colors. i.e lets assume we have 4 different Animals (R color names or hex codes both work)
#' Plot2DEmbedding(US, transitionmatrices = "kmeans.25", colorby = "Animal", colors = c("red","blue","#228B22","orange))
#' 
#' #if we have a US$meta (i.e US$meta$Dosage) variable that is continuous we can create a continuous coloring from low to high
#' Plot2DEmbedding(US, transitionmatrices = "kmeans.25", colorby = "Dosage")
#' 
#' #if we have a US$meta variable that is discrete (i.e.US$meta$Group) we can create a group based coloring and also add more details such as density and the standard error of mean (sem) to the plot.
#' Plot2DEmbedding(US, transitionmatrices = "kmeans.25", colorby = "Group", plot.density = TRUE, plot.sem = TRUE)
#' 
#' #lets assume the US$meta$Group has two discrete values, "control" and "test". for the previous plot we could specify the colors by
#' Plot2DEmbedding(US, transitionmatrices = "kmeans.25", colorby = "Group", colors = c("black","#CD2990"), plot.density = TRUE, plot.sem = TRUE)
Plot2DEmbedding <- function(us, reports = NULL, transitionmatrices = NULL, transitionmatrices_stabilized = NULL, normalize = TRUE, method = "umap", colorby = NULL, colors = NULL, plot.density = FALSE, plot.sem = FALSE, seed = 123){
  require(M3C)
  require(imputeTS)
  if(is.null(reports) & is.null(transitionmatrices) & is.null(transitionmatrices_stabilized)){
    print("neiter reports nor transitionmatrices specified. using all data for 2D embedding")
    if(us$has.updated.analysis$Report){
      reports <- names(us$Report)
    }
    if(us$has.updated.analysis$Transmat){
      transitionmatrices <- us$transitionmatrix_names
    }
  }
  
  #checks if data indicated to be added to the embedding data actually exists, otherwise a warning for each dicrepancy will be produced and the function returns NULL
  check <- FALSE
  if(sum(!reports %in% names(us$Report))){
    warning("selected report(s) for embedding that is not present in US_data:\n", paste(reports[!reports %in% names(us$Report)]))
    check <- TRUE
  }
  if(sum(!transitionmatrices %in% us$transitionmatrix_names)){
    warning("selected transitionmatrix for embedding that is not present in US_data:\n", paste(transitionmatrices[!transitionmatrices %in% us$transitionmatrix_names]))
    check <- TRUE
  }
  for(i in us$file_names){
    if(sum(!transitionmatrices_stabilized %in% names(us$files[[i]]$transmatrix_stabilized))){
      warning(paste("selected stabilized transitionmatrix embedding that is not present in US_data of file:", i), "\n",
              paste(transitionmatrices_stabilized[!transitionmatrices_stabilized %in% names(us$files[[i]]$transmatrix_stabilized)], collapse = "\n"))
      check <- TRUE
    }
  }
  if(check){
    return(NULL)
  }

  # Build up data to be used for 2D Embedding by flattening all selected reports, transitionmatrices and normalized transitionmatices and building up a matrix with samples as rows and features as columns
  mat <- data.frame()
  if(!is.null(reports)){
    print(paste("including report:", reports, sep = " "))
    for(r in reports){
      mat <- rbind(mat,t(us$Report[[r]]))
    }
  }
  if(!is.null(transitionmatrices)){
    print(paste("including transitionmatrix:", transitionmatrices, sep = " "))
    tmat <- list()
    for(i in names(us$files)){
      rdat <- NULL
      for(r in transitionmatrices){
        rdat <- append(rdat,c(us$files[[i]]$transmatrix[[r]]))
      }
      tmat[[i]] <- rdat
    }
    tmat <- do.call(cbind,tmat)
    mat <- rbind(mat,tmat)
  }
  if(!is.null(transitionmatrices_stabilized)){
    print(paste("including stabilized transitionmatrix:", transitionmatrices_stabilized, sep = " "))
    tmat <- list()
    for(i in names(us$files)){
      rdat <- NULL
      for(r in transitionmatrices_stabilized){
        rdat <- append(rdat,c(us$files[[i]]$transmatrix_stabilized[[r]]))
      }
      tmat[[i]] <- rdat
    }
    tmat <- do.call(cbind,tmat)
    mat <- rbind(mat,tmat)
  }
  
  # fill NA and Normalize data
  mat <- as.matrix(imputeTS::na_replace(mat,fill = 0))
  if(normalize){
    mat <- t(apply(mat,MARGIN = 1, FUN = function(x){if(sum(x) > 0){(x - mean(x)) / (sd(x))}else{x}})) 
  }
  
  # Embedding
  if(method == "umap"){
    res <- M3C::umap(mat, seed = seed)
  }else if(method == "tsne"){
    res <- M3C::tsne(mat, seed = seed)
  }
  
  # Build up final plot
  if(!is.null(colorby)){
    res$data <- cbind(res$data,col = us$meta[,colorby])
    res$data$sz <- 1
    res$data$sz[is.na(res$data$col)] <- 0
    p1 <- ggplot(data.frame(res$data),aes(X1,X2, color = col)) + 
      geom_point(size = (1 + res$data$sz),alpha = (0.5 + 0.5*res$data$sz)) + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    if(plot.density){
      p1 <- p1 + geom_density2d()
    }
    if(plot.sem){
      aggdat <- aggregate(res$data[,1:2], FUN = mean, by = list(res$data$col))
      aggdat_sem <- aggregate(res$data[,1:2], FUN = function(x){sd(x) / sqrt(length(x))}, by = list(res$data$col))
      aggdat <- cbind(aggdat,X1.sem = aggdat_sem$X1,X2.sem = aggdat_sem$X2)
      
      semdim <- mean(c(aggdat_sem$X1, aggdat_sem$X2))
      p1 <- p1 + geom_errorbar(data = aggdat,aes(ymin = X2 - X2.sem,ymax =  X2 + X2.sem, col = Group.1), width = semdim) + 
        geom_errorbarh(data = aggdat,aes(xmin = X1 - X1.sem,xmax = X1 + X1.sem, col = Group.1), height = semdim)
    }
    if(!is.null(colors)){
      p1 <- p1 + scale_color_manual(values = colors)
    }
  }else{
    p1 <- ggplot(data.frame(res$data),aes(X1,X2)) + 
      geom_point() + theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }
  return(p1)
}


#' Plots two group transition matrix distance stats that were generated with the function TwoGroupAnalysis()
#' 
#' @param us an object of type USData
#' @param analysis character or vector of characters. name of the US$Results element(s) for which the plot should be generated. defaults to NULL = all analyses will be plotted
#' @param labels character or vector of characters. name of the label group(s) for which the plot should be generated. defaults to NULL = all labels groups will be plotted
#' @return an object of type USData
#' @examples
#' PlotTransitionsStats(US)
#' 
#' #selecting a specific analysis to plot. needs to correspond existing elements(s) with these name(s) in names(US$Results)
#' PlotTransitionsStats(US, analysis = "Ctrl-vs-Test")
#' 
#' #selecting specific label group(s) to plot
#' PlotTransitionsStats(US, labels = c("kmeans.25","kmeans.50"))
#' 
#' #selecting a specific label group in a specific analysis
#' PlotTransitionsStats(US, analysis = "Ctrl-vs-Test", labels = "kmeans.25")
PlotTransitionsStats <- function(us, analysis = NULL, labels = NULL){
  require(ggplot2)
  require(cowplot)
  if(is.null(analysis)){
    analysis <- names(us$Results)
  }
  if(is.null(labels)){
      for(i in analysis){
        labels <- unique(c(labels, names(us$Results[[i]]$TransitionStats)))
      }
  }
  
  ps <- list()
  for(i in analysis){
    for(j in labels[labels %in% names(us$Results[[i]]$TransitionStats)]){
      bdat <- data.frame(distance = us$Results[[i]]$TransitionStats[[j]]$bootstraps)
      dis <- us$Results[[i]]$TransitionStats[[j]]$distance
      sigma <- us$Results[[i]]$TransitionStats[[j]]$sigma
      percentile <- us$Results[[i]]$TransitionStats[[j]]$percentile
      p1 <- ggplot(bdat,aes(distance)) + 
          geom_histogram(color = "black", fill = "gray70", bins = 30) + 
          geom_vline(xintercept = dis, color = "red") +
          theme_bw() + 
          ggtitle(paste(i,j, sep = " ")) 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      p2 <- ggplot() + 
        ggtitle(paste("percentile:", signif(percentile,digits = 4), "\n", "sigma:", round(sigma,digits = 2), "\n", "p.value:", signif((1 - pracma::erf(sigma / sqrt(2))) / 2,digits = 3) , sep = " ")) + 
        theme_minimal() +
        theme(plot.title = element_text(size = 10))
      ps[[paste(i,j, sep = "-")]] <- cowplot::plot_grid(p1,p2, nrow = 1,rel_widths = c(1,0.3))
    }
  }
  cowplot::plot_grid(plotlist = ps, ncol = 1)
}



