#' Reads DLC Tracking data from a csv file and returns a TrackingData object
#' 
#' @param file path to a DLC tracking .csv file
#' @param fps frames per second of the recording. required to enable time resolved metrics
#' @return An object of type TrackingData
#' @examples
#' ReadDLCDataFromCSV("DLCData/Data.csv", fps = 25)
ReadDLCDataFromCSV <- function(file,fps = 1){
  out <- list()
  out$data <- list()
  data.header <- read.table(file, sep =",", header = T, nrows = 1)
  data.header <- data.frame(sapply(data.header, as.character), stringsAsFactors=FALSE)
  raw.data <- read.table(file, sep =",", header = T, skip = 2)
  for(i in seq(from = 2, to = nrow(data.header), by = 3)){
    out$data[[paste(data.header[i,])]] <- data.frame(frame = raw.data$coords, x = raw.data[,i], y = raw.data[,(i+1)], likelihood = raw.data[,(i+2)])
  }
  
  if(fps == 1){
    warning("no fps set. setting fps to 1. keep in mind that time based analyses are resolved in frames / second")
  }
  out$frames <- raw.data$coords
  out$fps <- fps
  out$seconds <- out$frames / fps
  
  out$median.data <- NULL
  for(i in names(out$data)){
    out$median.data <- rbind(out$median.data, data.frame(PointName = i, x = median(out$data[[i]]$x), y = median(out$data[[i]]$y)))
  }
  rownames(out$median.data) <- out$median.data$PointName
  
  out$point.info <- data.frame(PointName = names(out$data), PointType = "NotDefined")
  out$distance.units <- "pixel"
  out$labels <- list()
  out$filename <- last(strsplit(file,split = "/")[[1]])
  out$object.type = "TrackingData"
  
  return(out)
}


#' Reads DLC Tracking data from a multianimal pickle file
#' 
#' @param file path to a mutianimal pickle file
#' @param fps frames per second of the recording. required to enable time resolved metrics
#' @return An object of type TrackingData
#' @examples
#' ReadDLCDataFromCSV("DLCData/Data.csv", fps = 25)
ReadDLCMultianimalFromPickle <- function(file, pointnames, arenanames ,fps = 1){
  library(reticulate)
  out <- list()
  out$data <- list()
  
  pd <- import("pandas")
  pickle_data <- pd$read_pickle(file)
  
  singledat <- pickle_data$single
  pickle_data$single <- NULL
  
  flat_data <- NULL
  for(i in (1:length(pickle_data))[-1]){
    flat_data <- rbind(flat_data, c(unlist(pickle_data[[i]]), unlist(singledat[i])))
  }
  flat_data[is.nan(flat_data)] <- NA
  
  
  
   #load data
    np <- length(pointnames)
    for(j in 1:n_animals){
      for(k in 1:length(pointnames)){
        index_x <- 4 * (j - 1) * np + k
        index_y <- 4 * (j - 1) * np + k +  np
        index_likelihood <- 4 * (j - 1) * np + k +  2* np
        out$data[[paste("animal",j,pointnames[k], sep = ".")]] <- data.frame(frame = 1:nrow(flat_data), x = flat_data[,index_x], 
                                                                             y = flat_data[,index_y], 
                                                                             likelihood =  flat_data[,index_likelihood])
      }
    }
    nap <- length(arenanames)
    for(i in 1:nap){
      index_x <- 4 * n_animals * np + i
      index_y <- 4 * n_animals * np +  i + nap
      index_likelihood <- 4 * n_animals * np + i + 2 * nap
      out$data[[paste("arena",arenanames[i], sep = ".")]] <- data.frame(frame = 1:nrow(flat_data), x = flat_data[,index_x], 
                                                                           y = flat_data[,index_y], 
                                                                           likelihood =  flat_data[,index_likelihood])
    }
    
    
    
  if(fps == 1){
    warning("no fps set. setting fps to 1. keep in mind that time based analyses are resolved in frames / second")
  }
  out$frames <- 1:nrow(out$data[[1]])
  out$fps <- fps
  out$seconds <- out$frames / fps
  
  out$median.data <- NULL
  for(i in names(out$data)){
    out$median.data <- rbind(out$median.data, data.frame(PointName = i, x = median(out$data[[i]]$x,na.rm = TRUE), y = median(out$data[[i]]$y,na.rm = TRUE)))
  }
  rownames(out$median.data) <- out$median.data$PointName
  
  out$point.info <- data.frame(PointName = names(out$data), PointType = "NotDefined")
  out$distance.units <- "pixel"
  out$labels <- list()
  out$filename <- last(strsplit(file,split = "/")[[1]])
  out$object.type = "TrackingData"
  
  return(out)
}




#' Reads DLC Tracking data from .h5 file(s) and returns a TrackingData object
#' 
#' @param file path to a SLEAP tracking .h5 file
#' @param fps frames per second of the recording. required to enable time resolved metrics
#' @return An object of type TrackingData
#' @examples
#' ReadSLEAPDataFromh5("SLEAPData/Data.h5", fps = 25)
ReadSLEAPDataFromh5 <- function(files,fps = 1){
  library(rhdf5)
  library(data.table)
  out <- list()
  out$data <- list()
  for(j in files){
    h5f = H5Fopen(j)
    for(i in 1:length(h5f$node_names)){
      out$data[[paste(h5f$node_names[i])]] <- data.frame(frame = 1:length(h5f$tracking_scores),
                                                       x = h5f$tracks[,i,1,], 
                                                       y = h5f$tracks[,i,2,], 
                                                       likelihood = h5f$point_scores[,i,])
      process <- is.nan(out$data[[paste(h5f$node_names[i])]]$x)
      out$data[[paste(h5f$node_names[i])]]$x[process] <- NA
      out$data[[paste(h5f$node_names[i])]]$y[process] <- NA
    }
  }
  
  if(fps == 1){
    warning("no fps set. setting fps to 1. keep in mind that time based analyses are resolved in frames / second")
  }
  out$frames <- 1:length(h5f$tracking_scores)
  out$fps <- fps
  out$seconds <- out$frames / fps
  
  out$median.data <- NULL
  for(i in names(out$data)){
    out$median.data <- rbind(out$median.data, data.frame(PointName = i, x = median(out$data[[i]]$x, na.rm = TRUE), y = median(out$data[[i]]$y, na.rm = TRUE)))
  }
  rownames(out$median.data) <- out$median.data$PointName
  
  out$point.info <- data.frame(PointName = names(out$data), PointType = "NotDefined")
  out$distance.units <- "pixel"
  out$labels <- list()
  out$filename <- last(strsplit(files[1],split = "/")[[1]])
  out$object.type = "TrackingData"
  
  return(out)
}

#' Reads SLEAP Tracking data from a Folder and returns a list of TrackingData objects
#' 
#' @param folder path to the folder with SLEAP .h5 files
#' @param fps frames per second of the recording. required to enable time resolved metrics
#' @return An object of type TrackingData
#' @examples
#' ReadSLEAPExperiment("SLEAPData/", fps = 25)
ReadSLEAPExperiment <- function(folder, fps = 1){
  f <- list.files(folder)
  f <- f[grep(f,pattern = ".h5")]
  f <- sapply(f, FUN = function(x){unlist(strsplit(x, "[.]"))[1]})
  Tracking <- list()
  for(i in unique(f)){
    Tracking[[i]] <-  ReadSLEAPDataFromh5(paste(folder,names(f)[i == f], sep = ""),fps)
  }
  return(Tracking)
}

#' Adds a dataframe with point info to the tracking data and checks its validity
#' 
#' @param t an object of type TrackingData
#' @param pointinfo A data frame with additional point info. Requires variables: 'PointName' and `PointTyp'
#' @return An object of type TrackingData
#' @examples
#' AddPointInfo(Data, PointInfoDataFrame)
#' AddPointInfo(Data, data.frame(PointName = c("nose","tail","top","bottom), PointTyp = c("Mouse","Mouse","Maze","Maze)))
AddPointInfo <- function(t,pointinfo){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  out <- t
  if(is.null(pointinfo$PointName)){
    warning("point info does not contain required variable PointName. could not add point info")
    return(t)
  }
  if(is.null(pointinfo$PointType)){
    warning("point info does not contain required variable PointType could not add point info")
    return(t)
  }
  
  if(length(setdiff(names(t$data),pointinfo$PointName)) != 0){
    warning(paste("point info missing for following points:", setdiff(names(t$data),pointinfo$PointName), " "))
  }
  out$point.info <- pointinfo
  
  return(out)
}

#' Checks if object is of type TrackingData
#' 
#' @param t an object of type TrackingData
#' @return a boolean
#' @examples
#' IsTrackingData(Tracking)
IsTrackingData <- function(t){
  if(!is.null(t$object.type)){
    if(t$object.type == "TrackingData"){
      return(TRUE)
    }
  }
  return(FALSE)
}

#' Cuts an object of type TrackingData into a shorter object of the same type
#' 
#' @param start if not null the first start frames will be removed
#' @param end if not null the last end frames will be removed
#' @param keep.frames a numeric vector. if not null, frames that intersect with this vector will be kept
#' @param remove.frames a numeric vector. if not null, frames that intersect with this vector will be removed
#' @return An object of type TrackingData
#' @examples
#' CutTrackingData(Data, start = 100, end = 100)
#' CutTrackingData(Data, keep.frames = c(10,11,12,13))
#' CutTrackingData(Data, remove.frames = c(21,22,23,24))
CutTrackingData <- function(t,start = NULL, end = NULL, remove.frames = NULL, keep.frames = NULL){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  keep <- t$frames
  if(!is.null(start)){
    keep <- keep[-(1:start)]
  }
  if(!is.null(end)){
    keep <- keep[-((length(keep) - end):length(keep))]
  }
  if(!is.null(remove.frames)){
    keep <- setdiff(keep, remove.frames)
  }
  if(!is.null(keep.frames)){
    keep <- intersect(keep, keep.frames)
  }

  if(!is.null(t$seconds)){
    t$seconds <- t$seconds[t$frames %in% keep]
  }
  if(length(names(t$labels)) > 0){
    for(i in (names(t$labels))){
      t$labels[[i]] <- t$labels[[i]][t$frames %in% keep]
    }
  }
  if(!is.null(t$features)){
    t$features <- t$features[t$frames %in% keep,]
  }
  
  for(i in 1:length(t$data)){
    t$data[[i]] <- t$data[[i]][t$data[[i]]$frame %in% keep,]
  }
  
  t$frames <- keep
  
  return(t)
}

#' Cleanes up an object of type Tracking data. missing data is replaced by data interpolation
#' 
#' @param t an object of type TrackingData
#' @param likelihoodcutoff points below this likelihoodcutoff from DLC will be interpolated
#' @param existence.pol points outside of the polygon existence.pol will be interpolated
#' @param maxdelta points that move more than maxdelta in one frame (cm or px, depending on calibrated data or not) will be interpolated
#' @return An object of type TrackingData
#' @examples
#' CleanTrackingData(t)
#' CleanTrackingData(t, likelihoodcutoff = 0.9)
#' CleanTrackingData(t, existence.pol = data.frame(x = c(0,100,100,0), y = c(100,100,0,0)))
#' CleanTrackingData(t, likelihoodcutoff = 1, maxdelta = 5)
CleanTrackingData <- function(t, likelihoodcutoff = 0.95, existence.pol = NULL, maxdelta = NULL){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  print(paste("interpolating points with likelihood < ", likelihoodcutoff, sep = ""))
  if(!is.null(existence.pol)){
    print("interpolating points which are outside of the existence area")
  }
  if(!is.null(maxdelta)){
    print(paste("interpolating points with a maximum delta of", maxdelta, t$distance.units,"per frame", sep = " "))
  }
  
  for(i in 1:length(t$data)){
    process <- t$data[[i]]$likelihood < likelihoodcutoff
    if(!is.null(existence.pol)){
      if(is.null(existence.pol$x) | is.null(existence.pol$y)){
        warning("warning. existence.pol is invalid. polygon data needs to include variable x and variable")
      }else if(length(existence.pol$x) != length(existence.pol$y)){
        warning("invalid polygon entered. polygon data needs to include variable x and variable y of equal length")
      }else{
        process <- process | !point.in.polygon(t$data[[i]]$x,t$data[[i]]$y,existence.pol$x, existence.pol$y)
      }
    }
    if(!is.null(maxdelta)){
      process <- process | (sqrt(integratevector(t$data[[i]]$x) ^2 + integratevector(t$data[[i]]$y)^2) > maxdelta)
    }
    t$data[[i]]$x[process] <- NA
    t$data[[i]]$y[process] <- NA
    t$data[[i]] <- na_interpolation(t$data[[i]])
  }
  return(t)
}

#' Calibrates an object of type TrackingData from pixel to metric space
#' 
#' @param t an object of type TrackingData to be calibrated
#' @param method use "distance" or "area". distance requires 2 points, area a polygon with > 2 points
#' @param in.metric the measured distance or area in metric units
#' @param points a vector of tracked points that are used for calibration
#' @return An object of type TrackingData
#' @examples
#' CalibrateTrackingData(t, method = "distance", in.metric = 40, points = c("top","bottom"))
#' CalibrateTrackingData(t, method = "area", in.metric = 1600, points = c("top.left","top.right","bottom.right","bottom.left"))
CalibrateTrackingData <- function(t, method, in.metric = NULL, points = NULL, ratio = NULL, new.units = NULL){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(!method %in% c("distance","area","ratio")){
    warning("invalid method: valid methods are distance, area or ratio, can not calibrate")
    return(t)
  }
  if(method == "ratio"){
    if(is.numeric(ratio)){
      t$px.to.cm <- ratio
    }else{
      warning("method ratio needs a valid ratio to be entered. can not calibrate")
      return(t)
    }
  }else{
    if(is.numeric(in.metric) & is.null(points)){
      warning("method requires both points and a in.metric (numeric!) measurement")
      return(t)
    }
    if(sum(!(points %in% t$median.data$PointName))){
      warning("invalid points entered, can not calibrate")
      return(t)
    }
    
    if(method == "distance"){
      if(length(points) == 2){
        t$px.to.cm <- in.metric / Distance2d(t$median.data[points[1],],t$median.data[points[2],]) 
      }else{
        warning("invalid number of points: distance needs 2 points, can not calibrate")
        return(t)
      }
    }
    if(method == "area"){
      if(length(points) > 2){
        t$px.to.cm <- sqrt(in.metric / AreaPolygon2d(t$median.data[points,])) 
      }else{
        warning("invalid number of points: area need polygon of > 2 points, can not calibrate")
        return(t)
      }
    }
  }
  
  for(i in 1:length(t$data)){
    t$data[[i]]$x <- t$data[[i]]$x * t$px.to.cm
    t$data[[i]]$y <- t$data[[i]]$y * t$px.to.cm
  }
  if(!is.null(t$median.data)){
    t$median.data[,c("x","y")] <- t$median.data[,c("x","y")] * t$px.to.cm
  }
  if(!is.null(t$zones)){
    for(i in 1:length(t$zones)){
      t$zones[[i]]$x <- t$zones[[i]]$x * t$px.to.cm
      t$zones[[i]]$y <- t$zones[[i]]$y * t$px.to.cm
    }
  }
  if(!is.null(new.units)){
    t$distance.units <- new.units
  }else{
    t$distance.units <- "cm"
  }
  
  return(t)
} 

#' Scales a polygon by a factor
#' 
#' @param p a polygon (type list() or data.frame() with elements x and y are required)
#' @param factor scaling facor
#' @return a polygon
#' @examples
#' ScalePolygon(p = data.frame(x = c(0,0,1,1), y = c(0,1,1,0)), factor = 1.3)
ScalePolygon <-function(p, factor){
  if(is.null(p$x) | is.null(p$y)){
    stop("invalid input. Needs polygon of type list() or data.frame() with 2 variables, x = and y =")
  }
  center_x <- mean(p$x)
  center_y <- mean(p$y)
  return(data.frame(x = center_x + factor * (p$x - center_x), y = center_y + factor * (p$y - center_y)))
}

#' Recenters a polygon to a new center point
#' 
#' @param p a polygon (type list() or data.frame() with elements x and y are required)
#' @param new_center the new center point (type list() or data.frame() with elements x and y are required)
#' @return a polygon
#' @examples
#' RecenterPolygon(p = data.frame(x = c(0,0,1,1), y = c(0,1,1,0)), new_center = data.frame(x=1,y=1))
RecenterPolygon <-function(p, new_center){
  center_x <- mean(p$x)
  center_y <- mean(p$y)
  return(data.frame(x = p$x + new_center$x - center_x , y = p$y + new_center$y - center_y))
}

#' Adds Zones required for OFT analysis to an object of type TrackingData
#' 
#' @param t an object of type TrackingData
#' @param points a vector of point names describing the corners of the arena (!THE ORDER IS IMPORTANT)
#' @param scale_center a factor by which the arena is scaled to define the center
#' @param scale_periphery a factor by which the arena is scaled to define the periphery
#' @param scale_corners a factor by which the arena is scaled before recentering to each corner
#' @return a polygon
#' @examples
#' AddOFTZones(t, c("tl","tr","br","bl"), 0.5,0.4,0.8)
AddOFTZones <- function(t, points = c("tl","tr","br","bl"), scale_center = 0.5, scale_corners = 0.4, scale_periphery = 0.8){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(length(intersect(points,names(t$data) != 4))){
    warning("invalid number or type of points entered. exactly 4 existing points needed for OFT Zones")
    return(t)
  }
  zones <- list()
  zones.invert <- list()
  zones[["arena"]] <- t$median.data[points,c("x","y")]
  zones.invert[["arena"]] <- FALSE
  zones[["center"]] <- ScalePolygon(t$median.data[points,c("x","y")], scale_center)
  zones.invert[["center"]] <- FALSE
  zones[["periphery"]] <- ScalePolygon(t$median.data[points,c("x","y")], scale_periphery)
  zones.invert[["periphery"]] <- TRUE
  t$corner.names <- NULL
  for(i in points){
    zones[[paste("corner",i, sep = ".")]] <- RecenterPolygon(ScalePolygon(t$median.data[points,c("x","y")], scale_corners), t$median.data[i,c("x","y")])
    zones.invert[[paste("corner",i, sep = ".")]] <- FALSE
    t$corner.names <- append(t$corner.names,paste("corner",i, sep = "."))
  }
  t$zones <- zones
  t$zones.invert <- zones.invert
  return(t)
}

#' Performs an OFT analysis on an object of type TrackingData
#' 
#' @param t an object of type TrackingData
#' @param movement_cutoff a numeric value that denotes the cutoff point above which an animal is considered moving
#' @param integration_period a numeric value that denotes the duration over which metrics are smoothed.
#' @param points a string or vector of strings that denotes the name of points which will be analysed
#' @return a TrackingData object
#' @examples
#' OFTAnalysis(t, movement_cutoff = 5,integration_period = 5,points = "bodycentre")
OFTAnalysis <- function(t, movement_cutoff,integration_period, points){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(t$zones)){
    warning("no zones defined for OFT analysis. Returing simple analysis only")
  }
  
  t <- CalculateMovement(t,movement_cutoff,integration_period)
  
  t$Report <- list()
  if(!is.null(t$labels)){
    t$Report <- append(t$Report,LabelReport(t, integration_period))
  }
  
  for(k in points){
    dat <- t$data[[k]]
    t$Report[[paste(k, "raw.distance", sep = ".")]] <- sum(dat[,"speed"], na.rm = T)
    t$Report[[paste(k, "distance.moving", sep = ".")]] <- sum(dat[dat$is.moving,"speed"], na.rm = T)
    t$Report[[paste(k, "raw.speed", sep = ".")]] <- mean(dat[,"speed"], na.rm = T) * t$fps
    t$Report[[paste(k, "speed.moving", sep = ".")]] <- mean(dat[dat$is.moving,"speed"], na.rm = T) * t$fps
    t$Report[[paste(k, "time.moving", sep = ".")]] <- sum(dat[,"is.moving"], na.rm = T) / t$fps
    t$Report[[paste(k, "total.time", sep = ".")]] <- length(dat[,"is.moving"]) / t$fps
    t$Report[[paste(k, "time.stationary", sep = ".")]] <- t$Report[[paste(k, "total.time", sep = ".")]] - t$Report[[paste(k, "time.moving", sep = ".")]]
    t$Report[[paste(k, "percentage.moving", sep = ".")]] <- t$Report[[paste(k, "time.moving", sep = ".")]] / t$Report[[paste(k, "total.time", sep = ".")]] * 100
    
    if(!is.null(t$zones)){
      t$Report <- append(t$Report, ZoneReport(t,k,"center", zone.name = paste(k,"center", sep = ".")))
      t$Report <- append(t$Report, ZoneReport(t,k,"periphery" , zone.name = paste(k,"periphery", sep = "."), invert = TRUE))
      t$Report <- append(t$Report, ZoneReport(t,k,t$corner.names, zone.name = paste(k,"corners", sep = ".")))
    }
  }
  return(t)
}

#' Performs an EPM analysis on an object of type TrackingData
#' 
#' @param t an object of type TrackingData
#' @param movement_cutoff a numeric value that denotes the cutoff point above which an animal is considered moving
#' @param integration_period a numeric value that denotes the duration over which metrics are smoothed.
#' @param points a string or vector of strings that denotes the name of points which will be analysed
#' @return a TrackingData object
#' @examples
#' EPMAnalysis(t, 5,5,"bodycentre")
EPMAnalysis <- function(t, movement_cutoff,integration_period, points,nosedips = FALSE){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(t$zones)){
    warning("no zones defined for EPM analysis. Returing simple analysis only")
  }
  
  t <- CalculateMovement(t, movement_cutoff,integration_period)
  t$Report <- list()
  
  #THIS IS A VERA ARBITRARY SECTION FOR THIS TYPE OF ANALYSIS. IT WILL ONLY WORK WITH THE CORRECT POINTS AND ZONES
  if(nosedips){
    if((length(setdiff(c("headcentre","bodycentre","neck"),names(t$data))) != 0) | (length(setdiff(c("closed.top","closed.bottom","arena"),names(t$zones))) != 0)){
      warning("Not all points or zones needed for nosedip analysis. Requires points : headcentre,bodycentre,neck and zones closed.top, closed.bottom, arena")
    }else{
      t$labels$automatic.nosedip <- avgbool(!IsInZone(t,"headcentre","arena") & IsInZone(t,"bodycentre","arena") &!IsInZone(t,"neck",c("closed.top","closed.bottom")),integration_period)
      t$Report[["nose.dip"]] <- CalculateTransitions(t$labels$automatic.nosedip,integration_period) / 2
      t$labels$automatic.nosedip <- ifelse(t$labels$automatic.nosedip == 1,"Nosedip","None")
    }
  }
  
  for(k in points){
    dat <- t$data[[k]]
    t$Report[[paste(k, "raw.distance", sep = ".")]] <- sum(dat[,"speed"], na.rm = T)
    t$Report[[paste(k, "distance.moving", sep = ".")]] <- sum(dat[dat$is.moving,"speed"], na.rm = T)
    t$Report[[paste(k, "raw.speed", sep = ".")]] <- mean(dat[,"speed"], na.rm = T) * t$fps
    t$Report[[paste(k, "speed.moving", sep = ".")]] <- mean(dat[dat$is.moving,"speed"], na.rm = T) * t$fps
    t$Report[[paste(k, "time.moving", sep = ".")]] <- sum(dat[,"is.moving"], na.rm = T) / t$fps
    t$Report[[paste(k, "total.time", sep = ".")]] <- length(dat[,"is.moving"]) / t$fps
    t$Report[[paste(k, "time.stationary", sep = ".")]] <- t$Report[[paste(k, "total.time", sep = ".")]] - t$Report[[paste(k, "time.moving", sep = ".")]]
    t$Report[[paste(k, "percentage.moving", sep = ".")]] <- t$Report[[paste(k, "time.moving", sep = ".")]] / t$Report[[paste(k, "total.time", sep = ".")]] * 100
    
    t$Report <- append(t$Report, ZoneReport(t,k,"center", zone.name = paste(k,"center", sep = ".")))
    t$Report <- append(t$Report, ZoneReport(t,k,c("open.right","open.left"), zone.name = paste(k,"open", sep = ".")))
    t$Report <- append(t$Report, ZoneReport(t,k,c("closed.top","closed.bottom"), zone.name = paste(k,"closed", sep = ".")))
    t$Report <- append(t$Report, ZoneReport(t,k,"closed.top",  zone.name = paste(k,"closed.top", sep = ".")))
    t$Report <- append(t$Report, ZoneReport(t,k,"closed.bottom",  zone.name = paste(k,"closed.bottom", sep = ".")))
    t$Report <- append(t$Report, ZoneReport(t,k,"open.left", zone.name = paste(k,"open.left", sep = ".")))
    t$Report <- append(t$Report, ZoneReport(t,k,"open.right", zone.name = paste(k,"open.right", sep = ".")))
  }
  return(t)
}

#' Performs an EPM analysis on an object of type TrackingData
#' 
#' @param t an object of type TrackingData
#' @param movement_cutoff a numeric value that denotes the cutoff point above which an animal is considered moving
#' @param integration_period a numeric value that denotes the duration over which metrics are smoothed.
#' @param points a string or vector of strings that denotes the name of points which will be analysed
#' @return a TrackingData object
#' @examples
#' EPMAnalysis(t, 5,5,"bodycentre")
LDBAnalysis <- function(t, movement_cutoff,integration_period, points){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(t$zones)){
    warning("no zones defined for LDB analysis. Returing simple analysis only")
  }
  
  t <- CalculateMovement(t, movement_cutoff,integration_period)
  t$Report <- list()
  if(!is.null(t$labels)){
    t$Report <- append(t$Report,LabelReport(t, integration_period))
  }
  
  for(k in points){
    dat <- t$data[[k]]
    t$Report[[paste(k, "raw.distance", sep = ".")]] <- sum(dat[,"speed"], na.rm = T)
    t$Report[[paste(k, "distance.moving", sep = ".")]] <- sum(dat[dat$is.moving,"speed"], na.rm = T)
    t$Report[[paste(k, "raw.speed", sep = ".")]] <- mean(dat[,"speed"], na.rm = T) * t$fps
    t$Report[[paste(k, "speed.moving", sep = ".")]] <- mean(dat[dat$is.moving,"speed"], na.rm = T) * t$fps
    t$Report[[paste(k, "time.moving", sep = ".")]] <- sum(dat[,"is.moving"], na.rm = T) / t$fps
    t$Report[[paste(k, "total.time", sep = ".")]] <- length(dat[,"is.moving"]) / t$fps
    t$Report[[paste(k, "time.stationary", sep = ".")]] <- t$Report[[paste(k, "total.time", sep = ".")]] - t$Report[[paste(k, "time.moving", sep = ".")]]
    t$Report[[paste(k, "percentage.moving", sep = ".")]] <- t$Report[[paste(k, "time.moving", sep = ".")]] / t$Report[[paste(k, "total.time", sep = ".")]] * 100
    
    t$Report <- append(t$Report, ZoneReport(t,k,"light", zone.name = paste(k,"light", sep = ".")))
    t$Report <- append(t$Report, ZoneReport(t,k,"dark", zone.name = paste(k,"dark", sep = ".")))
    t$Report[["latency.to.entry.dark"]] <- LatencyToEntry(t,k,"dark")
  }
  return(t)
}

#' Checks if point p is in zone(s) z
#' 
#' @param t an object of type TrackingData
#' @param p string value of a point
#' @param z a string or vector of strings naming the zones to be checked
#' @param invert if TRUE instead it will be checked if the point is outside the zone
#' @return a boolean vector for the assessment at each frame
#' @examples
#' IsInZone(t, "bodycentre","center")
IsInZone <- function(t,p,z,invert = FALSE){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(length(intersect(p,names(t$data))) != 1){
    warning(paste("Points not available in Tracking data:",paste(setdiff(points, names(t$data)),collapse = " "), sep = " "))
    return(NULL)
  }
  if(length(intersect(z,names(t$zones))) != length(z)){
    warning(paste("Zones in Tracking data:",paste(setdiff(z, names(t$zones)),collapse = " "), sep = " "))
    return(NULL)
  }
  
  
  zones <- t$zones[z]
  in.zone <- rep(FALSE,nrow(t$data[[p]]))
  for(i in zones){
    in.zone <- in.zone | (point.in.polygon(t$data[[p]]$x,t$data[[p]]$y,i$x,i$y) == 1)
  }
  if(invert){
    in.zone <- !in.zone
  }
  return(in.zone)
}

LatencyToEntry <- function(t,p,z,invert = FALSE){
  in.zone <- IsInZone(t,p,z,invert)
  return(min(which(in.zone == TRUE)) / t$fps)
}

#' Creates a report for each behavior which includes total time of behavior and number of onset/offsets
#' 
#' @param t an object of type TrackingData
#' @param integration_period string value of a point
#' @return a list of metrics for each behavior
#' @examples
#' ClassificationReport(t, integration_period = 5)
LabelReport <- function(t, integration_period = 0){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(length(t$labels) == 0){
    warning("no label data present")
    return(NULL)
  }
  Report <- list()
  for(j in names(t$labels)){
    c <- SmoothLabel(t$labels[[j]],integration_period)
    c <- na.omit(c)
    for(i in unique(c)){
      Report[[paste(j,i,"time", sep = ".")]] <- sum(c == i) / t$fps
      Report[[paste(j,i,"count", sep = ".")]] <- sum(CalculateTransitions(c == i, 0)) / 2
    }
  }
  return(Report)
}

CalculateMovement <- function(t, movement_cutoff, integration_period){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  for(i in 1:length(t$data)){
    t$data[[i]]$delta_x <- integratevector(t$data[[i]]$x)
    t$data[[i]]$delta_y <- integratevector(t$data[[i]]$y)
    t$data[[i]]$speed <- sqrt(t$data[[i]]$delta_x ^2 + t$data[[i]]$delta_y^2)
    t$data[[i]]$acceleration <- integratevector(t$data[[i]]$speed)
    t$data[[i]]$is.moving <- as.logical(avgbool(t$data[[i]]$speed > (movement_cutoff / t$fps), integration_period))
  }
  t$integration_period <- integration_period
  return(t)
}

CalculateAccelerations <- function(t){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  for(i in 1:length(t$data)){
    t$data[[i]]$delta_x <- integratevector(t$data[[i]]$x)
    t$data[[i]]$delta_y <- integratevector(t$data[[i]]$y)
    t$data[[i]]$speed <- sqrt(t$data[[i]]$delta_x ^2 + t$data[[i]]$delta_y^2)
    t$data[[i]]$acceleration <- integratevector(t$data[[i]]$speed)
  }
  return(t)
}

ZoneReport <- function(t,point,zones, zone.name = NULL, invert = FALSE){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  Report <- list()
  if(is.null(zone.name)){
    zone.name <- paste(zones, collapse = ".")
  }
  if(!point %in% names(t$data)){
    warning("Invalid point")
    return(NULL)
  }
  if(!sum(zones %in% names(t$zones))){
    warning("Invalid zone(s)")
    return(NULL)
  }
  
  dat <- t$data[[point]]
  in.zone <- rep(FALSE,nrow(dat))
  for(i in t$zones[zones]){
    in.zone <- in.zone | (point.in.polygon(dat$x,dat$y,i$x,i$y) == 1)
  }
  if(invert){
    in.zone <- !in.zone
  }
  Report[[paste(zone.name, "raw.distance", sep = ".")]] <- sum(dat[in.zone,"speed"], na.rm = T)
  Report[[paste(zone.name, "distance.moving", sep = ".")]] <- sum(dat[dat$is.moving & in.zone,"speed"], na.rm = T)
  Report[[paste(zone.name, "raw.speed", sep = ".")]] <- mean(dat[in.zone,"speed"], na.rm = T) * t$fps
  Report[[paste(zone.name, "speed.moving", sep = ".")]] <- mean(dat[dat$is.moving & in.zone,"speed"], na.rm = T) * t$fps
  Report[[paste(zone.name, "time.moving", sep = ".")]] <- sum(dat[in.zone,"is.moving"], na.rm = T) / t$fps
  Report[[paste(zone.name, "total.time", sep = ".")]] <- length(dat[in.zone,"is.moving"]) / t$fps
  Report[[paste(zone.name, "time.stationary", sep = ".")]] <- Report[[paste(zone.name, "total.time", sep = ".")]] - Report[[paste(zone.name, "time.moving", sep = ".")]]
  Report[[paste(zone.name, "percentage.moving", sep = ".")]] <- Report[[paste(zone.name, "time.moving", sep = ".")]] / Report[[paste(zone.name, "total.time", sep = ".")]] * 100
  Report[[paste(zone.name, "transitions", sep = ".")]] <- CalculateTransitions(in.zone, t$integration_period) 
  
  return(Report)
}

Distance2d <- function(a,b){
  sqrt((a$x - b$x)^2 + (a$y - b$y)^2)
}

MedianMouseArea <- function(t,points = c("nose","earr","bcr","hipr","tailbase","hipl","bcl","earl")){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  median(GetPolygonAreas(t,points), na.rm = T)
}

MedianMouseLength <- function(t, front  = "nose", back = "tailbase"){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  median(GetDistances(t,front,back), na.rm = T)
}

AreaPolygon2d <- function(p){
  area <- 0
  for(i in 1:nrow(p)){
    if(i != nrow(p)){
      area <- area + p$x[i]*p$y[i+1] - p$y[i]*p$x[i+1]
    }else{
      area <- area + p$x[i]*p$y[1] - p$y[i]*p$x[1]
    }
  }
  return(abs(area/2))
}

GetPolygonAreas <-function(t,points){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(length(intersect(points,names(t$data))) != length(unique(points))){
    warning(paste("Points for distance measurement not available in Tracking data:",paste(setdiff(points, names(t$data)),collapse = " "), sep = " "))
    return(NULL)
  }
  x <- NULL
  y <- NULL
  for(i in points){
    x <- cbind(x,t$data[[i]]$x )
    y <- cbind(y,t$data[[i]]$y )
  }
  AreaPolygon3d(x,y)
}

GetDistances <- function(t,f,b){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(length(intersect(c(f,b),names(t$data))) != 2){
    warning(paste("Points for distance measurement not available in Tracking data:",paste(setdiff(c(f,b), names(t$data)),collapse = " "), sep = " "))
  }
  Distance2d(t$data[[f]],t$data[[b]])
}

AreaPolygon3d <- function(x,y){
  area <- 0
  for(i in 1:dim(x)[2]){
    if(i != dim(x)[2]){
      area <- area + x[,i]*y[,i+1] - y[,i]*x[,i+1]
    }else{
      area <- area + x[,i]*y[,1] - y[,i]*x[,1]
    }
  }
  return(abs(area/2))
}

AddLabelingData <- function(t, lab){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  t$labels$manual <- rep("None",length(t$seconds))
  for(i in 1:nrow(lab)){
    t$labels$manual[t$seconds >= lab[i,"from"] & t$seconds <= lab[i,"to"]] <- as.character(lab[i,"type"])
  }
  return(t)
}

CreateSkeletonData <- function(t){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  dat <- data.frame(S1 = GetDistances(t,"nose","headcentre"))
  dat$S2 <- GetDistances(t,"headcentre","neck")
  dat$S3 <- GetDistances(t,"neck","bodycentre")
  dat$S4 <- GetDistances(t,"bodycentre","bcr")
  dat$S5 <- GetDistances(t,"bodycentre","bcl")
  dat$S6 <- GetDistances(t,"bodycentre","tailbase")
  dat$S7 <- GetDistances(t,"tailbase","hipr")
  dat$S8 <- GetDistances(t,"tailbase","hipl")
  dat$S9 <- GetDistances(t,"tailbase","tailcentre")
  dat$S10 <- GetDistances(t,"tailcentre","tailtip")
  dat$A1 <- GetAngleTotal(t,"tailbase","tailcentre","tailcentre","tailtip")
  dat$A2 <- GetAngleTotal(t,"hipr","tailbase","tailbase","hipl")
  dat$A3 <- GetAngleTotal(t,"tailbase","bodycentre","bodycentre","neck")
  dat$A4 <- GetAngleTotal(t,"bcr","bodycentre","bodycentre","bcl")
  dat$A5 <- GetAngleTotal(t,"bodycentre","neck","neck","headcentre")
  dat$A6 <- GetAngleTotal(t,"tailbase","bodycentre","neck","headcentre")
  dat$Ar1 <- GetPolygonAreas(t,c("tailbase","hipr","hipl"))
  dat$Ar2 <- GetPolygonAreas(t,c("hipr","hipl","bcl","bcr"))
  dat$Ar3 <- GetPolygonAreas(t,c("bcr","earr","earl","bcl"))
  dat$Ar4 <- GetPolygonAreas(t,c("earr","nose","earl"))
  dat <- as.data.frame(na_replace(dat)) 
  t$features <- dat
  return(t)
}

CreateSkeletonData_OFT <- function(t){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  dat <- data.frame(S1 = GetDistances(t,"nose","headcentre"))
  dat$S2 <- GetDistances(t,"headcentre","neck")
  dat$S3 <- GetDistances(t,"neck","bodycentre")
  dat$S4 <- GetDistances(t,"bodycentre","bcr")
  dat$S5 <- GetDistances(t,"bodycentre","bcl")
  dat$S6 <- GetDistances(t,"bodycentre","tailbase")
  dat$S7 <- GetDistances(t,"tailbase","hipr")
  dat$S8 <- GetDistances(t,"tailbase","hipl")
  dat$S9 <- GetDistances(t,"tailbase","tailcentre")
  dat$S10 <- GetDistances(t,"tailcentre","tailtip")
  dat$A1 <- GetAngleTotal(t,"tailbase","tailcentre","tailcentre","tailtip")
  dat$A2 <- GetAngleTotal(t,"hipr","tailbase","tailbase","hipl")
  dat$A3 <- GetAngleTotal(t,"tailbase","bodycentre","bodycentre","neck")
  dat$A4 <- GetAngleTotal(t,"bcr","bodycentre","bodycentre","bcl")
  dat$A5 <- GetAngleTotal(t,"bodycentre","neck","neck","headcentre")
  dat$A6 <- GetAngleTotal(t,"tailbase","bodycentre","neck","headcentre")
  dat$Ar1 <- GetPolygonAreas(t,c("tailbase","hipr","hipl"))
  dat$Ar2 <- GetPolygonAreas(t,c("hipr","hipl","bcl","bcr"))
  dat$Ar3 <- GetPolygonAreas(t,c("bcr","earr","earl","bcl"))
  dat$Ar4 <- GetPolygonAreas(t,c("earr","nose","earl"))
  dat$P1 <- as.integer(IsInZone(t,"nose","arena"))
  dat$P2 <- as.integer(IsInZone(t,"headcentre","arena"))
  dat <- as.data.frame(na_replace(dat)) 
  t$features <- dat
  
  return(t)
}

CreateSkeletonData_OFT_v2 <- function(t){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  dat <- data.frame(Ac1 = t$data[["nose"]]$acceleration)
  dat$Ac2 <- t$data[["headcentre"]]$acceleration
  dat$Ac3 <- t$data[["neck"]]$acceleration
  dat$Ac4 <- t$data[["earr"]]$acceleration
  dat$Ac5 <- t$data[["earl"]]$acceleration
  dat$Ac6 <- t$data[["bodycentre"]]$acceleration
  dat$Ac7 <- t$data[["bcl"]]$acceleration
  dat$Ac8 <- t$data[["bcr"]]$acceleration
  dat$Ac9 <- t$data[["hipl"]]$acceleration
  dat$Ac10 <- t$data[["hipr"]]$acceleration
  dat$Ac11 <- t$data[["tailbase"]]$acceleration
  dat$S1 <- GetDistances(t,"nose","headcentre")
  dat$S2 <- GetDistances(t,"headcentre","neck")
  dat$S3 <- GetDistances(t,"neck","bodycentre")
  dat$S4 <- GetDistances(t,"bodycentre","bcr")
  dat$S5 <- GetDistances(t,"bodycentre","bcl")
  dat$S6 <- GetDistances(t,"bodycentre","tailbase")
  dat$S7 <- GetDistances(t,"tailbase","hipr")
  dat$S8 <- GetDistances(t,"tailbase","hipl")
  dat$S9 <- GetDistances(t,"tailbase","tailcentre")
  dat$S10 <- GetDistances(t,"tailcentre","tailtip")
  dat$A1 <- GetAngleTotal(t,"tailbase","tailcentre","tailcentre","tailtip")
  dat$A2 <- GetAngleTotal(t,"hipr","tailbase","tailbase","hipl")
  dat$A3 <- GetAngleTotal(t,"tailbase","bodycentre","bodycentre","neck")
  dat$A4 <- GetAngleTotal(t,"bcr","bodycentre","bodycentre","bcl")
  dat$A5 <- GetAngleTotal(t,"bodycentre","neck","neck","headcentre")
  dat$A6 <- GetAngleTotal(t,"tailbase","bodycentre","neck","headcentre")
  dat$Ar1 <- GetPolygonAreas(t,c("tailbase","hipr","hipl"))
  dat$Ar2 <- GetPolygonAreas(t,c("hipr","hipl","bcl","bcr"))
  dat$Ar3 <- GetPolygonAreas(t,c("bcr","earr","earl","bcl"))
  dat$Ar4 <- GetPolygonAreas(t,c("earr","nose","earl"))
  dat$P1 <- as.integer(IsInZone(t,"nose","arena"))
  dat$P2 <- as.integer(IsInZone(t,"headcentre","arena"))
  dat <- as.data.frame(na_replace(dat)) 
  t$features <- dat
  
  return(t)
}

CreateSkeletonData_OFT_v3 <- function(t){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  dat <- data.frame(Ac1 = t$data[["nose"]]$acceleration)
  dat$Ac2 <- t$data[["headcentre"]]$acceleration
  dat$Ac3 <- t$data[["neck"]]$acceleration
  dat$Ac4 <- t$data[["earr"]]$acceleration
  dat$Ac5 <- t$data[["earl"]]$acceleration
  dat$Ac6 <- t$data[["bodycentre"]]$acceleration
  dat$Ac7 <- t$data[["bcl"]]$acceleration
  dat$Ac8 <- t$data[["bcr"]]$acceleration
  dat$Ac9 <- t$data[["hipl"]]$acceleration
  dat$Ac10 <- t$data[["hipr"]]$acceleration
  dat$Ac11 <- t$data[["tailbase"]]$acceleration
  dat$S1 <- GetDistances(t,"nose","headcentre")
  dat$S2 <- GetDistances(t,"headcentre","neck")
  dat$S3 <- GetDistances(t,"neck","bodycentre")
  dat$S4 <- GetDistances(t,"bodycentre","bcr")
  dat$S5 <- GetDistances(t,"bodycentre","bcl")
  dat$S6 <- GetDistances(t,"bodycentre","tailbase")
  dat$S7 <- GetDistances(t,"tailbase","hipr")
  dat$S8 <- GetDistances(t,"tailbase","hipl")
  dat$S9 <- GetDistances(t,"tailbase","tailcentre")
  dat$S10 <- GetDistances(t,"tailcentre","tailtip")
  dat$A1 <- GetAngleTotal(t,"tailbase","tailcentre","tailcentre","tailtip")
  dat$A2 <- GetAngleTotal(t,"hipr","tailbase","tailbase","hipl")
  dat$A3 <- GetAngleTotal(t,"tailbase","bodycentre","bodycentre","neck")
  dat$A4 <- GetAngleTotal(t,"bcr","bodycentre","bodycentre","bcl")
  dat$A5 <- GetAngleTotal(t,"bodycentre","neck","neck","headcentre")
  dat$A6 <- GetAngleTotal(t,"tailbase","bodycentre","neck","headcentre")
  dat$Ar1 <- GetPolygonAreas(t,c("tailbase","hipr","hipl"))
  dat$Ar2 <- GetPolygonAreas(t,c("hipr","hipl","bcl","bcr"))
  dat$Ar3 <- GetPolygonAreas(t,c("bcr","earr","earl","bcl"))
  dat$Ar4 <- GetPolygonAreas(t,c("earr","nose","earl"))
  dat$D1 <- GetDistanceToZoneBorder(t,"arena","nose")
  dat$D2 <- GetDistanceToZoneBorder(t,"arena","neck")
  dat$D3 <- GetDistanceToZoneBorder(t,"arena","bodycentre")
  dat$D4 <- GetDistanceToZoneBorder(t,"arena","tailbase")
  dat <- as.data.frame(na_replace(dat)) 
  t$features <- dat
  
  return(t)
}

CreateSkeletonData_OFT_v4 <- function(t){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  dat <- data.frame(Ac1 = t$data[["nose"]]$acceleration)
  dat$Ac2 <- t$data[["headcentre"]]$acceleration
  dat$Ac3 <- t$data[["neck"]]$acceleration
  dat$Ac4 <- t$data[["earr"]]$acceleration
  dat$Ac5 <- t$data[["earl"]]$acceleration
  dat$Ac6 <- t$data[["bodycentre"]]$acceleration
  dat$Ac7 <- t$data[["bcl"]]$acceleration
  dat$Ac8 <- t$data[["bcr"]]$acceleration
  dat$Ac9 <- t$data[["hipl"]]$acceleration
  dat$Ac10 <- t$data[["hipr"]]$acceleration
  dat$Ac11 <- t$data[["tailbase"]]$acceleration
  dat$A1 <- GetAngleTotal(t,"hipr","tailbase","tailbase","hipl")
  dat$A2 <- GetAngleTotal(t,"tailbase","bodycentre","bodycentre","neck")
  dat$A3 <- GetAngleTotal(t,"bcr","bodycentre","bodycentre","bcl")
  dat$A4 <- GetAngleTotal(t,"bodycentre","neck","neck","headcentre")
  dat$A5 <- GetAngleTotal(t,"tailbase","bodycentre","neck","headcentre")
  dat$A6 <- GetAngleTotal(t,"bcl","hipl","bcl","earl")
  dat$A7 <- GetAngleTotal(t,"bcr","hipr","bcr","earr")
  dat$A8 <- GetAngleTotal(t,"nose","earr","nose","earl")
  dat$D1 <- GetDistanceToZoneBorder(t,"arena","nose")
  dat$D2 <- GetDistanceToZoneBorder(t,"arena","neck")
  dat$D3 <- GetDistanceToZoneBorder(t,"arena","bodycentre")
  dat$D4 <- GetDistanceToZoneBorder(t,"arena","tailbase")
  dat$S1 <- GetDistances(t,"nose","headcentre")
  dat$S2 <- GetDistances(t,"headcentre","neck")
  dat$S3 <- GetDistances(t,"neck","bodycentre")
  dat$S4 <- GetDistances(t,"bodycentre","bcr")
  dat$S5 <- GetDistances(t,"bodycentre","bcl")
  dat$S6 <- GetDistances(t,"bodycentre","tailbase")
  dat$S7 <- GetDistances(t,"tailbase","hipr")
  dat$S8 <- GetDistances(t,"tailbase","hipl")
  dat$S9 <- GetDistances(t,"bcr","hipr")
  dat$S10 <- GetDistances(t,"bcl","hipl")
  dat$S11 <- GetDistances(t,"bcl","earl")
  dat$S12 <- GetDistances(t,"bcr","earr")
  dat$S13 <- GetDistances(t,"nose","earr")
  dat$S14 <- GetDistances(t,"nose","earl")
  dat$Ar1 <- GetPolygonAreas(t,c("tailbase","hipr","hipl"))
  dat$Ar2 <- GetPolygonAreas(t,c("hipr","hipl","bcl","bcr"))
  dat$Ar3 <- GetPolygonAreas(t,c("bcr","earr","earl","bcl"))
  dat$Ar4 <- GetPolygonAreas(t,c("earr","nose","earl"))
  dat <- as.data.frame(na_replace(dat)) 
  t$features <- dat
  
  return(t)
}

CreateSkeletonData_LDB<- function(t){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  dat <- data.frame(Ac1 = t$data[["nose"]]$acceleration)
  dat$Ac2 <- t$data[["headcentre"]]$acceleration
  dat$Ac3 <- t$data[["neck"]]$acceleration
  dat$Ac4 <- t$data[["earr"]]$acceleration
  dat$Ac5 <- t$data[["earl"]]$acceleration
  dat$Ac6 <- t$data[["bodycentre"]]$acceleration
  dat$Ac7 <- t$data[["bcl"]]$acceleration
  dat$Ac8 <- t$data[["bcr"]]$acceleration
  dat$Ac9 <- t$data[["hipl"]]$acceleration
  dat$Ac10 <- t$data[["hipr"]]$acceleration
  dat$Ac11 <- t$data[["tailbase"]]$acceleration
  dat$S1 <- GetDistances(t,"nose","headcentre")
  dat$S2 <- GetDistances(t,"headcentre","neck")
  dat$S3 <- GetDistances(t,"neck","bodycentre")
  dat$S4 <- GetDistances(t,"bodycentre","bcr")
  dat$S5 <- GetDistances(t,"bodycentre","bcl")
  dat$S6 <- GetDistances(t,"bodycentre","tailbase")
  dat$S7 <- GetDistances(t,"tailbase","hipr")
  dat$S8 <- GetDistances(t,"tailbase","hipl")
  dat$S9 <- GetDistances(t,"tailbase","tailcentre")
  dat$S10 <- GetDistances(t,"tailcentre","tailtip")
  dat$A1 <- GetAngleTotal(t,"tailbase","tailcentre","tailcentre","tailtip")
  dat$A2 <- GetAngleTotal(t,"hipr","tailbase","tailbase","hipl")
  dat$A3 <- GetAngleTotal(t,"tailbase","bodycentre","bodycentre","neck")
  dat$A4 <- GetAngleTotal(t,"bcr","bodycentre","bodycentre","bcl")
  dat$A5 <- GetAngleTotal(t,"bodycentre","neck","neck","headcentre")
  dat$A6 <- GetAngleTotal(t,"tailbase","bodycentre","neck","headcentre")
  dat$Ar1 <- GetPolygonAreas(t,c("tailbase","hipr","hipl"))
  dat$Ar2 <- GetPolygonAreas(t,c("hipr","hipl","bcl","bcr"))
  dat$Ar3 <- GetPolygonAreas(t,c("bcr","earr","earl","bcl"))
  dat$Ar4 <- GetPolygonAreas(t,c("earr","nose","earl"))
  dat$D1 <- apply(data.frame(a = GetDistanceToZoneBorder(t,"light","nose"), b = GetDistanceToZoneBorder(t,"dark","nose")), 1,FUN = function(x){min(x)})
  dat$D2 <- apply(data.frame(a = GetDistanceToZoneBorder(t,"light","neck"), b = GetDistanceToZoneBorder(t,"dark","neck")), 1,FUN = function(x){min(x)})
  dat$D3 <- apply(data.frame(a = GetDistanceToZoneBorder(t,"light","bodycentre"), b = GetDistanceToZoneBorder(t,"dark","bodycentre")), 1,FUN = function(x){min(x)})
  dat$D4 <- apply(data.frame(a = GetDistanceToZoneBorder(t,"light","tailbase"), b = GetDistanceToZoneBorder(t,"dark","tailbase")), 1,FUN = function(x){min(x)})
  dat <- as.data.frame(na_replace(dat)) 
  t$features <- dat
  
  return(t)
}

CreateSkeletonData_EPM <- function(t){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  dat <- data.frame(Ac1 = t$data[["nose"]]$acceleration)
  dat$Ac2 <- t$data[["headcentre"]]$acceleration
  dat$Ac3 <- t$data[["neck"]]$acceleration
  dat$Ac4 <- t$data[["earr"]]$acceleration
  dat$Ac5 <- t$data[["earl"]]$acceleration
  dat$Ac6 <- t$data[["bodycentre"]]$acceleration
  dat$Ac7 <- t$data[["bcl"]]$acceleration
  dat$Ac8 <- t$data[["bcr"]]$acceleration
  dat$Ac9 <- t$data[["hipl"]]$acceleration
  dat$Ac10 <- t$data[["hipr"]]$acceleration
  dat$Ac11 <- t$data[["tailbase"]]$acceleration
  dat$S1 <- GetDistances(t,"nose","headcentre")
  dat$S2 <- GetDistances(t,"headcentre","neck")
  dat$S3 <- GetDistances(t,"neck","bodycentre")
  dat$S4 <- GetDistances(t,"bodycentre","bcr")
  dat$S5 <- GetDistances(t,"bodycentre","bcl")
  dat$S6 <- GetDistances(t,"bodycentre","tailbase")
  dat$S7 <- GetDistances(t,"tailbase","hipr")
  dat$S8 <- GetDistances(t,"tailbase","hipl")
  dat$S9 <- GetDistances(t,"tailbase","tailcentre")
  dat$S10 <- GetDistances(t,"tailcentre","tailtip")
  dat$A1 <- GetAngleTotal(t,"tailbase","tailcentre","tailcentre","tailtip")
  dat$A2 <- GetAngleTotal(t,"hipr","tailbase","tailbase","hipl")
  dat$A3 <- GetAngleTotal(t,"tailbase","bodycentre","bodycentre","neck")
  dat$A4 <- GetAngleTotal(t,"bcr","bodycentre","bodycentre","bcl")
  dat$A5 <- GetAngleTotal(t,"bodycentre","neck","neck","headcentre")
  dat$A6 <- GetAngleTotal(t,"tailbase","bodycentre","neck","headcentre")
  dat$Ar1 <- GetPolygonAreas(t,c("tailbase","hipr","hipl"))
  dat$Ar2 <- GetPolygonAreas(t,c("hipr","hipl","bcl","bcr"))
  dat$Ar3 <- GetPolygonAreas(t,c("bcr","earr","earl","bcl"))
  dat$Ar4 <- GetPolygonAreas(t,c("earr","nose","earl"))
  dat <- as.data.frame(na_replace(dat)) 
  t$features <- dat
  
  return(t)
}


CreateSkeletonData_FST_v2 <- function(t){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  dat <- data.frame(Ac1 = t$data[["nose"]]$acceleration)
  dat$Ac2 <- t$data[["headcentre"]]$acceleration
  dat$Ac3 <- t$data[["neck"]]$acceleration
  dat$Ac4 <- t$data[["earr"]]$acceleration
  dat$Ac5 <- t$data[["earl"]]$acceleration
  dat$Ac6 <- t$data[["bodycentre"]]$acceleration
  dat$Ac7 <- t$data[["bcl"]]$acceleration
  dat$Ac8 <- t$data[["bcr"]]$acceleration
  dat$Ac9 <- t$data[["hipl"]]$acceleration
  dat$Ac10 <- t$data[["hipr"]]$acceleration
  dat$Ac11 <- t$data[["tailbase"]]$acceleration
  dat$S1 <- GetDistances(t,"nose","headcentre")
  dat$S2 <- GetDistances(t,"headcentre","neck")
  dat$S3 <- GetDistances(t,"neck","bodycentre")
  dat$S4 <- GetDistances(t,"bodycentre","bcr")
  dat$S5 <- GetDistances(t,"bodycentre","bcl")
  dat$S6 <- GetDistances(t,"bodycentre","tailbase")
  dat$S7 <- GetDistances(t,"tailbase","hipr")
  dat$S8 <- GetDistances(t,"tailbase","hipl")
  dat$S9 <- GetDistances(t,"tailbase","tailcentre")
  dat$S10 <- GetDistances(t,"tailcentre","tailtip")
  dat$A1 <- GetAngleTotal(t,"tailbase","tailcentre","tailcentre","tailtip")
  dat$A2 <- GetAngleTotal(t,"hipr","tailbase","tailbase","hipl")
  dat$A3 <- GetAngleTotal(t,"tailbase","bodycentre","bodycentre","neck")
  dat$A4 <- GetAngleTotal(t,"bcr","bodycentre","bodycentre","bcl")
  dat$A5 <- GetAngleTotal(t,"bodycentre","neck","neck","headcentre")
  dat$A6 <- GetAngleTotal(t,"tailbase","bodycentre","neck","headcentre")
  dat$Ar1 <- GetPolygonAreas(t,c("tailbase","hipr","hipl"))
  dat$Ar2 <- GetPolygonAreas(t,c("hipr","hipl","bcl","bcr"))
  dat$Ar3 <- GetPolygonAreas(t,c("bcr","earr","earl","bcl"))
  dat$Ar4 <- GetPolygonAreas(t,c("earr","nose","earl"))
  dat <- as.data.frame(na_replace(dat)) 
  t$features <- abs(dat)
  
  return(t)
}

CreateAccelerationFeatures<- function(t){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  dat <- data.frame(Ac1 = t$data[["nose"]]$acceleration)
  dat$Ac2 <- t$data[["headcentre"]]$acceleration
  dat$Ac3 <- t$data[["neck"]]$acceleration
  dat$Ac4 <- t$data[["earr"]]$acceleration
  dat$Ac5 <- t$data[["earl"]]$acceleration
  dat$Ac6 <- t$data[["bodycentre"]]$acceleration
  dat$Ac7 <- t$data[["bcl"]]$acceleration
  dat$Ac8 <- t$data[["bcr"]]$acceleration
  dat$Ac9 <- t$data[["hipl"]]$acceleration
  dat$Ac10 <- t$data[["hipr"]]$acceleration
  dat$Ac11 <- t$data[["tailbase"]]$acceleration
  dat <- as.data.frame(dat) 
  t$features <- abs(dat)
  
  return(t)
}

ZscoreNormalizeFeatures <- function(t, omit = NULL, include = NULL, type = "mean"){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(t$features)){
    stop("Object has no feature data")
  }
  change <- setdiff(names(t$features),omit)
  if(!is.null(include)){
    change <- intersect(change, include)
  }
  for(i in change){
    if(type == "median"){
      t$features[i] <- NormalizeZscore_median(t$features[i])
    }else if(type == "mean"){
      t$features[i] <- NormalizeZscore(t$features[i])
    }else{
      warning("invalid normalization method")
    }
  }
  return(t)
}

PlotDensityPaths <- function(t,points,SDcutoff = 4, Title = "density path"){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  out <- NULL
  
  SDPlot <- function(x,nSD){
    ifelse((mean(x) - x) / sd(x) > -nSD, x, mean(x) + nSD * sd(x))
  }
  
  for(i in points){
    
    data_plot <- t$data[[i]]
    xbreaks <- seq(floor(min(data_plot$x)), ceiling(max(data_plot$x)), by = 0.1)
    ybreaks <- seq(floor(min(data_plot$y)), ceiling(max(data_plot$y)), by = 0.1)
    data_plot$latbin <- xbreaks[cut(data_plot$x, breaks = xbreaks, labels=F)]
    data_plot$longbin <- ybreaks[cut(data_plot$y, breaks = ybreaks, labels=F)]
    out[[i]] <- ggplot(data = data_plot, aes(x,y)) + 
      stat_density_2d(data = data_plot, aes(latbin,longbin, fill=..density..), geom = "raster", contour = FALSE) + 
      scale_fill_gradient(name = "Time Density", low = "blue", high = "yellow") +
      geom_path(data = data_plot, aes(x,y, color = SDPlot((speed * t$fps),SDcutoff))) + 
      theme_bw() + 
      ggtitle(paste(i,Title, sep = " ")) + 
      scale_color_gradient2(name = paste("speed (",t$distance.units,"/s)",sep = ""), high = "white", low="black", mid = "black")
  }
  return(out)
}


SmoothLabel <- function(x, integration_period){
  types <-unique(x)
  mat <- NULL
  for (i in types){
    mat <- cbind(mat,periodsum(x == i, integration_period))
  }
  c <- apply(mat,1,FUN = which.max)
  return(types[c])
}

AddZonesToPlots <- function(p,z){
  for(i in 1:length(p)){
    for(j in z){
      p[[i]] <- p[[i]] + geom_path(data=j[c(1:nrow(j),1),],aes(x,y))
    }
  }
  return(p)
}

PlotZones <- function(t, zones = NULL){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(t$zones)){
    warning("Zones not defined")
    return(NULL)
  }
  if(is.null(zones)){
    zones <- names(t$zones)
  }
  p <- ggplot()
  for(i in zones){
    dat <- t$zones[[i]]
    p <- p + geom_path(data=dat[c(1:nrow(dat),1),],aes(x,y))
  }
  return(p + theme_bw())
}

AddZones <- function(t,z){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(t$zones)){
    t$zones <- list()
    t$zones.invert <- list()
  }
  
  for(i in names(z)){
    t$zones[[i]] <- t$median.data[as.character(z[z[,i]!= "",i]),c("x","y")]
    t$zones.invert[[i]] <- FALSE
  }
  return(t)
}

AddBinData <- function(t, bindat = NULL, unit = "frame", binlength = NULL){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(!is.null(bindat)){
    if(unit == "second"){
      bindat$from <- bindat$from * t$fps
      bindat$to <- bindat$to * t$fps
    }
    if(unit == "minute"){
      bindat$from <- bindat$from * t$fps * 60
      bindat$to <- bindat$to * t$fps * 60
    }
    t$bins <- bindat
  }else if(!is.null(binlength)){
    binlength <- binlength * ifelse(unit != "frame",t$fps,1) * ifelse(unit == "minute",60,1)
    t$bins <- NULL
    for(i in 1:ceiling(length(t$frames) / binlength)){
      t$bins <- rbind(t$bins, data.frame(bin = paste("bin",i,sep="."), 
                                         from = t$frames[(i-1)*binlength + 1], 
                                         to = t$frames[min(i*binlength,length(t$frames))]))
    }
  }else{
    warning("To add bin data either a the bindata or binlenght has to be added")
  }
  return(t)
}

BinAnalysis <- function(t,FUN, ...){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(t$bins)){
    warning("No bins defined, can not perform bin analysis")
    return(NULL)
  }
  
  Report <- NULL
  for(i in 1:nrow(t$bins)){
    tb <- CutTrackingData(t, keep.frames = t$bins[i,"from"]:t$bins[i,"to"])
    Report <- rbind.fill(Report,data.frame(bin = t$bins[i,"bin"], FUN(tb,...)$Report))
  }
  return(Report)
}

FSTAnalysis <- function(t, cutoff_floating, integration_period = 0, points, Object){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(sum(t$point.info$PointType == Object) == 0){
    stop("Object to be tracked not available in point info")
  }
  for(i in 1:length(t$data)){
    t$data[[i]]$delta_x <- integratevector(t$data[[i]]$x)
    t$data[[i]]$delta_y <- integratevector(t$data[[i]]$y)
    t$data[[i]]$speed <- sqrt(t$data[[i]]$delta_x ^2 + t$data[[i]]$delta_y^2)
    t$data[[i]]$acceleration <- integratevector(t$data[[i]]$speed)
  }
  
  temp <- t$data[as.character(t$point.info[t$point.info$PointType == Object,"PointName"])]
  acc <- NULL
  for(i in 1:length(temp)){
    acc <- cbind(acc, i = temp[[i]]$acceleration)
  }
  
  t$object <- list()
  t$object$movement <- abs(apply(acc,1,FUN = mean))
  t$object$is.floating <- avgmean(t$object$movement,integration_period) < cutoff_floating
  t$labels$cutoff.floating <- ifelse(t$object$is.floating == 1, "Floating","None")
  
  t$Report <- list()
  t$Report[["time.floating"]] <- sum(t$object$is.floating) / t$fps
  t$Report[["total.time"]] <- length(t$object$is.floating) / t$fps
  t$Report[["percentage.floating"]] <- sum(t$object$is.floating) / length(t$object$is.floating) * 100
  
  for(k in points){
    if(!k %in% names(t$data)){
      stop(paste("point", k, "not vaild", sep = " "))
    }
    t$Report[[paste(k, "raw.distance", sep = ".")]] <- sum(t$data[[k]]$speed, na.rm = T)
    t$Report[[paste(k, "raw.speed", sep = ".")]] <- mean(t$data[[k]]$speed, na.rm = T) * t$fps
    t$Report[[paste(k, "distance.swiming", sep = ".")]] <- sum(t$data[[k]]$speed[!t$object$is.floating], na.rm = T)
    t$Report[[paste(k, "speed.swiming", sep = ".")]] <- mean(t$data[[k]]$speed[!t$object$is.floating], na.rm = T) * t$fps
  }
  
  return(t)
}  

HeadAngleAnalysis <- function(t, points = c("tailbase","neck","neck","nose"), angle_cutoff, integration_period){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  t$object <- list()
  t$object$head.angle.CW <- GetAngleClockwise(t,points[1],points[2],points[3],points[4])
  t$object$head.angle.total <- GetAngleTotal(t,points[1],points[2],points[3],points[4])
  t$object$head.angle.CW.degree <- (t$object$head.angle.CW + pi) * 180 / pi - 180
  t$object$head.angle.total.degree <- t$object$head.angle.total * 180 / pi
  t$object$head.tilted.CW <- t$object$head.angle.CW.degree > angle_cutoff
  t$object$head.tilted.CCW <- t$object$head.angle.CW.degree < -angle_cutoff
  
  t$Report <- list()
  t$Report[["average.head.angle.CW"]] <- mean(t$object$head.angle.CW.degree)
  t$Report[["time.head.tilted.CW"]] <- sum(avgbool(t$object$head.tilted.CW, integration_period)) / t$fps
  t$Report[["time.head.tilted.CCW"]] <- sum(avgbool(t$object$head.tilted.CCW, integration_period)) / t$fps
  t$Report[["time.total"]] <- length(t$object$head.tilted.CW) / t$fps
  
  
  return(t)
}

GetAngleClockwise <- function(t,p1,p2,p3,p4){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(sum(c(p1,p2,p3,p4) %in% names(t$data)) != 4){
    warning(paste("Points for angle measurement not available in Tracking data:",paste(setdiff(c(p1,p2,p3,p4), names(t$data)),collapse = " "), sep = " "))
    return(NULL)
  }
  
  a1 <- t$data[[p1]]
  a2 <- t$data[[p2]]
  b1 <- t$data[[p3]]
  b2 <- t$data[[p4]]
  ax <- a1$x - a2$x
  ay <- a1$y - a2$y
  bx <- b1$x - b2$x
  by <- b1$y - b2$y
  dot <- ax*bx + ay*by      # dot product between [x1, y1] and [x2, y2]
  det <- ax*by - ay*bx      # determinant
  res <- atan2(det, dot)  # atan2(y, x) or atan2(sin, cos)
  return(res)
}

GetAngleTotal <- function(t,p1,p2,p3,p4){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(sum(c(p1,p2,p3,p4) %in% names(t$data)) != 4){
    warning(paste("Points for angle measurement not available in Tracking data:",paste(setdiff(c(p1,p2,p3,p4), names(t$data)),collapse = " "), sep = " "))
    return(NULL)
  }
  
  a1 <- t$data[[p1]]
  a2 <- t$data[[p2]]
  b1 <- t$data[[p3]]
  b2 <- t$data[[p4]]
  ax <- a1$x - a2$x
  ay <- a1$y - a2$y
  bx <- b1$x - b2$x
  by <- b1$y - b2$y
  
  res <- rep(0,nrow(a1))
  for(i in 1:nrow(a1)){
    res[i] <- acos(sum(c(ax[i],ay[i])*c(bx[i],by[i])) / ( sqrt(sum(c(ax[i],ay[i]) * c(ax[i],ay[i]))) * sqrt(sum(c(bx[i],by[i]) * c(bx[i],by[i]))) ) )
  }
  return(res)
}

CalculateTransitions <- function(x,integration_period){
  x <- avgbool(x,integration_period)
  sum(append(0, (x[2:length(x)]!= x[1:length(x)-1])))
}

avgmean <- function(x, window){
  res <- rep(0, length(x))
  for(i in 1:length(x)){
    res[i] <- mean(x[max(0,i-window):min(length(x), i + window)], na.rm = T)
  }
  return(res)
}

avgmedian <- function(x, window){
  res <- rep(0, length(x))
  for(i in 1:length(x)){
    res[i] <- median(x[max(0,i-window):min(length(x), i + window)], na.rm = T)
  }
  return(res)
}

periodsum <- function(x, window){
  res <- rep(0, length(x))
  for(i in 1:length(x)){
    res[i] <- sum(x[max(0,i-window):min(length(x), i + window)], na.rm = T)
  }
  return(res)
}

CreateTrainingSet <- function(t, integration_period, label.group = "manual"){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(t$labels[[label.group]])){
    stop("manual labels or specified label group does not exist")
  }
  if(is.null(t$features)){
    stop("no feature data available. can not create training set")
  }
  x <- t$features
  x_window <- x[1:(nrow(x) - 2*integration_period),]
  labs <- paste(names(x),-integration_period, sep = ",")
  if(integration_period > 0){
    for(i in (-integration_period + 1):integration_period){
      labs <- append(labs, paste(names(x),i, sep = ",") )
      x_window <- cbind(x_window, x[(integration_period + i + 1):(nrow(x) - integration_period + i),])
    }
  }
  
  colnames(x_window) <- labs
  t$train_x <- as.matrix(x_window)
  t$train_y <- t$labels[[label.group]][(integration_period + 1):(length(t$labels[[label.group]]) - integration_period)]
  t$ml_integration <- integration_period
  return(t)
}


periodogram <- function(x, mean.x = mean(x)) { # simple periodogram
  n <- length(x)
  x <- unclass(x) - mean.x
  Mod(fft(x))[2:(n%/%2 + 1)]^2 / (2*pi*n) # drop I(0)
}

CreateTrainingSetFFT <- function(t, integration_period, label.group = "manual"){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(t$labels[[label.group]])){
    stop("manual labels or specified label group does not exist")
  }
  if(is.null(t$features)){
    stop("no feature data available. can not create training set")
  }
  if(integration_period < 2){
    stop("to low temporal resolution for fast fourier transform, integration period of at least 2")
  }
  
  x <- t$features
  x_window <- matrix(0, nrow = nrow(x) - 2 * integration_period, ncol = ncol(x)* integration_period)
  for(i in 1:nrow(x_window)){
    x_window[i,] <- c(apply(x[i:(i+2*integration_period),], 2, FUN = periodogram))
  }
  
  labs <- NULL
  for(i in names(t$features)){
    labs <- append(labs,paste(rep(i,integration_period), 1:integration_period, sep = ","))
  }
  # if(integration_period > 0){
  #   for(i in (-integration_period):integration_period){
  #     labs <- append(labs, paste(names(x),i, sep = ",") )
  #     x_window <- cbind(x_window, x[(integration_period + i + 1):(nrow(x) - integration_period + i),])
  #   }
  # }
  colnames(x_window) <- labs
  
  t$train_x <- as.matrix(x_window)
  t$train_y <- t$labels[[label.group]][(integration_period + 1):(length(t$labels[[label.group]]) - integration_period)]
  t$ml_integration <- integration_period
  return(t)
}


CreateTrainingSetLSTM <- function(t, integration_period, label.group = "manual"){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(t$labels[[label.group]])){
    stop("manual labels or specified label group does not exist")
  }
  if(is.null(t$features)){
    stop("no feature data available. can not create training set")
  }
  if(integration_period < 2){
    stop("to low temporal resolution for fast fourier transform, integration period of at least 2")
  }
  x <- t$features
  x_window <- array(rep(0,(nrow(x) - 2 * integration_period) * (2 * integration_period + 1) * ncol(x)), c(nrow(x) - 2 * integration_period, 2 * integration_period + 1, ncol(x))); 
  
  for(i in 1:(nrow(x) - 2 * integration_period)){
    x_window[i,,] <- as.matrix(x[(i):(i + 2*integration_period),])
  }
    
  dimnames(x_window)[[3]] <- colnames(x)
  t$train_x <- x_window
  t$train_y <- t$labels[[label.group]][(integration_period + 1):(length(t$labels[[label.group]]) - integration_period)]
  t$ml_integration <- integration_period
  return(t)
}


CreateTrainingSetConv1d <- function(t, integration_period, label.group = "manual"){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(t$labels[[label.group]])){
    stop("manual labels or specified label group does not exist")
  }
  if(is.null(t$features)){
    stop("no feature data available. can not create training set")
  }
  if(integration_period < 2){
    stop("to low temporal resolution for fast fourier transform, integration period of at least 2")
  }
  x <- t$features
  x_window <- array(rep(0,(nrow(x) - 2 * integration_period) * (2 * integration_period + 1) * ncol(x)), c(nrow(x) - 2 * integration_period, 2 * integration_period + 1, ncol(x))); 
  
  for(i in 1:(nrow(x) - 2 * integration_period)){
    x_window[i,,] <- as.matrix(x[(i):(i + 2*integration_period),])
  }
  
  dimnames(x_window)[[3]] <- colnames(x)
  t$train_y <- t$labels[[label.group]][(integration_period + 1):(length(t$labels[[label.group]]) - integration_period)]
  t$ml_integration <- integration_period
  return(t)
}


CreateTestSetLSTM <- function(t, integration_period, label.group = "manual"){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(t$labels[[label.group]])){
    stop("manual labels or specified label group does not exist")
  }
  if(is.null(t$features)){
    stop("no feature data available. can not create training set")
  }
  if(integration_period < 2){
    stop("to low temporal resolution for fast fourier transform, integration period of at least 2")
  }
  x <- t$features
  x_window <- array(rep(0,(nrow(x) - 2 * integration_period) * (2 * integration_period + 1) * ncol(x)), c(nrow(x) - 2 * integration_period, 2 * integration_period + 1, ncol(x))); 
  
  for(i in 1:(nrow(x) - 2 * integration_period)){
    x_window[i,,] <- as.matrix(x[(i):(i + 2*integration_period),])
  }
  
  dimnames(x_window)[[3]] <- colnames(x)
  t$train_x <- x_window
  t$ml_integration <- integration_period
  return(t)
}

CreateTestSetConv1d <- function(t, integration_period, label.group = "manual"){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(t$labels[[label.group]])){
    stop("manual labels or specified label group does not exist")
  }
  if(is.null(t$features)){
    stop("no feature data available. can not create training set")
  }
  if(integration_period < 2){
    stop("to low temporal resolution for fast fourier transform, integration period of at least 2")
  }
  x <- t$features
  x_window <- array(rep(0,(nrow(x) - 2 * integration_period) * (2 * integration_period + 1) * ncol(x)), c(nrow(x) - 2 * integration_period, 2 * integration_period + 1, ncol(x))); 
  
  for(i in 1:(nrow(x) - 2 * integration_period)){
    x_window[i,,] <- as.matrix(x[(i):(i + 2*integration_period),])
  }
  
  dimnames(x_window)[[3]] <- colnames(x)
  t$ml_integration <- integration_period
  return(t)
}

CreateTestSetFFT <- function(t, integration_period){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(integration_period < 2){
    stop("to low temporal resolution for fast fourier transform, integration period of at least 2")
  }
  
  x <- t$features
  x_window <- matrix(0, nrow = nrow(x) - 2 * integration_period, ncol = ncol(x)* integration_period)
  for(i in 1:nrow(x_window)){
    x_window[i,] <- c(apply(x[i:(i+2*integration_period),], 2, FUN = periodogram))
  }
  
  labs <- NULL
  for(i in names(t$features)){
    labs <- append(labs,paste(rep(i,integration_period), 1:integration_period, sep = ","))
  }
  # if(integration_period > 0){
  #   for(i in (-integration_period):integration_period){
  #     labs <- append(labs, paste(names(x),i, sep = ",") )
  #     x_window <- cbind(x_window, x[(integration_period + i + 1):(nrow(x) - integration_period + i),])
  #   }
  # }
  colnames(x_window) <- labs
  
  
  t$train_x <- as.matrix(x_window)
  t$ml_integration <- integration_period
  return(t)
}

CreateTestSet <- function(t, integration_period){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  x <- t$features
  x_window <- x[1:(nrow(x) - 2*integration_period),]
  
  labs <- paste(names(x),-integration_period, sep = ",")
  if(integration_period > 0){
    for(i in (-integration_period + 1):integration_period){
      labs <- append(labs, paste(names(x),i, sep = ",") )
      x_window <- cbind(x_window, x[(integration_period + i + 1):(nrow(x) - integration_period + i),])
    }
  }
  colnames(x_window) <- labs
  t$train_x <- as.matrix(x_window)
  t$ml_integration <- integration_period
  return(t)
}

ClassifyBehaviors <- function(t,model, model_parameters){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  t <- CreateTestSet(t, model_parameters$integration_period)
  t$labels$classifications <- model %>% predict(t$train_x)
  t$labels$classifications <- c(rep(NA,model_parameters$integration_period), 
                                model_parameters$Feature_names[apply(t$labels$classifications,MARGIN = 1, FUN = which.max)],
                                rep(NA,model_parameters$integration_period))  
  return(t)
}

ClassifyBehaviorsFFT <- function(t,model, model_parameters){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  t <- CreateTestSetFFT(t, model_parameters$integration_period)
  t$labels$classifications <- model %>% predict(t$train_x)
  t$labels$classifications <- c(rep(NA,model_parameters$integration_period), 
                                model_parameters$Feature_names[apply(t$labels$classifications,MARGIN = 1, FUN = which.max)],
                                rep(NA,model_parameters$integration_period))  
  return(t)
}

ClassifyBehaviorsLSTM <- function(t,model, model_parameters){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  t <- CreateTestSetLSTM(t, model_parameters$integration_period)
  temp <- t$train_x[1:(floor(nrow(t$train_x) / model_parameters$batchsize) * model_parameters$batchsize),,]
  t$labels$classifications <- model %>% predict(temp)
  t$labels$classifications <- c(rep(NA,model_parameters$integration_period), 
                                model_parameters$Feature_names[apply(t$labels$classifications,MARGIN = 1, FUN = which.max)],
                                rep(NA,nrow(t$train_x) %% model_parameters$batchsize),
                                rep(NA,model_parameters$integration_period))
  return(t)
}

ClassifyBehaviorsConv1d <- function(t,model, model_parameters){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  t <- CreateTestSetConv1d(t, model_parameters$integration_period)
  temp <- t$train_x[1:(floor(nrow(t$train_x) / model_parameters$batchsize) * model_parameters$batchsize),,]
  t$labels$classifications <- model %>% predict(temp)
  t$labels$classifications <- c(rep(NA,model_parameters$integration_period), 
                                model_parameters$Feature_names[apply(t$labels$classifications,MARGIN = 1, FUN = which.max)],
                                rep(NA,nrow(t$train_x) %% model_parameters$batchsize),
                                rep(NA,model_parameters$integration_period))
  return(t)
}


PlotLabels <- function(t, p.size = 2){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  dat <- NULL
  
  if(length(t$labels) == 0){
    warning("No labeling data available. can not plot")
    return(NULL)
  }
  
  for(i in names(t$labels)){
    dat <- rbind(dat,data.frame(seconds = t$seconds, behavior =  t$labels[[i]], type = i))
  }
  ggplot(data = na.omit(dat),aes(seconds, behavior, color = behavior)) + 
    geom_point(size = p.size, shape = 124) + 
    facet_grid(type~., scales = "free_y") +
    theme_bw()
}

PlotZoneVisits <- function(t, points, zones = NULL, p.size = 2){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(t$zones)){
    warning("no zones defined")
    return(NULL)
  }
  if(is.null(zones)){
    zones <- names(t$zones)
  }
  if(length(setdiff(zones,names(t$zones))) > 0){
    warning("invalid zones")
    return(NULL)
  }
  if(length(setdiff(points,names(t$data))) > 0){
    warning("invalid points")
    return(NULL)
  }
  
  dat <- NULL
  for(j in points){
    for(i in zones){
      dat <- rbind(dat,data.frame(seconds = t$seconds, zone = ifelse(IsInZone(t,j,i,t$zones.invert[[i]]),i,NA), type = "automatic", points = j))
    }
  }
  
  ggplot(data = na.omit(dat),aes(seconds, zone, color = zone)) + 
    geom_point(size = p.size, shape = 124) + 
    facet_grid(points~.) + theme_bw()
}

PlotPointData <- function(t, points = NULL, from = NULL, to = NULL, unit = "frame", type = NULL){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(points) & is.null(type)){
    points <- names(t$data)
  }
  
  if(unit == "second"){
    if(!is.null(from)){
      from <- t$frames[which(t$seconds >= from)[1]]
    }
    if(!is.null(to)){
      to <- t$frames[which(t$seconds >= to)[1]]
    }
  }
  if(is.null(from)){
    from = min(t$frames)
  }
  if(is.null(to)){
    to = max(t$frames)
  }
  
  if(!is.null(type)){
    points <- t$point.info[t$point.info$PointType == type,"PointName"]
  }
  
  range <- from:to
  
  p <- ggdraw()
  dim <- ceiling(sqrt(length(points)))
  nplot <- 0
  
  for(i in points){
    p <- p + draw_plot(ggplot(data = t$data[[i]][t$data[[i]]$frame %in% range,], aes(x,y, color = likelihood)) + geom_path() + ggtitle(i) + xlab(paste("x /",t$distance.units,sep = " ")) + ylab(paste("y /",t$distance.units,sep = " ")) + theme_bw(), 
                       x = (nplot %% dim / dim),
                       y = ((dim - 1)/ dim) - floor(nplot / dim) / dim, 
                       width = 1/dim, 
                       height = 1/dim)
    nplot <- nplot + 1
  }
  
  return(p)
}

RunPipeline <- function(files, path, FUN){
  out <- list()
  for(j in files){
    out[[paste(j)]] <- FUN(paste(path,j,sep = ""))
  }
  return(out)
}

CombineTrainingsData <- function(ts, shuffle =TRUE){
  if(IsTrackingData(ts)){
    temp <- ts
    ts <- list()
    ts[[paste(temp$filename)]] <- temp
  }
  
  train_x <- list()
  train_y <- NULL
  for(i in names(ts)){
    if(is.null(ts[[i]]$train_x) | is.null(ts[[i]]$train_y)){
      warning(paste("File:",i,"is missing trainings data. data for this file not included"))
    }else{
      train_x[[i]] <- ts[[i]]$train_x
      train_y <- append(train_y,ts[[i]]$train_y)
    }
  }
  train_x <- do.call(rbind, train_x)
  out <- PrepareMLData(train_x,train_y, shuffle)
  out$parameters$integration_period <- ts[[1]]$ml_integration
  return(out)
}

PrepareMLData <- function(x_train, y_train, shuffle = TRUE){
  if(nrow(x_train) != length(y_train)){
    warning("Unequal size of x_train and y_train")
    return(NULL)
  }
  out <- list()
  parameters <- list()
  x_train <- as.matrix(x_train)
  out$parameters$N_input <- ncol(x_train)
  out$parameters$N_features <- length(unique(y_train))
  out$parameters$Feature_names <- levels(as.factor(y_train))
  
  y_train_cat <- to_categorical(-1 + as.integer(as.factor(y_train)))
  
  if(shuffle){
    new_order <- sample(1:nrow(y_train_cat))
    x_train <- x_train[new_order,]
    y_train_cat <- y_train_cat[new_order,]
  }
  out$train_x <- x_train
  out$train_y <- y_train_cat
  
  return(out)
}

CombineTrainingsDataLSTM <- function(ts, shuffle =TRUE){
  if(IsTrackingData(ts)){
    temp <- ts
    ts <- list()
    ts[[paste(temp$filename)]] <- temp
  }
  
  train_x <- NULL
  train_y <- NULL
  for(i in names(ts)){
    if(is.null(ts[[i]]$train_x) | is.null(ts[[i]]$train_y)){
      warning(paste("File:",i,"is missing trainings data. data for this file not included"))
    }else{
      train_x <- abind(train_x,ts[[i]]$train_x, along = 1)
      train_y <- append(train_y,ts[[i]]$train_y)
    }
  }
  out <- PrepareMLDataLSTM(train_x,train_y, shuffle)
  out$parameters$integration_period <- ts[[1]]$ml_integration
  return(out)
}

PrepareMLDataLSTM <- function(x_train, y_train, shuffle = TRUE){
  if(nrow(x_train) != length(y_train)){
    warning("Unequal size of x_train and y_train")
    return(NULL)
  }
  out <- list()
  parameters <- list()
  out$parameters$N_input <- ncol(x_train)
  out$parameters$N_features <- length(unique(y_train))
  out$parameters$Feature_names <- levels(as.factor(y_train))
  y_train_cat <- to_categorical(-1 + as.integer(as.factor(y_train)))
  
  if(shuffle){
    new_order <- sample(1:nrow(y_train_cat))
    x_train <- x_train[new_order,,]
    y_train_cat <- y_train_cat[new_order,]
  }
  out$train_x <- x_train
  out$train_y <- y_train_cat
  return(out)
}

SmoothLabels <- function(t, integration_period){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(length(t$labels) == 0){
    warning("No labels present. Returning original object")
    return(t)
  }
  for(i in names(t$labels)){
    t$labels[[i]] <- SmoothLabel(t$labels[[i]], integration_period)
  }
  return(t)
}

UnsupervisedClusteringKmeans <- function(ts, N_clusters = 20, Z_score_Normalize = TRUE, dimensions = NULL){
  if(IsTrackingData(ts)){
    print("single file detected. Runing kmeans in single file mode")
    if(is.null(ts$train_x)){
      warning("No training or testing data present. Returning original data")
      return(ts)
    }
    if(Z_score_Normalize){
      test <- kmeans(NormalizeZscore(ts$train_x),centers = N_clusters)
    }
    else{
      test <- kmeans(ts$train_x,centers = N_clusters)
    }
    ts$labels$unsupervised <- c(rep(NA,ts$ml_integration),as.character(test$cluster),rep(NA,ts$ml_integration))
    return(ts)
  }
  
  allx <- list()
  id <- NULL
  print("multiple files detected. Runing kmeans in multi file mode")
  for(j in names(ts)){
    if(is.null(ts[[j]]$train_x)){
      warning(paste("File:",j, "does not cotain any training data. Returning original data"))
      return(ts)
    }
    allx[[j]] <- ts[[j]]$train_x
    id <- append(id, rep(paste(j),nrow(ts[[j]]$train_x)))
    #ts[[j]]$train_x <- NULL
  }
  
  print("fusing data")
  allx <- do.call(rbind, allx)
  
  if(Z_score_Normalize){
    print("Zscore Normalize data")
    for(i in 1:ncol(allx)){
      x <- allx[,i]
      allx[,i] <- (x - mean(x))/(sd(x))
    }
  }
  
  if(!is.null(dimensions)){
    print("Run dimension reduction using singular value decomposition")
    allx <- svd(allx, nu = dimensions, nv = 0)$u
  }
  
  print("Runing kmeans clustering")
  #allx <- as.big.matrix(allx)
  test <- bigkmeans(allx,centers = N_clusters)
  
  for(j in names(ts)){
    clust <- test$cluster[id == j]
    ts[[j]]$labels$unsupervised <- c(rep(NA,ts[[j]]$ml_integration),as.character(clust),rep(NA,ts[[j]]$ml_integration))
  }
  return(ts)
}


#' Plots the zone selection of a zone in an object of type TrackingData
#' 
#' @param t a objects of type TrackingData
#' @param point string. name of the point to be plotted
#' @param zones string or vector of strings. name of zones to be selected
#' @param invert boolean. a checkt to determine if the zone or a inversion of the zone should be used. defaults to FALSE
#' @return a plot
#' @examples
#' MultiFileReport(TrackingAll)
#'
PlotZoneSelection <- function(t,point,zones, invert = FALSE){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(!point %in% names(t$data)){
    warning("Invalid point")
    return(NULL)
  }
  if(!sum(zones %in% names(t$zones))){
    warning("Invalid zone(s)")
    return(NULL)
  }
  
  dat <- t$data[[point]]
  in.zone <- rep(FALSE,nrow(dat))
  for(i in t$zones[zones]){
    in.zone <- in.zone | (point.in.polygon(dat$x,dat$y,i$x,i$y) == 1)
  }
  if(invert){
    in.zone <- !in.zone
  }
  p <- ggplot(dat,aes(x,y, color = in.zone)) + geom_point() + theme_bw()
  return(p)
}

#' Creates a combined report of all files in a list of TrackingData objects
#' 
#' @param ts list of objects of type TrackingData
#' @return a data.frame
#' @examples
#' MultiFileReport(TrackingAll)
#'
MultiFileReport <- function(ts){
  if(IsTrackingData(ts)){
    stop("Expected a list() of TrackingData objects. You entered a single object")
  }
  out <- list()
  for(i in names(ts)){
    if(IsTrackingData(ts[[i]])){
      if(is.null(ts[[i]]$Report)){
        warning(paste("Object",i,"Does not contain any Report. omitting", sep = " "))
      }else{
        out <- rbindlist(list(out,append(c(file = ts[[i]]$filename), ts[[i]]$Report)),use.names = TRUE, fill = TRUE,idcol = F)
      }
    }else{
      warning(paste("List contains an element that is not of type TrackingData:",i,".No report produced for these", sep = " "))
    }
  }
  return(data.frame(out))
}

#' Performs a specified analysis on a list() of TrackingData objects and produces a final report
#' 
#' @param ts list of objects of type TrackingData
#' @param FUN An analysis function that should be applied to bins
#' @param ... any paramteres that are required for the function FUN
#' @return a data.frame
#' @examples
#' MultiFileBinanalysis(TrackingAll, FUN = OFTAnalysis, movement_cutoff = 5,integration_period = 5,points = "bodycentre")
#'
MultiFileBinanalysis <- function(ts, FUN, ...){
  if(IsTrackingData(ts)){
    stop("Expected a list() of TrackingData objects. You entered a single object")
  }
  out <- NULL
  for(i in names(ts)){
    if(IsTrackingData(ts[[i]])){
      binrep <- BinAnalysis(ts[[i]], FUN = FUN, ...)
      out <- rbind.fill(out, data.frame(file = i, binrep))
    }else{
      warning(paste("Object",i,"is not of type Tracking data. omiting from Analysis", sep = " "))
    }
  }
  return(out)
}

#' Creates an density path plot pdf for one or multiple objects of type TrackingData
#' 
#' @param ts list or single object of type TrackingData
#' @param points A string or vector of strings. name of the point(s) to be plotted
#' @param filename A string. the name of the pdf file. defaults to "DensityPathMulti"
#' @param width a numeric value. width of the plot in inch. defaults to 10
#' @param height a numeric value. height of the plot in inch. defaults to 8
#' @param add_zones a boolean check. should zones be plotted to the density path? requires all objects to have zones. default to FALSE
#' @param selected_zones a string or vector of strings. zone names of the zones to be plotted. default to NULL (= plot all)
#' @param ... any parameter(s) that is used for the function PlotZoneVisits()
#' @return NULL. creates a pdf file in the working directory
#' @examples
#' PlotZoneVisits.Multi.PDF(TrackingAll, filename = "MyZoneVisits")
#'
PlotDensityPaths.Multi.PDF <- function(ts, points, filename = "DensityPathMulti",width = 10, height = 8, add_zones = FALSE, selected_zones = NULL, ...){
  if(IsTrackingData(ts)){
    x <- ts
    ts <- list()
    ts[[paste(x$filename)]] <- x
  }
  out <- list()
  for(i in names(ts)){
    if(IsTrackingData(ts[[i]])){
      ps <- PlotDensityPaths(ts[[i]],points = points, Title = i, ...)
      
      if(add_zones){
        if(!is.null(ts[[i]]$zones)){
          if(is.null(selected_zones)){
            zones <- names(ts[[i]]$zones)
          }else{
            zones <- intersect(names(ts[[i]]$zones), selected_zones)
          }
          if(length(zones) == 0){
            warning(paste("in file:", i, "none of the indicated zones found in data", sep = " "))
          }
          else{
            ps <- AddZonesToPlots(ps,ts[[i]]$zones[zones])
          }
        }else{
          warning(paste("can not add zones to",i, ".Does not have zones defined", sep = " "))
        }
      }
      out <- append(out,ps)
    }else{
      warning(paste("List contains an element that is not of type TrackingData:",i,".No plot produced for these", sep = " "))
    }
  }
  pdf(paste(filename,".pdf",sep = ""), width = width, height = height)
  for(i in out){
    print(i)
  }
  dev.off()
  return(NULL)
}

#' Creates an zone visit plot pdf for one or multiple objects of type TrackingData
#' 
#' @param ts list or single object of type TrackingData
#' @param filename A string. the name of the pdf file. defaults to "ZoneVisitsMulti"
#' @param width a numeric value. width of the plot in inch. defaults to 10
#' @param height a numeric value. height of the plot in inch. defaults to 8
#' @param ... any parameter(s) that is used for the function PlotZoneVisits()
#' @return NULL. creates a pdf file in the working directory
#' @examples
#' PlotZoneVisits.Multi.PDF(TrackingAll, filename = "MyZoneVisits")
#'
PlotZoneVisits.Multi.PDF <- function(ts, points,filename = "ZoneVisitsMulti", width = 10, height = 8,...){
  if(IsTrackingData(ts)){
    x <- ts
    ts <- list()
    ts[[paste(x$filename)]] <- x
  }
  
  pdf(paste(filename,".pdf",sep = ""), width = width, height = height)
  for(i in names(ts)){
    print(PlotZoneVisits(ts[[i]], points,...) + ggtitle(paste(i)))
  }
  dev.off()
  return(NULL)
}

#' Creates an labels plot pdf for one or multiple objects of type TrackingData
#' 
#' @param ts list or single object of type TrackingData
#' @param filename A string. the name of the pdf file. defaults to "LabelsMulti"
#' @param width a numeric value. width of the plot in inch. defaults to 10
#' @param height a numeric value. height of the plot in inch. defaults to 8
#' @param ... any parameter(s) that is used for the function PlotLabels()
#' @return NULL. creates a pdf file in the working directory
#' @examples
#' PlotLabels.Multi.PDF(TrackingAll, filename = "MyLabels")
#'
PlotLabels.Multi.PDF <- function(ts, filename = "LabelsMulti", width = 10, height = 8, ...){
  if(IsTrackingData(ts)){
    x <- ts
    ts <- list()
    ts[[paste(x$filename)]] <- x
  }
  
  pdf(paste(filename,".pdf",sep = ""), width = width, height = height)
  for(i in names(ts)){
    print(PlotLabels(ts[[i]],...) + ggtitle(paste(i)))
  }
  dev.off()
  return(NULL)
}

#' Creates an overview plot for a object of type TrackingData
#' 
#' @param t object of type TrackingData
#' @param point a string indicating the name of the point to be plotted
#' @return a plot
#' @examples
#' NormalizeZscore(Trackingdata, point = "bodycentre")
#'
OverviewPlot <- function(t, point){
  if(!IsTrackingData(t)){
    stop("Input needs to be a single object of type Tracking")
  }
  if(!point %in% names(t$data)){
    warning("point does not exist in Trackingobject")
  }
  
  rel_h <- c(0.2,1.5)
  title <- ggdraw() + draw_label(paste("Overview file",t$filename, sep = " "))
  p1 <- PlotDensityPaths(t,point)
  if(!is.null(t$zones)){
    p1 <- AddZonesToPlots(p1, t$zones)
  }
  p1 <- p1[[1]] + scale_y_reverse()
  
  if(length(t$labels) > 0){
    p2 <- PlotLabels(t) +  theme(legend.position = "none")
    rel_h <- append(rel_h, 0.5)
  }
  if(!is.null(t$zones)){
    p3 <- PlotZoneVisits(t,points = point) +  theme(legend.position = "none")
    rel_h <- append(rel_h, 0.5)
  }
  
  if((length(t$labels) > 0) & !is.null(t$zones)){
    return(plot_grid(title,p1,p2,p3,rel_heights = rel_h, ncol = 1))
  }else if(length(t$labels) > 0){
    return(plot_grid(title,p1,p2,rel_heights = rel_h, ncol = 1))
  }else if(!is.null(t$zones)){
    return(plot_grid(title,p1,p3,rel_heights = rel_h, ncol = 1))
  }
  else{
    return(plot_grid(title,p1,rel_heights = rel_h, ncol = 1))
  }
}

#' mean Zscore normalization of a numeric matrix across columns
#' 
#' @param x a numeric matrix
#' @return a numeric matrix
#' @examples
#' NormalizeZscore(mymatrix)
#'
NormalizeZscore <- function(x){
  apply(x, 2, FUN = function(x){(x - mean(x)) / (sd(x))})
}

#' median Zscore normalization of a numeric matrix across columns
#' 
#' @param x a numeric matrix
#' @return a numeric matrix
#' @examples
#' NormalizeZscore_median(mymatrix)
#'
NormalizeZscore_median <- function(x){
  apply(x, 2, FUN = function(x){(x - median(x)) / (sd(x))})
}

#' integration of a numeric vector
#' 
#' @param x a numeric vector
#' @return a numeric vector
#' @examples
#' integratevector(myvector)
#'
integratevector <- function(x){
  if(length(x) < 2){
    stop("can  not integrate a vector of length < 2")
  }
  append(0, x[2:length(x)] - x[1:(length(x)-1)])
}

#' Boolean smoothing over an integration period
#' 
#' @param x a boolean vector
#' @param window a integration window length (in +- window entries)
#' @return a boolean vector
#' @examples
#' avgbool(booleanvector, window = 10)
#' 
avgbool <- function(x, window){
  res <- rep(0, length(x))
  for(i in 1:length(x)){
    res[i] <- mean(x[max(0,i-window):min(length(x), i + window)], na.rm = T) > 0.5
  }
  return(res)
}

#' Equalizes a training set so every group is equaly represented
#' 
#' @param x the training data
#' @param y the categorical labeling data (as produced by the library(keras) funciton to_categorical)
#' @return a list with 2 elements, adjusted x and adjusted y
#' @examples
#' EqualizeTrainingSet(x_train,y_train)
#' 
EqualizeTrainingSet <- function(x,y){
  N_obs <- NULL
  for(i in 1:ncol(y)){
    N_obs <- append(N_obs,sum(y[,i]))
  }
  
  keep <- rep(FALSE,nrow(y))
  
  for(i in 1:ncol(y)){
    keep[sample(which(y[,i] == 1))[1:min(N_obs)]] <- TRUE
  }
  out <- list()
  out$x <- x[keep,]
  out$y <- y[keep,]
  return(out)
}


#' Evaluates the performance of the classifier across one or mutlipe TrackingData objects by comparing a ground truth to any label
#' 
#' @param ts a list of objects of type TrackingData or a single object of type TrackingData
#' @param truth name of the label group that is considered the ground truth. by default "manual"
#' @param compare name of the label group that is considered the label to compare. by default "classifications"
#' @return a Report object with two sub-object, each a data.frame for the evaluation across or within files
#' @examples
#' EvaluateClassification(ts)
#' EvaluateClassification(ts, truth = "manual", compare = "classifications")
#' 
EvaluateClassification <- function(ts, truth = "manual", compare = "classifications"){
  if(IsTrackingData(ts)){
    x <- ts
    ts <- list()
    ts[[paste(x$filename)]] <- x
  }
  Report <- list()
  Report$files <- NULL
  for(i in names(ts)){
    if(is.null(ts[[i]]$labels)){
      stop(paste("file",i,"Does not have label data", sep = " "))
    }
    if(is.null(ts[[i]]$labels[[truth]])){
      stop(paste("file",i,"Does not have label type", truth, sep = " "))
    }
    if(is.null(ts[[i]]$labels[[compare]])){
      stop(paste("file",i,"Does not have label type", compare, sep = " "))
    }
    for(j in na.omit(unique(ts[[i]]$labels[[compare]]))){
      precision <- sum(ts[[i]]$labels[[compare]] == j & ts[[i]]$labels[[truth]] == j, na.rm = T) / sum(ts[[i]]$labels[[compare]] == j, na.rm = T)
      N_truth = sum(ts[[i]]$labels[[truth]] == j, na.rm = T)
      N_compare = sum(ts[[i]]$labels[[compare]] == j, na.rm = T)
      correct <- sum(ts[[i]]$labels[[compare]] == j & ts[[i]]$labels[[truth]] == j, na.rm = T)
      wrong <- sum(ts[[i]]$labels[[compare]] != j & ts[[i]]$labels[[truth]] == j, na.rm = T)
      recall <- sum(ts[[i]]$labels[[compare]] == j & ts[[i]]$labels[[truth]] == j, na.rm = T) / sum(ts[[i]]$labels[[truth]] == j, na.rm = T)
      Report$files <- rbind(Report$files,data.frame(file = i, 
                                                    label = j, 
                                                    accuracy = correct/(correct + wrong),
                                                    precision = correct / N_compare, 
                                                    recall = correct / N_truth,
                                                    correct = correct, 
                                                    wrong = wrong, 
                                                    N_truth = N_truth, 
                                                    N_compare = N_compare))
    }
    
  }
  
  Report$overall <- NULL
  for(i in unique(Report$files$label)){
    s <- apply(Report$files[Report$files$label == i,-c(1:5)],2,FUN = sum)
    Report$overall <- rbind(Report$overall, data.frame(label = i,
                                                       accuracy = s[1]/(s[1] +s[2]), 
                                                       precision =  s[1] / s[4], 
                                                       recall =   s[1] / s[3],
                                                       correct = s[1],
                                                       wrong = s[2],
                                                       N_truth = s[3],
                                                       N_compare = s[4]
    ))
  }
  return(Report)
}


#' Creates correlation plots between different labels over multiple files
#' 
#' @param ts a list of objects of type TrackingData
#' @param include a character vector that describes which labels to include. defaults to all
#' @param smooth a integer that describes over how many frames labels should be smoothed. defaults to NULL = no smoothing.
#' @param hclust a boolean that decides if plots should be ordered by hclust
#' @return a correlation plot
#' @examples
#' CorrelationPlotLabels(ts)
#' CorrelationPlotLabels(ts, include = c("manual.Unsupported.count","classification.Unsupported.count","manual.Supported.count","classification.Supported.count"), hclust = TRUE)
#' 
CorrelationPlotLabels <- function(ts, include = NULL, smooth = NULL, hclust = FALSE){
  compare <- list()
  for(i in ts){
    if(!is.null(smooth)){
      i <- SmoothLabels(i, smooth)
    }
    compare <- rbindlist(list(compare,LabelReport(i)),use.names = TRUE, fill = TRUE,idcol = F)
  }
  compare <-as.data.frame(compare)
  
  if(is.null(include)){
    include <- names(compare)
  }
  
  if(hclust){
    corrplot(cor(as.matrix(na_replace(compare[,include]))),title = "Correlation Plot", 
             method = "square", 
             outline = T, 
             addgrid.col = "darkgray", 
             order="hclust")
  }else{
  corrplot(cor(as.matrix(na_replace(compare[,include]))),
           method="color",
           type = "upper",
           addCoef.col = "black")
  }
}

#' Extracts labels from a long format data.frame
#' 
#' @param lab a long format data frame with labeling data
#' @param Experimenter a character or vector of character. the name of the experimenter(s) that should be extracted
#' @param type a character or vector of character. the name of type(s) that should be extracted
#' @param DLCFile a character or vector of character. the name of DLCFile(s) that should be extracted
#' @param ID a character or vector of character. the name of ID(s) that should be extracted
#' @return a data.frame extracted labels
#' @examples
#' ExtractLabels(lab, Experimenter = "Oliver", DLCFile = "FST_1.csv")
#' ExtractLabels(lab, Experimenter = "Oliver", type = c("Supported","Unsupported"))
#' 
ExtractLabels <- function(lab, Experimenter=NULL, type = NULL, DLCFile = NULL, ID = NULL){
  if(!is.null(Experimenter)){
    lab <- lab[lab$Experimenter %in% Experimenter,]
  }
  if(!is.null(type)){
    lab <- lab[lab$type %in% type,]
  }
  if(!is.null(DLCFile)){
    lab <- lab[lab$DLCFile %in% DLCFile,]
  }
  if(!is.null(ID)){
    lab <- lab[lab$ID %in% ID,]
  }
  return(lab)
}

#' Linearly scales a number of features by a set scaling factor
#' 
#' @param feat a numeric data.frame or matrix
#' @param select if set, a vector of column names that should be scaled. otherwise all will be scaled
#' @param factor vector of column names that should be scaled
#' @return a numeric data.frame or matrix
#' @examples
#' ScaleFeatures(feat, factor = 0.3)
#' ScaleFeatures(feat, select = c("feat1","feat3"), factor = 2)
#'
ScaleFeatures <- function(feat, select = NULL, factor){
  if(is.null(select)){
    feat <- names(feat)
  }
  feat[,select] <- feat[,select] * factor
  return(feat)
}

#' Combines multiple specified labels into a new, seperate label
#' 
#' @param ts a list of objects of type TrackingData or a single object of type TrackingData
#' @param which.lab the name of the label group that should be used, defaults to "unsupervised"
#' @param by a vector of label names that should be aggregated
#' @param name a string that contains the name of the newly created lable group
#' @return a list of objects of type TrackingData or a single object of type TrackingData
#' @examples
#' CombineLabels(ts, which.lab = "manual", by = c("Unsupported","Supported"), name = "AllRears")
#' CombineLabels(ts, which.lab = "unsupervised", by = list(c("cluster1","cluster2"),c("cluster5","cluster12","cluster6"),c("cluster8")), name = "CombinedClusters")
#' 
CombineLabels <- function(ts, which.lab = "unsupervised", by, name = "combined"){
  if(IsTrackingData(ts)){
    x <- ts
    ts <- list()
    ts[[paste(x$filename)]] <- x
  }
  
  for(i in names(ts)){
    if(!IsTrackingData(ts[[i]])){
      stop(paste("element",i,"is not of type TrackingData", sep = " "))
    }
    if(is.null(ts[[i]]$labels)){
      warning(paste("element",i,"has no labeling data", sep = " "))
      return(ts)
    }
    if(is.null(ts[[i]]$labels[[which.lab]])){
      warning(paste("element",i,"has no labeling data of type", which.lab, sep = " "))
      return(ts)
    }
    cdat <- ts[[i]]$labels[[which.lab]]
    out <- rep("None", length(ts[[i]]$labels[[which.lab]]))
    for(j in by){
      out[cdat %in% j] <- paste(j, collapse = ".")
    }
    ts[[i]]$labels[[name]] <- out
  }
  if(length(ts) > 1){
  return(ts)
  }
  return(ts[[1]])
}


#' Calculates the minimum distance of a point to a zone
#' 
#' @param t an object of type TrackingData
#' @param zone a name of zone (has to be present in t)
#' @param point a name of point (has to be present in t)
#' @return an array of distances of point p to edge of the zone at each frame
#' @examples
#' GetDistanceToZoneBorder(t = Tracking, zone = "arena", point = "bodycentre")
#' 
GetDistanceToZoneBorder <- function(t,zone,point){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(is.null(t$zones)){
    stop("TrackingData object has no zones")
  }
  if(!(zone %in% names(t$zones))){
    stop("invalid zone")
  }
  if(!(point %in% names(t$data))){
    stop("invalid point")
  }
  return(DistanceToPolygon(t$data[[point]][,c("x","y")], t$zones[[zone]][,c("x","y")]))
}

#' Calculates the minimum distance of a point p to any edge of a polygon pol
#' 
#' @param p a list() or data.frame() of a point. needs to contain p$x and p$y
#' @param pol a data.frame() of a polygon. needs to contain at least to columns, pol$x and pol$y
#' @return an object of type TrackingData 
#' @examples
#' DistanceToPolygon(p = data.frame(x = 0.5, y = 0.7), pol = data.frame(x = c(1,1,0,0), y = c(1,0,0,1)))
#' 
DistanceToPolygon <- function(p,pol){
  dist2d <- function(a,b,c){
    v1x <- b$x - c$x
    v1y <- b$y - c$y
    v2x <- a$x - b$x
    v2y <- a$y - b$y
    det <- v1x*v2y - v1y*v2x 
    
    d <- abs(det)/sqrt(v1x*v1x + v1y*v1y)
    return(d)
  } 
  
  pol <- pol[c(1:nrow(pol),1),]
  dist <- NULL
  
  for(i in 1:(nrow(pol) - 1)){
    dist <-  cbind(dist,dist2d(p, pol[i,], pol[i + 1,]))
  }
  return(apply(dist, 1, FUN=min))
}

#' Rotates an object of type Tracking data around a center of gravity defined by a number of points. rotation can be defined in degree, or alternatively a line of to points can be set parallel to the x-axis
#' 
#' @param t an object of type TrackingData
#' @param theta an angle in degree
#' @param center.of a vector of point names used to define a rotation center (mean value of x and y coordinats of the specified points)
#' @param set.straight.to.x a vector with the names of two points that define a line that should be set parallel to the x-axis
#' @return an object of type TrackingData 
#' @examples
#' RotateTrackingData(t = Tracking,theta = 45, center.of = c("tl","tr","bl","br"))
#' RotateTrackingData(t = Tracking,set.straight.to.x = c("br","tr"), center.of = c("tl","tr","bl","br"))
#' 
RotateTrackingData <- function(t, theta = 0, center.of = c("tl","tr","bl","br"), set.straight.to.x = NULL){
  theta <- theta * pi / 180
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }
  if(sum(center.of %in% names(t$data)) != length(center.of)){
    stop("invalid points selected for center.of")
  }
  if(!is.null(set.straight.to.x)){
    if(length(set.straight.to.x) != 2 | sum(set.straight.to.x %in% names(t$data)) != 2){
      stop("invalid set.straigth.to.x: needs to be a vector with valid point names of exactly 2 points that define a line")
    }
    ax <- t$median.data[set.straight.to.x[2],"x"] - t$median.data[set.straight.to.x[1],"x"] 
    ay <- t$median.data[set.straight.to.x[2],"y"] - t$median.data[set.straight.to.x[1],"y"] 
    print(ax)
    print(ay)
    theta = -acos(abs(ax) / sqrt(ax*ax + ay*ay))
  }
  center_x <- mean(t$median.data[center.of,"x"])
  center_y <- mean(t$median.data[center.of,"y"])
  
  if(!is.null(t$zones)){
    for(i in names(t$zones)){
      old_x <- t$zones[[i]]$x - center_x
      old_y <- t$zones[[i]]$y - center_y
      t$zones[[i]]$x <- old_x * cos(theta) - old_y * sin(theta) + center_x
      t$zones[[i]]$y <- old_x * sin(theta) + old_y * cos(theta) + center_y
    }
  }
  for(i in names(t$data)){
    old_x <- t$data[[i]]$x - center_x
    old_y <- t$data[[i]]$y - center_y
    t$data[[i]]$x <- old_x * cos(theta) - old_y * sin(theta) + center_x
    t$data[[i]]$y <- old_x * sin(theta) + old_y * cos(theta) + center_y
  }
  
  old_x <- t$median.data$x - center_x
  old_y <- t$median.data$y - center_y
  t$median.data$x <- old_x * cos(theta) - old_y * sin(theta) + center_x
  t$median.data$y <- old_x * sin(theta) + old_y * cos(theta) + center_y
  return(t)
}

CreateExampleVideos <- function(t, labels.group, label.type, video, n = NULL, random = TRUE, folder = "Examples", name = "Example", min.length = 0, lag = 0){
  if(!IsTrackingData(t)){
    stop("Object is not of type TrackingData")
  }  
  if(is.null(t$labels)){
    stop("Object has no labeling data")
  }
  if(is.null(t$labels[[labels.group]])){
    stop("Specified labels.group invalid")
  }
  
  dir.create(folder)
  
  examples <- data.frame()
  active = FALSE
  start = NULL
  end = NULL
  for(i in 1:length(t$labels[[labels.group]])){
    frame <- t$frames[i]
    if(is.na(t$labels[[labels.group]][i])){
      
    }
    else if(!active){
      if(t$labels[[labels.group]][i] == label.type){
        active = TRUE
        start = frame
      }
    }
    else{
      if(t$labels[[labels.group]][i] != label.type){
        active = FALSE
        end = frame
        examples <- rbind(examples,data.frame(start = start, end = end, name = paste(name,frame / t$fps,sep ="_")))
      }
    }
  }
  
  examples$valid <- (examples$end - examples$start) / t$fps >= min.length
  examples$start <- examples$start - lag
  examples$end <- examples$end + lag
  examples <- examples[examples$valid,]
  n <- min(n, nrow(examples))
  if(random){
    examples <- examples[sample(1:nrow(examples)),]
  }
  
  #ffmpeg -ss 300 -i vid\\OFT_1.mp4 -t 20 -an vid\\TEST2.mp4
  for(i in 1:n){
    out <- paste(folder, "\\",examples[i,"name"],".mp4", sep = "")
    command <- paste("ffmpeg -ss", examples[i,"start"] / t$fps, 
                     "-t", (examples[i,"end"] - examples[i,"start"])/ t$fps,
                     "-i", video, 
                     out,
                     sep = " ")
    system(command)
  }
  return(NULL)
}

TwoGroupComparisonReport <- function(ts, by, FDR.cutoff = 0.5){
  if((length(unique(by)) != 2) | (length(by) != length(ts))){
    stop("Incorect grouping variable by. needs to be same length as trackingdata list and have exactly 2 unique values")
  }

  Report <- MultiFileReport(ts)
  Report <- na.replace(Report,fill = 0)
  res <- NULL
  for(i in names(Report)[-1]){
    a <- Report[by == unique(by)[1],i]
    b <- Report[by == unique(by)[2],i]
    res <- rbind(res,data.frame(name = i, p = t.test(nafill(a, fill = 0),nafill(b, fill = 0))$p.value))
  }
  res$FDR <- p.adjust(res$p, method = "BY")
  
  ForPlot <- data.frame(Report, Group = by)
  row.names(res) <- res$name
  ForPlot <- gather(ForPlot,key = metric,value = value,2:(ncol(ForPlot)-1))
  significant <- res[ForPlot$metric,"FDR"] < FDR.cutoff

  p1 <- ggplot(data = ForPlot, aes(Group,value, color = significant)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_point(position = position_jitter(height = 0, width = 0.2)) + 
    facet_wrap(~metric, scales = "free") + 
    scale_color_manual(values = c("black","red")) + 
    theme_bw()
  p2 <- ggplot(data = ForPlot[significant,], aes(Group,value, color = Group)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_point(position = position_jitter(height = 0, width = 0.2)) + 
    facet_wrap(~metric, scales = "free") + 
    scale_color_manual(values = c("black","red")) + 
    theme_bw()
  
  out <- list()
  out$Results <- res
  out$MatrixData <- Report
  out$PlotData <- ForPlot
  out$PlotsAll <- p1
  out$PlotsSignificant <- p2
  
  return(out)
}

ANOVATwoWayReport <- function(ts, design, FDR.cutoff = 0.5){
  
  Report <- MultiFileReport(ts)
  
  coeffa <- NULL
  coeffb <- NULL
  coeffc <- NULL
  
  Report <- na_replace(Report)
  for(i in names(Report)[-1]){
    AOVdata <- data.frame(a = Report[,i], b = design[,1], c = design[,2])
    AOV <- summary(aov(AOVdata$a ~ AOVdata$b * AOVdata$c))
    
    coeffa <- rbind(coeffa,data.frame(i,AOV[[1]][1,5]))
    coeffb <- rbind(coeffb,data.frame(i,AOV[[1]][2,5]))
    coeffc <- rbind(coeffc,data.frame(i,AOV[[1]][3,5]))
  }
  coeffa$FDR <- p.adjust(coeffa[,2], method = "BY")
  coeffb$FDR <- p.adjust(coeffb[,2], method = "BY")
  coeffc$FDR <- p.adjust(coeffc[,2], method = "BY")
  
  Results <- cbind(coeffa, coeffb[,-1],coeffc[,-1])
  
  names(Results) <- c("name", 
                      paste(names(design)[1],"p.value", sep = "."), 
                      paste(names(design)[1],"FDR", sep = "."),
                      paste(names(design)[2],"p.value", sep = "."),
                      paste(names(design)[2],"FDR", sep = "."),
                      paste("interaction","p.value", sep = "."),
                      paste("interaction","FDR", sep = "."))
  rownames(Results) <- Results$name
   
  ForPlot <- data.frame(Report, design)
  ForPlot <- gather(ForPlot,key = metric,value = value,2:(ncol(ForPlot)-2))

  significant <- (Results[ForPlot$metric,3] < FDR.cutoff | Results[ForPlot$metric,5] < FDR.cutoff | Results[ForPlot$metric,7] < FDR.cutoff)
  
  p1 <- ggplot(data = ForPlot, aes(ForPlot[,3],value, color = ForPlot[,2])) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_point(position = position_jitterdodge()) + 
    facet_wrap(~metric, scales = "free") + 
    scale_color_manual(values = c("black","red")) + 
    theme_bw()
  p2 <- ggplot(data = ForPlot[significant,], aes(ForPlot[significant,3],value, color = ForPlot[significant,2])) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_point(position = position_jitterdodge()) + 
    facet_wrap(~metric, scales = "free") + 
    scale_color_manual(values = c("black","red")) + 
    theme_bw()
  
  
  out <- list()
  out$Results <- Results
  out$sig <- significant
  out$MatrixData <- Report
  out$PlotData <- ForPlot
  out$PlotsAll <- p1
  out$PlotsSignificant <- p2
  
  return(out)
}

LMReport <- function(ts, by, FDR.cutoff = 0.5){
  
  Report <- MultiFileReport(ts)
  Report <- na.replace(Report,fill = 0)
  
  res <- NULL
  for(i in names(Report)[-1]){
    mod <- lm(Report[,i]~by)
    mod <- summary(mod)

    res <- rbind(res,data.frame(name = i, p = mod$coefficients[2,4]))
  }
  res$FDR <- p.adjust(res$p, method = "BY")
  
  ForPlot <- data.frame(Report, Group = by)
  row.names(res) <- res$name
  ForPlot <- gather(ForPlot,key = metric,value = value,2:(ncol(ForPlot)-1))
  significant <- res[ForPlot$metric,"FDR"] < FDR.cutoff
  
  p1 <- ggplot(data = ForPlot, aes(Group,value, color = significant)) + 
    geom_point() + 
    facet_wrap(~metric, scales = "free") + 
    scale_color_manual(values = c("black","red")) + 
    theme_bw()
  p2 <- ggplot(data = ForPlot[significant,], aes(Group,value)) + 
    geom_point() + 
    facet_wrap(~metric, scales = "free") + 
    theme_bw()
  
  out <- list()
  out$Results <- res
  out$MatrixData <- Report
  out$PlotData <- ForPlot
  out$PlotsAll <- p1
  out$PlotsSignificant <- p2
  
  return(out)
}

multitsne <- function(ts, color.by = "unsupervised", perplex = 500, Z_score_Normalize = FALSE, dimensions = NULL){
  allx <- list()
  ally <- list()
  id <- NULL
  
  for(j in names(ts)){
    allx[[j]] <- ts[[j]]$train_x
    ally[[j]] <- na.omit(ts[[j]]$labels[[color.by]])
    id <- append(id, rep(paste(j),nrow(ts[[j]]$train_x)))
    ts[[j]]$train_x <- NULL
  }
  
  allx <- do.call(rbind, allx)
  ally <- do.call(rbind, ally)
  
  if(Z_score_Normalize){
    for(i in 1:ncol(allx)){
      x <- allx[,i]
      allx[,i] <- (x - mean(x))/(sd(x))
    }
  }
  
  if(!is.null(dimensions)){
    allx <- svd(allx, nu = dimensions, nv = 0)$u
  }
  
  tsneres <- tsne(t(allx),labels = ally,dotsize = 0.5, perplex = perplex, seed = 123)
  
  for(j in names(ts)){
    ts[[j]]$tsnecoords <- data.frame(tsneres$data[id == j,], group = ally[id == j])
  }
  return(ts)
}


Tracking2DEmbedding <- function(ts, color.by = "unsupervised", perplex = 500, Z_score_Normalize = FALSE, method = "umap", n_points = NULL, seed = 123, dotsize = 0.5){
  library(M3C)
  allx <- list()
  ally <- NULL
  id <- NULL
  
  for(j in names(ts)){
    allx[[j]] <- ts[[j]]$train_x
    ally <- append(ally,na.omit(ts[[j]]$labels[[color.by]]))
  }
  
  allx <- do.call(rbind, allx)
  
  if(!is.null(n_points) & (n_points < nrow(allx))){
    idx <- sample(1:nrow(allx))[1:n_points]
    allx <- allx[idx,]
    ally <- ally[idx]
  }
  
  
  if(Z_score_Normalize){
    for(i in 1:ncol(allx)){
      x <- allx[,i]
      allx[,i] <- (x - mean(x))/(sd(x))
    }
  }
  
  if(method == "tsne"){
    M3C::tsne(t(allx),labels = ally,dotsize = dotsize, perplex = perplex, seed = seed)
  }else
  {
    M3C::umap(t(allx),labels = ally,dotsize = dotsize, seed = seed)
  }
}



FeatureImportance_Permutation <- function(MLData, model, depth, type, shuffle = FALSE){
  if(sum(type == MLData$parameters$Feature_names) < 1){
    warning("Invalid type, needs to be contained in MLData$paramteres$Feature_names")
  }
  if(depth > nrow(MLData$train_x)){
    warning("Less data in MLData than selected depth")
  }
  
  if(shuffle){
    index <- sample(1:nrow(MLData$train_x),depth)
  }else{
    index <- 1:depth
  }
  select <- which(MLData$parameters$Feature_names == type)
  
  pred_wrapper <- function(object, newdata) {
    predict(object, x = as.matrix(newdata))[,select] %>%
      as.vector()
  }
  
  train_x <- MLData$train_x[index,]
  train_y <- MLData$train_y[index,]
  
  p1 <- vip(
    object = model,                     # fitted model
    method = "permute",                 # permutation-based VI scores
    num_features = ncol(train_x),       # default only plots top 10 features
    pred_wrapper = pred_wrapper,            # user-defined prediction function
    train = as.data.frame(train_x) ,    # training data
    target = train_y[,select],                   # response values used for training
    metric = "rsquared",                # evaluation metric
    progress = "text"                 # request a text-based progress bar
  )
  dat <- p1$data
  dat$Var <- str_split_fixed(dat$Variable,",",2)[,1]
  dat$Index <- as.numeric(str_split_fixed(dat$Variable,",",2)[,2])
  
  print(ggplot(dat, aes(fill=as.factor(Index), y=Importance , x=Var)) + 
          geom_bar(position="stack", stat="identity") + theme_bw())
  
  return(dat)
}

FeatureImportance_1DALE <- function(MLData, model, depth, type, feature, shuffle = FALSE){
  if(sum(type == MLData$parameters$Feature_names) < 1){
    warning("Invalid type, needs to be contained in MLData$paramteres$Feature_names")
  }
  if(depth > nrow(MLData$train_x)){
    warning("Less data in MLData than selected depth")
  }
  
  if(shuffle){
    index <- sample(1:nrow(MLData$train_x),depth)
  }else{
    index <- 1:depth
  }
  select <- which(MLData$parameters$Feature_names == type)
  
  pred_wrapper_ALE <- function(X.model, newdata) {
    predict(X.model, x = as.matrix(newdata))[,select] %>%
      as.vector()
  }
  
  train_x <- MLData$train_x[index,]
  dat <- NULL
  for(i in -MLData$parameters$integration_period:MLData$parameters$integration_period){
    a <- ALEPlot(X = as.data.frame(train_x), X.model = model,pred.fun = pred_wrapper_ALE, J = paste(feature,i, sep = ","))
    dat <- rbind(dat, data.frame(index = i, value = a$x.values, prediction = a$f.values))
  }
  
  return(dat)
}

PlotFeatureImportance_1DALE <- function(data, filename, lower = NULL, upper = NULL, xaxis = NULL, scale = NULL){
  if(!is.null(lower)){
    data <- data[data$value > lower,]
  }
  if(!is.null(upper)){
    data <- data[data$value < upper,]
  }
  if(!is.null(scale)){
    data$value <- data$value * scale
  }
  
  pdf(filename)
    p1 <- ggplot(data, aes(value, prediction, color = as.factor(index))) + 
      geom_path() + 
      geom_hline(yintercept = 0) + 
      ylab("Effect on prediction") + 
      theme_bw() +
      labs(color = "Frame of Sequence")
    if(!is.null(xaxis)){
      p1 <- p1 + xlab(xaxis)
    }
    print(p1)
    p2 <- ggplot(data, aes(value, prediction, color = as.factor(index))) + 
      geom_path() + 
      geom_hline(yintercept = 0) + 
      ylab("Effect on prediction") + 
      theme_bw() +
      labs(color = "Frame of Sequence") +
      facet_wrap(~as.factor(index))
    if(!is.null(xaxis)){
      p2 <- p2 + xlab(xaxis) 
    }
    print(p2)
  dev.off()
  return(NULL)
}

prettyConfused<-function(Actual,Predict,colors=c("white","red4","dodgerblue3"),text.scl=5){
  actual = as.data.frame(table(Actual))
  names(actual) = c("Actual","ActualFreq")
  
  #build confusion matrix
  confusion = as.data.frame(table(Actual, Predict))
  names(confusion) = c("Actual","Predicted","Freq")
  
  #calculate percentage of test cases based on actual frequency
  
  confusion = merge(confusion, actual, by=c('Actual','Actual'))
  confusion$Percent = confusion$Freq/confusion$ActualFreq*100
  confusion$ColorScale<-confusion$Percent*-1
  confusion[which(confusion$Actual==confusion$Predicted),]$ColorScale<-confusion[which(confusion$Actual==confusion$Predicted),]$ColorScale*-1
  confusion$Label<-paste(round(confusion$Percent,0),"%, n=",confusion$Freq,sep="")
  tile <- ggplot() +
    geom_tile(aes(x=Actual, y=Predicted,fill=ColorScale),data=confusion, color="black",size=0.1) +
    labs(x="Actual",y="Predicted")
  
  tile = tile +
    geom_text(aes(x=Actual,y=Predicted, label=Label),data=confusion, size=text.scl, colour="black") +
    scale_fill_gradient2(low=colors[2],high=colors[3],mid=colors[1],midpoint = 0,guide='none')
}

CalculateConfusionMatrix <- function(ts, truth = "manual", compare = "classifications",text.scl = 3){
  tr <- NULL
  comp <- NULL
  for(i in ts){
    tr <- append(tr,i$labels[[truth]])
    comp <- append(comp,i$labels[[compare]])
  }
  out <- list()
  out$plot <- prettyConfused(tr,comp,text.scl = text.scl)
  out$Results <- confusionMatrix(as.factor(tr),as.factor(comp))
  return(out)
}

CalculateTransitionMatrix <- function(ts, label = "unsupervised", smooth = NULL,text.scl = 3){
  lab <- NULL
  out <- list()
  for(i in ts){
    lab <- append(lab,i$labels[[label]])
  }
  if(!is.null(smooth)){
    lab <- SmoothLabel(lab,smooth)
  }
  out$labels <- lab
  out$transvector <- na.omit(lab[lab[2:length(lab)] != lab[1:(length(lab)-1)]])
  out$transmatrix <- table(out$transvector[1:(length(out$transvector)-1)],out$transvector[2:length(out$transvector)])
  out$plot <- prettyConfused(out$transvector[1:(length(out$transvector)-1)],out$transvector[2:length(out$transvector)],text.scl = text.scl)
  out$plot <- out$plot + xlab("From") + ylab("To")
  return(out)
}

CalculateTransitionMatrices <- function(ts, label = "unsupervised", smooth = NULL){
  lab <- NULL
  for(i in ts){
    lab <- append(lab,i$labels[[label]])
  }
  lev <- unique(na.omit(lab))
  
  out <- list()
  out$levels <- lev
  out$transitionmatrix <- list()
  out$norm_transitionmatrix <- list()
  
  #create ordered and complete transition matrix for each files
  for(i in names(ts)){
    lab <- ts[[i]]$labels[[label]]
    if(!is.null(smooth)){
      lab <- SmoothLabel(lab,smooth)
    }
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
    out$transitionmatrix[[paste(i)]] <- transmatrix
    out$norm_transitionmatrix[[paste(i)]] <- t(apply(transmatrix, MARGIN = 1, FUN = function(x){if(sum(x) > 0){x / sum(x)}else{x}}))
  }
  return(out)
}
 
CalculateOverlapMatrix <- function(ts, from = "unsupervised",to = "classifications", smooth = NULL,text.scl = 3){
  A <- NULL
  B <- NULL
  out <- list()
  for(i in ts){
    A <- append(A,i$labels[[from]])
    B <- append(B,i$labels[[to]])
  }
  if(!is.null(smooth)){
    A <- SmoothLabel(A,smooth)
    B <- SmoothLabel(B,smooth)
  }
  out$from <- A
  out$to <- B
  out$transmatrix <- table(from = A, to = B)
  out$transmatrix_norm <- t(apply(out$transmatrix,MARGIN = 1, FUN = function(x){x / sum(x)}))
  return(out)
}
  

EvaluateGroupedTransitions <- function(ts, groups, label = "unsupervised", smooth = NULL, n_bootstraps){
  if(length(unique(groups)) != 2){
    warning("grouping variable does have more or less unique values than 2!")
    return(NULL)
  }
  Results <- CalculateTransitionMatrices(ts, label = label, smooth = smooth)
  
  A <- Results$transitionmatrix[groups == unique(groups)[1]]
  B  <- Results$transitionmatrix[groups == unique(groups)[2]]
  
  a <- do.call(cbind, A)
  a <- array(a, dim=c(dim(A[[1]]), length(A)))
  a <- apply(a, c(1, 2), mean, na.rm = TRUE)
  
  b <- do.call(cbind, B)
  b <- array(b, dim=c(dim(B[[1]]), length(B)))
  b <- apply(b, c(1, 2), mean, na.rm = TRUE)
  
  distance <- sum(abs(a-b))
  
  bootstraps <- NULL
  for(i in 1:n_bootstraps){
    grp <- sample(groups)
    A <- Results$transitionmatrix[grp == unique(groups)[1]]
    B  <- Results$transitionmatrix[grp == unique(groups)[2]]
    
    a <- do.call(cbind, A)
    a <- array(a, dim=c(dim(A[[1]]), length(A)))
    a <- apply(a, c(1, 2), mean, na.rm = TRUE)
    
    b <- do.call(cbind, B)
    b <- array(b, dim=c(dim(B[[1]]), length(B)))
    b <- apply(b, c(1, 2), mean, na.rm = TRUE)
    
    bootstraps <- append(bootstraps,sum(abs(a-b)))
  }
  
  
  transitions <- NULL
  for(i in Results$levels){
    for(j in Results$levels){
      A <- lapply(Results$transitionmatrix[groups == unique(groups)[1]],FUN = function(x){x[i,j]})
      B <- lapply(Results$transitionmatrix[groups == unique(groups)[2]],FUN = function(x){x[i,j]})
      test <- t.test(unlist(A),unlist(B))
      transitions <- rbind(transitions, data.frame(from = j, to = i, p.value = test$p.value, A = mean(unlist(A)), B = mean(unlist(B)) ,FC = mean(unlist(A)) / mean(unlist(B))))
    }
  }
  names(transitions)[4:5] <- c(paste(unique(groups)[1]),paste(unique(groups)[2]))
  transitions$fdr <- p.adjust(transitions$p.value)
  
  transitions_norm <- NULL
  for(i in Results$levels){
    for(j in Results$levels){
      A <- lapply(Results$norm_transitionmatrix[groups == unique(groups)[1]],FUN = function(x){x[i,j]})
      B <- lapply(Results$norm_transitionmatrix[groups == unique(groups)[2]],FUN = function(x){x[i,j]})
      test <- t.test(unlist(A),unlist(B))
      transitions_norm <- rbind(transitions_norm, data.frame(from = i, to = j, p.value = test$p.value, A = mean(unlist(A)), B = mean(unlist(B)) ,FC = mean(unlist(A)) / mean(unlist(B))))
    }
  }
  names(transitions_norm)[4:5] <- c(paste(unique(groups)[1]),paste(unique(groups)[2]))
  transitions_norm$fdr <- p.adjust(transitions_norm$p.value)
  
  out <- list()
  out$bootstrap <- ggplot(data.frame(distance = bootstraps), aes(distance)) + geom_histogram(color = "black", fill = "grey50") + geom_vline(xintercept = distance, color = "red") + theme_bw()
  out$distance <- distance
  out$bootstraps <- data.frame(distance = bootstraps)
  out$transitions <- transitions
  out$transitions_norm <- transitions_norm
  out$transmatrices <- Results$transitionmatrix
  out$transmatrices_norm <- Results$norm_transitionmatrix
  out$groups <- groups
  out$n_sigma <- (distance - mean(bootstraps)) / sd(bootstraps)
  
  return(out)
}

EvaluateLMTransitions <- function(ts, group, label = "unsupervised", smooth = NULL, n_bootstraps){

  Results <- CalculateTransitionMatrices(ts, label = label, smooth = smooth)
  
  transitions <- NULL
  for(i in Results$levels){
    for(j in Results$levels){
      A <- lapply(Results$transitionmatrix,FUN = function(x){x[i,j]})
      test <- lm(unlist(A)~group)
      test <- summary(test)
      transitions <- rbind(transitions, data.frame(from = j, to = i, p.value =test$coefficient[2,4], slope = test$coefficient[2,1]))
    }
  }
  transitions$fdr <- p.adjust(transitions$p.value)
  

  out <- list()
  out$transitions <- transitions
  out$transmatrices <- Results$transitionmatrix
  out$groups <- group
  
  return(out)
}
