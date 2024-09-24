library(tidyverse)
library(plotly)

#For analyzing output of DLC models to calculate (x,y) animal positions 
#Can use either filtered or unfiltered analysis from DLC but prefers unfiltered data for statistics

#==Define Constants==
#DLC annotation confidence threshold
threshold = 0.6
fps = 10
verbosePlotting = FALSE
removeOutliers = TRUE
#set this to be 1/2 the resolution of the video in whichever dimension is smallest
#(i.e for a recording that is 640x480, the threshold should be ~240)
outlierThreshold = 200
binning = FALSE


#==Functions==
#calculates euclidean distance of two (x,y) coordinates
euc.dist <- function(x1, x2) {
  sqrt(sum((x1 - x2) ^ 2)) 
}

process.set <- function(dfData, animalID) {
  #calculate euclidean distance between each point, getting the distance and acceleration
  dfDist <- data.frame("dist" = numeric(), "acc" = numeric())
  tmp = data.frame(0, 0)
  names(tmp) = c("dist", "acc")
  dfDist <- rbind(dfDist, tmp)
  
  outliers = 0
  excluded = 0
  distTraveled = 0
  for(a in 2:nrow(dfData)) {
    #if it falls below threshold, mark as NA for distance. Cant compute distance if any one is NA
    if (dfData$likelihood[a] < threshold) {
      dfData$x[a] = NA
      dfData$y[a] = NA
      excluded <- excluded + 1
    }
    if ((dfData$likelihood[a] < threshold) || (dfData$likelihood[a-1] < threshold)) {
      tmp = data.frame(NA, NA)
      names(tmp) <- c( "dist", "acc")
      dfDist = rbind(dfDist, tmp)
    } else {
      pt1 = data.frame(x = dfData$x[a-1], y = dfData$y[a-1])
      pt2 = data.frame(x = dfData$x[a], y = dfData$y[a])
      dist = euc.dist(pt1, pt2)
      if ((dist >= outlierThreshold) && (!is.na(dist))) {
        if (removeOutliers == TRUE) {
          #if true, a is an outlier and should be removed
          dfData$x[a] = NA
          dfData$y[a] = NA
          dist <- NA
          acc <- NA
        }
        outliers <- outliers + 1
      } else {
        acc <- dist - dfDist$dist[a-1]
        outliers <- outliers + 1
      }
      tmp = data.frame(dist, acc)
      names(tmp) <- c("dist", "acc")
      dfDist = rbind(dfDist, tmp)
      if (!is.na(dist)) {
        distTraveled <- distTraveled + dist
      }
    }
  }
  
  #calculate statistics
  part_percent <- (excluded / nrow(dfData)) * 100
  outlier_percent <- (outliers / nrow(dfData)) * 100
  dfStats = c(part_percent, outliers, outlier_percent, distTraveled)
  ret = c(dfData, dfDist, dfStats)
  return(ret)
}

process.bin <- function(dfData, interval) {
  df <- data.frame("time" = numeric(), "cum_dist" = numeric(), "avg_acc" = numeric())
  it = 0
  time = 1
  sumAcc = 0
  sumDist = 0
  for (a in 1:nrow(dfData)) {
    if (!is.na(dfData$acc[a])) {
      sumAcc <- sumAcc + dfData$acc[a]
    } 
    if (!is.na(dfData$dist[a])) {
      sumDist <- sumDist + dfData$dist[a]
    }
    if (it == interval || a+1 >= nrow(dfData)) {
      avgAcc = sumAcc / it
      tmp <- c(time, sumDist, avgAcc)
      names(tmp) <- c("time", "cum_dist", "avg_acc")
      dfDist = rbind(df, tmp)
      time <- time + 1
      it = 0
      sumAcc = 0
      sumDist = 0
    } else {
      it <- it + 1
    }
  }
  return(df)
}

#==Main==
print("Starting Analysis...")
print("Select files to analyze")
files <- choose.files()
print("Select save location")
saveDir <- choose.dir(default= "", caption = "Select folder to save graphs in")

#loop through all selected files and perform batched analysis
#expects the passed DLC csv files to be formatted with filenames separated by underscores 
for (a in 1:length(files)) {
    tmp <- str_split_fixed(basename(files[a]), "_", 4)
    animalID <- tmp[3]
    print(paste("Processing DLC Data for Animal: ", animalID))
    
    
    #==cleanup original DLC data for whole body mouse tracking==
    dataOrig <- read.csv(files[a], skip = 1, header = TRUE)
    colNames <- names(dataOrig)
    tmp <- c()
    for (x in 1:ncol(dataOrig)) {
      tmp <- c(tmp, paste(colNames[x], dataOrig[1,x], sep = "_"))
    }
    colnames(dataOrig) <- tmp
    dataOrig <- dataOrig[-1,]
    dataOrig <- dataOrig %>% mutate_all(.funs = readr::parse_number)
    
    #get names of parts: loop through data, separate out column names
    #skip first column which is the frame number
    dataProcessed <- data.frame("frame"= dataOrig$bodyparts_coords)
    dataStats = data.frame("part" = character(), "percent_labeled" = numeric(), "num_outlers" = numeric(), "percent_outlers" = numeric(), "dist_traveled" = numeric())
    dfBinned <- data.frame("time" = numeric())
    
    #take x,y, and likelihood values for each part and analyze
    for (b in seq(from=2, to=ncol(dataOrig), by=3))  {
      part = colNames[b]
      print(paste("Analyzing part:", part, sep= " "))
      part_x = paste(part, "_x", sep = "")
      part_y = paste(part, "_y", sep= "")
      part_dist = paste(part, "_dist", sep="")
      part_acc = paste(part, "_acc", sep="")
      part_cumdist = paste(part, "_cumdist", sep="")
      part_avgacc = paste(part, "_avgacc", sep="")
      
      tmp = data.frame("x"= dataOrig[[b]], "y" = dataOrig[[b+1]], "likelihood" = dataOrig[[b+2]])
      ret = process.set(tmp, animalID)
      dataProcessed[, part_x] <- ret[[1]]
      dataProcessed[, part_y] <- ret[[2]]
      dataProcessed[, part_dist] <- ret[[4]]
      dataProcessed[, part_acc] <- ret[[5]]
      tmp2 = data.frame(part, ret[[6]], ret[[7]], ret[[8]], ret[[9]])
      names(tmp2) = c("part", "percent_labeled", "num_outlers", "percent_outlers", "dist_traveled")
      dataStats = rbind(dataStats, tmp2)
      
      #bin data into seconds 
      if (binning == TRUE) {
        tmp3 = data.frame(dataProcessed[, part_dist], dataProcessed[, part_acc])
        names(tmp3) = c("dist", "acc")
        ret = process.bin(tmp3, fps)
        if (nrow(dfBinned) < 1) {
          dfBinned[, time] <- ret[[1]]
          dfBinned[, part_cumdist] <- ret[[2]]
          dfBinned[, part_avgacc] <- ret[[3]]
        } else {
          dfBinned[, part_cumdist] <- ret[[2]]
          dfBinned[, part_avgacc] <- ret[[3]]
        }
      }
      
      strName = paste(animalID, part, sep="-")
      print(paste("Plotting", part, "Data...", sep=" "))
      if (verbosePlotting == TRUE) {
        #(x,y) part plot
        p1 <- ggplot(dataProcessed, aes(x=!!sym(part_x), y=!!sym(part_y))) + geom_point(color="darkblue") + scale_y_reverse() + 
          ggtitle(strName)
        print(p1)
        #plot distance
        strName= paste("Distance(Px): ", animalID, " - ", part, sep="")
        p2 <- ggplot(dataProcessed, aes(x=frame, y=!!sym(part_dist))) + geom_line(color="darkblue", alpha=0.9, size = 1) + 
          ggtitle(strName) + 
          theme_light()
        print(p2)
      }
      #plot acceleration
      strName= paste("Accelleration(Px): ", animalID, " - ", part, sep="")
      p3 <- ggplot(dataProcessed, aes(x=frame, y=!!sym(part_acc))) + geom_line(color="darkblue", alpha=0.9, size = 1) + 
        ggtitle(strName) + 
        theme_light()
      print(p3)
    }
    #clear unnecessary data
    rm(p1)
    rm(p2)
    rm(p3)
    rm(ret)
    rm(tmp)
    rm(tmp2)
    
    #save all plots
    plotsPath <- list.files(tempdir(), pattern="rs-graphics", full.names=TRUE)
    plotspngPath <- list.files(plotsPath, pattern=".png", full.names=TRUE)
    file.copy(from=plotspngPath, to=saveDir)
}

