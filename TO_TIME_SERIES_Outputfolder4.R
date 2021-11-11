library(rstudioapi) #foldering
library(tidyr)      #data wrangling
library(signal)     #signal processing
library(pracma)     #signal processing
library(kza)        #smoothing


############################
#Author Script: Wim Pouw (wimpouw@gmail.com)

############################

#folders and datafiles
basefolder <- dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path))) #get path of current R doc
MTfolder3 <- paste0(basefolder, "/3_SYNCED_KINECT_TIMESTAMPED/")
MTS <- list.files(MTfolder3, pattern = "*.csv")
wrifolder4 <- paste0(basefolder, "/4_PROCESSED_FOR_ANALAYSIS/")
wrifolder4_1 <- paste0(wrifolder4, "4_1_AUTO_ANNOTATIONS/")
#############################FUNCTIONS

# we are using a Kolmogorov-Zurbenko filter with an order of 2 and a span of 5
smooth.it <- function(x) 
{x <- kz(x, 2, 3)}

#this function loads into a time series annotations dataframes with columns (begintime, endtime, annotation1, annotation2, ...)
#the first argument is the timeseries vector of your timeseries, the second argument is your annotation dataframe
load.in.event <- function(time_ms_rec, g_d)
{
  output <- character(length = length(time_ms_rec))
  output <- NA
  for(i in g_d[,1])
  {
    output <- ifelse((time_ms_rec >= g_d[,1][g_d[,1] == i] & time_ms_rec <= g_d[,2][g_d[,1] == i]), as.character(g_d[,3][g_d[,1]==i]), output)
  }
  return(output)
}

#THis function allows you to make event identifiers, so that each run of observations in a time series
# gets an id.
make.event.identifiers <- function(runsval, originaltime)
{ #MAKE SPEECH EVENTS
  rle_x <- rle(runsval)
  values = rle_x$values
  endtime <- cumsum(rle_x$lengths)
  begintime <- c(1, diff(endtime))
  begintime <- (endtime - begintime)+1
  runs <- data.frame(originaltime[begintime], originaltime[endtime], values)
  runs <- subset(runs, values == 1)
  runs <- subset(runs, select =-c(values))
  runs$s_events <- seq.int(nrow(runs))
  events <<- load.in.event(originaltime, runs)
}


#detect gesture function 
#this is the main auto gesture detection funciton
#based on speed (cm/s) treshold, a mindistformerge treshold stating when gestures are merged (millisecond gap),
#a minum duration (in ms) before a gesture is kept, and a speed vector (e.g., left hand speed), and the
#the time vector
#trials are also given as input as we only want codings within  relevant trials
get.autogesture.annotations <- function(verticalmov, speed, gesturespace, minduration, mindistformerge, speedvector, timevector, trial)
{
  ges <- ifelse(speedvector > speed & timevector > 0 & !is.na(trial) & verticalmov > (mean(verticalmov)-(gesturespace*sd(verticalmov))), 1, 0) #check for gestures/movements higher then speed and mark them in the timeseries
  
  if(1 %in% ges) #check if there are any movements higher then speedvector is detected, and if so, give it an id
  {ges <- make.event.identifiers( ges, timevector)}
  
  #auto coding
  begintime <- endtime <- annotation <- vector() #initialize annotation variables
  for(i in unique(ges[!is.na(ges)]))  #save the begin and endtimes of the observed gestures
  {
    begintime   <- c(begintime, min(timevector[ges==i], na.rm = TRUE))
    endtime     <- c(endtime, max(timevector[ges==i], na.rm = TRUE))
    annotation  <- c(annotation, i)
  }
  annotation <- cbind.data.frame(begintime,endtime,annotation)
  
  #MERGE GESTURES based on mindistformerge
  if(nrow(annotation) > 1) #only if there are more than 2 gestures
    #check close gestures
  {
    for(ge in c(1:(nrow(annotation)-1)) )
    {
      #endtime gesture  - begintime next gesture
      if( (annotation[ge+1,1]-annotation[ge,2]) < mindistformerge )
      {
        newbegintime <- annotation[ge,1]
        annotation[ge,] <- NA
        #delete current gesture
        annotation[ge+1,1] <- newbegintime
      }
    }
    annotation <- annotation[!is.na(annotation[,3]),]
  }
  
  #check length and set duration (if gestures are very short then delete)
  annotation <-     annotation[(annotation[,2]-annotation[,1]) > minduration,]
  return(annotation)
}


#FUNCTION for computing 2D and 3D speed (and apply butterworth filter)
velocity.it <- function(x, y, z, dimension, t) 
{
  x <- as.numeric(x)
  y <- as.numeric(y)
  t <- as.numeric(t)
  z <- as.numeric(z)
  if(dimension == "3D")
  {
    xyz <- cbind.data.frame(x,y,z,t)
    colnames(xyz) <- cbind("x", "y", "z", "t")
    #euclidean distance calc 3D: ((x1-xn) (y1-yn) (z1-zn))^2
    xyzdiff <- as.data.frame(apply( xyz , 2 , diff ))
    xyzdiff$v <- sqrt(rowSums(xyzdiff[,1:3]^2))
    xyzdiff$v <- xyzdiff$v/xyzdiff$t*1000
    velocity <- c(0, xyzdiff$v)
  }
  if(dimension == "2D")
  {
    xy <- cbind.data.frame(x,y,t)
    colnames(xy) <- cbind("x", "y", "t")
    #euclidean distance calc 2D: ((x1-xn) (y1-yn)^2
    xydiff <- as.data.frame(apply( xy , 2 , diff ))
    xydiff$v <- sqrt(rowSums(xydiff[,1:2]^2))
    xydiff$v <- xydiff$v/xydiff$t*1000
    velocity <- c(0, xydiff$v)
  }
  velocity <- smooth.it(velocity)
  return(velocity)
}


############################Initializing variable names
kinjs <- c("handtip_left", "handtip_right")
kinvars <- vector()
for(i in kinjs)
  for(j in c("_x", "_y", "_z"))
  {kinvars <- c(kinvars, paste0(i,j))}
selectionvariables <- c("time_s", "round", "trial", "trial_nr", "trial_id", kinvars)

############################Initializing settings for the autocoder
#parameters setting for gesture detection (if you want to compare performance of other autocoding settings you can change these)
speed <- c(15)      #faster than 15 cm per second
duration <- c(200)  #minimal duration
mindist <- c(250)   #connect when less then 250ms between gestures
gesturespace <- c(1) #gesturespace*standard deviation under mean
############################################MAIN ROUTINE
  
#go through the motion tracking files, extract only a subset of relevant variables, calculate speed, and then add autoannotations
for(par in MTS)
{
  print(paste0("currently proccesing: ", par))
  id <- substring(par, 1, 3)
  MT <- read.csv(paste0(MTfolder3, par))
  
  #############make a selection
  MT <- MT[,colnames(MT) %in% selectionvariables]
  MT$time_ms <- round(MT$time_s*1000)
  ############SMOOTH TIMESERIES
  MT$handtip_left_x <- smooth.it(MT$handtip_left_x)
  MT$handtip_left_y <- smooth.it(MT$handtip_left_y)
  MT$handtip_left_z <- smooth.it(MT$handtip_left_z)
  MT$handtip_right_x <- smooth.it(MT$handtip_right_x)
  MT$handtip_right_y <- smooth.it(MT$handtip_right_y)
  MT$handtip_right_z <- smooth.it(MT$handtip_right_z)
  
  ############COMPUTE SPEED
  #compute velocities A
  MT$v_handtip_left     <- velocity.it(MT$handtip_left_x, MT$handtip_left_y,MT$handtip_left_z,"3D", MT$time_ms)
  MT$v_handtip_right    <- velocity.it(MT$handtip_right_x, MT$handtip_right_y, MT$handtip_right_z,"3D", MT$time_ms)
  
  
  #################################################AUTO ANNOTATIONS
  MT$autogesturel <- MT$autogesturer <- NA
  
  #compute auto annotations based on finger tip data (see info 'get.autogesture.annotations' at function section)
  print(paste0("working on autoannotations: xx"))
  annotationleftAUTO   <- get.autogesture.annotations(MT$handtip_left_y, speed, gesturespace, duration, mindist, MT$v_handtip_left, MT$time_ms, MT$trial)
  annotationrightAUTO  <- get.autogesture.annotations(MT$handtip_right_y, speed, gesturespace,duration, mindist, MT$v_handtip_right, MT$time_ms,MT$trial)
  
  #if there are gesture detected then load the auto annotations into a the time series
  print(paste0("working on autoannotations: .x"))
  if(nrow(    annotationleftAUTO) > 1) {MT$autogesturel <- load.in.event(MT$time_ms,annotationleftAUTO)}
  if(nrow(    annotationrightAUTO) > 1){MT$autogesturer <- load.in.event(MT$time_ms,annotationrightAUTO  )}
  print(paste0("working on autoannotations: DONE & WRITING FILES"))
  
  
  unique(ave(MT$autogesturer, MT$autogesturer, FUN=length))
  ########################################################WRITING AUTOCODED TIME SERIES FILES
  write.csv(MT, paste0(wrifolder4, id, "_MT_processed",".csv")) #Here we save the final output to be used for full analysis

  #save auto_annotations
  annor <- cbind(annotationrightAUTO, "right_hand")
  annol <-cbind(annotationleftAUTO, "left_hand")
  colnames(annor) <- colnames(annol) <- c("begintime", "endtime", "gesture_id", "hand")
  c_annotations <- rbind.data.frame(annor, annol)
  write.csv(c_annotations, paste0(wrifolder4_1, id, "_auto_annotations.csv")) #save the annotation files
  
}