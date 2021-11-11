library(rstudioapi) #foldering
library(tidyr)      #data wrangling
library(signal)     #signal processing
library(pracma)     #signal processing
library(kza)        #smoothing
library(zoo)        #interpolation

############################
#Author Script: Wim Pouw (wimpouw@gmail.com)

############################

#folders and datafiles
basefolder <- dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path))) #get path of current R doc
rawfolder  <- paste0(basefolder, "/2_SYNCED_RAW_KINECT/")
rawfiles   <- list.files(rawfolder, pattern = "*_WIDE.txt")
annofolder <- paste0(basefolder, "/0_ANNOTATIONDATA/")
annofiles  <- list.files(annofolder, pattern = "*accuracy.txt")
wrifolder3 <- paste0(basefolder, "/3_SYNCED_KINECT_TIMESTAMPED/")

#############################FUNCTIONS

# we are using a Kolmogorov-Zurbenko  filter with an order of 2 and a span of 3
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

############################Initializing variable names
kinectjoints <-
  c("spinebase",
    "spinemid",
    "neck",
    "head",
    "shoulder_left",
    "elbow_left",
    "wrist_left",
    "hand_left",
    "shoulder_right",
    "elbow_right",
    "wrist_right",
    "hand_right",
    "hip_left",
    "knee_left",
    "ankle_left",
    "foot_left",
    "hip_right",
    "knee_right",
    "ankle_right",
    "foot_right",
    "spineshoulder",
    "handtip_left",
    "thumbleft",
    "handtip_right",
    "thumb_right")

#make a variable list taking into account that we have x,y,z data
variablelist <- c("time_s") #the first column is synchronized time in seconds
for(i in kinectjoints)
  for(j in c("_x", "_y", "_z"))
  {variablelist <- c(variablelist, paste0(i, j))}
print(variablelist)

#the timings are sometimes in seconds and sometimes in milliseconds (for later pairs), so check what it is it and adjust accordingly
correct.samplingrate <- function(time_s)
{
time_sa <- time_s  
if(mean(unique(diff(time_s)))<10) #check sampling rate
{
  return(time_sa)
}
  if( mean(unique(diff(time_s)))>10) #if the sampling intervals exceed 10 units, it cant be seconds, and it must be milliseconds
  {
    print(paste0("pair", par, "is actually in milliseconds; converted to seconds"))
    time_sa <- round(time_s/1000,3)
  }
}

############################################MAIN ROUTINE
  
#go through the files
for(par in rawfiles)
{
  print(paste0("currently proccesing: ", par))
  MT <- read.delim(paste0(rawfolder, par), sep ="\t", header= FALSE)
  colnames(MT) <- variablelist
  MT$time_s <- correct.samplingrate(MT$time_s) #check timing variable and adjust accordingly to seconds
  
  #there are duplicate time stamps in the data (we remove these) #this si also giving a warning that the interpolation only takes unique values
  MT <- MT[!duplicated(MT$time_s),]
  
  #regularize sampling at exactly 30Hz using interpolation
  MTnewtime <- as.data.frame(round(seq( min(MT$time_s), max(MT$time_s), by = 1/30),3)) #make a new time variable
  colnames(MTnewtime) <- c("time_s")
  MT <- merge(MTnewtime,MT, by.x="time_s", by.y = "time_s", all.x = TRUE, all.y = TRUE) #merge time variable with data
  MT[,2:ncol(MT)]<- apply(MT[,2:ncol(MT)], 2, FUN = function(x) na.approx(x, x=MT$time_s,na.rm=FALSE)) #interpolate na's
  MT <- MT[which(MT$time_s%in%MTnewtime$time_s),] #only keep new time variable with interpolated data
  MT[,2:ncol(MT)]<- apply(MT[,2:ncol(MT)], 2, FUN = function(x) smooth.it(x)) #smooth interpolated data

  #add trial info
      #select the relevant annotationfile
  id <- substring(par, 1, 3)
  annotationfile <- annofiles[grepl( id,  substring(annofiles, 1,3), fixed = TRUE)|grepl( id,  substring(annofiles, 4,6), fixed = TRUE)] #check if id is in the filename of the annotation
  annot <- read.delim(paste0(annofolder, annotationfile), sep =",", header=TRUE)
  print(paste0("loading in anntations: xxxx"))
  #save the key annotation info to the time series (note you can add more variables here if you like)
  MT$pair     <- unique(annot$pair)
  MT$trial    <- load.in.event(MT$time_s*1000, cbind.data.frame(annot$onset_msec, annot$offset_msec, annot$trial))
  print(paste0("loading in anntations: .xxx"))
  MT$round    <- load.in.event(MT$time_s*1000, cbind.data.frame(annot$onset_msec, annot$offset_msec, annot$round))
  print(paste0("loading in anntations: ..xx"))
  MT$trial_nr <- load.in.event(MT$time_s*1000, cbind.data.frame(annot$onset_msec, annot$offset_msec, annot$trial_nr))
  print(paste0("loading in anntations: ...x"))
  MT$trial_id <- load.in.event(MT$time_s*1000, cbind.data.frame(annot$onset_msec, annot$offset_msec, annot$trial_id))
  print(paste0("loading in anntations: DONE!"))
  
  #write the timestamped file to folder 3 
  print(paste0("writing file to folder 3: x"))
  write.csv(MT, paste0(wrifolder3, id, "_TIMESTAMPED", ".csv"),row.names = FALSE)  
  print(paste0("writing file to folder 3: DONE"))
}