#A script to clean NextFlow outputs for final processing
#Developed by Dr. Jake Beierle (don't forget the Dr., it's important)

#----Documentation----
#See comprehensive documentation on the github repository
#https://github.com/jacobbeierle/JABS_nextflow_postprocess/tree/main

#Import libraries---------------------------------------------------------------
library(tidyverse)
library(writexl)
library(janitor)
library(data.table)
options(error = NULL) #helps with error handling in functions checking for directories and filenames
#########################################################Set QC and other Values------------------------------------------------------------------

##Set your working directory
#working.directory <- "C:\\Users\\beierj\\Desktop\\2025-04-09_NTG_C1-C5_Analysis"
working.directory <- "Y:\\kumarlab-new\\Jake\\Expts_Comp\\2025-10-29_OW_Pilot_V1-4_WS1-4"

##Set the QC values you will use to screen the Nextflow QC files
#Expected length of clip in frames (video clipping duration plus 5 seconds)
expected.length <- 150*60*30 + 5*30
#Max number of tracklets/hour recording
max.tracklets.per.hour <- 6
#Max frames (as percent of all frames) missing pose
max.percent.segmentation.missing <- 0.2
#Max percentage of Frames missing pose
max.percent.pose.missing <- 0.005
#Max percentage of Frames missing pose
max.percent.kp.missing <- 0.01

#Proportion of Highest and lowest fecal boli mice to plot separately for QC
fecal_boli_percent_threshold <- 0.05


#Create functions---------------------------------------------------------------
#Creates a directory if it does not already exist
directory_check_creation <- function(x){
  x <- paste0(x)
  if(!dir.exists(x)){
    dir.create(x)
  }
  else{print(paste0("'/", substitute(x), "' already exists"))}
}


#Set working directory and create qc directory----------------------------------

#If you have defined a working directory above, set it here
if(exists("working.directory")){
  setwd(working.directory)
}

#Check to ensure there is a nexflow ouput folder
if(!dir.exists("Nextflow_Output")){
  stop("YOU DO NOT HAVE A 'Nextflow_Output' DIRECTORY IN YOUR WORKING DIRECTORY /n
       YOU NEED TO CREATE A 'Nextflow_Output' DIRECTORY WITH YOUR RESULTS") 
}else{
  directory_check_creation(file.path("Nextflow_Output", "final_nextflow_feature_data"))
  final.dataframes.dir <- file.path("Nextflow_Output", "final_nextflow_feature_data")
}


#create QC directories if not made, programatically define these directories for publishing
directory_check_creation("qc")
#directory for QC logs
directory_check_creation(file.path("qc", "nextflow_qc_logs"))
qc.log.dir <- file.path("qc", "nextflow_qc_logs")
#directory for missing or duplicated data logs
directory_check_creation(file.path("qc", "missing_or_dup_data"))
qc.missing_dup.dir <- file.path("qc", "missing_or_dup_data")
#directory for QC figures
directory_check_creation(file.path("qc", "qc_figs"))
qc.figs.dir <- file.path("qc", "qc_figs")




#Process and Publish QC logs with success or failure annotated in a CSV-------------------------

#read QC files in NextFlow_Output directory
qc_log <- list.files(
              path = "NextFlow_Output/",
              pattern = "qc_batch_",
              full.names = TRUE,
              recursive = TRUE
              ) |> 
          read_csv(id = "QC_file")

#Create subset of data failing QC measures

#Record why QC failed for each video
qc_log$passed_duration_QC <- qc_log$video_duration == expected.length
qc_log$passed_tracklet_QC <- qc_log$pose_tracklets < max.tracklets.per.hour*(expected.length/108000)
qc_log$passed_segmentation_QC <- qc_log$seg_counts > (1-max.percent.segmentation.missing) * expected.length
qc_log$passed_pose_QC <- qc_log$pose_counts > (1-max.percent.pose.missing) * expected.length
qc_log$passed_kp_QC <- qc_log$missing_keypoint_frames < max.percent.kp.missing * expected.length

#Apply thresholds defined above to create a seperate 'failed QC' data frame
qc_log.failed <- subset(qc_log, video_duration != expected.length |
                      pose_tracklets > max.tracklets.per.hour*(expected.length/108000) |
                      seg_counts < (1-max.percent.segmentation.missing) * expected.length |
                      pose_counts < (1-max.percent.pose.missing) * expected.length|
                      missing_keypoint_frames > max.percent.kp.missing * expected.length)

#Write final Nextflow QC files for review by a human
write.csv(qc_log, file.path(qc.log.dir , "qc_all.csv"), row.names = FALSE)
write.csv(qc_log.failed,  file.path(qc.log.dir,"qc_failed.csv"), row.names = FALSE)


#Create a list of expected videos that would be in each Nextflow output---------------------------------------------------------
#Create a dataframe of videos in QC log
expected_vidoes <- as.data.frame(gsub("_with_fecal_boli", ".avi", qc_log$video_name))
expected_vidoes <- as.data.frame(gsub("_filtered", "", expected_vidoes[,1]))
colnames(expected_vidoes) <- "NetworkFilename"
#Make sure there are no duplicates
expected_vidoes <- distinct(expected_vidoes, NetworkFilename)


#Process Fecal Boli Data and publish QC plots--------------------------------------------------------
fecal_boli.raw <- list.files(
            path = "NextFlow_Output/",
            pattern = "fecal_boli.csv",
            recursive = TRUE,
            full.names = TRUE) |> 
  read_csv() |> 
  select(-nextflow_version)

colnames(fecal_boli.raw)[1] <- "NetworkFilename"
#Correct some discrepancies that may have arose in file names due to corner correction workflow
fecal_boli.raw$NetworkFilename <- gsub("_corrected", "", fecal_boli.raw$NetworkFilename)
fecal_boli.raw$NetworkFilename <- gsub("_filtered", "", fecal_boli.raw$NetworkFilename)

#Check for missing data in fecal boli data based on QC files
videos_with_missing_fecal_boli <- as.data.frame(setdiff(expected_vidoes$NetworkFilename, fecal_boli.raw$NetworkFilename))
colnames(videos_with_missing_fecal_boli) <- "NetworkFilename"
fecal_boli_videos_missing_in_qc <- as.data.frame(setdiff(fecal_boli.raw$NetworkFilename, expected_vidoes$NetworkFilename))
colnames(fecal_boli_videos_missing_in_qc) <- "NetworkFilename"

#Output warning and prepare to report report CSV if files are missing in QC log, but present in Gait analysis
if(length(fecal_boli_videos_missing_in_qc!=0)){
  fecal_boli_videos_missing_in_qc$fboli <- 1
}

#Output warning and report CSV if files are missing in QC log, but present in Gait analysis
if(length(videos_with_missing_fecal_boli!=0)){
  videos_with_missing_fecal_boli$missing_fboli <- 1
}

#Check for fecal boli rows with identical data (i.e. something went wrong in video recording)
#I remove the network file name col because data may have been misslabeled
#This is used to output QC measures in the final portion of this script
#There are many more duplicate rows here, and not nessisarily cause for alarm
duplicate_fboli_rows <- fecal_boli.raw[duplicated(fecal_boli.raw[-1]) | duplicated(fecal_boli.raw[-1], fromLast = TRUE), ]

#Pivot longer to facilitate plotting for QC
fecal_boli.plot <- fecal_boli.raw |> 
  pivot_longer(
    cols = !NetworkFilename,
    names_to = "min", 
    values_to = "fecal_boli",
    values_drop_na = TRUE) |> 
  mutate(min = parse_number(min))

#Plot fecal boli QC measures 
outFileNamePDF <- file.path(qc.figs.dir, "fecal_boli_qc_figs.pdf") 
pdf(outFileNamePDF, 6, 6)

#Growth curve for all mice
p1 <- ggplot(fecal_boli.plot, aes(min, fecal_boli, group = NetworkFilename, colour = NetworkFilename))+
  geom_line() +
  labs(title = "Fecal boli growth, all mice") +
  theme(legend.position = "none")
print(p1)

#Plot 10% mice with lowest fecal boli
p1 <- fecal_boli.plot  |> 
  summarise(across(fecal_boli, max), .by = NetworkFilename) |> 
  slice_min(fecal_boli, prop = fecal_boli_percent_threshold) |> 
  select(NetworkFilename) |> 
  merge(fecal_boli.plot, by.x = "NetworkFilename") |> 
  ggplot(aes(min, fecal_boli, group = NetworkFilename, colour = NetworkFilename))+
    geom_line() +
    labs(title = paste("Lowest ", fecal_boli_percent_threshold*100, "% of fecal boli mice", sep = "")) +
    theme(legend.position = "none")
print(p1)

#Plot 10% mice with highest fecal boli
p1 <- fecal_boli.plot  |> 
  summarise(across(fecal_boli, max), .by = NetworkFilename) |> 
  slice_max(fecal_boli, prop = fecal_boli_percent_threshold) |> 
  select(NetworkFilename) |> 
  merge(fecal_boli.plot, by.x = "NetworkFilename") |> 
  ggplot(aes(min, fecal_boli, group = NetworkFilename, colour = NetworkFilename))+
    geom_line() +
    labs(title = paste("Highest ", fecal_boli_percent_threshold*100, "% of fecal boli mice", sep = "")) +
    theme(legend.position = "none")
print(p1)

#Histogram of final fecal boli count
p1 <- fecal_boli.plot  |> 
  arrange(desc(min)) |> 
  distinct(NetworkFilename, .keep_all = TRUE) |> 
  ggplot(aes(fecal_boli, ifelse(after_stat(count) > 0, after_stat(count), NA)))+
    geom_histogram(binwidth = 1, boundary = 0)+
    labs(title = "Fecal boli highest bin, all mice") +
    ylab("count")
print(p1)

dev.off()

#Write out all raw, merged fecal boli counts
write.csv(fecal_boli.raw, file.path(final.dataframes.dir, "fecal_boli_raw.csv"), row.names = FALSE)

#Process Gait Data--------------------------------------------------------------

#Import Gait Data
gait.raw <- list.files(
  path = "NextFlow_Output/",
  pattern = "gait.csv",
  full.names = TRUE,
  recursive = TRUE) |> 
  read_csv() |> 
  select(-nextflow_version)
colnames(gait.raw)[1] <- "NetworkFilename"


#Amend video path information to match other Nextflow outputs, making merging easier later
gait.raw$NetworkFilename <- sub(".", "", gait.raw$NetworkFilename)

#Check for missing data in Gait and QC files
videos_with_all_gait_missing <- as.data.frame(setdiff(expected_vidoes$NetworkFilename, gait.raw$NetworkFilename))
colnames(videos_with_all_gait_missing) <- "NetworkFilename"
gait_videos_missing_in_qc <- as.data.frame(setdiff(gait.raw$NetworkFilename, expected_vidoes$NetworkFilename))
colnames(gait_videos_missing_in_qc) <- "NetworkFilename"

#Output warning and report CSV if files are missing in QC log, but present in Gait analysis
if(length(gait_videos_missing_in_qc!=0)){
  gait_videos_missing_in_qc$gait <- 1
}

if(length(videos_with_all_gait_missing!=0)){
  videos_with_all_gait_missing$gait <- 1
}

#Extract values repeated identically for each time bin (i.e. the cols identified below) and reduce to one per video
gait.duplicated_data <- gait.raw |> 
  select(c(NetworkFilename, `Distance Traveled`, `Body Length`, Speed, `Speed Variance`)) |>
  distinct(NetworkFilename, .keep_all = TRUE)

#Unmelt the data and add speed bin to col names
#Remove duplicated measures handled above
gait.speed_bin_data <- gait.raw |> 
  select(!c(`Distance Traveled`, `Body Length`, Speed, `Speed Variance`))

#Remove variance measures from speed bins with fewer than 3 strides
#These do not make statistical sense to report
gait.speed_bin_data[gait.speed_bin_data$`Stride Count` < 2, colnames(gait.speed_bin_data) %like% "Variance"] <- NA

#"Unmelt" non-duplicated data by speed bin
gait.speed_bin_data <-  dcast(as.data.table(gait.speed_bin_data),
                              formula = NetworkFilename ~ gait.speed_bin_data$'Speed Bin', 
                              value.var = colnames(gait.speed_bin_data)[3:(ncol(gait.speed_bin_data)-1)],
                              sep = ".")

#Merge binned and duplicated measur1es
gait.merged <- merge(gait.duplicated_data, gait.speed_bin_data, by = "NetworkFilename")

#Select all unique NetworkFilenames in the qc Log, and create empty rows unrepresented in the gait features
gait.final <- merge(expected_vidoes, gait.merged, by = "NetworkFilename", all = TRUE)
gait.final[is.na(gait.final$`Stride Count.10`) , 'Stride Count.10'] <- 0
gait.final[is.na(gait.final$`Stride Count.15`) , 'Stride Count.15'] <- 0
gait.final[is.na(gait.final$`Stride Count.20`) , 'Stride Count.20'] <- 0
gait.final[is.na(gait.final$`Stride Count.25`) , 'Stride Count.25'] <- 0

#Remove Speed bin cols, which are no longer useful
gait.final <- select(gait.final, !contains('Speed Bin'))

#Check for gait rows with identical data (i.e. something went wrong in video recording)
#This is used to output QC measures in the final portion of this script
duplicate_gait_rows <- as_tibble(gait.final[duplicated(gait.final[-1]) | duplicated(gait.final[-1], fromLast = TRUE), ])

#output to final CSV
write.csv(gait.final, file.path(final.dataframes.dir,"gait_final.csv"), row.names = FALSE)


#Process JABS Feature Data NEEDS DEVELOPMENT------------------------------------
JABS.features <- list.files(
  path = "NextFlow_Output/",
  pattern = "features.csv",
  recursive = TRUE,
  full.names = TRUE) |> 
  read_csv() |> 
  select(-nextflow_version)
colnames(JABS.features)[1] <- "NetworkFilename"

#Adjust networkfile names

JABS.features$NetworkFilename <- paste0("/", JABS.features$NetworkFilename, ".avi")
JABS.features$NetworkFilename <- gsub("_corrected", "", JABS.features$NetworkFilename)
JABS.features$NetworkFilename <- gsub("_filtered", "", JABS.features$NetworkFilename)


#Check for missing data in Gait and QC files
videos_with_JABS_features_missing <- as.data.frame(setdiff(expected_vidoes$NetworkFilename, JABS.features$NetworkFilename))
colnames(videos_with_JABS_features_missing) <- "NetworkFilename"
JABS_features_videos_missing_in_qc <- as.data.frame(setdiff(JABS.features$NetworkFilename, expected_vidoes$NetworkFilename))
colnames(JABS_features_videos_missing_in_qc) <- "NetworkFilename"

#Output warning and prepare to report report CSV if files are missing in QC log, but present in Gait analysis
if(length(JABS_features_videos_missing_in_qc!=0)){JABS_features_videos_missing_in_qc$JABS_features <- 1}

#Output warning and report CSV if files are missing in QC log, but present in Gait analysis
if(length(videos_with_JABS_features_missing!=0)){videos_with_JABS_features_missing$JABS_features <- 1}


#Check for JABs feature rows with identical data (i.e. something went wrong in video recording)
#This is used to output QC measures in the final portion of this script
#JABS features are rounded to the nearest tenth because of minor differences in estimates
duplicate_JABS_feature_rows <- JABS.features[duplicated(sapply(JABS.features[-1], \(x) round(x, digits = 0))) |
                                               duplicated(sapply(JABS.features[-1], \(x) round(x, digits = 0)), fromLast = TRUE), ]
#Write the final csv
write.csv(JABS.features, file.path(final.dataframes.dir, "JABS_features_final.csv"), row.names = FALSE)


#Process morphometrics feature data---------------------------------------------
morpho.raw <- list.files(
  path = "NextFlow_Output/",
  pattern = "morphometrics.csv",
  full.names = TRUE, 
  recursive = TRUE) |> 
  read_csv() |> 
  relocate(NetworkFilename) |> 
  select(!c(nextflow_version))

#reformat NetworkFilename
morpho.raw$NetworkFilename <- sub(".", "", morpho.raw$NetworkFilename)

#Check for missing data in Gait and QC files
videos_with_morphometrics_features_missing <- as.data.frame(setdiff(expected_vidoes$NetworkFilename, morpho.raw$NetworkFilename))
colnames(videos_with_morphometrics_features_missing) <- "NetworkFilename"
morphometrics_videos_missing_in_qc <- as.data.frame(setdiff(morpho.raw$NetworkFilename, expected_vidoes$NetworkFilename))
colnames(morphometrics_videos_missing_in_qc) <- "NetworkFilename"

#Output warning and prepare to report report CSV if files are missing in QC log, but present in Gait analysis
if(length(morphometrics_videos_missing_in_qc!=0)){morphometrics_videos_missing_in_qc$morpho_features <- 1}

#Output warning and report CSV if files are missing in QC log, but present in Gait analysis
if(length(videos_with_morphometrics_features_missing!=0)){videos_with_morphometrics_features_missing$morpho_features <- 1}


#Check for JABs feature rows with identical data (i.e. something went wrong in video recording)
#This is used to output QC measures in the final portion of this script
#JABS features are rounded to the nearest tenth because of minor differences in estimates
duplicate_morphometrics_rows <- morpho.raw[duplicated(morpho.raw[-1]) | duplicated(morpho.raw[-1], fromLast = TRUE), ]


write.csv(morpho.raw, file.path(final.dataframes.dir, "morphometrics_final.csv"), row.names = FALSE)



#Merge and output all data types missing for videos-----

#Publish all vids with missing data
#Combine into a single list
all_missing_data <- list(
  "videos_with_missing_fecal_boli" = videos_with_missing_fecal_boli,
  "videos_with_all_gait_missing" = videos_with_all_gait_missing,
  "videos_with_JABS_features_missing" = videos_with_JABS_features_missing,
  "videos_with_morphometrics_features_missing" = videos_with_morphometrics_features_missing
)

#Select the dfs with actual missing data represented, dropping the rest
publish_missing_data <- NULL
for(i in seq_along(all_missing_data)){
  if(length(all_missing_data[[i]]) > 0){
    if(length(publish_missing_data) == 0){
      publish_missing_data <- all_missing_data[[i]]
    }else{publish_missing_data <- full_join(publish_missing_data, all_missing_data[[i]])}
  }
}

#Fill the csv with something if all QC passes
if(length(publish_missing_data) == 0){
  publish_missing_data <- "NO DATA MISSING"
}

#Write the csv
write.csv(publish_missing_data, file.path(qc.missing_dup.dir, "missing_data.csv"), row.names = FALSE)


#Publish all vids in some feature csv, but not in the qc dataframe
#Combine into a single list
videos_not_in_qc_report <- list(
  "fecal_boli_videos_missing_in_qc" = fecal_boli_videos_missing_in_qc,
  "gait_videos_missing_in_qc" = gait_videos_missing_in_qc,
  "JABS_features_videos_missing_in_qc" = JABS_features_videos_missing_in_qc,
  "morphometrics_videos_missing_in_qc" = morphometrics_videos_missing_in_qc
)

#Select the dfs with actual missing data represented, dropping the rest
publish_videos_not_in_qc_report <- NULL
for(i in seq_along(videos_not_in_qc_report)){
  if(length(videos_not_in_qc_report[[i]]) > 0){
    if(length(publish_videos_not_in_qc_report) == 0){
      publish_missing_data <- videos_not_in_qc_report[[i]]
    }else{publish_missing_data <- full_join(publish_missing_data, videos_not_in_qc_report[[i]])}
      
  }else{}
}

#Fill the csv with something if all QC passes
if(length(publish_videos_not_in_qc_report) == 0){
  publish_videos_not_in_qc_report <- "NO DATA MISSING"
}
write.csv(publish_videos_not_in_qc_report, file.path(qc.missing_dup.dir, "videos_not_in_qc_report.csv"), row.names = FALSE)



#Output data that is duplicated in the data frames
duplicated_data <- list(
  "dup_gait" = duplicate_gait_rows,
  "dup_JABS" = duplicate_JABS_feature_rows,
  "dup_morpho" = duplicate_morphometrics_rows,
  "dup_fboli" = duplicate_fboli_rows
)

#If all objecst in the list have a length of 0 (i.e. empty), report no dupli
if(all(sapply(duplicated_data, function(x) nrow(x)==0))){
  duplicated_data <- "NO DUPLICATED DATA"
  write_xlsx(as.data.frame(duplicated_data), path = file.path(qc.missing_dup.dir, "duplicated_data.xlsx"))
}else{
  write_xlsx(duplicated_data, path = file.path(qc.missing_dup.dir, "duplicated_data.xlsx"))
  
}


#Write warnings for failed QC--------------------------------------------------------------
error.reporting <- NULL

#Report to terminal if data is missing from QC log but in feature tables
if(!is.character(publish_videos_not_in_qc_report)){
  if(length(fecal_boli_videos_missing_in_qc)){ error.reporting <- c(error.reporting,"FECAL BOLI DATA PRESENT FOR VIDEOS NOT IN QC LOG") }
  if(length(gait_videos_missing_in_qc)){ error.reporting <- c(error.reporting,"GAIT DATA PRESENT FOR VIDEOS NOT IN QC LOG") }
  if(length(JABS_features_videos_missing_in_qc)){ error.reporting <- c(error.reporting,"JABS FEATURE DATA PRESENT FOR VIDEOS NOT IN QC LOG") }
  if(length(morphometrics_videos_missing_in_qc)){ error.reporting <- c(error.reporting,"MORPHOMETRIC DATA PRESENT FOR VIDEOS NOT IN QC LOG") }
}

#Report to terminal if data is missing feature tables but in QC log
if(!is.character(publish_missing_data)){
  if(length(videos_with_missing_fecal_boli)){error.reporting <- c(error.reporting,"YOU ARE MISSING FECAL BOLI DATA") }
  if(length(videos_with_JABS_features_missing)){ error.reporting <- c(error.reporting,"YOU ARE MISSING JABS FEATURE DATA") }
  if(length(videos_with_morphometrics_features_missing)){ error.reporting <- c(error.reporting,"YOU ARE MISSING MORPHOMETRIC FEATURE DATA") }
  if(length(videos_with_all_gait_missing)){ error.reporting <- c(error.reporting,"YOU HAVE VIDEOS WITH NO GAIT DATA, MAY BE INACTIVE MICE") }
}

#Report to terminal if there is duplicated data
if(!is.character(duplicated_data)){ error.reporting <- c(error.reporting,"YOU HAVE DUPLICATED DATA!") }


#Print out all errors after code done running
if(length(error.reporting) == 0){
  print("FINAL ERROR REPORT: NO ERRORS TO REPORT")
}else{
  print("FINAL ERROR REPORT:", )
  paste(error.reporting)
}
