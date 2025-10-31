# JABS_nextflow_postprocess

This repository is R code developed by Jake to take raw single mouse nextflow outputs from the JABS pipeline and conduct a series of qc screening, data visualization, and post processing to prepare the data for its final analysis. It currently consists of 5 scripts that should be run in order, denoted by the number at the suffix of the file name. They are described below.

## NextFlow_Output_QC_Postprocess_1.R
The goal of this script is to read in the disperate data from the nextflow pipeline, ensure all data from nextflow QC reports are present, post process the gait (descirbed below), screen for duplicate data, format them for easier merging, and produce QC figures for the raw fecal boli data.

### This script makes the following assumptions:
1. A project folder containing a directory with NextFlow Output inside titled **NextFlow_Output/**, see below for example file structure. The content in this directory can be nested, as files are searched for recursively. Files in this directory should include one or more: QC report, gait, unprocessed fecal boli, morphometrics, JABS features. These files can be in nested directories, but need the following words in their file names: "qc_batch_", "gait", "fecal_boli", "feature".
2. Missing gait videos/bins are a consequence of no predicted gaits, not because gait was not predicted on these videos.

### Before running this code:
1. Set the QC values to your chosen thresholds in the 'Set QC and other Values' section
2. Set your working directory in the in the 'Set QC and other Values' section

### This code will generate:
1. A new **/qc** directory, where many of the qc files will be placed. This includes:
	1. A **qc_all.csv** and **qc_failed.csv** file that include metrics about what videos passed and failed QC in **/nextflow_qc_logs**
	2. A **missing_data.csv** file listing video data that exists in the qc files output by nextflow but not in one or more of your feature files
	3. A **videos_not_in_qc_report.csv** listing video data that exists in your feature csv files but not your qc files output by nextflow
	4. A **duplicated_data.xlsx** spreadsheet listing videos/rows from feature csv files that have identical data (i.e. duplicate rows that may, or may not, have the same NetworkFilename). *Note, there will probably be some for fecal boli that are not cause for alarm*
	5. A **fecal_boli_qc_figs.pdf** file with figures from the raw fecal boli data to be manually reviewed for problematic fecal boli data.
2. Final csv files that will be used for coalating your final dataframe in the next script. They will be in '/WORKING_DIRECTORY/Nextflow_Output/final_nextflow_feature_data/' and include:
	1. **gait_final.csv** -- with edited NetworkFilenames for easy merging in next script, widened so that each feature in each speed bin is its own col, with features repeated (i.e. identical) across speed bins removed, with 'Stride count' for missing speed bins padded with 0s (not NAs), and with variances set to NA for speed bins with less than 3 strides.
	2. **morphometrics_final.csv** -- with edited NetworkFilenames for easy merging in next script
	3. **JABS_features_final.csv** -- with edited NetworkFilenames for easy merging in next script
	4. **fecal_boli_raw.csv** -- with edited NetworkFilenames for easy merging in next script, note that you, the user, will have to review the fecal boli QC plots, check for outliers and odd patterns, and manually correct those counts within '/qc/fecal_boli_raw.csv', relabel this file as "fecal_boli_final.csv" and move it to '/NEXTFLOW_PROJECT_FOLDER'

#### Before moving to *merge_all_dataframes_2.R* you need to:
1. Review the **qc_failed.csv**, and manually overlay pose on problematic files, decide which videos to include and exclude.
2. Create a text file, **videos_to_exclude.txt**, in your /WORKING_DIRECTORY including the 'NetworkFilename' of the videos to exclude based on your manual QC. A note, I usually manually create a file, **qc_failed_omitted.csv** recording the reason I excluded videos.
3. Review **fecal_boli_qc_figs.pdf**, check for outliers and odd patterns, and manually correct those counts that have errors. Relabel this file as **fecal_boli_final.csv** and move to '/WORKING_DIRECTORY/Nextflow_Output/final_nextflow_feature_data/'

## merge_all_dataframes_2.R
The goal of this script is to merge all final datasets from nextflow, ensure all videos are represetned in all files, and that all mice are represented in the metadata (and that all mice in the metadata are represetned in the data). Its final output is a **merged_nextflow_dataset.csv** that still needs additional QC from *final_dataframe_QC_3.R*

### This script makes the following assumptions:
1. You have generated the following csv files from the previous script *NextFlow_Output_QC_Postprocess_1.R*, which sould be in '/NEXTFLOW_PROJECT_FOLDER/Nextflow_Output/final_nextflow_feature_data/':
	1. **gait_final.csv**
	2. **morphometrics_final.csv**
	3. **JABS_features_final.csv**
 	4. **fecal_boli_final.csv** -- modified from 'qc/fecal_boli_raw.csv' based on the fecal boli QC plots
2.  You have a **metadata.csv** file within the '/WORKING_DIRECTORY' with a col 'MouseID' that matches a substring in your NetworkFilename cols

### Before running this code you need to:
1. Run *NextFlow_Output_QC_Postprocess_1.R* and generate all it's outputs.
2  Create a text file, **videos_to_exclude.txt**, in your /WORKING_DIRECTORY including the 'NetworkFilename' of the videos to exclude based on your manual QC
3. Review **fecal_boli_qc_figs.pdf**, check for outliers and odd patterns, and manually correct those counts that have errors. Relabel this file as **fecal_boli_final.csv** and move to '/WORKING_DIRECTORY/Nextflow_Output/final_nextflow_feature_data/'
4.  Set your working directory in the in the 'Set QC and other Values' section

### This code will generate:
1. A **NetworkFilenames_missing_in_data.csv** file comparing all NetworkFilenames in all data csv files and reporting any discrepancies.
2. A **mice_missing_in_metadata.csv** file comparing all mice (i.e. MouseID), in metadata and final merged dataset and report values missing from one or the other. Failure to match all mice will not stop merging of final dataset, as sometimes mice are included in metadata/runsheets and not excluded if they are exited before testing begins.
3. A **merged_nextflow_dataset.csv** that includes metadata, fecal boli, JABS features, gait, and morphometrics. This file will contain a row for each video processed, and further carpentry will take place in subsequent code. This csv file will be output in the parent/working directory.

### Before moving to *final_dataframe_QC_3.R*:
1. Review the generated qc tables and errors generated at the end of the script to ensure dataframes are congruent across final nextflow dataset

## final_dataframe_QC_3.R
The goal of this script is to prepare your dataframe for final analysis, and therefore will have to be edited and amended for each new project, or for each new paradigm. This script should serve as a robust template for your own needs and highlight some of the final post processing you want to conduct in your own analysis. It includes removing phenotypes with zero variance, phenotypes you manually choose, qc figures to identify any problematic mice/phenotypes you may have missed in the prior steps. Good luck, and godspeed.

### This script makes the following assumptions:
1. You have generated the 'merged_nextflow_dataset.csv' from *merge_all_dataframes_2.R* in the /WORKING_DIRECTORY/Nextflow_Output/final_nextflow_feature_data/ directory
2. You have used the previous 2 scripts in this pipeline and have therefore generated the appropriate directories within your working/parent directory.

### Before running this code you need to:
1. Within the 'Set QC and other Values' section:
   1. Set your working directory
   2. Define a project name to prefix your final dataset
   3. Define the strings or substrings of features you wish to manually exclude with the ```features.removed.manually``` variable. This could include features you do not trust, duplicate/redundant locomotor features, etc.
   4. Define a threshold for automated z-score outlier detection and figure generation
   5. Define values to determine the number of highest z-score outliers to plot and number of phenotyes with the highest number of outliers to plot
2. Within the 'Adjust your data frame to align with your analysis' section write code to adjust your final data into a format for more meaningful analysis. This could include:
   	1. Creating a date tested, or day of testing, from have multiple days of testing. This may be encoded deep within the mess that is the ```NetworkFilenames``` col.
   	2. Perhaps you have doses you wish to exclude from the final analysis in your treatment col.
3. Within the 'Summary information reporting' create code to output bespoke summaries of data points for your review.
4. Within the 'Preliminary QC figures: outlines from linear model of 2 phenotypes' section, decide if you have other comparisions you wish to add. The function in this section creates a scatter plot from 2 phenotypes and labels the 5 points with the largest orthoginal distance from the linear line of best fit. This is another opportunity for manual review.

### If you follow my general format, this code will generate:
1. A ***YOUR_PROJECT_HERE_final_nextflow_dataset.csv'***.
2. A **features_removed_from_curated_dataset.csv** file including all the features you excluded from the final dataset, and why, in your working dir.
2. Several qc figures and one table, which include:
   1. Scatter plots of 2 phenotypes, with the 5 points furthest from the linear line of best fit labled.
   2. A series of box plots visualizing outliers, as defined by a z-score more extreme than the threshold you provided. This includes
		1. 3 pdf files with boxplots of phenotypes with outliers, subset by the nextflow output csv file the phenotypes came from (i.e. gait, JABS features, morphometrics). This excludes fecal boli.
		2. 2 pdf files with boxplots of phenotypes with the most extreme z-score values and the phenotypes with the largest number of outlier z-score values.
		3. A  **outlier_videos_summary.csv** that summarizes the frequency of z-score outliers across all mice
3. Whatever bespoke summaries you may have created below.

### Before moving to *corr_heatmaps_4.R* or *phenotype_exploration_5.R* you need to:
1. Review the QC figures described above and decide if you need to exclude any additional videos or mice. If you decide you should, you will have to re-run *merge_all_dataframes_2.R* or *final_dataframe_QC_3.R* to exclude these phenotypes/mice.

## corr_heatmaps_4.R
The goal of this script is to make heatmaps of phenotypes across several observations. This script is currently under development, and is the least stable.

### This script makes the following assumptions:

### Before running this code you need to:
1. UNDER CONSTRUCTION

### This code will generate
1. Heatmaps, dummy. What, do I have to spell it out for you?

### Before moving to *the afterlife*:
1. Reflect on whether this was really worth it

## phenotype_exploration_5.R
Dawg, I don't even know. What are you doing down here?

## Expected file structure from scripts
Below is an example of how data might be organized and still work with this code:
```
/WORKING_DIRECTORY 
	|NextFlow_Output_QC_Postprocess_1.R
	|merge_all_dataframes_2.R
	|final_dataframe_QC_3.R
	|other_Rscripts.R
	|
	|YOUR_PROJECT_NAME_final_nextflow_dataset.csv ***    <--- the final output for downstream analysis
	|final_n_per_timepoint.csv *** <--- an example of bespoke summaries created in 'final_dataframe_QC_3.R'
	|
	|metadata.csv  &&
	|videos_to_exclude.txt &&
	|features_removed_from_curated_dataset.csv ***
	|
	|/qc *
	|	  |/nextflow_qc_logs
	|     |     |qc_all.csv *
	|     |     |qc_failed.csv *
	|	  |
	|	  |/missing_or_dup_data
	|	  |		 |missing_data.csv *
	|	  |		 |videos_not_in_qc_report.csv *
	|	  |		 |duplicated_data.xlsx *
	|	  |		 |NetworkFilenames_missing_in_data.csv **
	|	  |		 |mice_missing_in_metadata.csv **
	|	  |
	|	  |/qc_figs
	|	  		 |fecal_boli_qc_figs.pdf *
	|	  		 |scatter_plot_lm_figs.pdf ***
	|	  		 |
	|	  		 |/zscore_boxplots
	|	  		  		|gait_boxplot_outlier_figs.pdf ***
	|	  		  		|JABS_boxplot_outlier_figs.pdf ***
	|	  		  		|morpho_boxplot_outlier_figs.pdf ***
	|	  		  		|highest_z_boxplot_outlier_figs.pdf ***
	|	  		  		|most_freq_outliers_boxplot_outlier_figs.pdf ***
	|	  		  		|outlier_videos_summary.csv ***
	|
	|/Nextflow_Output
	|/final_nextflow_feature_data *
	|		 |fecal_boli_raw.csv *
	|		 |fecal_boli_final.csv &&
	|		 |gait_final.csv *
	|		 |JABS_features_final.csv *
	|		 |morphometrics_final.csv *
	|		 |merged_nextflow_dataset.csv **
	|		 
	|Some or all .csv files containing raw nextflow output
	|/CornerCorrection
	|	|Maybe some more .csv files containing raw nextflow output
	|	
	|/Nextflow_cohort1
		|Maybe some more .csv files containing raw nextflow output
		|/CornerCorrection
			|Maybe some more of the .csv files

 * files created by 'NextFlow_Output_QC_Postprocess_1.R'
 ** files created by 'merge_all_dataframes_2.R'
 *** files created by 'final_dataframe_QC_3.R'
 
 && files that must be manually created
```
