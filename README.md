# JABS_nextflow_postprocess

This repository is R code developed by Jake to take raw nextflow outputs from the JABS pipeline and conduct a series of qc screening, data visualization, and post processing to prepare the data for its final analysis. It currently consists of 5 scripts that should be run in order, denoted by the number at the suffix of the file name. They are described below.

## 





## Expected file structure from scripts
#3. Below is an example of how data might be organized and still work with this code:
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
			|	  |		 |qc_all.csv *
			|	  |		 |qc_failed.csv *
			|	  |
			|	  |/missing_or_dup_data
			|	  |		 |missing_data.csv *
			|	  |		 |videos_not_in_qc_report *
			|	  |		 |duplicated_data.xlsx *
			|	  |		 |NetworkFilenames_missing_in_data.csv **
			|	  |		 |mice_missing_in_metadata.csv **
			|	  |
			|	  |/qc_figs
			|	  		 |fecal_boli_qc_figs.pdf *
			|	  		 |scatter_plot_lm_figs.pdf **
			|	  		 |
			|	  		 |/zscore_boxplots
			|	  		  		|gait_boxplot_outlier_figs.pdf **
			|	  		  		|JABS_boxplot_outlier_figs.pdf **
			|	  		  		|morpho_boxplot_outlier_figs.pdf **
			|	  		  		|highest_z_boxplot_outlier_figs.pdf **
			|	  		  		|most_freq_outliers_boxplot_outlier_figs.pdf **
			|	  		  		|outlier_videos_summary.csv **
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
                 |Maybe some more .csv files containing raw nextflow output
                 |/Nextflow_cohort1
                     |Maybe some more .csv files containing raw nextflow output
                     |/CornerCorrection
                         |Maybe some more of the .csv files

 * files created by 'NextFlow_Output_QC_Postprocess_1.R'
 ** files created by 'merge_all_dataframes_2.R'
 *** files created by 'final_dataframe_QC_3.R'
 
 && files that must be manually created








QC plot ideas:

