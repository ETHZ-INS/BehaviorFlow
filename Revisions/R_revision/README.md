This folder contains all the scripts used for the revised manuscript as well as the updated code collection required for BehaviorFlow analysis. Below, you find more information about the scripts. Keep in mind that you need to download the pose estimation files from zenodo and adjust the path to these files accordingly in the corresponding scripts to reproduce all steps. You find more instructions in the main tutorial file. <br><br>

**BFA_BFF.R** <br>
Code package for all BehaviorFlow functions, including BFA, BFL and BFF. This file replaces the file "UnsupervisedAnalysis_Functions.R". You only need this file if you want to analyze your own clustering data with BehaviorFlow.

**DLCAnalyzer.R** <br>
An updated version of the DLCAnalyzer (see Tutorial) package, replacing the "DLCAnalyzer_Functions_v4.R". We recommend to use this file for any DLCAnalyzer analyses going forward since it is more integrated with the downstream applications in the BehaviorFlow package.

**VAME_Clustering_AllData.Rmd** <br>
Script containing the code to train a cluster classifier imitating the clustering produced by VAME on 5 datasets (CSI, AS, yohimbine, CRS, DREADD). Important for Suppl. Figure 5.

**VAME_Analysis_25Clusters.Rmd** <br>
Follow up script of the previous one, containing all analysis steps for the VAME clustering of several datasets. Findings are shown in Suppl. Figure 5.

**MBTY_Clustering.Rmd** <br>
Script containing the code used to cluster the marble burying test recordings after yohimbine injections (MBTY). The results of this are shown in the updated Figure 6.

**MBTY_Analysis.Rmd** <br>
Follow up script of the previous one, containing all the analysis of the MBTY experiment for Figure 6.

**Roche_Clustering_Diaz_Yohimbine.Rmd** <br>
Script containing the code to cluster the experimental data received from Roche. Important for Figure 6 and Suppl. Figure 7.

**Roche_Analysis_Diaz_Yohimbine.Rmd** <br>
Follow up script of the previous one, including all analysis steps of the yohimbine and diazepam experiments presented in Figure 6 and Suppl. Figure 7.

**CRS1_LDB_Clustering.Rmd** <br>
Script containing the code to cluster the light-dark box (LDB) recordings after chronic restraint stress. Important for Suppl. Figure 7. 

**CRS1_LDB_Analysis.Rmd** <br>
Follow up script of the previous one containing all analysis steps of the LDB data. Results are shown in Suppl. Figure 7.

**Clustering_IFS_Shockbox.Rmd** <br>
Code that clusters the top down recordings from the fear conditioning box during the inescapable footshock (IFS) session and during the 6 extinction sessions. Important for Suppl. Figure 7.

**Analysis_IFS_Shockbox.Rmd** <br>
Follow up script of the previous one analysing the fear conditioning recordings. Results are shown in Suppl. Figure 7.

**Clustering_diffClusters_MBT_LDB_IFS_Roche.Rmd** <br>
Script in which we added clusterings for different numbers of clusters (10, 50, 70, 100) for all the new datasets (MBTY, Roche data, LDB, fear conditioning box). These clusterings were used for the power analyses, shown in Figure 6 and Suppl. Figure 7.  

**Figure Folder**
Folder containing the scripts that were used to create the revised and new figures for the updated manuscript.
