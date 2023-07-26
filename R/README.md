This folder contains all the scripts used for the original publication (LINK), as well as the code collection required for BehaviorFlow. below you find more information about the scripts and what order they were used with. keep in mind that you will need to download the raw data from zenodo and place it correctly if you want to recreate these steps. you will find instructions for this on the main tutorial file<br><br>

**UnsupervisedAnalysis_Functions.R** <br>
Code package for all BehaviorFlow functions, including BFA and BFF. if you want to analyze your own clustering data with behavior flow you will only need this file

**DLCAnalyzer_Functions_v4.R** <br>
An updated version of the DLCAnalyzer (see Tutorial) code. we recommend you use this for any DLCAnalyzer analyses going forward since it is more integrated with the downstream applications in the BehaviorFlow package

**KmeansClustering_CSI_70Clusters.Rmd** <br>
script that explains how the clustering for Figure 1 of the manuscript was performed from pose estimation data. this only covered the CSI dataset with 70 kmeans clusters

**SensitivityAssays_CSI.Rmd** <br>
Script that contains the code to produce the data for the sensitivity curves (in silico reduced group size analysis) for Figure 1

**Comparison_ClusteringAlgorithms.Rmd** <br>
Script that contains the code to produce the data for the comparison between kmeans, BSOiD and VAME. important for Figure 1

**KmeansClustering_AllData_25Clusters.Rmd** <br>
Code that performs a big clustering across 3 experiments (for each we take 20 random files). it then runs multiple kmeans clusterings across these 60 files with a varying number of centers (10,25,50,100) once with dimensionality reduction prior (svd), and once without. next, for each label group it uses a neural network to train a clustering classifier. Note, for the manuscript we went ahead with 25 clusters later and no dimensionality reduction!, so most steps can be omitted to reprodcue the mansucript data. important for Figures 2,3,4,5

**ClusterClassifier_AllData_25Clusters.Rmd** <br>
Code that uses the clustering classifiers for kmeans 25 trained with the script KmeansClustering_AllData_25Clusters.Rmd on all data and with this generates clustering results for all data across all experiments. important for Figures 2,3,4,5

**SensitivityAssays_AllData_25Clusters.Rmd** <br>
Script that generates sensitivity curves data for all kmeans 25 BFA results across all experiments. important for Figures 2,3,4,5

**ClusteringAnalysis_IFS_25Clusters.Rmd** <br>
Code that uses the kmeans 25 classifier (see KmeansClustering_AllData_25Clusters.Rmd) on the IFS (=traumatic footshock or TFS in the manuscript) data. important for Figure 5

**ResiliencePhenotype_IFS.Rmd** <br>
Code that describes how we use the 2D embedding generated in ClusteringAnalysis_IFS_25Clusters.Rmd to determine which animals in the TFS groups are considered responders vs non-responders. important for Figure 5

**KmeansClustering_YohimbineRoche_25Clusters.Rmd** <br>
Code that describes how we use the external data from Roche to perform an independent clustering with 25 kmeans clusters on their data. important for Figure 6

**Figures Folder**
Folder that contains the scripts that were used to create the final figures for the manuscript
