For pose estimation data go to https://zenodo.org/deposit/8186065, download the data.zip file and extract its contents to this directory. This adds multiple folders that are named after each experiment. Each folder contains the pose estimation data from DeepLabCut for the corresponding experiment. For your convenience we reposited the pre-processed USData objects that we created from the raw data here .These do not containe any pose esimation data, but only label data(i.e the cluster number obtained with kmeans) at each frame for each recording. they are much smaller than the raw data but can be loaded quickly and analyzed with the code of the BehaviorFlow package (see Tutorial)

**US_AllData_25Clusters.rds**<br><br>
This file contains all the kmeans 25 results on a per-frame basis for the CSI, AS, Yohimbine, DREADD and CRS experiments. it is already pre-processed and annotated with metadata

**US_CSI_SensitivityAssay_processed.rds**<br><br>
This file contains all the kmeans (with 10,25,50,70,100 centers) as well as BSOiD and VAME clustering results on a per-frame basis for the CSI dataset. it is already pre-processed and annotated with metadata

**US_IFS_25Clusters_withResilience.rds**<br><br>
This file contains all the kmeans 25 results on a per-frame basis for the IFS (=traumatic footshock, TFS in the manuscript)  dataset. it is already pre-processed and annotated with metadata

**US_YohimbineRoche_25Clusters_processed.rds**<br><br>
This file contains all the kmeans 25 results on a per-frame basis for Yohimbine dataset produced by an external lab (Roche). it is already pre-processed and annotated with metadata

**BSOiD**<br><br>
This folder contains the clustering results for the CSI dataset as obtained by the BSOiD algorithm

**VAME**<br><br>
This folder contains the clustering results for the CSI dataset as obtained by the VAME algorithm
