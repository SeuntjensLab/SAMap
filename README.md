# Styfhals_2022

## Scripts used for data analysis
```
Styfhals, R., Zolotarov, G., Hulselmans, G., Spanier, K. I., Poovathingal, S., Elagoz, A. M., Deryckere, A., Rajewsky, N., Ponte,
G., Fiorito, G. Aerts, S., Seuntjens, E. (2022). Cell type diversity in a developing octopus brain. BioRxiv, 1–34.
```

#### Cell type tree_script.R

Rscript to build a neighbor-joining cell type tree based on the distance matrix of averaged gene expression per cell type.

#### Data_analysis.R

Rscript to analyze, subset, visualize the sn/sc data.

#### DEseq2.R

Code used to analyze Bulk RNA seq data and to identify differentially expressed genes between the optic lobe and central brain.

#### iCytoTRACE.R

Code to assess transcriptional diversity

#### preprocessing_and_clusterstability.R

Pre-processing steps, data integration, cluster stability and clustering of the sn/sc data.

#### tau_octopusBrain_PLD1_snscRNAseq.ipynb

Jupyter notebook to calculate the enrichment of different transcription factor families within the rank of tau (cell type specicity index).

#### tissue specificity.R

Rscript for data preparation for tau calculations.

#### tSNE_RGB_to_CMY.ipynb

Jupyter notebook to change the color scale from RGB to CMY.

#### scopeloomr.R

Script to generate loom files for data visualization in SCope.


## Scripts used for cross-species comparisons with SAMap

For more information: https://github.com/atarashansky/SAMap

Please cite their most recent publication when using SAMap: https://elifesciences.org/articles/66747

#### SAMapMakeH5ADfiles.sh 

#### SAMap_Octopus_Fly

#### Cross_species_analysis.R
Rscript used to visualize cross-species mappings.

### Installation 

#### cloned from github:
```
git clone https://github.com/atarashansky/SAMap.git ${staging}kspan/progs/SAMap/SAMap_v0p1p7/
git checkout d560861
```
#### create conda environment
```
conda create -n SAMapv0.1.6 -c conda-forge python=3.7 pip pybind11 h5py=2.10.0 leidenalg python-igraph texttable
conda activate SAMapv0.1.6
```
#### install SAMap:
```
pip install .
```



