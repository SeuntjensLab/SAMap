# Styfhals_2022

All code was used to generate and analyze the data presented in the following paper:

```
Styfhals, R., Zolotarov, G., Hulselmans, G., Spanier, K. I., Poovathingal, S., Elagoz, A. M., Deryckere, A., Rajewsky, N., Ponte,
G., Fiorito, G. Aerts, S., Seuntjens, E. (2022). Cell type diversity in a developing octopus brain. BioRxiv, 1–34.
```

## Scripts used for data analysis

#### Cell type tree_script.R

Rscript to build a neighbor-joining cell type tree based on the distance matrix of averaged gene expression per cell type.

#### Data_analysis.R

Rscript to analyze, subset, visualize the sn/sc data.

#### DEseq2.R

Code used to analyze Bulk RNA seq data and to identify differentially expressed genes between the optic lobe and central brain.

#### iCytoTRACE.R

Code to assess transcriptional diversity.

#### preprocessing_and_clusterstability.R

Pre-processing steps, data integration, cluster stability and clustering of the sn/sc data.

#### tau_octopusBrain_PLD1_snscRNAseq.ipynb

Jupyter notebook to calculate the enrichment of different transcription factor families within the rank of tau (cell type specicity index).

#### tissue specificity.R

Rscript for data preparation for tau calculations.

#### gene_set_enrichment_test_Fisher.R

Rscript for a two-sided Fisher's Exact Test to analyze the overrepresentation of certain gene sets versus the background. 
Contingency tables should be constructed accordingly so that Odds Ratio's < 1 indicate enrichment, Odds Ratio's > indicate an underrepresentation.  

```
matrix(c(total number of other genes in cell type of interest, total number of genes of interest in cell type of interest, total number of other genes in all other cell types, total number of genes of interest in all other cell types), nrow=2)
```

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



