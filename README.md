# Styfhals_2022
## Scripts used for cross-species comparisons with SAMap
For more information: https://github.com/atarashansky/SAMap

Please cite their most recent publication when using SAMap: https://elifesciences.org/articles/66747

## Installation 

### cloned from github:
```
git clone https://github.com/atarashansky/SAMap.git ${staging}kspan/progs/SAMap/SAMap_v0p1p7/
git checkout d560861
```
### create conda environment
```
conda create -n SAMapv0.1.6 -c conda-forge python=3.7 pip pybind11 h5py=2.10.0 leidenalg python-igraph texttable
conda activate SAMapv0.1.6
```
### install SAMap:
```
pip install .
```
