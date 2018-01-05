[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/ASamarkRoth/gammalab_analysis/master?filepath=DataAnalysis_GammaSpectroscopy.ipynb)

# Gamma Lab Analysis

A software framework and jupyter-notebook to support the laboratory exercise KF6 - Gamma Spectroscopy.

## Running the Jupyter-Notebook

You can run it in the web browser on mybinder (without installing anything) by clicking the link [here](https://mybinder.org/v2/gh/ASamarkRoth/gammalab_analysis/master?filepath=DataAnalysis_GammaSpectroscopy.ipynb) (ignore the following in that case). 

It is possible to run the notebook on your local computer as follows:

1. Install [miniconda3](https://conda.io/miniconda.html) alternatively the full [anaconda3](https://www.anaconda.com/download) environment on your laptop (the latter is **much** larger).
2. [Download](https://github.com/mlund/jupyter-course/archive/master.zip) this repository.
3. Install and activate the `ChainCongruency` environment described by the file [`environment.yml`](/environment.yml)  by running the following in a terminal:

```bash
conda env create -f environment.yml
source activate GammaLabAnalysis
./postBuild
```
4. Run the notebook via `jupyter-notebook`

