<a href="https://mybinder.org/v2/gh/ASamarkRoth/gammalab_analysis/master?filepath=DataAnalysis_GammaSpectroscopy.ipynb" target="_blank" title="Run jupyter-notebook in binder.">
  <img src="https://mybinder.org/badge.svg">
</a>

# Gamma Lab Analysis

A software framework and jupyter-notebook to support the laboratory exercise KF6 - Gamma Spectroscopy.

If you find any errors or if you have any questions, make a pull request, post an issue or write me an [email](mailto:anton.samark-roth@nuclear.lu.se). 

## Running the Jupyter-Notebook

You can run it in the web browser on mybinder (without installing anything) by clicking the link [here](https://mybinder.org/v2/gh/ASamarkRoth/gammalab_analysis/master?filepath=DataAnalysis_GammaSpectroscopy.ipynb) (ignore the following in that case). However, it is recommended to instead install it locally.

It is possible to run the notebook on your local computer as follows:

1. Install [miniconda3](https://conda.io/miniconda.html) alternatively the full [anaconda3](https://www.anaconda.com/download) environment on your laptop (the latter is **much** larger).
2. [Download](https://github.com/ASamarkRoth/gammalab_analysis/archive/master.zip) this repository.
3. Install and activate the `GammaLabAnalysis` environment described by the file [`environment.yml`](/environment.yml)  by running the following in a terminal:

```bash
conda env create -f environment.yml
source activate GammaLabAnalysis
./postBuild
```
4. Run the notebook via `jupyter-notebook` or if you prefer with `jupyter-lab`.

If you want to try out jupyter-lab you can do so here: <a href="https://mybinder.org/v2/gh/ASamarkRoth/gammalab_analysis/master?urlpath=lab" target="_blank">jupyter-lab binder</a>
