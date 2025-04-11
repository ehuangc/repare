**repare** is a Python package for (ancient) pedigree reconstruction.

## Installation

### Recommended
```
conda create -n "repare" -c conda-forge python=3.13 pygraphviz
conda activate repare
pip install repare
```
repare uses PyGraphviz to plot reconstructed pedigrees. Since PyGraphviz relies on Graphviz which cannot be installed using `pip`, we recommend installing repare and its dependencies in a fresh conda environment.

If you don't need to plot reconstructed pedigrees, you can install repare directly with `pip install repare`. If you need to plot reconstructed pedigrees and have your own Graphviz installation, you can install repare and Pygraphviz with `pip install repare[plot]`.

To install conda, see [this page](https://www.anaconda.com/docs/getting-started/miniconda/install). To install PyGraphviz and Graphviz (yourself), see [this page](https://pygraphviz.github.io/documentation/stable/install.html).
