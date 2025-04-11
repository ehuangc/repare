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


## Usage
We recommend using repare as a command-line interface.
```
repare -n NODES -r RELATIONS [-o OUTPUT] [-m MAX_CANDIDATE_PEDIGREES] [-e EPSILON] [-s SEED] [-d] [-w] [-v]
```

**Nodes (-n, required)**: path to a CSV file that contains information about the individuals to be analyzed by repare.

**Relations (-r, required)**: path to a CSV file that contains information about inferred pairwise kinship relations. All individuals included in this file must be specified in the nodes CSV.

**Output (-o, optional)**: path to directory for saving repare outputs. Defaults to the current working directory.

**Max Candidate Pedigrees (-m, optional)**: maximum number of candidate pedigrees to keep after each algorithm iteration. Defaults to 1000.

**Epsilon (-e, optional)**: parameter for adapted epsilon-greedy sampling at the end of each algorithm iteration. Defaults to 0.2.

**Seed (-s, optional)**: random seed for reproducibility. Defaults to 42.

**Do Not Plot (-d, flag)**: do not plot reconstructed pedigree(s).

**Write Alternate Pedigrees (-w, flag)**: write outputs for alternate reconstructed pedigrees to disk.

**Verbose (-v, flag)**: enable verbose output (INFO-level logging).
