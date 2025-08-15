:evergreen_tree: **repare** is a Python package for (ancient) pedigree reconstruction.

<p align="center">
  <img src="examples/algorithm_diagram.png" alt="Reconstruction Process Diagram" width="600" />
  <br>
  <em>Diagram of repare's pedigree reconstruction process</em>
</p>

## Installation

### Recommended
```
conda create -n "repare" -c conda-forge python=3.13 pygraphviz matplotlib networkx pandas tqdm
conda activate repare
pip install repare
```
repare uses PyGraphviz to plot reconstructed pedigrees. Since PyGraphviz relies on Graphviz which cannot be installed using `pip`, we recommend installing repare and its dependencies in a fresh conda environment, as shown above.

If you don't need to plot reconstructed pedigrees, you can install repare directly with `pip install repare`. If you need to plot reconstructed pedigrees and have your own Graphviz installation, you can install repare and Pygraphviz with `pip install repare[plot]`.

To install conda, see [this page](https://www.anaconda.com/docs/getting-started/miniconda/install). To install PyGraphviz and Graphviz (yourself), see [this page](https://pygraphviz.github.io/documentation/stable/install.html).


## Usage

We recommend running repare through its command-line interface.
```
repare -n NODES -r RELATIONS [-o OUTPUT] [-m MAX_CANDIDATE_PEDIGREES] [-e EPSILON] [-s SEED] [-d] [-w] [-v]
```

> [!NOTE]
> Minimal command:
> ```
> repare -n nodes.csv -r relations.csv
> ```
> For example data inputs, see [examples/nodes.csv](examples/nodes.csv) and [examples/relations.csv](examples/relations.csv).

### Inputs
**Nodes** (-n) (*required*): Path to a CSV file that contains information about the individuals to be analyzed by repare. 

<blockquote><details open>
  <summary><ins>Nodes CSV file columns</ins></summary>

  - **id** *(required)*: ID of individual. Cannot be fully numeric, as numeric IDs are reserved for placeholder nodes.
  - **sex** *(required)*: Genetic sex of individual.
  - **y_haplogroup** *(required)*: Y chromosome haplogroup of individual. Can include "*" as a wildcard expansion character at the end if haplogroup is not fully inferred.
  - **mt_haplogroup** *(required)*: Mitochondrial haplogroup of individual. Can include "*" as a wildcard expansion character at the end if haplogroup is not fully inferred.
  - **can_have_children** *(optional)*: Whether the individual *can* have offspring (e.g., as indicated by age of death). Defaults to "True".
  - **can_be_inbred** *(optional)*: Whether the individual *can* have parents related at the 3rd-degree or closer (e.g., as indicated by ROH). Defaults to "True".
  - **years_before_present** *(optional)*: (Approximate) date of birth of individual, in years before present. If provided, will be used to prune temporally invalid pedigrees. *This column should only be used when backed by strong dating evidence.*
</details></blockquote>

**Relations** (-r) (*required*): Path to a CSV file that contains information about inferred pairwise kinship relations. All individuals included in this file must be specified in the nodes CSV.

<blockquote><details open>
  <summary><ins>Relations CSV file columns</ins></summary>

  - **id1** *(required)*: ID of individual 1.
  - **id2** *(required)*: ID of individual 2.
  - **degree** *(required)*: Degree of (inferred) kinship relation between individual 1 and individual 2. Must be "1", "2", or "3". Higher-degree relatives are considered unrelated.
  - **constraints** *(optional)*: Semicolon-delimited list of possible configurations of kinship relation. For example, a parental 1st-degree relation can be constrained with "parent-child;child-parent". Many kinship inference methods will classify 1st-degree relation types, which can be used as relation constraints.
</details></blockquote>

**Output** (-o) (*optional*): Path to directory for saving repare outputs. Defaults to the current working directory.

**Max Candidate Pedigrees** (-m) (*optional*): Maximum number of candidate pedigrees to keep after each algorithm iteration. Defaults to 1000.

**Epsilon** (-e) (*optional*): Parameter for adapted epsilon-greedy sampling at the end of each algorithm iteration. Defaults to 0.2.

**Seed** (-s) (*optional*): Random seed for reproducibility. Defaults to 42.

**Do Not Plot** (-d) (*flag*): If set, do not plot reconstructed pedigree(s).

**Write Alternate Pedigrees** (-w) (*flag*): If set, write outputs for alternate reconstructed pedigrees to disk.

**Verbose** (-v) (*flag*): If set, enable verbose output (INFO-level logging).

## Reproducibility
We recommend using [pixi](https://pixi.sh/) to reproduce the results in this repo.
```
git clone https://github.com/ehuangc/repare.git
cd repare
pixi shell
```

Once in the pixi shell, you can run the script(s) corresponding to the results you'd like to reproduce. For example:
```
python benchmarks/published/run_parameter_experiment.py
exit
```
To install pixi, see [this page](https://pixi.sh/latest/installation/).
