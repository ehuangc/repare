import tempfile

import matplotlib.pyplot as plt

from evaluator.pedigree_evaluator import PedigreeEvaluator
from repare.pedigree import Pedigree
from repare.pedigree_reconstructor import PedigreeReconstructor


def get_published_pedigree(nodes_path: str, relations_path: str) -> Pedigree:
    # Write outputs other than the plot to a temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        pedigree_reconstructor = PedigreeReconstructor(
            relations_path, nodes_path, outputs_dir=temp_dir, max_candidate_pedigrees=1000, plot=False
        )
        published_pedigree = pedigree_reconstructor.find_best_pedigree()
    return published_pedigree


def get_mt_colormap(
    inferred_pedigree: Pedigree, published_pedigree: Pedigree
) -> dict[str, tuple[float, float, float, float]]:
    # Build mt_haplogroup color mapping so both plots can use the same colormap
    inferred_pedigree_mt_haplogroups = set(
        [
            inferred_pedigree.get_data(node)["mt_haplogroup"].replace("*", "")
            for node in inferred_pedigree.node_to_data
            if not node.isnumeric()
        ]
    )
    published_pedigree_mt_haplogroups = set(
        [
            published_pedigree.get_data(node)["mt_haplogroup"].replace("*", "")
            for node in published_pedigree.node_to_data
            if not node.isnumeric()
        ]
    )
    mt_haplogroups = sorted(inferred_pedigree_mt_haplogroups | published_pedigree_mt_haplogroups)
    cmap = plt.get_cmap("tab20")
    mt_haplogroup_to_color = {haplogroup: cmap(i / len(mt_haplogroups)) for i, haplogroup in enumerate(mt_haplogroups)}
    return mt_haplogroup_to_color


def plot_inferred_pedigree(inferred_pedigree: Pedigree, plot_path: str, mt_haplogroup_to_color: dict[str, str]) -> None:
    inferred_pedigree.plot(path=plot_path, mt_haplogroup_to_color=mt_haplogroup_to_color)


def plot_published_pedigree(
    published_pedigree: Pedigree,
    plot_path: str,
    mt_haplogroup_to_color: dict[str, str] | None = None,
    nodes_to_remove: list[str] | None = None,
    edges_to_remove: list[tuple[str, str]] | None = None,
    dotted_edges_to_add: list[tuple[str, str]] | None = None,
) -> None:
    published_pedigree.plot(
        path=plot_path,
        mt_haplogroup_to_color=mt_haplogroup_to_color,
        nodes_to_remove=nodes_to_remove,
        edges_to_remove=edges_to_remove,
        dotted_edges_to_add=dotted_edges_to_add,
    )


def write_relation_differences(evaluator: PedigreeEvaluator, path: str):
    evaluator.write_relation_differences(path=path)
