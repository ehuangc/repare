import os
import tempfile

import matplotlib.pyplot as plt
from comparator.relation_comparison import RelationComparison

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


def write_relation_differences(relation_comparison: RelationComparison, path: str):
    relation_comparison.write_relation_differences(path=path)


def main():
    """
    Compare the Gurgy inferred and published pedigrees by plotting and writing relation differences.
    """
    data_dir = os.path.join(os.path.dirname(__file__), "data", "gurgy")
    results_dir = os.path.join(os.path.dirname(__file__), "results", "gurgy_analysis")

    relation_comparison = RelationComparison(
        published_relations_path=os.path.join(data_dir, "published_exact_relations.csv"),
        algorithm_nodes_path=os.path.join(data_dir, "nodes.csv"),
        algorithm_relations_path=os.path.join(data_dir, "inferred_relations_READv2.csv"),
    )
    inferred_pedigree = relation_comparison.algorithm_pedigree
    published_pedigree = get_published_pedigree(
        nodes_path=os.path.join(data_dir, "nodes.csv"),
        relations_path=os.path.join(data_dir, "published_exact_relations.csv"),
    )

    # Build mt_haplogroup color mapping so both plots use the same colormap
    inferred_pedigree_mt_haplogroups = set(
        [
            inferred_pedigree.node_to_data[node]["mt_haplogroup"].replace("*", "")
            for node in inferred_pedigree.node_to_data
            if not node.isnumeric()
        ]
    )
    published_pedigree_mt_haplogroups = set(
        [
            published_pedigree.node_to_data[node]["mt_haplogroup"].replace("*", "")
            for node in published_pedigree.node_to_data
            if not node.isnumeric()
        ]
    )
    mt_haplogroups = sorted(inferred_pedigree_mt_haplogroups | published_pedigree_mt_haplogroups)
    cmap = plt.get_cmap("tab20")
    mt_haplogroup_to_color = {haplogroup: cmap(i / len(mt_haplogroups)) for i, haplogroup in enumerate(mt_haplogroups)}

    plot_inferred_pedigree(
        inferred_pedigree,
        plot_path=os.path.join(results_dir, "inferred_pedigree.png"),
        mt_haplogroup_to_color=mt_haplogroup_to_color,
    )
    # Define dotted edges to add to the published pedigree plot
    nodes_to_remove = [
        "127",  # Placeholder mother of 27  (disconnect GLN276 and GLN215A/GLN215B)
        "26",  # Placeholder father of GLN215A and GLN215B (disconnect GLN276 and GLN215A/GLN215B)
        "27",  # Placeholder mother of GLN215A and GLN215B (disconnect GLN276 and GLN215A/GLN215B)
        "173",  # Placeholder father of 170 (disconnect GLN319 and GLN232C)
        "170",  # Placeholder father of GLN232C (disconnect GLN319 and GLN232C)
        "171",  # Placeholder mother of GLN232C (disconnect GLN319 and GLN232C)
        "189",  # Placeholder mother of 187 (disconnect GLN291 and GLN288)
        "186",  # Placeholder father of GLN288 (disconnect GLN291 and GLN288)
        "187",  # Placeholder mother of GLN288 (disconnect GLN291 and GLN288)
    ]
    edges_to_remove = [
        ("124", "122"),  # Disconnect GLN320 from GLN231A and GLN270B
        ("125", "122"),  # Disconnect GLN320 from GLN231A and GLN270B
        ("GLN319", "170"),  # Disconnect GLN319 from GLN232C
        ("GLN291", "187"),  # Disconnect GLN291 from GLN288
    ]
    dotted_edges_to_add = [
        ("GLN320", "GLN231A"),
        ("GLN320", "GLN270B"),
        ("GLN276", "GLN215A"),
        ("GLN276", "GLN215B"),
        ("GLN319", "GLN232C"),
        ("GLN291", "GLN288"),
    ]
    plot_published_pedigree(
        published_pedigree=published_pedigree,
        plot_path=os.path.join(results_dir, "published_pedigree.png"),
        mt_haplogroup_to_color=mt_haplogroup_to_color,
        nodes_to_remove=nodes_to_remove,
        edges_to_remove=edges_to_remove,
        dotted_edges_to_add=dotted_edges_to_add,
    )
    write_relation_differences(
        relation_comparison=relation_comparison,
        path=os.path.join(results_dir, "relation_differences.csv"),
    )


if __name__ == "__main__":
    main()
