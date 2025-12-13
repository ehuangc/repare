from pathlib import Path

from evaluator.comparison_utils import (
    get_mt_colormap,
    get_published_pedigree,
    plot_inferred_pedigree,
    plot_published_pedigree,
    write_relation_differences,
)
from evaluator.pedigree_evaluator import PedigreeEvaluator


def main():
    """
    Compare the Gurgy subset inferred and published pedigrees by plotting and writing relation differences.
    """
    script_dir = Path(__file__).resolve().parent
    data_dir = script_dir / "data" / "gurgy"
    results_dir = script_dir / "results" / "gurgy_subset_comparison"
    results_dir.mkdir(parents=True, exist_ok=True)
    reconstructor_outputs_dir = results_dir / "reconstructor_outputs"
    reconstructor_outputs_dir.mkdir(parents=True, exist_ok=True)

    evaluator = PedigreeEvaluator(
        published_relations_path=data_dir / "published_exact_relations_subset.csv",
        algorithm_nodes_path=data_dir / "nodes_subset.csv",
        algorithm_relations_path=data_dir / "inferred_relations_READv2_subset.csv",
        outputs_dir=reconstructor_outputs_dir,
    )
    inferred_pedigree = evaluator.algorithm_pedigree
    published_pedigree = get_published_pedigree(
        nodes_path=data_dir / "nodes_subset.csv",
        relations_path=data_dir / "published_exact_relations_subset.csv",
    )

    mt_haplogroup_to_color = get_mt_colormap(inferred_pedigree, published_pedigree)
    plot_inferred_pedigree(
        inferred_pedigree,
        plot_path=results_dir / "inferred_pedigree.svg",
        mt_haplogroup_to_color=mt_haplogroup_to_color,
        plot_haplogroups=False,
    )
    # Define dotted edges to add to the published pedigree plot
    nodes_to_remove = [
        "25",  # Placeholder mother of 23 (disconnect GLN291 and GLN288)
        "22",  # Placeholder father of GLN288 (disconnect GLN291 and GLN288)
        "23",  # Placeholder mother of GLN288 (disconnect GLN291 and GLN288)
    ]
    edges_to_remove = [
        ("GLN291", "23"),  # Disconnect GLN291 from GLN288
    ]
    dotted_edges_to_add = [
        ("GLN291", "GLN288"),
    ]
    plot_published_pedigree(
        published_pedigree=published_pedigree,
        plot_path=results_dir / "published_pedigree.svg",
        mt_haplogroup_to_color=mt_haplogroup_to_color,
        nodes_to_remove=nodes_to_remove,
        edges_to_remove=edges_to_remove,
        dotted_edges_to_add=dotted_edges_to_add,
        plot_haplogroups=False,
    )
    write_relation_differences(
        evaluator=evaluator,
        path=results_dir / "relation_differences.csv",
    )


if __name__ == "__main__":
    main()
