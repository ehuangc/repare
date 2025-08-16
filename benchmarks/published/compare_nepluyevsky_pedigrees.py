import os

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
    Compare the Nepluyevsky inferred and published pedigrees by plotting and writing relation differences.
    """
    data_dir = os.path.join(os.path.dirname(__file__), "data", "nepluyevsky")
    results_dir = os.path.join(os.path.dirname(__file__), "results", "nepluyevsky_comparison")
    os.makedirs(results_dir, exist_ok=True)

    evaluator = PedigreeEvaluator(
        published_relations_path=os.path.join(data_dir, "published_exact_relations.csv"),
        algorithm_nodes_path=os.path.join(data_dir, "nodes.csv"),
        algorithm_relations_path=os.path.join(data_dir, "inferred_relations_KIN.csv"),
    )
    inferred_pedigree = evaluator.algorithm_pedigree
    published_pedigree = get_published_pedigree(
        nodes_path=os.path.join(data_dir, "nodes.csv"),
        relations_path=os.path.join(data_dir, "published_exact_relations.csv"),
    )

    mt_haplogroup_to_color = get_mt_colormap(inferred_pedigree, published_pedigree)
    plot_inferred_pedigree(
        inferred_pedigree,
        # Plot as SVG because we will need to crop and shift nodes
        plot_path=os.path.join(results_dir, "inferred_pedigree.svg"),
        mt_haplogroup_to_color=mt_haplogroup_to_color,
    )
    plot_published_pedigree(
        published_pedigree=published_pedigree,
        # Plot as SVG because we will need to crop and shift nodes
        plot_path=os.path.join(results_dir, "published_pedigree.svg"),
        mt_haplogroup_to_color=mt_haplogroup_to_color,
    )
    write_relation_differences(
        evaluator=evaluator,
        path=os.path.join(results_dir, "relation_differences.csv"),
    )


if __name__ == "__main__":
    main()
