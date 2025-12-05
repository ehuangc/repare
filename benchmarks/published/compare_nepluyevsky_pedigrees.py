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
    Compare the Nepluyevsky inferred and published pedigrees by plotting and writing relation differences.
    """
    script_dir = Path(__file__).resolve().parent
    data_dir = script_dir / "data" / "nepluyevsky"
    results_dir = script_dir / "results" / "nepluyevsky_comparison"
    results_dir.mkdir(parents=True, exist_ok=True)

    evaluator = PedigreeEvaluator(
        published_relations_path=data_dir / "published_exact_relations.csv",
        algorithm_nodes_path=data_dir / "nodes.csv",
        algorithm_relations_path=data_dir / "inferred_relations_KIN.csv",
    )
    inferred_pedigree = evaluator.algorithm_pedigree
    published_pedigree = get_published_pedigree(
        nodes_path=data_dir / "nodes.csv",
        relations_path=data_dir / "published_exact_relations.csv",
    )

    mt_haplogroup_to_color = get_mt_colormap(inferred_pedigree, published_pedigree)
    plot_inferred_pedigree(
        inferred_pedigree,
        # Plot as SVG because we will need to crop and shift nodes
        plot_path=results_dir / "inferred_pedigree.svg",
        mt_haplogroup_to_color=mt_haplogroup_to_color,
    )
    plot_published_pedigree(
        published_pedigree=published_pedigree,
        # Plot as SVG because we will need to crop and shift nodes
        plot_path=results_dir / "published_pedigree.svg",
        mt_haplogroup_to_color=mt_haplogroup_to_color,
    )
    write_relation_differences(
        evaluator=evaluator,
        path=results_dir / "relation_differences.csv",
    )


if __name__ == "__main__":
    main()
