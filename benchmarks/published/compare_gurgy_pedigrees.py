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
    Compare the Gurgy inferred and published pedigrees by plotting and writing relation differences.
    """
    data_dir = os.path.join(os.path.dirname(__file__), "data", "gurgy")
    results_dir = os.path.join(os.path.dirname(__file__), "results", "gurgy_comparison")
    os.makedirs(results_dir, exist_ok=True)

    evaluator = PedigreeEvaluator(
        published_relations_path=os.path.join(data_dir, "published_exact_relations.csv"),
        algorithm_nodes_path=os.path.join(data_dir, "nodes.csv"),
        algorithm_relations_path=os.path.join(data_dir, "inferred_relations_READv2.csv"),
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
        # Plot as SVG because we will need to crop and shift nodes
        plot_path=os.path.join(results_dir, "published_pedigree.svg"),
        mt_haplogroup_to_color=mt_haplogroup_to_color,
        nodes_to_remove=nodes_to_remove,
        edges_to_remove=edges_to_remove,
        dotted_edges_to_add=dotted_edges_to_add,
    )
    write_relation_differences(
        evaluator=evaluator,
        path=os.path.join(results_dir, "relation_differences.csv"),
    )


if __name__ == "__main__":
    main()
