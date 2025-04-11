import logging
import os

from comparator.relation_comparison import RelationComparison
from tqdm.contrib.logging import logging_redirect_tqdm


def main():
    data_dir = os.path.join(os.path.dirname(__file__), "data", "fowler")
    published_relations_path = os.path.join(data_dir, "published_exact_relations.csv")
    algorithm_nodes_path = os.path.join(data_dir, "nodes.csv")
    algorithm_relations_path = os.path.join(data_dir, "inferred_relations_coeffs.csv")

    logging.basicConfig(level=logging.WARNING)  # Set to logging.INFO for more detailed output
    with logging_redirect_tqdm():
        relation_comparison = RelationComparison(
            published_relations_path=published_relations_path,
            algorithm_nodes_path=algorithm_nodes_path,
            algorithm_relations_path=algorithm_relations_path,
        )
        for relation, count in relation_comparison.get_metrics().items():
            print(f"{relation}: {round(count, 2)}")


if __name__ == "__main__":
    main()
