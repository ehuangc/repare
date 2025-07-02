import logging
import os

from evaluator.pedigree_evaluator import PedigreeEvaluator
from tqdm.contrib.logging import logging_redirect_tqdm


def main():
    results_dir = os.path.join(os.path.dirname(__file__), "results")
    os.makedirs(results_dir, exist_ok=True)
    results_path = os.path.join(results_dir, "results.csv")

    for idx, (site, relations_file_name) in enumerate(
        [
            ("hazleton_north", "inferred_relations_coeffs.csv"),
            ("hazleton_north", "inferred_relations_custom.csv"),
            ("nepluyevsky", "inferred_relations_KIN.csv"),
            ("nepluyevsky", "inferred_relations_custom.csv"),
            ("gurgy", "inferred_relations_READv2.csv"),
        ]
    ):
        print(f"Reconstructing pedigree: site={site}, relation_data={relations_file_name}")
        data_dir = os.path.join(os.path.dirname(__file__), "data", site)
        algorithm_nodes_path = os.path.join(data_dir, "nodes.csv")
        algorithm_relations_path = os.path.join(data_dir, relations_file_name)
        published_relations_path = os.path.join(data_dir, "published_exact_relations.csv")

        logging.basicConfig(level=logging.WARNING)  # Set to logging.INFO for more detailed output
        with logging_redirect_tqdm():
            evaluator = PedigreeEvaluator(
                published_relations_path=published_relations_path,
                algorithm_nodes_path=algorithm_nodes_path,
                algorithm_relations_path=algorithm_relations_path,
            )

            with open(results_path, "a") as file:
                metrics_values = evaluator.get_metrics()
                metrics = list(metrics_values.keys())
                values = list(metrics_values.values())
                if idx == 0:
                    file.truncate(0)
                    file.write(",".join(["site", "relation_data"] + list(metrics)) + "\n")
                file.write(f"{site},{relations_file_name}," + ",".join([str(value) for value in values]) + "\n")


if __name__ == "__main__":
    main()
