import logging
from pathlib import Path

from evaluator.pedigree_evaluator import PedigreeEvaluator
from tqdm.contrib.logging import logging_redirect_tqdm


def main():
    script_dir = Path(__file__).resolve().parent
    evaluation_outputs_dir = script_dir / "results" / "evaluation_outputs"
    evaluation_outputs_dir.mkdir(parents=True, exist_ok=True)
    results_path = evaluation_outputs_dir / "results.csv"

    for idx, (site, relations_file_name) in enumerate(
        [
            ("hazleton_north", "inferred_relations_coeffs.csv"),
            ("hazleton_north", "inferred_relations_custom.csv"),
            ("nepluyevsky", "inferred_relations_KIN.csv"),
            ("nepluyevsky", "inferred_relations_custom.csv"),
            ("koszyce", "inferred_relations_custom.csv"),
            ("gurgy", "inferred_relations_READv2.csv"),
        ]
    ):
        print(f"Reconstructing pedigree: site={site}, relation_data={relations_file_name}")
        data_dir = script_dir / "data" / site
        algorithm_nodes_path = data_dir / "nodes.csv"
        algorithm_relations_path = data_dir / relations_file_name
        published_relations_path = data_dir / "published_exact_relations.csv"
        outputs_dir = evaluation_outputs_dir / site / Path(relations_file_name).stem
        outputs_dir.mkdir(parents=True, exist_ok=True)

        logging.basicConfig(level=logging.WARNING)  # Set to logging.INFO for more detailed output
        with logging_redirect_tqdm():
            evaluator = PedigreeEvaluator(
                published_relations_path=published_relations_path,
                algorithm_nodes_path=algorithm_nodes_path,
                algorithm_relations_path=algorithm_relations_path,
                outputs_dir=outputs_dir,
            )

            with results_path.open("a") as file:
                metrics_values = evaluator.get_metrics()
                metrics = list(metrics_values.keys())
                values = list(metrics_values.values())
                if idx == 0:
                    file.truncate(0)
                    file.write(",".join(["site", "relation_data"] + list(metrics)) + "\n")
                file.write(f"{site},{relations_file_name}," + ",".join([str(value) for value in values]) + "\n")


if __name__ == "__main__":
    main()
