import logging
from pathlib import Path

from tqdm.contrib.logging import logging_redirect_tqdm

from repare.pedigree_reconstructor import PedigreeReconstructor


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    data_dir = script_dir / "data" / "koszyce"
    nodes_path = data_dir / "nodes.csv"
    relations_path = data_dir / "inferred_relations_custom.csv"
    outputs_dir = script_dir / "results" / "koszyce_reconstruction"
    outputs_dir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(level=logging.WARNING)
    with logging_redirect_tqdm():
        print("Beginning pedigree reconstruction for Koszyce...")
        pedigree_reconstructor = PedigreeReconstructor(
            relations_path=relations_path,
            nodes_path=nodes_path,
            outputs_dir=outputs_dir,
            plot_haplogroups=False,
        )
        pedigree_reconstructor.find_best_pedigree()

    print(f"Finished pedigree reconstruction. Outputs written to {outputs_dir.resolve()}.")


if __name__ == "__main__":
    main()
