import os
import argparse
import logging
from tqdm.contrib.logging import logging_redirect_tqdm
from . import PedigreeEnsemble


def parse_arguments():
    parser = argparse.ArgumentParser(description="Build and analyze pedigrees.")
    parser.add_argument("-n", "--nodes", type=str, required=True, help="Path to the nodes CSV file.")
    parser.add_argument("-r", "--relations", type=str, required=True, help="Path to the relations CSV file.")
    parser.add_argument("-o", "--output", type=str, default=".", help="Directory to save the output files. Defaults to the current directory.")
    parser.add_argument("-s", "--sample_count", type=int, default=1000, help="Number of pedigrees to keep at each relation step. Default is 100.")
    parser.add_argument("-e", "--epsilon", type=float, default=0.2, help="Epsilon value for the simulation. Default is 0.2.")
    parser.add_argument("-d", "--seed", type=int, default=42, help="Random seed for reproducibility. Default is 42.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output (INFO-level logging).")
    return parser.parse_args()

def main():
    args = parse_arguments()
    logging_level = logging.INFO if args.verbose else logging.WARNING
    logging.basicConfig(level=logging_level, format="%(levelname)s - %(message)s")

    output_dir = args.output
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    with logging_redirect_tqdm():
        pedigree_ensemble = PedigreeEnsemble(
            relations_path=args.relations,
            nodes_path=args.nodes,
            outputs_dir=output_dir,
            sample_count=args.sample_count,
            epsilon=args.epsilon,
            random_seed=args.seed
        )
        pedigree = pedigree_ensemble.find_best_pedigree()
        pedigree.plot(path=os.path.join(output_dir, "pedigree.png"))

if __name__ == "__main__":
    main()
