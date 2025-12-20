import argparse
import logging
from datetime import datetime
from pathlib import Path

from tqdm.contrib.logging import logging_redirect_tqdm

from repare.pedigree_reconstructor import PedigreeReconstructor


def parse_arguments():
    base_parser = argparse.ArgumentParser(add_help=False)
    base_parser.add_argument(
        "--print-allowed-constraints",
        action="store_true",
        help="Print the list of allowed constraint strings and exit.",
    )

    # Handle the informational flag before enforcing required arguments
    base_args, remaining_args = base_parser.parse_known_args()
    if base_args.print_allowed_constraints:
        for constraint in sorted(PedigreeReconstructor.get_allowed_constraints()):
            print(constraint)
        base_parser.exit()

    parser = argparse.ArgumentParser(
        parents=[base_parser],
        description="Reconstruct (ancient) pedigrees from pairwise kinship relations.",
    )
    parser.add_argument("-n", "--nodes", type=str, required=True, help="Path to the nodes CSV file.")
    parser.add_argument("-r", "--relations", type=str, required=True, help="Path to the relations CSV file.")
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default=datetime.now().strftime("repare_output_%Y%m%d_%H%M%S"),
        help=(
            "Directory to save the output files. Defaults to repare_output_<YYYYMMDD_HHMMSS> in the current directory."
        ),
    )
    parser.add_argument(
        "-m",
        "--max_candidate_pedigrees",
        type=int,
        default=1000,
        help="Number of pedigrees to keep after each iteration. Defaults to 1000.",
    )
    parser.add_argument(
        "-e",
        "--epsilon",
        type=float,
        default=0.2,
        help="Parameter for adapted epsilon-greedy sampling at the end of each algorithm iteration. Defaults to 0.2.",
    )
    parser.add_argument("-s", "--seed", type=int, default=42, help="Random seed for reproducibility. Defaults to 42.")
    parser.add_argument("-d", "--do_not_plot", action="store_false", help="Do not plot reconstructed pedigree(s).")
    parser.add_argument(
        "-w", "--write_alternates", action="store_true", help="Write outputs of alternate pedigrees to disk."
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose output (INFO-level logging).")
    return parser.parse_args(remaining_args)


def main():
    args = parse_arguments()
    logging_level = logging.INFO if args.verbose else logging.WARNING
    logging.basicConfig(level=logging_level, format="%(levelname)s - %(message)s")

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    with logging_redirect_tqdm():
        print("Beginning pedigree reconstruction...")
        pedigree_reconstructor = PedigreeReconstructor(
            relations_path=args.relations,
            nodes_path=args.nodes,
            outputs_dir=output_dir,
            max_candidate_pedigrees=args.max_candidate_pedigrees,
            epsilon=args.epsilon,
            plot=args.do_not_plot,
            write_alternate_pedigrees=args.write_alternates,
            random_seed=args.seed,
        )
        pedigree_reconstructor.find_best_pedigree()
    print(f"Finished pedigree reconstruction. Outputs written to {output_dir.resolve()}.")


if __name__ == "__main__":
    main()
