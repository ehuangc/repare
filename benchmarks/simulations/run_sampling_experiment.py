from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
from pathlib import Path

import pandas as pd
from simulator.simulated_pedigree import SimulatedPedigree


def simulate(
    script_dir: Path,
    p_mask_node: float,
    coverage_level: float,
    max_candidate_pedigrees: int,
    epsilon: float,
    random_seed: int,
) -> tuple[dict[str, int | float], dict[str, float]]:
    pedigree_data_dir = (
        script_dir
        / "results"
        / "sampling_experiment"
        / "pedigree_data"
        / f"max_candidate_pedigrees={max_candidate_pedigrees}_epsilon={epsilon}"
        / f"pedigree{random_seed}"
    )
    pedigree_data_dir.mkdir(parents=True, exist_ok=True)
    simulated_pedigree = SimulatedPedigree(
        pedigree_data_dir=pedigree_data_dir,
        p_mask_node=p_mask_node,
        coverage_level=coverage_level,
        max_candidate_pedigrees=max_candidate_pedigrees,
        epsilon=epsilon,
        random_seed=random_seed,
    )
    simulated_pedigree.create_pedigree()
    simulated_pedigree.mask_and_corrupt_data()
    simulated_pedigree.run_algorithm()
    pedigree_statistics = simulated_pedigree.get_pedigree_statistics()
    metrics = simulated_pedigree.get_metrics()
    return pedigree_statistics, metrics


def run_experiment(
    script_dir: Path,
    p_mask_node: float,
    coverage_level: float,
    max_candidate_pedigrees: int,
    epsilon: float,
    num_simulations: int = 100,
) -> None:
    print(
        f"Running {num_simulations} simulations: "
        f"max_candidate_pedigrees={max_candidate_pedigrees}, "
        f"epsilon={epsilon}, "
        f"p_mask_node={p_mask_node}, "
        f"coverage_level={coverage_level}x"
    )

    # Parallelize simulations across CPU cores
    seeds = list(range(num_simulations))
    experiment_pedigree_statistics = defaultdict(list)
    experiment_metrics = defaultdict(list)

    # Use as many workers as (logical) CPU cores by default
    with ProcessPoolExecutor() as ex:
        for pedigree_statistics, metrics in ex.map(
            simulate,
            repeat(script_dir),
            repeat(p_mask_node),
            repeat(coverage_level),
            repeat(max_candidate_pedigrees),
            repeat(epsilon),
            seeds,
        ):
            for k, v in pedigree_statistics.items():
                experiment_pedigree_statistics[k].append(v)
            for k, v in metrics.items():
                experiment_metrics[k].append(v)

    results_df = pd.concat(
        [pd.DataFrame.from_dict(experiment_pedigree_statistics), pd.DataFrame.from_dict(experiment_metrics)], axis=1
    )
    results_df["Max Candidate Pedigrees"] = max_candidate_pedigrees
    results_df["Epsilon"] = epsilon
    results_dir = script_dir / "results" / "sampling_experiment" / "data"
    results_dir.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(
        results_dir / f"max_candidate_pedigrees={max_candidate_pedigrees}_epsilon={epsilon}.csv",
        index=False,
    )


def main():
    script_dir = Path(__file__).resolve().parent
    for max_candidate_pedigrees in [10, 100, 1000, 10000]:
        for epsilon in [0.0, 0.2, 0.4]:
            run_experiment(
                script_dir=script_dir,
                p_mask_node=0.4,
                coverage_level=0.5,
                max_candidate_pedigrees=max_candidate_pedigrees,
                epsilon=epsilon,
                num_simulations=100,
            )


if __name__ == "__main__":
    main()
