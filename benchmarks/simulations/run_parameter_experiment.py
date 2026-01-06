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
    random_seed: int,
) -> tuple[dict[str, int | float], dict[str, float]]:
    pedigree_data_dir = (
        script_dir
        / "results"
        / "parameter_experiment"
        / "pedigree_data"
        / f"p_mask_node={p_mask_node}_coverage_level={coverage_level}"
        / f"pedigree{random_seed}"
    )
    pedigree_data_dir.mkdir(parents=True, exist_ok=True)
    simulated_pedigree = SimulatedPedigree(
        pedigree_data_dir=pedigree_data_dir,
        p_mask_node=p_mask_node,
        coverage_level=coverage_level,
        random_seed=random_seed,
    )
    simulated_pedigree.create_pedigree()
    simulated_pedigree.mask_and_corrupt_data()
    simulated_pedigree.run_algorithm()
    pedigree_statistics = simulated_pedigree.get_pedigree_statistics()
    metrics = simulated_pedigree.get_metrics()
    return pedigree_statistics, metrics


def run_experiment(script_dir: Path, p_mask_node: float, coverage_level: float, num_simulations: int = 100) -> None:
    print(f"Running {num_simulations} simulations: p_mask_node={p_mask_node}, coverage_level={coverage_level}x")

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
            seeds,
        ):
            for k, v in pedigree_statistics.items():
                experiment_pedigree_statistics[k].append(v)
            for k, v in metrics.items():
                experiment_metrics[k].append(v)

    results_df = pd.concat(
        [pd.DataFrame.from_dict(experiment_pedigree_statistics), pd.DataFrame.from_dict(experiment_metrics)], axis=1
    )
    results_df["p(Mask Node)"] = p_mask_node
    results_df["Coverage Level"] = coverage_level
    results_dir = script_dir / "results" / "parameter_experiment" / "data"
    results_dir.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(
        results_dir / f"p_mask_node={p_mask_node}_coverage_level={coverage_level}.csv",
        index=False,
    )


def main():
    script_dir = Path(__file__).resolve().parent
    for p_mask_node in [0.0, 0.2, 0.4, 0.6]:
        for coverage_level in [0.1, 0.2, 0.5, 1.0]:
            run_experiment(
                script_dir=script_dir, p_mask_node=p_mask_node, coverage_level=coverage_level, num_simulations=100
            )


if __name__ == "__main__":
    main()
