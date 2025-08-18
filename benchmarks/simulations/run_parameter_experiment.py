import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat

import pandas as pd
from simulator.simulated_pedigree import SimulatedPedigree


def simulate(
    p_mask_node: float, error_rate_scale: float, random_seed: int
) -> tuple[dict[str, int | float], dict[str, float]]:
    simulated_pedigree = SimulatedPedigree(
        p_mask_node=p_mask_node, error_rate_scale=error_rate_scale, random_seed=random_seed
    )
    simulated_pedigree.create_pedigree()
    simulated_pedigree.mask_and_corrupt_data()
    simulated_pedigree.run_algorithm()
    pedigree_statistics = simulated_pedigree.get_pedigree_statistics()
    metrics = simulated_pedigree.get_metrics()
    return pedigree_statistics, metrics


def run_experiment(p_mask_node: float, error_rate_scale: float, num_simulations: int = 100) -> None:
    print(f"Running {num_simulations} simulations: p_mask_node={p_mask_node}, error_rate_scale={error_rate_scale}")

    # Parallelize simulations across CPU cores
    seeds = list(range(num_simulations))
    experiment_pedigree_statistics = defaultdict(list)
    experiment_metrics = defaultdict(list)

    # Use as many workers as (logical) CPU cores by default
    with ProcessPoolExecutor() as ex:
        for pedigree_statistics, metrics in ex.map(
            simulate,
            repeat(p_mask_node),
            repeat(error_rate_scale),
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
    results_df["Error Rate Scale"] = error_rate_scale
    os.makedirs("results/parameter_experiment/data", exist_ok=True)
    results_df.to_csv(
        f"results/parameter_experiment/data/p_mask_node={p_mask_node}_error_rate_scale={error_rate_scale}.csv",
        index=False,
    )


def main():
    for p_mask_node in [0.0, 0.2, 0.4, 0.6]:
        for error_rate_scale in [0.0, 0.5, 1.0, 2.0]:
            run_experiment(p_mask_node=p_mask_node, error_rate_scale=error_rate_scale, num_simulations=100)


if __name__ == "__main__":
    main()
