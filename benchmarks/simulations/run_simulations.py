import os
from collections import defaultdict

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


def run_experiment(p_mask_node: float, error_rate_scale: float, num_simulations: int = 100) -> dict[str, float]:
    print(f"Running {num_simulations} simulations: p_mask_node={p_mask_node}, error_rate_scale={error_rate_scale}")
    experiment_pedigree_statistics = defaultdict(list)
    experiment_metrics = defaultdict(list)

    for idx in range(num_simulations):
        pedigree_statistics, metrics = simulate(
            p_mask_node=p_mask_node, error_rate_scale=error_rate_scale, random_seed=idx
        )
        for statistic, value in pedigree_statistics.items():
            experiment_pedigree_statistics[statistic].append(value)
        for metric, value in metrics.items():
            experiment_metrics[metric].append(value)

    results_df = pd.concat(
        [pd.DataFrame.from_dict(experiment_pedigree_statistics), pd.DataFrame.from_dict(experiment_metrics)], axis=1
    )
    results_df["p(Mask Node)"] = p_mask_node
    results_df["Error Rate Scale"] = error_rate_scale
    os.makedirs("results/data", exist_ok=True)
    results_df.to_csv(f"results/data/p_mask_node={p_mask_node}_error_rate_scale={error_rate_scale}.csv", index=False)


def main():
    for p_mask_node in [0.0, 0.2, 0.4, 0.6]:
        for error_rate_scale in [0.0, 0.5, 1.0, 2]:
            run_experiment(p_mask_node=p_mask_node, error_rate_scale=error_rate_scale, num_simulations=100)


if __name__ == "__main__":
    main()
