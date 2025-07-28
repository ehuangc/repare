import os
from collections import defaultdict

import pandas as pd
from simulator.simulated_pedigree import SimulatedPedigree


def simulate(
    p_mask_node: float, error_rate_scale: float, max_candidate_pedigrees: int, epsilon: float, random_seed: int
) -> tuple[dict[str, int | float], dict[str, float]]:
    simulated_pedigree = SimulatedPedigree(
        p_mask_node=p_mask_node,
        error_rate_scale=error_rate_scale,
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
    p_mask_node: float,
    error_rate_scale: float,
    max_candidate_pedigrees: int,
    epsilon: float,
    num_simulations: int = 100,
) -> None:
    print(
        f"Running {num_simulations} simulations: "
        f"max_candidate_pedigrees={max_candidate_pedigrees}, "
        f"epsilon={epsilon}, "
        f"p_mask_node={p_mask_node}, "
        f"error_rate_scale={error_rate_scale}"
    )
    experiment_pedigree_statistics = defaultdict(list)
    experiment_metrics = defaultdict(list)

    for idx in range(num_simulations):
        pedigree_statistics, metrics = simulate(
            p_mask_node=p_mask_node,
            error_rate_scale=error_rate_scale,
            max_candidate_pedigrees=max_candidate_pedigrees,
            epsilon=epsilon,
            random_seed=idx,
        )
        for statistic, value in pedigree_statistics.items():
            experiment_pedigree_statistics[statistic].append(value)
        for metric, value in metrics.items():
            experiment_metrics[metric].append(value)

    results_df = pd.concat(
        [pd.DataFrame.from_dict(experiment_pedigree_statistics), pd.DataFrame.from_dict(experiment_metrics)], axis=1
    )
    results_df["Max Candidate Pedigrees"] = max_candidate_pedigrees
    results_df["Epsilon"] = epsilon
    os.makedirs("results/sampling_experiment/data", exist_ok=True)
    results_df.to_csv(
        f"results/sampling_experiment/data/max_candidate_pedigrees={max_candidate_pedigrees}_epsilon={epsilon}.csv",
        index=False,
    )


def main():
    for max_candidate_pedigrees in [10, 100, 1000, 10000]:
        for epsilon in [0.0, 0.2, 0.4]:
            run_experiment(
                p_mask_node=0.4,
                error_rate_scale=1,
                max_candidate_pedigrees=max_candidate_pedigrees,
                epsilon=epsilon,
                num_simulations=100,
            )


if __name__ == "__main__":
    main()
