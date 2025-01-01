from statistics import mean, stdev
from collections import defaultdict
from simulated_pedigree import SimulatedPedigree


def simulate(p_mask_node: float, error_rate_scale: float, random_seed: int) -> tuple[dict[str, int | float], dict[str, float]]:
    simulated_pedigree = SimulatedPedigree(p_mask_node=p_mask_node, 
                                           error_rate_scale=error_rate_scale, 
                                           random_seed=random_seed
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
        pedigree_statistics, metrics = simulate(p_mask_node=p_mask_node, error_rate_scale=error_rate_scale, random_seed=idx)
        for statistic, value in pedigree_statistics.items():
            experiment_pedigree_statistics[statistic].append(value)
        for metric, value in metrics.items():
            experiment_metrics[metric].append(value)

    mean_pedigree_statistics = {statistic: mean(values) for statistic, values in experiment_pedigree_statistics.items()}
    stdev_pedigree_statistics = {statistic: stdev(values) for statistic, values in experiment_pedigree_statistics.items()}
    for statistic in experiment_pedigree_statistics:
        print(f"{statistic}: {round(mean_pedigree_statistics[statistic], 2)} ± {round(stdev_pedigree_statistics[statistic], 2)}")

    mean_metrics = {metric: mean(values) for metric, values in experiment_metrics.items()}
    stdev_metrics = {metric: stdev(values) for metric, values in experiment_metrics.items()}
    for metric in experiment_metrics:
        print(f"{metric}: {round(mean_metrics[metric], 2)} ± {round(stdev_metrics[metric], 2)}")

def main():
    run_experiment(p_mask_node=0.4, error_rate_scale=1, num_simulations=10)

if __name__ == "__main__":
    main()
