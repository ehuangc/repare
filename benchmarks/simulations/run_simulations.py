from statistics import mean, stdev
from collections import defaultdict
from simulated_pedigree import SimulatedPedigree


def simulate(p_mask_node: float, random_seed: int) -> dict[str, float]:
    simulated_pedigree = SimulatedPedigree(p_mask_node=p_mask_node, random_seed=random_seed)
    simulated_pedigree.create_pedigree()
    simulated_pedigree.mask_and_corrupt_data()
    simulated_pedigree.run_algorithm()
    metrics = simulated_pedigree.get_metrics()
    return metrics

def run_experiment(p_mask_node: float, num_simulations: int = 100) -> dict[str, float]:
    experiment_metrics = defaultdict(list)
    for idx in range(num_simulations):
        metrics = simulate(p_mask_node=p_mask_node, random_seed=idx)
        for metric, value in metrics.items():
            experiment_metrics[metric].append(value)
    
    mean_metrics = {metric: mean(values) for metric, values in experiment_metrics.items()}
    stdev_metrics = {metric: stdev(values) for metric, values in experiment_metrics.items()}
    for metric in experiment_metrics:
        print(f"{metric}: {round(mean_metrics[metric], 2)} ± {round(stdev_metrics[metric], 2)}")

def main():
    run_experiment(p_mask_node=0.4, num_simulations=10)

if __name__ == "__main__":
    main()
