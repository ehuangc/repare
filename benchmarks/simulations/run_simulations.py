from statistics import mean, stdev
from collections import defaultdict
from simulated_pedigree import SimulatedPedigree


def simulate(random_seed: int) -> dict[str, float]:
    simulated_pedigree = SimulatedPedigree(random_seed=random_seed)
    simulated_pedigree.create_pedigree()
    simulated_pedigree.mask_data()
    simulated_pedigree.run_algorithm()
    metrics = simulated_pedigree.get_metrics()
    return metrics

def run_experiment(num_simulations: int = 100) -> dict[str, float]:
    experiment_metrics = defaultdict(list)
    for idx in range(num_simulations):
        metrics = simulate(random_seed=idx)
        for metric, value in metrics.items():
            experiment_metrics[metric].append(value)
    
    mean_metrics = {metric: mean(values) for metric, values in experiment_metrics.items()}
    stdev_metrics = {metric: stdev(values) for metric, values in experiment_metrics.items()}
    for metric in experiment_metrics:
        print(f"{metric}: {round(mean_metrics[metric], 2)} ± {round(stdev_metrics[metric], 2)}")

def main():
    run_experiment(num_simulations=10)

if __name__ == "__main__":
    main()
