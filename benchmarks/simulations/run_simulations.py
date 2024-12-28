from collections import defaultdict
from simulated_pedigree import SimulatedPedigree


def simulate() -> dict[str, float]:
    simulated_pedigree = SimulatedPedigree()
    simulated_pedigree.create_pedigree()
    simulated_pedigree.mask_data()
    simulated_pedigree.run_algorithm()
    metrics = simulated_pedigree.get_metrics()
    return metrics

def run_experiment(num_simulations: int = 100) -> dict[str, float]:
    experiment_metrics = defaultdict(list)
    for _ in range(num_simulations):
        metrics = simulate()
        for metric, value in metrics.items():
            experiment_metrics[metric].append(value)
    
    mean_metrics = {metric: sum(values) / len(values) for metric, values in experiment_metrics.items()}
    print(mean_metrics)

def main():
    run_experiment(num_simulations=10)

if __name__ == "__main__":
    main()
