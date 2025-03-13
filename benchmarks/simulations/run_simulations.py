import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
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
    
    results_df = pd.concat([pd.DataFrame.from_dict(experiment_pedigree_statistics), pd.DataFrame.from_dict(experiment_metrics)], axis=1)
    results_df["p(Mask Node)"] = p_mask_node
    results_df["Error Rate Scale"] = error_rate_scale
    os.makedirs("simulation_results/data", exist_ok=True)
    results_df.to_csv(f"simulation_results/data/p_mask_node={p_mask_node}_error_rate_scale={error_rate_scale}.csv", index=False)

def plot_results():
    results_dir = "simulation_results/data"
    p_mask_nodes = []
    error_rate_scales = []
    mean_relation_f1s = []
    mean_degree_f1s = []

    for file in os.listdir(results_dir):
        results_df = pd.read_csv(os.path.join(results_dir, file))
        p_mask_node = results_df["p(Mask Node)"].iloc[0]
        error_rate_scale = results_df["Error Rate Scale"].iloc[0]
        mean_relation_f1 = mean(results_df["Relation F1"])
        mean_degree_f1 = mean(results_df["Degree F1"])

        p_mask_nodes.append(p_mask_node)
        error_rate_scales.append(error_rate_scale)
        mean_relation_f1s.append(mean_relation_f1)
        mean_degree_f1s.append(mean_degree_f1)

    results_df = pd.DataFrame({
        "p_mask_node": p_mask_nodes,
        "error_rate_scale": error_rate_scales,
        "mean_relation_f1": mean_relation_f1s,
        "mean_degree_f1": mean_degree_f1s
    })
    relation_f1_heatmap_data = results_df.pivot(index="p_mask_node", columns="error_rate_scale", values="mean_relation_f1")
    degree_f1_heatmap_data = results_df.pivot(index="p_mask_node", columns="error_rate_scale", values="mean_degree_f1")

    for heatmap_data, metric in zip([relation_f1_heatmap_data, degree_f1_heatmap_data], ["Relation F1", "Degree F1"]):
        heatmap_data = heatmap_data.sort_index(ascending=False)  # Error rate scale increases from left to right
        heatmap_data = heatmap_data.sort_index(axis=1, ascending=True)  # p(Mask Node) increases from top to bottom
        heatmap_data.rename(columns={1.0: "1.0\n(~0.5x coverage)"}, inplace=True)
        
        plt.figure(figsize=(8, 6), dpi=1200)
        # Set vmin and vmax so relation and degree F1 scores are on the same color scale
        ax = sns.heatmap(heatmap_data, annot=True, fmt=".2f", cmap="Blues", cbar_kws={"label": metric}, vmin=0.5, vmax=1.0)
        ax.figure.axes[-1].yaxis.label.set_size(12)  # Set colorbar label size
        plt.title(f"{metric} Heatmap", fontsize=14)
        plt.xlabel("Error Rate Scale", fontsize=12)
        plt.ylabel("p(Mask Node)", fontsize=12)
        plt.savefig(f"simulation_results/{metric.lower().replace(' ', '_')}_heatmap.png")

def main():
    for p_mask_node in [0.0, 0.2, 0.4, 0.6]:
        for error_rate_scale in [0.0, 0.5, 1.0, 2]:
            run_experiment(p_mask_node=p_mask_node, error_rate_scale=error_rate_scale, num_simulations=10)
    plot_results()

if __name__ == "__main__":
    main()
