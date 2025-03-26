import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statistics import mean, stdev


def plot_pedigree_summary_statistics(results_dir: str) -> None:
    # We can use any results file to get pedigree statistics because 
    # all experiments are run on same simulated pedigrees since seed is fixed
    results_path = os.listdir(results_dir)[0]
    results_df = pd.read_csv(os.path.join(results_dir, results_path))

    pedigree_sizes = results_df["Total Node Count"].values
    inbred_proportions = results_df["Proportion of Inbred Nodes"].values

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    sns.histplot(pedigree_sizes, ax=axes[0])
    axes[0].set_title("Pedigree Size Distribution")
    axes[0].set_xlabel("# of Individuals in Full Simulated Pedigree")
    
    sns.histplot(inbred_proportions, ax=axes[1])
    axes[1].set_title("Inbreeding Proportion Distrbution")
    axes[1].set_xlabel("Proportion of Inbred Individuals in Full Simulated Pedigree")
    plt.savefig("simulation_results/plots/pedigree_summary_statistics.png", dpi=600)

def plot_results(results_dir: str) -> None:
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
        
        plt.figure(figsize=(8, 6))
        # Set vmin and vmax so relation and degree F1 scores are on the same color scale
        ax = sns.heatmap(heatmap_data, annot=True, fmt=".2f", cmap="Blues", cbar_kws={"label": metric}, vmin=0.5, vmax=1.0)
        ax.figure.axes[-1].yaxis.label.set_size(12)  # Set colorbar label size
        plt.title(f"{metric} Heatmap", fontsize=14)
        plt.xlabel("Kinship Relation Error Rate Scale", fontsize=12)
        plt.ylabel("p(Mask Node)", fontsize=12)
        plt.savefig(f"simulation_results/plots/{metric.lower().replace(' ', '_')}_heatmap.png", dpi=600)

def main():
    os.makedirs("simulation_results/plots", exist_ok=True)
    results_dir = "simulation_results/data"
    plot_pedigree_summary_statistics(results_dir)
    plot_results(results_dir)

if __name__ == "__main__":
    main()
