import os
from statistics import mean

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot_pedigree_summary_statistics(results_dir: str) -> None:
    # We can use any results file to get pedigree statistics because
    # all experiments are run on same simulated pedigrees since seed is fixed
    results_path = os.listdir(results_dir)[0]
    results_df = pd.read_csv(os.path.join(results_dir, results_path))
    pedigree_sizes = results_df["Total Node Count"].values
    inbred_proportions = results_df["Proportion of Inbred Nodes"].values
    has_children_proportions = results_df["Proportion of Non-Final-Generation Nodes with Children"].values
    mean_children_count = results_df["Mean Children Count per Parent"].values

    with mpl.rc_context(
        {
            "figure.constrained_layout.h_pad": 0.1,
            "figure.constrained_layout.w_pad": 0.15,
        }
    ):
        fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)
        axes = axes.flatten()

        plt.suptitle("Pedigree Summary Statistics (Before Masking Nodes)", fontsize=16)

        for ax in axes:
            ax.set_ylabel("Pedigree Count", fontsize=14)
            ax.tick_params(axis="x", labelsize=12)
            ax.tick_params(axis="y", labelsize=12)

        sns.histplot(pedigree_sizes, ax=axes[0])
        axes[0].set_title("Pedigree Size Distribution", fontsize=14)
        axes[0].set_xlabel("# of Individuals", fontsize=14)

        sns.histplot(inbred_proportions, ax=axes[1])
        axes[1].set_title("Inbreeding Proportion Distribution", fontsize=14)
        axes[1].set_xlabel("Proportion of Inbred Individuals", fontsize=14)

        sns.histplot(has_children_proportions, ax=axes[2])
        axes[2].set_title("Has Children Proportion Distribution", fontsize=14)
        axes[2].set_xlabel("Proportion of Non-Final-Generation\nIndividuals with Children", fontsize=14)

        sns.histplot(mean_children_count, ax=axes[3])
        axes[3].set_title("Mean Children Count Distribution", fontsize=14)
        axes[3].set_xlabel("Mean # of Children per Parent", fontsize=14)

        plt.savefig(
            "results/parameter_experiment/plots/pedigree_summary_statistics.pdf",
            bbox_inches="tight",
        )


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

    results_df = pd.DataFrame(
        {
            "p_mask_node": p_mask_nodes,
            "error_rate_scale": error_rate_scales,
            "mean_relation_f1": mean_relation_f1s,
            "mean_degree_f1": mean_degree_f1s,
        }
    )
    relation_f1_heatmap_data = results_df.pivot(
        index="p_mask_node", columns="error_rate_scale", values="mean_relation_f1"
    )
    degree_f1_heatmap_data = results_df.pivot(index="p_mask_node", columns="error_rate_scale", values="mean_degree_f1")

    for heatmap_data, metric in zip(
        [relation_f1_heatmap_data, degree_f1_heatmap_data], ["Relation F1", "Degree F1"], strict=True
    ):
        # p(Mask Node) increases from bottom to top
        heatmap_data = heatmap_data.sort_index(ascending=False)
        # Error rate scale increases from left to right
        heatmap_data = heatmap_data.sort_index(axis=1, ascending=True)
        heatmap_data.rename(columns={1.0: "1.0\n(~0.5x coverage)"}, inplace=True)

        plt.figure(figsize=(8, 6))
        # Set vmin and vmax so relation and degree F1 scores are on the same color scale
        ax = sns.heatmap(
            heatmap_data,
            annot=True,
            fmt=".2f",
            cmap="Blues",
            vmin=0.5,
            vmax=1.0,
            annot_kws={"size": 14},
        )
        # Set colorbar label padding
        ax.figure.axes[-1].yaxis.labelpad = 10
        # Set colorbar tick label size
        ax.figure.axes[-1].tick_params(labelsize=14)
        plt.title(f"{metric} Scores", fontsize=18, pad=10)
        plt.xlabel("Kinship Relation Error Rate Scale", fontsize=16, labelpad=10)
        plt.ylabel("p(Mask Node)", fontsize=16, labelpad=10)
        ax.tick_params(axis="x", labelsize=14)
        ax.tick_params(axis="y", labelsize=14)

        plt.savefig(
            f"results/parameter_experiment/plots/{metric.lower().replace(' ', '_')}_heatmap.pdf",
            bbox_inches="tight",
        )


def main():
    os.makedirs("results/parameter_experiment/plots", exist_ok=True)
    results_dir = "results/parameter_experiment/data"
    plot_pedigree_summary_statistics(results_dir)
    plot_results(results_dir)


if __name__ == "__main__":
    main()
