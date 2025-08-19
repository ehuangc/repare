import os
from statistics import mean

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot_results(results_dir: str) -> None:
    max_candidate_pedigrees_values = []
    epsilons = []
    mean_relation_f1s = []
    mean_degree_f1s = []

    for file in os.listdir(results_dir):
        results_df = pd.read_csv(os.path.join(results_dir, file))
        max_candidate_pedigrees = results_df["Max Candidate Pedigrees"].iloc[0]
        epsilon = results_df["Epsilon"].iloc[0]
        mean_relation_f1 = mean(results_df["Relation F1"])
        mean_degree_f1 = mean(results_df["Degree F1"])

        max_candidate_pedigrees_values.append(max_candidate_pedigrees)
        epsilons.append(epsilon)
        mean_relation_f1s.append(mean_relation_f1)
        mean_degree_f1s.append(mean_degree_f1)

    results_df = pd.DataFrame(
        {
            "max_candidate_pedigrees": max_candidate_pedigrees_values,
            "epsilon": epsilons,
            "mean_relation_f1": mean_relation_f1s,
            "mean_degree_f1": mean_degree_f1s,
        }
    )
    relation_f1_heatmap_data = results_df.pivot(
        index="epsilon", columns="max_candidate_pedigrees", values="mean_relation_f1"
    )
    degree_f1_heatmap_data = results_df.pivot(
        index="epsilon", columns="max_candidate_pedigrees", values="mean_degree_f1"
    )

    for heatmap_data, metric in zip(
        [relation_f1_heatmap_data, degree_f1_heatmap_data], ["Relation F1", "Degree F1"], strict=True
    ):
        # Epsilon increases from bottom to top
        heatmap_data = heatmap_data.sort_index(ascending=False)
        # Max candidate pedigrees increases from left to right
        heatmap_data = heatmap_data.sort_index(axis=1, ascending=True)

        plt.figure(figsize=(8, 6))
        # Set vmin and vmax so relation and degree F1 scores are on the same color scale
        ax = sns.heatmap(
            heatmap_data,
            annot=True,
            fmt=".2f",
            cmap="Greens",
            vmin=0.5,
            vmax=1.0,
            annot_kws={"size": 14},
        )
        # Set colorbar label padding
        ax.figure.axes[-1].yaxis.labelpad = 10
        # Set colorbar tick label size
        ax.figure.axes[-1].tick_params(labelsize=14)
        plt.title(f"{metric} Scores", fontsize=18, pad=10)
        plt.xlabel("Max Candidate Pedigrees", fontsize=16, labelpad=10)
        plt.ylabel("Epsilon", fontsize=16, labelpad=10)
        ax.tick_params(axis="x", labelsize=14)
        ax.tick_params(axis="y", labelsize=14)

        plt.savefig(
            f"results/sampling_experiment/plots/{metric.lower().replace(' ', '_')}_heatmap.pdf",
            bbox_inches="tight",
        )


def main():
    os.makedirs("results/sampling_experiment/plots", exist_ok=True)
    results_dir = "results/sampling_experiment/data"
    plot_results(results_dir)


if __name__ == "__main__":
    main()
