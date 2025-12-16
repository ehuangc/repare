from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def load_results(results_path: Path | str) -> pd.DataFrame:
    results_path = Path(results_path)
    results_df = pd.read_csv(results_path)
    numeric_columns = [
        "Total Node Count",
        "Proportion of Inbred Nodes",
        "Pairwise Relation Accuracy",
        "Relation Precision",
        "Relation Recall",
        "Relation F1",
        "Pairwise Degree Accuracy",
        "Degree Precision",
        "Degree Recall",
        "Degree F1",
        "Connectivity R-squared",
        "p(Mask Node)",
        "Coverage Level",
    ]
    for column in numeric_columns:
        if column in results_df.columns:
            results_df[column] = pd.to_numeric(results_df[column], errors="coerce")
    return results_df.dropna(subset=["Relation F1", "Degree F1"])


def plot_pedigree_level_results(
    results_df: pd.DataFrame, *, results_path: Path | str, output_name: str | None = None
) -> None:
    plots_dir = Path(results_path).resolve().parent.parent / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    metric_columns: dict[str, str] = {
        "Relation F1": "Relation F1",
        "Degree F1": "Degree F1",
    }
    feature_columns: dict[str, str] = {
        "Total Node Count": "Pedigree Size (# of Individuals)",
        "Proportion of Inbred Nodes": "Inbreeding Proportion",
    }

    p_mask_node = results_df["p(Mask Node)"].iloc[0]
    coverage_level = results_df["Coverage Level"].iloc[0]
    if output_name is None:
        output_name = f"pedigree_level_results_p_mask_node={p_mask_node}_coverage_level={coverage_level}x.pdf"

    n_rows = len(metric_columns)
    n_cols = len(feature_columns)
    with mpl.rc_context(
        {
            "figure.constrained_layout.h_pad": 0.1,
            "figure.constrained_layout.w_pad": 0.15,
        }
    ):
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(10, 6.5), constrained_layout=True)
        axes = axes.reshape(n_rows, n_cols)

        for row_idx, (metric_col, metric_label) in enumerate(metric_columns.items()):
            for col_idx, (feature_col, feature_label) in enumerate(feature_columns.items()):
                ax = axes[row_idx][col_idx]
                sns.regplot(
                    data=results_df,
                    x=feature_col,
                    y=metric_col,
                    ax=ax,
                    scatter_kws={"s": 40, "alpha": 0.65},
                    line_kws={"lw": 2.0},
                    lowess=True,
                )
                if row_idx == n_rows - 1:
                    ax.set_xlabel(feature_label, fontsize=16)
                else:
                    ax.set_xlabel("")
                if col_idx == 0:
                    ax.set_ylabel(metric_label, fontsize=16)
                else:
                    ax.set_ylabel("")
                if row_idx == 0:
                    ax.set_title(feature_label, fontsize=16)
                ax.set_ylim(0.0, 1.05)
                ax.tick_params(axis="both", labelsize=14)
                ax.grid(True, linewidth=0.8, alpha=0.4)

        plt.suptitle("Pedigree-Level Reconstruction Performance", fontsize=18, y=1.09)
        fig.text(
            0.5,
            1.045,
            f"p(Mask Node)={p_mask_node}, Simulated Sequence Coverage={coverage_level}x",
            ha="center",
            va="top",
            fontsize=16,
        )
        sns.despine()

        output_path = plots_dir / output_name
        fig.savefig(output_path, bbox_inches="tight")
        plt.close(fig)


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    paths: list[Path] = [
        script_dir / "results" / "parameter_experiment" / "data" / "p_mask_node=0.4_coverage_level=0.5.csv",
        script_dir / "results" / "parameter_experiment" / "data" / "p_mask_node=0.4_coverage_level=0.1.csv",
    ]

    for results_path in paths:
        results_df = load_results(results_path)
        plot_pedigree_level_results(results_df, results_path=results_path)


if __name__ == "__main__":
    main()
