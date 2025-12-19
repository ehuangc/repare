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


def plot_regressions(
    results_df: pd.DataFrame,
    *,
    regression_output_path: Path,
    p_mask_node: float,
    coverage_level: float,
) -> None:
    metric_columns: dict[str, str] = {
        "Relation F1": "Relation F1",
        "Degree F1": "Degree F1",
    }
    feature_columns: dict[str, str] = {
        "Total Node Count": "Pedigree Size (# of Individuals)",
        "Proportion of Inbred Nodes": "Inbreeding Proportion",
    }

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
                    ci=95,
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

        fig.savefig(regression_output_path, bbox_inches="tight")
        plt.close(fig)


def plot_histograms(results_df: pd.DataFrame, *, histogram_output_path: Path) -> None:
    histogram_metrics = [
        ("Relation F1", "Relation F1 Distribution"),
        ("Degree F1", "Degree F1 Distribution"),
    ]
    with mpl.rc_context(
        {
            "figure.constrained_layout.h_pad": 0.1,
            "figure.constrained_layout.w_pad": 0.15,
        }
    ):
        fig, axes = plt.subplots(1, len(histogram_metrics), figsize=(12, 4.5), constrained_layout=True)
        axes = axes.flatten()

        plt.suptitle("Pedigree-Level Reconstruction Performance Distribution", fontsize=18, y=1.14)
        fig.text(
            0.5,
            1.07,
            f"p(Mask Node)={results_df['p(Mask Node)'].iloc[0]}, "
            f"Simulated Sequence Coverage={results_df['Coverage Level'].iloc[0]}x",
            ha="center",
            va="top",
            fontsize=16,
        )

        for ax in axes:
            ax.set_ylabel("Pedigree Count", fontsize=16)
            ax.tick_params(axis="x", labelsize=14)
            ax.tick_params(axis="y", labelsize=14)

        for ax, (metric_col, metric_title) in zip(axes, histogram_metrics, strict=True):
            sns.histplot(results_df[metric_col], ax=ax)
            ax.set_title(metric_title, fontsize=16)
            ax.set_xlabel(metric_col, fontsize=16)

        fig.savefig(histogram_output_path, bbox_inches="tight")
        plt.close(fig)


def plot_pedigree_level_results(results_df: pd.DataFrame, *, results_path: Path | str) -> None:
    plots_dir = Path(results_path).resolve().parent.parent / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    p_mask_node = results_df["p(Mask Node)"].iloc[0]
    coverage_level = results_df["Coverage Level"].iloc[0]
    regression_output_filename = Path(
        f"pedigree_level_regressions_p_mask_node={p_mask_node}_coverage_level={coverage_level}x.pdf"
    )
    histogram_output_filename = Path(
        f"pedigree_level_histograms_p_mask_node={p_mask_node}_coverage_level={coverage_level}x.pdf"
    )
    regression_output_path = plots_dir / regression_output_filename
    histogram_output_path = plots_dir / histogram_output_filename

    plot_regressions(
        results_df,
        regression_output_path=regression_output_path,
        p_mask_node=p_mask_node,
        coverage_level=coverage_level,
    )
    plot_histograms(results_df, histogram_output_path=histogram_output_path)


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
