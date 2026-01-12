from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def load_results(results_path: Path | str) -> pd.DataFrame:
    results_path = Path(results_path)
    results_df = pd.read_csv(results_path)
    numeric_columns = [
        "num_first_and_second_degree_relations",
        "runtime_seconds",
        "peak_memory_mb",
    ]
    for column in numeric_columns:
        if column in results_df.columns:
            results_df[column] = pd.to_numeric(results_df[column], errors="coerce")

    return results_df.dropna(subset=numeric_columns)


def _aggregate_by_relations(results_df: pd.DataFrame) -> pd.DataFrame:
    grouped = results_df.groupby(["dataset", "num_first_and_second_degree_relations"], as_index=False).agg(
        runtime_mean=("runtime_seconds", "mean"),
        runtime_std=("runtime_seconds", "std"),
        runtime_count=("runtime_seconds", "count"),
        peak_mean=("peak_memory_mb", "mean"),
        peak_std=("peak_memory_mb", "std"),
        peak_count=("peak_memory_mb", "count"),
    )
    grouped["runtime_std"] = grouped["runtime_std"].fillna(0.0)
    grouped["peak_std"] = grouped["peak_std"].fillna(0.0)
    grouped["runtime_count"] = grouped["runtime_count"].clip(lower=1)
    grouped["peak_count"] = grouped["peak_count"].clip(lower=1)
    grouped["runtime_ci"] = 1.96 * grouped["runtime_std"] / np.sqrt(grouped["runtime_count"])
    grouped["peak_ci"] = 1.96 * grouped["peak_std"] / np.sqrt(grouped["peak_count"])
    return grouped.sort_values(["dataset", "num_first_and_second_degree_relations"])


def plot_relations_vs_runtime(results_df: pd.DataFrame, *, output_path: Path) -> None:
    grouped = _aggregate_by_relations(results_df)
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)

    for dataset, subset in grouped.groupby("dataset"):
        ax.errorbar(
            subset["num_first_and_second_degree_relations"],
            subset["runtime_mean"] / 60.0,
            yerr=subset["runtime_ci"] / 60.0,
            fmt="-o",
            label=dataset,
            capsize=10,
            markersize=10,
            elinewidth=8,
        )

    ax.set_xlabel("# of 1st- and 2nd-Degree Relations", fontsize=16)
    ax.set_ylabel("Runtime (minutes)", fontsize=16)
    ax.set_title("Runtime vs. Relations", fontsize=18)
    ax.tick_params(axis="both", labelsize=14)
    ax.grid(True, linewidth=0.8, alpha=0.4)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=-2)
    sns.despine(ax=ax)

    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def plot_rss_vs_relations(results_df: pd.DataFrame, *, output_path: Path) -> None:
    grouped = _aggregate_by_relations(results_df)
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)

    for dataset, subset in grouped.groupby("dataset"):
        ax.errorbar(
            subset["num_first_and_second_degree_relations"],
            subset["peak_mean"],
            yerr=subset["peak_ci"],
            fmt="-o",
            label=dataset,
            capsize=10,
            markersize=10,
            elinewidth=8,
        )

    ax.set_xlabel("# of 1st- and 2nd-Degree Relations", fontsize=16)
    ax.set_ylabel("Peak RSS (MB)", fontsize=16)
    ax.set_title("Peak Resident Set Size (RSS) vs. Relations", fontsize=18)
    ax.tick_params(axis="both", labelsize=14)
    ax.grid(True, linewidth=0.8, alpha=0.4)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    sns.despine(ax=ax)

    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def plot_results(results_path: Path | str) -> None:
    script_dir = Path(__file__).resolve().parent
    plots_dir = script_dir / "results" / "resource_benchmark" / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    relations_vs_runtime_path = plots_dir / "gurgy_runtime_vs_relations.pdf"
    rss_vs_relations_path = plots_dir / "gurgy_rss_vs_relations.pdf"

    results_df = load_results(results_path)
    plot_relations_vs_runtime(results_df, output_path=relations_vs_runtime_path)
    plot_rss_vs_relations(results_df, output_path=rss_vs_relations_path)


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    results_path = script_dir / "results" / "resource_benchmark" / "gurgy_results.csv"
    plot_results(results_path)


if __name__ == "__main__":
    main()
