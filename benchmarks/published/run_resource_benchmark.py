import math
import multiprocessing as mp
import random
import resource
import sys
import tempfile
import time
from pathlib import Path

import pandas as pd

from repare.pedigree_reconstructor import PedigreeReconstructor


def reconstruct(nodes_path, relations_path, outputs_dir, conn):
    Path(outputs_dir).mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    reconstructor = PedigreeReconstructor(
        relations_path=relations_path,
        nodes_path=nodes_path,
        outputs_dir=outputs_dir,
    )
    reconstructor.find_best_pedigree()
    runtime = time.perf_counter() - start
    usage = resource.getrusage(resource.RUSAGE_SELF)
    peak_rss = usage.ru_maxrss
    if sys.platform != "darwin":
        peak_rss *= 1024
    conn.send({"runtime": runtime, "peak_bytes": peak_rss})
    conn.close()


def run_trial(nodes_df, relations_df):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        nodes_file = tmpdir_path / "nodes.csv"
        relations_file = tmpdir_path / "relations.csv"
        nodes_df.to_csv(nodes_file, index=False)
        relations_df.to_csv(relations_file, index=False)
        outputs_dir = tmpdir_path / "outputs"

        # Use a fresh interpreter per run to avoid forked memory sharing skewing RSS
        parent_conn, child_conn = mp.get_context("spawn").Pipe(duplex=False)
        process = mp.get_context("spawn").Process(
            target=reconstruct,
            args=(str(nodes_file), str(relations_file), str(outputs_dir), child_conn),
        )
        process.start()
        child_conn.close()
        result = parent_conn.recv()
        parent_conn.close()
        process.join()
        if process.exitcode != 0:
            raise RuntimeError(f"Worker exited with code {process.exitcode}.")
        runtime = float(result["runtime"])
        peak_mb = float(result["peak_bytes"]) / (1024**2)
        return runtime, peak_mb


def node_counts(total_nodes):
    counts = set()
    for frac in [0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.5, 0.9, 0.95, 1.0]:
        counts.add(max(1, min(total_nodes, math.ceil(total_nodes * frac))))
    counts = sorted(counts)
    if counts[-1] != total_nodes:
        counts.append(total_nodes)
    return counts


def sample_node_subset(nodes_df, relations_df, target_node_count, rng):
    if target_node_count <= 0:
        raise ValueError("target_node_count must be positive.")

    node_ids = nodes_df["id"].tolist()
    target_node_count = min(target_node_count, len(node_ids))
    if target_node_count == len(node_ids):
        keep_ids = node_ids
    else:
        keep_ids = rng.sample(node_ids, target_node_count)

    keep_ids_set = set(keep_ids)
    filtered_nodes = nodes_df[nodes_df["id"].isin(keep_ids_set)].copy()
    filtered_relations = relations_df[
        relations_df["id1"].isin(keep_ids_set) & relations_df["id2"].isin(keep_ids_set)
    ].copy()
    return filtered_nodes, filtered_relations


def relation_degree_counts(relations_df):
    degree_counts = relations_df["degree"].value_counts()
    first_degree = int(degree_counts.get("1", 0))
    second_degree = int(degree_counts.get("2", 0))
    return first_degree, second_degree


def main():
    script_dir = Path(__file__).resolve().parent
    data_dir = script_dir / "data" / "gurgy"
    nodes_path = data_dir / "nodes.csv"
    relations_path = data_dir / "inferred_relations_READv2.csv"
    results_dir = script_dir / "results" / "resource_benchmark"
    results_dir.mkdir(parents=True, exist_ok=True)
    results_path = results_dir / "gurgy_results.csv"

    nodes_df = pd.read_csv(nodes_path, dtype=str, comment="#", keep_default_na=False)
    relations_df = pd.read_csv(relations_path, dtype=str, comment="#", keep_default_na=False)
    counts = node_counts(len(nodes_df))
    repeats = 3
    results = []
    rng = random.Random(42)

    print("Benchmarking resource usage...\n")
    for count in counts:
        nodes_subset, relations_subset = sample_node_subset(nodes_df, relations_df, count, rng)
        first_degree, second_degree = relation_degree_counts(relations_subset)
        for repeat_idx in range(repeats):
            runtime, peak_mb = run_trial(nodes_subset, relations_subset)
            results.append(
                {
                    "dataset": "gurgy",
                    "num_nodes": len(nodes_subset),
                    "num_first_degree_relations": first_degree,
                    "num_second_degree_relations": second_degree,
                    "num_first_and_second_degree_relations": first_degree + second_degree,
                    "repeat_idx": repeat_idx,
                    "runtime_seconds": runtime,
                    "peak_memory_mb": peak_mb,
                }
            )
            print(
                f"nodes={len(nodes_subset)} relations={first_degree + second_degree} "
                f"repeat_idx={repeat_idx} runtime={runtime:.2f}s peak_mem={peak_mb:.1f}MB"
                "\n"
            )

    pd.DataFrame(results).to_csv(results_path, index=False)
    print(f"\nBenchmark complete. Results saved to {results_path}.")


if __name__ == "__main__":
    main()
