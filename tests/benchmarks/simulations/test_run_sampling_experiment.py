from pathlib import Path

import pytest

from benchmarks.simulations import run_sampling_experiment


class DummySimulatedPedigree:
    created_instances: list["DummySimulatedPedigree"] = []

    def __init__(self, pedigree_data_dir, **kwargs):
        self.pedigree_data_dir = Path(pedigree_data_dir)
        self.kwargs = kwargs
        self.calls: list[str] = []
        DummySimulatedPedigree.created_instances.append(self)

    def create_pedigree(self):
        self.calls.append("create")

    def mask_and_corrupt_data(self):
        self.calls.append("mask")

    def run_algorithm(self):
        self.calls.append("run")

    def get_pedigree_statistics(self):
        self.calls.append("stats")
        return {"nodes": 4}

    def get_metrics(self):
        self.calls.append("metrics")
        return {"accuracy": 0.9}


@pytest.fixture(autouse=True)
def clear_dummy_instances():
    DummySimulatedPedigree.created_instances.clear()


def test_simulate_runs_full_pipeline(monkeypatch, tmp_path):
    params = dict(
        p_mask_node=0.3,
        coverage_level=0.5,
        max_candidate_pedigrees=50,
        epsilon=0.2,
        random_seed=7,
    )
    monkeypatch.setattr(run_sampling_experiment, "SimulatedPedigree", DummySimulatedPedigree)
    monkeypatch.chdir(tmp_path)

    stats, metrics = run_sampling_experiment.simulate(**params)

    relative_expected_dir = (
        Path("results")
        / "sampling_experiment"
        / "pedigree_data"
        / "max_candidate_pedigrees=50_epsilon=0.2"
        / "pedigree7"
    )
    absolute_expected_dir = tmp_path / relative_expected_dir
    assert absolute_expected_dir.exists()

    instance = DummySimulatedPedigree.created_instances[-1]
    assert instance.pedigree_data_dir == relative_expected_dir
    assert instance.kwargs == params
    assert instance.calls == ["create", "mask", "run", "stats", "metrics"]

    assert stats == {"nodes": 4}
    assert metrics == {"accuracy": 0.9}
