from benchmarks.published.evaluator import comparison_utils


class StubPedigree:
    def __init__(self, haplogroups: dict[str, str]):
        self.node_to_data = {node: {"mt_haplogroup": haplo} for node, haplo in haplogroups.items()}

    def get_data(self, node: str) -> dict[str, str]:
        return self.node_to_data[node]


def test_get_mt_colormap_uses_union_of_haplogroups():
    inferred = StubPedigree(
        {
            "I1": "H*",
            "I2": "L",
            "101": "ignored",  # Numeric IDs should be skipped
        }
    )
    published = StubPedigree(
        {
            "P1": "J",
            "P2": "H",
        }
    )

    color_map = comparison_utils.get_mt_colormap(inferred, published)

    assert set(color_map.keys()) == {"H", "J", "L"}
    assert "H*" not in color_map
    assert list(color_map.keys()) == sorted(color_map.keys())

    for color in color_map.values():
        assert len(color) == 4
        for channel in color:
            assert 0 <= channel <= 1
