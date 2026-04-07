from pathlib import Path

import svgutils.compose as sc
import svgutils.transform as st
from lxml import etree


def stitch_pair(top_path: Path, bottom_path: Path, output_path: Path) -> None:
    top = sc.SVG(str(top_path))
    bottom = sc.SVG(str(bottom_path))
    # st.fromfile preserves original pt dimensions; sc.SVG inflates them to px
    top_fig = st.fromfile(str(top_path))
    bottom_fig = st.fromfile(str(bottom_path))
    top_w = float(top_fig.width.replace("pt", ""))
    top_h = float(top_fig.height.replace("pt", ""))
    bottom_w = float(bottom_fig.width.replace("pt", ""))
    bottom_h = float(bottom_fig.height.replace("pt", ""))

    label_size = 20
    gap = 16
    bottom_y = top_h + gap

    fig = sc.Figure(
        f"{max(top_w, bottom_w)}pt",
        f"{bottom_y + bottom_h}pt",
        top.move(0, 0),
        bottom.move(0, bottom_y),
        sc.Text("A", 5, label_size, size=label_size, weight="bold", font="sans-serif"),
        sc.Text("B", 5, bottom_y + label_size, size=label_size, weight="bold", font="sans-serif"),
    )
    # Add new white background so the gap between panels isn't transparent
    bg = etree.Element("rect", width="100%", height="100%", fill="white")
    fig.root.insert(0, bg)
    fig.save(str(output_path))


def main():
    script_dir = Path(__file__).resolve().parent
    results_dir = script_dir / "results"

    pairs = [
        ("parameter_experiment", "combined_heatmaps.svg"),
        ("sampling_experiment", "combined_heatmaps.svg"),
    ]
    for experiment, output_name in pairs:
        plots_dir = results_dir / experiment / "plots"
        stitch_pair(
            plots_dir / "degree_f1_heatmap.svg",
            plots_dir / "relation_f1_heatmap.svg",
            plots_dir / output_name,
        )


if __name__ == "__main__":
    main()
