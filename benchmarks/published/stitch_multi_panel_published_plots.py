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
    results_dir = Path(__file__).resolve().parent / "results"

    nepluyesvky_dir = results_dir / "nepluyevsky_comparison"
    stitch_pair(
        nepluyesvky_dir / "published_pedigree_cropped.svg",
        nepluyesvky_dir / "inferred_pedigree_cropped.svg",
        nepluyesvky_dir / "combined_pedigrees_cropped.svg",
    )

    gurgy_dir = results_dir / "gurgy_ancibd"
    stitch_pair(
        gurgy_dir / "ancibd_published.svg",
        gurgy_dir / "ancibd_repare.svg",
        gurgy_dir / "ancibd_combined.svg",
    )


if __name__ == "__main__":
    main()
