from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def load_ibd_table(ibd_path: Path | str) -> pd.DataFrame:
    ibd_path = Path(ibd_path)
    df = pd.read_excel(ibd_path)
    df = df.rename(columns={"sum_IBD>12": "sum_ibd_gt12", "n_IBD>12": "n_ibd_gt12"})
    for col in ("sum_ibd_gt12", "n_ibd_gt12"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df[["iid1", "iid2", "sum_ibd_gt12", "n_ibd_gt12"]]


def categories_from_constraints_and_degrees(
    rel_df: pd.DataFrame, *, relation_to_category: dict[str, str], unknown_category: str
) -> pd.Series:
    constraints = rel_df["constraints"].fillna("").astype(str).str.strip()
    categories = constraints.map(relation_to_category)
    degrees = pd.to_numeric(rel_df["degree"], errors="coerce")
    degree_labels = degrees.map(lambda val: str(int(val)) if pd.notna(val) else "")
    degree_categories = degree_labels.map(relation_to_category)
    return categories.fillna(degree_categories).fillna(unknown_category)


def pair_categories(
    relations_path: Path | str, *, published: bool, relation_to_category: dict[str, str], unknown_category: str
) -> dict[tuple[str, str], str]:
    rel_df = pd.read_csv(relations_path, comment="#" if published else None, dtype=str)
    if published:
        rel_df["category"] = categories_from_constraints_and_degrees(
            rel_df, relation_to_category=relation_to_category, unknown_category=unknown_category
        )
    else:
        rel_df["category"] = rel_df["relation"].map(relation_to_category).fillna(unknown_category)

    rel_df["pair_key"] = rel_df[["id1", "id2"]].apply(lambda row: tuple(sorted(row)), axis=1)
    rel_df = rel_df.drop_duplicates(subset="pair_key")
    return dict(zip(rel_df["pair_key"], rel_df["category"], strict=False))


def label_ibd_pairs(
    ibd_df: pd.DataFrame, pair_to_category: dict[tuple[str, str], str], *, unknown_category: str
) -> pd.DataFrame:
    labeled = ibd_df.copy()
    labeled["pair_key"] = labeled[["iid1", "iid2"]].apply(lambda row: tuple(sorted(row)), axis=1)
    labeled["category"] = labeled["pair_key"].map(pair_to_category).fillna(unknown_category)
    return labeled


def plot_pairs(ibd_df: pd.DataFrame, output_path: Path, title: str, *, categories: list[str]) -> None:
    # Plot the 2nd-degree direct lineage last so those points sit on top
    target_category = "2nd-degree related, direct lineage"
    plot_df = (
        ibd_df.assign(_plot_priority=ibd_df["category"].eq(target_category))
        .sort_values("_plot_priority", kind="stable")
        .drop(columns="_plot_priority")
    )

    base_palette = sns.color_palette(n_colors=10)
    palette = {
        "1st-degree related, direct lineage": base_palette[0],
        "1st-degree related": base_palette[2],
        "2nd-degree related, direct lineage": base_palette[1],
        "2nd-degree related": base_palette[4],
        "3rd-degree+/unrelated": "lightgray",
    }

    with mpl.rc_context({"figure.constrained_layout.h_pad": 0.15, "figure.constrained_layout.w_pad": 0.15}):
        fig, ax = plt.subplots(figsize=(9, 7), constrained_layout=True)
        sns.scatterplot(
            data=plot_df,
            x="sum_ibd_gt12",
            y="n_ibd_gt12",
            hue="category",
            hue_order=categories,
            palette=palette,
            s=45,
            alpha=0.8,
            edgecolor="black",
            linewidth=0.4,
            ax=ax,
        )

        ax.set_xlabel("Total Length of IBD segments >12 cM (cM)", fontsize=16)
        ax.set_ylabel("# of IBD Segments >12 cM", fontsize=16)
        ax.set_title(title, fontsize=18)
        ax.tick_params(axis="both", labelsize=14)
        ax.grid(True, linewidth=0.8, alpha=0.4)
        ax.legend(title="Relatedness", loc="lower right", title_fontsize=14, fontsize=12)
        sns.despine(ax=ax)

        fig.savefig(output_path, dpi=600, bbox_inches="tight")
        plt.close(fig)


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    data_dir = script_dir / "data" / "gurgy"
    ibd_path = data_dir / "raw_data" / "rivollat_supplement" / "Supp_tables_zip" / "Supplementary_table_10_IBDs.xlsx"
    published_relations_path = data_dir / "published_exact_relations.csv"
    reconstructed_relations_path = (
        script_dir
        / "results"
        / "evaluation_outputs"
        / "gurgy"
        / "inferred_relations_READv2"
        / "reconstructed_exact_relations.csv"
    )
    output_dir = script_dir / "results" / "gurgy_ancibd"
    output_dir.mkdir(parents=True, exist_ok=True)
    ibd_df = load_ibd_table(ibd_path)
    categories = [
        "1st-degree related, direct lineage",
        "1st-degree related",
        "2nd-degree related, direct lineage",
        "2nd-degree related",
        "3rd-degree+/unrelated",
    ]
    relation_to_category: dict[str, str] = {
        "child-parent": "1st-degree related, direct lineage",
        "parent-child": "1st-degree related, direct lineage",
        "parent-child;child-parent": "1st-degree related, direct lineage",
        "child-parent;parent-child": "1st-degree related, direct lineage",
        "siblings": "1st-degree related",
        "1": "1st-degree related",
        "paternal half-siblings": "2nd-degree related",
        "maternal half-siblings": "2nd-degree related",
        "paternal grandchild-grandparent": "2nd-degree related, direct lineage",
        "paternal grandparent-grandchild": "2nd-degree related, direct lineage",
        "maternal grandchild-grandparent": "2nd-degree related, direct lineage",
        "maternal grandparent-grandchild": "2nd-degree related, direct lineage",
        "maternal grandparent-grandchild;paternal grandparent-grandchild": "2nd-degree related, direct lineage",
        "maternal grandchild-grandparent;paternal grandchild-grandparent": "2nd-degree related, direct lineage",
        "paternal nephew/niece-aunt/uncle": "2nd-degree related",
        "paternal aunt/uncle-nephew/niece": "2nd-degree related",
        "maternal nephew/niece-aunt/uncle": "2nd-degree related",
        "maternal aunt/uncle-nephew/niece": "2nd-degree related",
        "2": "2nd-degree related",
    }
    unknown_category = categories[-1]

    published_pairs = pair_categories(
        published_relations_path,
        published=True,
        relation_to_category=relation_to_category,
        unknown_category=unknown_category,
    )
    published_ibd = label_ibd_pairs(ibd_df, published_pairs, unknown_category=unknown_category)
    plot_pairs(
        published_ibd,
        output_dir / "ancibd_published.pdf",
        "Gurgy ancIBD Pairs (Published Pedigree)",
        categories=categories,
    )

    reconstructed_pairs = pair_categories(
        reconstructed_relations_path,
        published=False,
        relation_to_category=relation_to_category,
        unknown_category=unknown_category,
    )
    reconstructed_ibd = label_ibd_pairs(ibd_df, reconstructed_pairs, unknown_category=unknown_category)
    plot_pairs(
        reconstructed_ibd,
        output_dir / "ancibd_repare.pdf",
        "Gurgy ancIBD Pairs (repare-Reconstructed Pedigree)",
        categories=categories,
    )


if __name__ == "__main__":
    main()
