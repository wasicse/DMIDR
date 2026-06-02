#!/usr/bin/env python3
"""Generate four summary figures.  Exact plotting code ported from the old repo.

Figure 1 — label_transitions_by_block.png
    Average D→O and O→D counts per mutant, grouped by mutation block size.
    Source: AnalyzeResults.ipynb (MutationRandom) pattern applied to .caid files.

Figure 2 — avg_disorder_change_by_block.png
    Average disorder probability change per mutant, grouped by mutation block size.
    Same source.

Figure 3 — asa_delta_do_residues.png
    ΔASA (Å²) for residues that flipped Disordered→Ordered.
    Exact code: Final.ipynb Cell 1.

Figure 4 — do_ss_sequence_track.png
    Letter-level DO / SS / AA track for D→O residues.
    Exact code: Final.ipynb Cell 3 / comp.ipynb Cell 16.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import seaborn as sns

DISORDER_THRESHOLD = 0.5


# ─────────────────────────────────────────────────────────────────────────────
# Helpers — parse .caid files (DMIDR format)
# ─────────────────────────────────────────────────────────────────────────────

def _read_wt_caid(path: Path) -> pd.DataFrame:
    rows: list[list] = []
    with path.open(encoding="utf-8", errors="ignore") as fh:
        seen = False
        for line in fh:
            if line.startswith(">"):
                if seen:
                    break
                seen = True
                continue
            parts = line.strip().split()
            if len(parts) >= 3:
                prob = float(parts[2])
                flag = int(parts[3]) if len(parts) >= 4 else int(prob >= DISORDER_THRESHOLD)
                rows.append([int(parts[0]), parts[1], prob, flag])
    return pd.DataFrame(rows, columns=["SeqNo", "AA", "DisorderProb", "Disordered"])


def _parse_all_mutants_caid(path: Path) -> dict[str, pd.DataFrame]:
    """Parse every >header block from a multi-mutant .caid file."""
    mutants: dict[str, pd.DataFrame] = {}
    current: str | None = None
    rows: list[list] = []

    with path.open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                if current is not None and rows:
                    mutants[current] = pd.DataFrame(
                        rows, columns=["SeqNo", "AA", "DisorderProb", "Disordered"]
                    )
                current = line[1:].strip()
                rows = []
            elif line.strip():
                parts = line.strip().split()
                if len(parts) >= 3:
                    prob = float(parts[2])
                    flag = int(parts[3]) if len(parts) >= 4 else int(prob >= DISORDER_THRESHOLD)
                    rows.append([int(parts[0]), parts[1], prob, flag])
        if current is not None and rows:
            mutants[current] = pd.DataFrame(
                rows, columns=["SeqNo", "AA", "DisorderProb", "Disordered"]
            )
    return mutants


def _load_best_asa(asa_dir: Path) -> pd.DataFrame | None:
    """Return the rank_001 ASA CSV from a structure's ASA output directory."""
    hits = list(asa_dir.glob("*rank_001*_asa.csv")) or list(asa_dir.glob("*_asa.csv"))
    if not hits:
        return None
    df = pd.read_csv(hits[0])
    if "resid" in df.columns and "Residue" not in df.columns:
        df = df.rename(columns={"resid": "Residue"})
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Figure 1 — Label Transitions by Mutation Group
# ─────────────────────────────────────────────────────────────────────────────

def plot_label_transitions(summary_df: pd.DataFrame, outdir: Path) -> None:
    """Exact style from MutationRandom analysis (grouped bar, blue + gray)."""
    grouped = (
        summary_df
        .groupby("SourceFile")[["DisorderedToOrdered", "OrderedToDisordered"]]
        .mean()
        .sort_index()
    )

    groups = grouped.index.tolist()
    do_vals = grouped["DisorderedToOrdered"].tolist()
    od_vals = grouped["OrderedToDisordered"].tolist()
    x = np.arange(len(groups))
    width = 0.35

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")

    bars_do = ax.bar(x - width / 2, do_vals, width,
                     label="Disordered → Ordered", color="#1f77b4", edgecolor="black")
    bars_od = ax.bar(x + width / 2, od_vals, width,
                     label="Ordered → Disordered", color="#7f7f7f", edgecolor="black")

    for bar, val in zip(bars_do, do_vals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                f"{val:.1f}", ha="center", va="bottom", fontsize=10)
    for bar, val in zip(bars_od, od_vals):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.3,
                f"{val:.1f}", ha="center", va="bottom", fontsize=10)

    ax.set_xlabel("Mutation Group", fontsize=12)
    ax.set_ylabel("Average Count", fontsize=12)
    ax.set_title("Label Transitions by Mutation Group", fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(groups, rotation=30, ha="right")
    ax.legend()
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    fig.tight_layout()

    out = outdir / "label_transitions_by_block.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"  Saved: {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 2 — Average Disorder Change by Mutation Group
# ─────────────────────────────────────────────────────────────────────────────

def plot_avg_disorder_change(summary_df: pd.DataFrame, outdir: Path) -> None:
    """Exact style from MutationRandom analysis (single blue bar, neg values)."""
    grouped = (
        summary_df
        .groupby("SourceFile")["AvgDisorderChange"]
        .mean()
        .sort_index()
    )

    groups = grouped.index.tolist()
    vals   = grouped.tolist()

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")

    bars = ax.bar(groups, vals, color="#1f77b4", edgecolor="black")

    for bar, val in zip(bars, vals):
        label_y = bar.get_height() - 0.003 if val < 0 else bar.get_height() + 0.001
        va = "top" if val < 0 else "bottom"
        ax.text(bar.get_x() + bar.get_width() / 2, label_y,
                f"{val:.3f}", ha="center", va=va, fontsize=10)

    ax.set_xlabel("Mutation Group", fontsize=12)
    ax.set_ylabel("AvgDisorderChange", fontsize=12)
    ax.set_title("Average Disorder Change by Mutation Group", fontsize=14)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, rotation=30, ha="right")
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    fig.tight_layout()

    out = outdir / "avg_disorder_change_by_block.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"  Saved: {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 3 — ΔASA for D→O residues
# Exact code: Final.ipynb Cell 1
# ─────────────────────────────────────────────────────────────────────────────

def plot_asa_delta_do(result_df: pd.DataFrame, outdir: Path) -> None:
    df = result_df.sort_values("Residue").reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(20, 10))

    ax.plot(df["Residue"], df["ASA_Change"],
            marker="o", markersize=8, linestyle="-",
            color="darkblue", linewidth=3, label="ΔASA")

    ax.axhline(0, color="gray", linestyle="--", linewidth=1.5)

    top_n = 5
    highlighted = pd.concat([
        df.nlargest(top_n, "ASA_Change"),
        df.nsmallest(top_n, "ASA_Change"),
    ])
    for _, row in highlighted.iterrows():
        ax.annotate(f"{int(row['Residue'])}",
                    (row["Residue"], row["ASA_Change"]),
                    textcoords="offset points", xytext=(0, 12), ha="center",
                    fontsize=14, fontweight="bold", color="darkred")

    ax.set_title("Accessible Surface Area (ASA) Changes in Disordered→Ordered Residues",
                 fontsize=20, fontweight="bold")
    ax.set_xlabel("Residue Position", fontsize=16, fontweight="bold")
    ax.set_ylabel("ΔASA (Å²)", fontsize=16, fontweight="bold")
    ax.grid(True, linestyle=":", alpha=0.7)
    ax.tick_params(axis="x", labelrotation=45, labelsize=14)
    ax.tick_params(axis="y", labelsize=14)
    fig.tight_layout()

    out = outdir / "asa_delta_do_residues.png"
    fig.savefig(out, dpi=400)
    plt.close(fig)
    print(f"  Saved: {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 4 — DO / SS / AA sequence letter track
# Exact code: Final.ipynb Cell 3
# ─────────────────────────────────────────────────────────────────────────────

def plot_sequence_track(df: pd.DataFrame, outdir: Path) -> None:
    """df must have columns: Residue, AA, SS_Before, SS_After, From, To."""
    df = df.copy()
    df["DO_Before"]  = df["From"]
    df["DO_After"]   = df["To"]
    df["Transition"] = df["From"] + "→" + df["To"]

    # Color maps
    transition_palette = sns.color_palette("husl", df["Transition"].nunique())
    transition_color_map = dict(zip(df["Transition"].unique(), transition_palette))

    aa_palette = sns.color_palette("hls", df["AA"].nunique())
    aa_color_map = dict(zip(df["AA"].unique(), aa_palette))

    fig, ax = plt.subplots(figsize=(20, 6))
    ax.axis("off")

    rows = ["DO_After", "DO_Before", "SS_After", "SS_Before", "AA", "Residue"]
    row_height     = 1.0
    font_size      = 16
    res_font_size  = 11
    rotation_angle = 90

    y_positions = {r: (len(rows) - i) * row_height for i, r in enumerate(rows)}

    for col_idx, row in df.iterrows():
        for row_name in rows:
            value    = row[row_name]
            rotation = rotation_angle if row_name == "Residue" else 0

            if row_name in ("DO_After", "DO_Before"):
                color = transition_color_map[row["Transition"]]
            elif row_name == "AA":
                color = aa_color_map[row["AA"]]
            elif row_name in ("SS_Before", "SS_After"):
                color = "red" if row["SS_Before"] != row["SS_After"] else "black"
            else:
                color = "black"

            fs = res_font_size if row_name == "Residue" else font_size
            fw = "bold" if row_name == "Residue" else None
            effects = ([pe.withStroke(linewidth=2, foreground="white")]
                       if row_name == "Residue" else None)

            ax.text(col_idx + 1, y_positions[row_name], str(value),
                    ha="center", va="center", fontsize=fs, fontweight=fw,
                    color=color, rotation=rotation, family="monospace",
                    path_effects=effects)

    for row_name, y in y_positions.items():
        ax.text(0, y, row_name, ha="right", va="center",
                fontsize=font_size + 2, fontweight="bold", family="monospace")

    ax.set_xlim(0, len(df) + 1)
    ax.set_ylim(0, max(y_positions.values()) + 1)
    fig.tight_layout()

    out_svg = outdir / "do_ss_sequence_track.svg"
    out_png = outdir / "do_ss_sequence_track.png"
    fig.savefig(out_svg, format="svg", dpi=400)
    fig.savefig(out_png, dpi=400)
    plt.close(fig)
    print(f"  Saved: {out_png} + {out_svg}")


# ─────────────────────────────────────────────────────────────────────────────
# Build result_df (needed by figures 3 & 4) from pipeline ASA + label CSVs
# ─────────────────────────────────────────────────────────────────────────────

def _build_result_df(label_csv: Path, wt_asa: pd.DataFrame,
                     mut_asa: pd.DataFrame) -> pd.DataFrame:
    """Produce a result_df with the same columns as the old compare_structure_features()."""
    label_df = pd.read_csv(label_csv)

    # Keep only D→O and O→D residues (matching old keep_all=False behaviour)
    changed  = label_df[label_df["LabelChange"].isin([-1, 1])].copy()

    # ASA CSVs use 'resid' → renamed to 'Residue' by _load_best_asa.
    # Label CSV uses 'Position'. Normalise both to 'Position' for the merge.
    wt_m  = wt_asa.rename(columns={"Residue": "Position", "ss": "SS_Before",
                                    "asa": "ASA_Before", "aa": "AA"})
    mut_m = mut_asa.rename(columns={"Residue": "Position", "ss": "SS_After",
                                     "asa": "ASA_After"})

    merged = (
        changed
        .merge(wt_m[["Position", "AA", "SS_Before", "ASA_Before"]], on="Position", how="left")
        .merge(mut_m[["Position", "SS_After", "ASA_After"]],         on="Position", how="left")
    )

    merged["ASA_Change"] = merged["ASA_After"] - merged["ASA_Before"]
    merged["SS_Before"]  = merged["SS_Before"].fillna("-").str.strip().replace("", "-")
    merged["SS_After"]   = merged["SS_After"].fillna("-").str.strip().replace("", "-")
    # AA comes from the ASA merge; fall back to AA_orig from the label CSV
    if "AA" not in merged.columns or merged["AA"].isna().all():
        merged["AA"] = merged.get("AA_orig", pd.Series("-", index=merged.index))
    merged["AA"] = merged["AA"].fillna("-")
    merged["Residue"] = merged["Position"]  # alias for plot functions

    from_to_map = {-1: ("D", "O"), 1: ("O", "D")}
    merged["From"] = merged["LabelChange"].map(lambda x: from_to_map.get(x, ("-", "-"))[0])
    merged["To"]   = merged["LabelChange"].map(lambda x: from_to_map.get(x, ("-", "-"))[1])

    return merged.reset_index(drop=True)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate 4 summary figures (exact code from old repo)."
    )
    parser.add_argument("--seq-name",    required=True)
    parser.add_argument("--results-dir", required=True,
                        help="outputs/{seq_name}/ directory")
    parser.add_argument("--dispred-dir", required=True,
                        help="dispred/ folder with .caid files")
    parser.add_argument("--max-block",   type=int, default=5)
    parser.add_argument("--outdir",      required=True)
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    dispred_dir = Path(args.dispred_dir)
    outdir      = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    wt_caid = dispred_dir / f"{args.seq_name}_original.caid"
    if not wt_caid.exists():
        print(f"ERROR: WT .caid not found: {wt_caid}", file=sys.stderr)
        sys.exit(1)

    wt_df = _read_wt_caid(wt_caid)

    # ── Figures 1 & 2: per-mutant stats grouped by block size ─────────────────
    summary_rows: list[dict] = []
    for block in range(1, args.max_block + 1):
        key  = f"Mutants_{block}_Res"
        caid = dispred_dir / f"{args.seq_name}_mutants_{block}res.caid"
        if not caid.exists():
            continue
        all_mutants = _parse_all_mutants_caid(caid)
        print(f"  Block {block}res: {len(all_mutants)} mutants")

        for mut_name, mut_df in all_mutants.items():
            merged = pd.merge(wt_df, mut_df, on="SeqNo", suffixes=("_orig", "_mut"))
            merged["DeltaDisorder"] = merged["DisorderProb_mut"] - merged["DisorderProb_orig"]
            merged["LabelChange"]   = merged["Disordered_mut"] - merged["Disordered_orig"]

            summary_rows.append({
                "SourceFile":          key,
                "Mutant":              mut_name,
                "AvgDisorderChange":   round(float(merged["DeltaDisorder"].mean()), 4),
                "DisorderedToOrdered": int((merged["LabelChange"] == -1).sum()),
                "OrderedToDisordered": int((merged["LabelChange"] ==  1).sum()),
            })

    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_csv = outdir / "all_mutant_summaries.csv"
        summary_df.to_csv(summary_csv, index=False)
        print(f"  Saved: {summary_csv}")
        plot_label_transitions(summary_df, outdir)
        plot_avg_disorder_change(summary_df, outdir)
    else:
        print("  [WARN] No mutant .caid files — skipping figures 1 & 2.")

    # ── Figures 3 & 4: need ASA data ─────────────────────────────────────────
    label_csv = results_dir / f"disorder_1res/{args.seq_name}_1res_disorder_label_comparison.csv"
    if not label_csv.exists():
        print(f"  [WARN] {label_csv} missing — skipping figures 3 & 4.")
        return

    wt_asa_dir    = results_dir / "asa"
    cand1_dirs    = sorted(results_dir.glob("asa_candidate_1_*"))
    cand1_asa_dir = cand1_dirs[0] if cand1_dirs else None

    wt_asa  = _load_best_asa(wt_asa_dir)  if wt_asa_dir.exists()  else None
    mut_asa = _load_best_asa(cand1_asa_dir) if cand1_asa_dir else None

    if wt_asa is None or mut_asa is None:
        print("  [WARN] ASA CSVs missing — skipping figures 3 & 4.")
        return

    result_df = _build_result_df(label_csv, wt_asa, mut_asa)

    # Figure 3 — D→O only (matching old keep_all=False, From=="D")
    do_df = result_df[result_df["From"] == "D"].copy()
    if not do_df.empty:
        plot_asa_delta_do(do_df, outdir)
    else:
        print("  [WARN] No D→O residues — skipping figure 3.")

    # Figure 4 — both D→O and O→D (matching old code that used full result_df)
    if not result_df.empty:
        plot_sequence_track(result_df, outdir)
    else:
        print("  [WARN] Empty result_df — skipping figure 4.")


if __name__ == "__main__":
    main()
