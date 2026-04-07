#!/usr/bin/env python3
"""
Analyze MD trajectories (WT and Mutant) for:
1) Secondary structure fractions vs time (DSSP)  [mdtraj + DSSP]
2) Intramolecular contact count vs time (Cα-Cα cutoff) [MDAnalysis]
3) Disorder probability vs RMSF correlation (Spearman) [MDAnalysis + scipy]

Inputs expected:
- Topology file (e.g., .pdb, .gro) for WT and mutant (can be same if residue numbering matches)
- Trajectory file (e.g., .xtc, .dcd, .trr) for WT and mutant
- Per-residue disorder file (CSV/TSV) with columns:
    resid, disorder_prob
  where resid matches the topology residue IDs (NOT 0-based index).

Outputs:
- CSV files with time series
- Summary txt
- PNG plots
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import spearmanr

import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSF
from MDAnalysis.lib.distances import distance_array

# Secondary structure: mdtraj + DSSP
try:
    import mdtraj as md
    MDT_AVAILABLE = True
except Exception:
    MDT_AVAILABLE = False


def read_disorder_table(path: str) -> pd.DataFrame:
    # Auto-detect delimiter
    with open(path, "r", encoding="utf-8") as f:
        head = f.readline()
    sep = "," if "," in head else ("\t" if "\t" in head else None)
    df = pd.read_csv(path, sep=sep)
    required = {"resid", "disorder_prob"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Disorder file missing columns {missing}. Found: {list(df.columns)}")
    df = df.copy()
    df["resid"] = df["resid"].astype(int)
    df["disorder_prob"] = df["disorder_prob"].astype(float)
    return df.sort_values("resid")


def compute_rmsf_and_map(universe: mda.Universe, atomsel: str = "protein and name CA"):
    """
    RMSF on CA atoms, aligned to first frame implicitly by using positions relative to average.
    Returns a DataFrame with columns: resid, rmsf
    """
    ca = universe.select_atoms(atomsel)
    if len(ca) == 0:
        raise ValueError(f"No atoms selected with '{atomsel}'. Check topology atom names.")

    # RMSF in MDAnalysis uses average structure as reference by default
    rmsf = RMSF(ca).run().rmsf
    resids = ca.resids  # residue IDs from topology
    out = pd.DataFrame({"resid": resids.astype(int), "rmsf": rmsf.astype(float)})
    # If multiple CA per resid (rare), average
    out = out.groupby("resid", as_index=False)["rmsf"].mean()
    return out


def compute_contact_counts(universe: mda.Universe,
                           cutoff_nm: float = 0.8,
                           atomsel: str = "protein and name CA",
                           exclude_neighbors: int = 1):
    """
    Count intramolecular CA-CA contacts per frame.
    cutoff_nm: contact threshold in nm.
    exclude_neighbors: ignore |i-j| <= exclude_neighbors to avoid counting bonded neighbors.
    Returns DataFrame: time_ps, contacts
    """
    ca = universe.select_atoms(atomsel)
    n = len(ca)
    if n == 0:
        raise ValueError(f"No atoms selected with '{atomsel}'.")

    # precompute index pairs mask for i<j and neighbor exclusion
    idx = np.arange(n)
    I, J = np.meshgrid(idx, idx, indexing="ij")
    upper = I < J
    neighbor_ok = (np.abs(I - J) > exclude_neighbors)
    mask = upper & neighbor_ok

    times = []
    counts = []

    for ts in universe.trajectory:
        # positions in Angstrom in MDAnalysis; convert to nm
        pos = ca.positions * 0.1  # Å -> nm
        d = distance_array(pos, pos, box=ts.dimensions)  # nm if positions are nm
        # Count distances below cutoff using mask
        c = np.sum((d < cutoff_nm) & mask)
        times.append(ts.time)  # usually ps; depends on trajectory
        counts.append(int(c))

    return pd.DataFrame({"time": np.array(times, dtype=float),
                         "contacts": np.array(counts, dtype=int)})


def compute_ss_fractions_with_mdtraj(top_path: str, traj_path: str):
    """
    DSSP secondary structure using mdtraj.compute_dssp.
    Returns DataFrame: time, frac_H, frac_E, frac_C
    """
    if not MDT_AVAILABLE:
        raise RuntimeError("mdtraj not available. Install with: pip install mdtraj")

    t = md.load(traj_path, top=top_path)

    # DSSP requires a DSSP executable installed (mkdssp or dssp).
    # mdtraj will error if DSSP is missing.
    ss = md.compute_dssp(t, simplified=True)  # returns array [n_frames, n_res], values in {H,E,C}
    # Fractions per frame
    frac_H = np.mean(ss == "H", axis=1)
    frac_E = np.mean(ss == "E", axis=1)
    frac_C = np.mean(ss == "C", axis=1)

    # time in ps if present, else frames
    if t.time is not None and len(t.time) == t.n_frames:
        time = t.time
    else:
        time = np.arange(t.n_frames, dtype=float)

    return pd.DataFrame({"time": time, "frac_H": frac_H, "frac_E": frac_E, "frac_C": frac_C})


def spearman_report(x, y):
    rho, p = spearmanr(x, y, nan_policy="omit")
    return float(rho), float(p)


def plot_timeseries(df_wt, df_mut, ycol, title, outpng):
    plt.figure()
    plt.plot(df_wt["time"], df_wt[ycol], label="WT")
    plt.plot(df_mut["time"], df_mut[ycol], label="Mutant")
    plt.xlabel("Time (ps or frame index)")
    plt.ylabel(ycol)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpng, dpi=200)
    plt.close()


def plot_scatter(x_wt, y_wt, x_mut, y_mut, xlabel, ylabel, title, outpng):
    plt.figure()
    plt.scatter(x_wt, y_wt, s=18, alpha=0.7, label="WT")
    plt.scatter(x_mut, y_mut, s=18, alpha=0.7, label="Mutant")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpng, dpi=200)
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--wt_top", required=True, help="WT topology (pdb/gro)")
    ap.add_argument("--wt_traj", required=True, help="WT trajectory (xtc/dcd/trr)")
    ap.add_argument("--mut_top", required=True, help="Mutant topology (pdb/gro)")
    ap.add_argument("--mut_traj", required=True, help="Mutant trajectory (xtc/dcd/trr)")
    ap.add_argument("--disorder", required=True, help="Per-residue disorder CSV/TSV with columns resid,disorder_prob")
    ap.add_argument("--outdir", default="md_postproc_out", help="Output directory")
    ap.add_argument("--cutoff_nm", type=float, default=0.8, help="Contact cutoff in nm (default 0.8)")
    ap.add_argument("--exclude_neighbors", type=int, default=1, help="Exclude |i-j| <= N neighbors (default 1)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    disorder = read_disorder_table(args.disorder)

    # Load universes
    u_wt = mda.Universe(args.wt_top, args.wt_traj)
    u_mut = mda.Universe(args.mut_top, args.mut_traj)

    # 1) Secondary structure fractions vs time (DSSP)
    ss_wt = None
    ss_mut = None
    ss_ok = True
    try:
        ss_wt = compute_ss_fractions_with_mdtraj(args.wt_top, args.wt_traj)
        ss_mut = compute_ss_fractions_with_mdtraj(args.mut_top, args.mut_traj)
        ss_wt.to_csv(os.path.join(args.outdir, "wt_ss_fractions.csv"), index=False)
        ss_mut.to_csv(os.path.join(args.outdir, "mut_ss_fractions.csv"), index=False)

        plot_timeseries(ss_wt, ss_mut, "frac_H", "Helix fraction vs time", os.path.join(args.outdir, "ss_frac_H.png"))
        plot_timeseries(ss_wt, ss_mut, "frac_E", "Beta fraction vs time", os.path.join(args.outdir, "ss_frac_E.png"))
        plot_timeseries(ss_wt, ss_mut, "frac_C", "Coil fraction vs time", os.path.join(args.outdir, "ss_frac_C.png"))
    except Exception as e:
        ss_ok = False
        print(f"[WARN] Secondary structure (DSSP) failed: {e}", file=sys.stderr)
        print("[WARN] Install DSSP (mkdssp) and mdtraj to enable SS analysis.", file=sys.stderr)

    # 2) Contact count vs time
    contacts_wt = compute_contact_counts(u_wt, cutoff_nm=args.cutoff_nm, exclude_neighbors=args.exclude_neighbors)
    contacts_mut = compute_contact_counts(u_mut, cutoff_nm=args.cutoff_nm, exclude_neighbors=args.exclude_neighbors)
    contacts_wt.to_csv(os.path.join(args.outdir, "wt_contacts.csv"), index=False)
    contacts_mut.to_csv(os.path.join(args.outdir, "mut_contacts.csv"), index=False)
    plot_timeseries(contacts_wt, contacts_mut, "contacts",
                    f"CA-CA contacts (<{args.cutoff_nm} nm) vs time",
                    os.path.join(args.outdir, "contacts_vs_time.png"))

    # 3) Disorder vs RMSF correlation
    rmsf_wt = compute_rmsf_and_map(u_wt)
    rmsf_mut = compute_rmsf_and_map(u_mut)

    # Merge with disorder by resid (inner join)
    m_wt = pd.merge(disorder, rmsf_wt, on="resid", how="inner")
    m_mut = pd.merge(disorder, rmsf_mut, on="resid", how="inner")

    # Spearman correlation
    rho_wt, p_wt = spearman_report(m_wt["disorder_prob"].values, m_wt["rmsf"].values)
    rho_mut, p_mut = spearman_report(m_mut["disorder_prob"].values, m_mut["rmsf"].values)

    m_wt.to_csv(os.path.join(args.outdir, "wt_disorder_rmsf.csv"), index=False)
    m_mut.to_csv(os.path.join(args.outdir, "mut_disorder_rmsf.csv"), index=False)

    plot_scatter(
        m_wt["disorder_prob"].values, m_wt["rmsf"].values,
        m_mut["disorder_prob"].values, m_mut["rmsf"].values,
        xlabel="Disorder probability", ylabel="RMSF (Å)",
        title="Disorder probability vs RMSF (per residue)",
        outpng=os.path.join(args.outdir, "disorder_vs_rmsf.png")
    )

    # Summary report
    summary_path = os.path.join(args.outdir, "summary.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("=== Computational post-processing summary ===\n\n")
        f.write(f"WT:  Spearman(disorder, RMSF) rho={rho_wt:.3f}, p={p_wt:.3e}\n")
        f.write(f"Mut: Spearman(disorder, RMSF) rho={rho_mut:.3f}, p={p_mut:.3e}\n\n")
        f.write(f"Contacts cutoff: {args.cutoff_nm} nm; exclude_neighbors={args.exclude_neighbors}\n")
        f.write(f"WT contacts mean={contacts_wt['contacts'].mean():.2f}, sd={contacts_wt['contacts'].std():.2f}\n")
        f.write(f"Mut contacts mean={contacts_mut['contacts'].mean():.2f}, sd={contacts_mut['contacts'].std():.2f}\n\n")
        if ss_ok:
            f.write("Secondary structure fractions (mean over all frames):\n")
            f.write(f"WT:  H={ss_wt['frac_H'].mean():.3f}, E={ss_wt['frac_E'].mean():.3f}, C={ss_wt['frac_C'].mean():.3f}\n")
            f.write(f"Mut: H={ss_mut['frac_H'].mean():.3f}, E={ss_mut['frac_E'].mean():.3f}, C={ss_mut['frac_C'].mean():.3f}\n")
        else:
            f.write("Secondary structure: DSSP not computed (see warnings).\n")

    print(f"Done. Outputs written to: {args.outdir}")
    print(f"- {summary_path}")


if __name__ == "__main__":
    main()