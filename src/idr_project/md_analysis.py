from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis.rms import RMSF
    from MDAnalysis.lib.distances import distance_array
except Exception as exc:
    raise RuntimeError('MDAnalysis is required for md_analysis.py') from exc

try:
    import mdtraj as md
    MDT_AVAILABLE = True
except Exception:
    MDT_AVAILABLE = False

from .io_utils import load_disorder_probability_table


def compute_rmsf_and_map(universe: mda.Universe, atomsel: str = 'protein and name CA') -> pd.DataFrame:
    ca = universe.select_atoms(atomsel)
    if len(ca) == 0:
        raise ValueError(f"No atoms selected with '{atomsel}'.")
    rmsf = RMSF(ca).run().rmsf
    out = pd.DataFrame({'resid': ca.resids.astype(int), 'rmsf': rmsf.astype(float)})
    return out.groupby('resid', as_index=False)['rmsf'].mean()


def compute_contact_counts(
    universe: mda.Universe,
    cutoff_nm: float = 0.8,
    atomsel: str = 'protein and name CA',
    exclude_neighbors: int = 1,
) -> pd.DataFrame:
    ca = universe.select_atoms(atomsel)
    if len(ca) == 0:
        raise ValueError(f"No atoms selected with '{atomsel}'.")

    n = len(ca)
    idx = np.arange(n)
    i_grid, j_grid = np.meshgrid(idx, idx, indexing='ij')
    mask = (i_grid < j_grid) & (np.abs(i_grid - j_grid) > exclude_neighbors)

    times: list[float] = []
    counts: list[int] = []
    for ts in universe.trajectory:
        pos_nm = ca.positions * 0.1
        dist = distance_array(pos_nm, pos_nm, box=ts.dimensions)
        counts.append(int(np.sum((dist < cutoff_nm) & mask)))
        times.append(float(ts.time))

    return pd.DataFrame({'time': np.array(times), 'contacts': np.array(counts)})


def compute_ss_fractions_with_mdtraj(top_path: str, traj_path: str) -> pd.DataFrame:
    if not MDT_AVAILABLE:
        raise RuntimeError('mdtraj is not available. Install mdtraj to compute DSSP.')
    traj = md.load(traj_path, top=top_path)
    ss = md.compute_dssp(traj, simplified=True)
    time = traj.time if traj.time is not None and len(traj.time) == traj.n_frames else np.arange(traj.n_frames)
    return pd.DataFrame({
        'time': time,
        'frac_H': np.mean(ss == 'H', axis=1),
        'frac_E': np.mean(ss == 'E', axis=1),
        'frac_C': np.mean(ss == 'C', axis=1),
    })


def plot_timeseries(df_wt: pd.DataFrame, df_mut: pd.DataFrame, ycol: str, title: str, outpng: str | Path) -> None:
    plt.figure()
    plt.plot(df_wt['time'], df_wt[ycol], label='WT')
    plt.plot(df_mut['time'], df_mut[ycol], label='Mutant')
    plt.xlabel('Time (ps or frame index)')
    plt.ylabel(ycol)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpng, dpi=200)
    plt.close()


def plot_scatter(
    x_wt: np.ndarray,
    y_wt: np.ndarray,
    x_mut: np.ndarray,
    y_mut: np.ndarray,
    xlabel: str,
    ylabel: str,
    title: str,
    outpng: str | Path,
) -> None:
    plt.figure()
    plt.scatter(x_wt, y_wt, s=18, alpha=0.7, label='WT')
    plt.scatter(x_mut, y_mut, s=18, alpha=0.7, label='Mutant')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpng, dpi=200)
    plt.close()


def spearman_report(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    rho, pvalue = spearmanr(x, y, nan_policy='omit')
    return float(rho), float(pvalue)


def run_md_analysis(
    wt_top: str,
    wt_traj: str,
    mut_top: str,
    mut_traj: str,
    disorder_path: str,
    outdir: str,
    cutoff_nm: float = 0.8,
    exclude_neighbors: int = 1,
) -> None:
    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)

    disorder = load_disorder_probability_table(disorder_path)
    u_wt = mda.Universe(wt_top, wt_traj)
    u_mut = mda.Universe(mut_top, mut_traj)

    ss_ok = True
    try:
        ss_wt = compute_ss_fractions_with_mdtraj(wt_top, wt_traj)
        ss_mut = compute_ss_fractions_with_mdtraj(mut_top, mut_traj)
        ss_wt.to_csv(outdir_path / 'wt_ss_fractions.csv', index=False)
        ss_mut.to_csv(outdir_path / 'mut_ss_fractions.csv', index=False)
        plot_timeseries(ss_wt, ss_mut, 'frac_H', 'Helix fraction vs time', outdir_path / 'ss_frac_H.png')
        plot_timeseries(ss_wt, ss_mut, 'frac_E', 'Beta fraction vs time', outdir_path / 'ss_frac_E.png')
        plot_timeseries(ss_wt, ss_mut, 'frac_C', 'Coil fraction vs time', outdir_path / 'ss_frac_C.png')
    except Exception as exc:
        ss_ok = False
        print(f'[WARN] DSSP calculation failed: {exc}', file=sys.stderr)

    contacts_wt = compute_contact_counts(u_wt, cutoff_nm=cutoff_nm, exclude_neighbors=exclude_neighbors)
    contacts_mut = compute_contact_counts(u_mut, cutoff_nm=cutoff_nm, exclude_neighbors=exclude_neighbors)
    contacts_wt.to_csv(outdir_path / 'wt_contacts.csv', index=False)
    contacts_mut.to_csv(outdir_path / 'mut_contacts.csv', index=False)
    plot_timeseries(
        contacts_wt,
        contacts_mut,
        'contacts',
        f'CA-CA contacts (<{cutoff_nm} nm) vs time',
        outdir_path / 'contacts_vs_time.png',
    )

    rmsf_wt = compute_rmsf_and_map(u_wt)
    rmsf_mut = compute_rmsf_and_map(u_mut)
    merged_wt = pd.merge(disorder, rmsf_wt, on='resid', how='inner')
    merged_mut = pd.merge(disorder, rmsf_mut, on='resid', how='inner')

    rho_wt, p_wt = spearman_report(merged_wt['disorder_prob'].values, merged_wt['rmsf'].values)
    rho_mut, p_mut = spearman_report(merged_mut['disorder_prob'].values, merged_mut['rmsf'].values)

    merged_wt.to_csv(outdir_path / 'wt_disorder_rmsf.csv', index=False)
    merged_mut.to_csv(outdir_path / 'mut_disorder_rmsf.csv', index=False)
    plot_scatter(
        merged_wt['disorder_prob'].values,
        merged_wt['rmsf'].values,
        merged_mut['disorder_prob'].values,
        merged_mut['rmsf'].values,
        'Disorder probability',
        'RMSF (Å)',
        'Disorder probability vs RMSF',
        outdir_path / 'disorder_vs_rmsf.png',
    )

    with (outdir_path / 'summary.txt').open('w', encoding='utf-8') as handle:
        handle.write('=== MD post-processing summary ===\n\n')
        handle.write(f'WT: rho={rho_wt:.3f}, p={p_wt:.3e}\n')
        handle.write(f'Mutant: rho={rho_mut:.3f}, p={p_mut:.3e}\n\n')
        handle.write(f'Contacts cutoff (nm): {cutoff_nm}\n')
        handle.write(f'Exclude neighbors: {exclude_neighbors}\n')
        handle.write(f"WT mean contacts: {contacts_wt['contacts'].mean():.3f}\n")
        handle.write(f"Mutant mean contacts: {contacts_mut['contacts'].mean():.3f}\n")
        if ss_ok:
            handle.write('\nSecondary structure means:\n')
            handle.write(f"WT H={ss_wt['frac_H'].mean():.3f}, E={ss_wt['frac_E'].mean():.3f}, C={ss_wt['frac_C'].mean():.3f}\n")
            handle.write(f"Mut H={ss_mut['frac_H'].mean():.3f}, E={ss_mut['frac_E'].mean():.3f}, C={ss_mut['frac_C'].mean():.3f}\n")
        else:
            handle.write('\nSecondary structure was skipped or failed.\n')
