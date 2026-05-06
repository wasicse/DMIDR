#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import mdtraj as md
import pandas as pd

# Wilke scale (Å²)
MAX_ASA = {
    'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0, 'CYS': 167.0,
    'GLN': 225.0, 'GLU': 223.0, 'GLY': 104.0, 'HIS': 224.0, 'ILE': 197.0,
    'LEU': 201.0, 'LYS': 236.0, 'MET': 224.0, 'PHE': 240.0, 'PRO': 159.0,
    'SER': 155.0, 'THR': 172.0, 'TRP': 285.0, 'TYR': 263.0, 'VAL': 174.0,
}

THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}


def compute_asa_table(
    pdb_path: str,
    chain_id: str | None,
    exposed_threshold: float,
    model_index: int,
) -> pd.DataFrame:
    traj = md.load(pdb_path)
    if model_index != 0:
        traj = traj[model_index]

    # Shrake-Rupley SASA in nm² per residue → convert to Å²
    sasa_a2 = md.shrake_rupley(traj, mode='residue')[0] * 100.0

    # DSSP secondary structure (full 8-state)
    dssp = md.compute_dssp(traj, simplified=False)[0]

    rows: list[dict] = []
    for i, res in enumerate(traj.topology.residues):
        chain = res.chain.chain_id if res.chain.chain_id else str(res.chain.index)
        if chain_id is not None and chain != chain_id:
            continue
        name3 = res.name
        aa = THREE_TO_ONE.get(name3, 'X')
        max_asa = MAX_ASA.get(name3, 200.0)
        asa = float(sasa_a2[i])
        rasa = asa / max_asa if max_asa else 0.0
        rows.append({
            "chain": chain,
            "resid": res.resSeq,
            "aa": aa,
            "ss": dssp[i],
            "asa": round(asa, 4),
            "rasa": round(rasa, 4),
            "exposure": "Exposed" if rasa > exposed_threshold else "Buried",
        })

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df.sort_values(["chain", "resid"], kind="stable")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Compute per-residue ASA and rASA from a PDB using mdtraj (no external binary required).",
    )
    p.add_argument("--pdb", required=True, help="Path to input PDB file.")
    p.add_argument("--out", required=True, help="Output CSV path.")
    p.add_argument("--chain", default=None, help="Chain ID to filter (e.g., A).")
    p.add_argument("--threshold", type=float, default=0.25,
                   help="Relative ASA threshold for exposure (default: 0.25).")
    p.add_argument("--model-index", type=int, default=0,
                   help="Frame/model index to analyse (default: 0).")
    return p


def main() -> None:
    args = build_parser().parse_args()

    df = compute_asa_table(
        pdb_path=args.pdb,
        chain_id=args.chain,
        exposed_threshold=args.threshold,
        model_index=args.model_index,
    )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)

    if df.empty:
        print("No residues found for the requested selection.")
        return

    exposed = (df["exposure"] == "Exposed").sum()
    buried = (df["exposure"] == "Buried").sum()
    print(f"Wrote ASA table to {out_path}")
    print(f"Residues: {len(df)} | Exposed: {exposed} | Buried: {buried}")
    print(f"Mean rASA: {df['rasa'].mean():.4f}")


if __name__ == "__main__":
    main()
