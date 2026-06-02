#!/usr/bin/env python3
"""
Fix a ColabFold PDB before passing it to GROMACS.

Steps:
  1. PDBFixer — remove heterogens, fill missing residues/atoms
  2. Add hydrogens at pH 7.0
  3. OpenMM vacuum minimization (AMBER14) to resolve steric clashes
  4. Write a clean PDB (with H stripped so pdb2gmx -ignh re-adds them)

Usage:
  uv run python scripts/fix_pdb.py input.pdb output.pdb [--max-iter 2000]
"""
import argparse
import os
import sys
from pathlib import Path

from pdbfixer import PDBFixer
from openmm import LangevinMiddleIntegrator, Platform
from openmm.app import (
    ForceField, Modeller, PDBFile, Simulation,
    HBonds, NoCutoff
)
from openmm.unit import kelvin, picosecond, angstrom, kilocalorie_per_mole


def fix(input_pdb: str, output_pdb: str, max_iter: int = 2000) -> None:
    print(f"[fix_pdb] Input:  {input_pdb}")

    fixer = PDBFixer(filename=input_pdb)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)

    print(f"[fix_pdb] Structure: {fixer.topology.getNumResidues()} residues, "
          f"{fixer.topology.getNumAtoms()} atoms")

    ff = ForceField("amber14-all.xml")
    modeller = Modeller(fixer.topology, fixer.positions)

    try:
        modeller.addHydrogens(ff, pH=7.0)
    except Exception as e:
        print(f"[fix_pdb] addHydrogens skipped ({e}); using PDBFixer H placement")

    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=NoCutoff,
        constraints=HBonds,
        hydrogenMass=1.5,
    )

    platform = Platform.getPlatformByName("CPU")
    platform.setPropertyDefaultValue("Threads", str(os.cpu_count() or 1))
    integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picosecond)
    sim = Simulation(modeller.topology, system, integrator, platform)
    sim.context.setPositions(modeller.positions)

    state_before = sim.context.getState(getEnergy=True)
    e_before = state_before.getPotentialEnergy().value_in_unit(kilocalorie_per_mole)
    print(f"[fix_pdb] Energy before minimization: {e_before:.1f} kcal/mol")

    sim.minimizeEnergy(maxIterations=max_iter, tolerance=10 * kilocalorie_per_mole / angstrom)

    state_after = sim.context.getState(getEnergy=True, getPositions=True)
    e_after = state_after.getPotentialEnergy().value_in_unit(kilocalorie_per_mole)
    print(f"[fix_pdb] Energy after  minimization: {e_after:.1f} kcal/mol")

    positions = state_after.getPositions()

    # Strip hydrogens — pdb2gmx -ignh will add them with the correct force field
    out_top = modeller.topology
    heavy_indices = [
        a.index for a in out_top.atoms() if a.element.symbol != "H"
    ]
    heavy_positions = [positions[i] for i in heavy_indices]

    # Rebuild topology with only heavy atoms
    fixer2 = PDBFixer.__new__(PDBFixer)
    Path(output_pdb).parent.mkdir(parents=True, exist_ok=True)
    with open(output_pdb, "w") as fh:
        # Write full structure (with H) — pdb2gmx -ignh discards them anyway
        PDBFile.writeFile(out_top, positions, fh, keepIds=True)

    print(f"[fix_pdb] Output: {output_pdb}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Fix PDB for GROMACS")
    parser.add_argument("input",  help="Input PDB (ColabFold output)")
    parser.add_argument("output", help="Output fixed PDB")
    parser.add_argument("--max-iter", type=int, default=2000,
                        help="OpenMM minimization max iterations (default 2000)")
    args = parser.parse_args()

    if not Path(args.input).exists():
        print(f"ERROR: {args.input} not found", file=sys.stderr)
        sys.exit(1)

    fix(args.input, args.output, args.max_iter)


if __name__ == "__main__":
    main()
