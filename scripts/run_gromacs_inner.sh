#!/usr/bin/env bash
# Runs inside the gmx-gpu.sif Apptainer container.
# Usage: bash run_gromacs_inner.sh <wt.pdb> <mut.pdb> <outdir> <mdp_dir>
#
# Outputs per run:  <outdir>/A_orig/md.gro  +  <outdir>/A_orig/md.xtc
#                   <outdir>/B_mut/md.gro   +  <outdir>/B_mut/md.xtc
set -euo pipefail

unset OMP_NUM_THREADS

ORIG="${1:?arg1: path to WT PDB}"
MUT="${2:?arg2: path to mutant PDB}"
OUT="${3:?arg3: output directory}"
MDP_DIR="${4:?arg4: directory containing *.mdp files}"

TEMP="${TEMP:-280}"
SALT="${SALT:-0.15}"
PADDING="${PADDING:-2.5}"
STEPS="${STEPS:-3000000}"
FF="${FF:-amber99sb-ildn}"
WATER="${WATER:-tip3p}"
MD_GPU_OPTS="${MD_GPU_OPTS:--nb gpu -pme gpu -pin on -ntmpi 1 -ntomp $(nproc)}"

mkdir -p "$OUT/A_orig" "$OUT/B_mut"

# Copy mdp files into output dir so GROMACS can find them
cp -f "$MDP_DIR"/*.mdp "$OUT"/

run_one() {
  local tag=$1
  local pdb=$2
  local dir="$OUT/$tag"
  echo "=== [$tag] GROMACS workflow in $dir ==="
  pushd "$dir" >/dev/null

  if [[ ! -f processed.gro ]]; then
    echo ">>> pdb2gmx"
    gmx pdb2gmx -f "$pdb" -o processed.gro -p topol.top -i posre.itp \
      -ff "$FF" -water "$WATER" -ignh <<EOF
1
EOF
  fi

  if [[ ! -f boxed.gro ]]; then
    echo ">>> editconf (box)"
    gmx editconf -f processed.gro -o boxed.gro -c -d "$PADDING" -bt cubic
  fi

  if [[ ! -f solv.gro ]]; then
    echo ">>> solvate"
    gmx solvate -cp boxed.gro -cs spc216.gro -o solv.gro -p topol.top
  fi

  if [[ ! -f solv_ions.gro ]]; then
    echo ">>> genion"
    gmx grompp -f ../minim.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 1
    echo "SOL" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top \
      -pname NA -nname CL -neutral -conc "$SALT"
  fi

  if [[ ! -f em.gro ]]; then
    echo ">>> energy minimization (CPU)"
    gmx grompp -f ../minim.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 1
    gmx mdrun -deffnm em -ntmpi 1 -ntomp $(nproc) -pin on
  fi

  if [[ ! -f nvt.gro ]]; then
    echo ">>> NVT equilibration (GPU)"
    if [[ -f nvt.cpt ]]; then
      echo ">>> Resuming NVT from nvt.cpt"
      gmx mdrun -deffnm nvt -cpi nvt.cpt -append $MD_GPU_OPTS
    else
      gmx grompp -f ../nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1
      gmx mdrun -deffnm nvt $MD_GPU_OPTS
    fi
  fi

  if [[ ! -f npt.gro ]]; then
    echo ">>> NPT equilibration (GPU)"
    if [[ -f npt.cpt ]]; then
      echo ">>> Resuming NPT from npt.cpt"
      gmx mdrun -deffnm npt -cpi npt.cpt -append $MD_GPU_OPTS
    else
      gmx grompp -f ../npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -maxwarn 1
      gmx mdrun -deffnm npt $MD_GPU_OPTS
    fi
  fi

  echo ">>> Production MD (GPU)"
  if [[ ! -f md.tpr ]]; then
    gmx grompp -f ../md.mdp -c npt.gro -p topol.top -o md.tpr -maxwarn 1
  fi

  if [[ -f md.cpt ]]; then
    echo ">>> Resuming from md.cpt"
    gmx mdrun -deffnm md -cpi md.cpt -append $MD_GPU_OPTS
  else
    gmx mdrun -deffnm md $MD_GPU_OPTS
  fi

  popd >/dev/null
  echo "=== [$tag] done ==="
}

run_one "A_orig" "$ORIG"
run_one "B_mut"  "$MUT"

echo "GROMACS simulation complete. Outputs in: $OUT"
