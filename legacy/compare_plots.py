#!/usr/bin/env python3
import sys, pandas as pd
import matplotlib.pyplot as plt

if len(sys.argv) < 4:
    print("Usage: compare_plots.py analysis_A.csv analysis_B.csv outdir")
    raise SystemExit(1)

dfa = pd.read_csv(sys.argv[1])
dfb = pd.read_csv(sys.argv[2])
outdir = sys.argv[3]

plt.figure()
plt.plot(dfa["time_ps"], dfa["Rg_nm"], label="Rg A (orig)")
plt.plot(dfb["time_ps"], dfb["Rg_nm"], label="Rg B (mut)")
plt.xlabel("Time (ps)"); plt.ylabel("Rg (nm)"); plt.legend(); plt.tight_layout()
plt.savefig(f"{outdir}/Rg_overlay.png", dpi=180)

plt.figure()
plt.plot(dfa["time_ps"], dfa["end_to_end_nm"], label="E2E A (orig)")
plt.plot(dfb["time_ps"], dfb["end_to_end_nm"], label="E2E B (mut)")
plt.xlabel("Time (ps)"); plt.ylabel("E2E distance (nm)"); plt.legend(); plt.tight_layout()
plt.savefig(f"{outdir}/EndToEnd_overlay.png", dpi=180)

summary = pd.DataFrame({
    "metric": ["Rg_nm", "end_to_end_nm"],
    "A_mean": [dfa["Rg_nm"].mean(), dfa["end_to_end_nm"].mean()],
    "B_mean": [dfb["Rg_nm"].mean(), dfb["end_to_end_nm"].mean()],
    "A_std":  [dfa["Rg_nm"].std(),  dfa["end_to_end_nm"].std()],
    "B_std":  [dfb["Rg_nm"].std(),  dfb["end_to_end_nm"].std()],
})
summary.to_csv(f"{outdir}/summary_AB.csv", index=False)
print(f"Wrote {outdir}/Rg_overlay.png, EndToEnd_overlay.png, summary_AB.csv")
