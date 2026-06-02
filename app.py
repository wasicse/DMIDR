"""
DMIDR — Drosophila Matrisome Intrinsically Disordered Regions Pipeline Dashboard
Run:  streamlit run app.py
"""

import os
import re
import subprocess
import time
from datetime import datetime
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

PROJECT_ROOT = Path(__file__).parent
INPUT_DIR    = PROJECT_ROOT / "data" / "input"
OUTPUT_DIR   = PROJECT_ROOT / "outputs"

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBXZUOJ")

# ── Page config ───────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="DMIDR — Drosophila Matrisome IDR Pipeline",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
code { font-size: 0.78rem; }
.stProgress > div > div { border-radius: 4px; }
</style>
""", unsafe_allow_html=True)

# ── Session state ─────────────────────────────────────────────────────────────
for k, v in {"pid": None, "log_path": None, "fasta_saved": None}.items():
    if k not in st.session_state:
        st.session_state[k] = v

# ── Helpers ───────────────────────────────────────────────────────────────────
def sh(cmd: str) -> str:
    try:
        return subprocess.check_output(cmd, shell=True, stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:
        return ""

def is_running() -> bool:
    pid = st.session_state.pid
    return bool(pid and Path(f"/proc/{pid}").exists())

def load_env(path: Path) -> dict:
    cfg: dict = {}
    if path.exists():
        for line in path.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                k, _, v = line.partition("=")
                cfg[k.strip()] = v.strip().strip('"').strip("'")
    return cfg

def save_env(path: Path, cfg: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f'{k}="{v}"' for k, v in cfg.items() if v]
    path.write_text("\n".join(lines) + "\n")

# ── FASTA helpers ─────────────────────────────────────────────────────────────
def parse_fasta(text: str) -> list[dict]:
    """Return list of {id, seq} from FASTA text."""
    records = []
    current_id, current_seq = None, []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_id:
                records.append({"id": current_id, "seq": "".join(current_seq).upper()})
            current_id = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line)
    if current_id:
        records.append({"id": current_id, "seq": "".join(current_seq).upper()})
    return records

def validate_fasta(records: list[dict]) -> list[str]:
    errors = []
    if not records:
        errors.append("No sequences found.")
        return errors
    for r in records:
        seq = r["seq"]
        if not seq:
            errors.append(f"{r['id']}: empty sequence.")
            continue
        bad = set(seq) - VALID_AA
        if bad:
            errors.append(f"{r['id']}: invalid characters {bad}")
    return errors

def seq_stats(seq: str) -> dict:
    length = len(seq)
    charged = sum(seq.count(aa) for aa in "DEKR")
    hydro   = sum(seq.count(aa) for aa in "AVILMFW")
    return {"Length": length,
            "Charged (%)": f"{charged/length*100:.1f}" if length else "—",
            "Hydrophobic (%)": f"{hydro/length*100:.1f}" if length else "—"}

# ── Pipeline helpers ──────────────────────────────────────────────────────────
def launch(steps: str, config_path: Path) -> int:
    log = OUTPUT_DIR / f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, "PIPELINE_STEPS": steps}
    cmd = f'bash "{PROJECT_ROOT}/run_all.sh" "{config_path}"'
    with open(log, "w") as fh:
        proc = subprocess.Popen(cmd, shell=True, stdout=fh, stderr=fh,
                                env=env, cwd=str(PROJECT_ROOT))
    st.session_state.pid      = proc.pid
    st.session_state.log_path = str(log)
    return proc.pid

def tail_log(path: str, n: int = 100) -> str:
    p = Path(path)
    if not p.exists():
        return ""
    lines = p.read_text(errors="replace").splitlines()
    # Filter out low-value noise
    keep = [l for l in lines if not any(x in l for x in
            ("oneDNN", "absl::InitializeLog", "I0000", "port.cc",
             "TF_ENABLE_ONEDNN", "WARNING: All log messages"))]
    return "\n".join(keep[-n:])

@st.cache_data(ttl=5)
def outputs_summary() -> dict:
    result: dict = {}
    if not OUTPUT_DIR.exists():
        return result
    for seq_dir in sorted(OUTPUT_DIR.iterdir()):
        if not seq_dir.is_dir() or seq_dir.name in ("pssm_cache", "sequences"):
            continue
        name = seq_dir.name
        disp  = seq_dir / "dispred"
        caids = sorted(disp.glob("*.caid")) if disp.exists() else []
        result[name] = {
            "pssm":      (seq_dir / "pssm" / f"{name}.pssm").exists(),
            "wt_caid":   (disp / f"{name}_original.caid").exists(),
            "mut_caids": [c for c in caids if "mutants" in c.name],
            "all_caids": caids,
            "analysis":  list((seq_dir).glob("disorder_*res/*.csv")),
        }
    return result

# ── GPU / Docker helpers ──────────────────────────────────────────────────────
@st.cache_data(ttl=3)
def gpu_info() -> dict | None:
    raw = sh("nvidia-smi --query-gpu=name,memory.used,memory.total,"
             "utilization.gpu,temperature.gpu --format=csv,noheader,nounits")
    if not raw:
        return None
    p = [x.strip() for x in raw.split(",")]
    if len(p) >= 5:
        return dict(name=p[0], mem_used=int(p[1]), mem_total=int(p[2]),
                    util=int(p[3]), temp=int(p[4]))
    return None

@st.cache_data(ttl=3)
def gpu_procs() -> list[dict]:
    raw = sh("nvidia-smi --query-compute-apps=pid,used_memory,process_name "
             "--format=csv,noheader")
    rows = []
    for line in raw.splitlines():
        p = [x.strip() for x in line.split(",")]
        if len(p) == 3:
            rows.append({"PID": p[0], "VRAM": p[1], "Process": Path(p[2]).name})
    return rows

@st.cache_data(ttl=8)
def docker_containers() -> list[dict]:
    raw = sh('docker ps --format "{{.ID}}|{{.Image}}|{{.Status}}|{{.Names}}"')
    rows = []
    for line in raw.splitlines():
        p = line.split("|")
        if len(p) == 4:
            rows.append({"id": p[0], "image": p[1], "status": p[2], "name": p[3]})
    return rows

@st.cache_data(ttl=8)
def container_progress(cid: str) -> dict | None:
    raw = sh(
        f'docker exec {cid} sh -c \''
        f'feat=$(ls /opt/ESMDisPred/features/Dispredict3.0/features/ 2>/dev/null | wc -l); '
        f'total=$(grep -c ">" /opt/ESMDisPred/example/*.fasta 2>/dev/null || echo 0); '
        f'fasta=$(ls /opt/ESMDisPred/example/*.fasta 2>/dev/null | head -1 | xargs -r basename); '
        f'esm2=$(ls /opt/ESMDisPred/features/ESM2/Mutant_*.csv 2>/dev/null | wc -l); '
        f'echo "$feat|$total|$fasta|$esm2"\''
    )
    p = raw.split("|")
    if len(p) == 4:
        try:
            return {"feat": int(p[0].strip()), "total": int(p[1].strip()),
                    "fasta": p[2].strip(), "esm2": int(p[3].strip())}
        except Exception:
            pass
    return None

# ══════════════════════════════════════════════════════════════════════════════
# Sidebar
# ══════════════════════════════════════════════════════════════════════════════
cfg_path = PROJECT_ROOT / "configs" / "local.env"
cfg = load_env(cfg_path)

with st.sidebar:
    st.title("🧬 DMIDR")
    st.caption("Drosophila Matrisome Intrinsically Disordered Regions Pipeline")

    # ── FASTA Input ───────────────────────────────────────────────────────────
    st.subheader("Sequence Input")

    input_tab, upload_tab = st.tabs(["Paste sequence", "Upload FASTA"])

    fasta_text = ""

    with input_tab:
        seq_name_input = st.text_input("Sequence name", placeholder="e.g. MyProtein")
        raw_seq = st.text_area("Paste protein sequence (single letter code)",
                               height=120,
                               placeholder="MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLRWVQTLSEQVQEELLSSQVTQELRALMDETMKELKAYKSELEEQLTPVAEETRARLSKELQAAQARLGADMEDVCGRLVQYRGEVQAMLGQSTEELRVRLASHLRKLRKRLLRDADDLQKRLAVYQAGAREGAERGLSAIRERLGPLVEQGRVRAATVGSLAGQPLQERAQAWGERLRARMEEMGSRTRDRLDEVKEQVAEVRAKLEEQAQQIRLQAEAFQARLKSWFEPLVEDMQRQWAGLVEKVQAAVGTSAAPVPSDNH")
        if seq_name_input and raw_seq.strip():
            seq_clean = re.sub(r'\s+', '', raw_seq).upper()
            bad = set(seq_clean) - VALID_AA
            if bad:
                st.error(f"Invalid characters: {bad}")
            else:
                fasta_text = f">{seq_name_input}\n{seq_clean}\n"
                st.caption(f"Length: {len(seq_clean)} aa")

    with upload_tab:
        uploaded = st.file_uploader("Upload .fasta / .fa file",
                                    type=["fasta", "fa", "txt"])
        if uploaded:
            fasta_text = uploaded.read().decode("utf-8", errors="replace")

    # Preview + save
    if fasta_text:
        records = parse_fasta(fasta_text)
        errors  = validate_fasta(records)
        if errors:
            for e in errors:
                st.error(e)
        else:
            st.success(f"{len(records)} sequence(s) ready")
            for r in records[:3]:
                st.caption(f"▸ **{r['id']}** — {len(r['seq'])} aa")
            if len(records) > 3:
                st.caption(f"  … and {len(records)-3} more")

            save_name = (records[0]["id"].replace(" ", "_").replace("/", "_")
                         + ".fasta")
            INPUT_DIR.mkdir(parents=True, exist_ok=True)
            save_path = INPUT_DIR / save_name

            if st.button("💾 Use this sequence", use_container_width=True):
                save_path.write_text(fasta_text)
                cfg["INPUT_FASTA"] = str(save_path)
                save_env(cfg_path, cfg)
                st.session_state.fasta_saved = str(save_path)
                st.success(f"Saved → {save_name}")
                st.rerun()

    st.divider()

    # ── Config ────────────────────────────────────────────────────────────────
    st.subheader("Config")
    cfg["BLAST_DB"]          = st.text_input("BLAST_DB",       value=cfg.get("BLAST_DB", ""))
    cfg["ESMDISPRED_IMAGE"]  = st.text_input("Docker image",   value=cfg.get("ESMDISPRED_IMAGE", "wasicse/esmdispred:latest"))
    cfg["LARGE_MODELS_DIR"]  = st.text_input("Large models",   value=cfg.get("LARGE_MODELS_DIR", ""))

    # FASTA selector from saved files
    fastas = sorted(INPUT_DIR.glob("*.fasta")) if INPUT_DIR.exists() else []
    if fastas:
        current = cfg.get("INPUT_FASTA", "")
        idx = next((i for i, f in enumerate(fastas) if str(f) == current), 0)
        sel = st.selectbox("Active FASTA", fastas,
                           format_func=lambda p: p.name, index=idx)
        cfg["INPUT_FASTA"] = str(sel)

    cfg["MAX_MUTATIONS"] = str(st.slider("Max mutations", 1, 10,
                                          int(cfg.get("MAX_MUTATIONS", 5))))
    cfg["MODEL"] = str(st.selectbox("Model", [1, 2, 3], index=2))

    if st.button("💾 Save config", use_container_width=True):
        save_env(cfg_path, cfg)
        st.success("Saved.")

    st.divider()

    # ── Steps ─────────────────────────────────────────────────────────────────
    st.subheader("Steps")
    step_labels = [
        (1, "PSSM"), (2, "ESMDisPred WT"), (3, "Generate mutants"),
        (4, "ESMDisPred mutants"), (5, "Disorder analysis"),
        (6, "ColabFold"), (7, "ASA"), (8, "MD analysis"),
    ]
    selected = [str(i) for i, lbl in step_labels
                if st.checkbox(lbl, value=(i <= 5), key=f"s{i}")]

    st.divider()
    col_run, col_stop = st.columns(2)
    run_clicked  = col_run.button("▶ Run",  type="primary",
                                  disabled=(is_running() or not selected),
                                  use_container_width=True)
    stop_clicked = col_stop.button("⏹ Stop", disabled=not is_running(),
                                   use_container_width=True)

# ── Actions ───────────────────────────────────────────────────────────────────
if run_clicked and selected:
    save_env(cfg_path, cfg)
    launch(" ".join(selected), cfg_path)
    st.rerun()

if stop_clicked and st.session_state.pid:
    sh(f"kill -- -{st.session_state.pid}")
    st.session_state.pid = None
    st.rerun()

# ══════════════════════════════════════════════════════════════════════════════
# Main tabs
# ══════════════════════════════════════════════════════════════════════════════
tab_run, tab_gpu, tab_results = st.tabs(
    ["▶ Run & Monitor", "🖥 GPU & Containers", "📊 Results"])

# Compute once — shared across all three tabs
_summary = outputs_summary()

# ─────────────────────────────────────────────────────────────────────────────
# Tab 1 — Run & Monitor
# ─────────────────────────────────────────────────────────────────────────────
with tab_run:
    running = is_running()

    # Status banner
    if running:
        st.success(f"**Pipeline running** — PID {st.session_state.pid}")
    elif st.session_state.pid:
        st.info("Pipeline finished.")
    else:
        st.info("Configure a sequence and steps in the sidebar, then click **▶ Run**.")

    # Per-sequence output summary
    summary = _summary
    if summary:
        st.subheader("Outputs")
        for seq, info in summary.items():
            mut_done = len(info["mut_caids"])
            max_mut  = int(cfg.get("MAX_MUTATIONS", 5))
            with st.expander(f"**{seq}**  —  {mut_done}/{max_mut} mutation blocks",
                             expanded=True):
                c = st.columns(5)
                c[0].metric("PSSM",       "✅" if info["pssm"]    else "⬜")
                c[1].metric("WT pred",    "✅" if info["wt_caid"] else "⬜")
                c[2].metric("Mut preds",  mut_done)
                c[3].metric("Analysis",   len(info["analysis"]))
                sizes = [m.group(1) for caid in info["mut_caids"]
                         for m in [re.search(r'(\d+)res', caid.name)] if m]
                c[4].metric("Mut sizes", ", ".join(sizes) or "—")

    # Live log
    st.subheader("Log")
    log_path = st.session_state.log_path
    if not log_path:
        logs = sorted(OUTPUT_DIR.glob("run_*.log"), reverse=True) if OUTPUT_DIR.exists() else []
        log_path = str(logs[0]) if logs else None

    if log_path:
        st.caption(f"📄 {Path(log_path).name}")
        log_content = tail_log(log_path)
        st.code(log_content or "(empty)", language=None)
    else:
        st.caption("No log yet.")

    if running:
        st.html('<meta http-equiv="refresh" content="5">')

# ─────────────────────────────────────────────────────────────────────────────
# Tab 2 — GPU & Containers
# ─────────────────────────────────────────────────────────────────────────────
with tab_gpu:
    col_gpu, col_docker = st.columns(2)

    with col_gpu:
        st.subheader("GPU")
        info = gpu_info()
        if info:
            st.markdown(f"**{info['name']}**")
            mem_pct = info["mem_used"] / info["mem_total"]
            st.progress(mem_pct,
                        text=f"VRAM  {info['mem_used']:,} / {info['mem_total']:,} MiB  ({mem_pct*100:.0f}%)")
            st.progress(info["util"] / 100,
                        text=f"Utilization  {info['util']}%")
            col_t, col_f = st.columns(2)
            col_t.metric("Temperature", f"{info['temp']} °C")
            col_f.metric("Free VRAM", f"{info['mem_total']-info['mem_used']:,} MiB")

            procs = gpu_procs()
            if procs:
                st.caption("Active GPU processes")
                st.dataframe(pd.DataFrame(procs), hide_index=True,
                             use_container_width=True)
        else:
            st.warning("No GPU detected.")

    with col_docker:
        st.subheader("Containers")
        containers = docker_containers()
        esm   = [c for c in containers if "esmdispred" in c["image"].lower()]
        other = [c for c in containers if "esmdispred" not in c["image"].lower()]

        if esm:
            st.caption("**ESMDisPred**")
            for c in esm:
                with st.container(border=True):
                    st.markdown(f"🟢 **{c['name']}**  `{c['id'][:12]}`")
                    st.caption(c["status"])
                    prog = container_progress(c["id"])
                    if prog:
                        if prog["total"] > 0:
                            pct = prog["feat"] / prog["total"]
                            st.progress(pct,
                                text=f"DisPredict3.0  {prog['feat']}/{prog['total']}  ({pct*100:.0f}%)")
                        if prog["esm2"] > 0:
                            st.caption(f"ESM2 embeddings: {prog['esm2']}")
                        if prog["fasta"]:
                            st.caption(f"📂 {prog['fasta']}")

        if other:
            st.caption("**Other**")
            rows = [{"Name": c["name"], "Image": c["image"], "Status": c["status"]}
                    for c in other]
            st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)

        if not containers:
            st.info("No running containers.")

    if st.button("↻ Refresh", key="refresh_gpu"):
        st.rerun()

    if esm:
        st.html('<meta http-equiv="refresh" content="15">')

# ─────────────────────────────────────────────────────────────────────────────
# Tab 3 — Results
# ─────────────────────────────────────────────────────────────────────────────
def _wt_profile_fig(wt_df: pd.DataFrame, seq: str) -> go.Figure:
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=wt_df["pos"], y=wt_df["prob"],
        mode="lines", name="Wild-type",
        line=dict(color="#2196F3", width=2),
    ))
    fig.add_hline(y=0.5, line_dash="dash", line_color="gray",
                  annotation_text="Disorder threshold (0.5)", annotation_position="bottom right")
    fig.update_layout(
        title=f"{seq} — Wild-Type Disorder Profile",
        xaxis_title="Residue Position",
        yaxis_title="Disorder Probability",
        yaxis_range=[0, 1],
        height=340,
        margin=dict(l=40, r=20, t=50, b=40),
        legend=dict(x=0.01, y=0.99),
        hovermode="x unified",
    )
    return fig


def _comparison_figs(comp_df: pd.DataFrame, title: str) -> tuple[go.Figure, go.Figure]:
    # Overlay chart
    fig_ov = go.Figure()
    fig_ov.add_trace(go.Scatter(
        x=comp_df["Position"], y=comp_df["DisorderProb_orig"],
        mode="lines", name="Wild-type",
        line=dict(color="#2196F3", width=2),
    ))
    fig_ov.add_trace(go.Scatter(
        x=comp_df["Position"], y=comp_df["DisorderProb_mut"],
        mode="lines", name="Mutant",
        line=dict(color="#F44336", width=2),
    ))
    # Shaded region between WT and mutant
    fig_ov.add_trace(go.Scatter(
        x=list(comp_df["Position"]) + list(comp_df["Position"])[::-1],
        y=list(comp_df["DisorderProb_orig"]) + list(comp_df["DisorderProb_mut"])[::-1],
        fill="toself",
        fillcolor="rgba(76,175,80,0.12)",
        line=dict(color="rgba(0,0,0,0)"),
        showlegend=False,
        hoverinfo="skip",
    ))
    fig_ov.add_hline(y=0.5, line_dash="dash", line_color="gray",
                     annotation_text="Threshold (0.5)", annotation_position="bottom right")
    fig_ov.update_layout(
        title=title,
        xaxis_title="Residue Position",
        yaxis_title="Disorder Probability",
        yaxis_range=[0, 1],
        height=360,
        margin=dict(l=40, r=20, t=50, b=40),
        legend=dict(x=0.01, y=0.99),
        hovermode="x unified",
    )

    # Delta bar chart
    colors = ["#4CAF50" if v > 0 else "#F44336" for v in comp_df["DeltaDisorder"]]
    fig_delta = go.Figure()
    fig_delta.add_trace(go.Bar(
        x=comp_df["Position"],
        y=comp_df["DeltaDisorder"],
        marker_color=colors,
        name="ΔDisorder",
    ))
    fig_delta.update_layout(
        title="ΔDisorder per Residue  (green = reduction in mutant)",
        xaxis_title="Residue Position",
        yaxis_title="ΔDisorder (WT − Mut)",
        height=260,
        margin=dict(l=40, r=20, t=50, b=40),
        showlegend=False,
    )
    return fig_ov, fig_delta


with tab_results:
    summary = _summary
    if not summary:
        st.info("No outputs yet. Run the pipeline first.")
    else:
        seq_names = list(summary.keys())
        selected_seq = st.selectbox("Select sequence", seq_names, key="results_seq_sel")

        if selected_seq:
            seq  = selected_seq
            info = summary[seq]
            seq_dir   = OUTPUT_DIR / seq
            disp_dir  = seq_dir / "dispred"
            wt_caid   = disp_dir / f"{seq}_original.caid"
            cand_csv  = seq_dir / "candidates" / "candidate_scores.csv"

            # ── Summary metrics ───────────────────────────────────────────────
            wt_df = None
            if wt_caid.exists():
                try:
                    wt_df = pd.read_csv(wt_caid, sep="\t", header=None,
                                        names=["pos", "aa", "prob", "label"])
                except Exception:
                    pass

            m1, m2, m3, m4 = st.columns(4)
            if wt_df is not None:
                n_res = len(wt_df)
                n_dis = int(wt_df["label"].sum())
                m1.metric("Residues", f"{n_res:,}")
                m2.metric("WT Disordered", f"{n_dis} ({n_dis / n_res * 100:.1f}%)")
            else:
                m1.metric("Residues", "—")
                m2.metric("WT Disordered", "—")

            m3.metric("Mutation blocks", len(info["mut_caids"]))

            cand_df = None
            if cand_csv.exists():
                try:
                    cand_df = pd.read_csv(cand_csv)
                except Exception:
                    pass
            m4.metric("Candidates scored", len(cand_df) if cand_df is not None else 0)

            st.divider()

            # ── Candidate rankings ────────────────────────────────────────────
            if cand_df is not None:
                st.subheader("Candidate Rankings")
                try:
                    top_n = st.slider(
                        "Show top N", 5, min(50, len(cand_df)),
                        min(10, len(cand_df)), key="cand_top_n"
                    )
                    display_df = cand_df.head(top_n).copy()
                    display_df.insert(0, "Rank", range(1, len(display_df) + 1))
                    st.dataframe(
                        display_df,
                        use_container_width=True,
                        hide_index=True,
                        column_config={
                            "Rank": st.column_config.NumberColumn(width="small"),
                            "mutant_id": st.column_config.TextColumn("Mutant ID"),
                            "block_size": st.column_config.NumberColumn("Block Size", width="small"),
                            "avg_disorder_reduction": st.column_config.NumberColumn(
                                "Avg Reduction ↑", format="%.4f",
                                help="Mean (DisorderProb_WT − DisorderProb_mut) — higher is better"
                            ),
                            "total_disorder_reduction": st.column_config.NumberColumn(
                                "Total Reduction", format="%.2f"
                            ),
                            "disordered_to_ordered": st.column_config.NumberColumn("Dis→Ord ↑", width="small",
                                help="Residues that flipped from disordered (label=1) to ordered (label=0)"),
                            "ordered_to_disordered": st.column_config.NumberColumn("Ord→Dis ↓", width="small"),
                        },
                    )

                    dl_col, _ = st.columns([1, 3])
                    dl_col.download_button(
                        "⬇ Full rankings CSV",
                        data=cand_csv.read_bytes(),
                        file_name=f"{seq}_candidate_scores.csv",
                        mime="text/csv",
                        key=f"dl_scores_{seq}",
                    )

                    # Per-candidate FASTA downloads
                    cand_dir   = seq_dir / "candidates"
                    cand_fastas = sorted(cand_dir.glob("candidate_*.fasta"))
                    if cand_fastas:
                        st.caption("**Download candidate sequences:**")
                        btn_cols = st.columns(min(len(cand_fastas), 4))
                        for i, fp in enumerate(cand_fastas):
                            btn_cols[i % 4].download_button(
                                f"⬇ {fp.stem[:30]}",
                                data=fp.read_bytes(),
                                file_name=fp.name,
                                mime="text/plain",
                                key=f"dl_cand_{seq}_{i}",
                            )
                except Exception as e:
                    st.error(f"Could not load candidate scores: {e}")

                st.divider()

            # ── Disorder profiles ─────────────────────────────────────────────
            st.subheader("Disorder Profiles")

            if wt_df is not None:
                st.plotly_chart(_wt_profile_fig(wt_df, seq), use_container_width=True)

                disorder_dirs = sorted(seq_dir.glob("disorder_*res"))
                if disorder_dirs:
                    block_names  = [d.name for d in disorder_dirs]
                    sel_block    = st.selectbox("Block size comparison", block_names, key=f"blk_{seq}")
                    block_dir    = seq_dir / sel_block
                    comp_csv_path = block_dir / f"{seq}_{sel_block.replace('disorder_','')}_disorder_probability_comparison.csv"

                    if comp_csv_path.exists():
                        try:
                            comp_df = pd.read_csv(comp_csv_path)
                            # Summary bar under charts
                            avg_red = comp_df["DeltaDisorder"].mean()
                            n_improved = int((comp_df["DeltaDisorder"] > 0).sum())
                            ca, cb, cc = st.columns(3)
                            ca.metric("Avg ΔDisorder (WT−Mut)", f"{avg_red:.4f}",
                                      delta_color="normal" if avg_red > 0 else "inverse")
                            cb.metric("Residues improved", n_improved)
                            cc.metric("Net improvement",
                                      f"{n_improved / len(comp_df) * 100:.1f}%")

                            fig_ov, fig_delta = _comparison_figs(
                                comp_df, f"{seq} — WT vs Best {sel_block} Mutant"
                            )
                            st.plotly_chart(fig_ov,    use_container_width=True)
                            st.plotly_chart(fig_delta, use_container_width=True)
                        except Exception as e:
                            st.error(f"Could not load comparison CSV: {e}")
                    else:
                        st.caption("No disorder comparison CSV for this block size yet.")

                    # PNG plots
                    png_plots = sorted(block_dir.glob("*.png"))
                    if png_plots:
                        with st.expander("PNG disorder plots"):
                            for png in png_plots:
                                st.image(str(png), caption=png.name, use_container_width=True)
            else:
                st.caption("Run ESMDisPred (step 2) to see disorder profiles.")

            st.divider()

            # ── Structure files ───────────────────────────────────────────────
            structure_dir = seq_dir / "structure"
            pdb_files = sorted(structure_dir.glob("**/*.pdb")) if structure_dir.exists() else []
            if pdb_files:
                st.subheader("Structure Files")
                pdb_cols = st.columns(min(len(pdb_files), 3))
                for i, pdb in enumerate(pdb_files):
                    pdb_cols[i % 3].download_button(
                        f"⬇ {pdb.name}",
                        data=pdb.read_bytes(),
                        file_name=pdb.name,
                        mime="chemical/x-pdb",
                        key=f"dl_pdb_{seq}_{i}",
                    )

            # ── Raw .caid viewer ──────────────────────────────────────────────
            all_caids = info["all_caids"]
            if all_caids:
                with st.expander("Raw .caid file viewer"):
                    sel_caid = st.selectbox(
                        "Prediction file", [c.name for c in all_caids], key=f"caid_{seq}"
                    )
                    caid_path = next(c for c in all_caids if c.name == sel_caid)
                    try:
                        caid_df = pd.read_csv(caid_path, sep="\t", header=None,
                                              names=["pos", "aa", "score", "binary"])
                        col_chart, col_stats = st.columns([3, 1])
                        with col_chart:
                            st.line_chart(caid_df.set_index("pos")["score"],
                                          use_container_width=True, color="#e05c2a")
                        with col_stats:
                            n_dis = int(caid_df["binary"].sum())
                            st.metric("Residues",  len(caid_df))
                            st.metric("Disordered", n_dis)
                            st.metric("Fraction",  f"{n_dis / len(caid_df) * 100:.1f}%")
                        st.download_button(
                            "⬇ Download .caid",
                            data=caid_path.read_bytes(),
                            file_name=caid_path.name,
                            mime="text/plain",
                            key=f"dl_{seq}_{sel_caid}",
                        )
                        with st.expander("Data table"):
                            st.dataframe(caid_df, use_container_width=True, hide_index=True)
                    except Exception as e:
                        st.error(f"Could not parse {sel_caid}: {e}")
