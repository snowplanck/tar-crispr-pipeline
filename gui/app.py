"""Streamlit GUI for the TAR-CRISPR pipeline.

Thin UI layer on top of ``gui/run_pipeline.py``. No pipeline logic here:
widgets -> a dict of arguments -> run_pipeline() -> render results.
"""
from __future__ import annotations

import io
import tempfile
import zipfile
from pathlib import Path

import streamlit as st

# Ensure the package is importable when running `streamlit run gui/app.py`
# from the repository root.
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from gui.run_pipeline import run_pipeline, _enzyme_site  # noqa: E402


# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="TAR-CRISPR",
    page_icon="🧬",
    layout="wide",
)

st.title("🧬 TAR-CRISPR Pipeline")
st.caption(
    "Automated sgRNA + homology arm + tailed primer design for "
    "CATCH-style BGC capture in *Streptomyces*."
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _save_upload(uploaded, subdir: Path) -> Path:
    """Save a Streamlit UploadedFile to disk and return the path."""
    subdir.mkdir(parents=True, exist_ok=True)
    path = subdir / uploaded.name
    path.write_bytes(uploaded.getbuffer())
    return path


def _zip_dir(directory: Path) -> bytes:
    """Zip all files in *directory* and return the bytes."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in sorted(directory.iterdir()):
            if f.is_file():
                zf.write(f, arcname=f.name)
    buf.seek(0)
    return buf.getvalue()


# ---------------------------------------------------------------------------
# Section 1 — Inputs
# ---------------------------------------------------------------------------
st.header("1. Inputs")

col1, col2, col3 = st.columns(3)
with col1:
    bgc_file = st.file_uploader(
        "BGC (FASTA or GenBank)", type=["fasta", "fa", "fna", "gbk", "gb", "genbank"]
    )
with col2:
    genome_file = st.file_uploader(
        "Genome (FASTA)", type=["fasta", "fa", "fna"]
    )
with col3:
    vector_file = st.file_uploader(
        "Capture vector (FASTA or GenBank)",
        type=["fasta", "fa", "fna", "gbk", "gb", "genbank"],
    )


# ---------------------------------------------------------------------------
# Section 2 — Parameters
# ---------------------------------------------------------------------------
st.header("2. Parameters")

p1, p2, p3, p4 = st.columns(4)
with p1:
    start = st.number_input(
        "Cluster start (0-based)", min_value=0, value=0, step=1,
        help="Leave 0 and set end=0 to rely on GenBank features.",
    )
with p2:
    end = st.number_input(
        "Cluster end (0-based)", min_value=0, value=0, step=1,
    )
with p3:
    enzyme = st.selectbox(
        "Vector linearization enzyme",
        options=[
            "(none)", "SwaI", "PmeI", "NruI", "AatII", "PvuI",
            "KpnI", "NheI", "ClaI", "EcoRI", "BamHI", "HindIII",
        ],
    )
with p4:
    use_blast = st.checkbox(
        "Use BLAST for specificity", value=False,
        help="Slower. Leave off to use the in-process Hamming fallback.",
    )

with st.expander("Advanced parameters"):
    a1, a2, a3 = st.columns(3)
    with a1:
        pam_window = st.number_input("PAM window (bp)", min_value=50, value=500, step=50)
        arm_length = st.number_input("Homology arm length (bp)", min_value=20, value=50, step=5)
        top_n = st.number_input("Top N sgRNAs per end", min_value=1, value=5, step=1)
    with a2:
        tm_min = st.number_input("Primer Tm min (°C)", min_value=40.0, value=58.0, step=0.5)
        tm_max = st.number_input("Primer Tm max (°C)", min_value=40.0, value=62.0, step=0.5)
    with a3:
        max_mismatches = st.number_input("Max mismatches", min_value=0, value=3, step=1)


# ---------------------------------------------------------------------------
# Section 3 — Run
# ---------------------------------------------------------------------------
st.header("3. Run")

run_button = st.button(
    "▶ Run pipeline", type="primary",
    disabled=not (bgc_file and genome_file and vector_file),
)

if not (bgc_file and genome_file and vector_file):
    st.info("Upload BGC, genome, and vector to enable the run button.")


if run_button:
    # Persist uploads in a stable folder for this run
    workdir = Path(tempfile.mkdtemp(prefix="tar_crispr_gui_", dir="."))
    input_dir = workdir / "inputs"
    output_dir = workdir / "outputs"

    bgc_path = _save_upload(bgc_file, input_dir)
    genome_path = _save_upload(genome_file, input_dir)
    vector_path = _save_upload(vector_file, input_dir)

    progress_bar = st.progress(0.0, text="Starting...")
    status = st.empty()

    def _on_progress(msg: str, frac: float) -> None:
        progress_bar.progress(min(max(frac, 0.0), 1.0), text=msg)
        status.write(msg)

    kwargs = dict(
        bgc_path=bgc_path,
        genome_path=genome_path,
        vector_path=vector_path,
        output_dir=output_dir,
        vector_enzyme=None if enzyme == "(none)" else enzyme,
        use_blast=use_blast,
        pam_window=int(pam_window),
        arm_length=int(arm_length),
        top_n=int(top_n),
        tm_min=float(tm_min),
        tm_max=float(tm_max),
        max_mismatches=int(max_mismatches),
        progress=_on_progress,
    )
    if start > 0 and end > start:
        kwargs["start"] = int(start)
        kwargs["end"] = int(end)

    try:
        result = run_pipeline(**kwargs)
    except Exception as e:
        st.error(f"Pipeline failed: {e}")
        st.exception(e)
        st.stop()

    progress_bar.progress(1.0, text="Done.")
    st.success(f"Pipeline completed. Output: `{output_dir}`")
    st.session_state["result"] = result


# ---------------------------------------------------------------------------
# Section 4 — Results
# ---------------------------------------------------------------------------
result = st.session_state.get("result")

if result is not None:
    st.header("4. Results")

    # Download buttons
    dcols = st.columns(4)
    with dcols[0]:
        st.download_button(
            "⬇ report.md",
            data=result["report_md"].read_bytes(),
            file_name="report.md",
            mime="text/markdown",
        )
    with dcols[1]:
        st.download_button(
            "⬇ primers.csv",
            data=result["primers_csv"].read_bytes(),
            file_name="primers.csv",
            mime="text/csv",
        )
    with dcols[2]:
        st.download_button(
            "⬇ cluster_map.svg",
            data=result["cluster_map_svg"].read_bytes(),
            file_name="cluster_map.svg",
            mime="image/svg+xml",
        )
    with dcols[3]:
        st.download_button(
            "⬇ all (.zip)",
            data=_zip_dir(result["output_dir"]),
            file_name="tar_crispr_outputs.zip",
            mime="application/zip",
        )

    st.divider()

    # Cluster summary
    c = result["cluster"]
    st.subheader("Cluster")
    st.write(
        f"**{c.name}** — start `{c.start}`, end `{c.end}`, "
        f"length `{c.length}` bp, strand `{c.strand}`"
    )

    # Selected sgRNAs
    st.subheader("Selected sgRNAs")
    sel = result["selected"]
    for side in ("left", "right"):
        s = sel[side]
        st.markdown(
            f"- **{side}**: `{s.protospacer}`  "
            f"(PAM `{s.pam}`, strand `{s.strand}`, cut `{s.cut_position}`, "
            f"GC {s.gc_percent:.1f}%)"
        )

    # Primers
    st.subheader("Tailed primers")
    primers = result["primers"]
    rows = []
    for key in ("left", "right"):
        p = primers[key]
        rows.append({
            "name": p.name,
            "sequence": p.sequence,
            "tm": p.tm,
            "gc%": p.gc_percent,
            "hairpin": p.hairpin,
            "self_dimer": p.self_dimer,
            "cross_dimer": p.cross_dimer,
            "valid": p.valid,
            "issues": "; ".join(p.issues) if p.issues else "",
            "warnings": "; ".join(getattr(p, "warnings", []) or []),
        })
    st.dataframe(rows, use_container_width=True)

    # Report markdown
    st.subheader("Full report")
    md_text = result["report_md"].read_text(encoding="utf-8")
    with st.expander("Show report.md", expanded=False):
        st.markdown(md_text)

    # SVG inline
    svg_text = result["cluster_map_svg"].read_text(encoding="utf-8")
    st.subheader("Cluster map")
    st.components.v1.html(svg_text, height=250, scrolling=True)
