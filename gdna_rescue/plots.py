"""Per-transcript IGV-like coverage plots for surviving consensus transcripts.

Investigating a suspected novel transcript normally means opening IGV, loading
every sample, and eyeballing coverage over the locus. This module rebuilds that
view as a self-contained HTML page per transcript, so the whole cohort's coverage
plus the processed evidence sit on one page.

For each surviving locus we draw, per sample, an IGV-like stacked coverage track
over the locus +/- a shoulder: ``+`` strand above the axis, ``-`` below, and
within each strand the coverage is split into unique / duplicate / multimapper so
a reviewer can see how much of any pileup is genuine signal versus PCR-duplicate
or low-quality inflation. A feature track underneath shows the external annotation
(enhancer / promoter / CTCF / motif) and reference genes at their real
coordinates, and an evidence panel carries every metric already computed.

Reads coverage back from the per-sample ``*.coverage.h5`` stores (h5py) and the
overlay indexes merge_candidates already built -- NO pysam, so this runs anywhere.
plotly + h5py only, imported lazily so the rest of merge_candidates works without
them.
"""

from __future__ import annotations

import html
import os
from typing import Dict, List, Optional

import numpy as np

from .coverage_store import CHANNELS, CoverageStore
from .overlay import OverlayIndex, overlay_columns
from .utils import get_logger

# Channel display: strand is encoded by above/below the axis, so each read
# category gets one colour shared by its +/- channels. Colour-blind-safe trio.
_CHANNEL_COLORS = {
    "unique": "#2C6FBB",   # genuine unique signal (blue)
    "dup": "#E1943B",      # PCR/optical duplicates (amber = caution)
    "multi": "#8C8C8C",    # low-MAPQ / multimapper (grey = low quality)
}
_CATEGORY_ORDER = ("unique", "dup", "multi")  # stacking order, bottom -> top

_MAX_POINTS = 3500  # per-channel points before max-pool downsampling kicks in


# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #

def _parse_region_id(region_id: str):
    """('chrom:start-end', 0-based half-open) -> (chrom, start, end) or None."""
    try:
        chrom, span = region_id.rsplit(":", 1)
        start_s, end_s = span.split("-")
        return chrom, int(start_s), int(end_s)
    except (ValueError, AttributeError):
        return None


def _downsample_max(x: np.ndarray, series: Dict[str, np.ndarray], max_points: int):
    """Max-pool every channel into <= max_points bins, preserving peaks.

    Coverage peaks are the whole point of the plot, so we bin by max (not mean):
    a 1 bp spike survives downsampling instead of being averaged away.
    """
    n = x.shape[0]
    if n <= max_points:
        return x, series
    bin_size = int(np.ceil(n / max_points))
    n_bins = int(np.ceil(n / bin_size))
    pad = n_bins * bin_size - n
    x_out = x[:: bin_size][:n_bins]
    out: Dict[str, np.ndarray] = {}
    for name, arr in series.items():
        a = arr
        if pad:
            a = np.concatenate([a, np.zeros(pad, dtype=a.dtype)])
        out[name] = a.reshape(n_bins, bin_size).max(axis=1)
    return x_out, out


def _fmt(v, nd=3):
    if v is None:
        return "NA"
    if isinstance(v, float):
        return f"{v:.{nd}g}"
    return str(v)


# --------------------------------------------------------------------------- #
# Reference gene track (optional)
# --------------------------------------------------------------------------- #

def _load_reference_genes(reference_gtf: Optional[str], logger) -> Optional[OverlayIndex]:
    if not reference_gtf:
        return None
    try:
        from .overlay import parse_overlay
        idx = parse_overlay(reference_gtf, label="gene", feature_types=["gene"])
        logger.info("Loaded %d reference genes for the plot feature track.",
                    idx.n_features)
        return idx
    except Exception as exc:  # best-effort: a missing/odd GTF must not kill plots
        logger.warning("Could not load reference genes for plots: %s", exc)
        return None


# --------------------------------------------------------------------------- #
# Figure
# --------------------------------------------------------------------------- #

def _build_figure(
    region, chrom, win_start, win_end,
    stores: Dict[str, CoverageStore], sample_order: List[str],
    member_by_sample: Dict[str, tuple],
    indexes: List[OverlayIndex], gene_index: Optional[OverlayIndex],
):
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go

    n = win_end - win_start
    x_full = np.arange(win_start, win_end, dtype=np.int64) + 1  # 1-based display

    feature_row = len(sample_order) + 1
    n_rows = feature_row
    row_heights = [1.0] * len(sample_order) + [0.6]
    titles = [f"Sample: {s}" for s in sample_order] + ["Annotation / gene features"]
    fig = make_subplots(
        rows=n_rows, cols=1, shared_xaxes=True, vertical_spacing=0.02,
        row_heights=row_heights,
        subplot_titles=titles,
    )
    # Left-align the per-sample subplot titles so they read as row labels.
    for ann in fig.layout.annotations:
        ann.update(x=0, xanchor="left", font=dict(size=13))

    first_legend = True
    global_ymax = 0.0
    called_markers: List[tuple] = []  # (row, member_start, member_end)
    for si, sample in enumerate(sample_order, start=1):
        store = stores.get(sample)
        if store is not None:
            dense = store.window(chrom, win_start, win_end)  # zeros if chrom absent
        else:
            dense = {c: np.zeros(n, dtype=np.int32) for c in CHANNELS}
        x, series = _downsample_max(x_full, dense, _MAX_POINTS)

        ymax = 0
        for sign, suffix, mult in (("p", "_plus", 1), ("m", "_minus", -1)):
            stacked = np.zeros(x.shape[0])
            strand_sym = "+" if sign == "p" else "−"
            for cat in _CATEGORY_ORDER:
                vals = series[f"{cat}{suffix}"].astype(float)
                stacked = stacked + vals
                fig.add_trace(
                    go.Scatter(
                        x=x, y=mult * vals, mode="lines",
                        line=dict(width=0, color=_CHANNEL_COLORS[cat]),
                        fillcolor=_CHANNEL_COLORS[cat],
                        stackgroup=f"r{si}{sign}", name=cat,
                        legendgroup=cat, showlegend=first_legend and sign == "p",
                        customdata=vals,
                        # Unified hover (set on the layout) labels each row by its
                        # own category+strand and shows that channel's own depth,
                        # so a stacked pile-up reads as unique / dup / multi at once.
                        hovertemplate=(f"{cat} {strand_sym}: "
                                       "%{customdata:.0f}<extra></extra>"),
                    ),
                    row=si, col=1,
                )
            ymax = max(ymax, float(stacked.max()) if stacked.size else 0.0)
        first_legend = False
        global_ymax = max(global_ymax, ymax)

        mem = member_by_sample.get(sample)
        if mem is not None:
            called_markers.append((si, mem[0], mem[1]))
        fig.update_yaxes(title_text="depth", row=si, col=1, zeroline=True,
                         zerolinecolor="#BBBBBB")

    # Shared depth scale across every sample panel so equal depth = equal height
    # (IGV-style): a 14x sample must look shorter than a 30x one. Zero-centred so
    # + and - strands share the scale too.
    M = global_ymax * 1.05 if global_ymax > 0 else 1.0
    for si in range(1, len(sample_order) + 1):
        fig.update_yaxes(range=[-M, M], row=si, col=1)
        yid = "y" if si == 1 else f"y{si}"
        xid = "x" if si == 1 else f"x{si}"
        # Mark (near the top of the panel, in domain coords so it is independent
        # of the shared scale) where each sample independently called the region.
        for row, ms, me in called_markers:
            if row != si:
                continue
            fig.add_shape(
                type="rect", x0=ms + 1, x1=me, y0=0.9, y1=1.0,
                xref=xid, yref=f"{yid} domain",
                fillcolor="#2C6FBB", opacity=0.35, line_width=0,
            )

    # Feature track: reference genes on the top lane, each overlay source below.
    lanes: List[tuple] = []
    if gene_index is not None:
        lanes.append(("gene", gene_index, "#5B8C5A"))
    palette = ["#B5651D", "#7D3C98", "#C0392B", "#16A085", "#2E4053"]
    for k, idx in enumerate(indexes):
        lanes.append((idx.label, idx, palette[k % len(palette)]))

    n_lanes = max(1, len(lanes))
    for li, (label, idx, color) in enumerate(lanes):
        lane_top = n_lanes - li
        lane_bot = lane_top - 0.7
        feats = idx.overlap(chrom, win_start, win_end)
        for f in feats:
            fs = max(f.start, win_start)
            fe = min(f.end, win_end)
            fig.add_shape(
                type="rect", x0=fs + 1, x1=fe, y0=lane_bot, y1=lane_top,
                fillcolor=color, opacity=0.55, line=dict(color=color, width=1),
                row=feature_row, col=1,
            )
            fig.add_trace(
                go.Scatter(
                    x=[(fs + fe) / 2 + 0.5], y=[(lane_top + lane_bot) / 2],
                    mode="markers", marker=dict(size=4, color=color),
                    showlegend=False,
                    hovertext=[f"{label}: {f.ftype} {f.fid} "
                               f"[{f.start + 1}-{f.end}] {f.strand}"],
                    hoverinfo="text",
                ),
                row=feature_row, col=1,
            )
    fig.update_yaxes(
        row=feature_row, col=1, range=[0, n_lanes + 0.2],
        tickvals=[n_lanes - i - 0.35 for i in range(len(lanes))],
        ticktext=[lab for lab, _, _ in lanes] or ["(none)"],
    )
    fig.update_xaxes(title_text=f"{chrom} position (1-based)", row=feature_row, col=1)

    # Locus highlight across every row.
    fig.add_vrect(
        x0=region.start + 1, x1=region.end,
        fillcolor="#F2C14E", opacity=0.12, line_width=0, layer="below",
    )

    height = 170 * len(sample_order) + 240
    n_samples = len(sample_order)
    fig.update_layout(
        height=height, template="plotly_white",
        margin=dict(l=70, r=30, t=90, b=50),
        legend=dict(orientation="h", yanchor="bottom", y=1.015, x=0,
                    title_text="read category:"),
        title=dict(
            text=(f"{region.consensus_transcript_name or region.consensus_id}"
                  f"<br><span style='font-size:12px;color:#666'>"
                  f"one coverage panel per sample ({n_samples} samples) + a feature "
                  f"track; + strand above / − strand below the axis</span>"),
            y=0.985, yanchor="top",
        ),
        hovermode="x unified",
    )
    return fig


# --------------------------------------------------------------------------- #
# Evidence panel + page assembly
# --------------------------------------------------------------------------- #

def _evidence_html(region, labels, member_rows, strandedness) -> str:
    e = html.escape
    locus = f"{region.chrom}:{region.start + 1:,}-{region.end:,} ({region.strand})"
    summary_rows = [
        ("Transcript", region.consensus_transcript_name or "NA"),
        ("Consensus id", region.consensus_id),
        ("Locus", locus),
        ("Length (bp)", f"{region.end - region.start:,}"),
        ("Consensus class", region.consensus_class),
        ("Majority class", region.majority_class),
        ("Class agreement", _fmt(region.class_agreement)),
        ("Samples / members", f"{region.n_samples} / {region.n_members}"),
        ("Library strandedness", strandedness or "unknown"),
        ("Mean unique fraction", _fmt(region.mean_unique_fraction)),
        ("Mean dual-strand fraction", _fmt(region.mean_dual_strand_fraction)),
        ("Mean profile correlation", _fmt(region.mean_profile_correlation)),
        ("Mean avg depth", _fmt(region.mean_avg_depth)),
    ]
    # Overlay annotation columns (feature evidence) if present.
    for label in labels or []:
        for col in overlay_columns(label):
            val = region.annotations.get(col)
            if val not in (None, "", "NA", 0):
                summary_rows.append((col, _fmt(val)))

    summ = "".join(
        f"<tr><th>{e(str(k))}</th><td>{e(str(v))}</td></tr>" for k, v in summary_rows
    )

    # Per-sample evidence table.
    cols = [
        ("sample", "sample"), ("called", "called"), ("class", "class"),
        ("avg_depth", "avg_depth"), ("max_depth", "max_depth"),
        ("covered_fraction", "cov_frac"), ("unique_fraction", "uniq_frac"),
        ("dual_strand_fraction", "dual"), ("profile_correlation", "prof_corr"),
        ("dominant_strand", "strand"), ("context_label", "context"),
        ("nearest_feature_id", "nearest"), ("nearest_feature_distance", "dist"),
        ("reason_for_classification", "reason"),
    ]
    head = "".join(f"<th>{e(h)}</th>" for _, h in cols)
    body = []
    for row in member_rows:
        cells = []
        for key, _ in cols:
            cells.append(f"<td>{e(_fmt(row.get(key)))}</td>")
        body.append(f"<tr>{''.join(cells)}</tr>")
    per_sample = f"<table class='ev'><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table>"

    return (
        "<div class='evidence'>"
        "<h2>Evidence</h2>"
        f"<table class='ev summary'>{summ}</table>"
        "<h3>Per-sample metrics (samples that called this locus)</h3>"
        f"{per_sample}"
        "</div>"
    )


_PAGE_CSS = """
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
     margin:0;padding:16px 20px;color:#1b1b1b;background:#ffffff;}
h1{font-size:20px;margin:0 0 4px;} h2{font-size:16px;margin:18px 0 6px;}
h3{font-size:13px;margin:14px 0 4px;color:#444;}
.sub{color:#666;font-size:13px;margin-bottom:8px;}
table.ev{border-collapse:collapse;font-size:12px;margin:4px 0;}
table.ev th,table.ev td{border:1px solid #e2e2e2;padding:3px 7px;text-align:left;
     vertical-align:top;}
table.ev.summary th{background:#f6f7f9;white-space:nowrap;}
table.ev thead th{background:#f0f3f7;position:sticky;top:0;}
.evidence{overflow-x:auto;} .legend-note{font-size:12px;color:#555;margin:6px 0;}
a{color:#2C6FBB;}
"""


def _write_page(path, region, fig_html, evidence_html):
    name = html.escape(region.consensus_transcript_name or region.consensus_id)
    doc = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>{name}</title>
<style>{_PAGE_CSS}</style></head><body>
<h1>{name}</h1>
<div class="sub">{html.escape(region.chrom)}:{region.start + 1:,}-{region.end:,}
 ({html.escape(region.strand)}) &middot; {html.escape(region.consensus_class)}
 &middot; {region.n_samples} samples</div>
<div class="legend-note"><b>How to read this:</b> each stacked panel below is
 <b>one sample's</b> per-base coverage across this locus (&plusmn; shoulder); the last,
 shorter panel is the annotation / reference-gene track. Within a panel, coverage on
 the <b>+</b> strand is drawn above the axis and the <b>&minus;</b> strand below, and
 each is split by read category: <b>blue = unique</b> signal, <b>amber = PCR
 duplicates</b>, <b>grey = low-MAPQ / multimapper</b> &mdash; so a tall amber or grey
 stack means the pileup is inflated by duplicate or low-quality reads rather than
 genuine signal. The yellow band marks the consensus locus; a blue bar over a panel
 marks where that sample independently called the region.</div>
{fig_html}
{evidence_html}
</body></html>
"""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(doc)


def _write_index(out_dir, entries):
    rows = []
    for entry in entries:
        rows.append(
            "<tr>"
            f"<td><a href='{html.escape(entry['file'])}'>{html.escape(entry['name'])}</a></td>"
            f"<td>{html.escape(entry['locus'])}</td>"
            f"<td>{html.escape(entry['class'])}</td>"
            f"<td>{entry['n_samples']}</td>"
            f"<td>{html.escape(entry['features'])}</td>"
            "</tr>"
        )
    doc = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Consensus transcript coverage plots</title>
<style>{_PAGE_CSS}</style></head><body>
<h1>Consensus transcript coverage plots</h1>
<div class="sub">{len(entries)} transcript(s)</div>
<table class="ev"><thead><tr><th>transcript</th><th>locus</th><th>class</th>
<th>samples</th><th>features</th></tr></thead>
<tbody>{''.join(rows)}</tbody></table>
</body></html>
"""
    with open(os.path.join(out_dir, "index.html"), "w", encoding="utf-8") as fh:
        fh.write(doc)


# --------------------------------------------------------------------------- #
# Entry point
# --------------------------------------------------------------------------- #

def build_plots(
    targets, df, indexes, labels, cfg,
    store_by_sample: Dict[str, str], sample_order: List[str], out_dir: str,
) -> int:
    """Write one IGV-like HTML per target region into ``out_dir``.

    ``targets`` are the surviving ``ConsensusRegion`` objects to plot; ``df`` is the
    full per-sample candidate table (polars); ``indexes``/``labels`` are the overlay
    sources; ``store_by_sample`` maps sample name -> coverage-store path. Returns the
    number of plots written.
    """
    logger = get_logger(getattr(cfg, "verbose", False))
    try:
        import plotly  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "plotly is required for --emit-plots but is not installed. "
            "Install it with `pip install plotly`."
        ) from exc

    os.makedirs(out_dir, exist_ok=True)

    # By default pages load plotly.js from the CDN (small pages, needs network).
    # With --plot-offline we write plotly.min.js once into out_dir and every page
    # references it relatively (include_plotlyjs="directory" only *references* it,
    # so we must write the file ourselves), keeping the pages fully offline.
    offline = bool(getattr(cfg, "plot_offline", False))
    plotlyjs_mode = "directory" if offline else "cdn"
    if offline:
        try:
            from plotly.offline import get_plotlyjs
            with open(os.path.join(out_dir, "plotly.min.js"), "w", encoding="utf-8") as fh:
                fh.write(get_plotlyjs())
        except Exception as exc:  # pragma: no cover
            logger.warning("Could not write plotly.min.js (%s); pages will still open "
                           "if the library is cached.", exc)

    # Open every sample's coverage store once (cached slicing across all loci).
    stores: Dict[str, CoverageStore] = {}
    strandedness = "unknown"
    for sample in sample_order:
        path = store_by_sample.get(sample)
        if not path or not os.path.exists(path):
            logger.warning("No coverage store for sample %r; it will show as empty.",
                           sample)
            continue
        try:
            store = CoverageStore(path)
            stores[sample] = store
            if store.strandedness and store.strandedness != "unknown":
                strandedness = store.strandedness
        except Exception as exc:
            logger.warning("Could not open coverage store %r: %s", path, exc)

    gene_index = _load_reference_genes(getattr(cfg, "reference_gtf", None), logger)

    # Row lookup for per-sample member metrics: (sample, region_id) -> row dict.
    row_by_key: Dict[tuple, dict] = {}
    for r in df.to_dicts():
        row_by_key[(str(r.get("sample")), str(r.get("region_id")))] = r

    shoulder = int(getattr(cfg, "plot_shoulder", 1000))
    entries = []
    n_written = 0
    try:
        for region in targets:
            chrom = region.chrom
            win_start = max(0, region.start - shoulder)
            win_end = region.end + shoulder

            # Per-sample called interval (from the exact member region_id).
            member_by_sample: Dict[str, tuple] = {}
            member_rows = []
            for sample, rid in zip(region.member_samples, region.member_region_ids):
                parsed = _parse_region_id(rid)
                if parsed and parsed[0] == chrom:
                    member_by_sample[sample] = (parsed[1], parsed[2])
                row = row_by_key.get((str(sample), str(rid)))
                if row is not None:
                    member_rows.append(row)

            fig = _build_figure(
                region, chrom, win_start, win_end, stores, sample_order,
                member_by_sample, indexes, gene_index,
            )
            fig_html = fig.to_html(full_html=False, include_plotlyjs=plotlyjs_mode,
                                   default_width="100%")
            evidence_html = _evidence_html(region, labels, member_rows, strandedness)

            base = region.consensus_transcript_name or region.consensus_id.replace(":", "_")
            fname = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in base)
            fname = f"{fname}.html"
            _write_page(os.path.join(out_dir, fname), region, fig_html, evidence_html)

            feat_types = sorted({
                t for lbl in (labels or [])
                for t in (region.overlay_types.get(lbl) or {})
            })
            entries.append({
                "file": fname,
                "name": region.consensus_transcript_name or region.consensus_id,
                "locus": f"{region.chrom}:{region.start + 1}-{region.end} ({region.strand})",
                "class": region.consensus_class,
                "n_samples": region.n_samples,
                "features": ", ".join(feat_types) or "-",
            })
            n_written += 1

        _write_index(out_dir, entries)
    finally:
        for store in stores.values():
            store.close()

    return n_written
