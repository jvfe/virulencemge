#!/usr/bin/env python3
"""Aggregate per-sample results into summary matrices and an HTML dashboard.

MultiQC has no modules for ABricate/VFDB, PhiSpy or IslandPath-DIMOB, so this
script takes their place: it concatenates the per-sample tables produced by
``integrate_elements.py`` and renders a single self-contained HTML report with
the Plotly library inlined, so it can be opened without network access.
"""

from __future__ import annotations

import argparse
import datetime as dt
import html
import sys
from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from jinja2 import Template

CONTEXT_CORE = "Core genome"
CONTEXT_PAI = "PAI"
CONTEXT_PROPHAGE = "Prophage"
CONTEXT_BOTH = "PAI+Prophage"

CONTEXT_ORDER = [CONTEXT_CORE, CONTEXT_PAI, CONTEXT_PROPHAGE, CONTEXT_BOTH]
CONTEXT_CODES = {CONTEXT_CORE: 0, CONTEXT_PAI: 1, CONTEXT_PROPHAGE: 2, CONTEXT_BOTH: 3}
CONTEXT_COLOURS = {
    CONTEXT_CORE: "#9aa5b1",
    CONTEXT_PAI: "#e8973a",
    CONTEXT_PROPHAGE: "#3f8fd2",
    CONTEXT_BOTH: "#b5539c",
}

PLOTLY_LAYOUT = {
    "template": "plotly_white",
    "margin": {"l": 60, "r": 20, "t": 50, "b": 60},
    "font": {"family": "Inter, -apple-system, Segoe UI, Helvetica, Arial, sans-serif", "size": 13},
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--stats", nargs="+", required=True, type=Path, help="Per-sample *.stats.tsv files.")
    parser.add_argument("--annotations", nargs="+", required=True, type=Path, help="Per-sample *.virulence_annotation.tsv files.")
    parser.add_argument("--elements", nargs="*", default=[], type=Path, help="Per-sample *.mge_regions.tsv files.")
    parser.add_argument("--abricate-summary", type=Path, help="Matrix produced by `abricate --summary`.")
    parser.add_argument("--outdir", type=Path, default=Path("."), help="Output directory.")
    parser.add_argument("--prefix", default="virulencemge", help="Output file prefix.")
    parser.add_argument("--min-overlap", type=float, default=0.5, help="Containment threshold used upstream, reported in the HTML.")
    parser.add_argument("--pipeline-version", default="", help="Pipeline version shown in the report header.")
    parser.add_argument("--top-genes", type=int, default=30, help="Number of most frequent virulence genes to plot.")
    return parser.parse_args(argv)


def concat_tables(paths: list[Path]) -> pd.DataFrame:
    frames = []
    for path in paths:
        if not path.exists() or path.stat().st_size == 0:
            continue
        frame = pd.read_csv(path, sep="\t", dtype={"sample": str})
        if not frame.empty:
            frames.append(frame)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def figure_html(figure: go.Figure, include_js: bool) -> str:
    figure.update_layout(**PLOTLY_LAYOUT)
    return figure.to_html(
        full_html=False,
        include_plotlyjs=True if include_js else False,
        config={"displaylogo": False, "responsive": True},
    )


def plot_context_counts(annotations: pd.DataFrame) -> go.Figure:
    counts = (
        annotations.groupby(["sample", "mge_context"]).size().reset_index(name="count")
        if not annotations.empty
        else pd.DataFrame(columns=["sample", "mge_context", "count"])
    )
    figure = px.bar(
        counts,
        x="sample",
        y="count",
        color="mge_context",
        category_orders={"mge_context": CONTEXT_ORDER},
        color_discrete_map=CONTEXT_COLOURS,
        labels={"sample": "Genome", "count": "Virulence factor hits", "mge_context": "Genomic context"},
        title="Virulence factors per genome, by genomic context",
    )
    figure.update_layout(barmode="stack", xaxis_tickangle=-40)
    return figure


def plot_element_counts(stats: pd.DataFrame) -> go.Figure:
    melted = stats.melt(
        id_vars="sample",
        value_vars=["n_pai", "n_prophage"],
        var_name="element_type",
        value_name="count",
    )
    melted["element_type"] = melted["element_type"].map({"n_pai": "Genomic islands (PAI)", "n_prophage": "Prophages"})
    figure = px.bar(
        melted,
        x="sample",
        y="count",
        color="element_type",
        barmode="group",
        color_discrete_map={"Genomic islands (PAI)": CONTEXT_COLOURS[CONTEXT_PAI], "Prophages": CONTEXT_COLOURS[CONTEXT_PROPHAGE]},
        labels={"sample": "Genome", "count": "Predicted elements", "element_type": "Element type"},
        title="Predicted mobile genetic elements per genome",
    )
    figure.update_layout(xaxis_tickangle=-40)
    return figure


def plot_colocalisation_frequency(stats: pd.DataFrame) -> go.Figure:
    frame = stats.copy()
    frame["pct_in_mge"] = frame["frac_vf_in_mge"] * 100
    frame = frame.sort_values("pct_in_mge", ascending=False)
    figure = px.bar(
        frame,
        x="sample",
        y="pct_in_mge",
        hover_data=["n_virulence_factors", "n_vf_in_mge", "n_pai", "n_prophage"],
        labels={"sample": "Genome", "pct_in_mge": "Virulence factors inside MGEs (%)"},
        title="Co-localisation frequency: share of virulence factors within mobile genetic elements",
    )
    figure.update_traces(marker_color=CONTEXT_COLOURS[CONTEXT_BOTH])
    figure.update_layout(xaxis_tickangle=-40, yaxis_range=[0, 100])
    return figure


def plot_context_proportions(annotations: pd.DataFrame) -> go.Figure:
    if annotations.empty:
        return px.bar(title="Genomic context of virulence factors (proportions)")
    proportions = (
        annotations.groupby(["sample", "mge_context"]).size().rename("count").reset_index()
    )
    totals = proportions.groupby("sample")["count"].transform("sum")
    proportions["proportion"] = proportions["count"] / totals * 100
    figure = px.bar(
        proportions,
        x="sample",
        y="proportion",
        color="mge_context",
        category_orders={"mge_context": CONTEXT_ORDER},
        color_discrete_map=CONTEXT_COLOURS,
        labels={"sample": "Genome", "proportion": "Virulence factors (%)", "mge_context": "Genomic context"},
        title="Genomic context of virulence factors (proportions)",
    )
    figure.update_layout(barmode="stack", xaxis_tickangle=-40, yaxis_range=[0, 100])
    return figure


def build_localisation_matrix(annotations: pd.DataFrame, samples: list[str]) -> pd.DataFrame:
    """Gene x sample matrix whose cells state where the gene was found.

    A gene detected more than once in a genome is summarised by its most
    informative context, so that a single MGE-borne copy is never masked by
    core-genome copies. Genomes without any hit are kept as empty columns so the
    matrix always covers every sample in the run.
    """
    if annotations.empty:
        return pd.DataFrame()

    ranked = annotations.copy()
    ranked["context_rank"] = ranked["mge_context"].map(CONTEXT_CODES).fillna(0)
    best = ranked.sort_values("context_rank").groupby(["gene", "sample"], as_index=False).last()
    matrix = best.pivot(index="gene", columns="sample", values="mge_context")
    return matrix.reindex(columns=samples).sort_index()


def plot_localisation_heatmap(matrix: pd.DataFrame, top_genes: int) -> go.Figure:
    if matrix.empty:
        return px.imshow([[0]], title="Virulence factor localisation across genomes")

    presence = matrix.notna().sum(axis=1).sort_values(ascending=False)
    selected = matrix.loc[presence.head(top_genes).index].sort_index()
    coded = selected.apply(lambda column: column.map(CONTEXT_CODES)).astype(float)

    figure = go.Figure(
        go.Heatmap(
            z=coded.values,
            x=list(coded.columns),
            y=list(coded.index),
            text=selected.fillna("Absent").values,
            hovertemplate="Gene: %{y}<br>Genome: %{x}<br>Context: %{text}<extra></extra>",
            colorscale=[
                [0.00, CONTEXT_COLOURS[CONTEXT_CORE]],
                [0.33, CONTEXT_COLOURS[CONTEXT_PAI]],
                [0.66, CONTEXT_COLOURS[CONTEXT_PROPHAGE]],
                [1.00, CONTEXT_COLOURS[CONTEXT_BOTH]],
            ],
            zmin=0,
            zmax=3,
            colorbar={
                "tickmode": "array",
                "tickvals": [CONTEXT_CODES[label] for label in CONTEXT_ORDER],
                "ticktext": CONTEXT_ORDER,
                "title": "Context",
            },
            xgap=1,
            ygap=1,
        )
    )
    figure.update_layout(
        title=f"Virulence factor localisation across genomes (top {len(selected)} genes)",
        xaxis_title="Genome",
        yaxis_title="Virulence gene",
        height=max(400, 22 * len(selected) + 160),
        xaxis_tickangle=-40,
    )
    return figure


def plot_element_sizes(elements: pd.DataFrame) -> go.Figure:
    if elements.empty:
        return px.box(title="Size distribution of predicted mobile genetic elements")
    frame = elements.copy()
    frame["element_type"] = frame["element_type"].map({"pai": "Genomic islands (PAI)", "prophage": "Prophages"})
    figure = px.box(
        frame,
        x="element_type",
        y="length",
        color="element_type",
        points="all",
        hover_data=["sample", "element_id", "contig", "n_virulence_factors"],
        color_discrete_map={"Genomic islands (PAI)": CONTEXT_COLOURS[CONTEXT_PAI], "Prophages": CONTEXT_COLOURS[CONTEXT_PROPHAGE]},
        labels={"element_type": "Element type", "length": "Element length (bp)"},
        title="Size distribution of predicted mobile genetic elements",
    )
    figure.update_layout(showlegend=False)
    return figure


def plot_top_genes(annotations: pd.DataFrame, top_genes: int) -> go.Figure:
    if annotations.empty:
        return px.bar(title="Most frequently detected virulence genes")
    counts = (
        annotations.groupby(["gene", "mge_context"])["sample"]
        .nunique()
        .reset_index(name="n_genomes")
    )
    ordering = counts.groupby("gene")["n_genomes"].sum().sort_values(ascending=False).head(top_genes).index
    counts = counts[counts["gene"].isin(ordering)]
    figure = px.bar(
        counts,
        x="n_genomes",
        y="gene",
        color="mge_context",
        orientation="h",
        category_orders={"mge_context": CONTEXT_ORDER, "gene": list(reversed(list(ordering)))},
        color_discrete_map=CONTEXT_COLOURS,
        labels={"n_genomes": "Genomes", "gene": "Virulence gene", "mge_context": "Genomic context"},
        title=f"Most frequently detected virulence genes (top {len(ordering)})",
    )
    figure.update_layout(barmode="stack", height=max(400, 20 * len(ordering) + 160))
    return figure


REPORT_TEMPLATE = Template(
    """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{{ title }}</title>
<style>
  :root { --fg: #16212e; --muted: #5b6b7c; --line: #e2e8ef; --accent: #b5539c; }
  * { box-sizing: border-box; }
  body { margin: 0; background: #f6f8fa; color: var(--fg);
         font-family: Inter, -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif; line-height: 1.55; }
  header { background: linear-gradient(135deg, #16212e, #2c4a63); color: #fff; padding: 2.4rem 2rem 2rem; }
  header h1 { margin: 0 0 .35rem; font-size: 1.6rem; font-weight: 650; }
  header p { margin: 0; opacity: .82; font-size: .9rem; }
  main { max-width: 1180px; margin: 0 auto; padding: 1.6rem 1.2rem 4rem; }
  section { background: #fff; border: 1px solid var(--line); border-radius: 10px; padding: 1.3rem 1.4rem; margin-bottom: 1.2rem; }
  section > h2 { margin: 0 0 .3rem; font-size: 1.12rem; }
  section > p.hint { margin: 0 0 1rem; color: var(--muted); font-size: .875rem; }
  .cards { display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: .8rem; margin-bottom: 1.2rem; }
  .card { background: #fff; border: 1px solid var(--line); border-radius: 10px; padding: 1rem 1.1rem; }
  .card .value { font-size: 1.7rem; font-weight: 660; letter-spacing: -.02em; }
  .card .label { color: var(--muted); font-size: .78rem; text-transform: uppercase; letter-spacing: .05em; }
  .scroll { overflow-x: auto; }
  table { border-collapse: collapse; width: 100%; font-size: .84rem; }
  th, td { padding: .45rem .6rem; border-bottom: 1px solid var(--line); text-align: right; white-space: nowrap; }
  th { position: sticky; top: 0; background: #f0f3f7; text-align: right; font-weight: 600; }
  th:first-child, td:first-child { text-align: left; }
  tbody tr:hover { background: #fbfcfe; }
  .note { background: #fff8e6; border: 1px solid #f0d9a0; border-radius: 8px; padding: .8rem 1rem; font-size: .86rem; }
  footer { text-align: center; color: var(--muted); font-size: .8rem; padding-bottom: 2rem; }
  code { background: #eef2f6; padding: .1rem .3rem; border-radius: 4px; font-size: .85em; }
</style>
</head>
<body>
<header>
  <h1>{{ title }}</h1>
  <p>{{ n_samples }} genome(s) &middot; generated {{ generated_at }}{% if pipeline_version %} &middot; nf-core/virulencemge {{ pipeline_version }}{% endif %}</p>
</header>
<main>
  <div class="cards">
    {% for card in cards %}
    <div class="card"><div class="value">{{ card.value }}</div><div class="label">{{ card.label }}</div></div>
    {% endfor %}
  </div>

  {% if warnings %}
  <section>
    <h2>Quality warnings</h2>
    {% for warning in warnings %}<div class="note">{{ warning }}</div>{% endfor %}
  </section>
  {% endif %}

  <section>
    <h2>Methods</h2>
    <p class="hint">
      Virulence factors were detected with ABricate against VFDB, prophages with PhiSpy and genomic
      islands with IslandPath-DIMOB. A virulence factor is reported as contained in an element when at
      least <strong>{{ min_overlap_pct }}%</strong> of the gene falls inside the predicted element.
      Genes overlapping both a genomic island and a prophage are reported as <code>PAI+Prophage</code>.
    </p>
  </section>

  {% for block in sections %}
  <section>
    <h2>{{ block.title }}</h2>
    {% if block.hint %}<p class="hint">{{ block.hint }}</p>{% endif %}
    {{ block.body }}
  </section>
  {% endfor %}
</main>
<footer>Generated by nf-core/virulencemge &middot; tables for downstream analysis are written alongside this report.</footer>
</body>
</html>
"""
)


def dataframe_to_html(frame: pd.DataFrame, max_rows: int = 200) -> str:
    if frame.empty:
        return '<p class="hint">No rows to display.</p>'
    truncated = frame.head(max_rows)
    table = truncated.to_html(index=False, border=0, na_rep="", float_format=lambda value: f"{value:g}")
    suffix = (
        f'<p class="hint">Showing the first {max_rows} of {len(frame)} rows; the full table is in the accompanying TSV.</p>'
        if len(frame) > max_rows
        else ""
    )
    return f'<div class="scroll">{table}</div>{suffix}'


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    stats = concat_tables(args.stats)
    annotations = concat_tables(args.annotations)
    elements = concat_tables(list(args.elements))

    if stats.empty:
        sys.exit("ERROR: no per-sample statistics were found, cannot build a summary report")
    stats = stats.sort_values("sample").reset_index(drop=True)

    localisation = build_localisation_matrix(annotations, stats["sample"].tolist())
    presence = (localisation.notna().astype(int) if not localisation.empty else pd.DataFrame())

    stats.to_csv(args.outdir / f"{args.prefix}_sample_summary.tsv", sep="\t", index=False)
    if not localisation.empty:
        localisation.fillna("Absent").to_csv(args.outdir / f"{args.prefix}_vf_localisation_matrix.tsv", sep="\t")
        presence.to_csv(args.outdir / f"{args.prefix}_vf_presence_matrix.tsv", sep="\t")
    if not annotations.empty:
        annotations.to_csv(args.outdir / f"{args.prefix}_virulence_annotation.tsv", sep="\t", index=False)
    if not elements.empty:
        elements.to_csv(args.outdir / f"{args.prefix}_mge_regions.tsv", sep="\t", index=False)

    total_vf = int(stats["n_virulence_factors"].sum())
    total_in_mge = int(stats["n_vf_in_mge"].sum())
    cards = [
        {"label": "Genomes", "value": len(stats)},
        {"label": "Virulence factor hits", "value": f"{total_vf:,}"},
        {"label": "Unique VF genes", "value": f"{annotations['gene'].nunique() if not annotations.empty else 0:,}"},
        {"label": "Genomic islands", "value": f"{int(stats['n_pai'].sum()):,}"},
        {"label": "Prophages", "value": f"{int(stats['n_prophage'].sum()):,}"},
        {"label": "VFs inside MGEs", "value": f"{(total_in_mge / total_vf * 100):.1f}%" if total_vf else "n/a"},
    ]

    warnings = []
    unresolved = int(stats.get("n_unresolved_contig_ids", pd.Series([0])).sum())
    if unresolved:
        warnings.append(
            f"{unresolved} contig identifier(s) reported by PhiSpy or IslandPath-DIMOB could not be matched to "
            "the assembly FASTA. Elements on those contigs were not cross-referenced; check the "
            "integrate_elements logs in the pipeline work directory."
        )
    empty_samples = stats.loc[stats["n_virulence_factors"] == 0, "sample"].tolist()
    if empty_samples:
        warnings.append(f"No VFDB hits were found for: {', '.join(empty_samples)}.")

    # Plotly is inlined once, in the first figure, to keep the report standalone.
    figures = [
        ("Virulence factor burden", "Absolute counts of ABricate/VFDB hits, split by genomic context.", plot_context_counts(annotations)),
        ("Genomic context proportions", "The same data normalised per genome.", plot_context_proportions(annotations)),
        ("Mobile genetic elements", "Counts of PhiSpy prophages and IslandPath-DIMOB genomic islands.", plot_element_counts(stats)),
        ("Co-localisation frequency", "Share of virulence factors residing in a prophage or genomic island.", plot_colocalisation_frequency(stats)),
        ("Element sizes", "Length distribution of predicted elements, one point per element.", plot_element_sizes(elements)),
        ("Localisation matrix", "Where each virulence gene was found in each genome.", plot_localisation_heatmap(localisation, args.top_genes)),
        ("Most frequent virulence genes", "Number of genomes carrying each gene, split by context.", plot_top_genes(annotations, args.top_genes)),
    ]

    sections = []
    for index, (title, hint, figure) in enumerate(figures):
        sections.append({"title": title, "hint": hint, "body": figure_html(figure, include_js=index == 0)})

    sections.append(
        {
            "title": "Per-genome summary",
            "hint": "Also written to "
            f"{args.prefix}_sample_summary.tsv.",
            "body": dataframe_to_html(stats),
        }
    )
    if not elements.empty:
        sections.append(
            {
                "title": "Predicted elements carrying virulence factors",
                "hint": "Elements with at least one contained virulence gene.",
                "body": dataframe_to_html(
                    elements[elements["n_virulence_factors"] > 0].sort_values(
                        ["sample", "n_virulence_factors"], ascending=[True, False]
                    )
                ),
            }
        )
    if args.abricate_summary and args.abricate_summary.exists():
        abricate_summary = pd.read_csv(args.abricate_summary, sep="\t")
        abricate_summary.to_csv(args.outdir / f"{args.prefix}_abricate_summary.tsv", sep="\t", index=False)
        sections.append(
            {
                "title": "ABricate presence/absence matrix",
                "hint": "Raw output of `abricate --summary` across all genomes.",
                "body": dataframe_to_html(abricate_summary, max_rows=50),
            }
        )

    report = REPORT_TEMPLATE.render(
        title="Virulence factors and mobile genetic elements",
        n_samples=len(stats),
        generated_at=dt.datetime.now().strftime("%Y-%m-%d %H:%M"),
        pipeline_version=html.escape(args.pipeline_version),
        min_overlap_pct=f"{args.min_overlap * 100:g}",
        cards=cards,
        warnings=warnings,
        sections=sections,
    )

    report_path = args.outdir / f"{args.prefix}_report.html"
    report_path.write_text(report, encoding="utf-8")
    print(f"Wrote {report_path} covering {len(stats)} genomes", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
