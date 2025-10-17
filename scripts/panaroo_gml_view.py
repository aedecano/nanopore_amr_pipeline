#!/usr/bin/env python3
"""
panaroo_gml_view.py — visualize Panaroo .gml with communities & hubs

Outputs (into --outdir):
  - LCC_nodes.tsv, LCC_edges.tsv
  - LCC_plot.png (qualitative colors for top-N communities, hubs highlighted)
  - LCC.graphml (if possible), LCC.gexf (if possible)

Example:
  python panaroo_gml_view.py results/final_graph.gml \
    --outdir ./gml_view --max-nodes 800 --legend-top-n 10 \
    --layout spring --hub-quantile 0.95 --outline
"""

import argparse
import json
import os
import sys
from collections import Counter

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from networkx.algorithms import community


# --------------------------- helpers ---------------------------

def sanitize_graph_for_export(Gin: nx.Graph) -> nx.Graph:
    """Convert non-scalar attrs to JSON strings so GraphML/GEXF exports don't fail."""
    Gout = Gin.copy()

    def _clean(v):
        if isinstance(v, (str, int, float, bool)) or v is None:
            return v
        return json.dumps(v, ensure_ascii=False)

    for n, attrs in list(Gout.nodes(data=True)):
        Gout.nodes[n].update({k: _clean(v) for k, v in attrs.items()})
    for u, v, attrs in list(Gout.edges(data=True)):
        Gout.edges[u, v].update({k: _clean(vv) for k, vv in attrs.items()})
    Gout.graph.update({k: _clean(v) for k, v in Gout.graph.items()})
    return Gout


def compute_layout(G: nx.Graph, layout: str):
    if layout == "kamada":
        return nx.kamada_kawai_layout(G, weight=None)
    # default: spring (force-directed)
    return nx.spring_layout(G, seed=42, k=None, iterations=100)


def plot_lcc_with_enhancements(
    LCC: nx.Graph,
    node2c: dict,
    out_png: str,
    *,
    max_nodes: int,
    top_n: int,
    layout: str,
    hub_quantile: float,
    draw_outline: bool,
):
    """
    Plot LCC with:
      - qualitative colors for top-N largest communities (others grey)
      - spring/kamada layout
      - log-scale node sizes by degree
      - high-degree hubs highlighted
      - optional centroid outlines for top-N communities
    """
    # Subsample for readability
    nodes_to_plot = list(LCC.nodes())[:min(max_nodes, LCC.number_of_nodes())]
    LP = LCC.subgraph(nodes_to_plot).copy()

    fig, ax = plt.subplots(figsize=(11, 10))

    # ---- Layout & visuals
    pos = compute_layout(LP, layout)
    degs = dict(LP.degree())

    # Log-scale sizes (stable across sparse graphs)
    sizes = [20 + 100 * np.log1p(degs[n]) for n in LP.nodes()]

    # Community assignments & sizes from full LCC
    comms_full = {n: node2c.get(n, -1) for n in LCC.nodes()}
    counts = Counter(comms_full.values())
    unique_comms = sorted(set(node2c.values()))
    n_comms = len(unique_comms) if unique_comms else 1

    # Top-N communities by size (rest = Other)
    top_comms = [cid for cid, _ in counts.most_common(max(0, top_n))]
    comms_lp = np.array([node2c.get(n, -1) for n in LP.nodes()])
    top_mask = np.isin(comms_lp, top_comms)
    comms_mapped = np.where(top_mask, comms_lp, -1)

    # Qualitative palette for top-N; grey for Other
    cmap = mpl.colormaps.get_cmap("tab10")
    color_map_dict = {cid: cmap(i % 10) for i, cid in enumerate(top_comms)}
    other_color = (0.75, 0.75, 0.75, 0.65)
    colors = [color_map_dict.get(int(c), other_color) for c in comms_mapped]

    # ---- Draw network
    nx.draw_networkx_edges(LP, pos, alpha=0.25, width=0.5, ax=ax)
    nx.draw_networkx_nodes(LP, pos, node_size=sizes, node_color=colors, ax=ax)

    # ---- Highlight high-degree hubs (top quantile)
    if 0.0 < hub_quantile < 1.0:
        deg_values = np.array([degs[n] for n in LP.nodes()])
        thresh = float(np.quantile(deg_values, hub_quantile)) if len(deg_values) else np.inf
        hub_nodes = [n for n in LP.nodes() if degs[n] >= thresh and degs[n] > 0]
        if hub_nodes:
            nx.draw_networkx_nodes(
                LP, pos, nodelist=hub_nodes,
                node_size=[min(200, 40 + 20 * np.log1p(degs[n])) for n in hub_nodes],
                node_color="black", ax=ax, label="High-degree hubs"
            )

    # ---- Optional centroid outlines per top community
    if draw_outline and top_comms:
        for cid in top_comms:
            nodes_c = [n for n in LP.nodes() if node2c.get(n) == cid]
            if len(nodes_c) < 5:
                continue
            xs = [pos[n][0] for n in nodes_c]
            ys = [pos[n][1] for n in nodes_c]
            ax.scatter(np.mean(xs), np.mean(ys), s=160, facecolor='none',
                       edgecolor=color_map_dict[cid], linewidth=2, alpha=0.9)

    # ---- Legend: top-N + Other
    handles, labels = [], []
    for cid in top_comms:
        handles.append(plt.Line2D([0], [0], marker='o', linestyle='',
                                  markerfacecolor=color_map_dict[cid],
                                  markeredgecolor='none', markersize=8))
        labels.append(f"Community {cid} (n={counts[cid]})")
    if np.any(~top_mask):
        handles.append(plt.Line2D([0], [0], marker='o', linestyle='',
                                  markerfacecolor=other_color,
                                  markeredgecolor='none', markersize=8))
        labels.append(f"Other (sum n={sum(counts[cid] for cid in counts if cid not in top_comms)})")
    if handles:
        ax.legend(handles, labels, loc='lower left', frameon=False, fontsize=8, ncol=1)

    ax.set_axis_off()
    ax.set_title(
        f"Panaroo LCC (showing {LP.number_of_nodes()} nodes / {LCC.number_of_nodes()} total)\n"
        f"Communities detected: {n_comms} — legend shows top {min(top_n, n_comms)}",
        fontsize=11,
    )

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


# --------------------------- main ---------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gml", help="Path to Panaroo .gml file")
    ap.add_argument("--outdir", default="gml_view", help="Output directory")
    ap.add_argument("--max-nodes", type=int, default=1000,
                    help="Max nodes to plot for readability")
    ap.add_argument("--legend-top-n", type=int, default=10,
                    help="Number of largest communities to label in the legend")
    ap.add_argument("--layout", choices=["spring", "kamada"], default="spring",
                    help="Graph layout for plotting")
    ap.add_argument("--hub-quantile", type=float, default=0.95,
                    help="Quantile threshold to mark hubs (0-1); e.g. 0.95 marks top 5%%. Set 1.0 to disable.")
    ap.add_argument("--outline", action="store_true",
                    help="Draw centroid outlines for top-N communities")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"[i] Loading GML: {args.gml}")
    try:
        G = nx.read_gml(args.gml)
    except Exception as e:
        print(f"[x] Failed to read GML: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"[i] Graph loaded: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # Largest Connected Component (treat as undirected)
    H_und = G.to_undirected() if nx.is_directed(G) else G
    comps = sorted((H_und.subgraph(c).copy() for c in nx.connected_components(H_und)),
                   key=lambda g: g.number_of_nodes(), reverse=True)
    if not comps:
        print("[x] No connected components found.", file=sys.stderr)
        sys.exit(1)
    LCC = comps[0]
    print(f"[i] LCC: {LCC.number_of_nodes()} nodes, {LCC.number_of_edges()} edges")

    # Communities (greedy modularity)
    print("[i] Detecting communities on LCC (greedy modularity)…")
    cmtys = list(community.greedy_modularity_communities(LCC))
    node2c = {n: cid for cid, nodes in enumerate(cmtys) for n in nodes}
    nx.set_node_attributes(LCC, node2c, "community")

    # TSV exports
    nodes_tsv = os.path.join(args.outdir, "LCC_nodes.tsv")
    edges_tsv = os.path.join(args.outdir, "LCC_edges.tsv")
    with open(nodes_tsv, "w") as f:
        f.write("node\tdegree\tcommunity\n")
        for n in LCC.nodes():
            f.write(f"{n}\t{LCC.degree(n)}\t{node2c.get(n, -1)}\n")
    with open(edges_tsv, "w") as f:
        f.write("source\ttarget\n")
        for u, v in LCC.edges():
            f.write(f"{u}\t{v}\n")
    print(f"[i] Wrote: {nodes_tsv}")
    print(f"[i] Wrote: {edges_tsv}")

    # Plot with enhancements
    out_png = os.path.join(args.outdir, "LCC_plot.png")
    plot_lcc_with_enhancements(
        LCC, node2c, out_png,
        max_nodes=args.max_nodes,
        top_n=args.legend_top_n,
        layout=args.layout,
        hub_quantile=args.hub_quantile,
        draw_outline=args.outline,
    )
    print(f"[i] Plot saved: {out_png}")

    # Safe exports
    LCC_clean = sanitize_graph_for_export(LCC)
    try:
        nx.write_graphml(LCC_clean, os.path.join(args.outdir, "LCC.graphml"))
        print("[i] Exported GraphML: LCC.graphml")
    except Exception as e:
        print(f"[w] GraphML export skipped: {e}")
    try:
        nx.write_gexf(LCC_clean, os.path.join(args.outdir, "LCC.gexf"))
        print("[i] Exported GEXF: LCC.gexf")
    except Exception as e:
        print(f"[w] GEXF export skipped: {e}")


if __name__ == "__main__":
    main()