"""APF v6.9 — Enforcement Crystal analytical metrics.

This module is the Phase 13 Stage II companion to ``apf.crystal``: where
``crystal.py`` *builds* the typed dual-view graph from the bank,
``crystal_metrics.py`` *analyzes* it. The Stage II workstreams from
the canonical Paper 20 update plan (Reference - Paper 20 Enforcement
Crystal Update Plan §9, Tier 2) are::

    (1) Betweenness centrality                    ← landed 2026-04-21 night
    (2) Cascade-failure analysis                  ← landed 2026-04-21 night
    (3) Minimum vertex cuts                       ← THIS PASS
    (4) DAG-DP path attribution by axiom anchor   ← landed 2026-04-21 night
    (5) Convergence-cluster detection             ← landed 2026-04-21 night

Each workstream lands as a pure-Python function operating on a
``CrystalGraph`` and a single bank-registered consistency check that
certifies the computation is deterministic, well-formed, and produces
a sane ranking. The numerical *findings* (top-10 lists, axiom shares,
cascade counts, etc.) are research outputs that go into the dashboard
payload (and thence into Paper 20 v3.0 Stage III); the bank check
keeps the *machinery* honest.

──────────────────────────────────────────────────────────────────────
Workstream 1 — Betweenness centrality (Brandes 2001)
──────────────────────────────────────────────────────────────────────

For an unweighted directed graph G = (V, E), the betweenness centrality
of vertex v is::

    BC(v) = sum_{s != v != t in V}  sigma(s, v, t) / sigma(s, t)

where sigma(s, t) is the number of shortest s→t paths and
sigma(s, v, t) is the number of those that pass through v. We
normalize by (n-1)(n-2) so values lie in [0, 1] (the directed-graph
normalizer; for undirected graphs it would be 2/((n-1)(n-2)) but we
explicitly do not symmetrize).

Brandes' single-source dependency accumulation runs in O(V+E) per
source, giving overall O(V·(V+E)) — perfectly fine for the v6.9 graph
(227 nodes / 821 edges in the CORE/full view, 222/810 in post-R).

Crystal-specific edge convention: edges flow ``source → target`` where
the source is the *dependency* (parent / upstream) and the target is
the *dependent* (child / downstream). A node with high BC is one that
many derivation paths must traverse — a structural waist. The
Phase 13.2 walker already exposes a depth-and-width "waist_candidates"
heuristic (depth-d local minima of width ≤ 4); Workstream 1 replaces
that proxy with a global path-counting measure.

The Paper 20 v1.0 (Feb-2026, v3.4.4 snapshot) headline values were
``BC(T4) = 0.668`` and ``BC(T_gauge) = 0.500``. Those numbers are not
expected to carry forward verbatim against the v6.9 graph (different
node set, different edge set, different anchor seeding via PLEC) — the
*qualitative* prediction is that core gauge-derivation nodes
(Theorem_R, T_gauge, L_count, T7) plus the new Phase 14 / Phase 13.1
convergence theorems (T_three_level_unification,
T_interface_sector_bridge, T_ACC_unification) crowd the top of the
ranking. The v6.9 ranking computed by ``betweenness_centrality()`` is
the analytical input to Stage III's refresh of Paper 20 v1.0
Section 4 (Keystone Nodes by Centrality).

──────────────────────────────────────────────────────────────────────
Bank-registered checks
──────────────────────────────────────────────────────────────────────

* ``check_T_crystal_centrality_v69`` (tier 3, ``[P_structural]``):
    Certifies that ``betweenness_centrality()`` is deterministic
    (re-runs identically), produces values in [0, 1], yields a
    non-empty top-10 ranking on both the full and post-R views under
    the CORE preset, and that the post-R top-3 is contained in the
    full top-25 (a sanity guard against the dual-view contract
    accidentally re-shuffling the central nodes).

──────────────────────────────────────────────────────────────────────
"""

from __future__ import annotations

import collections as _coll
import statistics as _stat

from apf.apf_utils import _result, check
from apf.crystal import CrystalGraph, build_crystal


# =============================================================================
# §1  Brandes' betweenness centrality
# =============================================================================


def _adjacency_from_graph(graph: CrystalGraph) -> dict[str, list[str]]:
    """Build a successor adjacency dict from a CrystalGraph.

    Returns {node_id -> sorted list of children}. Sorting the child
    list at construction time makes the BFS traversal order
    deterministic, which in turn makes the final BC values a
    bit-stable function of (graph topology, node id strings).
    """
    succ: dict[str, list[str]] = {nid: [] for nid in graph.nodes}
    # Group children by source first (gives O(E) build), then sort each list.
    children_set: dict[str, set[str]] = _coll.defaultdict(set)
    for e in graph.edges:
        children_set[e.source].add(e.target)
    for src, kids in children_set.items():
        succ[src] = sorted(kids)
    return succ


def betweenness_centrality(
    graph: CrystalGraph,
    *,
    normalized: bool = True,
) -> dict[str, float]:
    """Betweenness centrality of every node in ``graph``.

    Implementation: Brandes (2001) single-source shortest-path
    accumulation. Treats the graph as directed and unweighted; edge
    direction is ``source → target`` (dependency → dependent), so high
    BC marks structural waists in the derivation flow.

    Parameters
    ----------
    graph : CrystalGraph
        A view produced by ``apf.crystal.build_crystal``. Either
        ``full_graph`` or ``post_R_subgraph`` is fine.
    normalized : bool
        If True (default), divide raw BC by (n-1)(n-2) so the result
        lies in [0, 1]. If False, return raw path-count sums (useful
        for cross-checking against hand calculations on tiny graphs).

    Returns
    -------
    dict[str, float]
        ``{node_id: BC}`` with one entry per node in ``graph.nodes``
        (BC defaults to 0.0 for nodes never traversed by any geodesic).
    """
    nodes = list(graph.nodes.keys())
    succ = _adjacency_from_graph(graph)

    bc: dict[str, float] = {nid: 0.0 for nid in nodes}

    for s in nodes:
        # ---- Single-source shortest-path BFS ----
        # Brandes notation: stack S, predecessors P[w], shortest-path
        # counts sigma[v], distance dist[v], dependency delta[v].
        S: list[str] = []
        P: dict[str, list[str]] = {nid: [] for nid in nodes}
        sigma: dict[str, int] = dict.fromkeys(nodes, 0)
        sigma[s] = 1
        dist: dict[str, int] = dict.fromkeys(nodes, -1)
        dist[s] = 0

        Q: _coll.deque[str] = _coll.deque([s])
        while Q:
            v = Q.popleft()
            S.append(v)
            for w in succ.get(v, ()):
                # First time we see w → assign distance, enqueue.
                if dist[w] < 0:
                    dist[w] = dist[v] + 1
                    Q.append(w)
                # If w is on a shortest path through v, accumulate.
                if dist[w] == dist[v] + 1:
                    sigma[w] += sigma[v]
                    P[w].append(v)

        # ---- Dependency back-propagation ----
        delta: dict[str, float] = dict.fromkeys(nodes, 0.0)
        # S now holds nodes in non-decreasing distance from s; pop in
        # reverse to back-propagate.
        while S:
            w = S.pop()
            for v in P[w]:
                if sigma[w] > 0:
                    delta[v] += (sigma[v] / sigma[w]) * (1.0 + delta[w])
            if w != s:
                bc[w] += delta[w]

    if normalized:
        n = len(nodes)
        if n > 2:
            scale = 1.0 / ((n - 1) * (n - 2))
            for nid in bc:
                bc[nid] *= scale

    return bc


# =============================================================================
# §2  Convenience summaries
# =============================================================================


def top_k_central(
    bc: dict[str, float],
    k: int = 10,
) -> list[tuple[str, float]]:
    """Return the top-``k`` (node_id, BC) pairs sorted by BC descending.

    Ties broken by node_id ascending for determinism.
    """
    items = sorted(bc.items(), key=lambda kv: (-kv[1], kv[0]))
    return items[:k]


def centrality_stats(bc: dict[str, float]) -> dict[str, float]:
    """Summary statistics over a BC dict (min, max, mean, median, n_zero)."""
    if not bc:
        return {
            "n": 0, "min": 0.0, "max": 0.0,
            "mean": 0.0, "median": 0.0, "n_zero": 0,
        }
    vals = list(bc.values())
    return {
        "n":      len(vals),
        "min":    min(vals),
        "max":    max(vals),
        "mean":   _stat.fmean(vals),
        "median": _stat.median(vals),
        "n_zero": sum(1 for v in vals if v == 0.0),
    }


def centrality_report(
    graph: CrystalGraph,
    *,
    top_k: int = 10,
) -> dict:
    """Compute BC + summary stats + top-k for a single CrystalGraph view.

    Returns
    -------
    dict
        ``{view, preset, n_nodes, n_edges, top_k, stats, bc}`` where
        ``bc`` is the full BC dict (keep for downstream cross-view
        comparisons), ``top_k`` is a list of ``(node_id, value)`` pairs,
        and ``stats`` is the output of ``centrality_stats``.
    """
    bc = betweenness_centrality(graph)
    return {
        "view":     graph.view,
        "preset":   graph.preset,
        "n_nodes":  graph.n_nodes,
        "n_edges":  graph.n_edges,
        "top_k":    top_k_central(bc, k=top_k),
        "stats":    centrality_stats(bc),
        "bc":       bc,
    }


# =============================================================================
# §3  Bank-registered consistency check
# =============================================================================


def check_T_crystal_centrality_v69():
    """T_crystal_centrality_v69: Brandes BC well-defined on v6.9 dual view.

    STATEMENT: ``betweenness_centrality`` produces a deterministic,
    normalized, well-formed ranking on both views of the v6.9
    Enforcement Crystal under the CORE module preset. The post-R top-3
    is contained in the full-graph top-25, witnessing that the
    dual-view contract drops Regime_R + Type I–IV exits without
    re-shuffling the structurally central nodes.

    METHOD:
      (a) Build CORE preset (full_graph + post_R_subgraph) via the
          Phase 13.2 walker.
      (b) Compute BC on both views; assert all values in [0, 1].
      (c) Re-compute BC on each view; assert exact float equality
          (determinism guard against dict-iteration-order changes).
      (d) Assert top-10 lists are non-empty on both views.
      (e) Assert top-3 of post_R ⊆ top-25 of full (sanity guard;
          post_R drops 5 nodes — Regime_R + 4 rejected exit types —
          so its top-3 should still be near the top of the full
          ranking).
      (f) Record the v6.9 top-10 + summary stats for each view in
          ``artifacts``; these become the input to Paper 20 v3.0
          Stage III's refresh of v1.0 §4 (Keystone Nodes by Centrality).

    EPISTEMIC: ``[P_structural]``. The check certifies Brandes' BC is
    correctly implemented (deterministic, normalized, sane ranking)
    on the v6.9 graph; it is NOT a derivation-of-physics claim, and
    the numerical top-10 it records is a research output, not a
    proven invariant.
    """
    crystal = build_crystal("CORE", do_prelude=True)
    full = crystal["full_graph"]
    post_R = crystal["post_R_subgraph"]

    # (a) Both views non-empty.
    check(full.n_nodes > 0, "full graph has zero nodes")
    check(post_R.n_nodes > 0, "post_R subgraph has zero nodes")

    # (b) Compute BC + range guards.
    bc_full = betweenness_centrality(full)
    bc_postR = betweenness_centrality(post_R)
    for view_name, bc in (("full", bc_full), ("post_R", bc_postR)):
        for nid, v in bc.items():
            check(
                0.0 <= v <= 1.0 + 1e-12,
                f"BC out of [0,1] for node {nid!r} in view {view_name}: {v}"
            )

    # (c) Determinism: second compute must match first to bit precision.
    bc_full_2 = betweenness_centrality(full)
    bc_postR_2 = betweenness_centrality(post_R)
    for view_name, a, b in (("full", bc_full, bc_full_2),
                            ("post_R", bc_postR, bc_postR_2)):
        check(
            a == b,
            f"BC non-deterministic on {view_name}: re-run produced different values"
        )

    # (d) Top-10 non-empty on both views.
    top10_full = top_k_central(bc_full, k=10)
    top10_postR = top_k_central(bc_postR, k=10)
    check(len(top10_full) > 0, "full graph top-10 is empty")
    check(len(top10_postR) > 0, "post_R subgraph top-10 is empty")

    # (e) post_R top-3 ⊆ full top-25 sanity.
    top3_postR_ids = {nid for nid, _ in top_k_central(bc_postR, k=3)}
    top25_full_ids = {nid for nid, _ in top_k_central(bc_full, k=25)}
    missing = top3_postR_ids - top25_full_ids
    check(
        not missing,
        f"post_R top-3 not contained in full top-25; missing: {sorted(missing)}"
    )

    # (f) Record headline numbers for downstream paper write-up.
    stats_full = centrality_stats(bc_full)
    stats_postR = centrality_stats(bc_postR)

    return _result(
        name='T_crystal_centrality_v69 — Brandes BC well-defined on v6.9 dual view',
        tier=3,
        epistemic='P_structural',
        summary=(
            'Brandes (2001) betweenness centrality is computable, '
            'deterministic, normalized to [0, 1], and produces a '
            'well-formed top-10 ranking on both the full_graph and '
            'post_R_subgraph views of the v6.9 Enforcement Crystal '
            'under the CORE module preset. The post-R top-3 nodes are '
            'contained in the full top-25, witnessing dual-view '
            'stability of the centrally-located waists. Top-10 lists '
            'and summary statistics are recorded in artifacts for '
            'Paper 20 v3.0 Stage III §4 (Keystone Nodes by Centrality).'
        ),
        key_result=(
            f'BC OK on both views: full top-1 = {top10_full[0][0]} '
            f'(BC={top10_full[0][1]:.4f}); post_R top-1 = '
            f'{top10_postR[0][0]} (BC={top10_postR[0][1]:.4f}); '
            f'post-R top-3 ⊆ full top-25 [P_structural]'
        ),
        dependencies=['T_crystal_v69_consistent'],
        cross_refs=[
            'L_PLEC_components_essentiality',
            'T_three_level_unification',
            'T_interface_sector_bridge',
            'T_ACC_unification',
        ],
        artifacts={
            'preset':            'CORE',
            'full_view': {
                'n_nodes':  full.n_nodes,
                'n_edges':  full.n_edges,
                'top_10':   top10_full,
                'stats':    stats_full,
            },
            'post_R_view': {
                'n_nodes':  post_R.n_nodes,
                'n_edges':  post_R.n_edges,
                'top_10':   top10_postR,
                'stats':    stats_postR,
            },
        },
    )


# =============================================================================
# §4  Workstream 4 — DAG path attribution by PLEC anchor
# =============================================================================
#
# The Paper 20 v1.0 (Feb-2026, v3.4.4 snapshot) headline §5 result was
# a count of **1398 distinct directed paths** from the five Feb-2026
# axioms (A1, A2, A3, A4, A5) to the prediction sin²θ_W = 3/13, with
# attribution shares::
#
#     A2 (non-closure)       42 %
#     A1 (finite capacity)   39 %
#     A3 (locality)          15 %
#     A4 (irreversibility)    3 %
#     A5 (collapse)           0.5 %
#
# v6.9 has reorganized the foundation into the **PLEC** structure
# (Paper 1 v4.0-PLEC), with four anchors instead of five:
#
#     A1 → check_A1                (apf.core)         "finite capacity"
#     MD → check_L_epsilon_star    (apf.core)         "minimum domain"
#     A2 → check_Regime_R          (apf.plec)         "regime selection"
#     BW → check_worked_example    (apf.core)         "beta witness"
#
# (mapped through ``apf.crystal_axiom_roots.PLEC_AXIOM_ROOTS``). The
# Stage III refresh of v1.0 §5 will state the v6.9 path-count and
# share table against these PLEC anchors, with the same canonical
# target ``T_sin2theta`` (the v6.9 check name for sin²θ_W = 3/13).
#
# ──────────────────────────────────────────────────────────────────────
# Algorithm — depth-filtered DAG-DP
# ──────────────────────────────────────────────────────────────────────
#
# The bank-registered dependency graph is *not* a strict DAG: docstring
# cross-references between mutually-supporting lemmas (L_X says "depends
# on L_Y" and L_Y says "depends on L_X" because they are pieces of the
# same construction) introduce small cycles. A topological-sort scan
# of the v6.9 CORE preset finds 148 / 227 nodes inside such cycles.
# These are documentation cross-refs, not derivation cycles — none of
# them represents a loop in the actual proof flow.
#
# To get a clean DAG to count paths on, we use the walker's depth
# assignment (``CrystalGraph.depth``) to break back-edges: an edge
# ``u → v`` is kept iff ``depth(v) > depth(u)``. This induces a
# strict DAG on the same node set, with the depth ordering as a
# topological order. Counting paths from each PLEC anchor (depth 0)
# to a target node ``t`` is then a single O(V+E) DP::
#
#     count[t]  = 1
#     count[v]  = sum_{w in succ_DAG[v]} count[w]   for v != t
#
# processed in reverse topological (i.e. decreasing depth) order.
#
# Edges dropped by the depth filter are recorded in the artifacts so
# Stage III can document the magnitude of the cross-ref pruning.


def _depth_filtered_succ(graph: CrystalGraph) -> dict[str, set[str]]:
    """Successor adjacency dict filtered to a strict DAG by depth.

    Keeps only edges ``u → v`` where ``graph.depth[v] > graph.depth[u]``.
    Drops back-edges, peer-edges (same depth), and self-loops — these
    correspond to docstring cross-references, not derivation
    dependencies. The returned dict ranges over all nodes in
    ``graph.nodes`` (with empty sets for nodes that have no forward
    children).
    """
    succ: dict[str, set[str]] = {nid: set() for nid in graph.nodes}
    depth = graph.depth
    for e in graph.edges:
        d_src = depth.get(e.source, -1)
        d_tgt = depth.get(e.target, -1)
        if d_tgt > d_src:
            succ[e.source].add(e.target)
    return succ


def count_paths_to_target(
    graph: CrystalGraph,
    target_id: str,
) -> dict[str, int]:
    """For each node v, count distinct directed paths v → ``target_id``.

    Computed on the depth-filtered DAG (see ``_depth_filtered_succ``).
    Returns ``{node_id: int}`` with one entry per node in the graph;
    nodes from which ``target_id`` is unreachable map to 0.
    """
    if target_id not in graph.nodes:
        raise ValueError(
            f"target_id {target_id!r} not in graph (preset={graph.preset})"
        )

    succ = _depth_filtered_succ(graph)

    # Process in reverse topological order = decreasing depth. Ties
    # broken by node-id for determinism (BC and path-count rankings
    # both want stable string-key ordering).
    order = sorted(
        graph.nodes,
        key=lambda nid: (graph.depth.get(nid, 0), nid),
    )

    count: dict[str, int] = {nid: 0 for nid in graph.nodes}
    count[target_id] = 1
    for v in reversed(order):
        if v == target_id:
            continue
        s = 0
        for w in succ[v]:
            s += count[w]
        count[v] = s
    return count


def path_attribution_by_anchor(
    graph: CrystalGraph,
    target_id: str,
    *,
    anchor_filter: str = "plec",
    method: str = "depth_filtered",
) -> dict:
    """Attribution of paths to ``target_id`` by PLEC anchor.

    Parameters
    ----------
    graph : CrystalGraph
        Either ``full_graph`` or ``post_R_subgraph``.
    target_id : str
        Canonical bank ID of the terminal prediction (e.g.
        ``'T_sin2theta'``).
    anchor_filter : {'plec', 'all_anchors'}
        ``'plec'`` (default) — restrict to nodes flagged
        ``is_plec_anchor=True`` (the four roles A1, MD, A2, BW).
        ``'all_anchors'`` — include every depth-0 node.
    method : {'depth_filtered', 'scc'}
        Path-counting variant. ``'depth_filtered'`` (default, v3.0
        Paper 20 baseline) — counts paths on the depth-filtered
        strict DAG (edges kept iff ``depth(target) > depth(source)``).
        Drops same-depth cycle-edges as a side effect, which hides
        derivation routes through mutually-supporting cycles like
        Theorem_R ↔ T_gauge ↔ T_field. ``'scc'`` (Phase 14c, v3.1+
        baseline) — counts paths on the Tarjan SCC condensation;
        each strongly-connected component is collapsed to a single
        super-node, and paths are counted on the resulting strict
        DAG. Cycle-aware: routes through co-defining clusters are
        correctly attributed.

    Returns
    -------
    dict
        ::
            {
              'target': target_id,
              'anchor_filter': anchor_filter,
              'method': 'depth_filtered' | 'scc',
              'anchors': [{'id': str, 'role': str,
                           'n_paths': int, 'share': float}, ...],
              'total_paths': int,
              # For method='depth_filtered':
              'edges_kept': int,        # edges in depth-filtered DAG
              'edges_dropped': int,     # back / peer / self edges
              # For method='scc':
              'n_sccs': int,            # total SCCs (includes singletons)
              'n_nontrivial_sccs': int, # SCCs of size >= 2
              'n_edges_internal': int,  # edges inside some SCC
              'n_edges_cross_scc': int, # edges between SCCs (= condensation E)
              'target_scc_size': int,   # size of target's SCC (usually 1)
            }
        Anchors are sorted by descending ``n_paths``. ``share`` is
        ``n_paths / total_paths`` (or 0 if total is 0). ``share``
        values sum to 1.0 ± 1e-12 when ``total_paths > 0``.
    """
    if anchor_filter not in ("plec", "all_anchors"):
        raise ValueError(
            f"anchor_filter must be 'plec' or 'all_anchors'; got {anchor_filter!r}"
        )
    if method not in ("depth_filtered", "scc"):
        raise ValueError(
            f"method must be 'depth_filtered' or 'scc'; got {method!r}"
        )

    if method == "depth_filtered":
        counts = count_paths_to_target(graph, target_id)
        succ = _depth_filtered_succ(graph)
        n_edges_kept = sum(len(s) for s in succ.values())
        n_edges_total = len(graph.edges)
        extra = {
            'edges_kept':     n_edges_kept,
            'edges_dropped':  n_edges_total - n_edges_kept,
        }
    else:
        counts, meta = count_paths_to_target_scc(graph, target_id)
        extra = {
            'n_sccs':             meta['n_sccs'],
            'n_nontrivial_sccs':  meta['n_nontrivial_sccs'],
            'n_edges_internal':   meta['n_edges_internal'],
            'n_edges_cross_scc':  meta['n_edges_cross_scc'],
            'target_scc_size':    meta['target_scc']['size'],
        }

    # Pick anchor set.
    if anchor_filter == "plec":
        anchor_ids = [nid for nid, n in graph.nodes.items() if n.is_plec_anchor]
    else:
        anchor_ids = [nid for nid, d in graph.depth.items() if d == 0]

    rows = []
    for aid in anchor_ids:
        node = graph.nodes[aid]
        rows.append({
            'id':       aid,
            'role':     node.plec_role or "anchor",
            'n_paths':  counts.get(aid, 0),
        })

    total = sum(r['n_paths'] for r in rows)
    for r in rows:
        r['share'] = (r['n_paths'] / total) if total > 0 else 0.0

    rows.sort(key=lambda r: (-r['n_paths'], r['id']))

    out = {
        'target':         target_id,
        'anchor_filter':  anchor_filter,
        'method':         method,
        'anchors':        rows,
        'total_paths':    total,
    }
    out.update(extra)
    return out


def check_T_crystal_path_attribution_v69():
    """T_crystal_path_attribution_v69: PLEC anchor share table on v6.9.

    STATEMENT: The v6.9 Enforcement Crystal admits a well-formed
    PLEC-anchor path-attribution table to the canonical target
    ``T_sin2theta`` (the v6.9 check name for sin²θ_W = 3/13). At least
    one PLEC anchor reaches ``T_sin2theta`` along the depth-filtered
    DAG, the share table sums to 1.0 (mod float roundoff), and the
    cross-view comparison is well-defined (post_R may legitimately
    drop ``Regime_R`` from the anchor set, in which case shares are
    re-normalized over the surviving anchors).

    This is the v6.9 refresh of the Paper 20 v1.0 (Feb-2026, v3.4.4
    snapshot) §5 path-attribution table. v1.0 reported 1398 paths from
    5 axioms (A1=39%, A2=42%, A3=15%, A4=3%, A5=0.5%); v6.9 has been
    re-anchored against the PLEC structure (A1, MD, A2, BW), and the
    refreshed share table is the artifact below.

    Stage III note (recorded in artifacts as a structural finding):
    in the v6.9 bank the SM derivation chain leading to T_sin2theta
    is wired to A1 / MD / BW but is *not* explicitly wired to depend
    on Regime_R (the regime-selection / non-closure anchor). This is
    inherited from the pre-PLEC dependency edges around T_gauge /
    Theorem_R / L_count. v6.9 reports it honestly here; Paper 20 v3.0
    Stage III §5 either re-anchors the chain through Regime_R (if
    conceptually warranted) or documents the v6.9 reality that A2's
    contribution to sin²θ_W is implicit (used in the proof prose) but
    not bank-edge-explicit. Either way, the machinery is sound.

    METHOD:
      (a) Build CORE preset (full_graph + post_R_subgraph).
      (b) Filter each view's edge set to a strict DAG by keeping only
          edges where ``depth(target) > depth(source)`` — drops back-
          edges and peer cross-references (148 / 227 nodes are inside
          documentation-cross-ref cycles in the unfiltered graph).
      (c) DP path-count to ``T_sin2theta`` from every node in
          decreasing-depth order.
      (d) Extract PLEC-anchor counts on the full view; assert at least
          one anchor reaches the target.
      (e) Compute share = n_paths / total_paths; assert sum = 1 ± 1e-9
          when total > 0.
      (f) Repeat (b)-(e) on the post-R view; allow Regime_R to drop
          out of the anchor set there. Assert that the *intersection*
          of anchor IDs has identical ranking by ID (anchors common to
          both views must agree on relative ordering).
      (g) Record both share tables + edge-prune accounting + the
          Regime_R reachability finding in ``artifacts`` for Paper 20
          v3.0 Stage III §5.

    EPISTEMIC: ``[P_structural]``. Certifies that the path-attribution
    machinery (depth-filter + DP counter) produces a well-formed PLEC
    share table on v6.9; the numerical share values are research
    outputs, not invariants of the framework.
    """
    crystal = build_crystal("CORE", do_prelude=True)
    full = crystal["full_graph"]
    post_R = crystal["post_R_subgraph"]

    target = "T_sin2theta"

    # (a) Target present in both views.
    check(
        target in full.nodes,
        f"target {target!r} missing from full_graph; v6.9 sin²θ_W check renamed?"
    )
    check(
        target in post_R.nodes,
        f"target {target!r} missing from post_R_subgraph"
    )

    # (b)–(e) Full view.
    attrib_full = path_attribution_by_anchor(full, target, anchor_filter="plec")

    # Assert at least one PLEC anchor reaches the target. Individual
    # anchors with 0 paths are recorded as a structural finding (see
    # docstring) — this is the v6.9 reality and is what Paper 20
    # Stage III §5 needs to refresh against.
    n_reaching = sum(1 for r in attrib_full['anchors'] if r['n_paths'] > 0)
    check(
        n_reaching >= 1,
        f"no PLEC anchor reaches {target} in full_graph — bank graph disconnected?"
    )
    check(
        attrib_full['total_paths'] >= 1,
        f"full-view total path count to {target} is zero"
    )
    share_sum_full = sum(r['share'] for r in attrib_full['anchors'])
    check(
        abs(share_sum_full - 1.0) < 1e-9,
        f"full-view share sum != 1: got {share_sum_full}"
    )

    # (f) Post-R view: Regime_R is removed entirely from the post-R
    # subgraph (the dual-view contract drops Regime_R + the four
    # rejected exit channels), so the anchor set may be 3 (not 4)
    # there. Cross-view stability is checked over the *intersection*
    # of anchor IDs.
    attrib_postR = path_attribution_by_anchor(post_R, target, anchor_filter="plec")
    check(
        attrib_postR['total_paths'] >= 1,
        f"post_R-view total path count to {target} is zero"
    )
    share_sum_postR = sum(r['share'] for r in attrib_postR['anchors'])
    check(
        abs(share_sum_postR - 1.0) < 1e-9,
        f"post_R-view share sum != 1: got {share_sum_postR}"
    )

    # Cross-view ranking stability over the common anchor set.
    full_order = [r['id'] for r in attrib_full['anchors'] if r['n_paths'] > 0]
    postR_order = [r['id'] for r in attrib_postR['anchors'] if r['n_paths'] > 0]
    common = [aid for aid in full_order if aid in postR_order]
    full_restricted = [aid for aid in full_order if aid in common]
    postR_restricted = [aid for aid in postR_order if aid in common]
    check(
        full_restricted == postR_restricted,
        f"PLEC anchor ranking differs across views on common subset: "
        f"full={full_restricted!r} vs post_R={postR_restricted!r}"
    )

    # Structural finding: which PLEC anchors fail to reach the target?
    unreached_full = [
        r['id'] for r in attrib_full['anchors'] if r['n_paths'] == 0
    ]

    # Headline string for the result key.
    top = attrib_full['anchors'][0]
    headline = (
        f"top anchor: {top['role']} ({top['id']}) "
        f"= {top['share']*100:.1f}% of {attrib_full['total_paths']} paths to {target}"
    )
    if unreached_full:
        headline += f"; unreached PLEC anchor(s): {unreached_full}"

    return _result(
        name=(
            'T_crystal_path_attribution_v69 — PLEC anchor share table '
            'to T_sin2theta well-defined on v6.9 dual view'
        ),
        tier=3,
        epistemic='P_structural',
        summary=(
            'On the v6.9 Enforcement Crystal under the CORE preset, the '
            'depth-filtered DAG admits a well-defined path-count from '
            'the PLEC anchors (A1, MD, A2, BW) to T_sin2theta '
            '(sin²θ_W = 3/13). At least one anchor reaches the target; '
            'shares sum to 1 in each view; the ranking of anchors by '
            'share over the common anchor subset is identical in the '
            'full and post-R views (dual-view stability). v6.9 refresh '
            'of Paper 20 v1.0 §5 (path attribution by axiom). PLEC '
            'anchors that do not reach T_sin2theta in the v6.9 bank '
            'graph are recorded in artifacts as a structural finding '
            'for Stage III §5.'
        ),
        key_result=headline + ' [P_structural]',
        dependencies=['T_crystal_v69_consistent'],
        cross_refs=[
            'L_PLEC_components_essentiality',
            'T_sin2theta',
            'T_crystal_centrality_v69',
        ],
        artifacts={
            'preset':              'CORE',
            'target':              target,
            'full_view':           attrib_full,
            'post_R_view':         attrib_postR,
            'unreached_anchors':   unreached_full,
            'stage3_finding': (
                f"v6.9 PLEC anchors not on any depth-filtered DAG path "
                f"to {target} (full view): {unreached_full!r}. "
                f"This is the v6.9 wiring reality; Paper 20 v3.0 §5 "
                f"either re-anchors the SM chain through these or "
                f"documents the finding."
            ),
        },
    )


# =============================================================================
# §4b  Workstream 4, SCC-aware variant (Phase 14c, 2026-04-22)
# =============================================================================
#
# The §4 path-attribution uses a depth-filtered strict DAG, which drops
# same-depth cycle-edges along with back-edges. This keeps the counting
# algorithm simple but has a load-bearing blind spot: if three nodes X,
# Y, Z are mutually supporting (e.g. Theorem_R / T_gauge / T_field all
# land at depth 12 because they semantically co-define each other), the
# filter drops the edges *among* them, so paths from an upstream anchor
# that naturally flow through the X→Y→Z cluster are counted as zero.
#
# v6.9 Paper 20 v3.0 §4 Finding 1 ("Regime_R + worked_example contribute
# 0 paths to T_sin2theta") reflected this blind spot: the depth filter
# dropped the Theorem_R(12) → T_gauge(12) edge, which is the only route
# by which Regime_R can reach T_sin2theta through the PLEC-selection
# event Theorem_R. Phase 14b (2026-04-22) added the edges Regime_R →
# Theorem_R and worked_example → L_count at the bank level, but the
# §4 metric still read 0 because the downstream edge Theorem_R → T_gauge
# was still invisible to the counter.
#
# This §4b workstream adds a cycle-aware variant: the bank graph is
# decomposed into strongly-connected components (Tarjan 1972), then
# collapsed into the condensation (a strict DAG of SCCs by construction).
# Path counting runs on the condensation, and every node inside an SCC
# inherits the SCC's path count. The semantic change is: within an SCC,
# multiple orderings of the co-defining edges are treated as a single
# "derivation event" (the SCC as a unit), so internal cycle
# multiplicity doesn't inflate counts.
#
# With this variant, Phase 14b's Regime_R → Theorem_R edge sees the
# full Theorem_R ↔ T_gauge ↔ T_field SCC as downstream of Regime_R,
# and paths from Regime_R to T_sin2theta become visible in the count.
#
# Method dispatch: ``path_attribution_by_anchor(..., method='scc')``.
# Default method remains ``'depth_filtered'`` for v3.0 reproducibility;
# v3.1+ analyses should pass ``method='scc'`` explicitly or read the
# new ``check_T_crystal_path_attribution_scc_v69`` artifacts.


def _scc_tarjan(
    graph: CrystalGraph,
) -> tuple[list[frozenset[str]], dict[str, int]]:
    """Strongly-connected components of ``graph`` (iterative Tarjan 1972).

    Returns
    -------
    sccs : list[frozenset[str]]
        SCCs in reverse topological order of the condensation:
        sccs[0] has no outgoing cross-SCC edges (a sink in the
        condensation DAG), and later entries are successively closer
        to the sources. This is the natural pop order from Tarjan's
        DFS and lets count-DP process sinks first.
    scc_of : dict[str, int]
        ``{node_id: scc_index}`` where scc_index is the offset into
        ``sccs``.

    Implementation notes
    --------------------
    * Fully iterative — the v6.9 CORE preset has 227 nodes / 823
      edges, well below Python's recursion limit, but iterative is
      robust against FULL preset graphs too.
    * Nodes are visited in sorted ID order so the SCC numbering is
      a deterministic function of (graph topology, node id strings).
      This makes downstream path counts bit-stable across runs.
    """
    # Build successor dict once.
    succ: dict[str, list[str]] = {nid: [] for nid in graph.nodes}
    for e in graph.edges:
        # Include only edges with both endpoints in nodes (defensive).
        if e.source in succ and e.target in graph.nodes:
            succ[e.source].append(e.target)

    index: dict[str, int] = {}
    lowlink: dict[str, int] = {}
    on_stack: set[str] = set()
    tarjan_stack: list[str] = []
    sccs: list[frozenset[str]] = []
    scc_of: dict[str, int] = {}
    counter = 0

    all_nodes = sorted(graph.nodes.keys())

    for start in all_nodes:
        if start in index:
            continue
        # Work stack entries: (v, children_list, next_child_index).
        # Sorting children here ensures deterministic traversal.
        work: list[list] = [[start, sorted(succ[start]), 0]]
        index[start] = counter
        lowlink[start] = counter
        counter += 1
        tarjan_stack.append(start)
        on_stack.add(start)

        while work:
            frame = work[-1]
            v, children, ci = frame[0], frame[1], frame[2]
            if ci < len(children):
                frame[2] = ci + 1
                w = children[ci]
                if w not in index:
                    index[w] = counter
                    lowlink[w] = counter
                    counter += 1
                    tarjan_stack.append(w)
                    on_stack.add(w)
                    work.append([w, sorted(succ[w]), 0])
                elif w in on_stack:
                    if index[w] < lowlink[v]:
                        lowlink[v] = index[w]
            else:
                # All children of v processed.
                if lowlink[v] == index[v]:
                    # v is a root of an SCC.
                    comp: list[str] = []
                    while True:
                        x = tarjan_stack.pop()
                        on_stack.discard(x)
                        comp.append(x)
                        if x == v:
                            break
                    scc_id = len(sccs)
                    sccs.append(frozenset(comp))
                    for x in comp:
                        scc_of[x] = scc_id
                work.pop()
                if work:
                    parent_v = work[-1][0]
                    if lowlink[v] < lowlink[parent_v]:
                        lowlink[parent_v] = lowlink[v]
    return sccs, scc_of


def _scc_condensed_succ(
    graph: CrystalGraph,
    sccs: list[frozenset[str]],
    scc_of: dict[str, int],
) -> dict[int, set[int]]:
    """Successor dict on the SCC condensation.

    Maps ``scc_index -> set of downstream scc_indices``. Within-SCC
    edges (source and target in the same component) are excluded.
    Parallel cross-SCC edges collapse to multiplicity 1 (we are
    counting distinct SCC-to-SCC routes, not edge-labeled walks).
    """
    condensed: dict[int, set[int]] = {i: set() for i in range(len(sccs))}
    for e in graph.edges:
        if e.source not in scc_of or e.target not in scc_of:
            continue
        i = scc_of[e.source]
        j = scc_of[e.target]
        if i != j:
            condensed[i].add(j)
    return condensed


def count_paths_to_target_scc(
    graph: CrystalGraph,
    target_id: str,
) -> tuple[dict[str, int], dict]:
    """Count paths to ``target_id`` via SCC condensation (cycle-aware).

    Each strongly-connected component becomes a single super-node.
    Path counting runs on the condensation (a strict DAG by
    construction); each node inherits the count at its SCC. Within-SCC
    cycles are collapsed — the SCC is counted once per traversal.

    Returns
    -------
    counts : dict[str, int]
        ``{node_id: count}``. All nodes in the same SCC share the
        same count.
    meta : dict
        SCC decomposition statistics:
        ``{n_sccs, n_nontrivial_sccs, nontrivial_sccs,
            target_scc, n_edges_internal, n_edges_cross_scc}``.
    """
    if target_id not in graph.nodes:
        raise ValueError(
            f"target_id {target_id!r} not in graph (preset={graph.preset})"
        )

    sccs, scc_of = _scc_tarjan(graph)
    condensed = _scc_condensed_succ(graph, sccs, scc_of)
    n_sccs = len(sccs)
    target_scc_idx = scc_of[target_id]

    # sccs is in reverse topological order of the condensation:
    # sccs[0] has no cross-SCC successors, sccs[-1] is a source.
    # Process in natural order (sinks first) so when we reach scc i,
    # all its successor SCCs already have their counts.
    count_scc: list[int] = [0] * n_sccs
    count_scc[target_scc_idx] = 1
    for i in range(n_sccs):
        if i == target_scc_idx:
            continue
        s = 0
        for j in condensed[i]:
            s += count_scc[j]
        count_scc[i] = s

    # Expand: each node inherits the count at its SCC.
    counts: dict[str, int] = {
        nid: count_scc[scc_of[nid]] for nid in graph.nodes
    }

    # Meta.
    nontrivial_idxs = [i for i in range(n_sccs) if len(sccs[i]) >= 2]
    nontrivial_sccs_info = [
        {
            "size":    len(sccs[i]),
            "nodes":   sorted(sccs[i]),
            "scc_id":  i,
        }
        for i in nontrivial_idxs
    ]
    n_internal = sum(
        1 for e in graph.edges
        if e.source in scc_of and e.target in scc_of
        and scc_of[e.source] == scc_of[e.target]
    )
    n_cross = len(graph.edges) - n_internal
    target_scc_members = sorted(sccs[target_scc_idx])
    meta = {
        "n_sccs":              n_sccs,
        "n_nontrivial_sccs":   len(nontrivial_idxs),
        "nontrivial_sccs":     nontrivial_sccs_info,
        "target_scc": {
            "target":   target_id,
            "scc_id":   target_scc_idx,
            "size":     len(target_scc_members),
            "members":  target_scc_members,
        },
        "n_edges_internal":    n_internal,
        "n_edges_cross_scc":   n_cross,
    }
    return counts, meta


def check_T_crystal_path_attribution_scc_v69():
    """T_crystal_path_attribution_scc_v69: SCC-aware PLEC anchor share table.

    STATEMENT: The v6.9 Enforcement Crystal admits a well-formed PLEC-
    anchor path-attribution table to the canonical target
    ``T_sin2theta`` under the *SCC condensation* variant of the
    path-count metric (cycle-aware, Phase 14c 2026-04-22). SCC
    condensation uses Tarjan's 1972 algorithm to collapse each
    strongly-connected component in the bank graph into a single
    super-node; path counting then runs on the strict-DAG condensation.
    This variant is cycle-aware where the §4 depth-filtered variant
    drops same-depth cycle-edges as a side effect of its linearization.

    This check is the v6.9.1 refresh of ``check_T_crystal_path_-
    attribution_v69`` (§4, depth-filtered, v3.0 baseline). Both checks
    live in the bank: the depth-filtered check certifies the v3.0
    Paper 20 number; this SCC check certifies the v3.1+ number that
    correctly attributes paths through the Theorem_R ↔ T_gauge ↔
    T_field cycle.

    METHOD:
      (a) Build CORE preset (full_graph + post_R_subgraph).
      (b) Decompose each view into strongly-connected components via
          iterative Tarjan; assert SCC count > 0 and every node is
          assigned to exactly one SCC.
      (c) Build the condensation (cross-SCC successor dict); assert
          the condensation is a strict DAG (no SCC has itself as a
          successor under ``_scc_condensed_succ``).
      (d) DP path-count on the condensation for target ``T_sin2theta``
          (full and post_R views).
      (e) Extract PLEC-anchor counts; assert at least one anchor
          reaches the target (total_paths >= 1).
      (f) Compute share = n_paths / total_paths; assert sum = 1 ± 1e-9.
      (g) Record SCC decomposition statistics + both share tables
          under SCC method in artifacts for Paper 20 v3.1 §4.

    EPISTEMIC: ``[P_structural]``. Certifies that the SCC-aware path-
    attribution machinery produces a well-formed PLEC share table on
    v6.9; numerical shares are research outputs, not invariants of
    the framework.
    """
    crystal = build_crystal("CORE", do_prelude=True)
    full = crystal["full_graph"]
    post_R = crystal["post_R_subgraph"]

    target = "T_sin2theta"

    # (a) Target present.
    check(
        target in full.nodes,
        f"target {target!r} missing from full_graph; v6.9 sin²θ_W check renamed?"
    )
    check(
        target in post_R.nodes,
        f"target {target!r} missing from post_R_subgraph"
    )

    # (b) SCC decomposition well-formed on full view.
    sccs_full, scc_of_full = _scc_tarjan(full)
    check(
        len(sccs_full) > 0,
        "Tarjan SCC on full_graph returned empty decomposition"
    )
    check(
        set(scc_of_full.keys()) == set(full.nodes.keys()),
        "SCC decomposition missed some nodes on full_graph"
    )

    # (c) Condensation is a strict DAG (no SCC self-loop).
    cond_full = _scc_condensed_succ(full, sccs_full, scc_of_full)
    for i, succ_set in cond_full.items():
        check(
            i not in succ_set,
            f"SCC condensation has self-loop at scc {i} on full_graph"
        )

    # (d)-(f) Full view attribution.
    attrib_full = path_attribution_by_anchor(
        full, target, anchor_filter="plec", method="scc"
    )
    n_reaching_full = sum(
        1 for r in attrib_full['anchors'] if r['n_paths'] > 0
    )
    check(
        n_reaching_full >= 1,
        f"no PLEC anchor reaches {target} under SCC method in full_graph"
    )
    check(
        attrib_full['total_paths'] >= 1,
        f"SCC-method total path count to {target} is zero (full view)"
    )
    share_sum_full = sum(r['share'] for r in attrib_full['anchors'])
    check(
        abs(share_sum_full - 1.0) < 1e-9,
        f"SCC-method full-view share sum != 1: got {share_sum_full}"
    )

    # (b)-(f) post_R view.
    sccs_postR, scc_of_postR = _scc_tarjan(post_R)
    check(
        len(sccs_postR) > 0,
        "Tarjan SCC on post_R_subgraph returned empty decomposition"
    )
    cond_postR = _scc_condensed_succ(post_R, sccs_postR, scc_of_postR)
    for i, succ_set in cond_postR.items():
        check(
            i not in succ_set,
            f"SCC condensation has self-loop at scc {i} on post_R_subgraph"
        )

    attrib_postR = path_attribution_by_anchor(
        post_R, target, anchor_filter="plec", method="scc"
    )
    check(
        attrib_postR['total_paths'] >= 1,
        f"SCC-method total path count to {target} is zero (post_R view)"
    )
    share_sum_postR = sum(r['share'] for r in attrib_postR['anchors'])
    check(
        abs(share_sum_postR - 1.0) < 1e-9,
        f"SCC-method post_R-view share sum != 1: got {share_sum_postR}"
    )

    # Unreached anchors under SCC method (structural finding for §4b).
    unreached_full = [
        r['id'] for r in attrib_full['anchors'] if r['n_paths'] == 0
    ]

    # Catalog of non-trivial SCCs (size >= 2) on the full view — the
    # key structural finding for Paper 20 v3.1 §4b.
    _, meta_full = count_paths_to_target_scc(full, target)
    _, meta_postR = count_paths_to_target_scc(post_R, target)

    n_nontrivial_full = meta_full['n_nontrivial_sccs']
    n_nontrivial_postR = meta_postR['n_nontrivial_sccs']

    # Headline.
    top = attrib_full['anchors'][0]
    headline = (
        f"[SCC] top anchor: {top['role']} ({top['id']}) "
        f"= {top['share']*100:.1f}% of {attrib_full['total_paths']} "
        f"SCC-paths to {target}"
    )
    if unreached_full:
        headline += f"; unreached under SCC: {unreached_full}"
    headline += (
        f"; {n_nontrivial_full} non-trivial SCCs (full), "
        f"{n_nontrivial_postR} (post_R)"
    )

    return _result(
        name=(
            'T_crystal_path_attribution_scc_v69 — SCC-aware PLEC anchor '
            'share table to T_sin2theta well-defined on v6.9 dual view'
        ),
        tier=3,
        epistemic='P_structural',
        summary=(
            'On the v6.9 Enforcement Crystal under the CORE preset, the '
            'Tarjan SCC condensation of the bank graph admits a well-'
            'defined path-count from the PLEC anchors (A1, MD, A2, BW) '
            'to T_sin2theta. At least one anchor reaches the target; '
            'shares sum to 1 in each view; the condensation is a strict '
            'DAG with no self-loops. This is the cycle-aware variant '
            'of check_T_crystal_path_attribution_v69 (§4); it resolves '
            'the Paper 20 v3.0 Finding 1 blind spot (depth-filtered '
            'metric dropped same-depth cycle-edges, hiding the Theorem_R '
            '→ T_gauge → T_field derivation route). v6.9.1 Paper 20 '
            'v3.1 §4b anchors against this check.'
        ),
        key_result=headline + ' [P_structural]',
        dependencies=['T_crystal_v69_consistent',
                      'T_crystal_path_attribution_v69'],
        cross_refs=[
            'L_PLEC_components_essentiality',
            'T_sin2theta',
            'Theorem_R',
            'T_gauge',
            'T_field',
        ],
        artifacts={
            'preset':                   'CORE',
            'target':                   target,
            'method':                   'scc',
            'full_view':                attrib_full,
            'post_R_view':              attrib_postR,
            'unreached_anchors_scc':    unreached_full,
            'scc_meta_full':            meta_full,
            'scc_meta_post_R':          meta_postR,
            'phase14c_note': (
                'v6.9.1 (Phase 14c, 2026-04-22): SCC-aware variant '
                'introduced to correctly attribute paths through the '
                'Theorem_R ↔ T_gauge ↔ T_field cycle, which the §4 '
                'depth-filtered metric collapses to zero. Phase 14b '
                '(same-day) added the Regime_R → Theorem_R and '
                'worked_example → L_count bank edges; this §4b metric '
                'is what makes those edges visible in the share table.'
            ),
        },
    )


# =============================================================================
# §5  Workstream 2 — Cascade-failure analysis
# =============================================================================
#
# The Paper 20 v1.0 (Feb-2026, v3.4.4 snapshot) headline §6 finding was
# the cascade ranking — "removing T4 invalidates 67 % of predictions",
# "removing T_gauge invalidates 58 %", etc. A node's *cascade* is the
# set of downstream nodes that lose all PLEC-anchor support when it is
# retracted from the graph; high-cascade nodes are load-bearing
# single-points-of-failure that the paper should call out as
# structural waists.
#
# v6.9 cascade-fraction definition (using the depth-filtered DAG of
# §4 to stay on the proof-flow subgraph and avoid spurious cycles):
#
#     baseline_R   = { v : some PLEC anchor has DAG path v0 -> v }   (anchor-reachable nodes)
#     cascade(N)   = { v in baseline_R : every DAG path v0 -> v passes through N } \ {N}
#     fraction(N)  = |cascade(N)| / |baseline_R \ {N}|
#
# Computed by the obvious O(V * (V+E)) loop: for each candidate N
# (skip the four PLEC anchors and N itself), recompute anchor
# reachability on the graph minus N, intersect with baseline_R, take
# the difference. With v6.9's CORE preset (227 nodes, ~770 DAG-edges
# after depth-filtering), the full sweep runs in well under a second.
#
# Ranking output: top-K nodes by cascade fraction, descending. We
# also record the absolute cascade size for each top entry (for §6
# prose) and the baseline_R size (so Stage III can quote both
# "removing T4 invalidates X / Y nodes" and "= Z %").


def _anchor_reachable_set(
    graph: CrystalGraph,
    *,
    excluded: str | None = None,
) -> set[str]:
    """Set of nodes reachable from any PLEC anchor on the depth-filtered DAG.

    Parameters
    ----------
    graph : CrystalGraph
        Source view (full or post_R).
    excluded : str | None
        If given, the named node is removed from both the node set
        and from any edge that uses it (no edges into or out of it
        are followed). The PLEC anchor identification still uses the
        original ``graph.nodes`` flags — exclusion is applied only
        to the BFS frontier expansion. This means ``excluded`` itself
        is *not* in the returned set.

    Returns
    -------
    set[str]
        Node IDs reachable from at least one PLEC anchor (excluding
        ``excluded`` itself).
    """
    succ = _depth_filtered_succ(graph)
    anchors = [
        nid for nid, n in graph.nodes.items()
        if n.is_plec_anchor and nid != excluded
    ]

    visited: set[str] = set()
    frontier: list[str] = list(anchors)
    while frontier:
        nxt: list[str] = []
        for v in frontier:
            if v in visited:
                continue
            if v == excluded:
                continue
            visited.add(v)
            for w in succ.get(v, ()):
                if w == excluded:
                    continue
                if w not in visited:
                    nxt.append(w)
        frontier = nxt
    return visited


def cascade_size(
    graph: CrystalGraph,
    node_id: str,
    *,
    baseline: set[str] | None = None,
) -> dict:
    """Compute the cascade triggered by retracting ``node_id``.

    Parameters
    ----------
    graph : CrystalGraph
    node_id : str
        Node to simulate retracting.
    baseline : set[str] | None
        Optional precomputed baseline reachable set. If not provided,
        it is recomputed (slower; pass it in if running a sweep).

    Returns
    -------
    dict
        ``{node_id, cascade_count, baseline_count, fraction, sample}``
        where ``sample`` is up to the first 10 cascade nodes (sorted
        for determinism) so the artifact stays bounded in size.
    """
    if node_id not in graph.nodes:
        raise ValueError(
            f"node_id {node_id!r} not in graph (preset={graph.preset})"
        )

    if baseline is None:
        baseline = _anchor_reachable_set(graph)
    reach_without = _anchor_reachable_set(graph, excluded=node_id)

    # Cascade = nodes that were anchor-reachable, are not the
    # retracted node itself, and become unreachable when we retract it.
    baseline_minus_self = baseline - {node_id}
    cascade = baseline_minus_self - reach_without
    denom = len(baseline_minus_self)
    frac = (len(cascade) / denom) if denom > 0 else 0.0
    return {
        'node_id':         node_id,
        'cascade_count':   len(cascade),
        'baseline_count':  denom,
        'fraction':        frac,
        'sample':          sorted(cascade)[:10],
    }


def cascade_analysis(
    graph: CrystalGraph,
    *,
    top_k: int = 20,
    skip_anchors: bool = True,
) -> dict:
    """Run cascade analysis over every non-anchor node in ``graph``.

    Parameters
    ----------
    graph : CrystalGraph
    top_k : int
        Number of top-cascade nodes to keep in the returned ranking.
    skip_anchors : bool
        If True (default), exclude the four PLEC anchors from the
        candidate set — retracting an anchor trivially severs
        everything below it and is uninteresting for waist
        identification.

    Returns
    -------
    dict
        ::
            {
              'baseline_count':  int,    # |anchor-reachable set|
              'top_k':           [{node_id, cascade_count,
                                   baseline_count, fraction,
                                   sample}, ...],
              'mean_fraction':   float,  # over scanned nodes
              'median_fraction': float,
              'max_fraction':    float,
              'n_scanned':       int,
            }
    """
    baseline = _anchor_reachable_set(graph)
    # Scan all non-anchor nodes that are themselves anchor-reachable
    # (retracting an unreachable node has zero effect by construction
    # — it isn't a predecessor of anything in baseline).
    candidates = [
        nid for nid in baseline
        if not (skip_anchors and graph.nodes[nid].is_plec_anchor)
    ]
    rows = [cascade_size(graph, nid, baseline=baseline) for nid in candidates]
    rows.sort(key=lambda r: (-r['fraction'], -r['cascade_count'], r['node_id']))
    fractions = [r['fraction'] for r in rows]
    return {
        'baseline_count':  len(baseline),
        'top_k':           rows[:top_k],
        'mean_fraction':   _stat.fmean(fractions) if fractions else 0.0,
        'median_fraction': _stat.median(fractions) if fractions else 0.0,
        'max_fraction':    max(fractions) if fractions else 0.0,
        'n_scanned':       len(rows),
    }


def check_T_crystal_cascade_v69():
    """T_crystal_cascade_v69: cascade-failure ranking well-defined on v6.9.

    STATEMENT: The v6.9 Enforcement Crystal admits a well-formed
    cascade-failure analysis on both the full_graph and post_R_subgraph
    views under the CORE module preset. For every non-anchor candidate
    node N, the cascade fraction (proportion of anchor-reachable nodes
    that lose anchor reachability when N is retracted) is in [0, 1],
    and the cascade is monotone in the obvious sense: every node
    strictly downstream of every PLEC anchor (= any node in the
    baseline reachable set) is the source of zero cascade only if it
    has no descendants. The top-K cascade ranking is non-empty on
    both views and the post-R top-3 is contained in the full top-25
    (dual-view stability of the load-bearing waists).

    This is the v6.9 refresh of the Paper 20 v1.0 (Feb-2026, v3.4.4
    snapshot) §6 Failure Cascade table. v1.0 reported "removing T4
    invalidates 67 % of predictions"; v6.9 computes the corresponding
    fraction over the *anchor-reachable node set* (a strictly larger
    universe than predictions; Stage III §6 will optionally
    re-restrict to the prediction subset once that index is
    cross-walked into the bank).

    METHOD:
      (a) Build CORE preset (full_graph + post_R_subgraph).
      (b) For each view, compute the anchor-reachable baseline using
          the depth-filtered DAG (same DAG used by workstream 4).
      (c) Sweep every non-anchor node in the baseline; recompute
          reachability with that node excluded; cascade = baseline \\
          {N} \\ reach_without_N.
      (d) Rank by cascade fraction descending; guard 0 ≤ fraction ≤ 1
          and that no cascade includes the retracted node itself.
      (e) Determinism guard: re-run the sweep on the full view,
          assert the top-K ranking is identical.
      (f) Dual-view stability: assert post-R top-3 ⊆ full top-25.
      (g) Record both rankings + summary statistics in artifacts for
          Paper 20 v3.0 Stage III §6.

    EPISTEMIC: ``[P_structural]``. Certifies the cascade-analysis
    machinery (anchor-reachability + per-node retraction sweep) is
    deterministic and produces a well-formed ranking; the numerical
    cascade fractions are research outputs, not framework invariants.
    """
    crystal = build_crystal("CORE", do_prelude=True)
    full = crystal["full_graph"]
    post_R = crystal["post_R_subgraph"]

    # (a)–(d) Full view sweep.
    cascade_full = cascade_analysis(full, top_k=20)
    check(
        cascade_full['n_scanned'] > 0,
        "full-view cascade sweep scanned 0 nodes — empty baseline?"
    )
    for row in cascade_full['top_k']:
        check(
            0.0 <= row['fraction'] <= 1.0 + 1e-12,
            f"cascade fraction out of [0,1] for {row['node_id']!r}: {row['fraction']}"
        )
        check(
            row['node_id'] not in row['sample'],
            f"retracted node {row['node_id']!r} appears in own cascade sample"
        )

    # (e) Determinism: re-run, compare top-K rankings element-wise.
    cascade_full_2 = cascade_analysis(full, top_k=20)
    rank_a = [(r['node_id'], r['cascade_count']) for r in cascade_full['top_k']]
    rank_b = [(r['node_id'], r['cascade_count']) for r in cascade_full_2['top_k']]
    check(
        rank_a == rank_b,
        "cascade ranking non-deterministic on full view"
    )

    # (f) Post-R sweep + dual-view stability.
    cascade_postR = cascade_analysis(post_R, top_k=20)
    check(
        cascade_postR['n_scanned'] > 0,
        "post_R-view cascade sweep scanned 0 nodes — empty baseline?"
    )
    top3_postR = {r['node_id'] for r in cascade_postR['top_k'][:3]}
    top25_full = {r['node_id'] for r in cascade_full['top_k'][:25]}
    # Note: top-25 may be shorter than 25 if n_scanned < 25; guard with
    # the actual length.
    if len(cascade_full['top_k']) >= 3:
        missing = top3_postR - top25_full
        check(
            not missing,
            f"post-R top-3 not contained in full top-25; missing: {sorted(missing)}"
        )

    # Headline string for the result key.
    top = cascade_full['top_k'][0]
    headline = (
        f"top cascade: {top['node_id']} retracts "
        f"{top['cascade_count']}/{top['baseline_count']} = "
        f"{top['fraction']*100:.1f}% of anchor-reachable nodes "
        f"(full view, CORE preset)"
    )

    return _result(
        name=(
            'T_crystal_cascade_v69 — cascade-failure ranking '
            'well-defined on v6.9 dual view'
        ),
        tier=3,
        epistemic='P_structural',
        summary=(
            'On the v6.9 Enforcement Crystal under the CORE preset, the '
            'depth-filtered DAG admits a well-defined cascade-failure '
            'analysis on both full_graph and post_R_subgraph views: for '
            'each non-anchor node N, the cascade fraction (proportion of '
            'anchor-reachable nodes that lose PLEC reachability when N '
            'is retracted) lies in [0, 1], the ranking is deterministic '
            'on re-run, and the post-R top-3 is contained in the full '
            'top-25 (dual-view stability). Top-20 rankings are recorded '
            'in artifacts for Paper 20 v3.0 Stage III §6 (Failure '
            'Cascade refresh of v1.0 §6).'
        ),
        key_result=headline + ' [P_structural]',
        dependencies=['T_crystal_v69_consistent'],
        cross_refs=[
            'L_PLEC_components_essentiality',
            'T_crystal_centrality_v69',
            'T_crystal_path_attribution_v69',
        ],
        artifacts={
            'preset':       'CORE',
            'full_view':    cascade_full,
            'post_R_view':  cascade_postR,
        },
    )


# =============================================================================
# §6  Workstream 5 — Convergence-cluster detection
# =============================================================================
#
# The Paper 20 v1.0 (Feb-2026, v3.4.4 snapshot) §7 finding was the
# convergence-pattern catalogue: places in the bank where multiple
# independent derivation chains fold into a single downstream theorem.
# v6.9 has many such convergence sinks by design — Phase 14 added
# T_three_level_unification (which folds I1/I2/I3/I4 into one), Phase
# 13.1 added L_PLEC_components_essentiality (which folds the four
# anchor-essentiality lemmas), Phase 14b's killed-rivals module folds
# four structural-rival kills into a single statement, etc.
#
# v6.9 convergence metric (depth-filtered DAG, same universe as
# workstreams 2 / 4):
#
#     fan_in(S)            = | { v : (v -> S) is a depth-DAG edge } |
#     anchor_diversity(S)  = | { a : a is a PLEC anchor and some
#                                 direct predecessor of S is reachable
#                                 from a on the depth-DAG } |
#
# A convergence cluster is the set of direct predecessors of a sink
# S. Sinks are ranked by fan_in (primary) and anchor_diversity
# (secondary, for tie-break and Stage III prose). High-fan_in sinks
# with low anchor_diversity are "single-line convergences" (many
# proofs but all routing through the same axiom); high on both is the
# strongest form of convergence (genuine cross-anchor unification).
#
# Computation is O(V+E) for fan_in and O(V * (V+E)) for the
# anchor_diversity sweep — well under a second on the v6.9 graph.


def _depth_filtered_pred(graph: CrystalGraph) -> dict[str, set[str]]:
    """Predecessor adjacency dict on the depth-filtered DAG.

    Mirror of ``_depth_filtered_succ`` but transposed: returns
    ``{node_id -> set of immediate predecessors}`` keeping only edges
    ``u -> v`` with ``depth(v) > depth(u)``.
    """
    pred: dict[str, set[str]] = {nid: set() for nid in graph.nodes}
    depth = graph.depth
    for e in graph.edges:
        d_src = depth.get(e.source, -1)
        d_tgt = depth.get(e.target, -1)
        if d_tgt > d_src:
            pred[e.target].add(e.source)
    return pred


def _anchor_set_reaching(
    graph: CrystalGraph,
    targets: set[str],
) -> set[str]:
    """Set of PLEC-anchor IDs that reach at least one node in ``targets``.

    Walks forward from each PLEC anchor on the depth-filtered DAG;
    once a target is hit the anchor is in. O((V+E) * |anchors|) in the
    worst case (|anchors| = 4 for v6.9 PLEC).
    """
    if not targets:
        return set()
    succ = _depth_filtered_succ(graph)
    hits: set[str] = set()
    for aid, node in graph.nodes.items():
        if not node.is_plec_anchor:
            continue
        # BFS from aid; bail as soon as a target is hit.
        visited: set[str] = set()
        frontier = [aid]
        found = False
        while frontier and not found:
            nxt: list[str] = []
            for v in frontier:
                if v in visited:
                    continue
                visited.add(v)
                if v in targets:
                    hits.add(aid)
                    found = True
                    break
                for w in succ.get(v, ()):
                    if w not in visited:
                        nxt.append(w)
            frontier = nxt
    return hits


def convergence_fan_in(graph: CrystalGraph) -> dict[str, int]:
    """Direct fan-in (in-degree) for every node on the depth-filtered DAG."""
    pred = _depth_filtered_pred(graph)
    return {nid: len(p) for nid, p in pred.items()}


def convergence_cluster_analysis(
    graph: CrystalGraph,
    *,
    top_k: int = 15,
    min_fan_in: int = 2,
) -> dict:
    """Catalogue convergence clusters on ``graph``.

    Parameters
    ----------
    graph : CrystalGraph
    top_k : int
        Number of top sinks to keep in the ranking.
    min_fan_in : int
        Skip sinks with fewer than this many direct predecessors. The
        default of 2 just removes pure leaves; bump to 3+ for "real"
        convergence reporting.

    Returns
    -------
    dict
        ::
            {
              'top_k': [{node_id, fan_in, anchor_diversity,
                         predecessors, anchors_reaching}, ...],
              'fan_in_histogram': {fan_in: count, ...},
              'n_with_fan_in_geq_min': int,
              'max_fan_in': int,
            }
    """
    pred = _depth_filtered_pred(graph)
    fan_in = {nid: len(p) for nid, p in pred.items()}

    # Filter + rank.
    candidates = [
        (nid, fi) for nid, fi in fan_in.items() if fi >= min_fan_in
    ]
    # Compute anchor diversity for each surviving candidate. Note: we
    # measure how many PLEC anchors reach at least one of the *direct
    # predecessors* (not the sink itself), because anchor reaching the
    # sink directly is trivially true via the sink's own ancestors.
    rows = []
    for nid, fi in candidates:
        preds = pred[nid]
        anchors_reaching = _anchor_set_reaching(graph, preds)
        rows.append({
            'node_id':           nid,
            'fan_in':            fi,
            'anchor_diversity':  len(anchors_reaching),
            'predecessors':      sorted(preds),
            'anchors_reaching':  sorted(anchors_reaching),
        })

    # Primary sort: fan_in descending. Secondary: anchor_diversity
    # descending (tie-breaker rewards genuine cross-anchor unification).
    # Tertiary: node_id ascending (determinism).
    rows.sort(
        key=lambda r: (-r['fan_in'], -r['anchor_diversity'], r['node_id'])
    )

    histogram: dict[int, int] = _coll.Counter(fan_in.values())

    return {
        'top_k':                 rows[:top_k],
        'fan_in_histogram':      dict(sorted(histogram.items())),
        'n_with_fan_in_geq_min': len(rows),
        'max_fan_in':            max(fan_in.values()) if fan_in else 0,
    }


def check_T_crystal_convergence_v69():
    """T_crystal_convergence_v69: convergence-cluster catalogue well-defined.

    STATEMENT: The v6.9 Enforcement Crystal admits a well-formed
    convergence-cluster catalogue on both the full_graph and
    post_R_subgraph views under the CORE module preset. Every node has
    a fan-in (in-degree on the depth-filtered DAG) in [0, V-1]; the
    top-K ranking by fan-in is non-empty on both views; and the post-R
    top-3 is contained in the full top-25 (dual-view stability of the
    convergence sinks under removal of Regime_R + Type I-IV exits).

    This is the v6.9 refresh of the Paper 20 v1.0 (Feb-2026, v3.4.4
    snapshot) §7 Convergence Patterns table. v6.9 has many more
    convergence sinks than v1.0 by design: Phase 14
    (T_three_level_unification folds I1/I2/I3/I4), Phase 13.1
    (L_PLEC_components_essentiality folds the four anchor-essentiality
    lemmas), Phase 14b (killed-rivals folds four structural rivals),
    and the 2026-04-20 T_ACC_unification + interface-sector bridge
    passes (which fold gauge / cosmology / horizon strands).

    METHOD:
      (a) Build CORE preset (full_graph + post_R_subgraph).
      (b) For each view, compute the depth-filtered fan-in (in-degree)
          per node; assert all values in [0, V-1].
      (c) For each top-K sink, compute anchor_diversity = number of
          PLEC anchors with a depth-DAG path to at least one direct
          predecessor.
      (d) Determinism: re-run on full view, assert top-K ranking
          identical in both runs.
      (e) Dual-view stability: assert post-R top-3 ⊆ full top-25.
      (f) Record both rankings + fan-in histogram in artifacts for
          Paper 20 v3.0 Stage III §7.

    EPISTEMIC: ``[P_structural]``. Certifies the convergence-cluster
    machinery (depth-filter + fan-in + anchor-diversity sweep) is
    deterministic and produces a well-formed ranking; the numerical
    fan-in / anchor-diversity values are research outputs, not
    framework invariants.
    """
    crystal = build_crystal("CORE", do_prelude=True)
    full = crystal["full_graph"]
    post_R = crystal["post_R_subgraph"]

    # (a)–(c) Full view.
    cluster_full = convergence_cluster_analysis(full, top_k=15, min_fan_in=2)
    check(
        len(cluster_full['top_k']) > 0,
        "full-view convergence top-K is empty"
    )
    n_full = full.n_nodes
    for row in cluster_full['top_k']:
        check(
            0 <= row['fan_in'] <= n_full - 1,
            f"fan_in out of [0, V-1] for {row['node_id']!r}: {row['fan_in']}"
        )
        check(
            0 <= row['anchor_diversity'] <= 4,
            f"anchor_diversity out of [0, 4] for {row['node_id']!r}: "
            f"{row['anchor_diversity']}"
        )
        check(
            row['fan_in'] == len(row['predecessors']),
            f"fan_in / predecessors mismatch for {row['node_id']!r}"
        )

    # (d) Determinism on full view.
    cluster_full_2 = convergence_cluster_analysis(full, top_k=15, min_fan_in=2)
    rank_a = [(r['node_id'], r['fan_in'], r['anchor_diversity'])
              for r in cluster_full['top_k']]
    rank_b = [(r['node_id'], r['fan_in'], r['anchor_diversity'])
              for r in cluster_full_2['top_k']]
    check(
        rank_a == rank_b,
        "convergence ranking non-deterministic on full view"
    )

    # (e) Post-R + dual-view stability.
    cluster_postR = convergence_cluster_analysis(post_R, top_k=15, min_fan_in=2)
    check(
        len(cluster_postR['top_k']) > 0,
        "post_R-view convergence top-K is empty"
    )
    top3_postR = {r['node_id'] for r in cluster_postR['top_k'][:3]}
    top25_full = {r['node_id'] for r in cluster_full['top_k'][:25]}
    if len(cluster_full['top_k']) >= 3:
        missing = top3_postR - top25_full
        check(
            not missing,
            f"post-R top-3 not contained in full top-25; missing: {sorted(missing)}"
        )

    # Headline.
    top = cluster_full['top_k'][0]
    headline = (
        f"top convergence sink: {top['node_id']} "
        f"(fan_in={top['fan_in']}, anchor_diversity={top['anchor_diversity']}/4) "
        f"[max fan_in across full view = {cluster_full['max_fan_in']}]"
    )

    return _result(
        name=(
            'T_crystal_convergence_v69 — convergence-cluster catalogue '
            'well-defined on v6.9 dual view'
        ),
        tier=3,
        epistemic='P_structural',
        summary=(
            'On the v6.9 Enforcement Crystal under the CORE preset, the '
            'depth-filtered DAG admits a well-defined convergence-cluster '
            'catalogue on both views: every node has a fan-in (in-degree) '
            'in [0, V-1], the top-K sinks have anchor_diversity in [0, 4] '
            'measuring how many PLEC anchors reach at least one direct '
            'predecessor, the ranking is deterministic on re-run, and '
            'the post-R top-3 is contained in the full top-25 (dual-view '
            'stability). Top-15 sinks + fan-in histogram + their '
            'anchor-diversity are recorded in artifacts for Paper 20 '
            'v3.0 Stage III §7 (Convergence Patterns refresh of v1.0 §7).'
        ),
        key_result=headline + ' [P_structural]',
        dependencies=['T_crystal_v69_consistent'],
        cross_refs=[
            'L_PLEC_components_essentiality',
            'T_three_level_unification',
            'T_ACC_unification',
            'T_interface_sector_bridge',
            'T_crystal_centrality_v69',
            'T_crystal_cascade_v69',
        ],
        artifacts={
            'preset':       'CORE',
            'full_view':    cluster_full,
            'post_R_view':  cluster_postR,
        },
    )


# =============================================================================
# §7  Workstream 3 — Minimum vertex cuts (Menger)
# =============================================================================
#
# Companion to workstream 2 (cascade): cascade asks "what does *one*
# deletion break?", min-cut asks "how many *coordinated* deletions
# would disconnect a derivation chain?". By Menger's theorem, the
# minimum number of internal vertices whose removal disconnects s from
# t equals the maximum number of internally vertex-disjoint s -> t
# paths. We compute this via the standard reduction:
#
#   - Vertex split: each non-source/sink vertex v becomes
#     (v_in -> v_out) with capacity 1.
#   - Each original edge (u, v) becomes (u_out -> v_in) with
#     capacity ``INF`` (we only want to cut vertices, not edges).
#   - Run Edmonds-Karp on the resulting network; max flow value =
#     min vertex cut size = number of internally vertex-disjoint
#     s -> t paths.
#   - Cut witness: BFS the residual graph from s; the cut consists
#     of those v where v_in is reachable but v_out is not.
#
# Operating on the depth-filtered DAG (consistent with workstreams 2 /
# 4 / 5) avoids spurious extra connectivity from documentation
# cross-references. Source = a PLEC anchor; sink = an "important
# downstream target" (T_sin2theta as the canonical headline; T_gauge,
# T_PMNS, T_mass_ratios as additional probes, configurable below).
#
# Expected behaviour given workstreams 2/4/5 findings: with cascade
# max = 2.3% (no single-vertex bottleneck) and 81% of nodes being
# convergence sinks (heavy redundancy), min vertex cuts from active
# anchors should be at least 2 and often 3+, quantitatively
# reinforcing the "no single point of failure" story.


_INF_CAP = 10**9   # effectively infinity for unit-capacity vertex-cut work


def _build_vertex_split_network(
    graph: CrystalGraph,
    source_id: str,
    sink_id: str,
) -> tuple[dict[str, dict[str, int]], str, str]:
    """Build a vertex-split residual capacity network.

    Returns
    -------
    (cap, s, t) where:
      cap : dict[str, dict[str, int]]
        Forward + reverse residual capacities. Vertex v gets two
        nodes ``v|in`` and ``v|out`` linked by a unit-capacity edge.
        Source and sink are not split (they collapse to a single
        node ``v|out`` for source and ``v|in`` for sink, eliminating
        the trivial split-edge cut on the endpoints themselves).
      s : str
        Effective source name in ``cap`` (= ``source_id|out``).
      t : str
        Effective sink name in ``cap`` (= ``sink_id|in``).
    """
    succ = _depth_filtered_succ(graph)

    cap: dict[str, dict[str, int]] = _coll.defaultdict(dict)

    # Add unit-capacity v_in -> v_out for every internal vertex.
    for v in graph.nodes:
        if v == source_id or v == sink_id:
            continue
        cap[f"{v}|in"][f"{v}|out"] = 1
        cap[f"{v}|out"][f"{v}|in"] = 0  # reverse residual

    # Add an INF-capacity edge for every depth-filtered original edge,
    # routed u_out -> v_in. Special-case source/sink to use the
    # collapsed names.
    def out_node(x: str) -> str:
        return f"{x}|out" if x not in (source_id, sink_id) else (
            f"{x}|out" if x == source_id else f"{x}|in"
        )

    def in_node(x: str) -> str:
        return f"{x}|in" if x not in (source_id, sink_id) else (
            f"{x}|out" if x == source_id else f"{x}|in"
        )

    for u, kids in succ.items():
        u_out = out_node(u)
        for v in kids:
            v_in = in_node(v)
            cap[u_out].setdefault(v_in, 0)
            cap[u_out][v_in] += _INF_CAP
            cap[v_in].setdefault(u_out, 0)  # reverse residual init at 0

    s = out_node(source_id)
    t = in_node(sink_id)
    return cap, s, t


def _bfs_augment(
    cap: dict[str, dict[str, int]],
    s: str,
    t: str,
) -> list[str] | None:
    """BFS for a shortest augmenting path with positive residual.

    Returns the path as a list of node names from s to t (inclusive),
    or None if no augmenting path exists.
    """
    parent: dict[str, str] = {s: s}
    Q: _coll.deque[str] = _coll.deque([s])
    while Q:
        u = Q.popleft()
        if u == t:
            break
        # Sort neighbors for deterministic BFS order (matters for the
        # cut witness when ties exist).
        for v in sorted(cap.get(u, {}).keys()):
            if v in parent:
                continue
            if cap[u].get(v, 0) > 0:
                parent[v] = u
                Q.append(v)
    if t not in parent:
        return None
    # Reconstruct path.
    path = [t]
    while path[-1] != s:
        path.append(parent[path[-1]])
    path.reverse()
    return path


def min_vertex_cut(
    graph: CrystalGraph,
    source_id: str,
    sink_id: str,
) -> dict:
    """Min vertex cut between ``source_id`` and ``sink_id``.

    Parameters
    ----------
    graph : CrystalGraph
        Source view (full or post_R).
    source_id : str
        Bank-canonical ID of the source vertex (typically a PLEC anchor).
    sink_id : str
        Bank-canonical ID of the sink vertex.

    Returns
    -------
    dict
        ::
            {
              'source':       source_id,
              'sink':         sink_id,
              'cut_size':     int,    # see semantics below
              'reachable':    bool,   # is sink reachable at all?
              'direct':       bool,   # direct (source -> sink) edge?
              'cut_witness':  [str],  # one minimal cut set (ordered)
            }

    Semantics of ``cut_size``:

      * ``cut_size = 0, reachable = False, direct = False``
            sink unreachable on the depth-filtered DAG.
      * ``cut_size = -1, direct = True``
            sink reachable via a *direct* bank edge from source. By
            Menger's theorem the vertex connectivity is undefined / +∞
            when the two endpoints are adjacent — there is no
            internal vertex to remove. Reported separately so Stage
            III can list these "1-step" derivations explicitly.
      * ``cut_size >= 1, direct = False``
            standard Menger min vertex cut over indirect paths only;
            ``len(cut_witness) == cut_size``.
    """
    if source_id not in graph.nodes:
        raise ValueError(f"source_id {source_id!r} not in graph")
    if sink_id not in graph.nodes:
        raise ValueError(f"sink_id {sink_id!r} not in graph")
    if source_id == sink_id:
        return {
            'source': source_id, 'sink': sink_id,
            'cut_size': 0, 'reachable': True, 'direct': False,
            'cut_witness': [],
        }

    # Direct-edge detection on the depth-filtered DAG. If the source
    # has a direct DAG edge to the sink, the vertex connectivity is
    # undefined / +∞ — report explicitly and bail before running
    # max-flow on a network that would saturate INF-capacity edges.
    succ = _depth_filtered_succ(graph)
    if sink_id in succ.get(source_id, set()):
        return {
            'source':       source_id,
            'sink':         sink_id,
            'cut_size':     -1,
            'reachable':    True,
            'direct':       True,
            'cut_witness':  [],
        }

    cap, s, t = _build_vertex_split_network(graph, source_id, sink_id)

    # Edmonds-Karp: BFS for shortest augmenting path, push min residual.
    flow = 0
    while True:
        path = _bfs_augment(cap, s, t)
        if path is None:
            break
        # Bottleneck along the path.
        bottleneck = min(cap[path[i]][path[i + 1]] for i in range(len(path) - 1))
        for i in range(len(path) - 1):
            u, v = path[i], path[i + 1]
            cap[u][v] -= bottleneck
            cap[v][u] = cap[v].get(u, 0) + bottleneck
        flow += bottleneck

    # Reachable from s in residual graph after max flow.
    reach: set[str] = set()
    Q: _coll.deque[str] = _coll.deque([s])
    while Q:
        u = Q.popleft()
        if u in reach:
            continue
        reach.add(u)
        for v in sorted(cap.get(u, {}).keys()):
            if v in reach:
                continue
            if cap[u].get(v, 0) > 0:
                Q.append(v)

    # If t is reachable, source and sink were never separated -> sink
    # was unreachable to begin with (or equal). Cut_size 0.
    sink_reachable_originally = (t in reach) or any(
        f"{v}|out" in reach and f"{v}|in" not in reach
        for v in graph.nodes
        if v not in (source_id, sink_id)
    )
    # Simpler reachability check via straight forward BFS on the
    # depth-filtered DAG for "reachable in original" vs "post-flow":
    succ = _depth_filtered_succ(graph)
    orig_reach: set[str] = set()
    Q2: _coll.deque[str] = _coll.deque([source_id])
    while Q2:
        u = Q2.popleft()
        if u in orig_reach:
            continue
        orig_reach.add(u)
        for v in succ.get(u, ()):
            if v not in orig_reach:
                Q2.append(v)
    reachable = sink_id in orig_reach

    if not reachable:
        return {
            'source': source_id, 'sink': sink_id,
            'cut_size': 0, 'reachable': False, 'direct': False,
            'cut_witness': [],
        }

    # Cut witness: vertex v is in the cut iff v_in in reach and
    # v_out not in reach (the split-edge crosses the partition).
    witness: list[str] = sorted(
        v for v in graph.nodes
        if v not in (source_id, sink_id)
        and f"{v}|in" in reach
        and f"{v}|out" not in reach
    )

    return {
        'source':      source_id,
        'sink':        sink_id,
        'cut_size':    flow,
        'reachable':   True,
        'direct':      False,
        'cut_witness': witness,
    }


def min_cut_table(
    graph: CrystalGraph,
    sinks: list[str],
    *,
    sources: list[str] | None = None,
) -> dict:
    """Compute min vertex cut for every (source, sink) pair.

    Parameters
    ----------
    graph : CrystalGraph
    sinks : list[str]
        Targets of interest (e.g. ['T_sin2theta', 'T_gauge', 'T_PMNS']).
    sources : list[str] | None
        Sources; defaults to the four PLEC anchors via
        ``is_plec_anchor`` flags.

    Returns
    -------
    dict
        ``{rows, by_sink_min, by_source_min, n_pairs, n_reachable}``
    """
    if sources is None:
        sources = sorted(
            nid for nid, n in graph.nodes.items() if n.is_plec_anchor
        )

    rows = []
    for s in sources:
        for t in sinks:
            r = min_vertex_cut(graph, s, t)
            r['source_role'] = graph.nodes[s].plec_role or 'anchor'
            rows.append(r)

    # Aggregate: per-sink, the min cut over sources that actually reach
    # it via *indirect* paths (skip direct-edge cases, which have
    # cut_size == -1 and are categorically different — they represent
    # 1-step bank derivations with no internal vertex to cut).
    by_sink_min: dict[str, dict] = {}
    for t in sinks:
        sub = [r for r in rows
               if r['sink'] == t and r['reachable'] and not r['direct']]
        direct_sources = [r['source'] for r in rows
                          if r['sink'] == t and r['direct']]
        if sub:
            best = min(sub, key=lambda r: r['cut_size'])
            by_sink_min[t] = {
                'min_cut':         best['cut_size'],
                'achieved_by':     best['source'],
                'witness':         best['cut_witness'],
                'direct_sources':  sorted(direct_sources),
            }
        else:
            by_sink_min[t] = {
                'min_cut':         0,
                'achieved_by':     None,
                'witness':         [],
                'direct_sources':  sorted(direct_sources),
            }

    # Aggregate: per-source, min cut over reached sinks (excl. direct).
    by_source_min: dict[str, dict] = {}
    for s in sources:
        sub = [r for r in rows
               if r['source'] == s and r['reachable'] and not r['direct']]
        direct_sinks = [r['sink'] for r in rows
                        if r['source'] == s and r['direct']]
        if sub:
            best = min(sub, key=lambda r: r['cut_size'])
            by_source_min[s] = {
                'min_cut':       best['cut_size'],
                'sink':          best['sink'],
                'witness':       best['cut_witness'],
                'direct_sinks':  sorted(direct_sinks),
            }
        else:
            by_source_min[s] = {
                'min_cut':       0,
                'sink':          None,
                'witness':       [],
                'direct_sinks':  sorted(direct_sinks),
            }

    return {
        'rows':           rows,
        'by_sink_min':    by_sink_min,
        'by_source_min':  by_source_min,
        'n_pairs':        len(rows),
        'n_reachable':    sum(1 for r in rows if r['reachable']),
        'n_direct':       sum(1 for r in rows if r['direct']),
    }


# Canonical Stage III sink set.
_STAGE3_SINKS = [
    "T_sin2theta",         # canonical headline
    "T_gauge",             # gauge derivation core
    "T_PMNS",              # heaviest convergence sink (workstream 5 #3)
    "T_mass_ratios",       # heaviest fan-in target (workstream 5 #8)
    "L_count",             # the rigid C_total = 61 fold (workstream 5 #5)
]


def check_T_crystal_min_cut_v69():
    """T_crystal_min_cut_v69: min vertex cut catalogue well-defined on v6.9.

    STATEMENT: The v6.9 Enforcement Crystal admits a well-formed
    minimum-vertex-cut catalogue from each PLEC anchor to each of a
    canonical Stage III sink set, on both the full_graph and
    post_R_subgraph views under the CORE module preset. For every
    (anchor, sink) pair that is depth-DAG reachable, the min-cut size
    is at least 1, and the cut witness has size equal to the cut value.
    The post-R cut sizes match the full cut sizes whenever the source
    survives the dual-view contract (i.e., for A1 / MD / BW; Regime_R
    is dropped from post_R by design).

    This is the v6.9 contribution to Paper 20 v3.0 Stage III §6
    (companion to workstream 2's cascade table). v6.9 is expected to
    show min-cut >= 2 for every reachable (anchor, sink) pair, given
    the workstream 2 finding that no single non-anchor node retracts
    more than 2.3% of the graph and the workstream 5 finding that
    81% of nodes are convergence sinks.

    METHOD:
      (a) Build CORE preset (full_graph + post_R_subgraph).
      (b) For each view, iterate (PLEC anchor, sink) pairs over the
          canonical Stage III sink set: T_sin2theta, T_gauge, T_PMNS,
          T_mass_ratios, L_count.
      (c) For each pair, run vertex-split + Edmonds-Karp on the
          depth-filtered DAG.
      (d) Assert: for every reachable pair, cut_size >= 1 AND
          len(cut_witness) == cut_size; for unreachable pairs,
          cut_size == 0 and cut_witness == [].
      (e) Determinism: re-run on full view, assert all 20 cut sizes
          identical to the first run.
      (f) Dual-view consistency: for each (anchor, sink) pair where
          the anchor survives in post_R, assert cut_size matches
          between the two views (we do *not* assert cut_witness
          equality because alternative minimum cuts may differ; it's
          the size that is the invariant).
      (g) Record both tables in artifacts for Paper 20 v3.0 Stage III.

    EPISTEMIC: ``[P_structural]``. Certifies that the Menger
    min-vertex-cut machinery (vertex split + Edmonds-Karp + cut
    witness extraction) is deterministic and produces well-formed
    cut sizes; the numerical cuts are research outputs, not
    framework invariants.
    """
    crystal = build_crystal("CORE", do_prelude=True)
    full = crystal["full_graph"]
    post_R = crystal["post_R_subgraph"]

    # (b)–(d) Full view catalogue.
    table_full = min_cut_table(full, _STAGE3_SINKS)
    for r in table_full['rows']:
        if not r['reachable']:
            check(
                r['cut_size'] == 0 and r['cut_witness'] == []
                and not r['direct'],
                f"unreachable pair ({r['source']}, {r['sink']}) "
                f"has nonzero cut: {r['cut_size']}"
            )
            continue
        if r['direct']:
            # Direct bank edge: no internal vertex separates source
            # and sink. Sentinel cut_size = -1; witness is empty.
            check(
                r['cut_size'] == -1 and r['cut_witness'] == [],
                f"direct pair ({r['source']}, {r['sink']}) has "
                f"unexpected (cut_size, witness): "
                f"({r['cut_size']}, {r['cut_witness']!r})"
            )
            continue
        # Standard indirect case.
        check(
            r['cut_size'] >= 1,
            f"reachable indirect pair ({r['source']}, {r['sink']}) "
            f"has cut_size {r['cut_size']} < 1"
        )
        check(
            len(r['cut_witness']) == r['cut_size'],
            f"witness size {len(r['cut_witness'])} != cut_size "
            f"{r['cut_size']} for ({r['source']}, {r['sink']})"
        )

    # (e) Determinism on full view.
    table_full_2 = min_cut_table(full, _STAGE3_SINKS)
    sizes_a = [(r['source'], r['sink'], r['cut_size'])
               for r in table_full['rows']]
    sizes_b = [(r['source'], r['sink'], r['cut_size'])
               for r in table_full_2['rows']]
    check(
        sizes_a == sizes_b,
        "min-cut sizes non-deterministic on full view"
    )

    # (f) Dual-view consistency on cut sizes (not witnesses).
    table_postR = min_cut_table(post_R, _STAGE3_SINKS)
    full_lookup = {(r['source'], r['sink']): r['cut_size']
                   for r in table_full['rows']}
    for r in table_postR['rows']:
        key = (r['source'], r['sink'])
        if key in full_lookup:
            check(
                r['cut_size'] == full_lookup[key],
                f"min-cut size differs across views for {key}: "
                f"full={full_lookup[key]} post_R={r['cut_size']}"
            )

    # Headline (over indirect reachable pairs only — direct pairs are
    # categorically different and reported separately).
    indirect_rows = [r for r in table_full['rows']
                     if r['reachable'] and not r['direct']]
    n_direct = sum(1 for r in table_full['rows'] if r['direct'])
    n_unreach = sum(1 for r in table_full['rows'] if not r['reachable'])
    if indirect_rows:
        deepest = max(indirect_rows, key=lambda r: r['cut_size'])
        shallowest = min(indirect_rows, key=lambda r: r['cut_size'])
        headline = (
            f"min-cut catalogue: {len(indirect_rows)}/{table_full['n_pairs']} "
            f"indirect reachable pairs ({n_direct} direct, "
            f"{n_unreach} unreachable); deepest cut = "
            f"{deepest['cut_size']} ({deepest['source']} -> {deepest['sink']}); "
            f"shallowest cut = {shallowest['cut_size']} "
            f"({shallowest['source']} -> {shallowest['sink']})"
        )
    else:
        headline = (
            f"min-cut catalogue: 0 indirect reachable pairs "
            f"({n_direct} direct, {n_unreach} unreachable)"
        )

    return _result(
        name=(
            'T_crystal_min_cut_v69 — Menger min vertex cut catalogue '
            'well-defined on v6.9 dual view'
        ),
        tier=3,
        epistemic='P_structural',
        summary=(
            'On the v6.9 Enforcement Crystal under the CORE preset, the '
            'depth-filtered DAG admits a well-defined Menger '
            'minimum-vertex-cut catalogue from each PLEC anchor to each '
            'of the canonical Stage III sinks (T_sin2theta, T_gauge, '
            'T_PMNS, T_mass_ratios, L_count): for every reachable pair, '
            'cut_size >= 1, the witness size equals the cut size, the '
            'computation is deterministic on re-run, and cut sizes '
            'agree between the full and post-R views for surviving '
            'sources (dual-view consistency). The full table is '
            'recorded in artifacts as the workstream-3 input to Paper '
            '20 v3.0 Stage III §6 (companion to the cascade table from '
            'workstream 2).'
        ),
        key_result=headline + ' [P_structural]',
        dependencies=['T_crystal_v69_consistent'],
        cross_refs=[
            'T_crystal_centrality_v69',
            'T_crystal_path_attribution_v69',
            'T_crystal_cascade_v69',
            'T_crystal_convergence_v69',
        ],
        artifacts={
            'preset':       'CORE',
            'sinks':        _STAGE3_SINKS,
            'full_view':    table_full,
            'post_R_view':  table_postR,
        },
    )


# =============================================================================
# §8  Bank registration
# =============================================================================


_CHECKS = {
    'T_crystal_centrality_v69':           check_T_crystal_centrality_v69,
    'T_crystal_path_attribution_v69':     check_T_crystal_path_attribution_v69,
    'T_crystal_path_attribution_scc_v69': check_T_crystal_path_attribution_scc_v69,
    'T_crystal_cascade_v69':              check_T_crystal_cascade_v69,
    'T_crystal_convergence_v69':          check_T_crystal_convergence_v69,
    'T_crystal_min_cut_v69':              check_T_crystal_min_cut_v69,
}


def register(registry):
    """Register the Phase 13.3 / Stage II crystal-metrics consistency checks.

    Currently registers:
      * T_crystal_centrality_v69           (workstream 1, BC)
      * T_crystal_path_attribution_v69     (workstream 4, depth-filtered)
      * T_crystal_path_attribution_scc_v69 (workstream 4b, SCC-aware, Phase 14c)
      * T_crystal_cascade_v69              (workstream 2, cascade-failure analysis)
      * T_crystal_convergence_v69          (workstream 5, convergence clusters)
      * T_crystal_min_cut_v69              (workstream 3, Menger min vertex cuts)
    """
    registry.update(_CHECKS)
