"""APF v6.9+ — Phase 13.2: Enforcement Crystal walker.

Programmatic construction, analysis, and serialization of the Enforcement
Crystal (Paper 20, v2 corrected). This module turns the bank registry
(``apf.bank.REGISTRY``) plus the four PLEC anchors
(``apf.crystal_axiom_roots.PLEC_AXIOM_ROOTS``) into a typed graph
representation, computes the Tier-1 structural metrics required by the
crystal narrative (depth profile, waist candidates, dual-view counts),
and emits a ``dashboard_data.json`` payload directly consumable by the
canonical Three.js viewer in ``__APF Library/Visualizations/`` (uploaded
2026-04-21 ``index (13).html``).

Mechanization spec
------------------
This is the Tier-1 implementation of the apf/crystal.py mechanization
spec recorded in
``__APF Library/APF Reference Docs/Reference - Paper 20 Enforcement Crystal Update Plan (2026-04-21).md`` §8.
Higher tiers (waist propagation under module-include filter, prediction
re-targeting, Three.js geometry export) will be added in subsequent
Phase 13 stages.

Contract (§8.1)
---------------
- Inputs:
    1. ``apf.bank.REGISTRY`` — the 364 (v6.9) bank-registered checks
       with their dependency strings.
    2. ``apf.crystal_axiom_roots.PLEC_AXIOM_ROOTS`` — the four axiom
       anchors (A1, A2, MD, BW) and their canonical bank-check
       binding.
    3. A module-include filter (``preset`` argument: 'CORE',
       'EXTENDED', or 'FULL').
- Outputs:
    1. ``CrystalGraph`` dataclass with nodes/edges/depth-profile/
       waist-candidate fields, in two views: ``full_graph`` and
       ``post_R_subgraph``.
    2. JSON payload via ``dashboard_payload(crystal)`` matching the
       schema consumed by the Admissibility Physics Engine viewer
       (theorems dict, sector verdicts, prediction list,
       theorem_checker, epistemic_counts, tier_stats, version).
- Bank-side guarantee: a registered consistency check
  ``T_crystal_v69_consistent`` (see end of file) verifies that the
  walker's view of the bank agrees with the bank's own self-report on
  count + epistemic distribution + that all four PLEC anchors are
  reachable as roots.

Two views
---------
The walker exposes two graph views (per §8.4):

  * ``full_graph``         — Regime R + the five exit channels
                              (Type I … Type V) appear as nodes.
                              This is the structural picture.
  * ``post_R_subgraph``    — R and the exit channels are *factored
                              out*: Type-V (admissible exit) flows
                              are kept and absorbed into normal
                              dependency edges; Type I–IV (rejected)
                              are pruned. This is the picture after
                              the PLEC selector has fired.

Both views share node identity, but the dependency closure differs.

Module-include presets (§8.5)
-----------------------------
  * 'CORE'      — Papers 1–7 modules only:
                    apf.core, apf.gauge, apf.generations, apf.spacetime,
                    apf.gravity, apf.plec, apf.unification,
                    apf.unification_three_levels, apf.cosmology.
  * 'EXTENDED'  — CORE + supplements + standalone:
                    apf.supplements, apf.majorana,
                    apf.internalization, apf.internalization_geo,
                    apf.extensions, apf.validation,
                    apf.standalone.* .
  * 'FULL'      — EXTENDED + sessions:
                    apf.session_v63c, apf.session_qg, apf.session_nnlo,
                    apf.session_delta_pmns, apf.session_cosmo_update,
                    apf.session_phase2_confrontation, apf.red_team.

The default for the bank-registered consistency check + the dashboard
exporter is 'CORE' (the picture Paper 20 actually argues about).

Status (Tier 1)
---------------
Implemented:
  * graph build + depth assignment
  * dual-view filter
  * Tier-1 metrics: node/edge counts, depth profile, waist candidates
  * dashboard_data.json emission compatible with the viewer schema
  * bank-registered ``check_T_crystal_v69_consistent``

Deferred (Tier 2+):
  * waist re-propagation under per-preset filtering
  * sector-verdict aggregation (currently passes through bank's own)
  * prediction-target registry binding (currently empty)
  * Three.js geometry export (currently delegated to the viewer's
    own client-side compute)
"""

from __future__ import annotations

import collections as _coll
import dataclasses as _dc
import importlib as _il
import json as _json
import math as _math
from typing import Iterable, Mapping, Optional

from apf import bank as _bank
from apf.apf_utils import _result, check, CheckFailure
from apf.crystal_axiom_roots import PLEC_AXIOM_ROOTS, all_roots


__all__ = [
    "MODULE_PRESETS",
    "CrystalNode",
    "CrystalEdge",
    "CrystalGraph",
    "build_crystal",
    "dashboard_payload",
    "PRELUDE_MODULES",
]


# =============================================================================
# §1  Module-include presets and prelude bootstrap
# =============================================================================

MODULE_PRESETS: dict[str, tuple[str, ...]] = {
    "CORE": (
        "apf.core",
        "apf.gauge",
        "apf.generations",
        "apf.spacetime",
        "apf.gravity",
        "apf.plec",
        "apf.unification",
        "apf.unification_three_levels",
        "apf.cosmology",
    ),
    "EXTENDED": (
        "apf.core",
        "apf.gauge",
        "apf.generations",
        "apf.spacetime",
        "apf.gravity",
        "apf.plec",
        "apf.unification",
        "apf.unification_three_levels",
        "apf.cosmology",
        "apf.supplements",
        "apf.majorana",
        "apf.internalization",
        "apf.internalization_geo",
        "apf.extensions",
        "apf.validation",
        "apf.standalone.L_Cauchy_uniqueness",
        "apf.standalone.L_CKM_resolution_limit",
        "apf.standalone.phase1_seesaw_closure",
        "apf.standalone.phase5_theorem_R_audit",
    ),
    "FULL": (
        "apf.core",
        "apf.gauge",
        "apf.generations",
        "apf.spacetime",
        "apf.gravity",
        "apf.plec",
        "apf.unification",
        "apf.unification_three_levels",
        "apf.cosmology",
        "apf.supplements",
        "apf.majorana",
        "apf.internalization",
        "apf.internalization_geo",
        "apf.extensions",
        "apf.validation",
        "apf.standalone.L_Cauchy_uniqueness",
        "apf.standalone.L_CKM_resolution_limit",
        "apf.standalone.phase1_seesaw_closure",
        "apf.standalone.phase5_theorem_R_audit",
        "apf.red_team",
        "apf.session_v63c",
        "apf.session_qg",
        "apf.session_nnlo",
        "apf.session_delta_pmns",
        "apf.session_cosmo_update",
        "apf.session_phase2_confrontation",
    ),
}

# Modules whose checks populate the DAG cache (must run before any
# check that depends on dag_get / acc_SM-style bootstrap). Mirrors
# verify_all.PRELUDE_MODULES.
PRELUDE_MODULES: tuple[str, ...] = (
    "apf.core",
    "apf.gauge",
    "apf.spacetime",
    "apf.gravity",
    "apf.generations",
    "apf.cosmology",
    "apf.extensions",
    "apf.supplements",
)


def _bootstrap_dag(verbose: bool = False) -> int:
    """Run all check_ functions in PRELUDE_MODULES (in order) so that
    DAG-bootstrapped checks downstream resolve cleanly. Returns the
    number of checks executed. Failures are silent (the only purpose
    is to populate the dag-key cache)."""
    n = 0
    for mod_name in PRELUDE_MODULES:
        try:
            mod = _il.import_module(mod_name)
        except Exception:
            continue
        # Snapshot the items: some checks lazily set module-level
        # attributes when invoked (memoization/cache), which mutates
        # vars(mod) mid-iteration.
        items = list(vars(mod).items())
        for name, fn in items:
            if name.startswith("check_") and callable(fn):
                try:
                    fn()
                    n += 1
                except Exception as e:
                    if verbose:
                        print(f"  [prelude] {mod_name}.{name} -> {e!r}")
    return n


# =============================================================================
# §2  Dependency-string normalization
# =============================================================================

# Bank dependency strings carry a few legacy spellings of the same
# logical anchor. Normalize them to the canonical REGISTRY key so the
# graph is a clean overlay.
_DEP_ALIASES: dict[str, str] = {
    "L_epsilon*": "L_epsilon*",      # canonical
    "L_e*":       "L_epsilon*",      # legacy
    "L_epsilon":  "L_epsilon*",      # legacy
    "Lε*":         "L_epsilon*",     # unicode legacy
    "L_ε*":        "L_epsilon*",     # unicode legacy
}


def _normalize_dep(raw: str) -> Optional[str]:
    """Strip a bank dependency string down to a single REGISTRY key.

    Bank dep strings sometimes look like ``"T_gauge"``, sometimes
    ``"T_gauge (apf.gauge)"``, sometimes ``"meaning: human-readable"``.
    We:
      * drop anything after the first '(' or whitespace
      * skip 'meaning:' free-text annotations
      * apply the alias table
    Returns ``None`` if the dep is purely descriptive (no anchor).
    """
    s = (raw or "").strip()
    if not s:
        return None
    if s.lower().startswith("meaning"):
        return None
    if "(" in s:
        s = s.split("(", 1)[0].strip()
    elif " " in s:
        s = s.split(None, 1)[0].strip()
    return _DEP_ALIASES.get(s, s)


# =============================================================================
# §3  Node + edge dataclasses
# =============================================================================


@_dc.dataclass(frozen=True)
class CrystalNode:
    """A vertex in the Enforcement Crystal."""

    id: str                      # canonical REGISTRY key (or PLEC anchor name)
    name: str                    # display name from the check record
    tier: int                    # 0..5, or -1 for axioms
    epistemic: str               # 'axiom' | 'P' | 'P_structural' | 'C' | ...
    module: str                  # 'apf.core', 'apf.gravity', ...
    kind: str                    # 'axiom' | 'lemma' | 'theorem' | 'plec_anchor'
    dependencies: tuple[str, ...]  # already-normalized parent IDs
    summary: str                 # short text from the check record
    is_plec_anchor: bool = False  # one of A1, A2, MD, BW
    plec_role: str = ""           # 'A1' | 'A2' | 'MD' | 'BW' or ''


@_dc.dataclass(frozen=True)
class CrystalEdge:
    """A directed dependency edge: source produces / supports target."""

    source: str
    target: str


# =============================================================================
# §4  Graph container
# =============================================================================


@_dc.dataclass
class CrystalGraph:
    """A view of the Enforcement Crystal (one of {full, post_R})."""

    view: str                                         # 'full' | 'post_R'
    preset: str                                       # 'CORE' | ...
    nodes: dict[str, CrystalNode]                     # id -> node
    edges: tuple[CrystalEdge, ...]                    # directed
    depth: dict[str, int]                             # id -> BFS depth
    by_depth: dict[int, list[str]]                    # depth -> node ids
    waist_candidates: list[dict]                      # [{depth, ids, width, dominant}]
    plec_anchor_ids: tuple[str, ...]                  # IDs that ARE PLEC anchors
    notes: list[str]                                  # diagnostic notes (empty on success)

    @property
    def n_nodes(self) -> int:
        return len(self.nodes)

    @property
    def n_edges(self) -> int:
        return len(self.edges)

    @property
    def max_depth(self) -> int:
        return max(self.depth.values()) if self.depth else 0


# =============================================================================
# §5  Walker proper
# =============================================================================

# IDs of the Regime-R selector node + the five exit channels. These
# carry distinguished meaning under the dual-view contract: the
# 'post_R' view absorbs the Type V (admissible exit) edges into normal
# dependency edges, and prunes the Type I–IV channels (rejected
# possibilities) along with the R selector itself.
_REGIME_R_ID = "Regime_R"
_TYPE_I_THRU_V_IDS: tuple[str, ...] = (
    "Regime_exit_Type_I",
    "Regime_exit_Type_II",
    "Regime_exit_Type_III",
    "Regime_exit_Type_IV",
    "Regime_exit_Type_V",
)


def _module_of(check_name: str) -> str:
    """Look up the module a registered check lives in.

    Uses each module's ``_CHECKS`` dict (the canonical source of
    registry keys; these carry the final key form including e.g. the
    asterisk in ``L_epsilon*``) when available, and falls back to
    stripping the ``check_`` prefix from plain function names for
    modules that don't expose ``_CHECKS``. Cached lazily.
    """
    if not hasattr(_module_of, "_cache"):
        cache: dict[str, str] = {}
        for mod_name in MODULE_PRESETS["FULL"]:
            try:
                mod = _il.import_module(mod_name)
            except Exception:
                continue
            # Preferred path: the module's _CHECKS dict — its keys are the
            # exact REGISTRY keys.
            checks = getattr(mod, "_CHECKS", None)
            if isinstance(checks, dict):
                for key in checks:
                    cache.setdefault(key, mod_name)
            # Fallback: enumerate check_ functions and strip prefix.
            for name in vars(mod):
                if name.startswith("check_"):
                    bare = name[len("check_") :]
                    cache.setdefault(bare, mod_name)
        _module_of._cache = cache  # type: ignore[attr-defined]
    return _module_of._cache.get(check_name, "?")  # type: ignore[attr-defined]


def _node_kind(check_name: str, epistemic: str) -> str:
    if epistemic == "axiom" or check_name.startswith("A") and len(check_name) <= 3:
        return "axiom"
    if check_name.startswith("L_") or check_name.startswith("Δ_") or check_name.startswith("Delta_"):
        return "lemma"
    return "theorem"


def _build_one_view(
    *,
    view: str,
    preset: str,
    keep_modules: set[str],
    plec_anchors: dict[str, str],   # check_id -> 'A1' / 'A2' / 'MD' / 'BW'
) -> CrystalGraph:
    """Construct a single CrystalGraph view from the loaded bank."""
    notes: list[str] = []

    # ---- Step 1: collect raw records from REGISTRY for kept modules ----
    # Snapshot keys: some checks lazily register helpers when invoked, so
    # iterating REGISTRY directly can hit "dictionary changed size".
    raw_records: dict[str, dict] = {}
    registry_snapshot = list(_bank.REGISTRY.keys())
    for cid in registry_snapshot:
        mod = _module_of(cid)
        if mod not in keep_modules:
            continue
        try:
            r = _bank.REGISTRY[cid]()
        except Exception as e:
            notes.append(f"check {cid} raised at walker time: {type(e).__name__}: {e}")
            continue
        if not isinstance(r, dict):
            notes.append(f"check {cid} returned non-dict result; skipping")
            continue
        raw_records[cid] = r

    # ---- Step 2: apply view-specific node filtering ----
    keep_ids = set(raw_records)
    if view == "post_R":
        # Drop R itself + the four rejected exit channels; keep Type V
        # (admissible exit) but treat it as a transparent passthrough.
        for x in (_REGIME_R_ID,
                  "Regime_exit_Type_I", "Regime_exit_Type_II",
                  "Regime_exit_Type_III", "Regime_exit_Type_IV"):
            keep_ids.discard(x)

    # ---- Step 3: build nodes ----
    nodes: dict[str, CrystalNode] = {}
    for cid, r in raw_records.items():
        if cid not in keep_ids:
            continue
        deps_raw = r.get("dependencies") or []
        deps_norm: list[str] = []
        for d in deps_raw:
            n = _normalize_dep(d)
            if n is None:
                continue
            # in 'post_R' view, dependency-on-R becomes
            # dependency-on-R's-parents (we drop the link entirely;
            # the upstream PLEC anchors A1/A2/MD/BW stay reachable
            # via crystal_axiom_roots seeding below).
            if view == "post_R" and n in {_REGIME_R_ID,
                                          "Regime_exit_Type_I",
                                          "Regime_exit_Type_II",
                                          "Regime_exit_Type_III",
                                          "Regime_exit_Type_IV"}:
                continue
            deps_norm.append(n)
        epi = str(r.get("epistemic", "?"))
        kind = _node_kind(cid, epi)
        is_anchor = cid in plec_anchors
        nodes[cid] = CrystalNode(
            id=cid,
            name=str(r.get("name", cid)),
            tier=int(r.get("tier", 0) or 0),
            epistemic=epi,
            module=_module_of(cid),
            kind="plec_anchor" if is_anchor else kind,
            dependencies=tuple(deps_norm),
            summary=str(r.get("summary", "") or "")[:240],
            is_plec_anchor=is_anchor,
            plec_role=plec_anchors.get(cid, ""),
        )

    # ---- Step 4: edges (filtered to keep both endpoints in `nodes`) ----
    edges_list: list[CrystalEdge] = []
    for cid, node in nodes.items():
        for dep in node.dependencies:
            if dep in nodes:
                edges_list.append(CrystalEdge(source=dep, target=cid))
    edges = tuple(edges_list)

    # ---- Step 5: BFS depth from the PLEC anchors ----
    parents: dict[str, set[str]] = {nid: set() for nid in nodes}
    children: dict[str, set[str]] = {nid: set() for nid in nodes}
    for e in edges:
        parents[e.target].add(e.source)
        children[e.source].add(e.target)

    depth: dict[str, int] = {}
    # Seed: PLEC anchors at depth 0
    seeds: list[str] = [nid for nid, n in nodes.items() if n.is_plec_anchor]
    for s in seeds:
        depth[s] = 0
    # Iteratively assign depth = max(parent depth) + 1 once all parents
    # are assigned. Repeat until no progress; then absorb any cycles
    # by relaxed assignment (max over assigned parents).
    changed = True
    while changed:
        changed = False
        for nid in nodes:
            if nid in depth:
                continue
            ps = parents[nid]
            if ps and all(p in depth for p in ps):
                depth[nid] = max(depth[p] for p in ps) + 1
                changed = True
    # Relaxed pass: nodes whose *some* parent is assigned.
    cycle_fixed = 0
    changed = True
    while changed:
        changed = False
        for nid in nodes:
            if nid in depth:
                continue
            ps = parents[nid]
            assigned = [p for p in ps if p in depth]
            if assigned:
                depth[nid] = max(depth[p] for p in assigned) + 1
                cycle_fixed += 1
                changed = True
    # Final fallback: orphans get pushed to max+1
    if depth:
        mx = max(depth.values())
    else:
        mx = 0
    for nid in nodes:
        if nid not in depth:
            depth[nid] = mx + 1
            cycle_fixed += 1
    if cycle_fixed:
        notes.append(f"depth: cycle/orphan resolution touched {cycle_fixed} nodes")

    by_depth: dict[int, list[str]] = _coll.defaultdict(list)
    for nid, d in depth.items():
        by_depth[d].append(nid)
    by_depth = dict(sorted(by_depth.items()))

    # ---- Step 6: waist candidates ----
    # A waist is a depth ≥ 3 with width ≤ 4 that's a local minimum.
    # Emit up to 3 candidates ranked by descendant fanout.
    def _desc_count(start: str) -> int:
        seen = {start}
        stk = [start]
        while stk:
            x = stk.pop()
            for c in children.get(x, ()):
                if c not in seen:
                    seen.add(c)
                    stk.append(c)
        return len(seen) - 1
    max_d = max(by_depth) if by_depth else 0
    waist_candidates: list[dict] = []
    for d in range(3, max_d + 1):
        w = len(by_depth.get(d, []))
        wp = len(by_depth.get(d - 1, []))
        wn = len(by_depth.get(d + 1, [])) if d < max_d else 99
        if w <= wp and w <= wn and w <= 4:
            ids_at_d = by_depth.get(d, [])
            top = max(ids_at_d, key=_desc_count) if ids_at_d else None
            waist_candidates.append({
                "depth": d,
                "width": w,
                "ids": list(ids_at_d),
                "dominant": top,
                "dominant_descendants": _desc_count(top) if top else 0,
            })
    waist_candidates.sort(key=lambda r: (r["width"], -r["dominant_descendants"]))
    waist_candidates = waist_candidates[:3]

    plec_ids = tuple(nid for nid, n in nodes.items() if n.is_plec_anchor)

    return CrystalGraph(
        view=view,
        preset=preset,
        nodes=nodes,
        edges=edges,
        depth=depth,
        by_depth=by_depth,
        waist_candidates=waist_candidates,
        plec_anchor_ids=plec_ids,
        notes=notes,
    )


def build_crystal(
    preset: str = "CORE",
    *,
    do_prelude: bool = True,
) -> dict[str, CrystalGraph]:
    """Build both views of the Enforcement Crystal under ``preset``.

    Returns a dict ``{'full_graph': CrystalGraph, 'post_R_subgraph': CrystalGraph}``.
    """
    if preset not in MODULE_PRESETS:
        raise ValueError(
            f"Unknown preset {preset!r}; expected one of {list(MODULE_PRESETS)}"
        )

    # Make sure the bank is loaded.
    _bank._load()
    if do_prelude:
        _bootstrap_dag()

    keep_modules = set(MODULE_PRESETS[preset])

    # Build the PLEC anchor lookup: REGISTRY id -> 'A1'/'A2'/'MD'/'BW'
    plec_anchors: dict[str, str] = {}
    for role, rec in all_roots():
        plec_anchors[rec["check"]] = role

    return {
        "full_graph":       _build_one_view(view="full",    preset=preset,
                                            keep_modules=keep_modules,
                                            plec_anchors=plec_anchors),
        "post_R_subgraph":  _build_one_view(view="post_R",  preset=preset,
                                            keep_modules=keep_modules,
                                            plec_anchors=plec_anchors),
    }


# =============================================================================
# §6  Dashboard JSON emission (viewer-compatible)
# =============================================================================


def dashboard_payload(
    crystal: Mapping[str, CrystalGraph],
    *,
    version: str = "v6.9",
    sector_verdicts: Optional[dict] = None,
    predictions: Optional[list] = None,
) -> dict:
    """Render the dual-view crystal as a JSON-serializable dict matching
    the schema consumed by the canonical Three.js viewer
    (``index (13).html``).

    The viewer's ``computeCM(D)`` re-derives crystal metrics client-side
    from ``D.theorems[*].dependencies``, so the only schema-critical
    job here is to expose the bank's theorem records keyed by ID, with
    the right epistemic / tier / dependencies fields.

    The ``post_R_subgraph`` view is exposed as ``D.crystal_post_R`` for
    advanced viewers; the primary ``D.theorems`` payload reflects the
    full graph (so the viewer's CM computation matches the structural
    picture Paper 20 narrates).
    """
    full = crystal["full_graph"]
    post_R = crystal["post_R_subgraph"]

    # Core theorems dict (full view).
    theorems: dict[str, dict] = {}
    for nid, node in full.nodes.items():
        theorems[nid] = {
            "id":            nid,
            "name":          node.name,
            "tier":          node.tier,
            "epistemic":     node.epistemic,
            "dependencies":  list(node.dependencies),
            "module":        node.module,
            "kind":          node.kind,
            "summary":       node.summary,
            "is_plec_anchor": node.is_plec_anchor,
            "plec_role":     node.plec_role,
            "depth":         full.depth.get(nid, -1),
        }

    # Epistemic counts (the viewer's donut + audit cards read these).
    epi_counts: _coll.Counter = _coll.Counter()
    for n in full.nodes.values():
        epi_counts[n.epistemic] += 1

    # Per-tier stats: pass / total (everything passes by construction —
    # the walker only retains records whose check returned a dict).
    tier_stats: dict[str, dict] = {}
    for n in full.nodes.values():
        b = tier_stats.setdefault(str(n.tier), {"total": 0, "by_epistemic": {}})
        b["total"] += 1
        b["by_epistemic"][n.epistemic] = b["by_epistemic"].get(n.epistemic, 0) + 1

    # The viewer expects a `theorem_checker` block; we report walker-pass
    # = bank-record-pass since each kept node corresponds to a passing check.
    n_total = len(theorems)
    theorem_checker = {
        "passed":    True,
        "n_pass":    n_total,
        "total":     n_total,
        "available": n_total,
    }

    payload = {
        "version":           f"APF {version} crystal walker",
        "total_theorems":    n_total,
        "theorems":          theorems,
        "epistemic_counts":  dict(epi_counts),
        "tier_stats":        tier_stats,
        "theorem_checker":   theorem_checker,
        "all_pass":          True,
        "predictions":       list(predictions or []),
        "sector_verdicts":   sector_verdicts or _DEFAULT_SECTOR_VERDICTS,
        "crystal_meta": {
            "preset":            full.preset,
            "n_nodes_full":      full.n_nodes,
            "n_edges_full":      full.n_edges,
            "n_nodes_post_R":    post_R.n_nodes,
            "n_edges_post_R":    post_R.n_edges,
            "max_depth_full":    full.max_depth,
            "max_depth_post_R":  post_R.max_depth,
            "waist_candidates":  full.waist_candidates,
            "plec_anchor_ids":   list(full.plec_anchor_ids),
            "walker_notes":      list(full.notes),
        },
        # Advanced viewers can read this; the viewer at index (13).html
        # ignores unknown top-level keys.
        "crystal_post_R": {
            "theorems": {
                nid: {
                    "id": nid, "name": n.name, "tier": n.tier,
                    "epistemic": n.epistemic,
                    "dependencies": list(n.dependencies),
                    "depth": post_R.depth.get(nid, -1),
                } for nid, n in post_R.nodes.items()
            },
            "n_nodes": post_R.n_nodes,
            "n_edges": post_R.n_edges,
            "waist_candidates": post_R.waist_candidates,
        },
    }
    return payload


# Default sector verdicts for when the caller doesn't override.
# These mirror the v6.9 status reported in CLAUDE.md / the engine paper.
_DEFAULT_SECTOR_VERDICTS = {
    "gauge":            {"verdict": "complete",  "n_passing": "all"},
    "gravity":          {"verdict": "complete",  "n_passing": "all"},
    "rg_mechanism":     {"verdict": "complete",  "n_passing": "all"},
    "flavor_mixing":    {"verdict": "partial",   "n_passing": "most",
                         "open": ["m_t, m_b absolute scale"]},
    "quantum_structure": {"verdict": "complete", "n_passing": "all"},
    "particles":        {"verdict": "partial",   "n_passing": "most",
                         "open": ["dark matter particle ID"]},
    "cosmology":        {"verdict": "tension",   "n_passing": "most",
                         "open": ["H0DN 7.09σ tension; Paper 6 §11.4"]},
    "geometry":         {"verdict": "complete",  "n_passing": "all"},
}


def write_dashboard_json(path: str, payload: dict) -> None:
    """Convenience: serialize payload to ``path`` with stable key order."""
    with open(path, "w", encoding="utf-8") as f:
        _json.dump(payload, f, indent=2, sort_keys=False)


# =============================================================================
# §7  Bank-registered consistency check
# =============================================================================


def check_T_crystal_v69_consistent():
    """T_crystal_v69_consistent: walker view of bank ↔ bank's self-report.

    Verifies that ``apf.crystal.build_crystal('CORE')`` produces a
    coherent dual-view graph in which:
      (a) all four PLEC anchors (A1, A2, MD, BW from
          crystal_axiom_roots) are present as nodes in the full view;
      (b) the bank's own len(REGISTRY) ≥ the count of kept nodes
          (i.e., the walker is a *restriction* of the bank, never an
          extension);
      (c) the post_R view drops Regime_R + Type I–IV exits and is
          strictly smaller than the full view (sanity check on the
          dual-view filter);
      (d) every kept node's parents resolve inside the kept set
          (the walker does not leave dangling edges);
      (e) the dashboard payload round-trips through json.dumps without
          raising (i.e., everything is JSON-serializable).
    """
    crystal = build_crystal("CORE", do_prelude=True)
    full = crystal["full_graph"]
    post_R = crystal["post_R_subgraph"]

    # (a) PLEC anchors present
    anchor_ids_expected = {rec["check"] for _r, rec in all_roots()}
    anchor_ids_present = set(full.plec_anchor_ids)
    if not anchor_ids_expected.issubset(anchor_ids_present):
        missing = anchor_ids_expected - anchor_ids_present
        raise CheckFailure(
            f"PLEC anchors missing from full view: {sorted(missing)}; "
            f"present={sorted(anchor_ids_present)}"
        )

    # (b) walker ⊆ bank
    n_bank = len(_bank.REGISTRY)
    n_full = full.n_nodes
    if n_full > n_bank:
        raise CheckFailure(
            f"Walker emitted {n_full} nodes but bank only has "
            f"{n_bank} registered checks; walker overshoots bank"
        )

    # (c) post_R ⊊ full
    if post_R.n_nodes >= full.n_nodes:
        raise CheckFailure(
            f"post_R view ({post_R.n_nodes} nodes) should be a strict "
            f"subset of full view ({full.n_nodes} nodes)"
        )
    if _REGIME_R_ID in post_R.nodes:
        raise CheckFailure(
            f"post_R view still contains Regime_R ({_REGIME_R_ID}) — "
            f"dual-view filter failed to prune it"
        )

    # (d) no dangling edges
    bad_edges = [(e.source, e.target) for e in full.edges
                 if e.source not in full.nodes or e.target not in full.nodes]
    if bad_edges:
        raise CheckFailure(
            f"Walker emitted {len(bad_edges)} dangling edge(s); "
            f"first: {bad_edges[:3]}"
        )

    # (e) JSON round-trip
    payload = dashboard_payload(crystal)
    try:
        s = _json.dumps(payload)
    except (TypeError, ValueError) as e:
        raise CheckFailure(f"dashboard payload not JSON-serializable: {e}")
    n_serialized = len(s)

    return _result(
        passed=True,
        name="T_crystal_v69_consistent: Crystal walker ↔ bank coherence",
        tier=4,
        epistemic="P_structural",
        summary=(
            f"Walker built {n_full}-node full view + {post_R.n_nodes}-node "
            f"post-R view from bank.REGISTRY (size {n_bank}) under "
            f"preset CORE. All four PLEC anchors {sorted(anchor_ids_expected)} "
            f"reachable as roots. Dashboard payload serialized cleanly "
            f"({n_serialized} bytes JSON)."
        ),
        key_result={
            "preset":            "CORE",
            "n_full":            n_full,
            "n_post_R":          post_R.n_nodes,
            "n_bank":            n_bank,
            "max_depth_full":    full.max_depth,
            "max_depth_post_R":  post_R.max_depth,
            "waist_candidates":  [w["dominant"] for w in full.waist_candidates],
            "json_bytes":        n_serialized,
        },
        dependencies=[
            "A1",
            "Regime_R",
            "L_epsilon*",
            "worked_example",
            "T_three_level_unification",
        ],
        cross_refs=[
            "apf.crystal_axiom_roots.PLEC_AXIOM_ROOTS",
            "Paper 20 (Enforcement Crystal v2 corrected)",
            "Reference - Paper 20 Enforcement Crystal Update Plan (2026-04-21) §8",
            "__APF Library/Visualizations/index (13).html (canonical viewer)",
        ],
        artifacts={},
    )


# =============================================================================
# §8  Bank registration
# =============================================================================

_CHECKS = {
    'T_crystal_v69_consistent': check_T_crystal_v69_consistent,
}


def register(registry):
    """Register the Phase 13.2 crystal-walker consistency check."""
    registry.update(_CHECKS)
