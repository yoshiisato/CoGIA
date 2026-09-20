"""Solve an ALDY-like copy number estimation problem with SCIP.

This implements a simplified version of ALDY's Copy Number Estimation Problem (CNEP)
with the three-term objective function:
  o1: coverage match error
  o2: gene-pseudogene interaction error  
  o3: parsimony (minimize number of configurations)

Install:
    python -m pip install pyscipopt

Run:
    python solve_aldy_like.py --gene-cn 50 50 50 50 --pseudo-cn 20 20 20 20 --output trace.json
"""

import argparse
import json
import time
import tracemalloc
from pathlib import Path

from pyscipopt import Eventhdlr, Model, SCIP_EVENTTYPE, SCIP_PARAMSETTING, quicksum


# ---------------------------------------------------------------------------
# 1. ALDY-like genomic data (simplified CYP2D6/CYP2D7 model)
# ---------------------------------------------------------------------------

# Gene is segmented into r regions. We'll use r=4 regions for simplicity.
# Each configuration vector g_i describes which regions are present in a copy.
# g_i = [g_i,1, ..., g_i,2r] where first r are gene regions, last r are pseudogene regions.
# g_i,j ∈ {0,1} indicates if region j is present in this configuration.

# Simplified configuration vectors for CYP2D6-like gene
# Format: [gene_region_1, gene_region_2, gene_region_3, gene_region_4,
#          pseudo_region_1, pseudo_region_2, pseudo_region_3, pseudo_region_4]
CONFIGURATION_VECTORS = {
    "CN1":  [1, 1, 1, 1, 0, 0, 0, 0],  # Normal copy (gene only)
    "CN2":  [1, 1, 0, 1, 0, 0, 0, 0],  # Deletion in region 3
    "CN3":  [1, 1, 1, 0, 0, 0, 0, 0],  # Deletion in region 4
    "CN4":  [1, 1, 1, 1, 1, 0, 0, 0],  # Fusion with pseudogene region 1
    "CN5":  [1, 1, 1, 1, 0, 1, 0, 0],  # Fusion with pseudogene region 2
    "CN6":  [1, 1, 1, 1, 1, 1, 0, 0],  # Fusion with pseudogene regions 1-2
    "CN7":  [0, 0, 0, 0, 1, 1, 1, 1],  # Pseudogene only (rare)
    "CN8":  [1, 0, 1, 1, 0, 0, 0, 0],  # Deletion in region 2
    "CN9":  [1, 1, 1, 1, 0, 0, 1, 0],  # Fusion with pseudogene region 3
    "CN10": [1, 1, 1, 1, 0, 0, 0, 1],  # Fusion with pseudogene region 4
}

CONFIG_NAMES = list(CONFIGURATION_VECTORS.keys())
CONFIGS = [CONFIGURATION_VECTORS[name] for name in CONFIG_NAMES]

# Parsimony scores: higher for unlikely configurations (e.g., left fusions)
PARSIMONY_SCORES = {
    "CN1": 1.0,   # Normal - most likely
    "CN2": 2.0,   # Deletion - moderately likely
    "CN3": 2.0,
    "CN4": 3.0,   # Fusion - less likely
    "CN5": 3.0,
    "CN6": 4.0,   # Complex fusion - even less likely
    "CN7": 5.0,   # Pseudogene only - very unlikely
    "CN8": 2.0,
    "CN9": 3.0,
    "CN10": 3.0,
}

NUM_REGIONS = 4  # r
BASE_COVERAGE = 50  # Expected coverage per region per copy


def make_coverage_vector(gene_cn, pseudo_cn, noise=0):
    """
    Generate observed coverage vector cn = [cn_1, ..., cn_2r].
    
    Args:
        gene_cn: List of coverage for gene regions [cn_1, ..., cn_r]
        pseudo_cn: List of coverage for pseudogene regions [cn_{r+1}, ..., cn_{2r}]
        noise: Optional noise to add to simulate sequencing bias
    """
    import random
    random.seed(42)  # For reproducibility
    
    cn = []
    for g_cov, p_cov in zip(gene_cn, pseudo_cn):
        cn.append(g_cov + random.randint(-noise, noise) if noise else g_cov)
    for p_cov in pseudo_cn:
        cn.append(p_cov + random.randint(-noise, noise) if noise else p_cov)
    
    return tuple(cn)


# ---------------------------------------------------------------------------
# 2. ALDY-like SCIP model with three-term objective
# ---------------------------------------------------------------------------

def build_aldy_model(observed_cn):
    """
    Create ALDY-like CNEP model with three-term objective.
    
    Minimize: o1 + o2 + o3
    where:
      o1 = sum_{j=1}^{2r} |cn_j - sum_i z_i * g_i,j|
      o2 = sum_{j=1}^{r} |(cn_j - cn_{j+r}) / v_j - sum_i z_i * (g_i,j - g_i,j+r) / v_j|
      o3 = sum_i mu_i * z_i (parsimony)
    
    where v_j = max(cn_j, cn_{j+r}) + 1 (normalization factor)
    """
    model = Model("aldy_cnep")
    model.hideOutput()

    # Keep search simple and reproducible
    model.setIntParam("parallel/maxnthreads", 1)
    model.setIntParam("randomization/randomseedshift", 0)
    model.setPresolve(SCIP_PARAMSETTING.OFF)
    model.setHeuristics(SCIP_PARAMSETTING.OFF)
    model.setSeparating(SCIP_PARAMSETTING.OFF)
    model.setIntParam("misc/usesymmetry", 0)
    model.setIntParam("branching/mostinf/priority", 1_000_000)

    # Binary variables z_i for each configuration
    z_vars = [
        model.addVar(name=f"z_{name}", vtype="B", lb=0, ub=1)
        for name in CONFIG_NAMES
    ]

    # Constraint: select at least 1, at most 4 configurations (typical copy number range)
    model.addCons(quicksum(z_vars) >= 1, name="min_copies")
    model.addCons(quicksum(z_vars) <= 4, name="max_copies")

    # Calculate normalization factors v_j
    v = []
    for j in range(NUM_REGIONS):
        v_j = max(observed_cn[j], observed_cn[j + NUM_REGIONS]) + 1
        v.append(v_j)

    # Term o1: coverage match error
    o1_terms = []
    for j in range(2 * NUM_REGIONS):
        predicted = quicksum(
            CONFIGS[i][j] * z_vars[i] * BASE_COVERAGE
            for i in range(len(CONFIGS))
        )
        observed = observed_cn[j]
        
        error = model.addVar(name=f"o1_error_{j}", vtype="C", lb=0)
        model.addCons(error >= predicted - observed)
        model.addCons(error >= observed - predicted)
        o1_terms.append(error)

    # Term o2: gene-pseudogene interaction error (normalized)
    o2_terms = []
    for j in range(NUM_REGIONS):
        # Normalized observed difference
        obs_diff = (observed_cn[j] - observed_cn[j + NUM_REGIONS]) / v[j]
        
        # Normalized predicted difference
        pred_diff = quicksum(
            (CONFIGS[i][j] - CONFIGS[i][j + NUM_REGIONS]) * z_vars[i] * BASE_COVERAGE / v[j]
            for i in range(len(CONFIGS))
        )
        
        error = model.addVar(name=f"o2_error_{j}", vtype="C", lb=0)
        model.addCons(error >= pred_diff - obs_diff)
        model.addCons(error >= obs_diff - pred_diff)
        o2_terms.append(error)

    # Term o3: parsimony (weighted sum of selected configurations)
    # o3 = sum_i mu_i * z_i where mu_i is the parsimony score (constant)
    o3 = quicksum(PARSIMONY_SCORES[name] * z_vars[i] for i, name in enumerate(CONFIG_NAMES))

    # Objective: minimize o1 + o2 + o3
    o1 = quicksum(o1_terms)
    o2 = quicksum(o2_terms)
    model.setObjective(o1 + o2 + o3, "minimize")

    return model, z_vars, (o1, o2, o3)


# ---------------------------------------------------------------------------
# 3. SCIP node recorder (same as before)
# ---------------------------------------------------------------------------

def describe_node(node):
    """Copy the node fields supplied by SCIP into ordinary Python data."""
    parent = node.getParent()
    decisions = []
    branching = node.getParentBranchings()

    if branching:
        variables, bounds, bound_types = branching
        for variable, bound, bound_type in zip(variables, bounds, bound_types):
            relation = ">=" if bound_type == 0 else "<="
            decisions.append(f"{variable.name} {relation} {bound:g}")

    return {
        "id": node.getNumber(),
        "parent_id": None if parent is None else parent.getNumber(),
        "depth": node.getDepth(),
        "branch_decision": "; ".join(decisions) or "root",
    }


class NodeRecorder(Eventhdlr):
    """Record which branch-and-bound nodes SCIP creates and visits."""

    def __init__(self):
        super().__init__()
        self.nodes = {}
        self.traversal = []
        self.event_mask = SCIP_EVENTTYPE.NODEFOCUSED | SCIP_EVENTTYPE.NODEBRANCHED
        self.memory_snapshots = []
        self.start_time = None

    def eventinit(self):
        self.model.catchEvent(self.event_mask, self)

    def eventexit(self):
        self.model.dropEvent(self.event_mask, self)

    def save_node(self, node):
        snapshot = describe_node(node)
        self.nodes[snapshot["id"]] = {
            **self.nodes.get(snapshot["id"], {}),
            **snapshot,
        }

    def eventexec(self, event):
        node = event.getNode()
        self.save_node(node)

        if event.getType() == SCIP_EVENTTYPE.NODEFOCUSED:
            open_nodes = sum(len(group) for group in self.model.getOpenNodes())
            focus_order = len(self.traversal)
            
            # Capture memory snapshot at each node focus
            current_mem = tracemalloc.get_traced_memory()
            elapsed_time = time.time() - self.start_time if self.start_time else 0
            
            self.nodes[node.getNumber()]["focus_order"] = focus_order
            self.memory_snapshots.append({
                "focus_order": focus_order,
                "node_id": node.getNumber(),
                "memory_mb": current_mem[0] / 1024 / 1024,
                "peak_memory_mb": current_mem[1] / 1024 / 1024,
                "elapsed_time_ms": elapsed_time * 1000,
            })
            
            self.traversal.append({
                "order": focus_order,
                "node_id": node.getNumber(),
                "open_nodes_excluding_focus": open_nodes,
                "frontier_size": open_nodes + 1,
                "lp_iterations_cumulative": self.model.getNLPIterations(),
            })

        if event.getType() == SCIP_EVENTTYPE.NODEBRANCHED:
            for child in self.model.getChildren():
                self.save_node(child)


# ---------------------------------------------------------------------------
# 4. Solve and export
# ---------------------------------------------------------------------------

def solve(gene_coverage, pseudo_coverage, output_path):
    observed_cn = make_coverage_vector(gene_coverage, pseudo_coverage)
    model, z_vars, (o1, o2, o3) = build_aldy_model(observed_cn)

    # Start memory and time tracking
    tracemalloc.start()
    recorder = NodeRecorder()
    recorder.start_time = time.time()
    
    model.includeEventhdlr(recorder, "node_recorder", "Record SCIP tree nodes")

    model.optimize()

    selected_configs = []
    objective = None
    o1_val = o2_val = o3_val = None
    
    if model.getNSols() > 0:
        for name, variable in zip(CONFIG_NAMES, z_vars):
            if model.getVal(variable) > 0.5:
                selected_configs.append(name)
        objective = model.getObjVal()
        o1_val = model.getVal(o1)
        o2_val = model.getVal(o2)
        o3_val = model.getVal(o3)

    # Get final memory statistics
    final_mem = tracemalloc.get_traced_memory()
    total_time = time.time() - recorder.start_time
    tracemalloc.stop()

    output = {
        "problem": {
            "observed_cn": list(observed_cn),
            "gene_coverage": list(gene_coverage),
            "pseudo_coverage": list(pseudo_coverage),
            "configurations": CONFIGURATION_VECTORS,
            "parsimony_scores": PARSIMONY_SCORES,
        },
        "solution": {
            "status": str(model.getStatus()),
            "selected_configurations": selected_configs,
            "objective": objective,
            "o1_coverage_error": o1_val,
            "o2_interaction_error": o2_val,
            "o3_parsimony": o3_val,
        },
        "scip_statistics": {
            "total_nodes": model.getNTotalNodes(),
            "total_lp_iterations": model.getNLPIterations(),
        },
        "system_metrics": {
            "total_time_ms": total_time * 1000,
            "peak_memory_mb": final_mem[1] / 1024 / 1024,
            "final_memory_mb": final_mem[0] / 1024 / 1024,
        },
        "trace": {
            "nodes": sorted(recorder.nodes.values(), key=lambda node: node["id"]),
            "traversal": recorder.traversal,
            "memory_snapshots": recorder.memory_snapshots,
        },
    }

    output_path = Path(output_path)
    output_path.write_text(json.dumps(output, indent=2))
    model.freeProb()

    print(f"Selected configurations: {' + '.join(selected_configs)}")
    print(f"Total objective: {objective:.2f}")
    print(f"  Coverage error (o1): {o1_val:.2f}")
    print(f"  Interaction error (o2): {o2_val:.2f}")
    print(f"  Parsimony (o3): {o3_val:.2f}")
    print(f"Processed nodes: {output['scip_statistics']['total_nodes']}")
    print(f"Total time: {total_time*1000:.2f}ms")
    print(f"Peak memory: {final_mem[1]/1024/1024:.2f}MB")
    print(f"Trace written to: {output_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-cn", type=int, nargs=4, default=[50, 50, 50, 50],
                        help="Coverage for 4 gene regions")
    parser.add_argument("--pseudo-cn", type=int, nargs=4, default=[20, 20, 20, 20],
                        help="Coverage for 4 pseudogene regions")
    parser.add_argument("--output", default="aldy_trace.json",
                        help="JSON file to create")
    args = parser.parse_args()
    solve(tuple(args.gene_cn), tuple(args.pseudo_cn), args.output)


if __name__ == "__main__":
    main()
