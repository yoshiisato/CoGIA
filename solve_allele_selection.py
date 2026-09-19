"""Solve a small allele-selection problem with SCIP.

Install:
    python -m pip install pyscipopt

Run:
    python solve_allele_selection.py --c1 40 --output c40_trace.json

The JSON output contains the solution plus the branch-and-bound nodes reported
by SCIP.  Use visualize_node_traversal.py to turn it into an HTML page.
"""

import argparse
import json
from pathlib import Path

from pyscipopt import Eventhdlr, Model, SCIP_EVENTTYPE, SCIP_PARAMSETTING, quicksum


# ---------------------------------------------------------------------------
# 1. Toy genomic data
# ---------------------------------------------------------------------------

# Each row is one fictional allele. Each bit describes one genomic position:
# 1 means the allele predicts C; 0 means it predicts T.
ALLELE_PATTERNS = [
    "00000000",  # A01
    "11111111",  # A02
    "11001010",  # A03
    "10110100",  # A04
    "01100011",  # A05
    "00011101",  # A06
    "10101010",  # A07
    "01110100",  # A08
    "11010001",  # A09
    "00101110",  # A10
    "10000111",  # A11
    "01011001",  # A12
]

ALLELES = [tuple(map(int, pattern)) for pattern in ALLELE_PATTERNS]
ALLELE_NAMES = [f"A{i + 1:02d}" for i in range(len(ALLELES))]

READS_PER_POSITION = 100
NUMBER_OF_COPIES = 2
FIXED_C_COUNTS = (48, 97, 52, 49, 51, 47, 53, 2)


def make_observations(c_at_position_1):
    """Return the eight observed C-read counts, changing only position 1."""
    if not 0 <= c_at_position_1 <= READS_PER_POSITION:
        raise ValueError("c1 must be between 0 and 100")
    return (c_at_position_1,) + FIXED_C_COUNTS[1:]


# ---------------------------------------------------------------------------
# 2. SCIP model
# ---------------------------------------------------------------------------

def build_model(observed_c):
    """Create the allele-selection optimization model."""
    model = Model("allele_selection")
    model.hideOutput()

    # Keep the search simple and reproducible for studying the tree.
    model.setIntParam("parallel/maxnthreads", 1)
    model.setIntParam("randomization/randomseedshift", 0)
    model.setPresolve(SCIP_PARAMSETTING.OFF)
    model.setHeuristics(SCIP_PARAMSETTING.OFF)
    model.setSeparating(SCIP_PARAMSETTING.OFF)
    model.setIntParam("misc/usesymmetry", 0)
    model.setIntParam("branching/mostinf/priority", 1_000_000)

    # copies[i] is 0, 1, or 2 copies of allele i.
    copies = [
        model.addVar(name=name, vtype="I", lb=0, ub=NUMBER_OF_COPIES)
        for name in ALLELE_NAMES
    ]
    model.addCons(quicksum(copies) == NUMBER_OF_COPIES, name="two_copies")

    absolute_errors = []
    for position, observed_c_count in enumerate(observed_c):
        predicted_c = 50 * quicksum(
            ALLELES[i][position] * copies[i] for i in range(len(ALLELES))
        )
        predicted_t = READS_PER_POSITION - predicted_c
        observed_t_count = READS_PER_POSITION - observed_c_count

        for base, predicted, observed in (
            ("C", predicted_c, observed_c_count),
            ("T", predicted_t, observed_t_count),
        ):
            error = model.addVar(
                name=f"error_{base}_{position + 1}", vtype="C", lb=0
            )
            model.addCons(error >= predicted - observed)
            model.addCons(error >= observed - predicted)
            absolute_errors.append(error)

    model.setObjective(quicksum(absolute_errors), "minimize")
    return model, copies


# ---------------------------------------------------------------------------
# 3. Minimal SCIP node recorder
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
            self.nodes[node.getNumber()]["focus_order"] = focus_order
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
# 4. Solve and export SCIP's result
# ---------------------------------------------------------------------------

def solve(c_at_position_1, output_path):
    observed_c = make_observations(c_at_position_1)
    model, copies = build_model(observed_c)

    recorder = NodeRecorder()
    model.includeEventhdlr(recorder, "node_recorder", "Record SCIP tree nodes")

    model.optimize()

    selected_alleles = []
    objective = None
    if model.getNSols() > 0:
        for name, variable in zip(ALLELE_NAMES, copies):
            selected_alleles.extend([name] * int(round(model.getVal(variable))))
        objective = model.getObjVal()

    output = {
        "problem": {
            "observed_c": list(observed_c),
            "observed_t": [READS_PER_POSITION - value for value in observed_c],
            "allele_patterns": dict(zip(ALLELE_NAMES, ALLELE_PATTERNS)),
        },
        "solution": {
            "status": str(model.getStatus()),
            "selected_alleles": selected_alleles,
            "objective": objective,
        },
        "scip_statistics": {
            "total_nodes": model.getNTotalNodes(),
            "total_lp_iterations": model.getNLPIterations(),
        },
        "trace": {
            "nodes": sorted(recorder.nodes.values(), key=lambda node: node["id"]),
            "traversal": recorder.traversal,
        },
    }

    output_path = Path(output_path)
    output_path.write_text(json.dumps(output, indent=2))
    model.freeProb()

    print(f"Selected alleles: {' + '.join(selected_alleles)}")
    print(f"Mismatch score: {objective}")
    print(f"Processed nodes: {output['scip_statistics']['total_nodes']}")
    print(f"Trace written to: {output_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--c1", type=int, default=48,
                        help="Observed C reads at genomic position 1")
    parser.add_argument("--output", default="allele_trace.json",
                        help="JSON file to create")
    args = parser.parse_args()
    solve(args.c1, args.output)


if __name__ == "__main__":
    main()
