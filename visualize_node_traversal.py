"""Create a self-contained HTML visualization from a solver trace.

Run:
    python visualize_node_traversal.py allele_trace.json traversal.html

This file does not import or run SCIP. It only reads the solver's JSON output.
"""

import argparse
import json
from pathlib import Path


HTML_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>SCIP Node Traversal</title>
  <style>
    :root { color-scheme: light; font-family: Inter, system-ui, sans-serif; }
    body { margin: 0; background: #f5f7fb; color: #172033; }
    header { padding: 24px 28px 14px; }
    h1 { margin: 0 0 8px; font-size: 25px; }
    .summary { display: flex; flex-wrap: wrap; gap: 10px; }
    .pill { background: white; border: 1px solid #dce2ee; border-radius: 10px;
            padding: 9px 12px; box-shadow: 0 2px 8px #1f293711; }
    main { display: grid; grid-template-columns: minmax(600px, 1fr) 310px;
           gap: 16px; padding: 10px 28px 28px; }
    .panel { background: white; border: 1px solid #dce2ee; border-radius: 12px;
             box-shadow: 0 3px 14px #1f29370c; }
    #tree-wrap { overflow: auto; min-height: 560px; }
    #tree { display: block; }
    #details { padding: 18px; }
    #details h2 { font-size: 18px; margin: 0 0 14px; }
    #details dl { display: grid; grid-template-columns: 90px 1fr; gap: 9px; }
    #details dt { color: #667085; }
    #details dd { margin: 0; overflow-wrap: anywhere; }
    .edge { stroke: #b9c2d0; stroke-width: 2; }
    .node { cursor: pointer; stroke: #344054; stroke-width: 2; }
    .node:hover { stroke-width: 4; }
    .node-label { pointer-events: none; text-anchor: middle; dominant-baseline: middle;
                  fill: white; font-size: 12px; font-weight: 700; }
    .order-label { pointer-events: none; text-anchor: middle; fill: #475467;
                   font-size: 11px; }
    .legend { padding-top: 18px; color: #667085; font-size: 13px; line-height: 1.5; }
    @media (max-width: 900px) { main { grid-template-columns: 1fr; padding: 10px; }
                                header { padding: 20px 10px 10px; } }
  </style>
</head>
<body>
  <header>
    <h1>SCIP branch-and-bound traversal</h1>
    <div class="summary" id="summary"></div>
  </header>
  <main>
    <section class="panel" id="tree-wrap"><svg id="tree"></svg></section>
    <aside class="panel" id="details">
      <h2>Node details</h2>
      <p>Click a node in the tree.</p>
      <div class="legend">
        Blue nodes were focused by SCIP. Gray nodes were created but were not
        observed as focused. The number below a node is its traversal order.
      </div>
    </aside>
  </main>
  <script>
    const data = __TRACE_DATA__;
    const nodes = data.trace.nodes;
    const byId = new Map(nodes.map(node => [node.id, node]));
    const children = new Map(nodes.map(node => [node.id, []]));
    nodes.forEach(node => {
      if (node.parent_id !== null && children.has(node.parent_id)) {
        children.get(node.parent_id).push(node);
      }
    });
    children.forEach(group => group.sort((a, b) => a.id - b.id));

    const roots = nodes.filter(node => node.parent_id === null || !byId.has(node.parent_id));
    let nextLeaf = 0;
    const positions = new Map();

    function place(node) {
      const kids = children.get(node.id) || [];
      let x;
      if (kids.length === 0) {
        x = nextLeaf++;
      } else {
        kids.forEach(place);
        x = kids.reduce((sum, kid) => sum + positions.get(kid.id).x, 0) / kids.length;
      }
      positions.set(node.id, {x, y: node.depth});
    }
    roots.forEach(place);

    const xGap = 105;
    const yGap = 105;
    const margin = 55;
    const width = Math.max(650, nextLeaf * xGap + margin * 2);
    const maxDepth = Math.max(0, ...nodes.map(node => node.depth));
    const height = Math.max(560, (maxDepth + 1) * yGap + margin * 2);
    const svg = document.getElementById("tree");
    svg.setAttribute("width", width);
    svg.setAttribute("height", height);
    svg.setAttribute("viewBox", `0 0 ${width} ${height}`);

    function point(node) {
      const p = positions.get(node.id);
      return {x: margin + p.x * xGap, y: margin + p.y * yGap};
    }
    function svgElement(name, attributes) {
      const element = document.createElementNS("http://www.w3.org/2000/svg", name);
      Object.entries(attributes).forEach(([key, value]) => element.setAttribute(key, value));
      return element;
    }

    nodes.forEach(node => {
      if (node.parent_id === null || !byId.has(node.parent_id)) return;
      const from = point(byId.get(node.parent_id));
      const to = point(node);
      svg.appendChild(svgElement("line", {
        x1: from.x, y1: from.y + 22, x2: to.x, y2: to.y - 22, class: "edge"
      }));
    });

    nodes.forEach(node => {
      const p = point(node);
      const visited = node.focus_order !== undefined;
      const circle = svgElement("circle", {
        cx: p.x, cy: p.y, r: 22, class: "node",
        fill: visited ? "#2563eb" : "#98a2b3"
      });
      circle.addEventListener("click", () => showDetails(node));
      svg.appendChild(circle);

      const id = svgElement("text", {x: p.x, y: p.y, class: "node-label"});
      id.textContent = node.id;
      svg.appendChild(id);

      if (visited) {
        const order = svgElement("text", {x: p.x, y: p.y + 40, class: "order-label"});
        order.textContent = `visit ${node.focus_order}`;
        svg.appendChild(order);
      }
    });

    const stats = data.scip_statistics;
    const solution = data.solution;
    const maxFrontier = Math.max(0, ...data.trace.traversal.map(step => step.frontier_size));
    document.getElementById("summary").innerHTML = [
      `Alleles: <strong>${solution.selected_alleles.join(" + ") || "none"}</strong>`,
      `Mismatch: <strong>${solution.objective ?? "n/a"}</strong>`,
      `Processed nodes: <strong>${stats.total_nodes}</strong>`,
      `LP iterations: <strong>${stats.total_lp_iterations}</strong>`,
      `Max sampled frontier: <strong>${maxFrontier}</strong>`
    ].map(text => `<div class="pill">${text}</div>`).join("");

    function showDetails(node) {
      const step = data.trace.traversal.find(item => item.node_id === node.id);
      document.getElementById("details").innerHTML = `
        <h2>Node ${node.id}</h2>
        <dl>
          <dt>Visit order</dt><dd>${node.focus_order ?? "not focused"}</dd>
          <dt>Parent</dt><dd>${node.parent_id ?? "none"}</dd>
          <dt>Depth</dt><dd>${node.depth}</dd>
          <dt>Decision</dt><dd>${node.branch_decision}</dd>
          <dt>Frontier</dt><dd>${step?.frontier_size ?? "n/a"}</dd>
          <dt>LP iterations</dt><dd>${step?.lp_iterations_cumulative ?? "n/a"}</dd>
        </dl>
        <div class="legend">Blue nodes were focused by SCIP. Gray nodes were
        created but were not observed as focused. The number below a node is
        its traversal order.</div>`;
    }
  </script>
</body>
</html>
"""


def create_html(trace_path, html_path):
    """Read one trace JSON file and write one standalone HTML file."""
    trace_data = json.loads(Path(trace_path).read_text())
    html = HTML_TEMPLATE.replace("__TRACE_DATA__", json.dumps(trace_data))
    Path(html_path).write_text(html)
    print(f"Visualization written to: {html_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("trace", help="JSON created by solve_allele_selection.py")
    parser.add_argument("output", nargs="?", default="node_traversal.html")
    args = parser.parse_args()
    create_html(args.trace, args.output)


if __name__ == "__main__":
    main()
