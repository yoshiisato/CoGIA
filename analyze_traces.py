"""Analyze SCIP solver traces to identify potential side-channel leakage.

This script compares traces from different genomic inputs to detect patterns
in execution time, memory usage, node traversal, and other system metrics that
could leak information about private genomic data.

Usage:
    python analyze_traces.py traces/*.json
"""

import argparse
import json
from pathlib import Path
from collections import defaultdict
import statistics


def load_trace(trace_path):
    """Load a JSON trace file."""
    with open(trace_path) as f:
        return json.load(f)


def extract_metrics(trace):
    """Extract key metrics from a trace for comparison."""
    solution = trace.get("solution", {})
    scip_stats = trace.get("scip_statistics", {})
    system_metrics = trace.get("system_metrics", {})
    trace_data = trace.get("trace", {})
    
    return {
        "selected_configs": tuple(sorted(solution.get("selected_configurations", []))),
        "objective": solution.get("objective"),
        "total_nodes": scip_stats.get("total_nodes"),
        "lp_iterations": scip_stats.get("total_lp_iterations"),
        "total_time_ms": system_metrics.get("total_time_ms"),
        "peak_memory_mb": system_metrics.get("peak_memory_mb"),
        "traversal_order": [t["node_id"] for t in trace_data.get("traversal", [])],
        "frontier_sizes": [t["frontier_size"] for t in trace_data.get("traversal", [])],
        "memory_snapshots": trace_data.get("memory_snapshots", []),
    }


def compare_traces(traces):
    """Compare traces and identify potential leakage patterns."""
    
    print("=" * 80)
    print("SIDE-CHANNEL LEAKAGE ANALYSIS")
    print("=" * 80)
    
    # Group by selected configurations
    by_config = defaultdict(list)
    metrics_by_trace = {}
    
    for trace_path, trace in traces.items():
        metrics = extract_metrics(trace)
        metrics_by_trace[trace_path] = metrics
        by_config[metrics["selected_configs"]].append((trace_path, metrics))
    
    print(f"\nAnalyzed {len(traces)} traces")
    print(f"Found {len(by_config)} distinct solution configurations\n")
    
    # Analysis 1: Same solution, different execution patterns?
    print("-" * 80)
    print("ANALYSIS 1: Same genomic input → Same solution (baseline)")
    print("-" * 80)
    
    normal_traces = [t for t in traces if "normal" in t]
    if len(normal_traces) >= 2:
        print(f"\nComparing {len(normal_traces)} 'normal' traces (same input):")
        for i, (path, metrics) in enumerate([(p, metrics_by_trace[p]) for p in normal_traces]):
            print(f"  {Path(path).name}:")
            print(f"    Nodes: {metrics['total_nodes']}, Time: {metrics['total_time_ms']:.2f}ms, "
                  f"Memory: {metrics['peak_memory_mb']:.2f}MB")
        
        # Check variance
        times = [metrics_by_trace[p]["total_time_ms"] for p in normal_traces]
        nodes = [metrics_by_trace[p]["total_nodes"] for p in normal_traces]
        mems = [metrics_by_trace[p]["peak_memory_mb"] for p in normal_traces]
        
        if len(times) > 1:
            print(f"\n  Time variance: {statistics.variance(times):.4f} (std: {statistics.stdev(times):.4f})")
            print(f"  Node variance: {statistics.variance(nodes):.4f} (std: {statistics.stdev(nodes):.4f})")
            print(f"  Memory variance: {statistics.variance(mems):.4f} (std: {statistics.stdev(mems):.4f})")
    
    # Analysis 2: Different inputs → Different execution patterns?
    print("\n" + "-" * 80)
    print("ANALYSIS 2: Different genomic inputs → Execution patterns")
    print("-" * 80)
    
    for config, group in by_config.items():
        print(f"\nSolution: {' + '.join(config) if config else 'none'}")
        print(f"  Traces with this solution: {len(group)}")
        
        if len(group) > 1:
            times = [m["total_time_ms"] for _, m in group if m["total_time_ms"] is not None]
            nodes = [m["total_nodes"] for _, m in group if m["total_nodes"] is not None]
            if times:
                print(f"  Time range: {min(times):.2f}ms - {max(times):.2f}ms")
            if nodes:
                print(f"  Node range: {min(nodes)} - {max(nodes)}")
    
    # Analysis 3: Can we distinguish inputs by execution metrics?
    print("\n" + "-" * 80)
    print("ANALYSIS 3: Input distinguishability via side channels")
    print("-" * 80)
    
    all_metrics = list(metrics_by_trace.values())
    
    # Check if different inputs have distinguishable node counts
    node_counts = [m["total_nodes"] for m in all_metrics if m["total_nodes"] is not None]
    unique_nodes = set(node_counts)
    print(f"\nUnique node counts: {len(unique_nodes)} (out of {len(node_counts)} valid traces)")
    print(f"Node counts: {sorted(node_counts)}")
    
    if len(unique_nodes) == len(node_counts):
        print("  ⚠️  HIGH RISK: Each input produces unique node count!")
    elif len(unique_nodes) > len(node_counts) / 2:
        print("  ⚠️  MEDIUM RISK: Most inputs have unique node counts")
    else:
        print("  ✓ LOW RISK: Many inputs share node counts")
    
    # Check time distinguishability
    times = [m["total_time_ms"] for m in all_metrics if m["total_time_ms"] is not None]
    unique_times = len(set(round(t, 2) for t in times))
    print(f"\nUnique execution times (2dp): {unique_times} (out of {len(times)} valid traces)")
    if times:
        print(f"Time range: {min(times):.2f}ms - {max(times):.2f}ms")
    
    if unique_times == len(times):
        print("  ⚠️  HIGH RISK: Each input has unique execution time!")
    elif unique_times > len(times) / 2:
        print("  ⚠️  MEDIUM RISK: Most inputs have unique times")
    else:
        print("  ✓ LOW RISK: Many inputs share similar times")
    
    # Analysis 4: Traversal pattern analysis
    print("\n" + "-" * 80)
    print("ANALYSIS 4: Branch-and-bound traversal patterns")
    print("-" * 80)
    
    traversal_patterns = {}
    for path, metrics in metrics_by_trace.items():
        traversal = tuple(metrics["traversal_order"])
        traversal_patterns[traversal] = traversal_patterns.get(traversal, []) + [path]
    
    print(f"\nUnique traversal patterns: {len(traversal_patterns)}")
    for pattern, paths in traversal_patterns.items():
        if len(paths) > 1:
            print(f"  Pattern shared by {len(paths)} traces: {[Path(p).name for p in paths]}")
        else:
            print(f"  Unique pattern: {Path(paths[0]).name}")
    
    if len(traversal_patterns) == len(all_metrics):
        print("  ⚠️  HIGH RISK: Each input has unique traversal pattern!")
    elif len(traversal_patterns) > len(all_metrics) / 2:
        print("  ⚠️  MEDIUM RISK: Most inputs have unique traversal patterns")
    else:
        print("  ✓ LOW RISK: Many inputs share traversal patterns")
    
    # Analysis 5: Memory pattern analysis
    print("\n" + "-" * 80)
    print("ANALYSIS 5: Memory access patterns")
    print("-" * 80)
    
    for path, metrics in metrics_by_trace.items():
        snapshots = metrics["memory_snapshots"]
        if snapshots:
            mem_values = [s["memory_mb"] for s in snapshots]
            print(f"\n{Path(path).name}:")
            print(f"  Memory snapshots: {len(snapshots)}")
            print(f"  Memory range: {min(mem_values):.4f}MB - {max(mem_values):.4f}MB")
            print(f"  Memory growth: {mem_values[-1] - mem_values[0]:.4f}MB")
    
    # Analysis 6: Frontier size patterns
    print("\n" + "-" * 80)
    print("ANALYSIS 6: Search frontier size patterns")
    print("-" * 80)
    
    frontier_patterns = {}
    for path, metrics in metrics_by_trace.items():
        frontier = tuple(metrics["frontier_sizes"])
        frontier_patterns[frontier] = frontier_patterns.get(frontier, []) + [path]
    
    print(f"\nUnique frontier patterns: {len(frontier_patterns)}")
    for pattern, paths in frontier_patterns.items():
        if len(paths) > 1:
            print(f"  Pattern shared by {len(paths)} traces")
        else:
            print(f"  Unique pattern: {Path(paths[0]).name}")
    
    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY: Potential Side-Channel Leakage")
    print("=" * 80)
    
    risk_factors = []
    
    if len(unique_nodes) > len(all_metrics) / 2:
        risk_factors.append("Node count distinguishability")
    if unique_times > len(all_metrics) / 2:
        risk_factors.append("Execution time distinguishability")
    if len(traversal_patterns) > len(all_metrics) / 2:
        risk_factors.append("Traversal pattern distinguishability")
    
    if risk_factors:
        print(f"\n⚠️  DETECTED POTENTIAL LEAKAGE:")
        for factor in risk_factors:
            print(f"  - {factor}")
        print(f"\nRecommendation: These components should be made data-oblivious")
    else:
        print("\n✓ No obvious leakage patterns detected in current traces")
        print("  (Note: This is a small-scale analysis; real-world scenarios may differ)")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("traces", nargs="+", help="Trace JSON files to analyze")
    args = parser.parse_args()
    
    traces = {}
    for trace_path in args.traces:
        try:
            traces[trace_path] = load_trace(trace_path)
        except Exception as e:
            print(f"Error loading {trace_path}: {e}")
    
    if not traces:
        print("No valid traces to analyze")
        return
    
    compare_traces(traces)


if __name__ == "__main__":
    main()
