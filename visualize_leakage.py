"""Visualize side-channel leakage patterns from SCIP solver traces.

Creates plots showing how different genomic inputs affect:
- Execution time
- Memory usage
- Node count
- Traversal patterns

Usage:
    python visualize_leakage.py traces/*.json
"""

import argparse
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np


def load_traces(trace_paths):
    """Load all trace files."""
    traces = {}
    for path in trace_paths:
        try:
            with open(path) as f:
                traces[path] = json.load(f)
        except Exception as e:
            print(f"Error loading {path}: {e}")
    return traces


def plot_execution_times(traces, output_dir):
    """Plot execution times for different inputs."""
    fig, ax = plt.subplots(figsize=(12, 6))
    
    names = []
    times = []
    nodes = []
    
    for path, trace in traces.items():
        name = Path(path).stem
        solution = trace.get("solution", {})
        scip_stats = trace.get("scip_statistics", {})
        system_metrics = trace.get("system_metrics", {})
        
        names.append(name)
        times.append(system_metrics.get("total_time_ms", 0))
        nodes.append(scip_stats.get("total_nodes", 0))
    
    x = np.arange(len(names))
    width = 0.35
    
    bars1 = ax.bar(x - width/2, times, width, label='Execution Time (ms)', alpha=0.8)
    ax2 = ax.twinx()
    bars2 = ax2.bar(x + width/2, nodes, width, label='Node Count', alpha=0.8, color='orange')
    
    ax.set_xlabel('Input Trace')
    ax.set_ylabel('Execution Time (ms)', color='blue')
    ax2.set_ylabel('Node Count', color='orange')
    ax.set_title('Side-Channel Leakage: Execution Time and Node Count by Input')
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=45, ha='right')
    ax.tick_params(axis='y', labelcolor='blue')
    ax2.tick_params(axis='y', labelcolor='orange')
    
    # Add value labels
    for bar in bars1:
        height = bar.get_height()
        ax.annotate(f'{height:.1f}',
                   xy=(bar.get_x() + bar.get_width() / 2, height),
                   xytext=(0, 3),
                   textcoords="offset points",
                   ha='center', va='bottom', fontsize=8)
    
    for bar in bars2:
        height = bar.get_height()
        ax2.annotate(f'{int(height)}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=8, color='orange')
    
    plt.tight_layout()
    output_path = Path(output_dir) / "execution_leakage.png"
    plt.savefig(output_path, dpi=150)
    print(f"Saved: {output_path}")
    plt.close()


def plot_memory_patterns(traces, output_dir):
    """Plot memory usage patterns during solver execution."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    
    # Select 4 representative traces
    trace_items = list(traces.items())[:4]
    
    for idx, (path, trace) in enumerate(trace_items):
        ax = axes[idx]
        name = Path(path).stem
        memory_snapshots = trace.get("trace", {}).get("memory_snapshots", [])
        
        if memory_snapshots:
            focus_orders = [s["focus_order"] for s in memory_snapshots]
            memory_mb = [s["memory_mb"] for s in memory_snapshots]
            
            ax.plot(focus_orders, memory_mb, marker='o', linewidth=2, markersize=6)
            ax.fill_between(focus_orders, memory_mb, alpha=0.3)
            ax.set_xlabel('Node Focus Order')
            ax.set_ylabel('Memory (MB)')
            ax.set_title(f'{name}')
            ax.grid(True, alpha=0.3)
    
    plt.suptitle('Memory Access Patterns During Branch-and-Bound', fontsize=14)
    plt.tight_layout()
    output_path = Path(output_dir) / "memory_patterns.png"
    plt.savefig(output_path, dpi=150)
    print(f"Saved: {output_path}")
    plt.close()


def plot_traversal_comparison(traces, output_dir):
    """Compare traversal patterns across inputs."""
    fig, ax = plt.subplots(figsize=(12, 8))
    
    for path, trace in traces.items():
        name = Path(path).stem
        traversal = trace.get("trace", {}).get("traversal", [])
        
        if traversal:
            frontier_sizes = [t["frontier_size"] for t in traversal]
            focus_orders = range(len(frontier_sizes))
            
            ax.plot(focus_orders, frontier_sizes, marker='o', label=name, alpha=0.7)
    
    ax.set_xlabel('Node Focus Order')
    ax.set_ylabel('Search Frontier Size')
    ax.set_title('Search Frontier Size Evolution (Traversal Patterns)')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_path = Path(output_dir) / "traversal_patterns.png"
    plt.savefig(output_path, dpi=150)
    print(f"Saved: {output_path}")
    plt.close()


def plot_solution_space(traces, output_dir):
    """Plot solution space showing objective components."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    for path, trace in traces.items():
        name = Path(path).stem
        solution = trace.get("solution", {})
        
        o1 = solution.get("o1_coverage_error", 0)
        o2 = solution.get("o2_interaction_error", 0)
        o3 = solution.get("o3_parsimony", 0)
        configs = solution.get("selected_configurations", [])
        
        ax.scatter(o1, o2, s=o3*50, alpha=0.6, label=name)
        ax.annotate(name, (o1, o2), fontsize=8, alpha=0.7)
    
    ax.set_xlabel('Coverage Error (o1)')
    ax.set_ylabel('Interaction Error (o2)')
    ax.set_title('Solution Space: Coverage vs Interaction Error (size = parsimony)')
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    
    plt.tight_layout()
    output_path = Path(output_dir) / "solution_space.png"
    plt.savefig(output_path, dpi=150)
    print(f"Saved: {output_path}")
    plt.close()


def plot_leakage_heatmap(traces, output_dir):
    """Create heatmap showing distinguishability of inputs."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Extract metrics
    trace_names = []
    metrics = []
    
    for path, trace in traces.items():
        name = Path(path).stem
        trace_names.append(name)
        
        scip_stats = trace.get("scip_statistics", {})
        system_metrics = trace.get("system_metrics", {})
        
        metrics.append([
            scip_stats.get("total_nodes", 0),
            system_metrics.get("total_time_ms", 0),
            system_metrics.get("peak_memory_mb", 0) * 1000,  # Scale for visibility
        ])
    
    metrics = np.array(metrics)
    
    # Normalize for comparison
    if metrics.size > 0:
        metrics_norm = (metrics - metrics.min(axis=0)) / (metrics.max(axis=0) - metrics.min(axis=0) + 1e-6)
        
        im = ax.imshow(metrics_norm, aspect='auto', cmap='viridis')
        ax.set_xticks(range(3))
        ax.set_xticklabels(['Node Count', 'Time (ms)', 'Memory (KB)'])
        ax.set_yticks(range(len(trace_names)))
        ax.set_yticklabels(trace_names)
        ax.set_title('Input Distinguishability Heatmap (normalized metrics)')
        
        plt.colorbar(im, ax=ax, label='Normalized Value')
        plt.tight_layout()
        output_path = Path(output_dir) / "leakage_heatmap.png"
        plt.savefig(output_path, dpi=150)
        print(f"Saved: {output_path}")
        plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("traces", nargs="+", help="Trace JSON files to visualize")
    parser.add_argument("--output", default="visualizations", help="Output directory")
    args = parser.parse_args()
    
    traces = load_traces(args.traces)
    
    if not traces:
        print("No valid traces to visualize")
        return
    
    output_dir = Path(args.output)
    output_dir.mkdir(exist_ok=True)
    
    print(f"Visualizing {len(traces)} traces...")
    plot_execution_times(traces, output_dir)
    plot_memory_patterns(traces, output_dir)
    plot_traversal_comparison(traces, output_dir)
    plot_solution_space(traces, output_dir)
    plot_leakage_heatmap(traces, output_dir)
    
    print(f"\nAll visualizations saved to: {output_dir}")


if __name__ == "__main__":
    main()
