#!/usr/bin/env python3
"""
Visualize BEAST tree using Baltic library
"""

import argparse
import sys
try:
    import baltic as bt
except ImportError:
    print("Error: Baltic library not found. Install with: pip install baltic", file=sys.stderr)
    sys.exit(1)

try:
    import matplotlib.pyplot as plt
    from matplotlib import patches
    from matplotlib.collections import LineCollection
except ImportError:
    print("Error: Matplotlib not found. Install with: pip install matplotlib", file=sys.stderr)
    sys.exit(1)

import numpy as np


def plot_timetree(tree_file, output_prefix, formats=['png', 'svg']):
    """
    Create a time-scaled tree visualization using Baltic
    """
    
    # Load tree
    print(f"Loading tree from {tree_file}")
    try:
        tree = bt.loadNexus(tree_file, absoluteTime=False)
    except:
        try:
            tree = bt.loadNewick(tree_file, absoluteTime=False)
        except Exception as e:
            print(f"Error loading tree: {e}", file=sys.stderr)
            sys.exit(1)
    
    # Set tree to absolute time
    tree.setAbsoluteTime(tree.treeStats()['height'])
    
    # Get tree statistics
    stats = tree.treeStats()
    print(f"Tree statistics:")
    print(f"  Tips: {len(tree.getExternal())}")
    print(f"  Height: {stats['height']:.4f}")
    print(f"  TMRCA: {tree.root.absoluteTime:.4f}")
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, max(8, len(tree.getExternal()) * 0.2)))
    
    # Set up tree traversal
    tree.traverse_tree()
    tree.drawTree()
    
    # Get y positions for tips
    tip_positions = {tip.name: tip.y for tip in tree.getExternal()}
    tip_heights = {tip.name: tip.absoluteTime for tip in tree.getExternal()}
    
    # Draw branches
    x_coords = []
    y_coords = []
    
    for node in tree.Objects:
        if node.parent:
            # Horizontal line (branch)
            x_coords.append([node.parent.absoluteTime, node.absoluteTime])
            y_coords.append([node.y, node.y])
            
            # Vertical line connecting to parent
            if node.parent.children[0] == node:
                # Draw vertical line from parent to all children
                child_ys = [c.y for c in node.parent.children]
                if len(child_ys) > 1:
                    x_coords.append([node.parent.absoluteTime, node.parent.absoluteTime])
                    y_coords.append([min(child_ys), max(child_ys)])
    
    # Plot all branches
    for x, y in zip(x_coords, y_coords):
        ax.plot(x, y, color='black', linewidth=1.5)
    
    # Draw tip labels
    tips = sorted(tree.getExternal(), key=lambda x: x.y)
    for tip in tips:
        ax.text(tip.absoluteTime, tip.y, f'  {tip.name}', 
                va='center', ha='left', fontsize=8)
    
    # Draw nodes with posterior support if available
    for node in tree.getInternal():
        # Check for posterior probability
        posterior = None
        if hasattr(node, 'traits'):
            if 'posterior' in node.traits:
                posterior = node.traits['posterior']
        
        if posterior and posterior >= 0.95:
            ax.plot(node.absoluteTime, node.y, 'o', 
                   color='red', markersize=4, zorder=10)
        elif posterior and posterior >= 0.75:
            ax.plot(node.absoluteTime, node.y, 'o', 
                   color='orange', markersize=3, zorder=10)
    
    # Get time range
    min_time = min([node.absoluteTime for node in tree.Objects])
    max_time = max([node.absoluteTime for node in tree.Objects])
    time_range = max_time - min_time
    
    # Set axis properties
    ax.set_xlabel('Time', fontsize=12, fontweight='bold')
    ax.set_xlim(min_time - time_range * 0.05, max_time + time_range * 0.3)
    ax.set_ylim(-0.5, len(tips) - 0.5)
    ax.set_yticks([])
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # Add title
    ax.set_title('Time-scaled Phylogenetic Tree', fontsize=14, fontweight='bold', pad=20)
    
    # Add legend for posterior support
    legend_elements = [
        patches.Patch(facecolor='red', label='Posterior ≥ 0.95'),
        patches.Patch(facecolor='orange', label='Posterior ≥ 0.75')
    ]
    ax.legend(handles=legend_elements, loc='upper left', frameon=False)
    
    plt.tight_layout()
    
    # Save figures
    for fmt in formats:
        output_file = f"{output_prefix}.{fmt}"
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved figure: {output_file}")
    
    plt.close()


def main():
    parser = argparse.ArgumentParser(description='Visualize BEAST tree with Baltic')
    parser.add_argument('--input', required=True, help='Input tree file (Nexus or Newick)')
    parser.add_argument('--output', required=True, help='Output file prefix')
    parser.add_argument('--format', default='png,svg', help='Output formats (comma-separated)')
    
    args = parser.parse_args()
    
    formats = args.format.split(',')
    
    plot_timetree(args.input, args.output, formats)


if __name__ == '__main__':
    main()
