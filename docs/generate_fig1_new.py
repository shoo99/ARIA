#!/usr/bin/env python3
"""Generate clean Figure 1: ARIA Architecture (no embedded caption)."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

fig, ax = plt.subplots(figsize=(10, 13))
ax.set_xlim(0, 10)
ax.set_ylim(0, 15)
ax.axis('off')

# Title
ax.text(5, 14.5, 'ARIA Architecture', fontsize=20, fontweight='bold',
        ha='center', va='center', color='#1a365d')

# Colors
layer_colors = {
    1: '#2C3E50',  # Reasoning Engine
    2: '#2980B9',  # Decision Protocol
    3: '#27AE60',  # Execution Engine
    4: '#E74C3C',  # Analysis Modules
}

def draw_layer(y_top, height, color, title, subtitle, items=None):
    box = FancyBboxPatch((0.5, y_top - height), 9, height,
                          boxstyle="round,pad=0.2", facecolor=color,
                          edgecolor='none', alpha=0.92)
    ax.add_patch(box)
    ax.text(5, y_top - 0.35, title, fontsize=14, fontweight='bold',
            ha='center', va='center', color='white')
    if subtitle:
        ax.text(5, y_top - 0.75, subtitle, fontsize=9,
                ha='center', va='center', color='#ddd', style='italic')
    if items:
        for i, item in enumerate(items):
            x = 1.0 + (i % 4) * 2.2
            y = y_top - height + 0.55 if i < 4 else y_top - height + 0.25
            ax.text(x, y, item, fontsize=7.5, ha='center', va='center',
                    color='white', fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.12',
                              facecolor='black', alpha=0.25, edgecolor='none'))

def draw_arrow(y_from, y_to):
    ax.annotate('', xy=(5, y_to + 0.1), xytext=(5, y_from - 0.1),
                arrowprops=dict(arrowstyle='->', color='#7f8c8d', lw=2.5))

# Layer 1: Reasoning Engine
draw_layer(13.5, 2.0, layer_colors[1],
           'Layer 1: Reasoning Engine (LLM)',
           'QC evaluation  ·  Strategy adaptation  ·  Literature interpretation')

draw_arrow(11.5, 11.1)

# Layer 2: Decision Protocol
draw_layer(11.0, 2.5, layer_colors[2],
           'Layer 2: Decision Protocol (8 DPs)',
           'Rule-based thresholds + LLM contextual reasoning',
           items=['DP1: QC', 'DP2: DE Strategy', 'DP3: Design', 'DP4: Signatures',
                  'DP5: Validation', 'DP6: Literature', 'DP7: Sensitivity', 'DP8: Report'])

draw_arrow(8.5, 8.1)

# Layer 3: Execution Engine
draw_layer(8.0, 1.8, layer_colors[3],
           'Layer 3: Execution Engine',
           'Docker containers  ·  R/Python script generation  ·  API calls  ·  File I/O')

draw_arrow(6.2, 5.8)

# Layer 4: Analysis Modules
draw_layer(5.7, 2.5, layer_colors[4],
           'Layer 4: Analysis Modules',
           None)

modules = [
    ('DESeq2\nedgeR\nlimma', 1.2),
    ('fgsea\nssGSEA', 3.0),
    ('STRING\nPPI', 4.8),
    ('WGCNA', 6.5),
    ('Cell Type\nDeconv', 8.0),
    ('Report\nGenerator', 9.3),
]
for text, x in modules:
    box = FancyBboxPatch((x - 0.65, 3.55), 1.3, 1.1,
                          boxstyle="round,pad=0.1",
                          facecolor='#c0392b', edgecolor='none', alpha=0.85)
    ax.add_patch(box)
    ax.text(x, 4.1, text, fontsize=7.5, ha='center', va='center',
            color='white', fontweight='bold', linespacing=1.2)

draw_arrow(3.2, 2.8)

# Output box
out_box = FancyBboxPatch((1.0, 1.2), 8, 1.4,
                          boxstyle="round,pad=0.2",
                          facecolor='#F39C12', edgecolor='none', alpha=0.92)
ax.add_patch(out_box)
ax.text(5, 2.2, 'Output', fontsize=14, fontweight='bold',
        ha='center', va='center', color='white')
ax.text(5, 1.7, 'HTML Report  ·  Excel Data  ·  PNG Figures  ·  Decision Log (JSON)',
        fontsize=9.5, ha='center', va='center', color='white')

# Layer labels on left side
labels = [
    (12.5, 'Reasoning', '#2C3E50'),
    (9.75, 'Decisions', '#2980B9'),
    (7.1, 'Execution', '#27AE60'),
    (4.45, 'Analysis', '#E74C3C'),
    (1.9, 'Output', '#F39C12'),
]

plt.tight_layout()
plt.savefig('/home/sysoft/ARIA/docs/figures/Fig1_architecture.png',
            dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig('/home/sysoft/ARIA/docs/figures/Fig1_architecture.pdf',
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print("Figure 1 regenerated (clean, no embedded caption)")
