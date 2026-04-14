#!/usr/bin/env python3
"""Figure 1 v2: Bold arrows, clear flow."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(10, 14))
ax.set_xlim(0, 10)
ax.set_ylim(0, 16)
ax.axis('off')
fig.patch.set_facecolor('white')

ax.text(5, 15.5, 'ARIA Architecture', fontsize=22, fontweight='bold',
        ha='center', va='center', color='#1a365d')

# ── Helper functions ──
def layer_box(y, h, color, title, subtitle=None):
    box = FancyBboxPatch((0.8, y), 8.4, h, boxstyle="round,pad=0.25",
                          facecolor=color, edgecolor='white', linewidth=2, alpha=0.95)
    ax.add_patch(box)
    ax.text(5, y + h - 0.4, title, fontsize=15, fontweight='bold',
            ha='center', va='center', color='white')
    if subtitle:
        ax.text(5, y + h - 0.85, subtitle, fontsize=9.5,
                ha='center', va='center', color='#e0e0e0')

def big_arrow(y_from, y_to):
    mid_x = 5
    ax.annotate('',
        xy=(mid_x, y_to), xytext=(mid_x, y_from),
        arrowprops=dict(
            arrowstyle='-|>',
            color='#2c3e50',
            lw=3.5,
            mutation_scale=25,
            shrinkA=0, shrinkB=0
        ))

def module_box(x, y, w, h, text, color):
    box = FancyBboxPatch((x - w/2, y - h/2), w, h,
                          boxstyle="round,pad=0.1",
                          facecolor=color, edgecolor='white', linewidth=1.5)
    ax.add_patch(box)
    ax.text(x, y, text, fontsize=8, ha='center', va='center',
            color='white', fontweight='bold', linespacing=1.2)

# ═══════════════════════════════════════════════════════════════
# Layer 1: Reasoning Engine
# ═══════════════════════════════════════════════════════════════
y1 = 13.0
layer_box(y1, 2.0, '#2C3E50',
          'Layer 1: Reasoning Engine (LLM)',
          'QC evaluation  ·  Strategy adaptation  ·  Literature interpretation')

# Arrow 1→2
big_arrow(y1 - 0.05, y1 - 0.55)

# ═══════════════════════════════════════════════════════════════
# Layer 2: Decision Protocol
# ═══════════════════════════════════════════════════════════════
y2 = 10.0
layer_box(y2, 2.5, '#2980B9',
          'Layer 2: Decision Protocol (8 DPs)',
          'Rule-based thresholds  +  LLM contextual reasoning')

# DP boxes inside Layer 2
dps_top = ['DP1: QC', 'DP2: DE Strategy', 'DP3: Design', 'DP4: Signatures']
dps_bot = ['DP5: Validation', 'DP6: Literature', 'DP7: Sensitivity', 'DP8: Report']

for i, dp in enumerate(dps_top):
    x = 1.8 + i * 2.15
    module_box(x, y2 + 1.15, 1.8, 0.5, dp, '#1a5276')

for i, dp in enumerate(dps_bot):
    x = 1.8 + i * 2.15
    module_box(x, y2 + 0.45, 1.8, 0.5, dp, '#1a5276')

# Arrow 2→3
big_arrow(y2 - 0.05, y2 - 0.55)

# ═══════════════════════════════════════════════════════════════
# Layer 3: Execution Engine
# ═══════════════════════════════════════════════════════════════
y3 = 7.5
layer_box(y3, 2.0, '#27AE60',
          'Layer 3: Execution Engine',
          'Docker containers  ·  R/Python script generation  ·  API calls  ·  File I/O')

# Arrow 3→4
big_arrow(y3 - 0.05, y3 - 0.55)

# ═══════════════════════════════════════════════════════════════
# Layer 4: Analysis Modules
# ═══════════════════════════════════════════════════════════════
y4 = 4.2
layer_box(y4, 2.8, '#E74C3C', 'Layer 4: Analysis Modules', None)

modules = [
    ('DESeq2\nedgeR\nlimma', 1.7),
    ('fgsea\nssGSEA', 3.4),
    ('STRING\nPPI', 5.0),
    ('WGCNA', 6.6),
    ('Cell Type\nDeconv', 8.1),
]
for text, x in modules:
    module_box(x, y4 + 0.9, 1.4, 1.2, text, '#c0392b')

# Arrow 4→Output
big_arrow(y4 - 0.05, y4 - 0.55)

# ═══════════════════════════════════════════════════════════════
# Output
# ═══════════════════════════════════════════════════════════════
y5 = 2.2
layer_box(y5, 1.5, '#F39C12', 'Output', None)
ax.text(5, y5 + 0.35, 'HTML Report  ·  Excel  ·  PNG Figures  ·  Decision Log (JSON)',
        fontsize=10, ha='center', va='center', color='white')

# ═══════════════════════════════════════════════════════════════
# Flow labels on arrows
# ═══════════════════════════════════════════════════════════════
arrow_labels = [
    (y1 - 0.3, 'Decisions & Prompts'),
    (y2 - 0.3, 'Analysis Commands'),
    (y3 - 0.3, 'R/Python Scripts'),
    (y4 - 0.3, 'Results & Figures'),
]
for y, label in arrow_labels:
    ax.text(8.2, y, label, fontsize=7.5, ha='left', va='center',
            color='#555', style='italic')

plt.tight_layout()
plt.savefig('/home/sysoft/ARIA/docs/figures/Fig1_architecture.png',
            dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig('/home/sysoft/ARIA/docs/figures/Fig1_architecture.pdf',
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print("Figure 1 v2: bold arrows + flow labels")
