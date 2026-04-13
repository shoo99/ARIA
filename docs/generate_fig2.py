#!/usr/bin/env python3
"""Generate improved Figure 2: ARIA Decision Flow with proper arrows."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(14, 10))
ax.set_xlim(0, 14)
ax.set_ylim(0, 10)
ax.axis('off')

# Title
ax.text(7, 9.7, 'ARIA Decision Flow', fontsize=18, fontweight='bold',
        ha='center', color='#2C3E50')

# Color scheme
C = {
    'input': '#95A5A6', 'dp1': '#2C3E50', 'dp3': '#2980B9', 'de': '#27AE60',
    'dp2': '#E74C3C', 'gsea': '#E74C3C', 'ora': '#27AE60',
    'dp5': '#8E44AD', 'dp4': '#D35400', 'dp7': '#16A085',
    'dp6': '#2C3E50', 'dp8': '#F39C12', 'arrow': '#555555'
}

def box(x, y, w, h, text, color, fontsize=9):
    b = FancyBboxPatch((x-w/2, y-h/2), w, h, boxstyle="round,pad=0.12",
                        facecolor=color, edgecolor='none', alpha=0.92)
    ax.add_patch(b)
    ax.text(x, y, text, fontsize=fontsize, ha='center', va='center',
            color='white', fontweight='bold', linespacing=1.3)

def arrow(x1, y1, x2, y2, color='#555555', lw=1.8):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw,
                                connectionstyle='arc3,rad=0'))

def arrow_curved(x1, y1, x2, y2, color='#555555', lw=1.5, rad=0.2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw,
                                connectionstyle=f'arc3,rad={rad}'))

# ── Row 1: Main pipeline (top) ──────────────────────────────────────────────
y1 = 8.5
box(1.5, y1, 1.8, 0.9, 'Input\nData', C['input'])
box(4.0, y1, 1.8, 0.9, 'DP1\nQC Check', C['dp1'])
box(6.5, y1, 1.8, 0.9, 'DP3\nDesign\nRecognition', C['dp3'])
box(9.5, y1, 1.8, 0.9, 'DE\nAnalysis', C['de'])

arrow(2.4, y1, 3.1, y1)
arrow(4.9, y1, 5.6, y1)
arrow(7.4, y1, 8.6, y1)

# ── Row 2: DP2 diamond + branching ──────────────────────────────────────────
y2 = 6.5
# Diamond for DP2
diamond = plt.Polygon([[9.5, y2+0.65], [10.4, y2], [9.5, y2-0.65], [8.6, y2]],
                       facecolor=C['dp2'], edgecolor='none', alpha=0.92)
ax.add_patch(diamond)
ax.text(9.5, y2, 'DP2\nDEGs?', fontsize=9, ha='center', va='center',
        color='white', fontweight='bold')

arrow(9.5, 8.05, 9.5, 7.15)  # DE → DP2

# ≥100 branch → ORA+GSEA
box(12.2, y2, 1.4, 0.8, 'ORA +\nGSEA', C['ora'])
arrow(10.4, y2, 11.5, y2, C['ora'])
ax.text(11.0, y2+0.25, '≥100', fontsize=9, color=C['ora'], fontweight='bold')

# <50 branch → GSEA Priority
box(6.5, y2, 1.5, 0.8, 'GSEA\nPriority', C['gsea'])
arrow(8.6, y2, 7.25, y2, C['gsea'])
ax.text(7.8, y2+0.25, '<50', fontsize=9, color=C['gsea'], fontweight='bold')

# ── Row 3: Secondary analyses ────────────────────────────────────────────────
y3 = 4.5

box(2.0, y3, 1.6, 0.8, 'DP5\nCross-\nvalidate', C['dp5'])
box(4.5, y3, 1.6, 0.8, 'DP4\nSignature\nDetection', C['dp4'])
box(7.0, y3, 1.6, 0.8, 'DP7\nSensitivity\nAnalysis', C['dp7'])
box(9.5, y3, 1.6, 0.8, 'DP6\nLiterature\nInterpretation', C['dp6'])

# Arrows from ORA+GSEA / GSEA Priority down to secondary analyses
# ORA+GSEA → DP5, DP4, DP7, DP6
arrow(12.2, y2-0.4, 12.2, y3+0.8, C['arrow'])
ax.annotate('', xy=(9.5, y3+0.4), xytext=(12.2, y3+0.8),
            arrowprops=dict(arrowstyle='->', color=C['arrow'], lw=1.3,
                            connectionstyle='arc3,rad=0.15'))
ax.annotate('', xy=(7.0, y3+0.4), xytext=(12.2, y3+0.8),
            arrowprops=dict(arrowstyle='->', color=C['arrow'], lw=1.3,
                            connectionstyle='arc3,rad=0.2'))
ax.annotate('', xy=(4.5, y3+0.4), xytext=(12.2, y3+0.8),
            arrowprops=dict(arrowstyle='->', color=C['arrow'], lw=1.3,
                            connectionstyle='arc3,rad=0.25'))
ax.annotate('', xy=(2.0, y3+0.4), xytext=(12.2, y3+0.8),
            arrowprops=dict(arrowstyle='->', color=C['arrow'], lw=1.3,
                            connectionstyle='arc3,rad=0.3'))

# GSEA Priority → DP7, DP6
arrow_curved(6.5, y2-0.4, 7.0, y3+0.4, C['arrow'], rad=-0.1)
arrow_curved(6.5, y2-0.4, 9.5, y3+0.4, C['arrow'], rad=-0.2)

# Label
ax.text(12.5, y3+1.1, 'Parallel\nanalyses', fontsize=7, ha='center',
        color='#888', fontstyle='italic')

# ── Row 4: DP8 Report ───────────────────────────────────────────────────────
y4 = 2.5
box(5.75, y4, 2.2, 1.0, 'DP8\nReport Generation', C['dp8'], fontsize=10)

# All secondary → DP8
arrow(2.0, y3-0.4, 4.8, y4+0.4, C['arrow'])
arrow(4.5, y3-0.4, 5.2, y4+0.4, C['arrow'])
arrow(7.0, y3-0.4, 6.3, y4+0.4, C['arrow'])
arrow(9.5, y3-0.4, 6.7, y4+0.4, C['arrow'])

# Output
box(5.75, 1.2, 3.5, 0.7, 'HTML Report  •  Excel  •  PNG  •  Decision Log', '#34495e', fontsize=8)
arrow(5.75, y4-0.5, 5.75, 1.55, C['dp8'])

# Legend
legend_y = 1.0
ax.text(0.3, legend_y+0.3, 'Decision types:', fontsize=7, color='#888', fontstyle='italic')
for i, (label, color) in enumerate([('Rule-based', '#2C3E50'), ('LLM-assisted', '#2980B9'),
                                       ('LLM-driven', '#8E44AD'), ('Adaptive', '#E74C3C')]):
    bx = FancyBboxPatch((0.3 + i*1.8, legend_y-0.3), 0.3, 0.3, boxstyle="round,pad=0.05",
                         facecolor=color, edgecolor='none', alpha=0.8)
    ax.add_patch(bx)
    ax.text(0.75 + i*1.8, legend_y-0.15, label, fontsize=6.5, va='center', color='#555')

plt.tight_layout()
plt.savefig('/home/sysoft/ARIA/docs/figures/Fig2_decision_flow.png', dpi=300,
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.savefig('/home/sysoft/ARIA/docs/figures/Fig2_decision_flow.pdf',
            bbox_inches='tight', facecolor='white', edgecolor='none')
plt.close()
print("Figure 2 regenerated with complete arrows")
