#!/usr/bin/env python3
"""Generate publication figures for the ARIA paper."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np

fig_dir = "/home/sysoft/ARIA/docs/figures"
import os
os.makedirs(fig_dir, exist_ok=True)

###############################################################################
# Figure 1: ARIA Architecture
###############################################################################
def fig1_architecture():
    fig, ax = plt.subplots(1, 1, figsize=(10, 12))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 14)
    ax.axis('off')

    colors = {
        'layer1': '#2C3E50',
        'layer2': '#2980B9',
        'layer3': '#27AE60',
        'layer4': '#E74C3C',
        'arrow': '#7F8C8D',
        'text': 'white'
    }

    # Title
    ax.text(5, 13.5, 'ARIA Architecture', fontsize=18, fontweight='bold',
            ha='center', va='center', color='#2C3E50')

    # Layer 1: Reasoning Engine
    box1 = FancyBboxPatch((0.5, 11), 9, 2, boxstyle="round,pad=0.2",
                           facecolor=colors['layer1'], edgecolor='none', alpha=0.9)
    ax.add_patch(box1)
    ax.text(5, 12.3, 'Layer 1: Reasoning Engine (LLM)', fontsize=13,
            fontweight='bold', ha='center', va='center', color='white')
    ax.text(5, 11.6, 'QC evaluation • Strategy adaptation • Literature interpretation',
            fontsize=9, ha='center', va='center', color='#BDC3C7')

    # Arrow
    ax.annotate('', xy=(5, 11), xytext=(5, 10.6),
                arrowprops=dict(arrowstyle='->', color=colors['arrow'], lw=2))

    # Layer 2: Decision Protocol
    box2 = FancyBboxPatch((0.5, 8.2), 9, 2.2, boxstyle="round,pad=0.2",
                           facecolor=colors['layer2'], edgecolor='none', alpha=0.9)
    ax.add_patch(box2)
    ax.text(5, 9.7, 'Layer 2: Decision Protocol (8 DPs)', fontsize=13,
            fontweight='bold', ha='center', va='center', color='white')

    dps = ['DP1:QC', 'DP2:DE Strategy', 'DP3:Design', 'DP4:Signatures',
           'DP5:Validation', 'DP6:Literature', 'DP7:Sensitivity', 'DP8:Report']
    for i, dp in enumerate(dps):
        x = 1.2 + (i % 4) * 2.2
        y = 9.0 if i < 4 else 8.5
        ax.text(x, y, dp, fontsize=7, ha='center', va='center', color='white',
                bbox=dict(boxstyle='round,pad=0.15', facecolor='#1A5276', edgecolor='none'))

    # Arrow
    ax.annotate('', xy=(5, 8.2), xytext=(5, 7.8),
                arrowprops=dict(arrowstyle='->', color=colors['arrow'], lw=2))

    # Layer 3: Execution Engine
    box3 = FancyBboxPatch((0.5, 5.8), 9, 1.8, boxstyle="round,pad=0.2",
                           facecolor=colors['layer3'], edgecolor='none', alpha=0.9)
    ax.add_patch(box3)
    ax.text(5, 7.0, 'Layer 3: Execution Engine', fontsize=13,
            fontweight='bold', ha='center', va='center', color='white')
    ax.text(5, 6.3, 'Docker containers • R/Python scripts • API calls • File I/O',
            fontsize=9, ha='center', va='center', color='#D5F5E3')

    # Arrow
    ax.annotate('', xy=(5, 5.8), xytext=(5, 5.4),
                arrowprops=dict(arrowstyle='->', color=colors['arrow'], lw=2))

    # Layer 4: Analysis Modules
    box4 = FancyBboxPatch((0.5, 3), 9, 2.2, boxstyle="round,pad=0.2",
                           facecolor=colors['layer4'], edgecolor='none', alpha=0.9)
    ax.add_patch(box4)
    ax.text(5, 4.6, 'Layer 4: Analysis Modules', fontsize=13,
            fontweight='bold', ha='center', va='center', color='white')

    modules = ['DESeq2\nedgeR\nlimma', 'fgsea\nssGSEA', 'STRING\nPPI', 'WGCNA', 'Cell Type\nDeconv', 'Report\nGen']
    for i, mod in enumerate(modules):
        x = 1.2 + i * 1.5
        ax.text(x, 3.5, mod, fontsize=7, ha='center', va='center', color='white',
                bbox=dict(boxstyle='round,pad=0.2', facecolor='#C0392B', edgecolor='none'))

    # Output
    box_out = FancyBboxPatch((1.5, 1), 7, 1.5, boxstyle="round,pad=0.2",
                              facecolor='#F39C12', edgecolor='none', alpha=0.9)
    ax.add_patch(box_out)
    ax.annotate('', xy=(5, 3), xytext=(5, 2.6),
                arrowprops=dict(arrowstyle='->', color=colors['arrow'], lw=2))
    ax.text(5, 2.0, 'Output', fontsize=13, fontweight='bold',
            ha='center', va='center', color='white')
    ax.text(5, 1.4, 'HTML Report • Excel • PNG Figures • Decision Log',
            fontsize=9, ha='center', va='center', color='white')

    plt.tight_layout()
    plt.savefig(f'{fig_dir}/Fig1_architecture.png', dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.savefig(f'{fig_dir}/Fig1_architecture.pdf', bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()
    print("Fig 1: Architecture saved")

###############################################################################
# Figure 2: Decision Flow
###############################################################################
def fig2_decision_flow():
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.axis('off')

    ax.text(6, 7.6, 'ARIA Decision Flow', fontsize=16, fontweight='bold',
            ha='center', color='#2C3E50')

    # Flow boxes
    steps = [
        (1.5, 6.5, 'Input\nData', '#95A5A6'),
        (4, 6.5, 'DP1\nQC Check', '#2C3E50'),
        (6.5, 6.5, 'DP3\nDesign\nRecognition', '#2980B9'),
        (9, 6.5, 'DE\nAnalysis', '#27AE60'),
    ]

    for x, y, txt, col in steps:
        box = FancyBboxPatch((x-0.8, y-0.5), 1.6, 1.0, boxstyle="round,pad=0.1",
                              facecolor=col, edgecolor='none', alpha=0.9)
        ax.add_patch(box)
        ax.text(x, y, txt, fontsize=8, ha='center', va='center', color='white', fontweight='bold')

    # Arrows between top row
    for x1, x2 in [(2.3, 3.2), (4.8, 5.7), (7.3, 8.2)]:
        ax.annotate('', xy=(x2, 6.5), xytext=(x1, 6.5),
                    arrowprops=dict(arrowstyle='->', color='#7F8C8D', lw=1.5))

    # DP2 Decision diamond
    diamond_x, diamond_y = 9, 4.5
    diamond = plt.Polygon([[diamond_x, diamond_y+0.7], [diamond_x+1, diamond_y],
                            [diamond_x, diamond_y-0.7], [diamond_x-1, diamond_y]],
                           facecolor='#E74C3C', edgecolor='none', alpha=0.9)
    ax.add_patch(diamond)
    ax.text(diamond_x, diamond_y, 'DP2\nDEGs?', fontsize=8, ha='center', va='center',
            color='white', fontweight='bold')
    ax.annotate('', xy=(9, 5.15), xytext=(9, 6.0),
                arrowprops=dict(arrowstyle='->', color='#7F8C8D', lw=1.5))

    # Many DEGs path
    ax.text(10.5, 4.5, '≥100', fontsize=8, color='#27AE60', fontweight='bold')
    box_ora = FancyBboxPatch((10.5, 3.5), 1.3, 0.8, boxstyle="round,pad=0.1",
                              facecolor='#27AE60', edgecolor='none')
    ax.add_patch(box_ora)
    ax.text(11.15, 3.9, 'ORA +\nGSEA', fontsize=7, ha='center', va='center', color='white')
    ax.annotate('', xy=(10.5, 4.3), xytext=(10.0, 4.5),
                arrowprops=dict(arrowstyle='->', color='#27AE60', lw=1.5))

    # Few DEGs path
    ax.text(7, 4.5, '<50', fontsize=8, color='#E74C3C', fontweight='bold')
    box_gsea = FancyBboxPatch((5.5, 3.5), 1.5, 0.8, boxstyle="round,pad=0.1",
                               facecolor='#E74C3C', edgecolor='none')
    ax.add_patch(box_gsea)
    ax.text(6.25, 3.9, 'GSEA\nPriority', fontsize=7, ha='center', va='center', color='white')
    ax.annotate('', xy=(7.5, 4.3), xytext=(8.0, 4.5),
                arrowprops=dict(arrowstyle='->', color='#E74C3C', lw=1.5))

    # Additional analyses
    add_steps = [
        (2, 3.5, 'DP5\nCross-validate', '#8E44AD'),
        (4, 3.5, 'DP4\nSignatures', '#D35400'),
        (2, 2, 'DP7\nSensitivity', '#16A085'),
        (4, 2, 'DP6\nLiterature', '#2C3E50'),
        (6.25, 2, 'DP8\nReport', '#F39C12'),
    ]

    for x, y, txt, col in add_steps:
        box = FancyBboxPatch((x-0.7, y-0.4), 1.4, 0.8, boxstyle="round,pad=0.1",
                              facecolor=col, edgecolor='none', alpha=0.85)
        ax.add_patch(box)
        ax.text(x, y, txt, fontsize=7, ha='center', va='center', color='white', fontweight='bold')

    # Arrows to additional steps
    ax.annotate('', xy=(6.25, 2.8), xytext=(6.25, 3.5),
                arrowprops=dict(arrowstyle='->', color='#7F8C8D', lw=1))

    plt.tight_layout()
    plt.savefig(f'{fig_dir}/Fig2_decision_flow.png', dpi=300, bbox_inches='tight',
                facecolor='white')
    plt.savefig(f'{fig_dir}/Fig2_decision_flow.pdf', bbox_inches='tight',
                facecolor='white')
    plt.close()
    print("Fig 2: Decision flow saved")

###############################################################################
# Figure 3: Airway Benchmark — Paired vs Unpaired
###############################################################################
def fig3_airway_benchmark():
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: DEG gain
    cutoffs = ['|LFC|>1.0', '|LFC|>0.5', 'No LFC']
    unpaired = [785, 1860, 2773]
    paired = [951, 2426, 4081]
    gain_pct = [21, 30, 47]

    x = np.arange(len(cutoffs))
    w = 0.35

    axes[0].bar(x - w/2, unpaired, w, label='Unpaired (naive)', color='#95A5A6', edgecolor='white')
    axes[0].bar(x + w/2, paired, w, label='Paired (ARIA)', color='#E74C3C', edgecolor='white')

    for i, (u, p, g) in enumerate(zip(unpaired, paired, gain_pct)):
        axes[0].text(i + w/2, p + 50, f'+{g}%', ha='center', fontsize=10,
                     fontweight='bold', color='#E74C3C')

    axes[0].set_xticks(x)
    axes[0].set_xticklabels(cutoffs)
    axes[0].set_ylabel('Number of DEGs')
    axes[0].set_title('A. DEG Detection: Paired vs Unpaired', fontweight='bold')
    axes[0].legend()
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)

    # Panel B: Known targets
    genes = ['ZBTB16', 'FKBP5', 'KLF15', 'PER1', 'TSC22D3', 'DUSP1', 'CRISPLD2']
    lfc = [7.35, 4.04, 4.46, 3.19, 3.19, 2.94, 2.63]
    colors = ['#E74C3C'] * len(genes)

    axes[1].barh(genes, lfc, color=colors, edgecolor='white', alpha=0.85)
    axes[1].set_xlabel('log2 Fold Change')
    axes[1].set_title('B. Known Dex Targets (all detected)', fontweight='bold')
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)

    for i, v in enumerate(lfc):
        axes[1].text(v + 0.1, i, f'{v:.2f}', va='center', fontsize=9)

    plt.tight_layout()
    plt.savefig(f'{fig_dir}/Fig3_airway_benchmark.png', dpi=300, bbox_inches='tight',
                facecolor='white')
    plt.savefig(f'{fig_dir}/Fig3_airway_benchmark.pdf', bbox_inches='tight',
                facecolor='white')
    plt.close()
    print("Fig 3: Airway benchmark saved")

###############################################################################
# Figure 4: Comparison Table
###############################################################################
def fig4_comparison():
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.axis('off')

    features = ['Execution automation', 'Result-based adaptation', 'Design recognition',
                'Cross-method validation', 'Biological interpretation',
                'Report generation', 'Decision transparency', 'Reproducibility']
    tools = ['nf-core', 'Galaxy', 'iDEP', 'ARIA']

    data = [
        [1, 1, 1, 1],  # Execution
        [0, 0, 0.5, 1],  # Adaptation
        [0, 0, 0.5, 1],  # Design
        [0, 0, 0, 1],  # Cross-method
        [0, 0, 0.5, 1],  # Interpretation
        [1, 1, 1, 1],  # Report
        [0, 0, 0, 1],  # Transparency
        [1, 0.5, 0.5, 1],  # Reproducibility
    ]

    cell_colors = []
    for row in data:
        row_colors = []
        for val in row:
            if val == 1:
                row_colors.append('#27AE60')
            elif val == 0.5:
                row_colors.append('#F39C12')
            else:
                row_colors.append('#E74C3C')
        cell_colors.append(row_colors)

    cell_text = []
    for row in data:
        cell_text.append(['✓' if v == 1 else '△' if v == 0.5 else '✗' for v in row])

    table = ax.table(cellText=cell_text, rowLabels=features, colLabels=tools,
                      cellColours=cell_colors, loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.8)

    for (i, j), cell in table.get_celld().items():
        if i == 0:
            cell.set_text_props(fontweight='bold', color='white')
            cell.set_facecolor('#2C3E50')
        if j == -1:
            cell.set_text_props(fontsize=9)
        if i > 0 and j >= 0:
            cell.set_text_props(color='white', fontweight='bold', fontsize=12)

    ax.set_title('Feature Comparison: ARIA vs Existing Tools', fontsize=14,
                  fontweight='bold', pad=20)

    plt.tight_layout()
    plt.savefig(f'{fig_dir}/Fig4_comparison.png', dpi=300, bbox_inches='tight',
                facecolor='white')
    plt.savefig(f'{fig_dir}/Fig4_comparison.pdf', bbox_inches='tight',
                facecolor='white')
    plt.close()
    print("Fig 4: Comparison table saved")

# Generate all figures
fig1_architecture()
fig2_decision_flow()
fig3_airway_benchmark()
fig4_comparison()
print(f"\nAll figures saved to {fig_dir}/")
