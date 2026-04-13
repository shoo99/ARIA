#!/usr/bin/env python3
"""Fix Figure 3 (sort by LFC) and Figure 4 (Galaxy reproducibility)."""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

fig_dir = "/home/sysoft/ARIA/docs/figures"

###############################################################################
# Figure 3: Airway Benchmark — Fix Panel B gene order
###############################################################################
def fig3():
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: DEG gain (unchanged)
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

    # Panel B: Known targets — SORTED by LFC (ascending for horizontal bar)
    genes = ['CRISPLD2', 'DUSP1', 'PER1', 'TSC22D3', 'FKBP5', 'KLF15', 'ZBTB16']
    lfc =   [2.63,       2.94,    3.19,   3.19,      4.04,    4.46,    7.35]

    # Sort by LFC
    sorted_pairs = sorted(zip(lfc, genes))
    lfc_sorted = [p[0] for p in sorted_pairs]
    genes_sorted = [p[1] for p in sorted_pairs]

    colors = ['#E74C3C'] * len(genes)
    axes[1].barh(genes_sorted, lfc_sorted, color=colors, edgecolor='white', alpha=0.85)
    axes[1].set_xlabel('log2 Fold Change')
    axes[1].set_title('B. Known Dex Targets (7/7 detected)', fontweight='bold')
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)

    for i, v in enumerate(lfc_sorted):
        axes[1].text(v + 0.1, i, f'{v:.2f}', va='center', fontsize=9)

    plt.tight_layout()
    plt.savefig(f'{fig_dir}/Fig3_airway_benchmark.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(f'{fig_dir}/Fig3_airway_benchmark.pdf', bbox_inches='tight', facecolor='white')
    plt.close()
    print("Fig 3 fixed: genes sorted by LFC")

###############################################################################
# Figure 4: Comparison Table — Fix Galaxy reproducibility
###############################################################################
def fig4():
    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.axis('off')

    features = ['Execution automation', 'Result-based adaptation', 'Design recognition',
                'Cross-method validation', 'Biological interpretation',
                'Report generation', 'Decision transparency', 'Reproducibility']
    tools = ['nf-core', 'Galaxy', 'iDEP', 'ARIA']

    # 1 = full, 0.5 = partial, 0 = none
    # Fixed: Galaxy reproducibility = 1 (Docker/Conda support)
    data = [
        [1, 1, 1, 1],       # Execution
        [0, 0, 0.5, 1],     # Adaptation (iDEP has some interactive choices)
        [0, 0, 0, 1],       # Design recognition (iDEP doesn't auto-detect)
        [0, 0, 0, 1],       # Cross-method
        [0, 0, 0.5, 1],     # Interpretation (iDEP has some annotation)
        [1, 1, 1, 1],       # Report
        [0, 0, 0, 1],       # Transparency
        [1, 1, 0.5, 1],     # Reproducibility (Galaxy has Docker; iDEP web-only)
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
        cell_text.append(['\u2713' if v == 1 else '\u25B3' if v == 0.5 else '\u2717' for v in row])

    table = ax.table(cellText=cell_text, rowLabels=features, colLabels=tools,
                      cellColours=cell_colors, loc='center', cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 1.8)

    for (i, j), cell in table.get_celld().items():
        if i == 0:  # Header row
            cell.set_text_props(fontweight='bold', color='white')
            cell.set_facecolor('#2C3E50')
        if j == -1:  # Row labels
            cell.set_text_props(fontsize=9.5, ha='left')
            cell.set_facecolor('#f8f8f8')
            cell.set_edgecolor('#ddd')
        if i > 0 and j >= 0:
            cell.set_text_props(color='white', fontweight='bold', fontsize=13)
            cell.set_edgecolor('#ffffff')

    ax.set_title('Feature Comparison: ARIA vs Existing Tools', fontsize=14,
                  fontweight='bold', pad=20)

    # Legend
    ax.text(0.15, -0.02, '\u2713 = Supported    \u25B3 = Partial    \u2717 = Not supported',
            transform=ax.transAxes, fontsize=9, color='#666')

    plt.tight_layout()
    plt.savefig(f'{fig_dir}/Fig4_comparison.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(f'{fig_dir}/Fig4_comparison.pdf', bbox_inches='tight', facecolor='white')
    plt.close()
    print("Fig 4 fixed: Galaxy reproducibility=✓, iDEP design=✗, legend added")

fig3()
fig4()
print("All figures updated")
