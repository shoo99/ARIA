#!/usr/bin/env python3
"""Generate bioRxiv-ready single PDF with embedded figures, line numbers, page numbers."""

import base64, os, re

DOCS = "/home/sysoft/ARIA/docs"
FIGS = f"{DOCS}/figures"
OUT = f"{DOCS}/ARIA_bioRxiv.html"

def img_b64(path, width="85%"):
    if not os.path.exists(path):
        return f"<p>[Figure not found: {path}]</p>"
    with open(path, "rb") as f:
        e = base64.b64encode(f.read()).decode()
    return f'<div class="figure"><img src="data:image/png;base64,{e}" style="width:{width};"></div>'

# Read manuscript markdown
with open(f"{DOCS}/manuscript_full.md") as f:
    md = f.read()

# Convert markdown tables to HTML tables
def md_table_to_html(match):
    lines = match.group(0).strip().split('\n')
    headers = [c.strip() for c in lines[0].split('|')[1:-1]]
    rows = []
    for line in lines[2:]:
        cells = [c.strip() for c in line.split('|')[1:-1]]
        rows.append(cells)

    html = '<table><thead><tr>'
    for h in headers:
        html += f'<th>{h}</th>'
    html += '</tr></thead><tbody>'
    for row in rows:
        html += '<tr>'
        for c in row:
            c = c.replace('**', '')  # remove bold markers
            html += f'<td>{c}</td>'
        html += '</tr>'
    html += '</tbody></table>'
    return html

# Process markdown to HTML manually (simple conversion)
def md_to_html(text):
    lines = text.split('\n')
    html_lines = []
    in_table = False
    table_lines = []
    line_num = 0

    for line in lines:
        # Skip front matter
        if line.strip() == '---':
            html_lines.append('<hr>')
            continue

        # Tables
        if '|' in line and line.strip().startswith('|'):
            table_lines.append(line)
            in_table = True
            continue
        elif in_table:
            html_lines.append(md_table_to_html(type('', (), {'group': lambda s, x=0: '\n'.join(table_lines)})()))
            table_lines = []
            in_table = False

        # Headers
        if line.startswith('# '):
            html_lines.append(f'<h1>{line[2:]}</h1>')
        elif line.startswith('## '):
            html_lines.append(f'<h2>{line[3:]}</h2>')
        elif line.startswith('### '):
            html_lines.append(f'<h3>{line[4:]}</h3>')
        elif line.startswith('### '):
            html_lines.append(f'<h4>{line[5:]}</h4>')
        elif line.strip().startswith('- **'):
            html_lines.append(f'<li>{line.strip()[2:]}</li>')
        elif line.strip().startswith('- '):
            html_lines.append(f'<li>{line.strip()[2:]}</li>')
        elif line.strip().startswith(('1.', '2.', '3.', '4.', '5.')):
            html_lines.append(f'<li>{line.strip()[3:]}</li>')
        elif line.strip() == '':
            html_lines.append('<br>')
        else:
            line_num += 1
            # Bold
            line = re.sub(r'\*\*(.+?)\*\*', r'<strong>\1</strong>', line)
            # Italic
            line = re.sub(r'\*(.+?)\*', r'<em>\1</em>', line)
            html_lines.append(f'<p><span class="linenum">{line_num}</span>{line}</p>')

    if in_table and table_lines:
        html_lines.append(md_table_to_html(type('', (), {'group': lambda s, x=0: '\n'.join(table_lines)})()))

    return '\n'.join(html_lines)

body = md_to_html(md)

html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<style>
@page {{
    size: letter;
    margin: 1in;
    @bottom-center {{
        content: counter(page);
        font-size: 10pt;
        color: #666;
    }}
}}
body {{
    font-family: "Times New Roman", Times, serif;
    font-size: 11pt;
    line-height: 2.0;
    color: #000;
    max-width: 7in;
    margin: 0 auto;
    counter-reset: line-number;
}}
h1 {{
    font-size: 16pt;
    text-align: center;
    margin: 0 0 5px 0;
    line-height: 1.3;
}}
h2 {{
    font-size: 13pt;
    margin: 20pt 0 8pt 0;
    page-break-after: avoid;
}}
h3 {{
    font-size: 11pt;
    margin: 15pt 0 5pt 0;
    page-break-after: avoid;
}}
p {{
    margin: 0 0 6pt 0;
    text-align: justify;
    position: relative;
    padding-left: 35px;
}}
.linenum {{
    position: absolute;
    left: 0;
    color: #999;
    font-size: 8pt;
    width: 30px;
    text-align: right;
    font-family: Arial, sans-serif;
}}
.author {{
    text-align: center;
    font-size: 11pt;
    margin: 5px 0;
}}
.affiliation {{
    text-align: center;
    font-size: 10pt;
    color: #555;
    margin: 2px 0;
}}
.email {{
    text-align: center;
    font-size: 10pt;
    margin: 5px 0 20px 0;
}}
.keywords {{
    font-size: 10pt;
    margin: 10px 0 20px 0;
    padding-left: 35px;
}}
table {{
    width: 100%;
    border-collapse: collapse;
    margin: 10pt 0;
    font-size: 9pt;
    line-height: 1.4;
    page-break-inside: avoid;
}}
th {{
    background: #f0f0f0;
    border-top: 2px solid #000;
    border-bottom: 1px solid #000;
    padding: 5px 6px;
    text-align: left;
    font-weight: bold;
}}
td {{
    padding: 4px 6px;
    border-bottom: 1px solid #ddd;
}}
tr:last-child td {{
    border-bottom: 2px solid #000;
}}
li {{
    margin: 3pt 0;
    padding-left: 35px;
}}
hr {{
    border: none;
    border-top: 1px solid #ccc;
    margin: 15pt 0;
}}
.figure {{
    text-align: center;
    margin: 15pt 0;
    page-break-inside: avoid;
}}
.figure img {{
    max-width: 100%;
    border: 1px solid #ddd;
}}
.fig-caption {{
    font-size: 9pt;
    text-align: left;
    margin: 5px 0 15px 0;
    padding-left: 35px;
    line-height: 1.4;
}}
.abstract-label {{
    font-weight: bold;
    font-size: 11pt;
}}
</style>
</head>
<body>

{body}

<h2>Figures</h2>

<div class="figure">
{img_b64(f"{FIGS}/Fig1_architecture.png", "80%")}
</div>
<p class="fig-caption"><strong>Figure 1.</strong> ARIA four-layer architecture. Layer 1: LLM-based Reasoning Engine evaluates results and makes decisions. Layer 2: Decision Protocol with 8 formalized Decision Points. Layer 3: Execution Engine manages Docker containers and script execution. Layer 4: Six analysis modules (DESeq2, fgsea, STRING DB, WGCNA, cell type deconvolution, ssGSEA).</p>

<div class="figure">
{img_b64(f"{FIGS}/Fig2_decision_flow.png", "90%")}
</div>
<p class="fig-caption"><strong>Figure 2.</strong> ARIA decision flow. After QC assessment (DP1) and design recognition (DP3), DE analysis results are evaluated at DP2, which determines the downstream strategy: standard (ORA + GSEA) for abundant DEGs, or GSEA-priority for datasets with few DEGs. Cross-method validation (DP5), sensitivity analysis (DP7), and report generation (DP8) follow.</p>

<div class="figure">
{img_b64(f"{FIGS}/Fig3_airway_benchmark.png", "95%")}
</div>
<p class="fig-caption"><strong>Figure 3.</strong> Airway benchmark results. (A) DEG detection comparison between unpaired (naive) and paired (ARIA) models. The paired model detected 21–47% more DEGs across all cutoff thresholds. (B) All 7 known dexamethasone target genes were recovered with strong fold changes (LFC 2.6–7.4).</p>

<div class="figure">
{img_b64(f"{FIGS}/Fig4_comparison.png", "80%")}
</div>
<p class="fig-caption"><strong>Figure 4.</strong> Feature comparison between ARIA and existing tools. Green (checkmark) = fully supported; orange (triangle) = partially supported; red (X) = not supported. ARIA's unique features include result-based adaptation, design recognition, cross-method validation, and decision transparency.</p>

</body>
</html>"""

with open(OUT, "w") as f:
    f.write(html)

# Convert to PDF
import subprocess
result = subprocess.run(["weasyprint", OUT, f"{DOCS}/ARIA_bioRxiv.pdf"], capture_output=True, text=True)
if result.returncode == 0:
    size = os.path.getsize(f"{DOCS}/ARIA_bioRxiv.pdf") / 1024
    print(f"bioRxiv PDF generated: {size:.0f} KB")
else:
    print(f"Error: {result.stderr[-200:]}")
