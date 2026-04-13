#!/usr/bin/env python3
"""Generate bioRxiv-ready PDF using pandoc for proper markdown conversion."""

import base64, os, subprocess

DOCS = "/home/sysoft/ARIA/docs"
FIGS = f"{DOCS}/figures"
PANDOC = "/home/sysoft/bin/pandoc"

def img_b64(path, width="85%"):
    if not os.path.exists(path):
        return ""
    with open(path, "rb") as f:
        e = base64.b64encode(f.read()).decode()
    return f'<div style="text-align:center;margin:20px 0;page-break-inside:avoid;"><img src="data:image/png;base64,{e}" style="width:{width};border:1px solid #ddd;"></div>'

# Step 1: Convert manuscript markdown to HTML body using pandoc
result = subprocess.run(
    [PANDOC, f"{DOCS}/manuscript_full.md", "-t", "html", "--no-highlight"],
    capture_output=True, text=True
)
body_html = result.stdout

# Step 2: Build figure section
figures_html = f"""
<h2>Figures</h2>

{img_b64(f"{FIGS}/Fig1_architecture.png", "85%")}
<p class="fig-cap"><strong>Figure 1.</strong> ARIA four-layer architecture. Layer 1: LLM-based Reasoning Engine evaluates results and makes decisions. Layer 2: Decision Protocol with 8 formalized Decision Points combining rule-based thresholds with LLM reasoning. Layer 3: Execution Engine manages Docker containers and R/Python script execution. Layer 4: Six analysis modules (DESeq2, fgsea, STRING DB, WGCNA, cell type deconvolution, ssGSEA).</p>

{img_b64(f"{FIGS}/Fig2_decision_flow.png", "90%")}
<p class="fig-cap"><strong>Figure 2.</strong> ARIA decision flow. After QC assessment (DP1) and design recognition (DP3), DE analysis results are evaluated at DP2, which determines the downstream strategy based on DEG count. Cross-method validation (DP5), sensitivity analysis (DP7), and report generation (DP8) follow adaptively.</p>

{img_b64(f"{FIGS}/Fig3_airway_benchmark.png", "95%")}
<p class="fig-cap"><strong>Figure 3.</strong> Airway benchmark results. (A) DEG detection comparison between unpaired (naive) and paired (ARIA) models. The paired model detected 21–47% more DEGs across all cutoff thresholds. (B) All 7 known dexamethasone target genes were recovered with strong fold changes (LFC 2.6–7.4).</p>

{img_b64(f"{FIGS}/Fig4_comparison.png", "80%")}
<p class="fig-cap"><strong>Figure 4.</strong> Feature comparison between ARIA and existing tools. Green (checkmark) = fully supported; orange (triangle) = partially supported; red (X) = not supported. ARIA uniquely provides result-based adaptation, experimental design recognition, cross-method validation, and decision transparency.</p>
"""

# Step 3: Wrap in full HTML with bioRxiv styling
full_html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<style>
@page {{
    size: letter;
    margin: 1in 1in 1in 1.2in;
    @bottom-center {{ content: counter(page); font-size: 10pt; color: #666; }}
}}
body {{
    font-family: "Times New Roman", Times, "DejaVu Serif", serif;
    font-size: 11pt;
    line-height: 2.0;
    color: #000;
    counter-reset: line-counter;
}}
h1 {{ font-size: 16pt; text-align: center; margin: 0 0 8px 0; line-height: 1.3; }}
h2 {{ font-size: 13pt; margin: 24pt 0 10pt 0; page-break-after: avoid; border-bottom: 1px solid #ccc; padding-bottom: 4px; }}
h3 {{ font-size: 11.5pt; margin: 18pt 0 6pt 0; page-break-after: avoid; }}
h4 {{ font-size: 11pt; font-style: italic; margin: 14pt 0 4pt 0; page-break-after: avoid; }}
p {{ margin: 0 0 8pt 0; text-align: justify; }}
ul, ol {{ margin: 6pt 0 6pt 20pt; }}
li {{ margin: 3pt 0; }}
hr {{ border: none; border-top: 1px solid #ccc; margin: 20pt 0; }}

/* Tables - academic style */
table {{
    width: 100%;
    border-collapse: collapse;
    margin: 12pt 0;
    font-size: 9.5pt;
    line-height: 1.5;
    page-break-inside: avoid;
}}
thead tr {{
    border-top: 2px solid #000;
    border-bottom: 1px solid #000;
}}
th {{
    padding: 6px 8px;
    text-align: left;
    font-weight: bold;
    background: #f8f8f8;
}}
td {{
    padding: 5px 8px;
    border-bottom: 1px solid #e0e0e0;
    vertical-align: top;
}}
tbody tr:last-child {{
    border-bottom: 2px solid #000;
}}

/* Figure captions */
.fig-cap {{
    font-size: 9.5pt;
    text-align: justify;
    margin: 5px 0 25px 0;
    line-height: 1.5;
}}

/* Code blocks */
code {{
    font-family: "Courier New", monospace;
    font-size: 9pt;
    background: #f5f5f5;
    padding: 1px 3px;
}}
pre {{
    font-size: 8.5pt;
    background: #f5f5f5;
    padding: 8px;
    border: 1px solid #ddd;
    overflow-x: auto;
    line-height: 1.4;
}}
blockquote {{
    border-left: 3px solid #ccc;
    margin: 8pt 0;
    padding: 4pt 12pt;
    color: #444;
}}
strong {{ font-weight: bold; }}
em {{ font-style: italic; }}
</style>
</head>
<body>

{body_html}

{figures_html}

</body>
</html>"""

out_html = f"{DOCS}/ARIA_bioRxiv.html"
out_pdf = f"{DOCS}/ARIA_bioRxiv.pdf"

with open(out_html, "w") as f:
    f.write(full_html)

# Step 4: Convert to PDF
result = subprocess.run(["weasyprint", out_html, out_pdf], capture_output=True, text=True)
size = os.path.getsize(out_pdf) / 1024
print(f"bioRxiv PDF generated: {size:.0f} KB")
if result.stderr:
    warnings = [l for l in result.stderr.split('\n') if 'WARNING' in l]
    if len(warnings) > 3:
        print(f"({len(warnings)} CSS warnings, non-critical)")
