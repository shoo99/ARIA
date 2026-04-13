#!/usr/bin/env python3
"""Generate Korean version PDF."""
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

result = subprocess.run(
    [PANDOC, f"{DOCS}/manuscript_korean.md", "-t", "html", "--no-highlight"],
    capture_output=True, text=True
)
body = result.stdout

figures = f"""
<h2>Figures</h2>
{img_b64(f"{FIGS}/Fig1_architecture.png", "85%")}
<p class="fig-cap"><strong>Figure 1.</strong> ARIA 4계층 아키텍처. Layer 1: LLM 기반 추론 엔진. Layer 2: 8개 의사결정 지점(DP)으로 구성된 의사결정 프로토콜. Layer 3: Docker 컨테이너 및 스크립트 관리 실행 엔진. Layer 4: 6개 분석 모듈.</p>

{img_b64(f"{FIGS}/Fig2_decision_flow.png", "90%")}
<p class="fig-cap"><strong>Figure 2.</strong> ARIA 의사결정 흐름. QC 평가(DP1)와 설계 인식(DP3) 후 DE 분석 결과가 DP2에서 평가되어 DEG 수에 따라 하류 전략이 결정된다. 교차 방법 검증(DP5), 민감도 분석(DP7), 보고서 생성(DP8)이 적응적으로 수행된다.</p>

{img_b64(f"{FIGS}/Fig3_airway_benchmark.png", "95%")}
<p class="fig-cap"><strong>Figure 3.</strong> Airway 벤치마크 결과. (A) Unpaired(naive) 모델과 Paired(ARIA) 모델의 DEG 검출 비교. Paired 모델이 모든 cutoff에서 21–47% 더 많은 DEG를 검출하였다. (B) 7개의 알려진 dexamethasone 표적 유전자가 모두 강한 fold change(LFC 2.6–7.4)로 복원되었다.</p>

{img_b64(f"{FIGS}/Fig4_comparison.png", "80%")}
<p class="fig-cap"><strong>Figure 4.</strong> ARIA와 기존 도구 간 기능 비교. 초록(✓) = 완전 지원; 주황(△) = 부분 지원; 빨강(✗) = 미지원. ARIA는 결과 기반 적응, 설계 인식, 교차 방법 검증, 의사결정 투명성을 고유하게 제공한다.</p>
"""

html = f"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="UTF-8">
<style>
@page {{
    size: letter;
    margin: 1in 1in 1in 1.2in;
    @bottom-center {{ content: counter(page); font-size: 10pt; color: #666; }}
}}
body {{
    font-family: "Malgun Gothic", "Apple SD Gothic Neo", "Noto Sans KR", sans-serif;
    font-size: 10.5pt;
    line-height: 2.0;
    color: #000;
}}
h1 {{ font-size: 15pt; text-align: center; margin: 0 0 8px 0; line-height: 1.4; }}
h2 {{ font-size: 12.5pt; margin: 22pt 0 10pt 0; page-break-after: avoid; border-bottom: 1px solid #ccc; padding-bottom: 4px; }}
h3 {{ font-size: 11pt; margin: 16pt 0 6pt 0; page-break-after: avoid; }}
h4 {{ font-size: 10.5pt; font-style: italic; margin: 12pt 0 4pt 0; }}
p {{ margin: 0 0 7pt 0; text-align: justify; }}
ul, ol {{ margin: 6pt 0 6pt 20pt; }}
li {{ margin: 3pt 0; }}
hr {{ border: none; border-top: 1px solid #ccc; margin: 18pt 0; }}
table {{
    width: 100%;
    border-collapse: collapse;
    margin: 10pt 0;
    font-size: 9pt;
    line-height: 1.5;
    page-break-inside: avoid;
}}
thead tr {{ border-top: 2px solid #000; border-bottom: 1px solid #000; }}
th {{ padding: 5px 6px; text-align: left; font-weight: bold; background: #f8f8f8; }}
td {{ padding: 4px 6px; border-bottom: 1px solid #e0e0e0; vertical-align: top; }}
tbody tr:last-child {{ border-bottom: 2px solid #000; }}
.fig-cap {{ font-size: 9pt; text-align: justify; margin: 5px 0 20px 0; line-height: 1.5; }}
code {{ font-family: monospace; font-size: 9pt; background: #f5f5f5; padding: 1px 3px; }}
strong {{ font-weight: bold; }}
em {{ font-style: italic; }}
</style>
</head>
<body>
{body}
{figures}
</body>
</html>"""

out_html = f"{DOCS}/ARIA_korean.html"
out_pdf = f"{DOCS}/ARIA_korean.pdf"

with open(out_html, "w") as f:
    f.write(html)

subprocess.run(["weasyprint", out_html, out_pdf], capture_output=True)
size = os.path.getsize(out_pdf) / 1024
print(f"Korean PDF: {size:.0f} KB")
