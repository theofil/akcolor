#!/usr/bin/env python3
import os
import base64
import webbrowser
import tempfile
from pathlib import Path

FIGS = Path(__file__).resolve().parent
pngs = sorted(FIGS.glob('*.png'))

if not pngs:
    print('No PNG files found in', FIGS)
    raise SystemExit

def to_data_uri(path, mime):
    data = base64.b64encode(path.read_bytes()).decode()
    return f'data:{mime};base64,{data}'

rows = ''.join(
    f'<div class="fig">'
    f'<p class="name">{p.stem}</p>'
    f'<a href="{to_data_uri(p.with_suffix(".pdf"), "application/pdf")}" target="_blank">'
    f'<img src="{to_data_uri(p, "image/png")}" />'
    f'</a>'
    f'</div>'
    for p in pngs
    if p.with_suffix('.pdf').exists()
) + ''.join(
    f'<div class="fig">'
    f'<p class="name">{p.stem}</p>'
    f'<img src="{to_data_uri(p, "image/png")}" />'
    f'</div>'
    for p in pngs
    if not p.with_suffix('.pdf').exists()
)

html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Figures</title>
<style>
  body   {{ background: #1a1a1a; margin: 0; padding: 16px;
            font-family: sans-serif; color: #ccc; }}
  h1     {{ text-align: center; font-size: 1.1em; color: #888;
            margin: 0 0 16px; }}
  .grid  {{ display: flex; flex-wrap: wrap; gap: 16px;
            justify-content: center; }}
  .fig   {{ background: #2a2a2a; border-radius: 6px;
            padding: 8px; text-align: center; }}
  .name  {{ font-size: 0.8em; color: #888; margin: 0 0 6px; }}
  img    {{ max-width: 600px; width: 100%; display: block; }}
</style>
</head>
<body>
<h1>{FIGS}</h1>
<div class="grid">{rows}</div>
</body>
</html>"""

tmp = Path(tempfile.mktemp(suffix='.html'))
tmp.write_text(html)
print(f'Opening {len(pngs)} figures...')
webbrowser.open(tmp.as_uri())
