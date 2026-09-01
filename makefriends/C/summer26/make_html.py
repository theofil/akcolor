#!/usr/bin/env python3
"""Regenerate figs/summer26/summer26.html from summer26.md."""
import pathlib
import re
import html as _html

_MD_PATH   = pathlib.Path(__file__).parent.parent / 'summer26.md'
_HTML_PATH = pathlib.Path(__file__).parent.parent / 'figs' / 'summer26' / 'summer26.html'

_CSS = '''
:root{--bg:#ffffff;--bg2:#f8fafc;--bg3:#f1f5f9;--border:#e2e8f0;--fg:#1e293b;--fg-dim:#64748b;--h1:#92400e;--h2:#1e40af;--h3:#0e7490;--code-fg:#be185d;--pre-fg:#334155;--link:#1d4ed8;}
*,*::before,*::after{box-sizing:border-box;margin:0;padding:0;}html{scroll-behavior:smooth;}
body{background:var(--bg);color:var(--fg);font-family:'JetBrains Mono','Fira Code','Cascadia Code','Menlo','Consolas',monospace;font-size:13.5px;line-height:1.75;display:flex;}
nav{width:220px;min-width:220px;background:var(--bg2);border-right:1px solid var(--border);height:100vh;position:sticky;top:0;overflow-y:auto;padding:1.4rem 0;font-size:11.5px;}
nav .nav-title{color:var(--h1);font-size:11px;letter-spacing:.1em;text-transform:uppercase;padding:0 1rem .8rem;border-bottom:1px solid var(--border);margin-bottom:.6rem;}
nav a{display:block;padding:.2rem 1rem;color:var(--fg-dim);text-decoration:none;border-left:2px solid transparent;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;}
nav a.h2{padding-left:1rem;}nav a.h3{padding-left:1.7rem;font-size:11px;}
nav a:hover{color:var(--h2);border-left-color:var(--h2);}
main{flex:1;min-width:0;padding:2.5rem 3rem 4rem;max-width:900px;}
h1{color:var(--h1);font-size:1.45rem;margin:0 0 .6rem;padding-bottom:.4rem;border-bottom:1px solid var(--border);line-height:1.3;}
h2{color:var(--h2);font-size:1.05rem;margin:2rem 0 .5rem;padding-bottom:.25rem;border-bottom:1px solid var(--border);}
h3{color:var(--h3);font-size:.95rem;margin:1.3rem 0 .3rem;}
h4{color:var(--fg);font-size:.9rem;margin:1rem 0 .2rem;}
p{margin:.45rem 0;}hr{border:none;border-top:1px solid var(--border);margin:1.8rem 0;}
a{color:var(--link);text-decoration:none;}a:hover{text-decoration:underline;}
code{background:var(--bg3);color:var(--code-fg);padding:.1em .35em;border-radius:3px;font-size:.88em;border:1px solid var(--border);}
pre{background:var(--bg2);border:1px solid var(--border);border-left:3px solid var(--h2);border-radius:4px;padding:.9rem 1.1rem;margin:.7rem 0;overflow-x:auto;}
pre code{background:none;border:none;color:var(--pre-fg);font-size:.84em;padding:0;}
.tbl-wrap{overflow-x:auto;margin:.7rem 0;}
table{border-collapse:collapse;width:100%;font-size:.82em;white-space:nowrap;}
th{background:var(--bg3);color:var(--h3);text-align:left;padding:.35rem .65rem;border:1px solid var(--border);font-weight:normal;}
td{padding:.28rem .65rem;border:1px solid var(--border);color:var(--fg);}
tr:nth-child(even) td{background:var(--bg2);}tr:hover td{background:var(--bg3);}
ul,ol{padding-left:1.6rem;margin:.4rem 0;}li{margin:.2rem 0;}
blockquote{border-left:3px solid var(--border);padding:.3rem .9rem;color:var(--fg-dim);margin:.5rem 0;font-style:italic;}
strong{color:#0f172a;font-weight:600;}
::-webkit-scrollbar{width:6px;height:6px;}::-webkit-scrollbar-track{background:var(--bg);}
::-webkit-scrollbar-thumb{background:var(--border);border-radius:3px;}
'''

_NAV = '''
<div class="nav-title">summer26</div>
<a class="h2" href="#campaign-details">Campaign details</a>
<a class="h2" href="#step-2-run-event-reconstruction">Step 2 — Reconstruction</a>
<a class="h2" href="#step-4-diagnostic-plots">Step 4 — Plots</a>
<a class="h2" href="#step-5-train-the-neural-networks">Step 5 — Train NNs</a>
<a class="h3" href="#nnkin-summer26-nnkin">NNkin</a>
<a class="h3" href="#nnj-summer26-nnj">NNj</a>
<a class="h3" href="#nnjb-summer26-nnjb">NNjB</a>
<a class="h3" href="#nnjj-summer26-nnjj">NNjj</a>
<a class="h3" href="#nnjjb-summer26-nnjjb">NNjjB</a>
<a class="h3" href="#nnjjbj-summer26-nnjjbj">NNjjBj</a>
<a class="h2" href="#step-6-run-inference">Step 6 — Inference</a>
<a class="h2" href="#step-7-head-to-head-comparison">Step 7 — Comparison</a>
<a class="h2" href="#step-8-publish-plots">Step 8 — Publish</a>
<a class="h2" href="#kweight-factors-and-run-3-statistical-uncertainty">kWeight &amp; Run 3</a>
<a class="h3" href="#where-evweight-lives">evweight source</a>
<a class="h3" href="#kweight-formula">kWeight formula</a>
<a class="h3" href="#per-sample-values">Per-sample values</a>
<a class="h3" href="#character-of-the-weights-by-generator">Generator weights</a>
<a class="h2" href="#signal-background-estimation-from-roc-curves">S/B estimation</a>
<a class="h2" href="#input-files">Input files</a>
<a class="h2" href="#file-map">File map</a>
'''


def _md_to_html_body(md):
    lines = md.split('\n')
    out = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.strip().startswith('```'):
            i += 1
            code_lines = []
            while i < len(lines) and not lines[i].strip().startswith('```'):
                code_lines.append(_html.escape(lines[i]))
                i += 1
            out.append('<pre><code>' + '\n'.join(code_lines) + '</code></pre>')
            i += 1; continue
        if '|' in line and i+1 < len(lines) and re.match(r'[\s|:-]+$', lines[i+1].replace('|', '')):
            headers = [_html.escape(c.strip()) for c in line.strip().strip('|').split('|')]
            i += 2
            out.append('<div class="tbl-wrap"><table>')
            out.append('<tr>' + ''.join(f'<th>{h}</th>' for h in headers) + '</tr>')
            while i < len(lines) and '|' in lines[i]:
                cells = [_html.escape(c.strip()) for c in lines[i].strip().strip('|').split('|')]
                cells = [re.sub(r'`([^`]+)`', r'<code>\1</code>', c) for c in cells]
                cells = [re.sub(r'\*\*([^*]+)\*\*', r'<strong>\1</strong>', c) for c in cells]
                out.append('<tr>' + ''.join(f'<td>{c}</td>' for c in cells) + '</tr>')
                i += 1
            out.append('</table></div>'); continue
        m = re.match(r'^(#{1,4})\s+(.*)', line)
        if m:
            lvl = len(m.group(1))
            text = _html.escape(m.group(2))
            text = re.sub(r'`([^`]+)`', r'<code>\1</code>', text)
            anchor = re.sub(r'[^a-z0-9]+', '-', m.group(2).lower()).strip('-')
            out.append(f'<h{lvl} id="{anchor}">{text}</h{lvl}>')
            i += 1; continue
        if re.match(r'^---+\s*$', line):
            out.append('<hr>'); i += 1; continue
        if line.startswith('> '):
            out.append(f'<blockquote>{_html.escape(line[2:])}</blockquote>'); i += 1; continue
        if re.match(r'^[-*]\s+', line):
            out.append('<ul>')
            while i < len(lines) and re.match(r'^[-*]\s+', lines[i]):
                t = _html.escape(lines[i][2:].strip())
                t = re.sub(r'`([^`]+)`', r'<code>\1</code>', t)
                t = re.sub(r'\*\*([^*]+)\*\*', r'<strong>\1</strong>', t)
                out.append(f'<li>{t}</li>'); i += 1
            out.append('</ul>'); continue
        if re.match(r'^\d+\.\s+', line):
            out.append('<ol>')
            while i < len(lines) and re.match(r'^\d+\.\s+', lines[i]):
                t = _html.escape(re.sub(r'^\d+\.\s+', '', lines[i]))
                t = re.sub(r'`([^`]+)`', r'<code>\1</code>', t)
                t = re.sub(r'\*\*([^*]+)\*\*', r'<strong>\1</strong>', t)
                out.append(f'<li>{t}</li>'); i += 1
            out.append('</ol>'); continue
        if not line.strip(): i += 1; continue
        text = _html.escape(line)
        text = re.sub(r'`([^`]+)`', r'<code>\1</code>', text)
        text = re.sub(r'\*\*([^*]+)\*\*', r'<strong>\1</strong>', text)
        text = re.sub(r'\[([^\]]+)\]\(([^)]+)\)', r'<a href="\2">\1</a>', text)
        if line.endswith('  '): text = text.rstrip() + '<br>'
        out.append(f'<p>{text}</p>'); i += 1
    return '\n'.join(out)


def main():
    if not _MD_PATH.exists():
        print(f'Warning: {_MD_PATH} not found, skipping HTML generation')
        return
    body = _md_to_html_body(_MD_PATH.read_text())
    html_doc = f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1.0">
<title>summer26 analysis</title>
<style>{_CSS}</style>
</head>
<body>
<nav>{_NAV}</nav>
<main>
{body}
</main>
</body>
</html>'''
    _HTML_PATH.parent.mkdir(parents=True, exist_ok=True)
    _HTML_PATH.write_text(html_doc)
    print('Saved', _HTML_PATH)


if __name__ == '__main__':
    main()
