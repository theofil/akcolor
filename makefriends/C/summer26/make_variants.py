"""
Creates the 9 missing (network, process) directories for W and H variants.
Run once from the summer26/ directory or with python summer26/make_variants.py.
"""
import os
import re
import sys
from pathlib import Path

BASE = Path(__file__).parent

# NNjZ_H and NNjjZ_H are deliberately NOT in VARIANTS: they keep bosonM,
# which the base (Z) and _W dirs dropped as a Herwig Z/W generation artifact
# (see commit 9ce67b6). Regenerating them from base would silently drop the
# feature and desync them from the trained best_model_H.pt (input dim).
# They are maintained manually; check_boson_invariant() below enforces this.
VARIANTS = [
    ('NNj',    'NNj_H',    'H'),
    ('NNj',    'NNj_W',    'W'),
    ('NNjZ',   'NNjZ_W',   'W'),
    ('NNjj',   'NNjj_H',   'H'),
    ('NNjj',   'NNjj_W',   'W'),
    ('NNjjZ',  'NNjjZ_W',  'W'),
    ('NNkin',  'NNkin_H',  'H'),
    ('NNkin',  'NNkin_W',  'W'),
]


def check_boson_invariant(base=BASE):
    """bosonM must stay in the _H boson nets and out of the base/W ones."""
    must_have    = ['NNjZ_H/dataset.py', 'NNjjZ_H/dataset.py']
    must_not_have = ['NNjZ/dataset.py',  'NNjZ_W/dataset.py',
                     'NNjjZ/dataset.py', 'NNjjZ_W/dataset.py']
    errors = []
    for rel in must_have:
        if "'bosonM'" not in (base / rel).read_text():
            errors.append(f"{rel}: bosonM MISSING — was it regenerated from the "
                          f"bosonM-less base dir? Restore it from git.")
    for rel in must_not_have:
        if "'bosonM'" in (base / rel).read_text():
            errors.append(f"{rel}: bosonM present — it must stay excluded "
                          f"(Herwig Z/W on-shell generation artifact).")
    if errors:
        print('\nbosonM invariant VIOLATED:', file=sys.stderr)
        for e in errors:
            print(f'  {e}', file=sys.stderr)
        sys.exit(1)
    print('bosonM invariant OK (kept in NNjZ_H/NNjjZ_H, excluded from base/W).')

SKIP_DIRS = {'logs', '__pycache__'}
SKIP_EXTS = {'.pt', '.pkl', '.npz', '.pdf'}
# save_scores runs once from the base (Z) dir and writes the NNj_jet0/jet1
# h5 columns — variants must not get a copy (they would clobber those columns)
SKIP_PREFIXES = ('save_scores',)


def apply_subs(text, source, dest, proc):
    # Functional: PROCESS variable and data file names
    text = re.sub(r'PROCESS\s*=\s*"Z"', f'PROCESS = "{proc}"', text)
    text = text.replace('QCDZjj', f'QCD{proc}jj')
    text = text.replace('VBFZ',   f'VBF{proc}')
    # Paths in .sh and .sub files
    text = text.replace(f'/summer26/{source}', f'/summer26/{dest}')
    # Cosmetic: print statements and comments
    text = text.replace('Herwig Z ',      f'Herwig {proc} ')
    text = text.replace('Herwig Z\n',     f'Herwig {proc}\n')
    text = text.replace('Herwig Z.',      f'Herwig {proc}.')
    text = text.replace('MG5_Pythia Z ',  f'MG5_Pythia {proc} ')
    text = text.replace('MG5_Pythia Z\n', f'MG5_Pythia {proc}\n')
    text = text.replace('MG5_Pythia Z.', f'MG5_Pythia {proc}.')
    text = text.replace('MG5+Py Z ',     f'MG5+Py {proc} ')
    text = text.replace('VBF Z ',        f'VBF {proc} ')
    text = text.replace('VBF Z\n',       f'VBF {proc}\n')
    text = text.replace('VBF Z.',        f'VBF {proc}.')
    text = text.replace('VBF Z vs',      f'VBF {proc} vs')
    text = text.replace('QCD Z+',        f'QCD {proc}+')
    text = text.replace('QCD Z ',        f'QCD {proc} ')
    return text


for source, dest, proc in VARIANTS:
    src_dir = BASE / source
    dst_dir = BASE / dest

    if dst_dir.exists():
        # Check if fully populated; if not, fill in missing files
        existing = {f.name for f in dst_dir.iterdir() if f.is_file()}
        source_files = {f.name for f in src_dir.iterdir()
                        if f.is_file() and f.name not in SKIP_DIRS
                        and f.suffix not in SKIP_EXTS
                        and not f.name.startswith(SKIP_PREFIXES)}
        missing = source_files - existing
        if not missing:
            print(f'SKIP {dest} — already complete')
            continue
        print(f'PATCH {dest}/ — filling in: {sorted(missing)}')
    else:
        print(f'Creating {dest}/ from {source}/ (PROCESS={proc})')
        dst_dir.mkdir()
        (dst_dir / 'logs').mkdir()

    for item in sorted(src_dir.iterdir()):
        if item.name in SKIP_DIRS or item.suffix in SKIP_EXTS:
            continue
        if item.name.startswith(SKIP_PREFIXES):
            continue
        if not item.is_file():
            continue

        dst_file = dst_dir / item.name
        try:
            text = item.read_text()
        except OSError:
            # EOS transient error: fall back to cp then read from destination
            import subprocess
            subprocess.run(['cp', str(item), str(dst_file)], check=True)
            text = dst_file.read_text()

        text = apply_subs(text, source, dest, proc)
        dst_file.write_text(text)
        if os.access(item, os.X_OK):
            dst_file.chmod(dst_file.stat().st_mode | 0o111)

check_boson_invariant()
print('\nAll done.')
