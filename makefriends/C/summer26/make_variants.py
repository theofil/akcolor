"""
Creates the 9 missing (network, process) directories for W and H variants.
Run once from the summer26/ directory or with python summer26/make_variants.py.
"""
import os
import re
from pathlib import Path

BASE = Path(__file__).parent

VARIANTS = [
    ('NNFullRaw', 'NNFullRaw_W', 'W'),
    ('NNkin',     'NNkin_H',     'H'),
    ('NNkin',     'NNkin_W',     'W'),
    ('NNPrePro',  'NNPrePro_H',  'H'),
    ('NNPrePro',  'NNPrePro_W',  'W'),
    ('NNpull',    'NNpull_H',    'H'),
    ('NNpull',    'NNpull_W',    'W'),
    ('NNRaw',     'NNRaw_H',     'H'),
    ('NNRaw',     'NNRaw_W',     'W'),
]

SKIP_DIRS = {'logs', '__pycache__'}
SKIP_EXTS = {'.pt', '.pkl', '.npz', '.pdf'}


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
                        and f.suffix not in SKIP_EXTS}
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

print('\nAll done.')
