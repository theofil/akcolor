import uproot
import numpy as np

samples = [
    {'label': 'QCD Hjj', 'file': '../QCDHtoInv.root'},
    {'label': 'QCD Zjj', 'file': '../QCDZtoInv.root'},
    {'label': 'QCD Wjj', 'file': '../QCDWtoInv.root'},
    {'label': 'VBF Hqq', 'file': '../VBFHtoInv.root'},
    {'label': 'VBF Zqq', 'file': '../VBFZtoInv.root'},
    {'label': 'VBF Wqq', 'file': '../VBFWtoInv.root'},
]


def load(file_path, branches, selection=None):
    with uproot.open(file_path) as f:
        arrays = f['events'].arrays(branches, library='np')
    if selection is not None:
        mask = selection(arrays)
        return {b: arrays[b][mask] for b in branches}
    return arrays


# ── total effective cross-sections ────────────────────────────────────────────
print('Effective cross-sections (sum of kWeight):')
for s in samples:
    data = load(s['file'], ['kWeight'])
    print(f"  {s['label']}: {data['kWeight'].sum():.6f}")

# ── cross-checks with published measurements ──────────────────────────────────

# EW Zqq: https://arxiv.org/abs/1712.09814
# mjj > 120, pTj > 25, mZ > 50 → 534 ± 20 (stat) ± 57 (syst) fb
data = load('../VBFZtoInv.root', ['kWeight', 'mjj'],
            selection=lambda a: a['mjj'] > 120)
print(f'\nVBF Zqq (mjj>120): {data["kWeight"].sum():.3f} pb')
print('  ref (arXiv:1712.09814): 534 ± 57 fb')

# EW Wqq: https://arxiv.org/abs/1903.04040
# mjj > 120, pTj > 25 → σEW(Wjj) = 6.23 ± 0.61 pb
data = load('../VBFWtoInv.root', ['kWeight', 'mjj'],
            selection=lambda a: a['mjj'] > 120)
print(f'\nVBF Wqq (mjj>120): {data["kWeight"].sum():.3f} pb')
print('  ref (arXiv:1903.04040): 6.23 ± 0.61 pb')

# EW Hqq: https://arxiv.org/pdf/2507.22574 p20
# σEW(Hjj) = 2.479 pb
data = load('../VBFHtoInv.root', ['kWeight', 'mjj', 'dYjj'],
            selection=lambda a: (a['mjj'] > 300) & (a['dYjj'] > 2))
print(f'\nVBF Hqq (mjj>300, dYjj>2): {data["kWeight"].sum():.3f} pb')
print('  ref (arXiv:2507.22574): 2.479 pb')
