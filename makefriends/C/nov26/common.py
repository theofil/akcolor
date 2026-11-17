"""Inputs and cross-fitting folds shared by the nov26 jet networks (NNjPF, ParTj).

One training example is one jet of friends/nov26/<sample>.h5 (Herwig only for training):
  role A (0): gate == 1, slot 0 (randomly the leading or subleading jet)   -> train / validation
  role B (1): gate == 1, slot 1 (the "unbiased 2nd jet")                   -> never trained on, scored only
  role C (2): gate == 0, every stored jet                                  -> train / validation
Every event gets one of 10 sub-folds (fixed seed), stratified within (sample, gate, nJets).
Fold model k infers on sub-folds {2k, 2k+1}, validates on (2k+2) % 10 and trains on the other
seven (70 / 10 / 20), so every event is scored once by a model that never saw it.
Details and the reasoning: nov26/training.md.
"""
import os
import time
import zlib

import h5py
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))     # makefriends/C
FRIENDS = os.path.join(BASE, 'friends', 'nov26')
LABEL = {'VBFZ_herwig': 1, 'QCDZjj_herwig': 0, 'VBFZ_mg5_pythia': 1, 'QCDZjj_mg5_pythia': 0}
TRAIN_SAMPLES = ('VBFZ_herwig', 'QCDZjj_herwig')
TRANSFER_SAMPLES = ('VBFZ_mg5_pythia', 'QCDZjj_mg5_pythia')

ROLE_A, ROLE_B, ROLE_C = 0, 1, 2
ROLE_NAMES = {ROLE_A: 'A (gate, slot 0)', ROLE_B: 'B (gate, slot 1)', ROLE_C: 'C (gate == 0)'}
NC_MAX = 80
N_SUB, N_FOLDS, SEED = 10, 5, 2026


def h5_path(sample):
    return os.path.join(FRIENDS, f'{sample}.h5')


# ── folds ─────────────────────────────────────────────────────────────────────

def assign_subfolds(sample, gate, n_jets, seed=SEED):
    """Sub-fold (0..9) of every event of one sample file, stratified within (gate, nJets).

    Within each (gate, nJets) group the events are shuffled with a seed derived from
    (seed, sample, group) and dealt out round-robin, so every sub-fold receives a tenth of
    each group (to within one event). Depends only on the file content, hence identical in
    every job and for both nets."""
    gate = np.asarray(gate).astype(np.int64)
    n_jets = np.asarray(n_jets).astype(np.int64)
    sub = np.full(len(gate), -1, dtype=np.int8)
    tag = zlib.crc32(sample.encode())
    for g in (0, 1):
        for nj in (0, 1, 2):
            idx = np.flatnonzero((gate == g) & (n_jets == nj))
            if len(idx):
                perm = np.random.default_rng([seed, tag, g, nj]).permutation(idx)
                sub[perm] = np.arange(len(perm)) % N_SUB
    assert (sub >= 0).all()
    return sub


def fold_sets(k):
    """(train, validation, inference) sub-folds of fold model k."""
    infer = (2 * k, 2 * k + 1)
    val = ((2 * k + 2) % N_SUB,)
    train = tuple(s for s in range(N_SUB) if s not in infer + val)
    return train, val, infer


# ── jets ──────────────────────────────────────────────────────────────────────

def load_jets(sample, subfolds=None, roles=(ROLE_A, ROLE_B, ROLE_C), with_p4=False, entry_stop=None):
    """Per-jet arrays for the stored jets of `sample` whose role is in `roles` and whose event
    is in one of `subfolds` (None = all). Jets come in (event, slot) order.

    x_jet  (n, 3)       |jetEta|, jetM, jetPt
    x_jcs  (n, 80, 3)   jcsDEta (signed, = dEtas), jcsDPhi, jcsPt / jetPt     [not standardised]
    jcs_pf (n, 80) int8 PF type 1 h+-, 2 e, 3 mu, 4 gamma, 5 h0; 0 = empty slot
    jcs_q  (n, 80) int8 charge
    mask   (n, 80) bool jcsPt > 0
    jcs_p4 (n, 80, 4)   with_p4 only: (px, py, pz, E) from pT, eta = |jetEta| + jcsDEta,
                        phi = jcsDPhi, m = jcsM (event reflected to put the jet at eta > 0 and
                        rotated to its phi: pair invariants unchanged); zero for empty slots
    kweight, ev, slot, role, sub, label  (n,)"""
    with h5py.File(h5_path(sample), 'r') as f:
        gate_all, nj_all = f['gate'][:], f['nJets'][:]
        sub_all = assign_subfolds(sample, gate_all, nj_all)          # always from the full file
        n = len(gate_all) if entry_stop is None else min(entry_stop, len(gate_all))
        gate, nj, sub = gate_all[:n], nj_all[:n], sub_all[:n]
        in_sub = np.ones(n, bool) if subfolds is None else np.isin(sub, subfolds)
        rows, slots = [], []
        for s in (0, 1):
            role_s = np.where(gate, ROLE_A if s == 0 else ROLE_B, ROLE_C)
            r = np.flatnonzero(in_sub & (s < nj) & np.isin(role_s, roles))
            rows.append(r)
            slots.append(np.full(len(r), s, np.int8))
        rows, slots = np.concatenate(rows), np.concatenate(slots)
        order = np.lexsort((slots, rows))
        rows, slots = rows[order], slots[order]

        def per_jet(name):          # (n_jets,) or (n_jets, 80) from an (N, 2) or (N, 2, 80) dataset
            return f[name][:n][rows, slots]

        eta, jet_m, jet_pt = per_jet('jetEta'), per_jet('jetM'), per_jet('jetPt')
        jcs_pt = per_jet('jcsPt')
        out = dict(
            x_jet=np.stack([np.abs(eta), jet_m, jet_pt], axis=1).astype(np.float32),
            x_jcs=np.stack([per_jet('jcsDEta'), per_jet('jcsDPhi'), jcs_pt / jet_pt[:, None]],
                           axis=2).astype(np.float32),
            jcs_pf=per_jet('jcsPFType').astype(np.int8),
            jcs_q=per_jet('jcsCharge').astype(np.int8),
            mask=jcs_pt > 0,
            kweight=f['kWeight'][:n][rows].astype(np.float32),
            ev=f['evIdx'][:n][rows].astype(np.int64),
            slot=slots,
            role=np.where(gate[rows], np.where(slots == 0, ROLE_A, ROLE_B), ROLE_C).astype(np.int8),
            sub=sub[rows],
            label=np.full(len(rows), LABEL[sample], np.int8),
        )
        if with_p4:
            deta, dphi, m = out['x_jcs'][..., 0], out['x_jcs'][..., 1], per_jet('jcsM')
            eta_c = np.abs(eta)[:, None] + deta
            px, py, pz = jcs_pt * np.cos(dphi), jcs_pt * np.sin(dphi), jcs_pt * np.sinh(eta_c)
            e = np.sqrt(px * px + py * py + pz * pz + m * m)
            out['jcs_p4'] = (np.stack([px, py, pz, e], axis=2) * out['mask'][..., None]).astype(np.float32)
    return out


def concat(parts):
    """Concatenate a list of load_jets() dicts."""
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}


def take(d, sel):
    """Rows `sel` (bool mask or indices) of every array in a load_jets() dict."""
    return {k: v[sel] for k, v in d.items()}


def class_balance_weights(label, role):
    """Per-jet training weight that removes the population-dependent class prior.

    The populations are the gate jets (roles A, B) and the gate == 0 jets (role C). Within
    population P each jet of class c gets w = N_P / (2 N_{P,c}): the two classes then carry
    equal weight inside P, and P keeps its total weight N_P (mean weight 1)."""
    pop = (np.asarray(role) == ROLE_C).astype(np.int64)
    label = np.asarray(label)
    w = np.zeros(len(label), np.float32)
    for p in (0, 1):
        in_p = pop == p
        for c in (0, 1):
            m = in_p & (label == c)
            if m.any():
                w[m] = in_p.sum() / (2.0 * m.sum())
    return w


# ── h5 output ─────────────────────────────────────────────────────────────────

def write_datasets(path, datasets, attempts=3):
    """Replace datasets in an h5 file, retrying: EOS FUSE occasionally fails HDF5 metadata
    operations with 'bad object header version number' (same as summer26/NNj/save_scores.py)."""
    for attempt in range(1, attempts + 1):
        try:
            with h5py.File(path, 'a') as f:
                for key, data in datasets.items():
                    if key in f:
                        del f[key]
                    f.create_dataset(key, data=data, compression='gzip', compression_opts=4)
            return
        except (KeyError, OSError, RuntimeError) as e:
            print(f'  write attempt {attempt}/{attempts} failed: {e}')
            if attempt == attempts:
                raise
            time.sleep(10)
