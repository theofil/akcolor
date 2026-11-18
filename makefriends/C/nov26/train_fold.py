#!/usr/bin/env python3
"""Train fold model k of a nov26 jet network, then score the jets of its inference sub-folds.

    python3 train_fold.py --net NNjBC --fold 0                           # condor GPU job
    python3 train_fold.py --net NNjBCdphi --fold 3 --smoke --outdir <dir>  # quick CPU check
    python3 train_fold.py --net NNj --fold 0 --load-only                 # data loading and memory only

--net is one of common.NETS: NNj (no beam connection), NNjBC (all four bc variables), NNjBC<v>
(NNj + bc variable v alone) and NNjBC_no<v> (NNjBC without v). All share the DeepSets of
model.py and the training recipe below; they differ only in the bc inputs appended to the
jet scalars |eta|, m, pT.

Phase 1 (training): Herwig jets of roles A and C (common.py) in the 7 training sub-folds and
the validation sub-fold of fold k. A per-fold StandardScaler is fitted on the training jets
(jet scalars incl. the bc inputs; continuous constituent features of valid constituents only).
The loss is BCEWithLogits weighted per jet with common.class_balance_weights; early stopping
on the validation loss.
Phase 2 (inference): the best model scores every stored jet (roles A, B, C) of the two
inference sub-folds, in the Herwig and the MG5+Pythia files -> <outdir>/scores_fold{k}.npz,
merged into the h5 files by save_scores.py. See training.md.
"""
import argparse
import os
import pickle
import resource
import sys
import time

import numpy as np
import torch
import torch.nn.functional as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402
from model import JetNN  # noqa: E402

CFG = dict(dropout=0.3, batch=512, infer_batch=4096, lr=1e-3, max_epochs=100, patience=15, sched_patience=5)
KEYS = ('x_jet', 'x_jcs', 'mask')
PER_CONSTITUENT = ('x_jcs', 'mask')
T0 = time.time()


def log(msg):
    print(f'[{time.time() - T0:7.0f}s] {msg}', flush=True)


def max_rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


# ── data ──────────────────────────────────────────────────────────────────────

def net_inputs(d, net):
    """Append the net's bc columns to the jet scalars, in place; x_bc is dropped."""
    d['x_jet'] = np.concatenate([d['x_jet'], d.pop('x_bc')[:, C.bc_columns(net)]], axis=1)
    return d


def load(samples, subfolds, roles, net, entry_stop):
    parts = []
    for s in samples:
        p = net_inputs(C.load_jets(s, subfolds=subfolds, roles=roles, entry_stop=entry_stop), net)
        log(f'  {s}: {len(p["label"])} jets, roles A/B/C = {np.bincount(p["role"], minlength=3).tolist()}')
        parts.append(p)
    return C.concat(parts)


def fit_scalers(d):
    from sklearn.preprocessing import StandardScaler
    return {'jet': StandardScaler().fit(d['x_jet']), 'jcs': StandardScaler().fit(d['x_jcs'][d['mask']])}


def apply_scalers(d, sc):
    """Standardise in place (padded constituents are standardised too; they are masked)."""
    d['x_jet'] = ((d['x_jet'] - sc['jet'].mean_) / sc['jet'].scale_).astype(np.float32)
    d['x_jcs'] -= sc['jcs'].mean_.astype(np.float32)
    d['x_jcs'] /= sc['jcs'].scale_.astype(np.float32)


def to_device(d, device):
    return {k: torch.from_numpy(np.ascontiguousarray(d[k])).to(device) for k in KEYS}


def get_batch(t, idx):
    """Rows idx, with the constituent axis trimmed to the largest constituent count in the batch
    (constituents are pT-ordered and zero-padded at the end, so this drops only padding)."""
    b = {k: t[k][idx] for k in KEYS}
    n = max(int(b['mask'].sum(1).max()), 1)
    for k in PER_CONSTITUENT:
        b[k] = b[k][:, :n]
    return b


# ── training / scoring ────────────────────────────────────────────────────────

def run_epoch(model, t, y, w, order, optimizer=None, max_batches=None):
    """Weighted mean BCE over `order`; trains when an optimizer is given."""
    training = optimizer is not None
    model.train(training)
    tot, wsum = 0.0, 0.0
    with torch.set_grad_enabled(training):
        for ib, i in enumerate(range(0, len(order), CFG['batch'])):
            if max_batches is not None and ib >= max_batches:
                break
            idx = order[i:i + CFG['batch']]
            if training and len(idx) < 2:                        # BatchNorm needs > 1 jet
                continue
            b = get_batch(t, idx)
            logit = model(*[b[k] for k in KEYS]).squeeze(1)
            wb = w[idx]
            loss = (F.binary_cross_entropy_with_logits(logit, y[idx], reduction='none') * wb).sum() / wb.sum()
            if training:
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                optimizer.step()
            tot += loss.item() * wb.sum().item()
            wsum += wb.sum().item()
    return tot / wsum


def score(model, d, device):
    t = to_device(d, device)
    model.eval()
    out = []
    with torch.no_grad():
        for i in range(0, len(d['label']), CFG['infer_batch']):
            idx = torch.arange(i, min(i + CFG['infer_batch'], len(d['label'])), device=device)
            b = get_batch(t, idx)
            out.append(torch.sigmoid(model(*[b[k] for k in KEYS]).squeeze(1)).float().cpu().numpy())
    return np.concatenate(out) if out else np.zeros(0, np.float32)


def auc(label, s):
    from sklearn.metrics import roc_auc_score
    return float(roc_auc_score(label, s)) if len(np.unique(label)) == 2 else float('nan')


def save_loss_curve(hist, path, title):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4))
    ep = np.arange(1, len(hist['train']) + 1)
    ax.plot(ep, hist['train'], color='#2a78d6', lw=2, label='train')
    ax.plot(ep, hist['val'], color='#eb6834', lw=2, ls='--', label='validation')
    ax.set_xlabel('epoch')
    ax.set_ylabel('class-balanced BCE')
    ax.set_title(title, fontsize=10)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--net', required=True, choices=list(C.NETS))
    ap.add_argument('--fold', required=True, type=int, choices=range(C.N_FOLDS))
    ap.add_argument('--outdir', default=None, help='default: nov26/<net>')
    ap.add_argument('--smoke', action='store_true', help='3000 entries per sample, 2 batches, 1 epoch')
    ap.add_argument('--load-only', action='store_true', help='load the training data, report memory, exit')
    args = ap.parse_args()

    cfg = dict(CFG)
    outdir = args.outdir or os.path.join(C.OUT, args.net)
    os.makedirs(outdir, exist_ok=True)
    entry_stop = 3000 if args.smoke else None
    max_batches = 2 if args.smoke else None
    if args.smoke:
        cfg.update(max_epochs=1)
    k = args.fold
    train_sub, val_sub, infer_sub = C.fold_sets(k)
    torch.manual_seed(C.SEED + k)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    gpu = torch.cuda.get_device_name(0) if device.type == 'cuda' else 'none'
    log(f'{args.net} fold {k} (bc inputs: {", ".join(C.NETS[args.net]) or "none"}): train sub-folds {train_sub}, '
        f'validation {val_sub}, inference {infer_sub}; device {device} ({gpu}); config {cfg}')

    # ── phase 1: training ─────────────────────────────────────────────────────
    log(f'loading Herwig training/validation jets (roles A, C): {", ".join(C.TRAIN_SAMPLES)}')
    d = load(C.TRAIN_SAMPLES, train_sub + val_sub, (C.ROLE_A, C.ROLE_C), args.net, entry_stop)
    is_train = np.isin(d['sub'], train_sub)
    tr, va = C.take(d, is_train), C.take(d, ~is_train)
    del d
    for name, part in (('train', tr), ('validation', va)):
        lab, rol = part['label'], part['role']
        log(f'  {name}: {len(lab)} jets | VBF A/C = {int(((lab == 1) & (rol == 0)).sum())}/'
            f'{int(((lab == 1) & (rol == 2)).sum())}, QCD A/C = {int(((lab == 0) & (rol == 0)).sum())}/'
            f'{int(((lab == 0) & (rol == 2)).sum())}')
    n_jet_feat = tr['x_jet'].shape[1]
    log(f'jet scalars: {n_jet_feat} (|eta|, m, pT + bc columns {C.bc_columns(args.net)}); '
        f'finite: {bool(np.isfinite(tr["x_jet"]).all() and np.isfinite(va["x_jet"]).all())}')
    scalers = fit_scalers(tr)
    apply_scalers(tr, scalers)
    apply_scalers(va, scalers)
    w_tr = C.class_balance_weights(tr['label'], tr['role'])
    w_va = C.class_balance_weights(va['label'], va['role'])
    log(f'data ready; max RSS {max_rss_gb():.2f} GB')
    if args.load_only:
        nbytes = sum(v.nbytes for v in tr.values()) + sum(v.nbytes for v in va.values())
        log(f'load-only: train+validation arrays {nbytes / 1e9:.2f} GB; max RSS {max_rss_gb():.2f} GB')
        return
    with open(os.path.join(outdir, f'scaler_fold{k}.pkl'), 'wb') as fh:
        pickle.dump(scalers, fh)

    t_tr, t_va = to_device(tr, device), to_device(va, device)
    y_tr = torch.from_numpy(tr['label'].astype(np.float32)).to(device)
    y_va = torch.from_numpy(va['label'].astype(np.float32)).to(device)
    ww_tr, ww_va = torch.from_numpy(w_tr).to(device), torch.from_numpy(w_va).to(device)
    n_tr, n_va = len(y_tr), len(y_va)
    del tr, va

    model = JetNN(n_jet_feat, dropout=cfg['dropout']).to(device)
    log(f'model JetNN: {sum(p.numel() for p in model.parameters())} parameters')
    optimizer = torch.optim.Adam(model.parameters(), lr=cfg['lr'])
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, factor=0.5, patience=cfg['sched_patience'])
    best_path = os.path.join(outdir, f'best_model_fold{k}.pt')
    best, wait, hist = float('inf'), 0, {'train': [], 'val': []}
    val_order = torch.arange(n_va, device=device)
    for epoch in range(1, cfg['max_epochs'] + 1):
        te = time.time()
        loss_tr = run_epoch(model, t_tr, y_tr, ww_tr, torch.randperm(n_tr, device=device), optimizer, max_batches)
        loss_va = run_epoch(model, t_va, y_va, ww_va, val_order, None, max_batches)
        hist['train'].append(loss_tr)
        hist['val'].append(loss_va)
        scheduler.step(loss_va)
        log(f'epoch {epoch:3d}  train {loss_tr:.5f}  val {loss_va:.5f}  '
            f'lr {optimizer.param_groups[0]["lr"]:.1e}  ({time.time() - te:.0f}s)')
        if loss_va < best:
            best, wait = loss_va, 0
            torch.save(model.state_dict(), best_path)
        else:
            wait += 1
            if wait >= cfg['patience']:
                log(f'early stopping at epoch {epoch}')
                break
    log(f'best validation loss {best:.5f}')
    np.savez(os.path.join(outdir, f'history_fold{k}.npz'), train=hist['train'], val=hist['val'])
    figdir = outdir if args.smoke else C.FIGS
    os.makedirs(figdir, exist_ok=True)
    save_loss_curve(hist, os.path.join(figdir, f'{args.net}_fold{k}_loss.pdf'), f'{args.net}, fold {k}')
    del t_tr, t_va, y_tr, y_va, ww_tr, ww_va
    if device.type == 'cuda':
        torch.cuda.empty_cache()

    # ── phase 2: score every stored jet of the inference sub-folds ───────────
    model.load_state_dict(torch.load(best_path, map_location=device))
    out = {'net': args.net, 'fold': k, 'infer_subfolds': np.array(infer_sub)}
    for sample in C.TRAIN_SAMPLES + C.TRANSFER_SAMPLES:
        d = net_inputs(C.load_jets(sample, subfolds=infer_sub, entry_stop=entry_stop), args.net)
        apply_scalers(d, scalers)
        s = score(model, d, device)
        out[f'{sample}__ev'], out[f'{sample}__slot'], out[f'{sample}__score'] = d['ev'], d['slot'], s
        out[f'{sample}__role'] = d['role']
        out[f'{sample}__kweight'] = d['kweight']
        log(f'scored {sample}: {len(s)} jets')
    for gen, (sig, bkg) in (('herwig', C.TRAIN_SAMPLES), ('pythia', C.TRANSFER_SAMPLES)):
        for name, roles in (('gate jets (A+B)', (C.ROLE_A, C.ROLE_B)), ('role B', (C.ROLE_B,)),
                            ('gate==0 jets (C)', (C.ROLE_C,))):
            ms = np.isin(out[f'{sig}__role'], roles)
            mb = np.isin(out[f'{bkg}__role'], roles)
            a = auc(np.r_[np.ones(ms.sum()), np.zeros(mb.sum())],
                    np.r_[out[f'{sig}__score'][ms], out[f'{bkg}__score'][mb]])
            out[f'auc_{gen}_{name.split()[0].replace("==", "eq")}'] = a
            log(f'AUC {gen:6s} {name:17s}: {a:.4f}')
    path = os.path.join(outdir, f'scores_fold{k}.npz')
    np.savez(path, **out)
    log(f'wrote {path}; max RSS {max_rss_gb():.2f} GB')


if __name__ == '__main__':
    main()
