#!/usr/bin/env python3
"""
Train DeepSetsNoFun: 4 constituent features (jcsDEta, jcsDPhi, jcsM, jcsPt),
no jcsW, no jet-level scalars.

Outputs (written next to this script):
  best_model.pt     — state dict at best val loss
  scaler.pkl        — {'jcs': StandardScaler}  (fitted on 4 features)
  split_indices.npz — train_idx, val_idx, test_idx  (same seed as other models)
  loss_curve.pdf
"""

import os
import sys
import pickle
import pathlib

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split

HERE = pathlib.Path(__file__).parent
sys.path.insert(0, str(HERE))

from dataset import load_features, DijetDataset, BKG_FILE, SIG_FILE, N_CONSTIT_FEAT
from model import DeepSetsNoFun
from style import apply_style, FIG_SIZE

SEED       = 42
BATCH_SIZE = 512
MAX_EPOCHS = 200
LR         = 1e-3
DROPOUT    = 0.3
PATIENCE   = 15
SCHED_PAT  = 5
SCHED_FAC  = 0.5
VAL_FRAC   = 0.15
TEST_FRAC  = 0.15


def scale_jcs(x_jcs, scaler):
    shape = x_jcs.shape
    return scaler.transform(x_jcs.reshape(-1, N_CONSTIT_FEAT)).reshape(shape).astype(np.float32)


def main():
    torch.manual_seed(SEED)

    print('Loading background ...')
    x_jet_bkg, x_jcs_bkg, _ = load_features(HERE / BKG_FILE)
    print(f'  QCD: {len(x_jet_bkg)} events')

    print('Loading signal ...')
    x_jet_sig, x_jcs_sig, _ = load_features(HERE / SIG_FILE)
    print(f'  VBF: {len(x_jet_sig)} events')

    x_jcs_all = np.concatenate([x_jcs_bkg, x_jcs_sig], axis=0)
    y_all = np.concatenate([np.zeros(len(x_jet_bkg), dtype=np.float32),
                            np.ones (len(x_jet_sig), dtype=np.float32)], axis=0)

    idx = np.arange(len(y_all))
    idx_trainval, idx_test = train_test_split(
        idx, test_size=TEST_FRAC, stratify=y_all, random_state=SEED)
    val_frac_of_trainval = VAL_FRAC / (1.0 - TEST_FRAC)
    idx_train, idx_val = train_test_split(
        idx_trainval, test_size=val_frac_of_trainval,
        stratify=y_all[idx_trainval], random_state=SEED)

    np.savez(HERE / 'split_indices.npz',
             train_idx=idx_train, val_idx=idx_val, test_idx=idx_test)
    print(f'Split: train={len(idx_train)}, val={len(idx_val)}, test={len(idx_test)}')

    mask_train = (x_jcs_all[idx_train][..., 3] > 0)   # jcsPt at index 3
    mask_val   = (x_jcs_all[idx_val]  [..., 3] > 0)

    flat_train = x_jcs_all[idx_train].reshape(-1, N_CONSTIT_FEAT)
    valid_rows = flat_train[:, 3] > 0
    scaler_jcs = StandardScaler()
    scaler_jcs.fit(flat_train[valid_rows])

    x_jcs_train = scale_jcs(x_jcs_all[idx_train], scaler_jcs)
    x_jcs_val   = scale_jcs(x_jcs_all[idx_val],   scaler_jcs)

    with open(HERE / 'scaler.pkl', 'wb') as fh:
        pickle.dump({'jcs': scaler_jcs}, fh)

    y_train = y_all[idx_train]
    y_val   = y_all[idx_val]

    train_ds = DijetDataset(x_jcs_train, mask_train, y_train)
    val_ds   = DijetDataset(x_jcs_val,   mask_val,   y_val)

    train_loader = DataLoader(train_ds, batch_size=BATCH_SIZE, shuffle=True,  num_workers=2, pin_memory=True)
    val_loader   = DataLoader(val_ds,   batch_size=BATCH_SIZE, shuffle=False, num_workers=2, pin_memory=True)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'Device: {device}')

    model = DeepSetsNoFun(dropout=DROPOUT).to(device)
    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f'Parameters: {n_params:,}')

    n_bkg = int((y_train == 0).sum())
    n_sig = int((y_train == 1).sum())
    pos_weight = torch.tensor([n_bkg / n_sig], dtype=torch.float32, device=device)
    criterion  = nn.BCEWithLogitsLoss(pos_weight=pos_weight)

    optimizer = torch.optim.Adam(model.parameters(), lr=LR)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
                    optimizer, factor=SCHED_FAC, patience=SCHED_PAT)

    best_val_loss = float('inf')
    patience_counter = 0
    train_losses, val_losses = [], []

    for epoch in range(1, MAX_EPOCHS + 1):
        model.train()
        running = 0.0
        for x_jcs_b, mask_b, yb in train_loader:
            x_jcs_b = x_jcs_b.to(device)
            mask_b  = mask_b.to(device)
            yb      = yb.to(device)
            optimizer.zero_grad()
            loss = criterion(model(x_jcs_b, mask_b), yb)
            loss.backward()
            optimizer.step()
            running += loss.item() * len(yb)
        train_loss = running / len(train_ds)

        model.eval()
        running = 0.0
        with torch.no_grad():
            for x_jcs_b, mask_b, yb in val_loader:
                x_jcs_b = x_jcs_b.to(device)
                mask_b  = mask_b.to(device)
                yb      = yb.to(device)
                running += criterion(model(x_jcs_b, mask_b), yb).item() * len(yb)
        val_loss = running / len(val_ds)

        train_losses.append(train_loss)
        val_losses.append(val_loss)
        scheduler.step(val_loss)

        print(f'Epoch {epoch:3d}  train={train_loss:.4f}  val={val_loss:.4f}  lr={optimizer.param_groups[0]["lr"]:.2e}')

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            torch.save(model.state_dict(), HERE / 'best_model.pt')
            patience_counter = 0
        else:
            patience_counter += 1
            if patience_counter >= PATIENCE:
                print(f'Early stopping at epoch {epoch}.')
                break

    print(f'Best val loss: {best_val_loss:.4f}')

    epochs = np.arange(1, len(train_losses) + 1)
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    ax.plot(epochs, train_losses, linewidth=2, color='black',   label='train')
    ax.plot(epochs, val_losses,   linewidth=2, color='#1f77b4', label='val', linestyle='--')
    apply_style(ax, xlabel='Epoch', ylabel='Loss', title='DeepSets (no fun)', legend_loc='upper right')
    plt.tight_layout()
    fig.savefig(HERE / 'loss_curve.pdf', bbox_inches='tight')
    plt.close(fig)
    print(f'Saved {HERE / "loss_curve.pdf"}')


if __name__ == '__main__':
    os.chdir(HERE)
    main()
