#!/usr/bin/env python3
"""Build MOD/friends/<sample>.h5 friend files from Zenodo MOD-format CMS Open Data
(https://zenodo.org/records/3340205 and its 8 Pythia6 QCD simulation companions).

Reproduces the same physics computations and output schema as makefriends.cpp +
root_to_hdf5.py (see MOD/README.md for the full derivation), so that the existing
summer26/NN*/dataset.py training code can read MOD/friends/*.h5 unchanged.

Key differences from makefriends.cpp, forced by the input format (see MOD/README.md
"Showstoppers" section for the full rationale):
  - Only 2 jets are ever available per event (MOD stores just the two hardest
    jets) -> no "slot 2" branches at all (jetPt2, jcsPt2, ... simply don't exist).
  - No boson truth in this dataset family -> bosonM/Y/Eta/Pt/Phi are always -99.
  - No re-clustering: jets are taken as pre-formed by MOD; only their own stored
    PF candidates are used for constituent/pull-vector features.
  - No opposite-eta-hemisphere dijet selection (that was a VBF-specific cut with
    no motivation for generic QCD dijets) -- kept selection is just the dataset's
    own recommended jet quality (>=2, "medium") and |eta|<1.9 cuts, on both jets.
  - kWeight = the jets_f 'weight' column directly (present in both DATA and SIM,
    despite the Zenodo docs implying it's SIM-only) -- no --xs/sumWtot rescaling.

Usage:
    python3 make_mod_friends.py [sample_name ...]   # default: all samples present in MOD/raw
"""

import glob
import os
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RAW_DIR = os.path.join(HERE, "raw")
OUT_DIR = os.path.join(HERE, "friends")

NC_MAX = 80
JET_ETA_MAX = 1.9
QUALITY_MIN = 2
PV_A, PV_B = 1.0, 1.0  # PV_C=0 in makefriends.cpp -> JET_R^0 = 1, dropped below

SAMPLES = {
    "CMS2011A_data": "data_3340205",
    "QCD_170-300":   "sim_170-300_3341500",
    "QCD_300-470":   "sim_300-470_3341498",
    "QCD_470-600":   "sim_470-600_3341419",
    "QCD_600-800":   "sim_600-800_3364139",
    "QCD_800-1000":  "sim_800-1000_3341413",
    "QCD_1000-1400": "sim_1000-1400_3341502",
    "QCD_1400-1800": "sim_1400-1800_3341770",
    "QCD_1800-inf":  "sim_1800-inf_3341772",
}

SCALAR_JETS = ["jetPt", "jetEta", "jetPhi", "jetM", "jetSPVA", "jetPVA", "jetPVM", "jetNC"]
CONSTIT_VARS = ["jcsPt", "jcsDEta", "jcsDEtaRaw", "jcsDPhi", "jcsM", "jcsW"]
SCALAR_EVT = ["kWeight", "mjj", "dYjj", "dPhijj", "ptjj",
              "bosonM", "bosonY", "bosonEta", "bosonPt", "bosonPhi"]


def delta_phi(p1, p2):
    return np.arctan2(np.sin(p1 - p2), np.cos(p1 - p2))


def derive_eta(pt, y, m):
    """Pseudorapidity from (pt, rapidity, mass), without needing px/py/pz/E."""
    mt = np.sqrt(m * m + pt * pt)
    pz = mt * np.sinh(y)
    return np.arcsinh(pz / pt)


def p4_from_ptymphi(pt, y, phi, m):
    mt = np.sqrt(m * m + pt * pt)
    px = pt * np.cos(phi)
    py = pt * np.sin(phi)
    pz = mt * np.sinh(y)
    e = mt * np.cosh(y)
    return px, py, pz, e


def group_events(jets_i, cols_i):
    """Return, for each event, the (>=1) row indices into jets_i/jets_f, grouped by
    (fn, rn, lbn, evn). Order within a group is arbitrary (caller re-sorts by pt)."""
    key_idx = [cols_i.index(c) for c in ("fn", "rn", "lbn", "evn")]
    keys = jets_i[:, key_idx]
    order = np.lexsort(keys.T[::-1])
    sorted_keys = keys[order]
    is_new = np.empty(len(sorted_keys), dtype=bool)
    is_new[0] = True
    is_new[1:] = np.any(sorted_keys[1:] != sorted_keys[:-1], axis=1)
    group_starts = np.where(is_new)[0]
    group_ends = np.concatenate([group_starts[1:], [len(sorted_keys)]])
    for s, e in zip(group_starts, group_ends):
        yield order[s:e]


def process_jet(pt, y, phi, m, jet_pt, jet_y, jet_phi, jet_eta):
    """Given a jet's own constituents (arrays) and the jet's own (pt,y,phi,eta),
    compute (pvm, pva, spva, nc_raw, jcs arrays[NC_MAX]) exactly mirroring
    fillPV()/the storage loop in makefriends.cpp, using derived eta for jcsDEta*."""
    nc_raw = len(pt)
    jcsPt = np.zeros(NC_MAX, dtype=np.float32)
    jcsDEta = np.zeros(NC_MAX, dtype=np.float32)
    jcsDEtaRaw = np.zeros(NC_MAX, dtype=np.float32)
    jcsDPhi = np.zeros(NC_MAX, dtype=np.float32)
    jcsM = np.zeros(NC_MAX, dtype=np.float32)
    jcsW = np.zeros(NC_MAX, dtype=np.float32)

    if nc_raw == 0:
        return -99.0, -99.0, -99.0, nc_raw, (jcsPt, jcsDEta, jcsDEtaRaw, jcsDPhi, jcsM, jcsW)

    order = np.argsort(-pt)
    if len(order) > NC_MAX:
        order = order[:NC_MAX]
    pt_s, y_s, phi_s, m_s = pt[order], y[order], phi[order], m[order]

    dY = y_s - jet_y
    dPhi = delta_phi(phi_s, jet_phi)
    r = np.sqrt(dY * dY + dPhi * dPhi)
    w = (pt_s / jet_pt) ** PV_A * r ** PV_B

    pv0 = float(np.sum(w * dY))
    pv1 = float(np.sum(w * dPhi))
    pvm = float(np.hypot(pv0, pv1))
    pva = spva = -99.0
    if pvm > 0:
        pva = float(np.arctan2(pv1 / pvm, pv0 / pvm))
        spva = float(np.arctan2(pv1 / pvm, -pv0 / pvm)) if jet_y < 0 else pva

    eta_s = derive_eta(pt_s, y_s, m_s)
    dEta_raw = eta_s - jet_eta
    dEtas = -dEta_raw if jet_eta < 0 else dEta_raw

    n = len(order)
    jcsPt[:n] = pt_s
    jcsDEta[:n] = dEtas
    jcsDEtaRaw[:n] = dEta_raw
    jcsDPhi[:n] = dPhi
    jcsM[:n] = m_s
    jcsW[:n] = w

    return pvm, pva, spva, nc_raw, (jcsPt, jcsDEta, jcsDEtaRaw, jcsDPhi, jcsM, jcsW)


def convert_file(path, event_counter_start):
    """Process one MOD h5 file -> dict of per-dijet-event numpy arrays, plus the
    updated running event counter (used for the deterministic slot-swap, mirroring
    makefriends.cpp's `swap = (iev % 2 == 0)`)."""
    with h5py.File(path, "r") as f:
        jets_i = f["jets_i"][:]
        jets_f = f["jets_f"][:]
        pfcs = f["pfcs"][:]
        pfcs_index = f["pfcs_index"][:]
        cols_i = [c.decode() for c in f["jets_i"].attrs["cols"]]
        cols_f = [c.decode() for c in f["jets_f"].attrs["cols"]]

    i_eta_j = cols_i.index("quality")
    f_pt, f_y, f_phi, f_m, f_eta = (cols_f.index(c) for c in ("jet_pt", "jet_y", "jet_phi", "jet_m", "jet_eta"))
    f_w = cols_f.index("weight")
    p_pt, p_y, p_phi, p_m = 0, 1, 2, 3  # pfcs cols are always [pt, y, phi, m, pid, vertex]

    out = {k: [] for k in SCALAR_JETS + CONSTIT_VARS + SCALAR_EVT + ["leadJetIndex"]}
    n1, n2 = 0, 0  # diagnostics: events with 1 vs 2 raw jet rows
    n_cutfail = 0

    iev = event_counter_start
    for idx in group_events(jets_i, cols_i):
        if len(idx) == 1:
            n1 += 1
            continue
        n2 += 1
        iev += 1
        # pt-order the two rows (leading first), matching makefriends.cpp's
        # pt-sorted jet list before any slot randomization.
        idx = idx[np.argsort(-jets_f[idx, f_pt])]
        eta_vals = jets_f[idx, f_eta]
        qual_vals = jets_i[idx, i_eta_j]
        if not (np.all(np.abs(eta_vals) < JET_ETA_MAX) and np.all(qual_vals >= QUALITY_MIN)):
            n_cutfail += 1
            continue

        jet_pt = jets_f[idx, f_pt]
        jet_y = jets_f[idx, f_y]
        jet_phi = jets_f[idx, f_phi]
        jet_m = jets_f[idx, f_m]
        jet_eta = eta_vals

        px, py, pz, e = p4_from_ptymphi(jet_pt, jet_y, jet_phi, jet_m)
        px_sum, py_sum, pz_sum, e_sum = px.sum(), py.sum(), pz.sum(), e.sum()
        mjj = float(np.sqrt(max(e_sum ** 2 - px_sum ** 2 - py_sum ** 2 - pz_sum ** 2, 0.0)))
        dYjj = float(jet_y[0] - jet_y[1])
        dPhijj = float(delta_phi(jet_phi[0], jet_phi[1]))
        ptjj = float(np.hypot(px_sum, py_sum))

        swap = (iev % 2 == 0)
        slot = [1, 0] if swap else [0, 1]

        pvm2, pva2, spva2, nc2, jcs2 = [None, None], [None, None], [None, None], [None, None], [None, None]
        for i in range(2):
            row = idx[i]
            c0, c1 = pfcs_index[row], pfcs_index[row + 1]
            pt_c = pfcs[c0:c1, p_pt]
            y_c = pfcs[c0:c1, p_y]
            phi_c = pfcs[c0:c1, p_phi]
            m_c = pfcs[c0:c1, p_m]
            pvm, pva, spva, nc, jcs = process_jet(pt_c, y_c, phi_c, m_c,
                                                   jet_pt[i], jet_y[i], jet_phi[i], jet_eta[i])
            si = slot[i]
            pvm2[si], pva2[si], spva2[si], nc2[si], jcs2[si] = pvm, pva, spva, nc, jcs

        jetPt2, jetEta2, jetPhi2, jetM2 = [None, None], [None, None], [None, None], [None, None]
        for i in range(2):
            si = slot[i]
            jetPt2[si] = jet_pt[i]
            jetEta2[si] = jet_eta[i]
            jetPhi2[si] = jet_phi[i]
            jetM2[si] = jet_m[i]

        out["jetPt"].append(jetPt2)
        out["jetEta"].append(jetEta2)
        out["jetPhi"].append(jetPhi2)
        out["jetM"].append(jetM2)
        out["jetSPVA"].append(spva2)
        out["jetPVA"].append(pva2)
        out["jetPVM"].append(pvm2)
        out["jetNC"].append(nc2)
        for vi, var in enumerate(CONSTIT_VARS):
            out[var].append([jcs2[0][vi], jcs2[1][vi]])
        out["kWeight"].append(float(jets_f[idx[0], f_w]))
        out["mjj"].append(mjj)
        out["dYjj"].append(dYjj)
        out["dPhijj"].append(dPhijj)
        out["ptjj"].append(ptjj)
        out["bosonM"].append(-99.0)
        out["bosonY"].append(-99.0)
        out["bosonEta"].append(-99.0)
        out["bosonPt"].append(-99.0)
        out["bosonPhi"].append(-99.0)
        out["leadJetIndex"].append(1 if swap else 0)

    return out, iev, n1, n2, n_cutfail


def finalize_arrays(out):
    """Turn one file's worth of python lists into numpy arrays."""
    arrays = {}
    for var in SCALAR_JETS:
        arrays[var] = np.asarray(out[var], dtype=np.float32 if var != "jetNC" else np.int32)
    for var in CONSTIT_VARS:
        arrays[var] = np.asarray(out[var], dtype=np.float32)
    for var in SCALAR_EVT:
        arrays[var] = np.asarray(out[var], dtype=np.float32)
    arrays["leadJetIndex"] = np.asarray(out["leadJetIndex"], dtype=np.int32)
    return arrays


ALL_VARS = SCALAR_JETS + CONSTIT_VARS + SCALAR_EVT + ["leadJetIndex"]


def append_h5(h, arrays, n_new):
    """Append one file's worth of arrays to (resizable) datasets in an open h5py
    File, creating them with maxshape=(None, ...) on the first call. Keeps peak
    memory bounded to a single input file's events, instead of a whole sample --
    a full sample can be several million events, which blew the container's
    memory cgroup limit (~35GB) when buffered whole in Python lists."""
    for var in ALL_VARS:
        data = arrays[var]
        if var not in h:
            maxshape = (None,) + data.shape[1:]
            h.create_dataset(var, data=data, maxshape=maxshape,
                              chunks=True, compression="gzip", compression_opts=4)
        else:
            dset = h[var]
            dset.resize(dset.shape[0] + n_new, axis=0)
            dset[-n_new:] = data


def convert_sample(name, subdir, limit_files=None):
    files = sorted(glob.glob(os.path.join(RAW_DIR, subdir, "*.h5")))
    # GEN*_compressed.h5 files (truth-only, no pfcs/quality/npv) are not usable by
    # this converter; only SIM*/CMS_Jet300* detector-level files are.
    files = [f for f in files if not os.path.basename(f).startswith("GEN")]
    if limit_files:
        files = files[:limit_files]
    if not files:
        print(f"[{name}] no input files found in {subdir}, skipping")
        return
    print(f"[{name}] {len(files)} input files")

    os.makedirs(OUT_DIR, exist_ok=True)
    out_path = os.path.join(OUT_DIR, f"{name}.h5")
    tmp_path = out_path + ".tmp"

    iev = 0
    tot1 = tot2 = totcutfail = totkept = 0
    with h5py.File(tmp_path, "w") as h:
        for i, path in enumerate(files):
            out, iev, n1, n2, ncf = convert_file(path, iev)
            n_new = len(out["mjj"])
            if n_new > 0:
                arrays = finalize_arrays(out)
                append_h5(h, arrays, n_new)
            del out
            tot1 += n1
            tot2 += n2
            totcutfail += ncf
            totkept += n_new
            print(f"  ({i+1}/{len(files)}) {os.path.basename(path)}: "
                  f"1-jet-events={n1} 2-jet-events={n2} cutfail={ncf} kept_so_far={totkept}")

    os.replace(tmp_path, out_path)
    print(f"[{name}] wrote {out_path} ({os.path.getsize(out_path)/1e6:.1f} MB), "
          f"{totkept} events "
          f"(1-jet dropped={tot1}, 2-jet raw={tot2}, eta/quality cutfail={totcutfail})")


def main():
    names = sys.argv[1:] if len(sys.argv) > 1 else list(SAMPLES.keys())
    for name in names:
        convert_sample(name, SAMPLES[name])


if __name__ == "__main__":
    main()
