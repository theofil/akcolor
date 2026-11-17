#!/usr/bin/env python3
"""Event displays of the nov26 beam-connection (bc) regions, and a Python reference
implementation of bc used to cross-check makefriends_nov26.

    python3 nov26/bc_display.py                         # -> figs/nov26/bc_regions.{pdf,png}
    python3 nov26/bc_display.py --friends fast_tmp_nov26 --sample QCDZjj_herwig

Draws two events in the (y, phi) plane: one with a forward jet, and one with a
negative-rapidity jet whose |dphi| < 0.4 strip crosses phi = 0/2pi. The events are
re-clustered here with FastJet (through PyROOT/cling: same library and settings as
makefriends_nov26) and bc is recomputed with numpy, independently of the C++ code.
Each panel prints the Python values next to those stored in the friend tree.

bc of a jet (y_j, phi_j), over the clustering particles (no H/nu/tau, pT > 0.5 GeV,
|eta| < 4.7), with s = sign(y_j) (+1 for y_j = 0):
    region R : s*(y_i - y_j) > 0.4  and  |y_i| < 4.7
    bc_in    = sum pT_i in R with |dphi(i, j)| < 0.4      bc_out = sum pT_i over the rest of R
    bc_met   = |sum of pT vectors in R|                   bc_dphi = Phi_mpi_pi(phi_j - phi(sum)),
                                                          -99 if R is empty
"""
import argparse
import os
import subprocess
import sys

import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # makefriends/C

CAMPAIGN23 = ('/eos/user/a/apapaefs/Projects/Pull/summer2026-hwsim-campaign-runner/'
              'campaigns/campaign-00023/merged')
PYDIPOLE = '/eos/user/a/apapaefs/Projects/Pull/pythia-dipole-recoil/merged'
INPUTS = {   # same inputs as run.nov26.sh
    'VBFZ_herwig': f'{CAMPAIGN23}/VBFZ_herwig.root',
    'VBFZ_mg5_pythia': f'{PYDIPOLE}/VBFZ/VBFZ_mg5_pythia_dipolerecoil.root',
    'QCDZjj_herwig': f'{CAMPAIGN23}/QCDZjj_herwig.root',
    'QCDZjj_mg5_pythia': f'{PYDIPOLE}/QCDZjj/QCDZjj_mg5_pythia_dipolerecoil.root',
}

# must match makefriends_nov26.cpp
EXCLUDED_PIDS = (25, 12, 14, 16, 15)
JC_PT_MIN, JC_ETA_MAX = 0.5, 4.7
JET_R, JET_PT_MIN, JET_ETA_MAX = 0.4, 30.0, 3.0
BC_DY, BC_DPHI, BC_YMAX = 0.4, 0.4, 4.7
TWOPI = 2 * np.pi


# ── physics: kinematics, particle selection, bc ──────────────────────────────

def kinematics(E, px, py, pz):
    """pt, eta, y, phi of each particle, with FastJet's PseudoJet formulas (phi in [0, 2pi))."""
    kt2 = px * px + py * py
    pt = np.sqrt(kt2)
    phi = np.where(kt2 > 0, np.arctan2(py, px), 0.0)
    phi = np.where(phi < 0, phi + TWOPI, phi)
    phi = np.where(phi >= TWOPI, phi - TWOPI, phi)
    m2 = np.maximum((E + pz) * (E - pz) - kt2, 0.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        y = 0.5 * np.log((kt2 + m2) / (E + np.abs(pz)) ** 2)
        y = np.where(pz > 0, -y, y)
        theta = np.arctan(pt / pz)
        theta = np.where(theta < 0, theta + np.pi, theta)
        eta = np.where(pz == 0, 0.0, -np.log(np.tan(theta / 2)))
    return pt, eta, y, phi


def select(pid, pt, eta):
    """Clustering particles: no H, nu, tau; pT > 0.5 GeV; |eta| < 4.7."""
    return ~np.isin(np.abs(pid), EXCLUDED_PIDS) & (pt > JC_PT_MIN) & (np.abs(eta) < JC_ETA_MAX)


def phi_mpi_pi(x):
    """ROOT's TVector2::Phi_mpi_pi: wrap an angle into [-pi, pi)."""
    return (np.asarray(x) + np.pi) % TWOPI - np.pi


def bc_masks(y, phi, yj, phij):
    """(in, out) particle masks: R = s*(y - y_j) > 0.4 and |y| < 4.7, split by |dphi| < 0.4."""
    s = -1.0 if yj < 0 else 1.0
    in_r = (s * (y - yj) > BC_DY) & (np.abs(y) < BC_YMAX)
    strip = np.abs(phi_mpi_pi(phi - phij)) < BC_DPHI
    return in_r & strip, in_r & ~strip


def bc_values(px, py, pt, y, phi, yj, phij):
    """(bc_in, bc_out, bc_met, bc_dphi, phi of the vector sum) for one jet."""
    m_in, m_out = bc_masks(y, phi, yj, phij)
    in_r = m_in | m_out
    sx, sy = px[in_r].sum(), py[in_r].sum()
    phisum = float(np.arctan2(sy, sx))
    dphi = float(phi_mpi_pi(phij - phisum)) if in_r.any() else -99.0
    return float(pt[m_in].sum()), float(pt[m_out].sum()), float(np.hypot(sx, sy)), dphi, phisum


# ── FastJet through PyROOT/cling ─────────────────────────────────────────────

_CXX = r'''
#include "fastjet/ClusterSequence.hh"
std::vector<double> nov26_antikt(int n, const double* px, const double* py, const double* pz,
                                 const double* e, double R, double ptmin, double etamax) {
    static bool quiet = (fastjet::ClusterSequence::set_fastjet_banner_stream(nullptr), true);
    (void)quiet;
    std::vector<fastjet::PseudoJet> in;
    in.reserve(n);
    for (int i = 0; i < n; ++i) in.emplace_back(px[i], py[i], pz[i], e[i]);
    fastjet::ClusterSequence cs(in, fastjet::JetDefinition(fastjet::antikt_algorithm, R));
    std::vector<double> out;
    for (const auto& j : fastjet::sorted_by_pt(cs.inclusive_jets(ptmin)))
        if (std::abs(j.eta()) < etamax)
            for (double v : {j.pt(), j.eta(), j.phi(), j.rapidity(), j.m()}) out.push_back(v);
    return out;
}
'''
_antikt = None


def antikt(px, py, pz, E):
    """anti-kT R=0.4 jets, pT > 30 GeV, |eta| < 3, pT-ordered: rows of (pt, eta, phi, y, m)."""
    global _antikt
    if _antikt is None:
        import ROOT
        prefix = subprocess.check_output(['fastjet-config', '--prefix']).decode().strip()
        ROOT.gInterpreter.AddIncludePath(prefix + '/include')
        ROOT.gSystem.Load(prefix + '/lib/libfastjet.so')
        if not ROOT.gInterpreter.Declare(_CXX):
            sys.exit('could not compile the FastJet helper')
        _antikt = ROOT.nov26_antikt
    a = [np.ascontiguousarray(v, dtype=np.float64) for v in (px, py, pz, E)]
    v = _antikt(len(a[0]), *a, JET_R, JET_PT_MIN, JET_ETA_MAX)
    return np.fromiter(v, dtype=float, count=v.size()).reshape(-1, 5)


# ── events ───────────────────────────────────────────────────────────────────

def iter_events(path, entry_start, entry_stop, step=200):
    """Yield (entry, E, px, py, pz, pid) from the input Data tree."""
    import uproot
    tree = uproot.open(path)['Data']
    entry = entry_start
    for arr in tree.iterate(['numparticles', 'objects'], entry_start=entry_start,
                            entry_stop=entry_stop, step_size=step, library='np'):
        for n, obj in zip(arr['numparticles'], arr['objects']):
            yield entry, obj[0, :n], obj[1, :n], obj[2, :n], obj[3, :n], obj[4, :n].astype(int)
            entry += 1


def reconstruct(E, px, py, pz, pid):
    """Clustering particles (dict of arrays) and the pT-ordered jets, as makefriends_nov26 does."""
    pt, eta, y, phi = kinematics(E, px, py, pz)
    sel = select(pid, pt, eta)
    parts = dict(px=px[sel], py=py[sel], pt=pt[sel], y=y[sel], phi=phi[sel])
    return parts, antikt(px[sel], py[sel], pz[sel], E[sel])


def pt_order_index(slot, n_jets, lead_slot):
    """pT-order index (0 leading, 1 subleading) of the jet stored in `slot`."""
    if n_jets == 1:
        return 0
    return 0 if slot == lead_slot else 1


# ── display ──────────────────────────────────────────────────────────────────

C_IN, C_OUT, C_SUM = '#2a78d6', '#eb6834', '#1baf7a'          # categorical slots 1-3
C_GRAY, INK, INK2, RULE, SURFACE = '#a8a7a1', '#0b0b0b', '#52514e', '#c3c2b7', '#fcfcfb'

FRIEND_COLS = (['evIdx', 'gate', 'nJets', 'leadJetIndex']
               + [f'{v}{s}' for v in ('jetPt', 'jetEta', 'jetPhi', 'jetBCin', 'jetBCout',
                                      'jetBCmet', 'jetBCdphi') for s in (0, 1)])


def pick(friend):
    """(evIdx, slot) of a forward-jet event and of a negative-rapidity-jet event whose
    |dphi| < 0.4 strip crosses phi = 0/2pi; both gate events with populated bc regions."""
    best = {}
    for slot in (0, 1):
        filled = friend['gate'] & (slot < friend['nJets'])
        eta, phi = friend[f'jetEta{slot}'], friend[f'jetPhi{slot}']
        bin_, bout = friend[f'jetBCin{slot}'], friend[f'jetBCout{slot}']
        busy = (bin_ > 3) & (bout > 10)
        wraps = (phi < 0.3) | (phi > TWOPI - 0.3)
        for key, mask in (('forward', filled & busy & (eta > 1.5)),
                          ('negative', filled & busy & (eta < -1.5) & wraps)):
            idx = np.flatnonzero(mask)
            if len(idx) and (key not in best or idx[0] < best[key][0]):
                best[key] = (int(idx[0]), slot)
    missing = {'forward', 'negative'} - set(best)
    if missing:
        sys.exit(f'no suitable event found for: {sorted(missing)}')
    return best['forward'], best['negative']


def periodic_bands(center, half):
    """[lo, hi] intervals of center +- half in [0, 2pi], split where they cross 0/2pi."""
    lo, hi = center - half, center + half
    out = [(max(lo, 0.0), min(hi, TWOPI))]
    if lo < 0:
        out.append((lo + TWOPI, TWOPI))
    if hi > TWOPI:
        out.append((0.0, hi - TWOPI))
    return out


def draw(ax, parts, jet, other, bc_py, bc_cpp, label):
    import matplotlib.patches as mpatches
    yj, phij = jet[3], jet[2]
    s = -1.0 if yj < 0 else 1.0
    y0, y1 = sorted((yj + s * BC_DY, s * BC_YMAX))           # region R in y
    # regions: R tinted as bc_out, the |dphi| < 0.4 strip as bc_in
    ax.axvspan(y0, y1, color=C_OUT, alpha=0.08, lw=0, zorder=0)
    for lo, hi in periodic_bands(phij, BC_DPHI):
        ax.fill_between([y0, y1], lo, hi, color=C_IN, alpha=0.16, lw=0, zorder=0)
        ax.hlines([v for v in (lo, hi) if 0 < v < TWOPI], y0, y1, colors=C_IN, lw=1, ls='--', zorder=1)
    ax.axvline(yj + s * BC_DY, color=INK2, lw=1, ls='--', zorder=1)
    for edge in (-BC_YMAX, BC_YMAX):
        ax.axvline(edge, color=INK2, lw=1, ls=':', zorder=1)
    # particles
    m_in, m_out = bc_masks(parts['y'], parts['phi'], yj, phij)
    rest = ~(m_in | m_out)
    size = np.clip(5 + 4 * parts['pt'], 5, 220)
    for mask, color, z in ((rest, C_GRAY, 2), (m_out, C_OUT, 3), (m_in, C_IN, 4)):
        ax.scatter(parts['y'][mask], parts['phi'][mask], s=size[mask], c=color,
                   edgecolors=SURFACE, linewidths=0.4, zorder=z)
    # jets: axis + R = 0.4 cone (drawn periodically)
    for j, color, ls, ms in ((jet, INK, '-', 15), (other, INK2, ':', 11)):
        if j is None:
            continue
        ax.plot(j[3], j[2], marker='*', ms=ms, color=color, mec=SURFACE, mew=0.8, zorder=6)
        for shift in (-TWOPI, 0.0, TWOPI):
            ax.add_patch(mpatches.Circle((j[3], j[2] + shift), JET_R, fill=False, ec=color,
                                         lw=1.3, ls=ls, zorder=5))
    # values
    b_in, b_out, b_met, b_dphi, phisum = bc_py
    c_in, c_out, c_met, c_dphi = bc_cpp
    left = s > 0                                              # region on the right -> text left
    ax.text(0.015 if left else 0.985, 0.975,
            f'{label}\njet:  y = {yj:+.3f},  φ = {phij:.3f},  pT = {jet[0]:.1f} GeV\n'
            f'bc_in = {b_in:.2f} GeV    bc_out = {b_out:.2f} GeV\n'
            f'bc_met = {b_met:.2f} GeV    bc_dphi = {b_dphi:+.3f}\n'
            f'makefriends_nov26: {c_in:.2f} / {c_out:.2f} / {c_met:.2f} / {c_dphi:+.3f}',
            transform=ax.transAxes, ha='left' if left else 'right', va='top', fontsize=8.5,
            color=INK, zorder=8, linespacing=1.35,
            bbox=dict(boxstyle='round,pad=0.4', fc=SURFACE, ec=RULE, alpha=0.93))
    # transverse-plane inset: jet direction and direction of the vector sum over R
    ins = ax.inset_axes([0.01 if left else 0.79, 0.03, 0.2, 0.3], projection='polar', zorder=9)
    ins.set_facecolor(SURFACE)
    ins.annotate('', xy=(phij, 1.0), xytext=(0, 0),
                 arrowprops=dict(arrowstyle='-|>', color=INK, lw=1.6))
    ins.annotate('', xy=(phisum, 0.8), xytext=(0, 0),
                 arrowprops=dict(arrowstyle='-|>', color=C_SUM, lw=2.2))
    ins.text(phij, 1.32, 'jet', ha='center', va='center', fontsize=7.5, color=INK)
    ins.text(phisum, 1.12, 'Σ pT(R)', ha='center', va='center', fontsize=7.5, color=INK)
    ins.set_ylim(0, 1.45)
    ins.set_yticks([])
    ins.set_xticks(np.arange(4) * np.pi / 2, ['0', 'π/2', 'π', '3π/2'], fontsize=6.5, color=INK2)
    ins.grid(color=RULE, lw=0.5)
    ins.set_title(f'transverse plane, Δφ = {b_dphi:+.2f}', fontsize=7, color=INK2, pad=11)
    # axes
    ax.set_xlim(-5, 5)
    ax.set_ylim(0, TWOPI)
    ax.set_aspect('equal')
    ax.set_xlabel('rapidity y', color=INK)
    ax.set_ylabel('φ [rad]', color=INK)
    ax.set_yticks(np.arange(5) * np.pi / 2, ['0', 'π/2', 'π', '3π/2', '2π'])
    ax.tick_params(colors=INK2, labelsize=9)
    for sp in ax.spines.values():
        sp.set_color(RULE)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--sample', default='VBFZ_herwig', choices=sorted(INPUTS))
    ap.add_argument('--friends', default=None,
                    help='folder with <sample>.friend.root (default friends/nov26, else fast_tmp_nov26)')
    ap.add_argument('--out', default=os.path.join(BASE, 'figs', 'nov26', 'bc_regions'))
    args = ap.parse_args()

    import uproot
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    folders = [args.friends] if args.friends else [os.path.join(BASE, 'friends', 'nov26'),
                                                    os.path.join(BASE, 'fast_tmp_nov26')]
    friend_path = next((p for p in (os.path.join(f, f'{args.sample}.friend.root') for f in folders)
                        if os.path.exists(p)), None)
    if friend_path is None:
        sys.exit(f'{args.sample}.friend.root not found in {folders}')
    friend = uproot.open(friend_path)['events'].arrays(FRIEND_COLS, library='np')

    fig, axes = plt.subplots(2, 1, figsize=(9.2, 12.2), facecolor=SURFACE)
    for ax, (row, slot), kind in zip(axes, pick(friend), ('Forward jet', 'Negative-rapidity jet')):
        ev = int(friend['evIdx'][row])
        _, E, px, py, pz, pid = next(iter_events(INPUTS[args.sample], ev, ev + 1))
        parts, jets = reconstruct(E, px, py, pz, pid)
        n_jets, lead = int(friend['nJets'][row]), int(friend['leadJetIndex'][row])
        k = pt_order_index(slot, n_jets, lead)
        jet = jets[k]
        if abs(round(jet[0], 2) - friend[f'jetPt{slot}'][row]) > 0.011:
            sys.exit(f'event {ev}: re-clustered jet pT {jet[0]:.3f} != friend {friend[f"jetPt{slot}"][row]}')
        other = jets[1 - k] if n_jets == 2 else None
        bc_py = bc_values(parts['px'], parts['py'], parts['pt'], parts['y'], parts['phi'], jet[3], jet[2])
        bc_cpp = tuple(float(friend[f'{v}{slot}'][row]) for v in ('jetBCin', 'jetBCout', 'jetBCmet', 'jetBCdphi'))
        draw(ax, parts, jet, other, bc_py, bc_cpp, f'{args.sample}, event {ev}')
        side = '+y' if jet[3] >= 0 else '−y'
        extra = ', the |Δφ| < 0.4 strip crosses φ = 0/2π' if kind.startswith('Negative') else ''
        ax.set_title(f'{kind}: region R runs from y_j {"+" if jet[3] >= 0 else "−"} 0.4 toward the {side} beam{extra}',
                     loc='left', fontsize=10, color=INK)
        print(f'{kind:22s} event {ev} slot {slot}: y_j={jet[3]:+.4f} phi_j={jet[2]:.4f}  '
              f'python bc = {bc_py[0]:.4f} {bc_py[1]:.4f} {bc_py[2]:.4f} {bc_py[3]:+.5f}  '
              f'friend = {bc_cpp[0]:.4f} {bc_cpp[1]:.4f} {bc_cpp[2]:.4f} {bc_cpp[3]:+.5f}')

    handles = [
        Line2D([], [], ls='', marker='o', ms=7, mfc=C_IN, mec=SURFACE, label='particle in bc_in (R, |Δφ| < 0.4)'),
        Line2D([], [], ls='', marker='o', ms=7, mfc=C_OUT, mec=SURFACE, label='particle in bc_out (R, |Δφ| ≥ 0.4)'),
        Line2D([], [], ls='', marker='o', ms=7, mfc=C_GRAY, mec=SURFACE, label='particle outside R'),
        Line2D([], [], ls='-', marker='*', ms=12, color=INK, label='jet axis and R = 0.4 cone'),
        Line2D([], [], ls=':', marker='*', ms=9, color=INK2, label='other leading jet'),
        Line2D([], [], ls='--', color=INK2, label='start of R: y = y_j ± 0.4'),
        Line2D([], [], ls=':', color=INK2, label='|y| = 4.7'),
        Line2D([], [], ls='--', color=C_IN, label='strip edges |Δφ| = 0.4'),
    ]
    fig.suptitle('nov26 beam connection regions. Particles: pT > 0.5 GeV, |η| < 4.7 '
                 '(marker area grows with pT)', y=0.995, fontsize=10, color=INK)
    fig.legend(handles=handles, loc='upper center', ncol=3, frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.5, 0.978))
    fig.subplots_adjust(left=0.08, right=0.98, top=0.895, bottom=0.05, hspace=0.18)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig(f'{args.out}.{ext}', dpi=150, facecolor=SURFACE)
    print(f'wrote {args.out}.pdf and .png')


if __name__ == '__main__':
    main()
