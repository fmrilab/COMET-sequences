"""Generate the 20261112comet sequences (collaborator request of 2026-11-11).

Protocol: 90 x 90 x 60 matrix, 2.4 mm isotropic, FOV 216 x 216 x 144 mm (ArbEPI-python's defaults),
GE_MR750, water excitation, human PNS weights, every file under the 80 % normal-mode PNS limit.
TE 30 ms (arm 1: 33.3 ms) is the time from the RF center to the kx = ky = kz = 0 sample (see "TE" below).
Each EPI scan is 5 min (volume TR is the fastest the arm allows: Nframes = ceil(300 s / volume TR)).

Folders (numbered in scan order, like 20260922xiaokai/seqs):
  1_laminar_1x_2.4mm_ETL90   fully sampled, laminar, 60 shots x ETL 90: every train is one full ky line
                             (one kz, ky = 0 .. 89 in unit steps, same direction every shot), kz ascending;
                             4 us dwell forced, TE 33.3 ms (user choices)
  2_radial_90x_2.4mm_ETL60     mask_60.mat (60 samples) in ONE echo train per frame (1 shot x ETL 60), radial
  3_radial_122.7x_2.4mm_ETL44  mask_44.mat (44 samples) in ONE echo train per frame (1 shot x ETL 44), radial
  4_meGRE_dwell4us           monopolar multi-echo 3D GRE (sequences/meGRE.py), same FOV / resolution as the
                             EPI and the receiver bandwidth (ADC dwell, 4 us) of the EPI. Echoes at multiples
                             of the fat period (1/f_fat = 2.236 ms), the first at 1 x, the last at the
                             latest multiple <= 50 ms (22 echoes).
Folders 1-3 hold ArbEPI, EPIcal, noise (.seq/.pge), scan_info.mat and diagnostic plots; 4-5 hold meGRE.

Things this script does that ArbEPI-python's main.py does not:
  * Laminar ordering: ArbEPI-python's mask2epi_laminar fills shots with ETL consecutive samples of a
    (kz outer, ky center-out inner) raster, which is not one monotonic ky line per train, so arm 1 swaps in
    `mask2epi_ky_lines` through sequences.ArbEPI._MASK2EPI (no repo edit). Arm 1 also forces the 4 us dwell
    over find_min_feasible_dwell's 6 us (acoustic-band margin) and prescribes TE 33.3 ms (30 ms is not
    reachable with ky = 0 on echo 45 at ETL 90 and 6 us; kept at the user's 33.3 ms with 4 us).
  * One train per frame (radial arms): the repo's mask2epi_radial is run once in its per-frame mode on the
    provided mask with Nshots = 1 (one spoke through the center; the repo's golden-angle chaining across
    frames is bypassed through sequences.ArbEPI._compute_schedules), and that schedule is played on every
    frame: same echo order every frame. Neither provided mask contains the k-space center (ky, kz) = (45, 30)
    (0-based; the index ArbEPI's prephasers/blips treat as k = 0), so `static_radial` replaces the train's
    sample nearest the center (echo (ETL - 1) // 2) by DC: the same one-sample move onto DC as the user
    chose on 2026-11-11; sample count, ETL and R unchanged.
  * TE: that DC echo is (ETL - 1) // 2. ArbEPI-python's TE anchors continuous echo index
    ETL/2 - 0.5, half an echo spacing away for even ETL, so the prescribed TE is chosen (`design`) such that
    the DC sample's own kx = 0 crossing is at 30.000 ms to within the 4 us block raster, and checked on a
    built snippet's k-space trajectory in the measured train (`verify_te`).
  * Timing is measured, not assumed: min TE/min TR come from lib/epi_shot.build_epi_shot_parts on the
    schedule, and volume TR = Nshots x (min TR per shot rounded up to the 4 us raster).

Run from an ArbEPI-python checkout on main with an unmodified params.py:
    cd <ArbEPI-python checkout, main at 28a28fa>
    .venv/bin/python <this folder>/generate.py measure [arm ...]   # timing only, writes nothing
    .venv/bin/python <this folder>/pns_search.py <radial arm>      # slew choice for SLEWS below
    .venv/bin/python <this folder>/generate.py build [arm ...]     # folders 1-3 (ArbEPI, EPIcal, noise, plots)
    .venv/bin/python <this folder>/generate.py megre               # folder 4
    .venv/bin/python -I <this folder>/verify.py                    # checks scan_info.mat
Arms: laminar_1x, radial_mask60, radial_mask44.
"""

import dataclasses
import glob
import json
import math
import os
import shutil
import sys
import tempfile

import numpy as np
import scipy.io as sio

sys.path.insert(0, os.getcwd())
import sequences.ArbEPI as arb  # noqa: E402
from params import load_params  # noqa: E402

ROOT = os.path.dirname(os.path.abspath(__file__))  # this folder: masks in, sequence folders out
DURATION = 300.0  # s per EPI run
TE_TARGET = 30e-3  # s, RF center -> kx = ky = kz = 0 sample (radial arms; per-arm 'te' in ARMS overrides)
# arm -> ADC dwell (s) forced over find_min_feasible_dwell's choice. Arm 1 at 4 us has a 528 us echo spacing, 18 us
# above the coil's 410-510 us forbidden acoustic band, inside the repo's 20 us margin (user choice).
FORCE_DWELL = {'laminar_1x': 4e-6}
RASTER = 4e-6  # s, GE block raster
T1 = 1.3  # s, load_params()'s `T1` (Ernst angle)
PNS_LIMIT = 80.0  # %, GE normal operating mode; build refuses to export above it
MEGRE_TE_MAX = 50e-3  # s, last meGRE echo
DC = (45, 30)  # 0-based (ky, kz) index of k = 0 in ArbEPI (Ny/2, Nz/2)

# arm -> folder, trajectory, ETL, shots per frame, mask file (None = fully sampled). The radial arms split the
# provided mask into one echo train per frame (1 shot, ETL = number of samples, `static_radial`); folder names
# carry R (Ny*Nz / samples per frame).
ARMS = {
    'laminar_1x': dict(folder='1_laminar_1x_2.4mm_ETL90', traj='laminar', ETL=90, Nshots=60, mask=None, te=33.3e-3),
    'radial_mask60': dict(folder='2_radial_90x_2.4mm_ETL60', traj='radial', ETL=60, Nshots=1, src_shots=1, mask='mask_60'),
    'radial_mask44': dict(folder='3_radial_122.7x_2.4mm_ETL44', traj='radial', ETL=44, Nshots=1, src_shots=1, mask='mask_44'),
}
# arm -> blip / readout-fall / readout-rise slew (T/m/s) overrides of params.py's defaults (200 / 190 / 155).
# At the defaults the single-train radial arms are over the 80 % normal-mode PNS limit or too close to it. The
# triples below are among the fastest per-shot TRs on pns_search.py's coarse grid with PNS <= 79 % (6 frames,
# every spoiler at its maximum), picked for margin among near-ties. Arm 1 (forced 4 us dwell): at the defaults
# its echo spacing is 524 us, 14 us above the coil's 510 us band edge, and the acoustics check reads 0.2995 of
# 0.3; a slower readout ramp-up (rise 120) moves the spacing to 552 us (acoustics 0.09) for +2 % TR.
SLEWS = {
    'laminar_1x': dict(blip_slew=200.0, ro_slew_fall=190.0, ro_slew_rise=120.0),
    'radial_mask60': dict(blip_slew=160.0, ro_slew_fall=160.0, ro_slew_rise=140.0),
    'radial_mask44': dict(blip_slew=180.0, ro_slew_fall=140.0, ro_slew_rise=140.0),
}


# ---------------------------------------------------------------------------------------------
# masks and trajectories
# ---------------------------------------------------------------------------------------------
def mask2epi_ky_lines(mask, ETL, Nshots, deltak=(1.0, 1.0)):
    """Fully sampled laminar: shot s is the full ky line at kz = s (ascending), ky = 0 .. Ny - 1 in unit
    steps. Same signature/returns as lib.mask2epi.mask2epi_laminar. ky always ascends, so the ky = 0
    sample (index Ny/2) is on echo Ny/2 of every train and has the same TE on every kz line."""
    Ny, Nz = mask.shape
    assert mask.all() and ETL == Ny and Nshots == Nz, 'ky-line laminar needs a full mask, ETL = Ny, Nshots = Nz'
    schedule = np.zeros((Nshots, ETL, 2), dtype=int)
    schedule[:, :, 0] = np.arange(Ny)[None, :]
    schedule[:, :, 1] = np.arange(Nz)[:, None]
    parts = np.repeat(np.arange(1, Nz + 1)[None, :], Ny, axis=0)
    return schedule, parts


def static_radial(mask, ETL, Nshots, deltak=(1 / 216e-3, 1 / 144e-3)):
    """Frame schedules from a provided radial mask, each with k = 0.

    The repo's `mask2epi_radial` (per-frame behavior, `tau0=None`) orders the provided mask into `Nshots`
    trains; train i is the whole of frame type i (Nshots = 1: the one train is the whole mask). Each train's
    sample nearest k-space center sits on echo `(ETL - 1) // 2`; it is replaced by the DC sample
    (ky, kz) = (45, 30), so every frame plays k = 0 once, on the same echo (same TE in every frame). Returns
    (schedule (Nshots, ETL, 2) 0-based, frame masks (Ny, Nz, Nshots) bool, list of the replaced (ky, kz)
    samples)."""
    from lib.mask2epi import mask2epi_radial

    schedule, _ = mask2epi_radial(mask, ETL, Nshots, deltak)  # tau0=None: no chaining across frames
    t = (ETL - 1) // 2
    removed = [tuple(int(v) for v in schedule[s, t]) for s in range(Nshots)]
    schedule = schedule.copy()
    schedule[:, t] = DC
    eff = np.zeros(mask.shape + (Nshots,), dtype=bool)
    for s in range(Nshots):
        eff[schedule[s, :, 0], schedule[s, :, 1], s] = True
        assert eff[..., s].sum() == ETL
    return schedule, eff, removed


_STATIC = {}  # arm -> (schedule, frame masks, replaced samples)


def load_arm_mask(arm):
    """(Ny, Nz, n_types) bool frame masks of an arm, repeated cyclically over the frames: all ones (arm 1, one
    type) or the two one-shot frames of the provided mask with DC in each (`static_radial`)."""
    a = ARMS[arm]
    if a['mask'] is None:
        return np.ones((90, 60, 1), bool)
    if arm not in _STATIC:
        src = sio.loadmat(os.path.join(ROOT, a['mask'] + '.mat'))['samp'].astype(bool)
        assert src.sum() == a['ETL'] * a['src_shots'] and not src[DC]
        _STATIC[arm] = static_radial(src, a['ETL'], a['src_shots'])
    return _STATIC[arm][1]


def use_arm(arm):
    """Install the arm's schedule builder in the repo's schedule computation (module-level patches, no repo
    edit): arm 1's ky-line laminar order and forced 4 us dwell, or the radial arm's two frame types alternating A, B, A, B, ..."""
    import lib.readout_from_params as rfp

    a = ARMS[arm]
    if not hasattr(rfp, '_orig_find_min_feasible_dwell'):
        rfp._orig_find_min_feasible_dwell = rfp.find_min_feasible_dwell
    if arm in FORCE_DWELL:  # the readout geometry is still checked by make_readout_grads
        rfp.find_min_feasible_dwell = lambda *args, **kw: FORCE_DWELL[arm]
    else:
        rfp.find_min_feasible_dwell = rfp._orig_find_min_feasible_dwell
    if a['traj'] == 'laminar':
        arb._MASK2EPI['laminar'] = mask2epi_ky_lines
        return
    load_arm_mask(arm)
    schedule = _STATIC[arm][0]  # (n_types, ETL, 2); one shot per frame

    def static_schedules(omegas, ETL, Nshots, trajectory, deltak):
        Ny, Nz, Nframes = omegas.shape
        assert Nshots == 1 and len(schedule) >= 1
        types = np.arange(Nframes) % len(schedule)
        parts = omegas.astype(int)  # 1-based shot label, 0 = unsampled: one shot per frame
        return schedule[types][:, None], parts

    arb._compute_schedules = static_schedules


# ---------------------------------------------------------------------------------------------
# parameters and timing
# ---------------------------------------------------------------------------------------------
def ernst(tr):
    return 180 / math.pi * math.acos(math.exp(-tr / T1))


def check_base(p):
    assert (p.Nx, p.Ny, p.Nz) == (90, 90, 60) and np.allclose(p.res, 2.4e-3)
    assert p.excitation == 'water' and p.discard_duration == 0 and p.spec.ge_coil == 'xrm'
    assert tuple(p.PNSwt) == (0.8, 1.0, 0.7) and p.custom_mask_path is None


def arm_params(arm, output_dir, tr_shot, te, Nframes, masks, slews=None):
    """Params of one EPI arm with per-shot TR `tr_shot`, prescribed TE `te` and `Nframes` frames."""
    a = ARMS[arm]
    p = load_params(output_dir=output_dir)
    check_base(p)
    volume_tr = round(a['Nshots'] * tr_shot, 9)
    # +0.1 ns: calc_te_tr_delays floors (TR - min_tr) / raster with no epsilon
    TR = volume_tr / a['Nshots'] + 1e-10
    assert len(set(int(masks[..., i].sum()) for i in range(masks.shape[2]))) == 1
    R = p.Ny * p.Nz / int(masks[..., 0].sum())  # per frame
    p = dataclasses.replace(
        p, R=R, ETL=a['ETL'], Nshots=a['Nshots'], volume_tr=volume_tr, TR=TR, Nframes=Nframes,
        discard_duration=0.0, fa=ernst(TR), TE=te, epi_trajectory=a['traj'], sampling_method=None,
        seed=None, custom_mask_path=None, custom_omegas=masks[:, :, np.arange(Nframes) % masks.shape[2]],
        **(SLEWS.get(arm, {}) if slews is None else slews),
    )
    return p


class Stats:
    """Blip steps and DC echo index of an arm's schedule. The schedule is periodic over the frames (period
    = number of frame types: 1, or 2 for A, B), so a few frames say it all."""

    def __init__(self, arm, masks, fov=(216e-3, 216e-3, 144e-3)):
        a = ARMS[arm]
        use_arm(arm)
        self.n_types = masks.shape[2]
        nf = 2 * self.n_types
        omegas = masks[:, :, np.arange(nf) % self.n_types]
        self.schedules, _ = arb._compute_schedules(
            omegas, a['ETL'], a['Nshots'], a['traj'], deltak=(1 / fov[1], 1 / fov[2]))
        s = self.schedules  # (frame, shot, echo, 2)
        assert (s == s[np.arange(nf) % self.n_types]).all(), 'the schedule must repeat with the frame types'
        self.ky_step = float(np.abs(np.diff(s[..., 0], axis=2)).max())
        self.kz_step = float(np.abs(np.diff(s[..., 1], axis=2)).max())
        is_dc = (s[..., 0] == DC[0]) & (s[..., 1] == DC[1])  # (frame, shot, echo)
        assert np.all(is_dc[:self.n_types].sum(axis=(1, 2)) == 1) or a['traj'] == 'laminar'
        assert is_dc[0].sum() == 1, 'k = 0 must be played exactly once per frame'
        self.dc_echo = np.unique(np.nonzero(is_dc)[2])
        assert len(self.dc_echo) == 1, f'DC lands on echoes {self.dc_echo} in different trains'
        self.dc_echo = int(self.dc_echo[0])
        # trains (frame * Nshots + shot) of the first n_types frames that play k = 0
        self.dc_shots = [int(f * a['Nshots'] + sh) for f, sh in zip(*np.nonzero(is_dc[:self.n_types].any(axis=2)))]

    def steps(self):
        return self.ky_step, self.kz_step


def design(arm, mask, stats, output_dir, slews=None):
    """Fixed point of (per-shot TR, Nframes, flip angle, TE) for an arm. Returns (params, info)."""
    import warnings

    from lib.epi_shot import build_epi_shot_parts
    import pypulseq as pp

    a = ARMS[arm]
    tr = 0.045 if a['Nshots'] <= 2 else 0.050
    te_target = a.get('te', TE_TARGET)
    te_p = te_target
    for it in range(12):
        n = math.ceil(DURATION / (a['Nshots'] * tr) - 1e-9)
        n = math.ceil(n / stats.n_types) * stats.n_types  # whole A, B cycles
        dky, dkz = stats.steps()
        p = arm_params(arm, output_dir, tr, te_p, n, mask, slews)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            parts = build_epi_shot_parts(p, dky, dkz)
        gro = pp.calc_duration(parts.rg.gro)
        offs = (stats.dc_echo - (a['ETL'] / 2 - 0.5)) * gro  # DC echo vs the nominal-TE echo position
        te_delay = max(0.0, round((te_target - offs - parts.min_te) / RASTER) * RASTER)
        te_new = parts.min_te + te_delay + RASTER / 2  # floors back to exactly te_delay
        # min TR depends on te_delay, so re-measure at the new TE before rounding TR
        p2 = dataclasses.replace(p, TE=te_new)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            parts2 = build_epi_shot_parts(p2, dky, dkz)
        assert abs(parts2.te_delay - te_delay) < 1e-12, (parts2.te_delay, te_delay)
        tr_new = math.ceil(parts2.min_tr / RASTER - 1e-6) * RASTER
        if abs(tr_new - tr) < 1e-12 and abs(te_new - te_p) < 1e-12:
            break
        tr, te_p = tr_new, te_new
    else:
        raise RuntimeError(f'{arm}: TE/TR did not converge')
    echo_time_dc = parts2.min_te + parts2.te_delay + offs
    info = dict(arm=arm, ETL=a['ETL'], Nshots=a['Nshots'], R=p2.R, Nframes=n, tr_shot=tr,
                volume_tr=a['Nshots'] * tr, fa=p2.fa, te_prescribed=te_new, min_te=parts2.min_te,
                min_tr=parts2.min_tr, dwell=parts2.rg.adc.dwell, esp=gro, dky=dky, dkz=dkz,
                dc_echo=stats.dc_echo, te_delay=parts2.te_delay, tr_delay=parts2.tr_delay,
                te_dc_nominal=echo_time_dc, scan_s=n * a['Nshots'] * tr)
    return p2, info


def measure(names=None):
    out = {}
    for arm in names or list(ARMS):
        mask = load_arm_mask(arm)
        a = ARMS[arm]
        stats = Stats(arm, mask)
        with tempfile.TemporaryDirectory() as tmpd:
            p, info = design(arm, mask, stats, tmpd)
        out[arm] = info
        print(json.dumps({k: (float(v) if isinstance(v, (np.floating, float)) else v)
                          for k, v in info.items()}, indent=None), flush=True)
    return out



# ---------------------------------------------------------------------------------------------
# meGRE
# ---------------------------------------------------------------------------------------------
def megre_kwargs(p, dwell):
    """design_megre/generate_megre keywords for a meGRE whose ADC dwell is `dwell` (so its receiver
    bandwidth is the EPI's) and whose last echo is the latest at or below MEGRE_TE_MAX: echoes sit
    at TE_n = n / f_fat (fat in phase with water), the first at the earliest multiple that fits the
    prephasers (1 x 2.236 ms)."""
    from sequences.meGRE import design_megre

    kw = dict(read_frac=dwell * p.Nx / (1 / p.fat_offres_freq))
    d1 = design_megre(p, n_echoes=1, **kw)
    n = int(math.floor((MEGRE_TE_MAX - d1.TE[0]) / d1.echo_spacing + 1e-9)) + 1
    d = design_megre(p, n_echoes=n, **kw)
    assert abs(d.dwell - dwell) < 1e-12, (d.dwell, dwell)
    assert d.TE[-1] <= MEGRE_TE_MAX + 1e-9 < d.TE[-1] + d.echo_spacing, d.TE
    return dict(kw, n_echoes=n), d


def pns_ok(name, rep):
    print(f'--- {name}\n{rep.summary()}', flush=True)
    assert rep.ok, f'{name} fails the GE feasibility check'
    # rep.ok gates PNS at 100 %; this is a human session and the scanner refuses > 80 % (normal mode)
    assert rep.peak_pns_percent < PNS_LIMIT, f'{name}: peak PNS {rep.peak_pns_percent:.1f}% >= {PNS_LIMIT}%'


def export(p, seq_path, base):
    from ge.ge_export import check_ge_feasibility, export_to_ge

    seq, rep = check_ge_feasibility(seq_path, p)
    pns_ok(os.path.basename(seq_path), rep)
    export_to_ge(seq_path, base, p, seq=seq, report=rep)
    return rep


def plot_megre_one_tr(seq, d, png):
    """one_tr.png of a meGRE: the first image excitation (iZ = 1, iY = 0, i.e. after the dummy z-loops and the
    receive-gain block; its phase-encode prephasers/rewinders are at their maximum), one TR from the start
    of its RF block. Dashed lines mark every echo center (TE_1, TE_n labelled), taken from the ADC times
    of the sequence and checked against the design's TE_n."""
    import matplotlib

    matplotlib.use('Agg')
    n_pre = (d.blocks['n_dummy_zloops'] + 1) * d.N[1]  # excitations before the first image one
    t_rf = np.asarray(seq.rf_times()[0]).ravel()
    t_blk = n_pre * d.TR  # RF blocks start every TR, the first one at t = 0
    assert abs(t_blk + t_rf[0] - t_rf[n_pre]) < 1e-9, 'excitations are not equally spaced by TR'
    splot = seq.plot(time_range=(t_blk, t_blk + d.TR), stacked=True, plot_now=False, time_disp='ms')
    fig = splot.fig1
    fig.set_size_inches(20, 10)
    t_adc, _ = seq.adc_times(time_range=(t_blk, t_blk + d.TR))
    n_fid = len(t_adc) // d.n_echoes
    assert n_fid * d.n_echoes == len(t_adc)
    centers = t_adc[n_fid // 2::n_fid]
    te = centers - t_rf[n_pre]
    assert np.abs(te - np.asarray(d.TE)).max() < 10e-6, (te, d.TE)
    for ax in splot.ax1:
        for c in centers:
            ax.axvline(c * 1e3, color='k', linestyle='--', linewidth=0.6, zorder=10)
    top = splot.ax1[0].get_xaxis_transform()
    for i in (0, d.n_echoes - 1):
        splot.ax1[0].text(centers[i] * 1e3, 1.02, f'TE{i + 1} {te[i] * 1e3:.2f} ms', transform=top,
                          ha='center', va='bottom', fontsize=8)
    fig.tight_layout()
    fig.savefig(png, dpi=300)


def build_megre(dwell, folder):
    """The meGRE of receiver dwell `dwell`, into ROOT/folder. Returns its design."""
    from sequences.meGRE import generate_megre

    out = os.path.join(ROOT, folder)
    os.makedirs(out, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix='build_megre_', dir=ROOT)
    p = load_params(output_dir=tmp)
    check_base(p)
    kw, _ = megre_kwargs(p, dwell)
    _, d = generate_megre(p, **kw)
    rep = export(p, os.path.join(tmp, 'meGRE.seq'), os.path.join(out, 'meGRE'))
    import pypulseq as pp

    seq = pp.Sequence(system=p.sys)
    seq.read(os.path.join(tmp, 'meGRE.seq'))
    plot_megre_one_tr(seq, d, os.path.join(out, 'one_tr.png'))
    shutil.move(os.path.join(tmp, 'meGRE.seq'), os.path.join(out, 'meGRE.seq'))
    shutil.rmtree(tmp)
    info = dict(n_echoes=d.n_echoes, TE_ms=[round(t * 1e3, 3) for t in d.TE], echo_spacing_us=d.echo_spacing * 1e6,
                dwell_us=d.dwell * 1e6, bandwidth_hz=1 / d.dwell, TR_ms=d.TR * 1e3, flip_deg=d.flip_deg,
                scan_s=d.scan_time, read_gradient_mT_m=d.g_read_mT_m, matrix=list(d.N),
                fov_mm=[f * 1e3 for f in d.fov], peak_pns_percent=rep.peak_pns_percent,
                acoustics=rep.acoustics_max_in_band)
    json.dump(info, open(os.path.join(out, 'meGRE_info.json'), 'w'), indent=1)
    return d


# ---------------------------------------------------------------------------------------------
# build + verification
# ---------------------------------------------------------------------------------------------
def verify_te(p, info, stats):
    """Build frame 0 of the arm (same timing as the full sequence: the largest blip steps of the whole
    run are forced) and measure the DC sample's kx = 0 crossing from the played k-space trajectory,
    relative to the RF center. Returns seconds."""
    import warnings

    from lib.epi_shot import build_epi_shot_parts

    tmp = tempfile.mkdtemp(prefix='verify_', dir=ROOT)
    full = (info['dky'], info['dkz'])
    nf = stats.n_types  # one cycle of frame types
    ps = dataclasses.replace(p, Nframes=nf, output_dir=tmp, custom_omegas=p.custom_omegas[..., :nf])
    orig = arb.max_blip_steps
    arb.max_blip_steps = lambda sch: full
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            seq = arb.generate_arbepi(ps.custom_omegas, ps, seqname='xcheck')
            rg = build_epi_shot_parts(ps, *full).rg
    finally:
        arb.max_blip_steps = orig
    k_traj_adc, _, t_exc, _, t_adc = seq.calculate_kspace()
    assert len(t_exc) == nf * ps.Nshots, (len(t_exc), nf, ps.Nshots)
    out = []
    for shot in stats.dc_shots:  # radial arms: every train; arm 1: the one train with ky = 0
        i0 = (shot * ps.ETL + stats.dc_echo) * rg.Nfid
        kx, t = k_traj_adc[0, i0:i0 + rg.Nfid], t_adc[i0:i0 + rg.Nfid]
        # DC in (ky, kz) too: the sample on that echo is the schedule's (45, 30)
        sc = np.nonzero(np.diff(np.sign(kx)))[0]
        assert len(sc) == 1
        j = sc[0]
        t_cross = t[j] + (0 - kx[j]) * (t[j + 1] - t[j]) / (kx[j + 1] - kx[j])
        ky, kz = k_traj_adc[1, i0 + j], k_traj_adc[2, i0 + j]
        out.append((shot, float(t_cross - t_exc[shot]), (float(ky), float(kz))))
    shutil.rmtree(tmp)
    return out


def build(arm, ge=True, plot=True):
    from sequences.EPIcal import generate_epical
    from sequences.meGRE import _patch_scan_info
    from sequences.noise import generate_noise

    a = ARMS[arm]
    mask = load_arm_mask(arm)
    stats = Stats(arm, mask)
    tmp = tempfile.mkdtemp(prefix=f'build_{arm}_', dir=ROOT)
    p, info = design(arm, mask, stats, tmp)
    out = os.path.join(ROOT, a['folder'])
    os.makedirs(out, exist_ok=True)
    print(json.dumps(info, default=float), flush=True)

    arb.generate_arbepi(p.custom_omegas, p)
    generate_epical(p)
    generate_noise(p)
    # scan_info.mat carries the realized meGRE geometry/TE/TR/flip for preprocessing, as generate_degre/
    # generate_megre would patch it: the meGRE whose dwell matches this EPI's.
    from sequences.meGRE import design_megre

    _patch_scan_info(p, design_megre(p, **megre_kwargs(p, info['dwell'])[0]))

    meas = verify_te(p, info, stats)
    info['dc_shots'] = [m[0] for m in meas]
    info['te_measured_dc_s'] = [m[1] for m in meas]
    info['dc_trajectory_ky_kz_at_crossing'] = [list(m[2]) for m in meas]
    for shot, te_meas, (ky, kz) in meas:
        print(f'verify: train {shot}: DC kx=0 crossing {te_meas * 1e3:.4f} ms after RF center '
              f'(ky, kz at that sample: {ky:.3g}, {kz:.3g})', flush=True)

    if plot:
        from plot.plot_last_run import plot_last_run

        plot_last_run(p)
        for f in glob.glob(os.path.join(tmp, '*.png')):
            shutil.move(f, os.path.join(out, os.path.basename(f)))
    shutil.move(os.path.join(tmp, 'scan_info.mat'), os.path.join(out, 'scan_info.mat'))
    for n in ('ArbEPI', 'EPIcal', 'noise'):
        sp = os.path.join(tmp, f'{n}.seq')
        if ge:
            rep = export(p, sp, os.path.join(out, n))
            info[f'peak_pns_percent_{n}'] = rep.peak_pns_percent
            info[f'acoustics_{n}'] = rep.acoustics_max_in_band
        shutil.move(sp, os.path.join(out, f'{n}.seq'))
    json.dump(info, open(os.path.join(out, 'build_info.json'), 'w'), indent=1, default=float)
    shutil.rmtree(tmp)
    return info


if __name__ == '__main__':
    args = sys.argv[1:]
    if args[:1] == ['measure']:
        measure(args[1:])
    elif args[:1] == ['build']:
        for name in args[1:] or list(ARMS):
            build(name)
    elif args[:1] == ['megre']:
        build_megre(4e-6, '4_meGRE_dwell4us')  # every EPI arm has a 4 us dwell
    else:
        raise SystemExit(__doc__)
