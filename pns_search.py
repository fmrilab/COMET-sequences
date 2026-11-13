"""Fastest-TR slews (blip / readout fall / readout rise, T/m/s) with peak PNS <= LIMIT for a radial arm.

    cd <ArbEPI-python>; .venv/bin/python <this folder>/pns_search.py <arm>

For each candidate triple, generate.design() gives the TE-anchored per-shot TR (cheap); candidates are
visited in ascending TR, in parallel batches, and the first passing one (in TR order) is the fastest on the
grid. Every frame has the same schedule, so PNS is evaluated on a short build (a few frames) with every
spoiler at its maximum (the worst case; the real build draws smaller ones at random), human weights. generate.build() still gates the real build at 80 %.
"""
import dataclasses
import itertools
import multiprocessing as mp
import os
import sys
import tempfile
import time
import warnings

import numpy as np

warnings.simplefilter('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.getcwd())
arm = sys.argv[1]
sys.argv = sys.argv[:1]
import generate as g  # noqa: E402
from ge.check import check_seq_feasibility  # noqa: E402

LIMIT = float(os.environ.get('LIMIT', 79.0))
BLIPS = [100, 120, 140, 160, 180, 200]
FALLS = [120, 140, 160, 175, 190]
RISES = [100, 120, 140, 155]

a = g.ARMS[arm]
mask = g.load_arm_mask(arm)
stats = g.Stats(arm, mask)  # also installs the arm's static schedule (every frame identical)
N_CHECK = 6  # frames (x Nshots trains) in the PNS check; the schedule is the same on every frame


def evaluate(triple):
    b, f, r = triple
    sl = dict(blip_slew=float(b), ro_slew_fall=float(f), ro_slew_rise=float(r))
    tmp = tempfile.TemporaryDirectory()  # removed when this candidate is done
    try:
        p, info = g.design(arm, mask, stats, tmp.name, slews=sl)
    except Exception as e:  # infeasible readout geometry etc.
        tmp.cleanup()
        return triple, float('inf'), None, repr(e)[:80]
    ps = dataclasses.replace(p, Nframes=N_CHECK, output_dir=tmp.name,
                             custom_omegas=p.custom_omegas[..., :N_CHECK])
    g.arb.draw_spoil_scales = lambda params, rng: (1.0, 1.0, 1.0)  # every spoiler at its maximum
    seq = g.arb.generate_arbepi(ps.custom_omegas, ps, seqname='pnscheck')
    rep = check_seq_feasibility(seq, ps.spec, pns_wt=tuple(ps.PNSwt))
    tmp.cleanup()
    return triple, info['tr_shot'], rep.peak_pns_percent, ''


if __name__ == '__main__':
    t0 = time.time()
    cands = list(itertools.product(BLIPS, FALLS, RISES))
    # cheap TR for every candidate first (design() only), then PNS in TR order
    with mp.get_context('fork').Pool(32) as pool:
        res = []
        for triple, tr, pns, err in pool.imap_unordered(evaluate, cands, chunksize=1):
            res.append((tr, triple, pns))
        res.sort(key=lambda x: (x[0], x[1]))
    # NOTE: evaluate() computes PNS for all candidates (the grid is small); report the passing frontier
    ok = [r for r in res if r[2] is not None and r[2] <= LIMIT]
    print(f'{arm}: max step ky {stats.steps()[0]:.0f} kz {stats.steps()[1]:.0f}, '
          f'{len(cands)} candidates in {time.time() - t0:.0f} s, LIMIT {LIMIT}')
    for tr, triple, pns in ok[:12]:
        print(f'  TR/shot {tr * 1e3:.3f} ms  blip/fall/rise {triple}  PNS {pns:.1f}%')
    if not ok:
        print('  none under the limit')
    print('  defaults (200/190/155):', [r for r in res if r[1] == (200, 190, 155)])
