"""Check each EPI folder's scan_info.mat against what was requested (run with ArbEPI-python's .venv python, -I):
  * the schedule is identical on every frame (same echo order every frame, same echo times);
  * radial arms (2, 3): 1 shot (one echo train) per frame holding all 60 / 44 samples; k = 0, i.e. (ky, kz) =
    (45, 30), is played once per frame on the same echo with a saved echo time of ~30 ms; the samples are the
    provided mask's with its one nearest-to-center sample moved onto DC;
  * arm 1: ETL 90, 60 shots, trains are ky lines (one kz, ky = 0..89), DC once per frame at ~33.3 ms;
  * scan length and the meGRE keys for preprocessing."""
import glob, os
import hdf5storage
import numpy as np
import scipy.io as sio

ROOT = os.path.dirname(os.path.abspath(__file__))
masks = {'2': 'mask_60.mat', '3': 'mask_44.mat'}
for d in sorted(glob.glob(os.path.join(ROOT, '[123]_*'))):
    k = os.path.basename(d)[0]
    si = hdf5storage.loadmat(os.path.join(d, 'scan_info.mat'))
    sch = si['schedules']  # (frames, shots, ETL, 3), 1-based ky, kz, echo time (s)
    nf, ns, etl, _ = sch.shape
    assert (sch == sch[:1]).all(), 'schedule (incl. echo times) differs between frames'
    iy, iz, te = sch[..., 0].astype(int) - 1, sch[..., 1].astype(int) - 1, sch[..., 2]  # (frame, shot, echo)
    dc = (iy == 45) & (iz == 30)
    assert (dc.sum(axis=(1, 2)) == 1).all(), 'k = 0 must be in every frame exactly once'
    if k == '1':
        assert (etl, ns) == (90, 60)
        assert (iz[0] == iz[0, :, :1]).all(), 'arm 1: kz must be constant within a train'
        assert (iy[0] == np.arange(90)).all(), 'arm 1: every train must be ky = 0..89'
        assert (iz[0, :, 0] == np.arange(60)).all()
    else:
        assert ns == 1
        orig = sio.loadmat(os.path.join(ROOT, masks[k]))['samp'].astype(bool)
        m = np.zeros_like(orig)
        m[iy[0].ravel(), iz[0].ravel()] = True
        assert m.sum() == etl == orig.sum(), 'a frame must hold every sample once'
        gone, new = orig & ~m, m & ~orig
        assert gone.sum() == 1 and new.sum() == 1 and new[45, 30], (gone.sum(), new.sum())
    e = np.unique(np.nonzero(dc)[2])
    assert len(e) == 1
    t_dc = te[dc]
    vt = float(np.squeeze(si['volume_tr']))
    print(f'{os.path.basename(d)}: {nf} frames x {ns} shot(s) x ETL {etl}, identical frames; DC on echo {e[0]} '
          f'(0-based) in {int(dc.any(axis=(1, 2)).sum())} of {nf} frames; saved echo time of DC '
          f'{t_dc.min() * 1e3:.4f}..{t_dc.max() * 1e3:.4f} ms; frame TR {vt * 1e3:.3f} ms; scan {nf * vt:.2f} s')
    if k != '1':
        print(f'   {etl} samples/frame (R {orig.size / etl:.1f}); sample moved onto DC (ky, kz, 0-based):',
              [tuple(int(v) for v in x) for x in np.argwhere(gone)])
    print('   meGRE keys:', {kk: (np.asarray(si[kk]).shape if np.asarray(si[kk]).size > 1 else float(np.squeeze(si[kk])))
                            for kk in ('n_echoes_degre', 'TE_degre', 'TR_degre', 'alpha_degre', 'dwell_degre',
                                       'Nx_degre', 'Ny_degre', 'Nz_degre') if kk in si})
print('OK')
