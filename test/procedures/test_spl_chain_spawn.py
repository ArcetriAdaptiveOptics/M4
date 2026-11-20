"""
Integration test: :meth:`m4.spl.SPLChain.analysis` with a REAL spawn-based
``ProcessPoolExecutor``.

Spawned workers re-import ``opticalib``, which builds its data tree from the
``AOCONF`` environment variable at import time, so monkeypatching in the test
process is invisible to them. The whole scenario therefore runs in a fresh
interpreter whose ``AOCONF`` points at a throw-away configuration inside
``tmp_path``; parent and workers then resolve the same folders.
"""

import json
import os
import subprocess
import sys
import textwrap

import numpy as np

from opticalib.core.config import CONFIG_SCHEMA_VERSION

N = 6
DELTA = [float(d) for d in np.arange(8) * 10.0]

_CONFIG = """SYSTEM:
  version: '{version}'
  data_path: '{data_path}'
PHASING:
  cameras: []
"""

_SCRIPT = textwrap.dedent(
    """
    import json, os
    import numpy as np
    from opticalib.ground import osutils as osu
    from opticalib.core.root import folders as fn
    from opticalib.procedures import phasing
    from m4.spl import SPLChain

    N = 6
    ALV = np.array([450, 500, 550, 600, 650])
    DELTA = np.arange(8) * 10.0


    class Filter:
        wl = None
        mode = None

        def move_to(self, wl):
            self.wl = wl

        def set_bandwidth_mode(self, n):
            self.mode = n


    class Cam:
        def __init__(self, k, flt):
            self.k, self.filter, self.name = k, flt, f"fakecam{k}"
            self.rng = np.random.default_rng(100 + k)

        def set_exptime(self, us):
            pass

        def acquire_frames(self, n):
            img = 10.0 + self.rng.normal(0, 1.0, (300, 300))
            if self.filter.mode != 1:
                yy, xx = np.mgrid[0:300, 0:300]
                period = self.filter.wl / 10.0 * (1 + 0.25 * self.k)
                env = np.exp(-((xx - 150) ** 2 + (yy - 150) ** 2) / (2 * 12.0**2))
                img = img + 1000.0 * env * (1 + 0.8 * np.cos(2 * np.pi * yy / period))
            return img


    def main():
        flt = Filter()
        cams = [Cam(k, flt) for k in range(N)]
        chain = SPLChain(
            cameras=cams, tunable_filter=flt, tnfringes="tnfr", angles=[0.0] * N
        )
        chain.set_filter_mode("narrow")
        chain.acquire_dark_frame(1.0)
        tn = chain.acquire(1.0, lambda_vector=ALV)

        fold = os.path.join(fn.SPL_FRINGES_ROOT_FOLDER, "tnfr")
        os.makedirs(fold)
        lam = np.arange(400, 705, 5, dtype=float)
        osu.save_fits(os.path.join(fold, "Differential_piston.fits"), DELTA)
        osu.save_fits(os.path.join(fold, "Lambda.fits"), lam)
        idx = np.isin(np.round(lam / 5).astype(int) * 5, ALV)
        rng = np.random.default_rng(0)
        mats = {}
        for k in range(N):
            ref = phasing.SPL(camera="none", tunable_filter="none")
            p = os.path.join(fn.SPL_DATA_ROOT_FOLDER, tn, f"cam{k}")
            ref.process_psfs(datapath=p, n_psfs=1, angles=[0.0], initial_half_size=100)
            cube = osu.load_fits(os.path.join(p, "psf0_cube.fits"))
            mats[k], _ = ref._compute_matrix(ALV, cube)
        for i in range(1, len(DELTA)):
            f = rng.normal(0, 1, (80, lam.size))
            k = i - 2
            if 0 <= k < N:
                f = rng.normal(0, 1e-3, (80, lam.size))
                f[:, idx] = np.asarray(mats[k])
            osu.save_fits(os.path.join(fold, "Fringe_%05d.fits" % i), f)

        seq = chain.analysis(tn=tn, max_workers=1, initial_half_size=100)
        par = chain.analysis(tn=tn, max_workers=2, initial_half_size=100)
        print("RESULT:" + json.dumps(
            {"seq": seq, "par": par, "root": fn.SPL_DATA_ROOT_FOLDER, "tn": tn}
        ))


    if __name__ == "__main__":
        main()
    """
)


def test_analysis_real_spawn_pool_matches_sequential(tmp_path):
    data = tmp_path / "data"
    conf = tmp_path / "configuration.yaml"
    conf.write_text(_CONFIG.format(
        data_path=data.as_posix(), version=CONFIG_SCHEMA_VERSION
    ))

    env = dict(os.environ, AOCONF=str(conf))
    proc = subprocess.run(
        [sys.executable, "-c", _SCRIPT],
        env=env,
        cwd=str(tmp_path),
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-4000:]
    assert "legacy" not in proc.stderr.lower(), proc.stderr[-2000:]
    line = [l for l in proc.stdout.splitlines() if l.startswith("RESULT:")]
    assert len(line) == 1, proc.stdout[-2000:]
    out = json.loads(line[0][len("RESULT:"):])

    # the scenario really ran inside the throw-away tree
    assert os.path.realpath(out["root"]).startswith(os.path.realpath(str(data)))
    for k in range(N):
        assert (data / "OPTData" / "SPL" / out["tn"] / f"cam{k}" / "psf0_fringes_result.fits").exists()

    expected = [DELTA[k + 1] for k in range(N)]
    for res in (out["seq"], out["par"]):
        assert len(res) == N
        assert all(isinstance(r, float) and np.isfinite(r) for r in res)
    assert out["par"] == out["seq"]
    assert out["par"] == expected
