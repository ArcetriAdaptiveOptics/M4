"""
Tests for :class:`m4.spl.SPLChain`: one shared tunable filter, several cameras.

Everything runs against fakes; no hardware nor real configuration is touched.
"""

import os
import glob
import time
import threading
import types
from concurrent import futures

import numpy as np
import pytest

from opticalib.procedures import phasing
from opticalib.ground import osutils as osu
from m4 import spl as m4spl
from m4.spl import SPLChain, SPLChainError

N = 6
MAIN = threading.get_ident()


# --------------------------------------------------------------------------
# Fakes
# --------------------------------------------------------------------------


class Recorder:
    def __init__(self):
        self.lock = threading.Lock()
        self.events = []  # (kind, arg, t0, t1, thread, who)

    def add(self, *ev):
        with self.lock:
            self.events.append(ev)

    def of(self, kind, who=None):
        return [
            e for e in self.events if e[0] == kind and (who is None or e[5] == who)
        ]


class FakeFilter:
    def __init__(self, rec, move_time=0.0, fail_move_at=None):
        self.rec = rec
        self.move_time = move_time
        self.current_wl = None
        self.mode = None
        self.calls = []  # ordered ("move_to", wl) / ("set_bandwidth_mode", n)
        self.fail_move_at = fail_move_at

    def move_to(self, wl):
        t0 = time.perf_counter()
        self.calls.append(("move_to", wl))
        self.current_wl = None
        time.sleep(self.move_time)
        if self.fail_move_at is not None and wl == self.fail_move_at:
            raise RuntimeError("move failed")
        self.current_wl = wl
        self.rec.add(
            "move", wl, t0, time.perf_counter(), threading.get_ident(), "filter"
        )

    def set_bandwidth_mode(self, n):
        t0 = time.perf_counter()
        self.calls.append(("set_bandwidth_mode", n))
        self.mode = n
        self.rec.add(
            "mode", n, t0, time.perf_counter(), threading.get_ident(), "filter"
        )


class FakeCamera:
    def __init__(
        self,
        k,
        rec,
        flt,
        exp_time=0.0,
        set_time=0.0,
        fail_at_wl=None,
        psf=False,
    ):
        self.k = k
        self.name = f"fakecam{k}"
        self.rec = rec
        self.filter = flt
        self.exp_time = exp_time
        self.set_time = set_time
        self.fail_at_wl = fail_at_wl
        self.psf = psf
        self.rng = np.random.default_rng(100 + k)
        self.closed = 0
        self.expected = {}  # wl -> returned frame

    def set_exptime(self, us):
        t0 = time.perf_counter()
        time.sleep(self.set_time)
        self.rec.add(
            "setexp", us, t0, time.perf_counter(), threading.get_ident(), self.k
        )

    def acquire_frames(self, nframes):
        t0 = time.perf_counter()
        wl = self.filter.current_wl
        if self.fail_at_wl is not None and wl == self.fail_at_wl:
            raise RuntimeError(f"cam{self.k} boom")
        time.sleep(self.exp_time)
        black = self.filter.mode == 1
        if self.psf:
            img = 10.0 + self.rng.normal(0, 1.0, (300, 300))
            if not black:
                yy, xx = np.mgrid[0:300, 0:300]
                period = (wl or 500) / 10.0 * (1 + 0.25 * self.k)
                env = np.exp(-((xx - 150) ** 2 + (yy - 150) ** 2) / (2 * 12.0**2))
                img = img + 1000.0 * env * (1 + 0.8 * np.cos(2 * np.pi * yy / period))
        elif black:
            img = np.full((8, 8), 10.0 * (self.k + 1))
        else:
            img = np.full((8, 8), 1000.0 * (self.k + 1) + wl)
            self.expected[wl] = img.copy()
        self.rec.add(
            "acq", wl, t0, time.perf_counter(), threading.get_ident(), self.k, black
        )
        return img

    def close(self):
        self.closed += 1


class _NoClose:
    name = "noclose"

    def set_exptime(self, us):
        pass

    def acquire_frames(self, n):
        return np.zeros((8, 8))


# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------


def _patch_env(mp, root):
    data = os.path.join(root, "SPL")
    fr = os.path.join(root, "SPLFringes")
    os.makedirs(data, exist_ok=True)
    os.makedirs(fr, exist_ok=True)
    mp.setattr(phasing._fn, "SPL_DATA_ROOT_FOLDER", data, raising=False)
    mp.setattr(phasing._fn, "SPL_FRINGES_ROOT_FOLDER", fr, raising=False)
    mp.setattr(m4spl, "_splconf", {})
    mp.setattr(phasing, "_splconf", {})
    mp.setattr(osu, "_OPTDATA", str(root))

    def _no_hw():
        raise AssertionError("_get_tunable_filter must not be called")

    mp.setattr(phasing, "_get_tunable_filter", _no_hw)
    return data, fr


@pytest.fixture
def env(tmp_path, monkeypatch):
    data, fr = _patch_env(monkeypatch, str(tmp_path))
    return types.SimpleNamespace(data=data, fringes=fr, tmp=str(tmp_path))


def make_chain(
    move=0.0, exp=0.0, setx=0.0, fail=None, psf=False, tnfringes=None, n=N, flt=None
):
    rec = Recorder()
    flt = flt or FakeFilter(rec, move_time=move)
    rec = flt.rec
    cams = [
        FakeCamera(
            k,
            rec,
            flt,
            exp_time=exp,
            set_time=setx,
            fail_at_wl=(fail[1] if fail and fail[0] == k else None),
            psf=psf,
        )
        for k in range(n)
    ]
    chain = SPLChain(
        cameras=cams, tunable_filter=flt, tnfringes=tnfringes, angles=[0.0] * n
    )
    return chain, cams, flt, rec


LV = np.array([450, 500, 540, 600, 660])


def _hdr(path):
    return osu.load_fits(path).header


def _first_moves(rec):
    """wl -> first move event (the final park at 600 nm is a later duplicate)."""
    out = {}
    for e in rec.of("move"):
        out.setdefault(e[1], e)
    return out


def _gain(wl):
    return float(phasing.SPL._expgain(np.array([wl]))[0])


# --------------------------------------------------------------------------
# Construction
# --------------------------------------------------------------------------


class TestConstruction:
    def test_members_share_filter(self, env):
        chain, cams, flt, _ = make_chain()
        assert len(chain.spls) == N
        assert all(s._filter is flt for s in chain.spls)
        assert [s._camera for s in chain.spls] == cams
        assert all(isinstance(s, phasing.SPL) for s in chain.spls)

    def test_angles_mismatch(self, env):
        rec = Recorder()
        flt = FakeFilter(rec)
        cams = [FakeCamera(k, rec, flt) for k in range(N)]
        with pytest.raises(ValueError):
            SPLChain(cameras=cams, tunable_filter=flt, angles=[0.0] * 5)

    def test_spls_returns_copy(self, env):
        chain, *_ = make_chain()
        chain.spls.clear()
        assert len(chain.spls) == N

    def test_config_cameras_and_angles(self, env, monkeypatch):
        monkeypatch.setattr(
            m4spl,
            "_splconf",
            {"cameras": ["cameras:none"] * 3, "psfs_angles": [1.0, 2.0, 3.0]},
        )
        c = SPLChain(tunable_filter="none")
        assert len(c.spls) == 3
        assert c._angles == [1.0, 2.0, 3.0]
        assert all(s._camera is None for s in c.spls)

    def test_filter_none_string_and_failed_creation(self, env):
        # "none" disables; None tries the (patched, raising) factory and degrades
        c = SPLChain(cameras=["none"], tunable_filter="none")
        with pytest.raises(RuntimeError):
            c.set_filter_mode("narrow")
        c2 = SPLChain(cameras=["none"], tunable_filter=None)
        with pytest.raises(RuntimeError):
            c2.set_filter_mode("narrow")

    def test_set_tn_fringes_propagates(self, env):
        chain, *_ = make_chain()
        chain.set_tn_fringes("abc")
        assert all(s._tnfringes == "abc" for s in chain.spls)
        chain.set_tn_fringes(None)
        assert all(s._tnfringes is None for s in chain.spls)


# --------------------------------------------------------------------------
# Filter mode
# --------------------------------------------------------------------------


class TestFilterMode:
    @pytest.mark.parametrize(
        "mode,n", [("narrow", 8), ("medium", 4), ("wide", 2), ("black", 1)]
    )
    def test_calls_filter_once(self, env, mode, n):
        chain, _, flt, _ = make_chain()
        chain.set_filter_mode(mode)
        assert flt.calls == [("set_bandwidth_mode", n)]

    @pytest.mark.parametrize("bad", ["bogus", "", None, 8, "NARROW"])
    def test_invalid_mode(self, env, bad):
        chain, _, flt, _ = make_chain()
        with pytest.raises((ValueError, TypeError)):
            chain.set_filter_mode(bad)
        assert flt.calls == []

    def test_black_not_remembered(self, env):
        chain, *_ = make_chain()
        chain.set_filter_mode("wide")
        chain.set_filter_mode("black")
        assert chain._filter_mode == "wide"

    def test_acquire_without_mode(self, env):
        chain, cams, flt, _ = make_chain()
        with pytest.raises(ValueError):
            chain.acquire(0.01, lambda_vector=LV)
        assert flt.calls == []
        assert os.listdir(env.data) == []
        assert all(c.expected == {} for c in cams)

    def test_acquire_black_rejected(self, env):
        chain, _, flt, _ = make_chain()
        with pytest.raises(ValueError):
            chain.acquire(0.01, filter_mode="black", lambda_vector=LV)
        assert flt.calls == []
        assert os.listdir(env.data) == []

    def test_acquire_invalid_args_leave_no_folder(self, env):
        chain, _, flt, _ = make_chain()
        chain.set_filter_mode("narrow")
        flt.calls.clear()
        for kw in (
            dict(exptime=0.01, lambda_vector=[300, 500]),
            dict(exptime=0.01, lambda_vector=[500, 800]),
            dict(exptime=[0.01] * 5, lambda_vector=LV),
            dict(exptime=[0.01] * 7, lambda_vector=LV),
            dict(exptime=0.01, lambda_vector=LV, masks=[None] * 5),
        ):
            with pytest.raises(ValueError):
                chain.acquire(**kw)
        assert flt.calls == []
        assert os.listdir(env.data) == []


# --------------------------------------------------------------------------
# Filter mode read from the filter
# --------------------------------------------------------------------------


class ReadableFilter(FakeFilter):
    """Fake filter exposing `get_bandwidth_mode` (settable code, or raising)."""

    def __init__(self, rec, code=None, **kw):
        super().__init__(rec, **kw)
        self.code = code
        self.raises = False

    def get_bandwidth_mode(self):
        if self.raises:
            raise RuntimeError("no answer")
        return self.code


def make_readable(code, n=2):
    flt = ReadableFilter(Recorder(), code=code)
    chain, cams, flt, rec = make_chain(n=n, flt=flt)
    return chain, flt


def _modes_before_first_move(flt):
    calls = flt.calls
    first = next(i for i, c in enumerate(calls) if c[0] == "move_to")
    return [c[1] for c in calls[:first] if c[0] == "set_bandwidth_mode"]


class TestFilterModeRead:
    def test_init_seeds_mode(self, env):
        chain, _ = make_readable(4)
        assert chain._filter_mode == "medium"

    def test_init_black_or_missing_filter(self, env):
        chain, _ = make_readable(1)
        assert chain._filter_mode is None
        c = SPLChain(cameras=["none"], tunable_filter="none")
        assert c._filter_mode is None

    def test_acquire_uses_current_mode(self, env):
        chain, flt = make_readable(4)
        chain.set_filter_mode("narrow")  # remembered: narrow
        flt.code = 2  # but the filter is now wide
        flt.calls.clear()
        chain.acquire(0.01, lambda_vector=LV)
        assert _modes_before_first_move(flt)[-1] == 2

    def test_black_falls_back_to_remembered(self, env):
        chain, flt = make_readable(4)
        flt.code = 1
        chain.acquire(0.01, lambda_vector=LV)
        assert _modes_before_first_move(flt)[-1] == 4

    def test_black_nothing_remembered(self, env):
        chain, flt = make_readable(1)
        with pytest.raises(ValueError):
            chain.acquire(0.01, lambda_vector=LV)
        assert flt.calls == []
        assert os.listdir(env.data) == []

    def test_getter_raising_falls_back(self, env):
        chain, flt = make_readable(8)
        flt.raises = True
        chain.acquire(0.01, lambda_vector=LV)
        assert _modes_before_first_move(flt)[-1] == 8

    def test_dark_restores_mode_filter_had(self, env):
        chain, flt = make_readable(4)
        flt.code = 2  # changed behind the chain's back
        chain.acquire_dark_frame(1.0)
        assert flt.calls == [("set_bandwidth_mode", 1), ("set_bandwidth_mode", 2)]

    def test_wavelength_upper_limit(self, env):
        chain, _, flt, _ = make_chain(n=2)
        chain.set_filter_mode("narrow")
        chain.acquire(0.01, lambda_vector=[720])
        flt.calls.clear()
        with pytest.raises(ValueError):
            chain.acquire(0.01, lambda_vector=[721])
        assert flt.calls == []


# --------------------------------------------------------------------------
# Dark
# --------------------------------------------------------------------------


class TestDark:
    def test_dark_restores_previous_mode(self, env):
        chain, cams, flt, _ = make_chain()
        chain.set_filter_mode("narrow")
        darks = chain.acquire_dark_frame(2.0)
        assert flt.calls == [
            ("set_bandwidth_mode", 8),
            ("set_bandwidth_mode", 1),
            ("set_bandwidth_mode", 8),
        ]
        assert len(darks) == N
        for k, d in enumerate(darks):
            np.testing.assert_allclose(np.asarray(d), 10.0 * (k + 1) / 2.0)
            assert d.header["ISDARK"] is True
            assert d.header["EXPTIME"] == 2.0
        assert all(len(c.rec.of("acq", k)) == 1 for k, c in enumerate(cams))

    def test_dark_without_previous_mode(self, env):
        chain, _, flt, _ = make_chain()
        chain.acquire_dark_frame(1.0)
        assert flt.calls == [("set_bandwidth_mode", 1)]

    def test_dark_returns_copies(self, env):
        chain, *_ = make_chain()
        d = chain.acquire_dark_frame(1.0)
        d[0][...] = -1
        assert np.all(np.asarray(chain.spls[0]._darkFrame1sec) == 10.0)

    def test_dark_camera_failure_still_restores_mode(self, env):
        chain, cams, flt, _ = make_chain()
        chain.set_filter_mode("wide")
        cams[3].acquire_frames = lambda n: (_ for _ in ()).throw(RuntimeError("x"))
        with pytest.raises(SPLChainError) as ei:
            chain.acquire_dark_frame(1.0)
        assert set(ei.value.errors) == {3}
        assert flt.calls[-1] == ("set_bandwidth_mode", 2)

    def test_acquire_takes_dark_before_sweep_mode(self, env):
        chain, cams, flt, rec = make_chain()
        chain.acquire(0.01, filter_mode="medium", lambda_vector=LV)
        calls = flt.calls
        first_move = next(i for i, c in enumerate(calls) if c[0] == "move_to")
        modes = [c[1] for c in calls[:first_move] if c[0] == "set_bandwidth_mode"]
        assert modes == [1, 4]
        # dark acquired (black) by all cameras, every sweep frame not black
        acqs = rec.of("acq")
        assert sum(1 for e in acqs if e[6]) == N
        assert sum(1 for e in acqs if not e[6]) == N * len(LV)
        first_move_t = rec.of("move")[0][2]
        assert all(e[3] <= first_move_t for e in acqs if e[6])

    def test_preacquired_dark_not_reacquired(self, env):
        chain, cams, flt, rec = make_chain()
        chain.set_filter_mode("narrow")
        chain.acquire_dark_frame(1.0)
        n_black = sum(1 for e in rec.of("acq") if e[6])
        flt.calls.clear()
        chain.acquire(0.01, lambda_vector=LV)
        assert sum(1 for e in rec.of("acq") if e[6]) == n_black
        assert [c for c in flt.calls if c[0] == "set_bandwidth_mode"] == [
            ("set_bandwidth_mode", 8)
        ]


# --------------------------------------------------------------------------
# Ordering and timing
# --------------------------------------------------------------------------


class TestOrdering:
    def test_no_move_overlaps_acquisition(self, env):
        chain, cams, flt, rec = make_chain(move=0.05, exp=0.03, setx=0.02)
        chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV)
        for m in rec.of("move"):
            for a in rec.of("acq"):
                assert m[3] <= a[2] or a[3] <= m[2], (m[:4], a[:4])

    def test_filter_only_on_main_thread(self, env):
        chain, cams, flt, rec = make_chain(move=0.01, exp=0.01)
        chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV)
        assert {e[4] for e in rec.of("move")} == {MAIN}
        assert {e[4] for e in rec.of("mode")} == {MAIN}
        # cameras are driven from worker threads
        assert MAIN not in {e[4] for e in rec.of("acq")}

    def test_barrier_per_wavelength(self, env):
        chain, cams, flt, rec = make_chain(move=0.03, exp=0.02)
        chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV)
        moves = _first_moves(rec)
        sweep = [e for e in rec.of("acq") if not e[6]]
        for i, wl in enumerate(LV.tolist()):
            mine = [e for e in sweep if e[1] == wl]
            assert sorted(e[5] for e in mine) == list(range(N))
            assert all(e[2] >= moves[wl][3] for e in mine)
            if i + 1 < len(LV):
                nxt = moves[LV[i + 1]]
                assert all(e[3] <= nxt[2] for e in mine)
        # every camera saw the wavelength it was told to acquire
        for c in cams:
            assert sorted(c.expected) == sorted(LV.tolist())

    def test_cameras_expose_in_parallel(self, env):
        # serial exposure would take N * exp = 3 s: robust to thread start jitter
        exp = 0.5
        chain, cams, flt, rec = make_chain(exp=exp)
        chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV[:1])
        sweep = [e for e in rec.of("acq") if not e[6]]
        assert len(sweep) == N
        # all exposures share a common instant (start jitter << exp)
        assert max(e[2] for e in sweep) < min(e[3] for e in sweep)
        wall = max(e[3] for e in sweep) - min(e[2] for e in sweep)
        assert wall < 0.5 * N * exp

    def test_exptime_update_overlaps_move(self, env):
        chain, cams, flt, rec = make_chain(move=0.1, exp=0.01, setx=0.05)
        base = 0.01
        chain.acquire(base, filter_mode="narrow", lambda_vector=LV)
        moves = _first_moves(rec)
        prev = 1.0  # exposure left by the dark
        checked = 0
        for wl in LV.tolist():
            target = base * _gain(wl)
            if np.isclose(target, prev):
                continue
            prev = target
            m = moves[wl]
            for k in range(N):
                sets = [
                    e
                    for e in rec.of("setexp", k)
                    if np.isclose(e[1], target * 1e6) and e[2] < m[3] + 0.5
                ]
                assert any(s[2] < m[3] and s[3] > m[2] for s in sets), (wl, k)
            checked += 1
        assert checked >= 3

    def test_unchanged_exptime_not_resent(self, env):
        chain, cams, flt, rec = make_chain()
        # LV 450 and 500 share gain 2 -> one exposure update for both
        chain.acquire(0.01, filter_mode="narrow", lambda_vector=[450, 500])
        vals = [e[1] for e in rec.of("setexp", 0)]
        assert vals.count(vals[-1]) == 1
        assert len(vals) == 2  # dark (1 s) then 0.02 s


class TestExposureCache:
    def test_resent_at_start_of_each_acquire(self, env):
        chain, cams, flt, rec = make_chain()
        chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV[:1])
        n1 = [len(rec.of("setexp", k)) for k in range(N)]
        assert n1 == [2] * N  # dark (1 s), then the sweep exposure
        # same sweep again: the cached value equals the target, yet it is
        # re-sent once per camera (the camera may have been power cycled)
        chain.acquire(0.01, lambda_vector=LV[:1])
        for k in range(N):
            vals = [e[1] for e in rec.of("setexp", k)]
            assert len(vals) == 3
            assert vals[-1] == vals[-2]

    def test_dark_resends_exposure(self, env):
        chain, cams, flt, rec = make_chain()
        chain.acquire_dark_frame(1.0)
        chain.acquire_dark_frame(1.0)
        assert [len(rec.of("setexp", k)) for k in range(N)] == [2] * N

    def test_resent_after_camera_failure(self, env):
        chain, cams, flt, rec = make_chain(fail=(2, int(LV[1])))
        with pytest.raises(SPLChainError):
            chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV)
        assert chain.spls[2]._curr_exptime is None
        assert all(chain.spls[k]._curr_exptime is not None for k in range(N) if k != 2)

    def test_resent_after_exptime_update_failure(self, env):
        chain, cams, flt, rec = make_chain()

        def boom(us):
            if us != 1e6:
                raise RuntimeError("expo")

        cams[5].set_exptime = boom
        with pytest.raises(SPLChainError):
            chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV)
        assert chain.spls[5]._curr_exptime is None

    def test_dark_failure_invalidates(self, env):
        chain, cams, flt, rec = make_chain()
        cams[3].acquire_frames = lambda n: (_ for _ in ()).throw(RuntimeError("x"))
        with pytest.raises(SPLChainError):
            chain.acquire_dark_frame(1.0)
        assert chain.spls[3]._curr_exptime is None
        assert chain.spls[0]._curr_exptime == 1.0


class TestTiming:
    def test_parallel_speedup(self, env):
        move, exp, n = 0.1, 0.2, 4
        lv = [450, 500, 600, 660]
        chain, *_ = make_chain(move=move, exp=exp)
        chain.set_filter_mode("narrow")
        chain.acquire_dark_frame(0.01)
        t0 = time.perf_counter()
        chain.acquire(0.01, lambda_vector=lv)
        dt = time.perf_counter() - t0
        ideal = n * (move + exp)
        assert dt < 0.6 * n * (move + 6 * exp)
        assert dt < 2.0 * ideal + move  # park move


# --------------------------------------------------------------------------
# Files and headers
# --------------------------------------------------------------------------


class TestFiles:
    def test_layout_and_data(self, env):
        chain, cams, flt, rec = make_chain()
        tn = chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV)
        assert os.listdir(env.data) == [tn]
        root = os.path.join(env.data, tn)
        assert "lambda_vector.fits" in os.listdir(root)
        np.testing.assert_array_equal(
            np.asarray(osu.load_fits(os.path.join(root, "lambda_vector.fits"))), LV
        )
        for k in range(N):
            d = os.path.join(root, f"cam{k}")
            files = set(os.listdir(d))
            assert files == {"lambda_vector.fits", "darkFrame.fits"} | {
                f"rawframe_{wl}nm.fits" for wl in LV.tolist()
            }
            for wl in LV.tolist():
                data = np.asarray(osu.load_fits(os.path.join(d, f"rawframe_{wl}nm.fits")))
                np.testing.assert_array_equal(data, cams[k].expected[wl])
            dark = np.asarray(osu.load_fits(os.path.join(d, "darkFrame.fits")))
            np.testing.assert_allclose(dark, 10.0 * (k + 1) / 1.0)
        assert chain._last_measure_tn == tn
        assert all(s._last_measure_tn == tn for s in chain.spls)

    def test_headers(self, env):
        chain, cams, *_ = make_chain()
        base = 0.01
        tn = chain.acquire(base, filter_mode="narrow", lambda_vector=LV)
        gains = phasing.SPL._expgain(LV)
        assert len(set(gains.tolist())) > 1
        for k in range(N):
            for wl, g in zip(LV.tolist(), gains):
                h = _hdr(os.path.join(env.data, tn, f"cam{k}", f"rawframe_{wl}nm.fits"))
                assert h["EXPTIME"] == pytest.approx(base * g)
                assert h["WAVELEN"] == wl
                assert h["EXPGAIN"] == pytest.approx(g)
                assert h["CAMNAME"] == f"fakecam{k}"

    def test_per_camera_exptime(self, env):
        chain, cams, flt, rec = make_chain()
        bases = [0.01 * (k + 1) for k in range(N)]
        tn = chain.acquire(bases, filter_mode="narrow", lambda_vector=LV)
        for k in range(N):
            for wl in LV.tolist():
                h = _hdr(os.path.join(env.data, tn, f"cam{k}", f"rawframe_{wl}nm.fits"))
                assert h["EXPTIME"] == pytest.approx(bases[k] * _gain(wl))
        # the camera was commanded in microseconds
        sent = [e[1] for e in rec.of("setexp", 3)]
        assert pytest.approx(bases[3] * _gain(LV[-1]) * 1e6) == sent[-1]

    def test_numpy_scalar_and_tuple_exptime(self, env):
        chain, *_ = make_chain()
        chain.acquire(np.float64(0.01), filter_mode="narrow", lambda_vector=LV[:2])
        chain.acquire(tuple([0.01] * N), lambda_vector=LV[:2])

    def test_masks_used(self, env):
        chain, *_ = make_chain()
        masks = [np.zeros((8, 8)) for _ in range(N)]
        masks[4][0, 0] = 1
        tn = chain.acquire(
            0.01, filter_mode="narrow", lambda_vector=LV[:1], masks=masks
        )
        wl = int(LV[0])
        for k in range(N):
            m = np.asarray(
                osu.load_fits(
                    os.path.join(env.data, tn, f"cam{k}", f"rawframe_{wl}nm.fits")
                ).mask
            )
            assert bool(m[0, 0]) == (k == 4)

    def test_masks_numpy_array(self, env):
        chain, *_ = make_chain()
        masks = np.zeros((N, 8, 8))
        masks[4, 0, 0] = 1
        tn = chain.acquire(
            0.01, filter_mode="narrow", lambda_vector=LV[:1], masks=masks
        )
        wl = int(LV[0])
        for k in range(N):
            m = np.asarray(
                osu.load_fits(
                    os.path.join(env.data, tn, f"cam{k}", f"rawframe_{wl}nm.fits")
                ).mask
            )
            assert bool(m[0, 0]) == (k == 4)
        with pytest.raises(ValueError):
            chain.acquire(0.01, lambda_vector=LV[:1], masks=np.zeros((N - 1, 8, 8)))

    def test_single_wavelength_and_boundaries(self, env):
        chain, flt = make_chain()[0], None
        chain.acquire(0.01, filter_mode="narrow", lambda_vector=[400])
        tn = chain.acquire(0.01, lambda_vector=[700])
        assert os.path.exists(
            os.path.join(env.data, tn, "cam5", "rawframe_700nm.fits")
        )

    def test_default_lambda_vector(self, env):
        chain, cams, flt, _ = make_chain()
        tn = chain.acquire(0.01, filter_mode="narrow")
        lv = np.arange(440, 721, 20)
        moves = [c[1] for c in flt.calls if c[0] == "move_to"]
        assert moves == lv.tolist() + [600]

    def test_unique_tn_per_acquire(self, env):
        chain, *_ = make_chain()
        a = chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV[:1])
        b = chain.acquire(0.01, lambda_vector=LV[:1])
        assert a != b


# --------------------------------------------------------------------------
# Failure policy
# --------------------------------------------------------------------------


def _nondaemon_threads():
    return {t for t in threading.enumerate() if not t.daemon}


class TestFailures:
    def test_camera_failure_midsweep(self, env):
        before = _nondaemon_threads()
        chain, cams, flt, rec = make_chain(exp=0.01, fail=(2, int(LV[2])))
        with pytest.raises(SPLChainError) as ei:
            chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV)
        err = ei.value
        assert set(err.errors) == {2}
        assert "cam2" in str(err) and "boom" in str(err)
        assert isinstance(err, RuntimeError)
        assert flt.calls[-1] == ("move_to", 600)
        assert [c[1] for c in flt.calls if c[0] == "move_to"] == LV[:3].tolist() + [600]

        (tn,) = os.listdir(env.data)
        for k in range(N):
            files = glob.glob(os.path.join(env.data, tn, f"cam{k}", "rawframe_*"))
            names = {os.path.basename(f) for f in files}
            exp = {f"rawframe_{wl}nm.fits" for wl in LV[:2].tolist()}
            if k != 2:
                exp.add(f"rawframe_{int(LV[2])}nm.fits")
            assert names == exp, k
        # no lingering worker threads
        time.sleep(0.05)
        assert _nondaemon_threads() <= before
        assert chain._last_measure_tn is None

    def test_multiple_camera_failures_reported(self, env):
        chain, cams, flt, _ = make_chain()
        for k in (1, 4):
            cams[k].fail_at_wl = int(LV[0])
        with pytest.raises(SPLChainError) as ei:
            chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV)
        assert set(ei.value.errors) == {1, 4}

    def test_exptime_update_failure(self, env):
        chain, cams, flt, _ = make_chain()

        def boom(us):
            if us != 1e6:  # let the dark (1 s) through
                raise RuntimeError("expo")

        cams[5].set_exptime = boom
        with pytest.raises(SPLChainError) as ei:
            chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV)
        assert set(ei.value.errors) == {5}
        assert flt.calls[-1] == ("move_to", 600)

    def test_filter_move_failure_propagates_and_parks(self, env):
        rec = Recorder()
        flt = FakeFilter(rec, fail_move_at=int(LV[1]))
        chain, cams, flt, rec = make_chain(flt=flt)
        with pytest.raises(RuntimeError, match="move failed"):
            chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV)
        assert flt.calls[-1] == ("move_to", 600)
        assert [c[1] for c in flt.calls if c[0] == "move_to"][-2:] == [int(LV[1]), 600]

    def test_writer_failure(self, env):
        before = _nondaemon_threads()
        chain, cams, flt, rec = make_chain()

        def bad(*a, **k):
            raise OSError("disk full")

        chain.spls[1]._save_frame = bad
        with pytest.raises(SPLChainError) as ei:
            chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV)
        assert set(ei.value.errors) == {1}
        assert "disk full" in str(ei.value)
        assert flt.calls[-1] == ("move_to", 600)
        # the sweep stops early (writes already finished are checked per
        # wavelength) and the filter is parked
        moves = [c[1] for c in flt.calls if c[0] == "move_to"]
        assert moves[-1] == 600 and len(moves) - 1 < len(LV)
        time.sleep(0.05)
        assert _nondaemon_threads() <= before

    def test_writer_failure_after_camera_failure_keeps_camera_error(self, env):
        # both fail at the first wavelength: neither error may be dropped
        chain, cams, flt, _ = make_chain(fail=(2, int(LV[0])))
        chain.spls[0]._save_frame = lambda *a, **k: (_ for _ in ()).throw(OSError("x"))
        with pytest.raises(SPLChainError) as ei:
            chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV)
        assert isinstance(ei.value.errors[2], RuntimeError)
        assert isinstance(ei.value.errors[0], OSError)
        assert flt.calls[-1] == ("move_to", 600)


    def test_error_tn_set_after_folder_creation(self, env):
        chain, cams, flt, _ = make_chain(fail=(2, int(LV[2])))
        with pytest.raises(SPLChainError) as ei:
            chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV)
        (tn,) = os.listdir(env.data)
        assert ei.value.tn == tn

    def test_error_tn_none_before_folder(self, env):
        chain, cams, flt, _ = make_chain()
        cams[3].acquire_frames = lambda n: (_ for _ in ()).throw(RuntimeError("x"))
        with pytest.raises(SPLChainError) as ei:  # failing dark inside acquire
            chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV)
        assert set(ei.value.errors) == {3}
        assert ei.value.tn is None
        assert os.listdir(env.data) == []
        assert SPLChainError("m", {}).tn is None

    def test_writer_failure_sets_tn(self, env):
        chain, cams, flt, rec = make_chain()

        def bad(*a, **k):
            raise OSError("disk full")

        chain.spls[1]._save_frame = bad
        with pytest.raises(SPLChainError) as ei:
            chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV)
        assert ei.value.tn == os.listdir(env.data)[0]
        assert chain._last_measure_tn is None

    def test_writer_error_aborts_sweep_early(self, env):
        chain, cams, flt, rec = make_chain(exp=0.05)

        def bad(*a, **k):
            raise OSError("disk full")

        chain.spls[1]._save_frame = bad
        with pytest.raises(SPLChainError) as ei:
            chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV)
        assert set(ei.value.errors) == {1}
        assert "disk full" in str(ei.value)
        moves = [c[1] for c in flt.calls if c[0] == "move_to"]
        # detected by the end of the second wavelength at the latest, not at
        # the end of the 5-wavelength sweep; the filter is parked
        assert moves[-1] == 600
        assert len(moves) - 1 <= 2 < len(LV)
        assert chain._last_measure_tn is None

    def test_writer_pool_has_one_worker_per_camera(self, env, monkeypatch):
        sizes = []
        real = futures.ThreadPoolExecutor

        class Spy(real):
            def __init__(self, max_workers=None, *a, **k):
                sizes.append(max_workers)
                super().__init__(max_workers, *a, **k)

        monkeypatch.setattr(
            m4spl, "_fut", types.SimpleNamespace(ThreadPoolExecutor=Spy)
        )
        chain, *_ = make_chain(n=4)
        chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV[:1])
        assert sizes and set(sizes) == {4}

    def test_park_failure_after_success(self, env):
        rec = Recorder()
        flt = FakeFilter(rec, fail_move_at=600)
        chain, cams, flt, rec = make_chain(flt=flt)
        with pytest.raises(RuntimeError) as ei:
            chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV[:3])
        assert not isinstance(ei.value, SPLChainError)
        (tn,) = os.listdir(env.data)
        assert tn in str(ei.value) and "parked" in str(ei.value)
        assert isinstance(ei.value.__cause__, RuntimeError)
        assert "move failed" in str(ei.value.__cause__)
        assert chain._last_measure_tn == tn
        assert all(s._last_measure_tn == tn for s in chain.spls)
        for k in range(N):
            assert os.path.exists(
                os.path.join(env.data, tn, f"cam{k}", f"rawframe_{int(LV[2])}nm.fits")
            )

    def test_park_failure_keeps_camera_error(self, env):
        rec = Recorder()
        flt = FakeFilter(rec, fail_move_at=600)
        chain, cams, flt, rec = make_chain(flt=flt, fail=(2, int(LV[1])))
        with pytest.raises(SPLChainError) as ei:
            chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV[:3])
        assert set(ei.value.errors) == {2}
        assert any("parked" in n for n in ei.value.__notes__)
        assert chain._last_measure_tn is None

    def test_park_failure_keeps_writer_error(self, env):
        rec = Recorder()
        flt = FakeFilter(rec, fail_move_at=600)
        chain, cams, flt, rec = make_chain(flt=flt)
        chain.spls[1]._save_frame = lambda *a, **k: (_ for _ in ()).throw(OSError("full"))
        with pytest.raises(SPLChainError) as ei:
            chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV[:3])
        assert set(ei.value.errors) == {1}
        assert any("parked" in n for n in ei.value.__notes__)


class TestConstructorCleanup:
    @staticmethod
    def _fake_gv(monkeypatch, fail_at):
        import opticalib.devices.cameras as cams_mod

        made = []

        class GV:
            def __init__(self, name):
                if len(made) == fail_at:
                    raise OSError("cannot open " + name)
                self.name = name
                self.closed = 0
                made.append(self)

            def close(self):
                self.closed += 1

        monkeypatch.setattr(cams_mod, "GigaVision", GV)
        return made

    def test_closes_opened_cameras_when_later_one_fails(self, env, monkeypatch):
        made = self._fake_gv(monkeypatch, fail_at=2)
        with pytest.raises(OSError, match="cannot open C"):
            SPLChain(cameras=["cameras:A", "cameras:B", "cameras:C"], tunable_filter="none")
        assert [c.closed for c in made] == [1, 1]

    def test_closes_opened_cameras_when_spl_fails(self, env, monkeypatch):
        made = self._fake_gv(monkeypatch, fail_at=99)
        real = m4spl.SPL
        calls = []

        def flaky(*a, **k):
            calls.append(1)
            if len(calls) == 2:
                raise RuntimeError("spl boom")
            return real(*a, **k)

        monkeypatch.setattr(m4spl, "SPL", flaky)
        with pytest.raises(RuntimeError, match="spl boom"):
            SPLChain(cameras=["cameras:A", "cameras:B", "cameras:C"], tunable_filter="none")
        assert [c.closed for c in made] == [1, 1, 1]

    def test_closes_opened_cameras_on_angles_mismatch(self, env, monkeypatch):
        made = self._fake_gv(monkeypatch, fail_at=99)
        with pytest.raises(ValueError):
            SPLChain(cameras=["cameras:A", "cameras:B"], tunable_filter="none", angles=[0.0])
        assert [c.closed for c in made] == [1, 1]

    def test_callers_cameras_are_not_closed(self, env):
        rec = Recorder()
        flt = FakeFilter(rec)
        cams = [FakeCamera(k, rec, flt) for k in range(2)]
        with pytest.raises(ValueError):
            SPLChain(cameras=cams, tunable_filter=flt, angles=[0.0])
        assert [c.closed for c in cams] == [0, 0]


# --------------------------------------------------------------------------
# Analysis
# --------------------------------------------------------------------------

ALV = np.array([450, 500, 550, 600, 650])
DELTA = np.arange(8) * 10.0


@pytest.fixture(scope="module")
def analysed(tmp_path_factory):
    root = str(tmp_path_factory.mktemp("splchain"))
    mp = pytest.MonkeyPatch()
    try:
        data, fr = _patch_env(mp, root)
        chain, cams, flt, rec = make_chain(psf=True, tnfringes="tnfr")
        chain.set_filter_mode("narrow")
        chain.acquire_dark_frame(1.0)
        tn = chain.acquire(1.0, lambda_vector=ALV)

        fold = os.path.join(fr, "tnfr")
        os.makedirs(fold)
        lam = np.arange(400, 705, 5, dtype=float)
        osu.save_fits(os.path.join(fold, "Differential_piston.fits"), DELTA)
        osu.save_fits(os.path.join(fold, "Lambda.fits"), lam)
        idx = np.isin(np.round(lam / 5).astype(int) * 5, ALV)
        rng = np.random.default_rng(0)
        mats = {}
        for k in range(N):
            ref = phasing.SPL(camera="none", tunable_filter="none")
            p = os.path.join(data, tn, f"cam{k}")
            ref.process_psfs(datapath=p, n_psfs=1, angles=[0.0], initial_half_size=100)
            cube = osu.load_fits(os.path.join(p, "psf0_cube.fits"))
            mats[k], _ = ref._compute_matrix(ALV, cube)
        for i in range(1, len(DELTA)):
            f = rng.normal(0, 1, (80, lam.size))
            k = i - 2  # Fringe_i matches camera i-2 (cameras 0..5 -> 2..7)
            if 0 <= k < N:
                f = rng.normal(0, 1e-3, (80, lam.size))
                f[:, idx] = np.asarray(mats[k])
            osu.save_fits(os.path.join(fold, "Fringe_%05d.fits" % i), f)
        yield types.SimpleNamespace(chain=chain, tn=tn, data=data)
    finally:
        mp.undo()


class TestAnalysis:
    expected = [DELTA[k + 1] for k in range(N)]  # Fringe_{k+2} -> delta[k+1]

    def test_sequential_in_camera_order(self, analysed):
        res = analysed.chain.analysis(
            tn=analysed.tn, max_workers=1, initial_half_size=100
        )
        assert len(res) == N
        assert all(isinstance(r, float) for r in res)
        assert res == self.expected
        for k in range(N):
            assert os.path.exists(
                os.path.join(analysed.data, analysed.tn, f"cam{k}", "psf0_fringes_result.fits")
            )

    def test_default_tn_is_last_acquired(self, analysed):
        res = analysed.chain.analysis(max_workers=1, initial_half_size=100)
        assert res == self.expected

    def test_parallel_path_uses_spawn_pool(self, analysed, monkeypatch):
        seen = {}

        class FakePool:
            def __init__(self, workers, mp_context=None, **kw):
                seen["workers"] = workers
                seen["ctx"] = mp_context

            def __enter__(self):
                return self

            def __exit__(self, *a):
                return False

            def submit(self, fn, *args, **kw):
                seen.setdefault("fns", []).append(fn)
                f = futures.Future()
                f.set_result(fn(*args, **kw))
                return f

        fake_fut = types.SimpleNamespace(
            ThreadPoolExecutor=futures.ThreadPoolExecutor, ProcessPoolExecutor=FakePool
        )
        monkeypatch.setattr(m4spl, "_fut", fake_fut)
        res = analysed.chain.analysis(
            tn=analysed.tn, max_workers=2, initial_half_size=100
        )
        assert seen["workers"] == 2
        assert seen["ctx"].get_start_method() == "spawn"
        assert set(seen["fns"]) == {phasing._analyze_datapath}  # picklable, module-level
        assert res == self.expected

    def test_default_workers_is_one_per_camera(self, analysed, monkeypatch):
        seen = {}

        class FakePool:
            def __init__(self, workers, mp_context=None):
                seen["w"] = workers

            def __enter__(self):
                return self

            def __exit__(self, *a):
                return False

            def submit(self, fn, *a, **k):
                f = futures.Future()
                f.set_result(fn(*a, **k))
                return f

        monkeypatch.setattr(
            m4spl, "_fut", types.SimpleNamespace(ProcessPoolExecutor=FakePool)
        )
        analysed.chain.analysis(tn=analysed.tn, initial_half_size=100)
        assert seen["w"] == N

    @staticmethod
    def _flaky(monkeypatch, bad_cams):
        orig = phasing._analyze_datapath

        def flaky(path, *a, **k):
            if any(path.endswith(f"cam{b}") for b in bad_cams):
                raise ValueError(f"bad {os.path.basename(path)}")
            return orig(path, *a, **k)

        monkeypatch.setattr(phasing, "_analyze_datapath", flaky)

    def test_failure_in_one_camera_in_process(self, analysed, monkeypatch):
        self._flaky(monkeypatch, [3])
        with pytest.raises(SPLChainError) as ei:
            analysed.chain.analysis(
                tn=analysed.tn, max_workers=1, initial_half_size=100
            )
        assert set(ei.value.errors) == {3}
        assert isinstance(ei.value.errors[3], ValueError)
        assert "cam3" in str(ei.value)

    def test_failures_all_reported_in_process(self, analysed, monkeypatch):
        self._flaky(monkeypatch, [0, 5])
        with pytest.raises(SPLChainError) as ei:
            analysed.chain.analysis(
                tn=analysed.tn, max_workers=1, initial_half_size=100
            )
        assert set(ei.value.errors) == {0, 5}

    def test_failure_in_one_camera_pool_path(self, analysed, monkeypatch):
        self._flaky(monkeypatch, [2])

        class FakePool:
            def __init__(self, workers, mp_context=None, **kw):
                pass

            def __enter__(self):
                return self

            def __exit__(self, *a):
                return False

            def submit(self, fn, *args, **kw):
                f = futures.Future()
                try:
                    f.set_result(fn(*args, **kw))
                except Exception as e:
                    f.set_exception(e)
                return f

        monkeypatch.setattr(
            m4spl, "_fut", types.SimpleNamespace(ProcessPoolExecutor=FakePool)
        )
        with pytest.raises(SPLChainError) as ei:
            analysed.chain.analysis(
                tn=analysed.tn, max_workers=2, initial_half_size=100
            )
        assert set(ei.value.errors) == {2}

    def test_tnfringes_none(self, env):
        chain, *_ = make_chain()
        chain.acquire(0.01, filter_mode="narrow", lambda_vector=LV)
        with pytest.raises(ValueError):
            chain.analysis()

    def test_tn_none_without_acquire(self, env):
        chain, *_ = make_chain(tnfringes="x")
        with pytest.raises(ValueError):
            chain.analysis()
        with pytest.raises(ValueError):
            chain.plot_comparison()


# --------------------------------------------------------------------------
# close / context manager
# --------------------------------------------------------------------------


class TestClose:
    def test_close_calls_all(self, env):
        chain, cams, *_ = make_chain()
        chain.close()
        assert [c.closed for c in cams] == [1] * N

    def test_close_swallows_errors_and_continues(self, env):
        chain, cams, *_ = make_chain()

        def boom():
            raise RuntimeError("nope")

        cams[1].close = boom
        chain.close()
        assert cams[2].closed == 1 and cams[5].closed == 1

    def test_close_without_close_attr_or_none_camera(self, env):
        rec = Recorder()
        flt = FakeFilter(rec)
        c0 = FakeCamera(0, rec, flt)
        chain = SPLChain(
            cameras=[c0, _NoClose(), "none"], tunable_filter=flt, angles=[0, 0, 0]
        )
        chain.close()
        assert c0.closed == 1

    def test_context_manager(self, env):
        chain, cams, *_ = make_chain()
        with chain as c:
            assert c is chain
            assert all(x.closed == 0 for x in cams)
        assert [x.closed for x in cams] == [1] * N

    def test_context_manager_closes_on_exception(self, env):
        chain, cams, *_ = make_chain()
        with pytest.raises(KeyError):
            with chain:
                raise KeyError("x")
        assert [x.closed for x in cams] == [1] * N
