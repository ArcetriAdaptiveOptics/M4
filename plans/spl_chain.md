# SPLChain: 1 tunable filter, 6 cameras (1 PSF each)

Status legend: `[ ]` todo · `[~]` in progress · `[x]` Done

## Context
`opticalib.procedures.phasing.SPL` assumes 1 camera + 1 tunable filter + `n` PSFs on the same frame.
M4 has 1 tunable filter shared by 6 cameras, each seeing 1 PSF. The wavelength is a shared, serial
resource, so 6 independent `SPL.acquire()` cannot run in parallel; everything between two filter
moves, however, is independent per camera.

**Architecture: composition.** `SPLChain` (M4) owns one shared filter + 6 single-camera `SPL`
instances. Per-instance state (`_camera`, `_curr_exptime`, `_darkFrame1sec`, `_Qt/_idp/_matrix`)
is what makes per-camera parallel acquisition/analysis safe. opticalib `SPL` stays single-camera.

Per-wavelength cost: sequential `T_move + 6·(T_set + T_exp·nframes + T_save)` →
target `T_move + max_cam(T_exp·nframes)`.

## Branches
- opticalib: `spl-generalization` (from `main`)
- M4: `calpygui-procedures` (current)
- No commits unless explicitly requested.

## Findings that shape the plan
- Nested tn (`"<tn>/cam0"`) is NOT supported by `find_tracknum`/`is_tn` → `datapath` override required.
- `SPL.analysis(tn=None)` is broken (joins `None` before resolving the last tn) → fixed in Step 1.
- No tests existed for `phasing.py` → baseline first (Step 0).
- Environment: run everything in conda env `m4` (`conda run -n m4 ...`, Python 3.13): opticalib and m4
  are installed editable from the two repos, all deps (incl. `zmq`) present.
- BUG: `SPL.acquire()` sets `filter_mode`, then (if no dark) `acquire_dark_frame` switches the filter to
  "black" and never restores it → the whole sweep is acquired in black mode. Fixed in Step 1.
- `plot_comparison` needs `self._Qt`, which only exists in the process that ran `analysis` → templates
  must be loadable on demand (needed once analysis runs in worker processes).
- Quirks found by Step 0 — user decisions:
  - default lambda vector 440–720 nm is intended → FIX: validation range becomes 400–720 nm;
  - `Fringe_k` → `delta[k-1]` in template comparison (possible off-by-one) → LEAVE as is;
  - `_heal_bad_pixels` with zero hot pixels re-medians the whole image (`[-0:]`) → FIX: no hot pixels = no healing;
  - `matrix == np.nan` is a no-op → FIX: `matrix[np.isnan(matrix)] = 0` (original intent).
- Filter mode (user): the filter always has a mode once powered on (set automatically), or set via
  `set_filter_mode`. `acquire(filter_mode=None)` uses the known mode; raises only if none is known.
  → plico client `get_bandwidth_mode() -> int` (1 black, 2 wide, 4 medium, 8 narrow). New
  `phasing._read_filter_mode(filt) -> str | None`: seeds `_filter_mode` at init; `acquire(None)` uses the filter's
  CURRENT mode, falling back to the remembered one; dark frames restore the mode the filter had just before.
- DECISION (user): `EXPGAIN` header must store the actual gain from `_expgain` (today it is always 1.0);
  `EXPTIME` stays the actual exposure used. Applied as a follow-up to Steps 1–2 (SPL.acquire, SPLChain, tests).

## File ownership
| Repo | Files | Step |
|---|---|---|
| opticalib | `tests/test_procedures_phasing.py` (new) | 0 |
| opticalib | `opticalib/procedures/phasing.py` | 1 |
| M4 | `m4/spl.py` | 2 |
| M4 | `test/procedures/test_spl_chain.py` (new) | 3 |

## Steps

### [x] Step 0: characterization tests for current SPL (opticalib) — Done
Result: `tests/test_procedures_phasing.py`, 57 tests: 55 passed + 2 strict xfail (black-mode sweep, `analysis(tn=None)`) on unmodified code.
Fake camera (synthetic PSF frames) + recording fake filter, data/fringes folders redirected to tmp.
Pins: files written, FITS headers (`EXPTIME`, `WAVELEN`, `EXPGAIN`), filter call order ending in
`move_to(600)`, `process_psfs`+`analysis` recover a known piston from synthetic templates.
Gate: green on the **unmodified** code.

### [x] Step 1: opticalib refactor (`phasing.py`) — Done
```python
@staticmethod
def _validate_lambda(lambda_vector) -> np.ndarray      # None -> arange(440, 721, 20); range check
@staticmethod
def _expgain(lambda_vector) -> np.ndarray
def _save_frame(self, datapath, img, wl, exptime, expgain, mask=None, header=None) -> None  # no instance state read (thread-safe); header merged (e.g. CAMNAME)
def _resolve_datapath(self, tn=None, datapath=None) -> str
def _load_templates(self, lambda_vector) -> tuple[np.ndarray, np.ndarray]     # (delta, Qt); sets self._Qt
# module level, picklable, used by SPLChain's process pool:
def _analyze_datapath(datapath, tnfringes, **analysis_kwargs) -> list[float]
def process_psfs(self, tn=None, ..., datapath: str | None = None)
def analysis(self, tn=None, n_psfs=None, datapath: str | None = None, **process_kwargs) -> list[float]
def plot_comparison(self, tn=None, psf_n="all", datapath: str | None = None)
# __init__: tunable_filter="none" -> None, no connection attempt
```
`acquire()` rebuilt on the helpers, identical output. Private helpers take the resolved folder and
list files via `get_file_list(fold=path, key=...)`; single-file `str` results normalized to list.
Dark-frame fix: acquire the dark BEFORE setting the sweep filter mode; `set_filter_mode` tracks the last
non-black mode, restored after the dark; if no mode is known and `filter_mode is None` → `ValueError`.
`plot_comparison` calls `_load_templates` when `_Qt` is missing.
Unchanged on purpose: `angles = angles or conf` fallback; camera/reconnect code.
Follow-up: `EXPGAIN` = actual gain (user decision).
Gate: Step 0 tests still green.
Progress: source refactor done (renames now basename-only, so nested `cam{k}` paths are safe). `acquire()`
without a known filter mode now raises `ValueError` (the old code silently swept in black mode): the
12 baseline tests relying on that get the mode set in the fixture.
Result: `EXPGAIN` = actual gain; phasing tests 72 passed (both former xfails fixed; new tests for
`"none"` filter, `_save_frame`, `_resolve_datapath`, nested `cam0` datapath, `plot_comparison` without `_Qt`,
`_analyze_datapath`). Full opticalib suite: 1 failure, likely unrelated
(`test_gui_procedures.py::TestDMCalibration::test_closed_loop_is_an_explicit_loop`) → checked in Step 4.
Note: conda env `m4` has no pytest; ran with pytest borrowed from the `opticalib` env via PYTHONPATH.

### [x] Step 2: M4 `SPLChain` (`m4/spl.py`) — Done
Result: implemented; acquisition smoke test with fakes (6 cams, 0.2 s exp, 0.1 s move): ~0.3 s/wl vs ~1.3 s/wl
sequential, no move/acquire overlap, per-camera failure → `SPLChainError`, others saved, filter back to 600.
Pending: run against the real Step 1 API (Step 4). Open risk: no timeout if a device call hangs
(`GigaVision.acquire_frames` has its own SDK timeout; `move_to` has none).
```python
class SPLChain:
    def __init__(self, cameras=None, tunable_filter=None, tnfringes=None, angles=None)
    spls: list[SPL]                                     # read-only
    def set_tn_fringes(self, tnfringes)
    def set_filter_mode(self, mode)                     # once, shared filter
    def acquire_dark_frame(self, exptime, nframes=1) -> list[ImageData]
    def preview_detection(self, exptime=None, filter_mode=None, wavelength=None, **kw)
    def acquire(self, exptime, filter_mode=None, lambda_vector=None, nframes=1, masks=None) -> str
    def analysis(self, tn=None, max_workers=None, **process_kwargs) -> list[float]
    def plot_comparison(self, tn=None, cam="all")
    def close(self) / __enter__ / __exit__
```
Config: `PHASING.cameras: ["cameras:Name0", ...]`, `PHASING.psfs_angles` (1 per camera).
Loop: `ThreadPoolExecutor(6)`; `set_exptime` ×6 concurrent with `filter.move_to(wl_k)` → barrier →
6× `acquire_frames` → writer pool (2) saves while the loop advances. `finally`: drain writers,
shutdown pools, `move_to(600)`. Per-camera errors aggregated after the step completes.
Layout: `SPL/<tn>/cam{k}/{lambda_vector,darkFrame,rawframe_*}.fits`, `CAMNAME` header.
Dark: chain-level (filter "black" once, 6 cameras in parallel, then restore the sweep mode).
Analysis: `ProcessPoolExecutor` (`spawn` context) mapping opticalib's `_analyze_datapath(datapath,
tnfringes, n_psfs=1, angles=[a_k], ...)` over the 6 `cam{k}` folders; `max_workers=1` → in-process.

### [x] Step 3: chain tests (M4) — Done
Result: `test/procedures/test_spl_chain.py`, 57 passed (21 s), no SPLChain bugs found. Real `spawn`
process pool only tested via a fake pool → real-spawn integration test added in Step 4.
Sleep-based fake cameras + recording fake filter: no `move_to` overlaps an acquisition; all cameras
acquire at the same wl; wall time ≈ `N·(T_move+T_exp)`; `set_exptime` overlaps the move; one failing
camera → others saved, aggregate error names it, filter ends at 600; parallel analysis == in-process.

### [x] Step 4: critic (full suites, both repos) + fix cycles (max 2 per failure) — Done
Run 1 (env `m4`): phasing 72 passed; SPLChain 57 passed ×3 (not flaky); real-spawn test 1 passed (10 s);
`m4/spl.py` pyflakes clean. Pre-existing failures, NOT related: opticalib
`test_gui_procedures.py::test_closed_loop_is_an_explicit_loop` (fails on `main` too); M4: 28 failed + 1 error
(OPC UA timeouts, hardware, stale tests: `M4Exapode`, `draw_mask`, `newtc`, packaging deps).
In progress: user-decided changes (720 nm range, `_heal_bad_pixels` zero-hot-pixel fix, NaN cleanup,
`_read_filter_mode` via `get_bandwidth_mode()` in SPL + SPLChain, 2 f-string lints).
Run 2 (final): opticalib 1347 passed, 1 failed (the pre-existing GUI test above); phasing 86 passed;
M4 `test/procedures` 66 passed (incl. real spawn); pyflakes clean on both files.

### [x] Step 5: reviewer on both diffs; Critical items → fixer — Done
Review: no Critical; concurrency design judged sound (tests also pass under heavy CPU load).
Fixing now — opticalib: `acquire(filter_mode="black")` accepted; dark-frame mode restore not in `finally`;
clear error for missing `tnfringes`; template cache keyed on lambda vector. M4: invalidate cached exposure
per sweep / after failures; writer pool = n cameras + early abort on write errors; constructor closes opened
cameras on failure; `SPLChainError.tn`; park-failure semantics; `masks is not None`; per-camera analysis
errors; test flake (unlocked event mutation) + fragile parallel test; spawn config `SYSTEM.version`.
Deferred to user / separate change:
  - M4 `pyproject.toml` pin → DONE: `opticalib>=2.5.5` (user: first release with this branch);
  - `GigaVision.reconnect` (cameras.py): lock around shared `VmbSystem` exit/enter; `_exptime` cache not
    refreshed after reconnect; no overall retry deadline;
  - no timeouts on filter `move_to` / `get_bandwidth_mode` (constructor/sweep could hang) → Step 7.
Result after fixes: phasing 91 passed; opticalib full suite 1353 passed, 0 failed; M4 `test/procedures` 87 passed
(flake-free ×3); pyflakes clean. M4 `test/packaging`: 3 pre-existing failures unrelated to this change
(`arte`/`tqdm` undeclared, `ruamel-yaml` missing, `m4.gui.ott_monitoring_GUI` import).
Note: on constructor failure SPLChain closes only cameras it opened itself, not caller-provided objects.

### [x] Step 6: documenter (opticalib new kwargs / `"none"` filter; `SPLChain` docstrings) → critic if sources edited — Done
Result: docstrings updated in `phasing.py` and `m4/spl.py` (M4: AST-verified docstring-only) + `docs/m4.rst`
(`m4.spl` automodule) + opticalib `docs/configuration.rst` (PHASING keys/defaults). Sphinx not built (not installed).
Found while documenting: the PHASING key `expected_psf_pos` was spelled `expected_psfs_pos` in the packaged
template config and docs stub (so it was silently ignored), and documented as `[y, x]` while the code uses
`[x, y]` (ordering only). Fixed key + docs; template default is now empty, and `detect_psf_centroids` treats
unfilled entries (`- []`) as "no expected positions" (would otherwise IndexError) + 5 tests.
Final run: opticalib 1358 passed, 0 failed; M4 `test/procedures` 87 passed (1 intermittent failure seen once
in 22 runs — name not captured; 21 later runs green, 8 of them under full CPU load); pyflakes clean.

### [ ] Step 7: bench validation (user) — open
Measure `T_move`, `T_exp`; time sequential vs `SPLChain.acquire`; force a camera disconnect mid-sweep
to check `GigaVision.reconnect()` vs the shared `VmbSystem` (if harmful: lock in `cameras.py`, separate change).
