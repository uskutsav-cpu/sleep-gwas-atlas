# Clean committed frailty test-suite verification — 2026-09-25 06:42 UTC

## Scope and result

Ran the package's `unittest discover -s tests -v` command in a temporary archive of committed `HEAD` `e0db1771c7de48e7f17a67bc43a390a38d08c9d8`, using the isolated Python 3.11.11 environment at `work/conda-envs/frailty-paper-py311-clean-2026-09-23`. All 130 tracked frailty tests passed (1.262 seconds). The test includes the incremental PubMed-window builder and its invariant that later-window records stay separate, unscreened, and do not mutate the frozen queue.

The archive contained only committed files; eight absolute/escaping external-volume symlinks were omitted during safe extraction. No repository or external analysis inputs were modified by this verification.

## Mixed-worktree distinction

The full test discovery in the shared working tree reports one import error in the untracked `frailty_paper/tests/test_audit_brain6_sleep_power_pilot.py`, which imports `pytest`; that package is absent from both the canonical Python 3.13.13 venv and isolated Python 3.11.11 venv. The other 130 tests pass in-place when that untracked Brain6-specific module is excluded. No Brain6 file or dependency was changed as part of this audit.
