# Full frailty package test suite — compatibility runtime

**Run:** 2026-09-26 04:04 UTC
**Command:** `PYTHONPATH="$PWD/frailty_paper/.venv/lib/python3.13/site-packages" python3 -m pytest frailty_paper/tests -q`
**Runtime:** system Python 3.13 with the project Python 3.13 environment's installed dependencies added to `PYTHONPATH`; the globally installed pytest entry point performed collection.

**Result:** 154 passed, 92 subtests passed in 4.43 seconds.

The standard integrated preflight remains nonzero in pinned Python 3.11.11 because the untracked Brain6-only pilot test imports pytest, which that isolated environment does not contain. This compatibility run confirms the whole discovered frailty test suite passes, but it is not a substitute for adding a reproducible test-runner dependency to the pinned environment. No Brain6 source or test file was changed.
