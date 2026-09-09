#!/usr/bin/env python3
"""Run via brain6 run-job; writes only to BRAIN6_JOB_DIR."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from brain6.native import execute_adapter
if __name__ == "__main__":
    if len(sys.argv)!=2:
        raise SystemExit("Usage: native_cli.py settings.json")
    execute_adapter(sys.argv[1])
