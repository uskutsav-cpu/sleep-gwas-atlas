#!/usr/bin/env python3
"""Execute only pairs admitted by the pre-outcome source eligibility lock."""

import argparse
import hashlib
import json
import pathlib
import subprocess
import sys


REPO = pathlib.Path(__file__).resolve().parents[2]
RESULTS = REPO / "brain6/results/brain6_alternative_local_validation_v1"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair", choices=["insomnia__adhd", "insomnia__mdd"])
    args = parser.parse_args()
    lock = json.loads((RESULTS / "source_eligibility_lock_v1.json").read_text())
    protocol = RESULTS / "protocol_freeze.json"
    if lock["parent_protocol_sha256"] != sha256(protocol):
        raise RuntimeError("Source eligibility lock is detached from its frozen protocol")
    if set(lock["pair_status"]) != set(json.loads(protocol.read_text())["pairs"]):
        raise RuntimeError("Eligibility lock does not cover exactly the frozen pair family")
    pairs = [args.pair] if args.pair else list(lock["pair_status"])
    for pair in pairs:
        decision = lock["pair_status"][pair]
        if decision == "METHOD_INAPPLICABLE":
            print(f"METHOD_INAPPLICABLE {pair} {lock['inapplicable_reason']}", flush=True)
            continue
        if decision != "ELIGIBLE_SECONDARY":
            raise RuntimeError(f"Unknown pair eligibility: {pair}={decision}")
        subprocess.run([sys.executable,
                        str(REPO / "brain6/scripts/run_supergnova_family_v1.py"),
                        "--stage", "run", "--pair", pair], check=True)


if __name__ == "__main__":
    main()
