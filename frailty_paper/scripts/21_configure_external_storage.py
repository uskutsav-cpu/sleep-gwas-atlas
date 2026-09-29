#!/usr/bin/env python3
"""Copy acquired data and optionally route empty derived-output dirs externally."""
from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import sys
from pathlib import Path


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def files(base: Path):
    # exFAT may create AppleDouble resource-fork sidecars while copying from macOS.
    return sorted(p for p in base.rglob("*") if p.is_file() and not p.is_symlink()
                  and not p.name.startswith("._"))


def verify(src: Path, dst: Path) -> int:
    src_files = files(src)
    dst_files = files(dst)
    src_rel = [p.relative_to(src) for p in src_files]
    dst_rel = [p.relative_to(dst) for p in dst_files]
    if src_rel != dst_rel:
        raise RuntimeError(f"file inventory mismatch: {src} -> {dst}")
    for rel in src_rel:
        a, b = src / rel, dst / rel
        if a.stat().st_size != b.stat().st_size or digest(a) != digest(b):
            raise RuntimeError(f"copy verification failed: {rel}")
    return len(src_rel)


def verify_preserved_subset(src: Path, dst: Path) -> int:
    """Verify the initial local snapshot remains intact within an expanded tree."""
    src_files = files(src)
    for source in src_files:
        rel = source.relative_to(src)
        target = dst / rel
        if not target.is_file() or source.stat().st_size != target.stat().st_size or digest(source) != digest(target):
            raise RuntimeError(f"preserved source differs or is missing in external tree: {rel}")
    return len(src_files)


def check_empty_placeholder_tree(path: Path) -> None:
    """Allow only the tracked empty .gitkeep in a derived-output directory."""
    if not path.is_dir() or path.is_symlink():
        raise RuntimeError(f"derived output path is not a plain directory: {path}")
    entries = list(path.iterdir())
    if any(entry.name != ".gitkeep" or not entry.is_file() or entry.stat().st_size != 0
           for entry in entries):
        raise RuntimeError(f"refusing to redirect non-empty derived output directory: {path}")


def prepare_output_links(root: Path, output_root: Path) -> list[tuple[Path, Path]]:
    specs = [(root / "data/harmonized", output_root / "data/harmonized"),
             (root / "data/munged", output_root / "data/munged")]
    if output_root == root or root in output_root.parents or output_root in root.parents:
        raise RuntimeError("repository and derived-output root must be separate")
    if not output_root.is_dir() or not os.access(output_root, os.W_OK):
        raise RuntimeError(f"derived-output root must be an existing writable directory: {output_root}")
    if shutil.disk_usage(output_root).free < 20 * 1024**3:
        raise RuntimeError(f"derived-output root has less than 20 GiB free: {output_root}")

    for src, dst in specs:
        if src.is_symlink():
            if src.resolve() != dst.resolve() or not dst.is_dir():
                raise RuntimeError(f"refusing unrelated or broken derived-output symlink: {src} -> {src.resolve()}")
            backup = src.with_name(src.name + ".local-preserved")
            if not backup.is_dir():
                raise RuntimeError(f"derived-output symlink lacks its preserved local directory: {backup}")
            continue
        check_empty_placeholder_tree(src)
        backup = src.with_name(src.name + ".local-preserved")
        if backup.exists() or backup.is_symlink():
            raise RuntimeError(f"derived-output preservation path already exists: {backup}")
        if dst.exists():
            check_empty_placeholder_tree(dst)
        else:
            dst.parent.mkdir(parents=True, exist_ok=True)
    return specs


def apply_output_links(specs: list[tuple[Path, Path]]) -> None:
    for src, dst in specs:
        if src.is_symlink() and src.resolve() == dst.resolve():
            print(f"derived output already linked: {src} -> {dst}")
            continue
        if not dst.exists():
            dst.mkdir(parents=True)
            placeholder = src / ".gitkeep"
            if placeholder.is_file():
                shutil.copyfile(placeholder, dst / ".gitkeep")
        backup = src.with_name(src.name + ".local-preserved")
        src.rename(backup)
        try:
            src.symlink_to(dst, target_is_directory=True)
        except Exception:
            backup.rename(src)
            raise
        print(f"linked derived output {src} -> {dst}; empty original retained at {backup}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo", required=True, type=Path)
    ap.add_argument("--storage-root", required=True, type=Path)
    ap.add_argument("--derived-output-root", type=Path,
                    help="Optional existing volume directory for empty data/harmonized and data/munged paths")
    ap.add_argument("--apply", action="store_true", help="Perform copy and replace repo data directories with symlinks")
    args = ap.parse_args()
    root = args.repo.resolve()
    storage = args.storage_root.expanduser().absolute()
    if not args.apply:
        ap.error("refusing to mutate without --apply")
    if not storage.is_dir() or not os.access(storage, os.W_OK):
        raise SystemExit(f"storage root is not an existing writable directory: {storage}")
    if storage == root or root in storage.parents or storage in root.parents:
        raise SystemExit("repository and external storage root must be separate")

    output_specs = None
    if args.derived_output_root:
        output_root = args.derived_output_root.expanduser().absolute()
        try:
            output_specs = prepare_output_links(root, output_root)
        except (OSError, RuntimeError) as exc:
            raise SystemExit(f"derived-output preflight failed; no output paths changed: {exc}")

    specs = [(root / "data/raw", storage / "data/raw"),
             (root / "frailty_paper/data", storage / "frailty_paper/data")]
    for src, dst in specs:
        if src.is_symlink():
            if src.resolve() == dst.resolve():
                backup = src.with_name(src.name + ".local-preserved")
                if backup.is_dir() and dst.is_dir():
                    n = verify_preserved_subset(backup, dst)
                    print(f"existing link and preserved source re-verified ({n} files): {src} -> {dst}")
                else:
                    raise SystemExit(f"linked path lacks its preserved source or storage target: {src}")
                continue
            raise SystemExit(f"refusing to replace unrelated symlink: {src} -> {src.resolve()}")
        if not src.is_dir():
            raise SystemExit(f"source data directory missing: {src}")
        if dst.exists():
            try:
                n = verify(src, dst)
                print(f"existing external copy verified ({n} files): {dst}")
            except Exception as exc:
                raise SystemExit(f"external target exists but does not match; preserving both: {exc}")
        else:
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copytree(src, dst, symlinks=True)
            n = verify(src, dst)
            print(f"copied and SHA-256 verified ({n} files): {dst}")

    # Only after both copies have passed full inventory and checksum verification.
    for src, dst in specs:
        if src.is_symlink() and src.resolve() == dst.resolve():
            continue
        backup = src.with_name(src.name + ".local-preserved")
        if backup.exists() or backup.is_symlink():
            raise SystemExit(f"preservation path already exists; no swap performed for {src}: {backup}")
    for src, dst in specs:
        if src.is_symlink() and src.resolve() == dst.resolve():
            continue
        backup = src.with_name(src.name + ".local-preserved")
        src.rename(backup)
        try:
            src.symlink_to(dst, target_is_directory=True)
        except Exception:
            backup.rename(src)
            raise
        print(f"linked {src} -> {dst}; original retained at {backup}")
    if output_specs is not None:
        apply_output_links(output_specs)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
