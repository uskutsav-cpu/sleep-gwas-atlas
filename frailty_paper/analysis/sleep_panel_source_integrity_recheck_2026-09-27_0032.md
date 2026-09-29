# Locked sleep-panel source integrity recheck — 2026-09-27 00:33 UTC

The current `all_acquired_resources.tsv` contains 12 locked sleep-panel summary-statistic resources. A scoped temporary manifest was generated from those rows and checked with the existing read-only `16_verify_acquired_resource_manifest.py` verifier using the authorized external storage root. Result: `RESOURCE_MANIFEST_OK rows=12 files_verified=12`. Each gzip stream then passed `gzip -t` to EOF, including CRC/trailer validation. No source data, runner state, receipts, or scientific settings were modified.

| Resource | File | Bytes | SHA-256 | Size/SHA | gzip integrity |
|---|---|---:|---|---|---|
| `sleep_panel_insomnia` | `data/raw/insomnia.txt.gz` | 314,778,585 | `32848cd92a6324c9048cad4a804cf1546299af1d0281981c6ca43c2046a9065b` | PASS | PASS |
| `sleep_panel_sleepdur` | `data/raw/sleepdur.txt.gz` | 378,539,229 | `decccd8e4ef5b9ae4d5c210bd7e0db695da7077c3385fae32c5d16c391bee559` | PASS | PASS |
| `sleep_panel_shortsleep` | `data/raw/shortsleep.txt.gz` | 383,448,891 | `4f754eefaa0ff9a415ef3fba20aac38f981b7d58ca198f04e3173903e561e194` | PASS | PASS |
| `sleep_panel_longsleep` | `data/raw/longsleep.txt.gz` | 384,522,668 | `93deae65cdb5a6422b50ebb15457f410a547d77e727cc32d20abb8fd6c16d63c` | PASS | PASS |
| `sleep_panel_chronotype` | `data/raw/chronotype.txt.gz` | 354,578,112 | `9fe6c797afbfad85e50625a6dcd3313af77c0528090c6bcea05959690b9c6327` | PASS | PASS |
| `sleep_panel_sleepiness` | `data/raw/sleepiness.txt.gz` | 382,692,981 | `447335097cb5872cddc28e95e9263749745c4735b038d69874a2e5c35a49bf93` | PASS | PASS |
| `sleep_panel_napping` | `data/raw/napping.txt.gz` | 356,906,540 | `222ada5749b5eab1bb265b691239cfe10d63b4a26bcbc156f1951525f52fda12` | PASS | PASS |
| `sleep_panel_snoring` | `data/raw/snoring.txt.gz` | 258,829,716 | `82eeb068648959afbbccf477745820a2118f88abdff8ff0076ebb47df21f2417` | PASS | PASS |
| `sleep_panel_sleep_apnea` | `data/raw/sleep_apnea.txt.gz` | 763,490,490 | `eb8cfc33febc044c3f608406f5f50c9be7f6973aba5effd89a0f94d1f58c1d4b` | PASS | PASS |
| `sleep_panel_sleep_efficiency` | `data/raw/sleep_efficiency.txt.gz` | 310,070,562 | `e9bcf2716a423ad99bcc1c39cbd2ba6fdff18a3b26f5a06ef6b75e29fd73a462` | PASS | PASS |
| `sleep_panel_accel_sleep_duration` | `data/raw/accel_sleep_duration.txt.gz` | 310,074,893 | `745bd3cd0c7b03e338df3ea7fcdf43cc5045f0688cc32a9aa8fca8da17f2dd0f` | PASS | PASS |
| `sleep_panel_sleep_timing` | `data/raw/sleep_timing.txt.gz` | 310,083,988 | `528c08ccf60c756d8fc9ca63cb108a837871511d0b3d29302e7a71aa7a06644a` | PASS | PASS |

Current acquired-resource registry SHA-256: `9312ba1a4862d4f76b2da45d0f99d0e118e0142a774907dfd8e8b2adc5f424f5` (registry contains 352 rows and 22 required metadata columns; metadata audit output: `analysis/acquisition_manifest_metadata_recheck_2026-09-27_0032.log`).
Verifier script SHA-256: `4be937ef1aca90c877d0aa93eb76649d5692b9f8545a0cdf58dc5875a84d51ed`.

Commands used (the temporary manifest selected every row whose `resource_id` starts with `sleep_panel_`):
```sh
python frailty_paper/scripts/45_audit_acquisition_manifest_metadata.py --manifest frailty_paper/manifests/all_acquired_resources.tsv
python frailty_paper/scripts/16_verify_acquired_resource_manifest.py --repo . --manifest /tmp/current_sleep_panel_manifest.tsv --external-storage-root '/Volumes/Extreme SSD/sleep-gwas-atlas-frailty-v1'
python - <<'PY'
import csv, subprocess
from pathlib import Path
with open("frailty_paper/manifests/all_acquired_resources.tsv", newline="") as stream:
    rows = csv.DictReader(stream, delimiter="\t")
    for row in rows:
        if row["resource_id"].startswith("sleep_panel_"):
            subprocess.run(["gzip", "-t", str(Path(row["file"]).resolve())], check=True)
PY
```
