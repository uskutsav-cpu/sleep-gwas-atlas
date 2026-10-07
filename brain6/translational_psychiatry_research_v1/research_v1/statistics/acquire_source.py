"""Acquire public pinned method source; bounded downloads, with manifest before fetch."""
import datetime, hashlib, json, pathlib, subprocess

HERE = pathlib.Path(__file__).resolve().parent
COMMIT = "319e84e114a4f954005a4592756c56cfee083667"
FILES = ["calculate.py", "heritability.py", "ldsc_thin.py", "pheno.py", "prep.py", "supergnova.py", "README.md", "LICENSE"]

def main():
    entries = [{"file": f, "url": f"https://raw.githubusercontent.com/qlu-lab/SUPERGNOVA/{COMMIT}/{f}", "max_bytes": 1000000} for f in FILES]
    (HERE / "source_acquisition_manifest.json").write_text(json.dumps({"created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(), "commit":COMMIT, "entries":entries}, indent=2)+"\n")
    receipts = []
    for item in entries:
        out = HERE / "source" / item["file"]
        p = subprocess.run(["curl", "--fail", "--location", "--max-time", "30", "--max-filesize", str(item["max_bytes"]), "--output", str(out), item["url"]], capture_output=True, text=True)
        record = {**item,"exit_code":p.returncode,"stderr":p.stderr,"retrieved_utc":datetime.datetime.now(datetime.timezone.utc).isoformat()}
        if p.returncode == 0:
            record.update(bytes=out.stat().st_size,sha256=hashlib.sha256(out.read_bytes()).hexdigest())
        receipts.append(record)
        (HERE / "source_acquisition_receipts.json").write_text(json.dumps(receipts, indent=2)+"\n")
        print(item["file"],p.returncode,record.get("sha256","UNAVAILABLE"),flush=True)

if __name__ == "__main__": main()
