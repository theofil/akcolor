#!/usr/bin/env python3
"""Idempotent downloader for the Zenodo MOD (CMS Open Data) record family used by
make_mod_friends.py. Queries each record's file list via the Zenodo API and
downloads any file missing or with a mismatched size into MOD/raw/<subdir>/.

Usage:
    python3 download_zenodo.py [sample_name ...]   # default: all samples
"""

import json
import os
import sys
import time
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
RAW_DIR = os.path.join(HERE, "raw")

# sample name -> (Zenodo record id, raw subfolder)
SAMPLES = {
    "CMS2011A_data": (3340205, "data_3340205"),
    "QCD_170-300":   (3341500, "sim_170-300_3341500"),
    "QCD_300-470":   (3341498, "sim_300-470_3341498"),
    "QCD_470-600":   (3341419, "sim_470-600_3341419"),
    "QCD_600-800":   (3364139, "sim_600-800_3364139"),
    "QCD_800-1000":  (3341413, "sim_800-1000_3341413"),
    "QCD_1000-1400": (3341502, "sim_1000-1400_3341502"),
    "QCD_1400-1800": (3341770, "sim_1400-1800_3341770"),
    "QCD_1800-inf":  (3341772, "sim_1800-inf_3341772"),
}


def record_files(record_id):
    url = f"https://zenodo.org/api/records/{record_id}"
    with urllib.request.urlopen(url, timeout=60) as resp:
        meta = json.load(resp)
    return [(f["key"], f["size"], f["links"]["self"]) for f in meta["files"]]


def download_one(key, size, link, dest_dir):
    dest = os.path.join(dest_dir, key)
    if os.path.exists(dest) and os.path.getsize(dest) == size:
        return "skip"
    tmp = dest + ".part"
    for attempt in range(3):
        try:
            urllib.request.urlretrieve(link, tmp)
            if os.path.getsize(tmp) == size:
                os.rename(tmp, dest)
                return "ok"
            print(f"    size mismatch ({os.path.getsize(tmp)} != {size}), retrying...")
        except Exception as e:
            print(f"    attempt {attempt+1} failed: {e}")
        time.sleep(5)
    return "FAILED"


def download_sample(name):
    record_id, subdir = SAMPLES[name]
    dest_dir = os.path.join(RAW_DIR, subdir)
    os.makedirs(dest_dir, exist_ok=True)
    files = record_files(record_id)
    # SIM records contain both GEN*_compressed.h5 (truth-only, no PF candidates --
    # not needed here) and SIM*_compressed.h5 (detector-level, has pfcs/quality/npv,
    # same schema make_mod_friends.py expects). Only fetch the latter; DATA record
    # files (CMS_Jet300_...) have no GEN counterpart and are always kept.
    files = [(k, sz, link) for k, sz, link in files if not k.startswith("GEN")]
    total_gb = sum(sz for _, sz, _ in files) / 1e9
    print(f"[{name}] record {record_id}: {len(files)} files needed (GEN-only files skipped), {total_gb:.1f} GB")
    for i, (key, size, link) in enumerate(files):
        status = download_one(key, size, link, dest_dir)
        print(f"  ({i+1}/{len(files)}) {key} [{size/1e6:.1f} MB] -> {status}")
        if status == "FAILED":
            print(f"  !! giving up on {key}")


def main():
    names = sys.argv[1:] if len(sys.argv) > 1 else list(SAMPLES.keys())
    for name in names:
        download_sample(name)
    print("all done")


if __name__ == "__main__":
    main()
