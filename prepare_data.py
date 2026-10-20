# Table: public coordinate datasets in README.md.
from pathlib import Path
import gzip
import hashlib
import json
import os
import shutil
import struct
import time
import urllib.request
import zipfile

import numpy as np


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
RAW = DATA / "raw"
RAW.mkdir(parents=True, exist_ok=True)
RADIUS = 6371008.8
GSHHG_HASH = "28600e8f7a08645aab43079326df6504212ec5ccb2b4bcf3b5f4f12ed60e82bc"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(2 ** 20), b""):
            digest.update(block)
    return digest.hexdigest()


def obtain(url, filename, expected=None, existing=None):
    target = RAW / filename
    if not target.exists():
        temporary = target.with_suffix(target.suffix + ".part")
        if existing and Path(existing).is_file():
            shutil.copyfile(existing, temporary)
        else:
            for attempt in range(3):
                try:
                    request = urllib.request.Request(url, headers={"User-Agent": "CGAL-Delaunay-reproduction/1.0"})
                    with urllib.request.urlopen(request, timeout=120) as response, temporary.open("wb") as output:
                        shutil.copyfileobj(response, output)
                    break
                except OSError:
                    if attempt == 2:
                        raise
                    time.sleep(5 * (attempt + 1))
        if expected and sha256(temporary) != expected:
            raise RuntimeError(f"Source checksum mismatch: {filename}")
        temporary.replace(target)
    digest = sha256(target)
    if expected and digest != expected:
        raise RuntimeError(f"Source checksum mismatch: {filename}")
    return target, {"url": url, "file": str(target.relative_to(ROOT)), "sha256": digest,
                    "bytes": target.stat().st_size}


def project(coordinates, center):
    longitude = np.deg2rad(coordinates[:, 0])
    latitude = np.deg2rad(coordinates[:, 1])
    longitude0, latitude0 = np.deg2rad(center)
    delta = longitude - longitude0
    sine, cosine = np.sin(latitude), np.cos(latitude)
    denominator = 1 + np.sin(latitude0) * sine + np.cos(latitude0) * cosine * np.cos(delta)
    if np.any(denominator <= 0) or not np.isfinite(coordinates).all():
        raise RuntimeError("Projection domain failure")
    scale = np.sqrt(2 / denominator)
    output = np.empty_like(coordinates, dtype="<f8")
    output[:, 0] = RADIUS * scale * cosine * np.sin(delta)
    output[:, 1] = RADIUS * scale * (np.cos(latitude0) * sine - np.sin(latitude0) * cosine * np.cos(delta))
    return output


def store(name, coordinates, center, metadata):
    points = project(coordinates, center)
    target = DATA / (name + ".xy")
    temporary = target.with_suffix(".xy.part")
    points.tofile(temporary)
    temporary.replace(target)
    return {"name": name, "n": len(points), "file": str(target.relative_to(ROOT)),
            "sha256": sha256(target), "bytes": target.stat().st_size,
            "projection": "Spherical Lambert azimuthal equal-area", "center_lon_lat_degrees": center,
            "radius_metres": RADIUS, "metric": "Euclidean distance in the projected plane, metres",
            "duplicates": "Preserved", "bounds_lon_lat_degrees": [coordinates.min(axis=0).tolist(),
                                                                     coordinates.max(axis=0).tolist()],
            **metadata}


def roads(code, count, center):
    filename = f"USA-road-d.{code}.co.gz"
    url = "https://www.diag.uniroma1.it/~challenge9/data/USA-road-d/" + filename
    source, record = obtain(url, filename)
    coordinates = np.empty((count, 2), dtype=np.float64)
    found = 0
    declared = None
    with gzip.open(source, "rt", encoding="ascii") as stream:
        for line in stream:
            fields = line.split()
            if not fields or fields[0] == "c":
                continue
            if fields[0] == "p":
                declared = int(fields[-1])
                continue
            if fields[0] != "v" or len(fields) != 4 or int(fields[1]) != found + 1:
                raise RuntimeError("Unexpected coordinate record or node identifier")
            coordinates[found] = [int(fields[2]) * 1e-6, int(fields[3]) * 1e-6]
            found += 1
    if found != count or declared != count:
        raise RuntimeError("DIMACS coordinate count mismatch")
    return store("roads_" + code.lower(), coordinates, center, {
        "source": record, "credit": "9th DIMACS Implementation Challenge, TIGER/Line US road networks",
        "source_description": "https://www.diag.uniroma1.it/~challenge9/download.shtml",
        "license": "Public domain, as stated on the DIMACS download page",
        "selection": "All coordinate records; graph edges are not used"})


def shorelines():
    name = "gshhg-bin-2.3.7.zip"
    archive, record = obtain("https://ftp.soest.hawaii.edu/gshhg/" + name, name,
                             GSHHG_HASH, os.environ.get("GSHHG_ARCHIVE"))
    selected = []
    global_count = 0
    with zipfile.ZipFile(archive) as source:
        for member in source.namelist():
            if member.upper().endswith("LICENSE.TXT"):
                (DATA / "GSHHG_LICENSE.txt").write_bytes(source.read(member))
        with source.open("gshhs_f.b") as stream:
            while header := stream.read(44):
                if len(header) != 44:
                    raise RuntimeError("Truncated GSHHG header")
                fields = struct.unpack(">11i", header)
                count, level = fields[1], fields[2] & 255
                raw = stream.read(8 * count)
                if len(raw) != 8 * count or count < 1:
                    raise RuntimeError("Truncated GSHHG coordinates")
                if level not in {1, 2, 3, 4, 5}:
                    continue
                global_count += count
                coordinates = np.frombuffer(raw, dtype=">i4").reshape(-1, 2).astype(np.float64) * 1e-6
                coordinates[:, 0] = (coordinates[:, 0] + 180) % 360 - 180
                keep = ((coordinates[:, 0] >= -25) & (coordinates[:, 0] <= 45)
                        & (coordinates[:, 1] >= 30) & (coordinates[:, 1] <= 72))
                if keep.any():
                    selected.append(coordinates[keep])
    if global_count != 10587854:
        raise RuntimeError("GSHHG source count mismatch")
    return store("shorelines_europe", np.concatenate(selected), [10.0, 52.0], {
        "source": record, "credit": "Paul Wessel and Walter H. F. Smith; GSHHG 2.3.7",
        "source_description": "https://ftp.soest.hawaii.edu/gshhg/",
        "license": "GNU Lesser General Public License version 3 or later",
        "selection": "gshhs_f.b, levels 1-5, longitude [-25,45], latitude [30,72], inclusive",
        "source_level_1_to_5_count": global_count})


def main():
    manifest = DATA / "datasets.json"
    if manifest.exists():
        data = json.loads(manifest.read_text())
        if all(sha256(ROOT / x["file"]) == x["sha256"] for x in data["datasets"]):
            print(json.dumps({"datasets_reused": len(data["datasets"])}), flush=True)
            return
        raise RuntimeError("Frozen dataset checksum mismatch")
    datasets = []
    for code, count, center in [("NY", 264346, [-74.0, 40.8]),
                                 ("CAL", 1890815, [-119.5, 37.25]),
                                 ("W", 6262104, [-115.0, 38.5])]:
        datasets.append(roads(code, count, center))
        print(json.dumps({"prepared": datasets[-1]["name"], "n": datasets[-1]["n"]}), flush=True)
    datasets.append(shorelines())
    data = {"format": "Little-endian float64 pairs x,y, no header", "datasets": datasets}
    temporary = manifest.with_suffix(".json.part")
    temporary.write_text(json.dumps(data, indent=2) + "\n")
    temporary.replace(manifest)
    print(json.dumps({"prepared": datasets[-1]["name"], "n": datasets[-1]["n"]}), flush=True)


if __name__ == "__main__":
    main()
