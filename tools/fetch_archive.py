"""
fetch_archive.py — pull-only docker-archive assembly for lightcone (lc), no podman.

Fetches an OCI/Docker image's manifest + layer blobs from a registry (ghcr.io,
using an anonymous pull token) in pure Python and writes them out as a
docker-archive tarball at the path `lc` expects for a pulled image
(`.datalad/environments/lc-env-<tag>/image`). Optionally grafts a second
image's `/uv` binary onto the base as an extra layer plus a
`UV_PYTHON_DOWNLOADS=never` env var, mimicking what `lc build`'s
`COPY --from=<uv image> /uv /usr/local/bin/uv` step does — this lets a
pull-only path stand in for `lc build` when `[tool.lightcone.image]` is
base-only (e.g. `ghcr.io/cosmostat/sp_validation`) and no podman/builder host
is available.

Usage: fetch_archive.py <base-repo> <base-ref> <out-path> [--uv <uv-repo> <uv-ref>]

Verified working invocation (candide e2e dogfood, 2026-08-27), pulling the uv
image as its own base to test the grafting logic before pointing base-repo at
the real sp_validation image:

    cd /tmp/lc-e2e/toy && python3 /tmp/lc-e2e/fetch_archive.py astral-sh/uv \
        "sha256:e5b65587bce7de595f299855d7385fe7fca39b8a74baa261ba1b7147afa78e58" \
        .datalad/environments/lc-env-b38c817c6793d767/image \
        --uv astral-sh/uv \
        "sha256:e85be844203885286c60ffad8a858d48afb6c5a5c237ca0e67f12e74b8f174b1" \
        2>&1 | tail -3

(An earlier attempt at the same command using a tagged ref
`0.12.5@sha256:...` for --uv failed with a 404 — ghcr wants a bare digest or
tag, not the `tag@digest` form; the retry above with a bare digest succeeded,
producing a 119 MB archive.)

Recovered verbatim (two Write-tool versions found; this is the later, final
one) from the Claude Code session transcript
3f074fd6-88dd-44c6-bdb1-fae0aabfc9c3.jsonl (2026-08-27/28), after the
original /tmp/lc-e2e/ working copy was lost to tmp cleanup. No Edit-tool
diffs were found on top of this version in that transcript.
"""
import json, tarfile, io, gzip, hashlib, urllib.request, sys, pathlib

def token(repo):
    return json.load(urllib.request.urlopen(
        f"https://ghcr.io/token?scope=repository:{repo}:pull"))["token"]

ACC = ("application/vnd.oci.image.index.v1+json,application/vnd.docker.distribution.manifest.list.v2+json,"
       "application/vnd.oci.image.manifest.v1+json,application/vnd.docker.distribution.manifest.v2+json")

def fetch(repo, ref):
    tok = token(repo)
    def get(path, accept=ACC):
        r = urllib.request.Request(f"https://ghcr.io/v2/{repo}/{path}",
                                   headers={"Authorization": f"Bearer {tok}", "Accept": accept})
        return urllib.request.urlopen(r).read()
    man = json.loads(get(f"manifests/{ref}"))
    if "manifests" in man:
        d = next(m["digest"] for m in man["manifests"]
                 if m.get("platform", {}).get("architecture") == "amd64"
                 and m["platform"]["os"] == "linux")
        man = json.loads(get(f"manifests/{d}"))
    cfg = json.loads(get(f"blobs/{man['config']['digest']}", "application/octet-stream"))
    layers = []
    for i, l in enumerate(man["layers"]):
        blob = get(f"blobs/{l['digest']}", "application/octet-stream")
        print(f"  {repo} layer {i+1}/{len(man['layers'])}: {len(blob)//2**20} MB", flush=True)
        layers.append(blob)
    return cfg, layers

def main():
    base_repo, base_ref, out = sys.argv[1], sys.argv[2], pathlib.Path(sys.argv[3])
    uv_spec = None
    if "--uv" in sys.argv:
        i = sys.argv.index("--uv"); uv_spec = (sys.argv[i+1], sys.argv[i+2])

    cfg, layers = fetch(base_repo, base_ref)
    extra = []  # (uncompressed_tar_bytes,)
    if uv_spec:
        ucfg, ulayers = fetch(*uv_spec)
        # find /uv in the uv image's layers (distroless: one layer, /uv at root)
        uv_bin = None
        for blob in ulayers:
            data = gzip.decompress(blob) if blob[:2] == b"\x1f\x8b" else blob
            with tarfile.open(fileobj=io.BytesIO(data)) as tf:
                for m in tf.getmembers():
                    if m.name.lstrip("./") == "uv" and m.isfile():
                        uv_bin = tf.extractfile(m).read()
        assert uv_bin, "no /uv found in uv image"
        buf = io.BytesIO()
        with tarfile.open(fileobj=buf, mode="w") as tf:
            for d in ("usr", "usr/local", "usr/local/bin"):
                ti = tarfile.TarInfo(d); ti.type = tarfile.DIRTYPE; ti.mode = 0o755
                tf.addfile(ti)
            ti = tarfile.TarInfo("usr/local/bin/uv"); ti.size = len(uv_bin); ti.mode = 0o755
            tf.addfile(ti, io.BytesIO(uv_bin))
        extra.append(buf.getvalue())
        cfg["history"] = cfg.get("history", []) + [{"created_by": "lc-pull-only: graft uv"}]
        cfg["rootfs"]["diff_ids"] += ["sha256:" + hashlib.sha256(extra[0]).hexdigest()]
        env = [e for e in cfg["config"].get("Env", []) if not e.startswith("UV_PYTHON_DOWNLOADS=")]
        cfg["config"]["Env"] = env + ["UV_PYTHON_DOWNLOADS=never"]

    cfg_bytes = json.dumps(cfg).encode()
    cfg_name = hashlib.sha256(cfg_bytes).hexdigest() + ".json"
    out.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(out, "w") as tf:
        def add(name, data):
            ti = tarfile.TarInfo(name); ti.size = len(data)
            tf.addfile(ti, io.BytesIO(data))
        add(cfg_name, cfg_bytes)
        names = []
        for i, blob in enumerate(layers):
            n = f"layer{i}.tar.gz"; names.append(n); add(n, blob)
        for j, ub in enumerate(extra):
            n = f"extra{j}.tar"; names.append(n); add(n, ub)
        add("manifest.json", json.dumps([{"Config": cfg_name, "Layers": names}]).encode())
    print(out, out.stat().st_size // 2**20, "MB; image id:", cfg_name[:-5])

main()
