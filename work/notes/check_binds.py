"""Check every lc task's binds against the host paths cat_config resolves for it."""
import re, sys, yaml
from pathlib import Path
from lightcone.engine import plan

ROOT = Path(".")
cc = yaml.safe_load(open("src/cosmo_val/cat_config.yaml"))
cfg = yaml.safe_load(open("src/papers/bmodes/config/config.yaml"))
SWEEP = [v for v in cfg["versions"] if "_ecut" not in v and v != cfg["fiducial"]["version"]]

def resolve(ver):
    key = ver if ver in cc else re.sub(r"_leak_corr$", "", ver)
    e = cc[key]; sub = Path(e["subdir"]); d = e.get("shear", {})
    out = {}
    for k in ("path", "redshift_path"):
        if k in d:
            p = Path(d[k]); out[k] = Path(p if p.is_absolute() else sub / p)
    return out

# which reads each script performs
CATALOG = ("run_2pcf.py", "generate_pseudo_cl.py", "generate_pseudo_cl_cov.py")
NZ = ("generate_pseudo_cl_cov.py", "precompute_pure_eb_chunk.py",
      "run_pure_eb_sweep.sh", "run_pure_eb_semianalytic.sh", "run_pure_eb_ptes_sweep.sh")

def cosmocov_nz(version, blind="A"):
    """generate_cosmocov_ini.build_redshift_path — CosmoCov composes its own
    n(z) path from the version string; it does NOT use cat_config's."""
    b = re.sub(r"_leak_corr$", "", version); b = re.sub(r"_ecut\d+", "", b)
    if "v1.4.11" in b: b = "SP_v1.4.6"
    return Path(f"/n17data/sguerrini/UNIONS/WL/nz/{b.replace('SP_','')}/nz_{b}_{blind}.txt")

def versions_of(recipe):
    m = re.findall(r"--ver(?:sion)?\s+(\S+)", recipe)
    if m: return sorted(set(m))
    if any(s in recipe for s in ("run_pure_eb_sweep.sh", "run_pure_eb_ptes_sweep.sh")):
        return SWEEP
    return []

g = plan.build(ROOT)
bad = []
checked = 0
for (u, o), t in sorted(g.tasks.items()):
    r = t.recipe
    binds = [Path(b).resolve() for p in t.inputs.values()
             if (b := plan.readable_source(p)) is not None]
    def granted(target):
        target = Path(target).resolve()
        return any(target == b or b in target.parents for b in binds)
    want = []
    for v in versions_of(r):
        res = resolve(v)
        if any(s in r for s in CATALOG) and "path" in res:
            want.append((v, "catalog", res["path"]))
        if any(s in r for s in NZ) and "redshift_path" in res:
            want.append((v, "n(z)", res["redshift_path"]))
        if "run_cosmocov_chain.sh" in r:
            want.append((v, "n(z) [cosmocov ini]", cosmocov_nz(v)))
    if "run_cosmocov_chain.sh" in r:
        want.append(("-", "cosmocov binary",
                     Path("/n23data1/n06data/lgoh/scratch/UNIONS/CosmoCov/covs/cov")))
        want.append(("-", "libgsl", Path("/lib64/libgsl.so.23.1.0")))
    for v, kind, path in want:
        checked += 1
        if not granted(path):
            bad.append((o, v, kind, str(path)))
print(f"tasks={len(g.tasks)} path-checks={checked} UNCOVERED={len(bad)}")
for b in bad:
    print("   ", b)
sys.exit(1 if bad else 0)
