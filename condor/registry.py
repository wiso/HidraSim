#!/usr/bin/env python3
"""Registry of HidraSim simulation runs (stdlib only).

Layout (RUNS_DIR defaults to <repo>/runs, override with $HIDRA_RUNS):
    RUNS_DIR/<tag>/run_info.json   full record of one run
    RUNS_DIR/registry.jsonl        one JSON line per run (summary, rewritten on update)
    RUNS_DIR/RUNS.md               human-readable table + one section per run

Commands:
    registry.py new <tag> <joblist> [--note TEXT] [--setup FILE] [--physlist PL]
        collect code/environment/configuration info, write run_info.json (status "created")
    registry.py set <tag> key=value [key=value ...]
        set top-level fields (e.g. cluster_id=10473 status=running submitted=...)
    registry.py close <tag> [--merged FILE ...] [--note TEXT]
        scan job outputs/logs, record per-energy events, timing, failed jobs; status "done"
    registry.py render
        regenerate registry.jsonl and RUNS.md from all run_info.json files

The joblist format is the one of condor/sampling.sub: "energy_GeV iE job nevents" per line;
seeds are seed1 = 1000*iE + job + 1, seed2 = 7*seed1 + 13 (sampling/make_macro.py).
"""
import argparse
import datetime as dt
import glob
import hashlib
import json
import os
import re
import shutil
import socket
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get("HIDRA_DIR", os.path.dirname(HERE))
RUNS = os.environ.get("HIDRA_RUNS", os.environ.get("OUTBASE", "/storage_tmp/atlas/turra/hidra/sim"))
EXE = os.environ.get("HIDRA_EXE", os.path.join(REPO, "build", "DREMTubes"))


def now():
    return dt.datetime.now().astimezone().isoformat(timespec="seconds")


def sh(cmd, env=None):
    try:
        return subprocess.run(cmd, shell=True, capture_output=True, text=True, env=env,
                              executable="/bin/bash").stdout.strip()
    except Exception as e:  # noqa
        return f"ERROR: {e}"


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_joblist(path):
    jobs = []
    for line in open(path):
        p = line.split()
        if len(p) >= 4 and not line.lstrip().startswith("#"):
            jobs.append(dict(energy=p[0], ie=int(p[1]), job=int(p[2]), nev=int(p[3])))
    return jobs


def geopar(path):
    """Active (uncommented) constants of DREMTubesGeoPar.hh."""
    txt = open(path).read()
    txt = re.sub(r"/\*.*?\*/", "", txt, flags=re.S)
    txt = re.sub(r"//[^\n]*", "", txt)
    out = {}
    for m in re.finditer(r"const\s+G4(?:int|double|bool)\s+(\w+)(\[[^\]]*\])?\s*=\s*([^;]+);", txt):
        name, arr, val = m.group(1), m.group(2), " ".join(m.group(3).split())
        if arr:
            continue  # module maps are long; recorded via the header md5
        out[name] = val
    return out


def git_info():
    g = f"git -C {REPO}"
    info = dict(commit=sh(f"{g} rev-parse HEAD"), branch=sh(f"{g} rev-parse --abbrev-ref HEAD"),
                describe=sh(f"{g} log -1 --format='%h %ci %s'"),
                upstream=sh(f"{g} rev-parse --abbrev-ref --symbolic-full-name @{{u}} 2>/dev/null"))
    diff = sh(f"{g} diff HEAD")
    info["dirty"] = bool(diff)
    return info, diff


def env_info(setup):
    """Environment as the jobs see it (after sourcing the setup file)."""
    pre = f"source {setup} >/dev/null 2>&1; " if setup and os.path.exists(setup) else ""
    e = dict(setup_file=setup if pre else None,
             setup_contents=open(setup).read() if pre else None)
    e["geant4_config_version"] = sh(pre + "geant4-config --version 2>/dev/null")
    e["geant4_config_prefix"] = sh(pre + "geant4-config --prefix 2>/dev/null")
    e["root_version"] = sh(pre + "root-config --version 2>/dev/null")
    e["python"] = sh(pre + "python3 --version 2>&1")
    e["g4_datasets"] = dict(l.split("=", 1) for l in sh(pre + "env | grep -E '^G4[A-Z]*DATA='").splitlines() if "=" in l)
    if os.path.exists(EXE):
        libs = sh(pre + f"ldd {EXE} | grep -E 'libG4run' | awk '{{print $3}}'")
        e["geant4_linked_lib"] = libs
        cache = os.path.join(os.path.dirname(EXE), "CMakeCache.txt")
        if os.path.exists(cache):
            for l in open(cache):
                for k in ("Geant4_DIR", "CMAKE_CXX_COMPILER", "CMAKE_BUILD_TYPE"):
                    if l.startswith(k + ":"):
                        e["cmake_" + k] = l.split("=", 1)[1].strip()
        if e.get("cmake_CMAKE_CXX_COMPILER"):
            e["compiler"] = sh(f"{pre}{e['cmake_CMAKE_CXX_COMPILER']} --version | head -1")
        # actual Geant4 version of the linked libraries
        m = re.search(r"Geant4-([\d.]+)", e.get("cmake_Geant4_DIR", "")) or None
        ver = sh(pre + f"strings $(ldd {EXE} | awk '/libG4run/{{print $3}}') | grep -m1 -oE 'geant4-[0-9]+-[0-9]+(-patch-[0-9]+)?'")
        e["geant4_linked_version"] = ver or (m.group(1) if m else "")
    e["os"] = sh("grep PRETTY_NAME /etc/os-release | cut -d= -f2 | tr -d '\"'")
    e["kernel"] = sh("uname -r")
    return e


def load(tag):
    with open(os.path.join(RUNS, tag, "run_info.json")) as f:
        return json.load(f)


def save(info):
    d = os.path.join(RUNS, info["tag"])
    tmp = os.path.join(d, ".run_info.json.tmp")
    with open(tmp, "w") as f:
        json.dump(info, f, indent=1)
    os.replace(tmp, os.path.join(d, "run_info.json"))
    render()


def cmd_new(a):
    d = os.path.join(RUNS, a.tag)
    if os.path.exists(os.path.join(d, "run_info.json")):
        sys.exit(f"tag {a.tag} already registered")
    os.makedirs(os.path.join(d, "logs"), exist_ok=True)
    jobs = read_joblist(a.joblist)
    shutil.copy(a.joblist, os.path.join(d, "joblist.txt"))
    git, diff = git_info()
    if diff:
        open(os.path.join(d, "git_diff.patch"), "w").write(diff)
    energies = {}
    for j in jobs:
        e = energies.setdefault(j["energy"], dict(njobs=0, events=0, nev_per_job=set(), seeds=[]))
        e["njobs"] += 1
        e["events"] += j["nev"]
        e["nev_per_job"].add(j["nev"])
        e["seeds"].append(1000 * j["ie"] + j["job"] + 1)
    os.makedirs(os.path.join(d, "macros"), exist_ok=True)
    mk = os.path.join(REPO, "sampling", "make_macro.py")
    extra = a.macro_args or ""
    for E, e in energies.items():
        e["nev_per_job"] = sorted(e["nev_per_job"])
        e["seed1_range"] = [min(e["seeds"]), max(e["seeds"])]
        del e["seeds"]
        j0 = next(j for j in jobs if j["energy"] == E)
        s = 1000 * j0["ie"] + j0["job"] + 1
        sh(f"python3 {mk} --energy {E} --nevents {j0['nev']} --seed {s} --particle {a.particle} {extra} "
           f"> {os.path.join(d, 'macros', f'E{E}.mac')}")
    mac0 = open(os.path.join(d, "macros", f"E{jobs[0]['energy']}.mac")).read()
    geo = dict(re.findall(r"^/tbgeo/(\w+)\s+(.+)$", mac0, flags=re.M))
    pos = re.search(r"^/gps/(?:position|pos/centre)\s+(.+)$", mac0, flags=re.M)
    gps_pos = dict(re.findall(r"^/gps/pos/(\w+)\s+(.+)$", mac0, flags=re.M))
    gp = os.path.join(REPO, "include", "DREMTubesGeoPar.hh")
    info = dict(
        tag=a.tag, status="created", created=now(), submitted=None, finished=None,
        submit_host=socket.getfqdn(), user=os.environ.get("USER"), cluster_id=None,
        note=a.note or "",
        code=dict(git=git, executable=EXE, executable_md5=md5(EXE) if os.path.exists(EXE) else None,
                  executable_mtime=dt.datetime.fromtimestamp(os.path.getmtime(EXE)).astimezone().isoformat(timespec="seconds")
                  if os.path.exists(EXE) else None),
        environment=env_info(a.setup),
        simulation=dict(physics_list=a.physlist, optical_flag="-opt not used (default optical physics)",
                        threads_per_job=1, particle=a.particle,
                        energies_GeV=sorted(energies, key=float), per_energy=energies,
                        njobs=len(jobs), events_total=sum(j["nev"] for j in jobs),
                        seed_scheme="seed1 = 1000*iE + job + 1, seed2 = 7*seed1 + 13 (iE, job from joblist)",
                        macro_generator="sampling/make_macro.py " + extra,
                        beam=dict(type=gps_pos.get("type", "Point"), shape=gps_pos.get("shape"), halfx=gps_pos.get("halfx"),
                                  halfy=gps_pos.get("halfy"), direction="0 0 1", position=pos.group(1) if pos else None),
                        tbgeo=geo, geopar=geopar(gp), geopar_md5=md5(gp)),
        files=dict(run_dir=d, joblist="joblist.txt", macros="macros/E<E>.mac",
                   job_outputs="out/sf_E<E>_<job>.root", job_info="jobinfo/E<E>_<job>.json",
                   logs="logs/", merged=[]),
        result=None,
    )
    save(info)
    print(d)


def parse_val(v):
    try:
        return json.loads(v)
    except ValueError:
        return v


def cmd_set(a):
    info = load(a.tag)
    for kv in a.kv:
        k, v = kv.split("=", 1)
        info[k] = parse_val(v)
    save(info)


def job_records(d):
    """Per-job records: from jobinfo/*.json when present, else parsed from logs/*.out."""
    recs = {}
    for f in glob.glob(os.path.join(d, "jobinfo", "*.json")):
        r = json.load(open(f))
        recs[(str(r["energy"]), int(r["job"]))] = r
    for f in glob.glob(os.path.join(d, "logs", "sf_E*_*.out")) + glob.glob(os.path.join(d, "logs", "E*_*.run.log")):
        m = re.search(r"E([\d.]+)_(\d+)\.(?:out|run\.log)$", f)
        key = (m.group(1), int(m.group(2)))
        if key in recs:
            continue
        txt = open(f, errors="replace").read()
        r = dict(energy=key[0], job=key[1])
        h = re.search(r"host=(\S+)", txt)
        el = re.search(r"elapsed_s=(\d+) nev=(\d+)", txt)
        g4 = re.search(r"Geant4 version Name: (\S+)", txt)
        r.update(host=h.group(1) if h else None, elapsed_s=int(el.group(1)) if el else None,
                 nev=int(el.group(2)) if el else None, geant4=g4.group(1) if g4 else None,
                 exit_code=0 if el else None)
        recs[key] = r
    return recs


def cmd_close(a):
    info = load(a.tag)
    d = info["files"]["run_dir"]
    jobs = read_joblist(os.path.join(d, "joblist.txt"))
    recs = job_records(d)
    per_e, failed, missing = {}, [], []
    for j in jobs:
        key = (j["energy"], j["job"])
        out = os.path.join(d, "out", f"sf_E{j['energy']}_{j['job']}.root")
        if not os.path.exists(out):  # old layout: outputs at the top of the run dir
            out = os.path.join(d, f"sf_E{j['energy']}_{j['job']}.root")
        r = recs.get(key, {})
        e = per_e.setdefault(j["energy"], dict(jobs_ok=0, events=0, cpu_s=0, hosts={}, geant4=set()))
        if os.path.exists(out) and r.get("exit_code", 0) == 0 and r.get("elapsed_s"):
            e["jobs_ok"] += 1
            e["events"] += j["nev"]
            e["cpu_s"] += r["elapsed_s"]
            if r.get("host"):
                e["hosts"][r["host"].split(".")[0]] = e["hosts"].get(r["host"].split(".")[0], 0) + 1
            if r.get("geant4"):
                e["geant4"].add(r["geant4"])
        else:
            (failed if r else missing).append(f"E{j['energy']}_{j['job']}")
    for e in per_e.values():
        e["s_per_event"] = round(e["cpu_s"] / e["events"], 2) if e["events"] else None
        e["geant4"] = sorted(e["geant4"])
    info["result"] = dict(per_energy=per_e, jobs_failed=failed, jobs_missing=missing,
                          events_total=sum(e["events"] for e in per_e.values()),
                          cpu_hours=round(sum(e["cpu_s"] for e in per_e.values()) / 3600, 1))
    info["files"]["merged"] = sorted(set(info["files"].get("merged", []) + (a.merged or [])))
    if a.note:
        info["note"] = (info.get("note", "") + "\n" + a.note).strip()
    info["status"] = "done" if not failed and not missing else "done-with-failures"
    info["finished"] = a.finished or now()
    save(info)
    print(json.dumps(info["result"], indent=1))


def summary(i):
    s = i["simulation"]
    env = i["environment"]
    return dict(tag=i["tag"], status=i["status"], submitted=i.get("submitted"), finished=i.get("finished"),
                cluster_id=i.get("cluster_id"), commit=i["code"]["git"]["commit"][:12],
                branch=i["code"]["git"]["branch"], dirty=i["code"]["git"]["dirty"],
                geant4=env.get("geant4_linked_version"), physics_list=s["physics_list"],
                particle=s["particle"], energies_GeV=s["energies_GeV"],
                events={E: v["events"] for E, v in s["per_energy"].items()},
                events_done=(i.get("result") or {}).get("events_total"),
                tilt=f"{s['tbgeo'].get('horizrot', '?')}/{s['tbgeo'].get('vertrot', '?')}",
                note=i.get("note", ""), path=i["files"]["run_dir"])


def render():
    infos = []
    for f in sorted(glob.glob(os.path.join(RUNS, "*", "run_info.json"))):
        try:
            infos.append(json.load(open(f)))
        except ValueError:
            pass
    infos.sort(key=lambda i: i.get("submitted") or i.get("created") or "")
    with open(os.path.join(RUNS, ".registry.jsonl.tmp"), "w") as f:
        for i in infos:
            f.write(json.dumps(summary(i)) + "\n")
    os.replace(os.path.join(RUNS, ".registry.jsonl.tmp"), os.path.join(RUNS, "registry.jsonl"))
    L = ["# HidraSim runs", "",
         "Generated by `condor/registry.py render` from `<tag>/run_info.json`; do not edit by hand.", "",
         "Seeds: seed1 = 1000*iE + job + 1, seed2 = 7*seed1 + 13. Reserved iE blocks: 1-999 sampling study (weekdays), "
         "2600-2699 tb26sym (TB2026 symmetry), 3000-3099 sampling study weekend runs. Never reuse an iE across runs.", "",
         "Status: `⚠ IN CORSO (non terminata)` = jobs still queued/running; `⚠ TERMINATA, DA CHIUDERE (merge mancante)` = no jobs left, "
         "not merged/closed yet; `done` = merged and closed with close_run.sh.", "",
         "| tag | status | submitted | cluster | commit | Geant4 | phys. list | particle | E [GeV] | events (planned) | tilt h/v | note |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for i in infos:
        s = summary(i)
        ev = sum(s["events"].values())
        L.append(f"| {s['tag']} | {s['status']} | {(s['submitted'] or '')[:16]} | {s['cluster_id'] or ''} | "
                 f"{s['commit'][:7]}{'+dirty' if s['dirty'] else ''} | {s['geant4']} | {s['physics_list']} | "
                 f"{s['particle']} | {', '.join(s['energies_GeV'])} | {ev} | {s['tilt']} | "
                 f"{s['note'].splitlines()[0] if s['note'] else ''} |")
    for i in infos:
        s, e, c, r = i["simulation"], i["environment"], i["code"], i.get("result")
        L += ["", f"## {i['tag']}", "",
              f"- status: **{i['status']}**; created {i.get('created')}, submitted {i.get('submitted')}, finished {i.get('finished')}"
              + (f"; jobs done {i['progress']['jobs_done']}/{i['progress']['jobs_total']} (running {i['progress']['jobs_running']}, idle {i['progress']['jobs_idle']}, held {i['progress']['jobs_held']}) as of {i['progress']['updated']}" if i.get('progress') else ""),
              f"- submitted from {i.get('submit_host')} by {i.get('user')}, HTCondor cluster {i.get('cluster_id')}",
              f"- code: `{c['git']['commit']}` on `{c['git']['branch']}` ({c['git']['describe']}), dirty={c['git']['dirty']}",
              f"- executable: `{c['executable']}` md5 `{c['executable_md5']}` built {c.get('executable_mtime')}",
              f"- Geant4 (linked): **{e.get('geant4_linked_version')}** via `{e.get('cmake_Geant4_DIR')}`; "
              f"geant4-config in PATH: {e.get('geant4_config_version')}; ROOT {e.get('root_version')}; {e.get('compiler')}; {e.get('os')}",
              f"- data sets: {', '.join(os.path.basename(v) for v in sorted(e.get('g4_datasets', {}).values()))}",
              f"- physics: {s['physics_list']}, {s['optical_flag']}, {s['threads_per_job']} thread/job",
              f"- beam: {s['particle']}, {s['beam']}; tbgeo {s['tbgeo']}",
              f"- geometry (DREMTubesGeoPar.hh, md5 {s['geopar_md5'][:8]}): "
              + ", ".join(f"{k}={v}" for k, v in s["geopar"].items()),
              f"- seeds: {s['seed_scheme']}",
              f"- jobs: {s['njobs']}, events planned: {s['events_total']}",
              "", "| E [GeV] | jobs | ev/job | events planned | seed1 range | jobs ok | events done | s/event |",
              "|---|---|---|---|---|---|---|---|"]
        for E in s["energies_GeV"]:
            p = s["per_energy"][E]
            q = (r or {}).get("per_energy", {}).get(E, {})
            L.append(f"| {E} | {p['njobs']} | {','.join(map(str, p['nev_per_job']))} | {p['events']} | "
                     f"{p['seed1_range'][0]}-{p['seed1_range'][1]} | {q.get('jobs_ok', '')} | {q.get('events', '')} | "
                     f"{q.get('s_per_event', '') or ''} |")
        if r:
            L.append("")
            L.append(f"- result: {r['events_total']} events, {r['cpu_hours']} CPU-h; failed jobs: "
                     f"{', '.join(r['jobs_failed']) or 'none'}; missing outputs: {', '.join(r['jobs_missing']) or 'none'}")
        if i["files"].get("merged"):
            L.append(f"- merged files: {', '.join(i['files']['merged'])}")
        L.append(f"- directory: `{i['files']['run_dir']}`")
        if i.get("note"):
            L.append(f"- note: {i['note']}")
    with open(os.path.join(RUNS, ".RUNS.md.tmp"), "w") as f:
        f.write("\n".join(L) + "\n")
    os.replace(os.path.join(RUNS, ".RUNS.md.tmp"), os.path.join(RUNS, "RUNS.md"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("new")
    p.add_argument("tag"); p.add_argument("joblist")
    p.add_argument("--note"); p.add_argument("--particle", default="e+")
    p.add_argument("--physlist", default="FTFP_BERT")
    p.add_argument("--macro-args", default="", help="extra make_macro.py args, e.g. '--horizrot 0 --vertrot 0'")
    p.add_argument("--setup", default=os.environ.get("HIDRA_SETUP", os.path.join(REPO, "setup_proof.sh")))
    p.set_defaults(f=cmd_new)
    p = sub.add_parser("set"); p.add_argument("tag"); p.add_argument("kv", nargs="+"); p.set_defaults(f=cmd_set)
    p = sub.add_parser("close"); p.add_argument("tag"); p.add_argument("--merged", nargs="*")
    p.add_argument("--note"); p.add_argument("--finished"); p.set_defaults(f=cmd_close)
    p = sub.add_parser("render"); p.set_defaults(f=lambda a: render())
    a = ap.parse_args()
    a.f(a)


if __name__ == "__main__":
    main()
