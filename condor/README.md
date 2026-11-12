# HTCondor scripts and run registry (proof cluster)

Submit from the HTCondor submit host (proof-11). Every run is registered so that results can be
compared later with the exact code, environment and configuration that produced them.

- `sampling_job.sh` – one job: generates the macro with `sampling/make_macro.py`, runs `DREMTubes`
  (1 thread) in `<runs>/<tag>/work/E<E>_<job>/`, moves the ntuple to `<runs>/<tag>/sf_E<E>_<job>.root`
  and writes `<runs>/<tag>/jobinfo/E<E>_<job>.json` (host, OS, timing, exit code, Geant4 version from
  the banner, linked libG4run, executable md5, data sets).
  The environment is `<builddir>/hidra_env.sh` (templates in `env/`), so each build carries its own.
- `sampling.sub` – submit description; one job per line of a joblist `energy iE job nevents`
  (seed1 = 1000*iE + job + 1). Variables: `tag joblist build hidra runsdir particle physlist`.
- `submit_run.sh <tag> <joblist> [note]` – refuses an existing tag, records code/environment/configuration
  in `<runs>/<tag>/run_info.json`, submits, stores the cluster id. Env: `BUILD`, `PARTICLE`, `PHYSLIST`, `HIDRA_RUNS`.
- `close_run.sh <tag> [note]` – after the cluster has left the queue: `hadd` per energy into
  `<runs>/<tag>/merged/`, records events, timing and failed jobs, status `done`.
- `registry.py` – the registry itself (`new`, `set`, `close`, `render`); writes `<runs>/registry.jsonl`
  (one line per run) and `<runs>/RUNS.md`.
- `joblists/` – joblists of the runs submitted so far; `env/` – environment of each build.

`<runs>` defaults to `<repo>/runs` (git-ignored by not being added; override with `HIDRA_RUNS`).
