# OmicVerse Submodule Policy

SCOOP keeps Fast modules as the primary implementation. OmicVerse is available
only as a hidden optional reference path.

Repository layout:

```text
third_party/omicverse/        # git submodule, pinned upstream OmicVerse
src/scsp_agent_sop/omicverse_facilities.py
```

`scsp_agent_sop.omicverse_facilities` is the only default bridge that may import
OmicVerse. It first tries an installed `omicverse` package; if unavailable, it
can temporarily add `third_party/omicverse` to `sys.path` and import the
submodule package.

Default behavior:

- FastCore uses `fastcore_cpu`, `fastcore_mixed`, and `fastcore_oom` names.
- `omicverse_*` backend names are accepted only as legacy aliases.
- `scoop_fast.registry` does not expose OmicVerse engines to agents.
- `configs/default_run.yaml` keeps OmicVerse reference disabled by default.
- The Fast environment does not install OmicVerse.

Reference use cases:

- cNMF validation when FastCNMF stability is insufficient.
- marker/CCC compatibility checks when explicitly requested.
- external FastCore mixed/OOM smoke tests in a separate reference environment.

This keeps OmicVerse reproducible and traceable without making it the visible
default engine.
