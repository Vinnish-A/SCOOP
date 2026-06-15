# SCOOP Skill Packages

SCOOP skills are versioned biological rule packages. They are not free-form
prompts and they must not mutate H5AD directly.

A future skill package should use this layout:

```text
skills/<skill_id>/
  skill.yaml
  markers.tsv
  anti_markers.tsv
  states.yaml
  programs.yaml
  naming_rules.yaml
  ontology_map.tsv
  conflict_rules.yaml
  examples/
```

The first architecture implementation only defines loader models and schemas.
It does not ship a large marker database.
