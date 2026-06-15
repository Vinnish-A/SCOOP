from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class NMFTask:
    """One cNMF replicate factorization task."""

    task_id: str
    k: int
    iteration: int
    seed: int
    status: str = "pending"


@dataclass(frozen=True)
class TaskManifest:
    """Serializable collection of NMF replicate tasks."""

    run_name: str
    tasks: tuple[NMFTask, ...]

    def to_json(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "run_name": self.run_name,
            "tasks": [asdict(task) for task in self.tasks],
        }
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    @classmethod
    def from_json(cls, path: Path) -> "TaskManifest":
        payload = json.loads(path.read_text(encoding="utf-8"))
        return cls(
            run_name=payload["run_name"],
            tasks=tuple(NMFTask(**task) for task in payload["tasks"]),
        )


def build_nmf_tasks(
    run_name: str,
    k_values: tuple[int, ...],
    n_iter: int,
    seed: int,
) -> TaskManifest:
    """Build deterministic NMF replicate tasks for dynamic scheduling."""

    tasks: list[NMFTask] = []
    for k in k_values:
        for iteration in range(n_iter):
            task_seed = seed + (k * 1_000_003) + iteration
            tasks.append(
                NMFTask(
                    task_id=f"k{k}_iter{iteration}",
                    k=k,
                    iteration=iteration,
                    seed=task_seed,
                )
            )
    return TaskManifest(run_name=run_name, tasks=tuple(tasks))

