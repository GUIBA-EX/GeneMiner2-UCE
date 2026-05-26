"""Console entry point for the legacy GeneMiner2 CLI script."""

from __future__ import annotations

import os
from pathlib import Path
import runpy
import sys


def _find_scripts_dir() -> Path:
    candidates = []
    env_path = os.environ.get("GENEMINER2_SCRIPTS")

    if env_path:
        candidates.append(Path(env_path))

    package_path = Path(__file__).resolve()
    candidates.extend([
        package_path.parents[2] / "scripts",
        Path.cwd() / "scripts",
    ])

    for candidate in candidates:
        script = candidate / "unix_command.py"

        if script.is_file():
            return candidate

    raise RuntimeError(
        "Unable to locate GeneMiner2 scripts. Set GENEMINER2_SCRIPTS to the "
        "directory containing unix_command.py."
    )


def main() -> None:
    scripts_dir = _find_scripts_dir()
    script = scripts_dir / "unix_command.py"
    sys.path.insert(0, str(scripts_dir))
    sys.argv[0] = "geneminer2"
    runpy.run_path(str(script), run_name="__main__")


if __name__ == "__main__":
    main()

