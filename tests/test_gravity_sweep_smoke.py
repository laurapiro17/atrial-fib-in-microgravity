import json
import os
import subprocess
import sys

ROOT = os.path.join(os.path.dirname(__file__), "..")


def test_gravity_sweep_smoke_runs_and_writes_results():
    proc = subprocess.run(
        [sys.executable, "experiments/gravity_sweep.py"],
        cwd=ROOT, capture_output=True, text=True, timeout=600,
    )
    assert proc.returncode == 0, proc.stderr
    with open(os.path.join(ROOT, "figures", "results_crn.json")) as f:
        data = json.load(f)
    assert "gravity_law" in data
    rows = data["gravity_law"]["rows"]
    gmap = {r["g"]: r["N_g"] for r in rows}
    # vulnerability rises as gravity falls: N_g(0) > N_g(1)
    assert gmap[0.0] > gmap[1.0]
