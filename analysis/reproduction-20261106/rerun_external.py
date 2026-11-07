#!/usr/bin/env python3
"""Rerun the existing isolated CRNT4SBML worker and retain its SBML inputs."""

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--synkit-root", type=Path, required=True)
    parser.add_argument("--legacy-python", type=Path, required=True)
    parser.add_argument("--native-library", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.synkit_root))
    from Experiment.CRN.benchmarks import BENCHMARK_NETWORKS
    from Experiment.CRN.studies.external_crnt4sbml import compare_payload
    from Experiment.CRN.tools.common import environment
    from synkit.CRN import crnt_summary, native_backend, write_sbml

    out = args.output_dir.resolve()
    (out / "sbml").mkdir(parents=True)
    baseline = {row["name"]: row for row in json.loads(args.baseline.read_text())["networks"]}
    ours = {}
    manifest = {"networks": []}
    with native_backend(args.native_library):
        for fixture in BENCHMARK_NETWORKS:
            model = fixture.build()
            row = crnt_summary(model).to_dict()
            row["weakly_reversible"] = row.pop("is_weakly_reversible")
            path = out / "sbml" / (fixture.name + ".xml")
            write_sbml(model, path, model_id=fixture.name, model_name=fixture.description)
            row["input_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
            assert row["input_sha256"] == baseline[fixture.name]["input_sha256"], fixture.name
            ours[fixture.name] = row
            manifest["networks"].append({"name": fixture.name,
                                         "path": str(path.relative_to(out))})
        env_record = environment()

    worker = args.synkit_root / "Experiment/CRN/tools/crnt4sbml_worker.py"
    # Keep the worker with the portable, relative-path input manifest.
    (out / "crnt4sbml_worker.py").write_bytes(worker.read_bytes())
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    command = [str(args.legacy_python), str(out / "crnt4sbml_worker.py"),
               "manifest.json", "crnt4sbml-raw.json"]
    (out / "command.json").write_text(json.dumps(command, indent=2) + "\n")
    env = dict(os.environ, MPLBACKEND="Agg", MPLCONFIGDIR=str(out / "matplotlib"))
    with (out / "worker.log").open("w") as log:
        result = subprocess.run(command, cwd=out, env=env,
                                stdout=log, stderr=subprocess.STDOUT, timeout=120)
    if result.returncode:
        raise RuntimeError((out / "worker.log").read_text())
    external = json.loads((out / "crnt4sbml-raw.json").read_text())
    comparison = compare_payload(ours, external)
    assert len(comparison["networks"]) == 12
    assert all(comparison["checks"].values()), comparison
    assert all(row["crnt4sbml"] == baseline[row["name"]]["crnt4sbml"]
               for row in comparison["networks"])
    report = {"schema": "synkit.crn-external-reproduction/1", "status": "PASS",
              "external_execution": {"interpreter_rerun": True,
                                     "python_executable": str(args.legacy_python),
                                     "all_inputs_match_paper_sha256": True,
                                     "worker_sha256": hashlib.sha256(worker.read_bytes()).hexdigest()},
              "environment": env_record, "external_environment": external["environment"],
              **comparison}
    (out / "comparison.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print("PASS: fresh isolated CRNT4SBML run agrees on 12 networks and 8 shared quantities")


if __name__ == "__main__":
    main()
