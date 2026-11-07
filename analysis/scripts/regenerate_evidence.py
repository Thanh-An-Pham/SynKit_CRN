#!/usr/bin/env python3
"""Regenerate every evidence file under analysis/ from the SynKit source tree.

Run with the Python environment whose SynKit checkout the paper describes:

    python analysis/scripts/regenerate_evidence.py --work-dir /tmp/synkit-crn-evidence

Stage 1 runs Experiment/CRN/run_all.py (all studies, cached BioModels inputs,
native C++ analysis and rule expansion). Stage 2 adds the CRNT4SBML comparison
on byte-identical SBML inputs, the 31-input structural cross-check against the
separate Python siphon search, and the glycolysis replay. Stage 3 installs the
reports under the names the manuscript scripts read. Sources are hashed before
and after; a change during the run aborts the installation.
"""

import argparse
import hashlib
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import date
from pathlib import Path
from time import perf_counter

PAPER = Path(__file__).resolve().parents[2]
ROOT = Path(os.environ.get("SYNKIT_ROOT", PAPER.parents[1] / "SynEco/SynKit")).resolve()
SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SCRIPTS))

#: run_all report name -> name read by the manuscript scripts.
INSTALL_NAMES = {
    "crn-representation.json": "crn-review-experiments.json",
    "crn-supplementary.json": "crn-supplement-audit.json",
    "crn-glycolysis-donors.json": "crn-si-donor-control.json",
}


def _read(out, name):
    return json.loads((out / name).read_text())


def _write(out, name, value):
    (out / name).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def run_studies(out, native_dir):
    """Stage 1: all run_all studies with cached BioModels inputs."""
    command = [
        sys.executable, str(ROOT / "Experiment/CRN/run_all.py"),
        "--output-dir", str(out), "--summary", str(out / "summary.json"),
        "--with-biomodels", "--offline", "--native-build-dir", str(native_dir),
    ]
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    subprocess.run(command, cwd=ROOT, env=env, check=True)


def external_comparison(out, library):
    """Recompute SynKit quantities and compare with recorded CRNT4SBML outputs."""
    from Experiment.CRN.benchmarks import BENCHMARK_NETWORKS
    from Experiment.CRN.studies.external_crnt4sbml import compare_payload
    from Experiment.CRN.tools.common import environment
    from synkit.CRN import crnt_summary, write_sbml
    import synkit

    archived_path = PAPER / "analysis/reference-1.6.3/crn-external-crnt4sbml.json"
    archived = json.loads(archived_path.read_text())
    archived_rows = {r["name"]: r for r in archived["networks"]}
    ours = {}
    started = perf_counter()
    with tempfile.TemporaryDirectory(prefix="crn-sbml-") as directory:
        for fixture in BENCHMARK_NETWORKS:
            model = fixture.build()
            row = crnt_summary(model, library_path=library).to_dict()
            row["weakly_reversible"] = row.pop("is_weakly_reversible")
            path = Path(directory) / (fixture.name + ".xml")
            write_sbml(model, path, model_id=fixture.name, model_name=fixture.description)
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            assert digest == archived_rows[fixture.name]["input_sha256"], fixture.name
            row["input_sha256"] = digest
            ours[fixture.name] = row
    external = {"networks": [{"name": r["name"], **r["crnt4sbml"]} for r in archived["networks"]]}
    comparison = compare_payload(ours, external)
    assert all(comparison["checks"].values())
    reference_dir = out / "reference-1.6.3"
    reference_dir.mkdir(exist_ok=True)
    (reference_dir / archived_path.name).write_bytes(archived_path.read_bytes())
    report = {
        "schema": archived["schema"], "status": "PASS",
        "environment": environment(),
        "external_environment": archived["external_environment"],
        "external_tool": archived["external_tool"],
        "interface": archived["interface"],
        "seconds": perf_counter() - started,
        "external_execution": {
            "mode": "recorded_outputs_on_identical_sbml_inputs",
            "interpreter_rerun": False,
            "reference_path": "reference-1.6.3/" + archived_path.name,
            "reference_sha256": hashlib.sha256(archived_path.read_bytes()).hexdigest(),
            "all_exported_input_hashes_match": True,
            "reference_environment": archived["environment"],
        },
        "claim_boundary": (
            f"SynKit {synkit.__version__} outputs are recomputed and compared with "
            "recorded CRNT4SBML outputs for byte-identical exported SBML inputs. The "
            "legacy CRNT4SBML interpreter was not rerun. This checks shared "
            "structural quantities, not the SynKit-only analyses or relative speed."
        ),
        **comparison,
    }
    _write(out, "crn-external-crnt4sbml.json", report)
    return report


def structural_audit(out, library):
    """Cross-check complete C++ enumeration with the separate Python search."""
    from siphon_oracle import enumerate_models
    from Experiment.CRN.benchmarks import BENCHMARK_NETWORKS, build_kegg_crn
    from Experiment.CRN.benchmarks.scaling import generate_network
    from Experiment.CRN.studies.formose_case_study import expand
    from Experiment.CRN.tools.common import environment, source_hashes
    from synkit.CRN import enumerate_semiflows, enumerate_structural_sets, read_sbml
    from synkit.CRN.Petrinet.net import PetriNet
    from synkit.CRN.Props.stoich import build_S

    sys.path.insert(0, str(PAPER))
    from verify_revision import structural_witnesses

    cases = [("reference/" + f.name, f.build(), None, None) for f in BENCHMARK_NETWORKS]
    for row in _read(out, "crn-scaling.json")["measurements"]:
        if row["family"] == "random_sparse" and row["task"] == "minimal_siphons":
            cases.append((f"random/{row['size']}", generate_network("random_sparse", row["size"], seed=0),
                          row["result"], None))
    for row in _read(out, "crn-biomodels.json")["models"]:
        if row["n_reactions"]:
            model = read_sbml(ROOT / "Experiment/CRN/data/biomodels" / (row["model_id"] + ".xml"))
            cases.append(("biomodels/" + row["model_id"], model, row["n_minimal_siphons"],
                          row["n_conserved_moieties"]))
    for row in _read(out, "crn-kegg.json")["modules"]:
        cases.append(("kegg/" + row["module_id"], build_kegg_crn(row["module_id"], drop_currency=True),
                      len(row["siphons"]), len(row["moieties"])))
    cases.append(("formose", expand(4), 2, 1))
    records = []
    for name, model, expected_siphons, expected_pools in cases:
        net = PetriNet.from_syncrn(model)
        native = enumerate_structural_sets(model, library_path=library, max_nodes=None,
                                           max_results=None, time_limit=None)
        assert native.complete and native.termination_reason == "complete", name
        actual = {frozenset(g) for g in native.value["sets"]}
        oracle, states = enumerate_models(model)
        assert actual == {frozenset(str(s) for s in group) for group in oracle}, name
        assert expected_siphons is None or len(actual) == expected_siphons, name
        species_order, reaction_order, matrix = build_S(model)
        pools = enumerate_semiflows(matrix, kind="p", library_path=library)
        assert pools.complete, name
        flows = pools.value["flows"]
        assert expected_pools is None or len(flows) == expected_pools, name
        record = {
            "name": name, "siphons_complete": True, "p_semiflows_complete": True,
            "independent_siphon_outputs_agree": True,
            "independent_search_states": states,
            "all_siphons_valid_and_minimal": True, "all_p_witnesses_exact": True,
            "species": dict(net.place_labels),
            "species_source_ids": dict(net.place_source_node_ids),
            "reaction_source_ids": dict(net.transition_source_node_ids),
            "reactions": {rid: {"pre": t.pre, "post": t.post} for rid, t in net.transitions.items()},
            "minimal_siphons": sorted(sorted(g) for g in actual),
            "matrix_species_order": list(species_order),
            "matrix_reaction_order": list(reaction_order),
            "matrix_shape": list(matrix.shape),
            "sparse_stoichiometry": [{str(j): int(v) for j, v in enumerate(row) if v} for row in matrix],
            "minimal_p_semiflows": flows,
        }
        structural_witnesses(record)
        records.append(record)
        print(name, len(actual), "siphons;", len(flows), "pools: PASS", flush=True)
    assert len(records) == 31
    _write(out, "structural-verification.json", {
        "status": "PASS", "complete": True, "cases": records,
        "source_sha256": source_hashes(), "environment": environment(),
        "native_backend": environment()["crn_backend"],
        "independent_search_source": (SCRIPTS / "siphon_oracle.py").read_text(),
        "verification_source": Path(__file__).read_text(),
        "limits": {"structural_nodes": None, "structural_time": None, "structural_results": None},
        "claim_boundary": (
            "Complete native enumeration is cross-checked against the separate "
            "Python implication search on 31 inputs. Exact witnesses and all "
            "returned siphons are checked against incidence records. The searches "
            "share inputs and a strategy, not independent scientific models."
        ),
    })


def complete(out, native_dir, rdkit_dir):
    """Stage 2: external comparison, structural audit and glycolysis replay."""
    from Experiment.CRN.tools.common import source_hashes
    from glycolysis_analysis import build_analysis, render_tex
    from synkit.CRN import native_backend
    from synkit.CRN.Construct import build_native_expander, native_expander

    summary = _read(out, "summary.json")
    assert summary["status"] == "PASS", "run_all did not pass"
    library = Path(summary["environment"]["crn_backend"]["library_path"])
    initial_sources = summary["source_sha256"]
    assert initial_sources == source_hashes(), "Sources changed during the run"
    expander = build_native_expander(native_dir / "expander", rdkit_dir=rdkit_dir)
    with native_backend(library), native_expander(expander):
        external = external_comparison(out, library)
        structural_audit(out, library)
        glycolysis = build_analysis()
        _write(out, "glycolysis_analysis.json", glycolysis)
        (out / "glycolysis_analysis.tex").write_text(render_tex(glycolysis))
    report_path = out / "crn-external-crnt4sbml.json"
    summary["studies"] = [s for s in summary["studies"] if s["name"] != "external-crnt4sbml"]
    summary["studies"].append({
        "name": "external-crnt4sbml", "path": report_path.name,
        "schema": external["schema"], "status": "PASS", "checks": external["checks"],
        "sha256": hashlib.sha256(report_path.read_bytes()).hexdigest(), "error": None,
    })
    summary["omitted_studies"] = {}
    summary["all_studies_included"] = True
    summary["external_execution"] = external["external_execution"]
    summary["run_date"] = date.today().isoformat()
    assert initial_sources == source_hashes(), "Sources changed during the audit"
    _write(out, "summary.json", summary)
    _write(out, "source-before-run.json", initial_sources)
    _write(out, "native-build.json", json.loads(library.with_suffix(".json").read_text()))


def install(out):
    """Stage 3: copy reports into analysis/ under the manuscript's names."""
    target = PAPER / "analysis"
    summary = _read(out, "summary.json")
    for row in summary["studies"]:
        name = row["path"]
        installed = INSTALL_NAMES.get(name, name)
        shutil.copy2(out / name, target / installed)
        row["path"] = installed
    for name in ("structural-verification.json", "source-before-run.json",
                 "native-build.json", "glycolysis_analysis.json"):
        shutil.copy2(out / name, target / name)
    shutil.copytree(out / "reference-1.6.3", target / "reference-1.6.3", dirs_exist_ok=True)
    _write(target, "summary.json", summary)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--native-build-dir", type=Path)
    parser.add_argument("--rdkit-dir", type=Path)
    parser.add_argument("--skip-studies", action="store_true",
                        help="reuse stage-1 reports already in --work-dir")
    arguments = parser.parse_args()
    logging.disable(logging.INFO)
    from Experiment.CRN.tools.common import _default_rdkit_sdk

    out = arguments.work_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    native_dir = (arguments.native_build_dir or out / "native").resolve()
    rdkit_dir = arguments.rdkit_dir or _default_rdkit_sdk()
    if rdkit_dir is None:
        raise SystemExit("An RDKit C++ SDK is required; pass --rdkit-dir")
    if not arguments.skip_studies:
        run_studies(out, native_dir)
    complete(out, native_dir, rdkit_dir)
    install(out)
    print("Evidence regenerated and installed under", PAPER / "analysis")


if __name__ == "__main__":
    main()
