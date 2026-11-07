#!/usr/bin/env python3
"""Check revision evidence, printed equations, and generated numerical assets."""

import os
from pathlib import Path
import hashlib
import json
import re
import sys
import tempfile
from fractions import Fraction
from itertools import product
from math import gcd
from functools import reduce

PAPER = Path(__file__).resolve().parent
ROOT = Path(os.environ.get("SYNKIT_ROOT", PAPER.parents[1] / "SynEco/SynKit")).resolve()


def complete_decisions(value):
    if isinstance(value, list):
        for item in value:
            complete_decisions(item)
    elif isinstance(value, dict):
        if "termination_reason" in value:
            assert value["termination_reason"] in {
                "complete",
                "certificate",
                "exhausted",
            }, value
        if "decision_complete" in value:
            assert value["decision_complete"] is True, value
            if value.get("realizable") is True:
                assert value["independent_replay_ok"] is True, value
            elif value.get("realizable") is False:
                assert value["exhaustive"] is True, value
        assert value.get("status") != "inconclusive", value
        for item in value.values():
            complete_decisions(item)


def exact_feasibility(certificate):
    assert certificate["exact"] is True
    rows, columns = certificate["shape"]
    matrix = certificate["problem_matrix"]
    assert len(matrix) == rows and all(len(row) == columns for row in matrix)
    assert all(type(x) is int for row in matrix for x in row)
    if certificate["holds"] is True:
        witness = list(map(Fraction, certificate["primal_witness"]))
        assert len(witness) == columns and all(x > 0 for x in witness)
        assert all(sum(a * x for a, x in zip(row, witness)) == 0 for row in matrix)
    else:
        assert certificate["holds"] is False
        witness = list(map(Fraction, certificate["dual_witness"]))
        assert len(witness) == rows
        image = [
            sum(matrix[i][j] * witness[i] for i in range(rows)) for j in range(columns)
        ]
        assert all(x >= 0 for x in image) and any(x > 0 for x in image)
        assert image == list(map(Fraction, certificate["dual_image"]))


def structural_witnesses(record):
    assert record["siphons_complete"] is True and record["p_semiflows_complete"] is True
    assert record["independent_siphon_outputs_agree"] is True
    species = set(record["species"])
    reactions = record["reactions"]
    transitions = [(set(r["post"]), set(r["pre"])) for r in reactions.values()]
    groups = [frozenset(s) for s in record["minimal_siphons"]]
    assert len(set(groups)) == len(groups)
    for group in groups:
        assert group and group <= species
        assert all(
            not group.intersection(t) or group.intersection(r) for t, r in transitions
        )
        for place in group:
            child = set(group) - {place}
            while True:
                before = set(child)
                for trigger, required in transitions:
                    if child.intersection(trigger) and not child.intersection(required):
                        child.difference_update(trigger)
                if child == before:
                    break
            assert not child, (record["name"], "non-minimal siphon")
    rows, columns = record["matrix_shape"]
    matrix = [
        {int(j): value for j, value in row.items()}
        for row in record["sparse_stoichiometry"]
    ]
    assert len(matrix) == rows
    assert all(
        0 <= j < columns and type(value) is int
        for row in matrix
        for j, value in row.items()
    )
    row_index = {str(s): i for i, s in enumerate(record["matrix_species_order"])}
    column_index = {str(r): i for i, r in enumerate(record["matrix_reaction_order"])}
    expected = [{} for _ in range(rows)]
    for rid, reaction in reactions.items():
        j = column_index[str(record["reaction_source_ids"][rid])]
        for side, sign in (("post", 1), ("pre", -1)):
            for sid, value in reaction[side].items():
                i = row_index[str(record["species_source_ids"][sid])]
                expected[i][j] = expected[i].get(j, 0) + sign * value
    expected = [{j: value for j, value in row.items() if value} for row in expected]
    assert matrix == expected, (record["name"], "matrix differs from reaction records")
    for witness in record["minimal_p_semiflows"]:
        assert len(witness) == rows and all(type(x) is int and x >= 0 for x in witness)
        assert any(witness) and reduce(gcd, witness, 0) == 1
        image = [0] * columns
        for i, row in enumerate(matrix):
            for j, value in row.items():
                image[j] += value * witness[i]
        assert not any(image), (record["name"], "invalid exact P-semiflow")


def main():
    source_manifest = json.loads(
        (PAPER / "analysis/source-before-run.json").read_text()
    )
    for filename, digest in source_manifest.items():
        assert (
            hashlib.sha256((ROOT / filename).read_bytes()).hexdigest() == digest
        ), filename
    sys.path.insert(0, str(ROOT))
    from Experiment.CRN.tools.common import source_hashes
    from Experiment.CRN.tools.check_reports import check_reports
    import synkit

    assert (
        source_manifest == source_hashes()
    ), "The evidence does not cover the current CRN source tree"
    summary = json.loads((PAPER / "analysis/summary.json").read_text())
    assert summary["status"] == "PASS" and summary["complete"] is True
    assert summary["schema"] == "synkit.crn-experiments/2"
    assert not check_reports(PAPER / "analysis/summary.json", require_all=True)
    assert summary["source_sha256"] == source_manifest
    assert synkit.__version__ == summary["environment"]["synkit"] == "1.6.5"
    assert len(summary["studies"]) == 10
    build = json.loads((PAPER / "analysis/native-build.json").read_text())
    for name in (
        "validation",
        "scaling",
        "kegg",
        "formose",
        "biomodels",
        "external-crnt4sbml",
        "review-experiments",
        "supplement-audit",
        "si-donor-control",
        "hydrogen",
    ):
        report = json.loads(
            (PAPER / "analysis" / ("crn-" + name + ".json")).read_text()
        )
        assert report["status"] == "PASS" and all(report["checks"].values()), name
        assert report["environment"]["synkit"] == "1.6.5", name
        assert (
            report["environment"]["git_commit"]
            == summary["environment"]["git_commit"]
        ), name
        complete_decisions(report)
        backend = report["environment"]["crn_backend"]
        assert backend["production_backend"] == "cpp", name
        assert backend["build_sha256"] == build["build_sha256"], name
        assert backend["library_sha256"] == build["library_sha256"], name
    external = json.loads((PAPER / "analysis/crn-external-crnt4sbml.json").read_text())
    execution = external["external_execution"]
    assert execution["mode"] == "recorded_outputs_on_identical_sbml_inputs"
    assert execution["interpreter_rerun"] is False
    assert execution["all_exported_input_hashes_match"] is True
    reference_path = PAPER / "analysis" / execution["reference_path"]
    assert hashlib.sha256(reference_path.read_bytes()).hexdigest() == execution["reference_sha256"]
    reference = json.loads(reference_path.read_text())
    old_rows = {row["name"]: row for row in reference["networks"]}
    current_rows = {row["name"]: row for row in external["networks"]}
    assert len(current_rows) == len(old_rows) == 12
    assert current_rows.keys() == old_rows.keys()
    from synkit.CRN import write_sbml
    from Experiment.CRN.benchmarks import BENCHMARK_NETWORKS

    with tempfile.TemporaryDirectory(prefix="synkit-crn-input-check-") as folder:
        for fixture in BENCHMARK_NETWORKS:
            row = current_rows[fixture.name]
            old = old_rows[fixture.name]
            path = Path(folder) / (fixture.name + ".xml")
            write_sbml(fixture.build(), path, model_id=fixture.name, model_name=fixture.description)
            assert hashlib.sha256(path.read_bytes()).hexdigest() == row["input_sha256"] == old["input_sha256"]
            assert row["crnt4sbml"] == old["crnt4sbml"] == row["synkit"]
    scaling = json.loads((PAPER / "analysis/crn-scaling.json").read_text())
    parameters = scaling["parameters"]
    expected_grid = set(
        product(parameters["families"], parameters["sizes"], parameters["tasks"])
    )
    rows = scaling["measurements"]
    actual_grid = [(r["family"], r["size"], r["task"]) for r in rows]
    assert (
        scaling["complete"] is True and len(actual_grid) == len(set(actual_grid)) == 60
    )
    assert set(actual_grid) == expected_grid
    assert (
        not scaling["errors"]
        and not scaling["timeouts"]
        and not scaling["resource_limits"]
    )
    for row in rows:
        assert row["backend"] == "cpp" and row["complete"] is True
        assert row["termination_reason"] == "complete" and row["result"] is not None
        assert row["seconds"] is not None and row["seconds"] >= 0
        assert row["error"] is None and not row["timed_out"]
    random_rows = sorted(
        (
            r
            for r in rows
            if r["family"] == "random_sparse" and r["task"] == "minimal_siphons"
        ),
        key=lambda r: r["size"],
    )
    assert [(r["n_species"], r["result"]) for r in random_rows] == [
        (49, 331),
        (99, 272),
        (197, 47),
        (395, 87),
    ]
    biomodels = json.loads((PAPER / "analysis/crn-biomodels.json").read_text())
    reaction_models = [r for r in biomodels["models"] if r["n_reactions"] > 0]
    assert len(reaction_models) == 10
    for model in reaction_models:
        assert model["roundtrip_ok"] is True
        assert isinstance(model["persistent"], bool)
        for key in ("n_conserved_moieties", "n_minimal_siphons", "deficiency"):
            assert isinstance(model[key], int) and model[key] >= 0
        for key in ("siphon_seconds", "moiety_seconds", "persistence_seconds"):
            assert model[key] is not None and model[key] >= 0
        assert not any(key.endswith("_error") for key in model)
        assert (
            hashlib.sha256(
                (
                    ROOT
                    / "Experiment/CRN/data/biomodels"
                    / (model["model_id"] + ".xml")
                ).read_bytes()
            ).hexdigest()
            == model["sha256"]
        )
    model19 = next(m for m in reaction_models if m["model_id"] == "BIOMD0000000019")
    assert model19["n_minimal_siphons"] == 20 and model19["n_conserved_moieties"] == 19
    kegg = json.loads((PAPER / "analysis/crn-kegg.json").read_text())
    assert len(kegg["modules"]) == 4
    assert all(
        m["siphons_complete"] is True
        and m["siphon_termination_reason"] == "complete"
        and isinstance(m["persistent"], bool)
        for m in kegg["modules"]
    )

    structures = json.loads(
        (PAPER / "analysis/structural-verification.json").read_text()
    )
    assert structures["status"] == "PASS" and structures["complete"] is True
    assert structures["source_sha256"] == source_manifest
    assert all(value is None for value in structures["limits"].values())
    assert structures["native_backend"]["build_sha256"] == build["build_sha256"]
    records = {row["name"]: row for row in structures["cases"]}
    assert len(records) == len(structures["cases"]) == 31
    for record in records.values():
        structural_witnesses(record)
    for row in random_rows:
        assert len(records[f"random/{row['size']}"]["minimal_siphons"]) == row["result"]
    for model in reaction_models:
        record = records["biomodels/" + model["model_id"]]
        assert len(record["minimal_siphons"]) == model["n_minimal_siphons"]
        assert len(record["minimal_p_semiflows"]) == model["n_conserved_moieties"]
    for model in kegg["modules"]:
        record = records["kegg/" + model["module_id"]]
        assert len(record["minimal_siphons"]) == len(model["siphons"])
        assert len(record["minimal_p_semiflows"]) == len(model["moieties"])

    formose = json.loads((PAPER / "analysis/crn-formose.json").read_text())
    assert all(
        row["rule_backend"] == "native" and row["hydrogen_model"] == "lewis"
        for row in formose["growth"]
    ), "Formose growth must use native expansion with hydrogen counts"
    hydrogen = json.loads((PAPER / "analysis/crn-hydrogen.json").read_text())
    assert hydrogen["schema"] == "synkit.crn-hydrogen/1"
    assert [row["depth"] for row in hydrogen["formose"]] == [4, 5, 6]
    assert all(row["identical"] for row in hydrogen["formose"])
    for case in hydrogen["balance_cases"]:
        lewis = case["lewis"]
        assert lewis["n_reactions"] > 0 and lewis["n_balanced"] == lewis["n_reactions"], case["name"]
        assert all(not r["imbalance"] for r in lewis["reactions"]), case["name"]
    expander = hydrogen["environment"]["crn_expander"]
    assert expander and expander["build_sha256"], "Native expander build is not recorded"
    assert (
        f"There are {formose['structure']['n_minimal_traps']} minimal traps"
        in (PAPER / "SI/si.tex").read_text()
    )
    for filename, schema in [
        ("crn-review-experiments.json", "synkit.crn-representation/1"),
        ("crn-supplement-audit.json", "synkit.crn-supplementary/1"),
        ("crn-si-donor-control.json", "synkit.crn-glycolysis-donors/1"),
    ]:
        report = json.loads((PAPER / "analysis" / filename).read_text())
        assert report["schema"] == schema, filename
        assert report["status"] == "PASS" and all(report["checks"].values()), filename
    review = json.loads((PAPER / "analysis/crn-review-experiments.json").read_text())
    for filename, digest in review["source_sha256"].items():
        assert (
            hashlib.sha256((ROOT / filename).read_bytes()).hexdigest() == digest
        ), f"Stale experiment source: {filename}"
    supplement = json.loads((PAPER / "analysis/crn-supplement-audit.json").read_text())

    for name, decisions in supplement["feasibility_certificates"].items():
        for decision, certificate in decisions.items():
            exact_feasibility(certificate)
            assert certificate["holds"] is supplement["feasibility"][name][decision]
        left = decisions["conservative"]["problem_matrix"]
        right = decisions["consistent"]["problem_matrix"]
        assert left == [list(column) for column in zip(*right)]
    from synkit.CRN import SynCRN
    from synkit.CRN.Props.stoich import stoichiometric_matrix

    glycolysis = SynCRN.from_reaction_strings(
        [
            supplement["equations"][key]
            for key in sorted(supplement["equations"], key=lambda key: int(key[1:]))
        ]
    )
    assert (
        stoichiometric_matrix(glycolysis).tolist()
        == supplement["feasibility_certificates"]["glycolysis"]["consistent"][
            "problem_matrix"
        ]
    )
    assert all(row["complete"] is True for row in supplement["symmetry"].values())
    assert all(
        row["exact_identity_holds"] is True
        for row in supplement["numerical_comparison"]
    )

    donors = json.loads((PAPER / "analysis/crn-si-donor-control.json").read_text())
    for filename, digest in donors["source_sha256"].items():
        assert (
            hashlib.sha256((ROOT / filename).read_bytes()).hexdigest() == digest
        ), filename
    assert donors["equations"] == supplement["equations"]
    assert (
        hashlib.sha256(
            (ROOT / "Experiment/CRN/studies/supplementary.py").read_bytes()
        ).hexdigest()
        == supplement["script_sha256"]
    )
    printed = dict(
        (f"r{n}", eq.replace(r"\(\rightarrow\)", ">>"))
        for n, eq in re.findall(
            r"^r(\d+) & (.*?) \\\\", (PAPER / "SI/si.tex").read_text(), re.M
        )
    )
    def reaction_sides(equations):
        result = {}
        for key, equation in equations.items():
            model = SynCRN.from_reaction_strings([equation])
            reaction, = model.reactions.values()
            result[key] = tuple(
                {model.species[sid].label: count for sid, count in side.items()}
                for side in (reaction.lhs, reaction.rhs)
            )
        return result

    assert (
        reaction_sides(printed) == reaction_sides(supplement["equations"])
    ), "SI reaction table differs from audited stoichiometric multisets"
    sys.path.insert(0, str(PAPER / "analysis/scripts"))
    from sync_manuscript import generated_blocks, block_pattern, marked_block

    for filename, blocks in generated_blocks(PAPER).items():
        source = (PAPER / filename).read_text()
        for name, body in blocks.items():
            matches = list(block_pattern(name).finditer(source))
            assert len(matches) == 1, (filename, name)
            assert matches[0].group() == marked_block(name, body), (
                "Stale inline numerical values", filename, name
            )
    for path in (PAPER / "Figure").iterdir():
        assert path.suffix == ".pdf" or (
            path.suffix == ".tex" and r"\begin{tikzpicture}" in path.read_text()
        ), f"Non-figure file in Figure/: {path.name}"
    glycolysis_figure = json.loads(
        (PAPER / "analysis/glycolysis_analysis.json").read_text()
    )
    backend = glycolysis_figure["environment"]["crn_backend"]
    assert glycolysis_figure["environment"]["synkit"] == "1.6.5"
    assert (
        backend["production_backend"] == "cpp"
        and backend["build_sha256"] == build["build_sha256"]
    )
    replay = glycolysis_figure["realizability"]
    assert (
        replay["status"] == "realizable"
        and replay["termination_reason"] == "certificate"
    )
    assert replay["certificate_length"] == len(replay["certificate"]) == 15
    assert replay["explored_states"] == 30
    assert all(
        row["nonnegative"] is True and row["backbone_invariant"] == 2
        for row in replay["trajectory"]
    )
    assert replay["trajectory"][-1]["marking"] == {"Pyruvate": 2}
    text = (PAPER / "main.tex").read_text()
    abstract = text.split("\\abstract{", 1)[1].split("\\keywords", 1)[0]
    assert len(abstract.split()) <= 350, "Abstract exceeds 350 words"
    assert "\\cite" not in abstract, "Abstract contains citations"
    assert "version~1.6.5 source checkout" in text
    assert summary["environment"]["git_commit"][:12] in text
    print(
        "Revision evidence, source hashes, SI equations, generated tables, and abstract: PASS"
    )


if __name__ == "__main__":
    main()
