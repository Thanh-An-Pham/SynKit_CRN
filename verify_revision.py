#!/usr/bin/env python3
"""Check revision evidence and manuscript consistency with saved numerical assets."""

import argparse
from collections import Counter
import os
from pathlib import Path
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


def verify_document_content():
    """Check manuscript values against saved reports without running experiments."""
    sys.path.insert(0, str(PAPER / "analysis/scripts"))
    from tables import table_blocks

    for filename, blocks in table_blocks(PAPER).items():
        source = (PAPER / filename).read_text()
        for name, body in blocks.items():
            rows = '\n'.join(line for line in body.splitlines() if not line.lstrip().startswith('%')).strip()
            assert rows in source, ('Table differs from saved results', filename, name)

    main_text = (PAPER / "main.tex").read_text()
    si_text = (PAPER / "SI/si.tex").read_text()
    sections = re.findall(r"\\section\{([^}]+)\}", si_text)
    assert sections == [
        "Python interface example",
        "Reference-network validation and numerical comparison",
        "Computational cost across network families",
        "BioModels: structural SBML interoperability",
        "KEGG: module structure and glycolysis controls",
        "Formose networks: expansion and precursor requirements",
    ]
    for obsolete in (
        "tab:glycolysis_full_reactions", "tab:glycolysis_species_map",
        "fig:glycolysis_donors",
        "fig:structural_landscape", "Automorphisms", "Jacobian",
    ):
        assert obsolete not in si_text, obsolete
    assert r"\input{../Figure/arch.tex}" not in si_text
    for source, target, prefix in (
        (main_text, si_text, "supp:"), (si_text, main_text, "main:")
    ):
        labels = set(re.findall(r"\\label\{([^}]+)\}", target))
        labels.update(re.findall(r"\blabel=\{([^}]+)\}", target))
        references = re.findall(r"\\(?:ref|eqref|pageref)\*?\{([^}]+)\}", source)
        for label in references:
            if label.startswith(prefix):
                assert label[len(prefix):] in labels, ("Unresolved reference", label)
    bibliography = set(re.findall(r"@\w+\{([^,]+),", (PAPER / "refs.bib").read_text()))
    for filename, source in (("main.tex", main_text), ("SI/si.tex", si_text)):
        local_labels = set(re.findall(r"\\label\{([^}]+)\}", source))
        local_labels.update(re.findall(r"\blabel=\{([^}]+)\}", source))
        for label in re.findall(r"\\(?:ref|eqref|pageref)\*?\{([^}]+)\}", source):
            if not label.startswith(("main:", "supp:", "#")):
                assert label in local_labels, ("Unresolved local reference", filename, label)
        for group in re.findall(r"\\cite\w*\*?(?:\[[^]]*\])*\{([^}]+)\}", source):
            for key in group.split(","):
                assert key.strip() in bibliography, ("Missing bibliography entry", filename, key)

    def read(name):
        return json.loads((PAPER / "analysis" / name).read_text())

    biomodels = read("crn-biomodels.json")
    assert len(biomodels["models"]) == 11
    reaction_models = [row for row in biomodels["models"] if row["n_reactions"]]
    assert len(reaction_models) == 10 and all(row["roundtrip_ok"] for row in reaction_models)
    kegg = read("crn-kegg.json")
    assert len(kegg["modules"]) == 4
    flux = kegg["glycolysis_flux"]
    assert len(flux["certificate"]) == sum(flux["flux"].values()) == 15
    assert Counter(flux["certificate"]) == flux["flux"]
    replay = read("glycolysis_analysis.json")["realizability"]
    assert replay["certificate"] == flux["certificate"]
    assert replay["explored_states"] == 30
    assert all(row["nonnegative"] and row["backbone_invariant"] == 2
               for row in replay["trajectory"])
    assert replay["trajectory"][-1]["marking"] == {"Pyruvate": 2}

    review = read("crn-review-experiments.json")
    assert review["status"] == "PASS" and all(review["checks"].values())
    complete_decisions(review["representation"])
    representation = review["representation"]
    assert representation["net_matrix_unchanged"] is True
    assert representation["species_order_unchanged"] is True
    assert representation["reaction_order_unchanged"] is True
    for arm in representation["currency_arms"]:
        outcome = arm["outcome"]
        assert arm["n_reactions"] == 15
        assert arm["n_species"] == (14 if arm["drop_currency"] else 24)
        assert outcome["initial_marking"] == {"alpha-D-Glucose": 1}
        assert outcome["requested_flux"] == flux["flux"]
        assert outcome["realizable"] is arm["drop_currency"]
    formose = read("crn-formose.json")
    assert formose["status"] == "PASS" and all(formose["checks"].values())
    structure = formose["structure"]
    assert (structure["n_species"], structure["n_reactions"], structure["rank"]) == (37, 44, 36)
    carbon = formose["carbon_conservation"]
    assert len(carbon["entries"]) == 37 and carbon["coefficient_equals_carbon_count"]
    assert all(row["carbons"] == row["coefficient"] for row in carbon["entries"])
    supplement = read("crn-supplement-audit.json")
    assert supplement["status"] == "PASS" and all(supplement["checks"].values())
    audit_table = si_text.split(r"\label{tab:reference_audit}", 1)[1].split(r"\end{table}", 1)[0]
    numerical_rows = re.findall(r"^([^\n]+?) & ([0-9.]+e[+-][0-9]+) & ([0-9]+) & ([0-9]+) \\", audit_table, re.M)
    assert len(numerical_rows) == len(supplement["numerical_comparison"]) == 12
    for actual, expected in zip(numerical_rows, supplement["numerical_comparison"]):
        assert actual[1:] == (f"{expected['svd_max_abs_residual']:.1e}",
                              str(expected["rounded_invalid_vectors"]), "0"), (
            "Stale SI numerical-audit row", expected["network"]
        )
    scaling = read("crn-scaling.json")
    plotted_tasks = {"minimal_siphons", "crnt_summary", "canonical_form", "conserved_moieties"}
    assert len(scaling["measurements"]) == 60
    assert sum(row["task"] in plotted_tasks for row in scaling["measurements"]) == 48
    assert sum(row["task"] == "rank" for row in scaling["measurements"]) == 12
    assert len(supplement["formose_routes"]) == 10
    complete_decisions(supplement["formose_routes"])
    complete_decisions(supplement["qualitative_routes"])
    for row in supplement["formose_routes"]:
        outcome = row["outcome"]
        if outcome["realizable"]:
            assert Counter(outcome["certificate"]) == outcome["requested_flux"]
    abstract = main_text.split(r"\begin{abstract}", 1)[1].split(r"\end{abstract}", 1)[0]
    assert len(abstract.split()) <= 350 and "\\cite" not in abstract
    if "glycolysis_inventory_controls" in review:
        inventory = review["glycolysis_inventory_controls"]
        complete_decisions(inventory)
        assert inventory["adenylate_siphon"]["valid_on_prescribed_transitions"] is True
        cases = {row["name"]: row for row in inventory["cases"]}
        assert cases["one_atp_with_adp_redox_and_phosphate"]["net_inventory_feasible"] is True
        assert cases["one_atp_with_adp_redox_and_phosphate"]["outcome"]["status"] == "unrealizable"
        assert cases["two_atp_with_adp_redox_and_phosphate"]["outcome"]["final_marking"]["Pyruvate"] == 2


def main():
    sys.path.insert(0, str(ROOT))
    from Experiment.CRN.tools.check_reports import check_reports
    import synkit

    summary = json.loads((PAPER / "analysis/summary.json").read_text())
    assert summary["status"] == "PASS" and summary["complete"] is True
    assert summary["schema"] == "synkit.crn-experiments/2"
    assert not check_reports(PAPER / "analysis/summary.json", require_all=True)
    assert synkit.__version__ == summary["environment"]["synkit"] == "1.7.0"
    assert len(summary["studies"]) == 10
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
        complete_decisions(report)
    external = json.loads((PAPER / "analysis/crn-external-crnt4sbml.json").read_text())
    execution = external["external_execution"]
    assert execution['mode'] == 'live_isolated_crnt4sbml'
    assert execution['interpreter_rerun'] is True
    current_rows = {row["name"]: row for row in external["networks"]}
    assert len(current_rows) == 12
    from synkit.CRN import crnt_summary
    from synkit.CRN.backend import python_backend
    from Experiment.CRN.studies import BENCHMARK_NETWORKS
    for fixture in BENCHMARK_NETWORKS:
        row = current_rows[fixture.name]
        assert row['reactions'] == list(fixture.reactions)
        assert row['crnt4sbml'] == row['synkit']
        with python_backend():
            computed = crnt_summary(fixture.build()).to_dict()
        computed['weakly_reversible'] = computed.pop('is_weakly_reversible')
        computed = json.loads(json.dumps(computed))
        assert {key: computed[key] for key in row['synkit']} == row['synkit']
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
    model19 = next(m for m in reaction_models if m["model_id"] == "BIOMD0000000019")
    assert model19["n_conserved_moieties"] == 19
    participants = json.loads((PAPER / "analysis/crn-biomodels-modifiers.json").read_text())
    assert participants["status"] == "PASS" and all(participants["checks"].values())
    for row in participants["models"]:
        assert row["agrees"] and row["net_stoichiometry_unchanged"]
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
    assert all(value is None for value in structures["limits"].values())
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
        row["rule_backend"] == "synreactor-native"
        and row["engine"] == "synkit_dpo_cpp"
        and row["controller"] == "cpp"
        and row["implicit_h"] is True
        for row in formose["growth"]
    ), "Formose growth must execute the genuine full C++ SynReactor DPO path"
    hydrogen = json.loads((PAPER / "analysis/crn-hydrogen.json").read_text())
    assert hydrogen["schema"] == "synkit.crn-hydrogen/1"
    assert [row["depth"] for row in hydrogen["formose"]] == [4, 5, 6]
    assert all(row["identical"] for row in hydrogen["formose"])
    for case in hydrogen["balance_cases"]:
        lewis = case["lewis"]
        assert lewis["n_reactions"] > 0 and lewis["n_balanced"] == lewis["n_reactions"], case["name"]
        assert all(not r["imbalance"] for r in lewis["reactions"]), case["name"]
    for filename, schema in [
        ("crn-review-experiments.json", "synkit.crn-representation/1"),
        ("crn-supplement-audit.json", "synkit.crn-supplementary/1"),
        ("crn-si-donor-control.json", "synkit.crn-glycolysis-donors/1"),
    ]:
        report = json.loads((PAPER / "analysis" / filename).read_text())
        assert report["schema"] == schema, filename
        assert report["status"] == "PASS" and all(report["checks"].values()), filename
    review = json.loads((PAPER / "analysis/crn-review-experiments.json").read_text())
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
    assert donors["equations"] == supplement["equations"]
    verify_document_content()
    sources = [PAPER / "main.tex", PAPER / "SI/si.tex", *(PAPER / "Figure").glob("*.tex")]
    used_pdfs = set()
    for path in sources:
        for name in re.findall(r"\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}", path.read_text()):
            for candidate in (path.parent / name, PAPER / name):
                if candidate.is_file():
                    used_pdfs.add(candidate.resolve())
                    break
    for path in (PAPER / "Figure").iterdir():
        if path.suffix == ".tex":
            source = path.read_text()
            assert r"\begin{tikzpicture}" in source, ("Not a TikZ fragment", path.name)
            assert r"\documentclass" not in source and r"\begin{document}" not in source, (
                "TikZ figures must be loaded inline", path.name
            )
            assert not path.with_suffix(".pdf").exists(), ("Compiled TikZ duplicate", path.name)
        else:
            assert path.suffix == ".pdf" and path.resolve() in used_pdfs, (
                "Unused or generated standalone asset in Figure", path.name
            )
    glycolysis_figure = json.loads(
        (PAPER / "analysis/glycolysis_analysis.json").read_text()
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
    abstract = text.split(r"\begin{abstract}", 1)[1].split(r"\end{abstract}", 1)[0]
    assert len(abstract.split()) <= 350, "Abstract exceeds 350 words"
    assert "\\cite" not in abstract, "Abstract contains citations"
    from formose_dpo_assets import verify as verify_dpo
    dpo = verify_dpo(render=False)
    assert len(dpo["completed_networks"]) == 18 and dpo["verified_expansions"] == 414
    route_comparison = json.loads((PAPER / "analysis/formose-expansion/route-comparison.json").read_text())
    assert route_comparison["all_ten_agree"] is True
    assert route_comparison["positive_certificates_replayed"] == 3
    route_rows = {(row['pool'], row['target']): row for row in route_comparison['rows']}
    for row in supplement['formose_routes']:
        compared = route_rows[(row['pool'], row['target'])]
        assert compared['initial_marking'] == row['outcome']['initial_marking']
        assert compared['requested_flux'] == row['outcome']['requested_flux']
        assert compared['synkit_realizable'] == compared['molfoundry_realizable'] == row['outcome']['realizable']
    print(
        "Scientific results, exact witnesses, SI tables, and abstract: PASS"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--documents-only", action="store_true",
        help="Check manuscripts against saved reports; omit live source and experiment checks.",
    )
    args = parser.parse_args()
    if args.documents_only:
        verify_document_content()
        print("Manuscript structure, cross-references, tables, and saved results: PASS")
    else:
        main()
