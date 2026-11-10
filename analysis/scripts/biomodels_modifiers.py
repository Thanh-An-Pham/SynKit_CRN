#!/usr/bin/env python3
"""Deficiency of the BioModels entries under three readings of catalysts.

* kept      : species listed on both sides of a reaction are kept (SynKit's reader)
* cancelled : a species with equal stoichiometry on both sides is removed, as in
              net-stoichiometry tools
* modifiers : every SBML modifier is added to both sides with coefficient 1

Reactions marked reversible="true" are written out as two directed reactions
before import, so the result does not depend on how a SynKit version expands them.
Input files are checked against the SHA-256 hashes in crn-biomodels.json.

Usage:  python analysis/scripts/biomodels_modifiers.py <dir with BIOMD*.xml>
Writes: analysis/crn-biomodels-modifiers.json
"""
import copy
import hashlib
import importlib.metadata as md
import json
import platform
import sys
import xml.etree.ElementTree as ET
from datetime import date
from pathlib import Path

import synkit.CRN as C

ANALYSIS = Path(__file__).resolve().parents[1]


def _loc(e):
    return e.tag.split("}")[-1]


def variant(root, mode):
    root = copy.deepcopy(root)
    ns = root.tag.split("}")[0] + "}" if root.tag.startswith("{") else ""
    model = next(e for e in root.iter() if _loc(e) == "model")
    lst = next((e for e in model if _loc(e) == "listOfReactions"), None)
    changed = 0
    for r in list(lst) if lst is not None else []:
        sides = {_loc(e): e for e in r if _loc(e) in ("listOfReactants", "listOfProducts", "listOfModifiers")}
        reac, prod = sides.get("listOfReactants"), sides.get("listOfProducts")
        if mode == "cancelled" and reac is not None and prod is not None:
            coef = lambda e: e.get("stoichiometry", "1")
            left = {e.get("species"): e for e in reac}
            for p in list(prod):
                q = left.get(p.get("species"))
                if q is not None and coef(q) == coef(p):
                    reac.remove(q)
                    prod.remove(p)
                    changed += 1
        if mode == "modifiers" and sides.get("listOfModifiers") is not None:
            for m in sides["listOfModifiers"]:
                for tag in ("listOfReactants", "listOfProducts"):
                    side = sides.get(tag)
                    if side is None:
                        side = sides[tag] = ET.SubElement(r, ns + tag)
                    side.append(ET.Element(ns + "speciesReference", {"species": m.get("species"), "stoichiometry": "1"}))
                changed += 1
        if r.get("reversible") == "true":
            r.set("reversible", "false")
            back = copy.deepcopy(r)
            back.set("id", r.get("id") + "_back")
            for side in back:
                if _loc(side) == "listOfReactants":
                    side.tag = side.tag.replace("listOfReactants", "listOfProducts")
                elif _loc(side) == "listOfProducts":
                    side.tag = side.tag.replace("listOfProducts", "listOfReactants")
            lst.append(back)
    return C.crn_from_sbml(root), changed


def main(folder):
    folder = Path(folder)
    evidence = json.loads((ANALYSIS / "crn-biomodels.json").read_text())
    models = []
    for entry in evidence["models"]:
        path = folder / f"{entry['model_id']}.xml"
        data = path.read_bytes()
        if hashlib.sha256(data).hexdigest() != entry["sha256"]:
            raise SystemExit(f"{path.name}: SHA-256 differs from crn-biomodels.json")
        root = ET.fromstring(data)
        n_modifier_refs = sum(_loc(e) == "modifierSpeciesReference" for e in root.iter())
        row = {"model_id": entry["model_id"], "n_modifier_refs": n_modifier_refs}
        for mode in ("kept", "cancelled", "modifiers"):
            crn, changed = variant(root, mode)
            if crn.n_reactions == 0:
                row = {"model_id": entry["model_id"], "n_modifier_refs": n_modifier_refs, "reaction_free": True}
                break
            row[f"deficiency_{mode}"] = C.crnt_summary(crn).deficiency
            if mode == "cancelled":
                row["n_cancelled_participants"] = changed
        models.append(row)
    out = {
        "schema": "crn-biomodels-modifiers/1",
        "date": date.today().isoformat(),
        "source": "analysis/scripts/biomodels_modifiers.py",
        "inputs": "BioModels files checked against the SHA-256 hashes in crn-biomodels.json",
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "synkit": md.version("synkit"),
            "numpy": md.version("numpy"),
            "networkx": md.version("networkx"),
        },
        "models": models,
    }
    (ANALYSIS / "crn-biomodels-modifiers.json").write_text(json.dumps(out, indent=2) + "\n")
    for row in models:
        print(row)


if __name__ == "__main__":
    main(sys.argv[1])
