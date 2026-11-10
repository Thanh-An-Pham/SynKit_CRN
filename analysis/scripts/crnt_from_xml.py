"""SynKit-free CRNT numbers from an SBML file: complexes, linkage classes, exact rank, deficiency.

Only the standard library is used. A reaction marked reversible="true" contributes its
reverse; an omitted attribute is read as irreversible (the paper's convention).

Usage:  python3 -I analysis/scripts/crnt_from_xml.py <model.xml>
"""
import sys
import xml.etree.ElementTree as ET
from fractions import Fraction


def loc(e):
    return e.tag.split("}")[-1]


def side(r, name):
    out = {}
    for lst in (c for c in r if loc(c) == name):
        for ref in lst:
            if loc(ref) == "speciesReference":
                out[ref.get("species")] = out.get(ref.get("species"), 0) + Fraction(ref.get("stoichiometry", "1"))
    return out


def main(path):
    root = ET.parse(path).getroot()
    reactions = []
    for r in (e for e in root.iter() if loc(e) == "reaction"):
        lhs, rhs = side(r, "listOfReactants"), side(r, "listOfProducts")
        reactions.append((lhs, rhs))
        if r.get("reversible") == "true":
            reactions.append((rhs, lhs))
    key = lambda d: tuple(sorted((s, c) for s, c in d.items() if c))
    complexes = sorted({key(l) for l, _ in reactions} | {key(p) for _, p in reactions})
    idx = {c: i for i, c in enumerate(complexes)}
    parent = list(range(len(complexes)))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    for l, p in reactions:
        parent[find(idx[key(l)])] = find(idx[key(p)])
    linkage = len({find(i) for i in range(len(complexes))})
    species = sorted({s for l, p in reactions for s in list(l) + list(p)})
    rows = [[p.get(s, 0) - l.get(s, 0) for l, p in reactions] for s in species]
    rank = 0
    for c in range(len(reactions)):
        piv = next((i for i in range(rank, len(rows)) if rows[i][c] != 0), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for i in range(len(rows)):
            if i != rank and rows[i][c] != 0:
                f = rows[i][c] / rows[rank][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[rank])]
        rank += 1
    print(f"species {len(species)}  reactions {len(reactions)}  complexes {len(complexes)}  "
          f"linkage classes {linkage}  rank {rank}  deficiency {len(complexes) - linkage - rank}  "
          f"conservation laws {len(species) - rank}")


if __name__ == "__main__":
    main(sys.argv[1])
