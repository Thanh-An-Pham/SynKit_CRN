SynKit 1.7.0 repair and reproduction plan
======================================

Target: prepare and validate a local 1.7.0 source snapshot. The user will push
and publish after the checks below pass. A local version number is not evidence
that a public release or DOI exists.

1. Correct the supported reaction model
-------------------------------------

* Give imported forward and reverse SBML reactions distinct source identities.
* Honor SBML Level 2's omitted reversible default; explicitly report unsupported
  or incomplete Level 3 input and fixed-integer stoichiometry restrictions.
* Raise on incomplete siphon enumeration by default; explicit partial-result
  requests must warn and cannot certify persistence.
* Regressions: reversible A <=> B survives digraph conversion and round trip;
  BIOMD0000000637 reproduces exact rank 9, deficiency 9, and three pools under
  the declared import convention; node limits never masquerade as completion.

2. Make the scientific controls explicit
---------------------------------------

* Use an unambiguous toy catalyst for the retained/cancelled participant control.
  Keep KEGG's generic polyphosphate identifier out of catalyst evidence.
* For the prescribed glycolytic flux, test glucose only, ATP-only additions,
  incomplete redox/phosphate additions, and a fully supplied inventory. Replay
  every positive certificate, and record final inventories and exact deficits.
* For formose, record the two fixed route definitions and their two/three
  glycolaldehyde requirements. Identify an initially empty initiator siphon for
  formaldehyde-only input without confusing it with the marked M3/M4 pools.
* Describe SVD residuals and naive scaled rounding honestly; do not infer failure
  of floating-point nullspace methods from an arbitrary rounding heuristic.

3. Deliver the complete reproduction code
----------------------------------------

* Put the DPO timing driver, rule/seed inputs, workers, and equality audits in
  Experiments/CRN, independent of manuscript files and archived executables.
* Verify the sdist contains the native controller, CMake files, experiments,
  fixtures, and environment constraints. Build native code from packaged source.
* Record the actual Python/packages/compiler/Boost/RDKit SDK, source hashes,
  thread and CPU policy, search limits, warm-ups, repeats, and failed runs.
* Provide a smoke profile and a full profile. A native-only run cannot validate
  equality or speed against an absent comparator. Archived comparator output
  must remain explicitly labeled as archived.

4. Regenerate evidence and the manuscript
----------------------------------------

* Run the full native study suite with checksum-verified cached BioModels input.
* Cross-check BioModels CRNT quantities directly from XML using rational algebra.
* Regenerate the participant-convention audit and glycolysis certificate replay.
* Rerun DPO expansion at the six stated depths with the declared repetitions;
  rerun external comparators when their isolated environments are available.
* Replace corrected historical rows with fresh reports; never relabel historical
  timings or comparator outputs as 1.7.0 measurements.
* Synchronize generated tables/macros, then render figures and both documents.
  Repair the standalone deficiency example, route scopes, siphon interpretation,
  numerical-comparison claims, figure description, and abstract counting.

5. Release checks
-----------------

* CRN regression tests pass, including native expansion and new import controls.
* Full saved-report validation checks checksums, source/version consistency,
  workload completeness, and certificate replay.
* The sdist builds and contains every required source and input. Installed-wheel
  analysis and a fresh source-build expansion smoke check pass where applicable.
* Manuscript checks pass and both PDFs build with resolved references.
* Remaining external comparisons, release/tag, DOI, and publication steps stay
  explicit pending items until actually completed; no commit/push is automatic.

Execution record
----------------

Work started 10 November 2026. Detailed commands, results, environment and any
remaining requirements will be recorded in analysis/release-1.7.0-status.json.

Completed locally on 10 November 2026: all ten full-profile studies, a live
CRNT4SBML comparison, the 31-case structural audit, all 30 participant-convention
variants, ten MolFoundry fixed-route decisions, and 414 matched DPO expansions.
The full CRN suite passed 1,048 tests and 487 subtests; the added real-model
search-budget regression passed separately. Both distributions built and passed
archive/native checks, including DPO builds from the unpacked sdist and wheel.
The regenerated manuscript and supplement build with resolved references; no
unfinished revision notes remain. Publication steps below remain with the user.

Reproduction commands
---------------------

Run the science from the SynKit source root, with a C++ toolchain, Boost headers
and the RDKit C++ SDK configured. The Python constraints and isolated comparator
environment specifications are in ``Experiments/CRN/environments``. BioModels
inputs are checksum-verified against the packaged manifest; omit ``--offline``
on the initial download. Keep CRNT4SBML in its separate legacy environment::

    python Experiments/CRN/run_all.py --with-biomodels --offline \
      --crnt4sbml-python /path/to/crnt4sbml-env/bin/python \
      --output-dir /tmp/crn-170 --summary /tmp/crn-170/summary.json
    python Experiments/CRN/tools/check_reports.py \
      /tmp/crn-170/summary.json --require-all
    python Experiments/CRN/studies/participant_conventions.py --offline \
      --output /tmp/crn-170/participant-conventions.json
    python Experiments/CRN/studies/formose_dpo.py --external \
      --mod-python /path/to/mod-env/bin/python \
      --molfoundry-python /path/to/molfoundry-env/bin/python \
      --depths 1 2 3 4 5 6 --warmups 3 --repeats 20 \
      --output-dir /tmp/formose-170

Set ``SYNKIT_ROOT=/path/to/SynKit`` when the two repositories are not in the
original directory layout. From the manuscript root, the orchestrator stages the full study reports and
31-case independent structural audit. ``--install`` copies successful reports
under the manuscript filenames; historical comparator reuse requires the
separate, explicit ``--allow-recorded-crnt4sbml`` option::

    python analysis/scripts/regenerate_evidence.py \
      --work-dir /tmp/paper-170 \
      --crnt4sbml-python /path/to/crnt4sbml-env/bin/python --install
    python analysis/scripts/sync_manuscript.py
    make all

The DPO output directory supplies ``analysis/formose-expansion/dpo``. Its complete
JSON reports, exact snapshots, warm-up and measured schedules must be copied as
a unit before regenerating the DPO figure. Independent fixed-route comparison
is ``analysis/scripts/compare_formose_routes.py`` in the MolFoundry interpreter.
The document-specific figure and LaTeX code belongs to this manuscript archive.

Release order after review
--------------------------

1. Review the source diffs in both repositories and the recorded status file.
2. Commit the reviewed SynKit sources, including the full ``Experiments/CRN``
   tree and native DPO sources; remove the old singular experiment paths and
   runtime benchmark package. Preserve unrelated work as separate changes.
3. Create and push the intended ``v1.7.0`` release with its tested source and
   package artifacts. Archive the release, manuscript reproduction assets and
   checksum-verified BioModels inputs separately from the Python distribution.
4. Replace the explicitly unpublished snapshot citation with the real release
   URL/DOI, verify source hashes against the tag, and rerun the document build.

No release, tag, commit, push or DOI creation is performed by this repair task.
