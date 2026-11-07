SYNKIT_ROOT ?= $(abspath ../../SynEco/SynKit)
export SYNKIT_ROOT
PYTHON ?= python
MPL_CACHE ?= /tmp/synkit-crn-matplotlib
# Figures use Latin Modern Roman (TikZ, via fontspec) and CMU Serif (matplotlib);
# only the visual abstract keeps Arial and colour.
LATEXMK_FLAGS := -xelatex -interaction=nonstopmode -halt-on-error
ANALYSIS_EVIDENCE := $(wildcard analysis/crn-*.json)
REVIEW_ASSETS := Figure/visual_abstract.pdf \
	Figure/representation_sensitivity.pdf Figure/scaling.pdf Figure/formose.pdf
REPLAY_ASSETS := Figure/glycolysis_replay.pdf
FIGURE_SOURCES := Figure/synstyle.tex Figure/arch.tex Figure/glycolysis.tex \
	Figure/formose_chemistry.tex Figure/deficiency.tex Figure/analysis_overview.pdf \
	Figure/performance_diagnostics.pdf Figure/formose_dpo.pdf $(REVIEW_ASSETS) $(REPLAY_ASSETS)

.PHONY: all final-main evidence-check document-check sync-values glycolysis-analysis review-experiments clean
.NOTPARALLEL: all
all: sync-values SI/si.pdf final-main evidence-check


Figure/formose_dpo.pdf: analysis/scripts/formose_dpo_assets.py analysis/scripts/result_style.py analysis/formose-expansion/dpo/comparison.json
	PYTHONPATH=analysis/scripts MPLCONFIGDIR=$(MPL_CACHE) $(PYTHON) analysis/scripts/formose_dpo_assets.py

sync-values:
	PYTHONPATH=analysis/scripts MPLCONFIGDIR=$(MPL_CACHE) $(PYTHON) analysis/scripts/sync_manuscript.py

review-experiments:
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 $(PYTHON) $(SYNKIT_ROOT)/Experiment/CRN/studies/representation.py --output analysis/crn-review-experiments.json

glycolysis-analysis:
	PYTHONPATH=$(SYNKIT_ROOT):analysis/scripts $(PYTHON) analysis/scripts/glycolysis_analysis.py --backend cpp --json analysis/glycolysis_analysis.json
	$(MAKE) sync-values

$(REVIEW_ASSETS) &: analysis/scripts/review_figures.py analysis/scripts/visual_abstract.py analysis/scripts/results_panels.py analysis/scripts/biochemical_schemes.py analysis/scripts/result_style.py $(ANALYSIS_EVIDENCE)
	PYTHONPATH=analysis/scripts MPLCONFIGDIR=$(MPL_CACHE) $(PYTHON) analysis/scripts/review_figures.py --evidence-dir analysis --output-dir Figure

Figure/analysis_overview.pdf: analysis/scripts/analysis_overview.py analysis/scripts/result_style.py $(ANALYSIS_EVIDENCE)
	PYTHONPATH=analysis/scripts MPLCONFIGDIR=$(MPL_CACHE) $(PYTHON) analysis/scripts/analysis_overview.py --evidence-dir analysis --output Figure/analysis_overview.pdf

$(REPLAY_ASSETS): analysis/scripts/glycolysis_replay.py analysis/scripts/result_style.py analysis/glycolysis_analysis.json
	PYTHONPATH=analysis/scripts MPLCONFIGDIR=$(MPL_CACHE) $(PYTHON) analysis/scripts/glycolysis_replay.py --output-dir Figure --input analysis/glycolysis_analysis.json

Figure/performance_diagnostics.pdf: analysis/scripts/performance_diagnostics.py analysis/crn-scaling.json
	PYTHONPATH=analysis/scripts MPLCONFIGDIR=$(MPL_CACHE) $(PYTHON) $< --evidence analysis/crn-scaling.json --output $@

Figure/structural_landscape.pdf: analysis/scripts/structural_landscape.py $(ANALYSIS_EVIDENCE)
	PYTHONPATH=analysis/scripts MPLCONFIGDIR=$(MPL_CACHE) $(PYTHON) $< --evidence-dir analysis --output $@

main.pdf: main.tex paper_layout.tex refs.bib $(FIGURE_SOURCES)
	latexmk $(LATEXMK_FLAGS) main.tex

SI/si.pdf: SI/si.tex paper_layout.tex main.pdf $(FIGURE_SOURCES)
	latexmk $(LATEXMK_FLAGS) -cd SI/si.tex

final-main: SI/si.pdf
	latexmk -g $(LATEXMK_FLAGS) main.tex

evidence-check:
	MPLCONFIGDIR=$(MPL_CACHE) $(PYTHON) verify_revision.py

document-check:
	MPLCONFIGDIR=$(MPL_CACHE) $(PYTHON) verify_revision.py --documents-only

clean:
	latexmk -c main.tex
	latexmk -c -cd SI/si.tex
