SYNKIT_ROOT ?= ../..
PYTHON ?= python
MPL_CACHE ?= /tmp/synkit-crn-matplotlib
ANALYSIS_EVIDENCE := $(wildcard analysis/crn-*.json)
REVIEW_ASSETS := Figure/visual_abstract.pdf \
	Figure/representation_sensitivity.pdf Figure/scaling.pdf Figure/formose.pdf \
	Figure/glycolysis_donors.pdf
REPLAY_ASSETS := Figure/glycolysis_replay.pdf
FIGURE_SOURCES := $(wildcard Figure/*.tex) Figure/analysis_overview.pdf \
	Figure/performance_diagnostics.pdf Figure/structural_landscape.pdf $(REVIEW_ASSETS) $(REPLAY_ASSETS)

.PHONY: all final-main evidence-check glycolysis-analysis review-experiments clean
all: SI/si.pdf final-main evidence-check

review-experiments:
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 $(PYTHON) $(SYNKIT_ROOT)/Experiment/CRN/review_experiments.py --output analysis/crn-review-experiments.json

glycolysis-analysis:
	PYTHONPATH=$(SYNKIT_ROOT):analysis/scripts $(PYTHON) analysis/scripts/glycolysis_analysis.py --json Figure/glycolysis_analysis.json --tex Figure/glycolysis_analysis.tex

$(REVIEW_ASSETS) &: analysis/scripts/review_figures.py analysis/scripts/visual_abstract.py analysis/scripts/results_panels.py analysis/scripts/biochemical_schemes.py analysis/scripts/result_style.py $(ANALYSIS_EVIDENCE)
	PYTHONPATH=analysis/scripts MPLCONFIGDIR=$(MPL_CACHE) $(PYTHON) analysis/scripts/review_figures.py --evidence-dir analysis --output-dir Figure

Figure/analysis_overview.pdf: analysis/scripts/analysis_overview.py analysis/scripts/result_style.py $(ANALYSIS_EVIDENCE)
	PYTHONPATH=analysis/scripts MPLCONFIGDIR=$(MPL_CACHE) $(PYTHON) analysis/scripts/analysis_overview.py --evidence-dir analysis --output Figure/analysis_overview.pdf

$(REPLAY_ASSETS): analysis/scripts/glycolysis_replay.py analysis/scripts/result_style.py Figure/glycolysis_analysis.json
	PYTHONPATH=analysis/scripts MPLCONFIGDIR=$(MPL_CACHE) $(PYTHON) analysis/scripts/glycolysis_replay.py --output-dir Figure --json Figure/glycolysis_analysis.json

Figure/performance_diagnostics.pdf: analysis/scripts/performance_diagnostics.py analysis/crn-scaling.json
	PYTHONPATH=analysis/scripts MPLCONFIGDIR=$(MPL_CACHE) $(PYTHON) $< --evidence analysis/crn-scaling.json --output $@

Figure/structural_landscape.pdf: analysis/scripts/structural_landscape.py $(ANALYSIS_EVIDENCE)
	PYTHONPATH=analysis/scripts MPLCONFIGDIR=$(MPL_CACHE) $(PYTHON) $< --evidence-dir analysis --output $@

main.pdf: main.tex paper_layout.tex refs.bib $(FIGURE_SOURCES)
	latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex

SI/si.pdf: SI/si.tex paper_layout.tex main.pdf $(FIGURE_SOURCES)
	latexmk -pdf -interaction=nonstopmode -halt-on-error -cd SI/si.tex

final-main: SI/si.pdf
	latexmk -g -pdf -interaction=nonstopmode -halt-on-error main.tex

evidence-check:
	MPLCONFIGDIR=$(MPL_CACHE) $(PYTHON) verify_revision.py

clean:
	latexmk -c main.tex
	latexmk -c -cd SI/si.tex
