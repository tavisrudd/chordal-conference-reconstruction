export SOURCE_DATE_EPOCH = 1767225600
export FORCE_SOURCE_DATE = 1

TEXSHELL ?= nix develop .\#manuscript --command
LATEXMK ?= $(TEXSHELL) latexmk
LATEXMK_FLAGS ?= -xelatex -interaction=nonstopmode -halt-on-error
SOURCE := golden_companion_reconstruction.tex

.PHONY: all check evidence warnings clean distclean

all: golden_companion_reconstruction.pdf

check: evidence golden_companion_reconstruction.pdf warnings

evidence:
	python3 verification/evidence/paper_ii_chordal_axis.py --check

golden_companion_reconstruction.pdf: $(SOURCE)
	$(LATEXMK) $(LATEXMK_FLAGS) $(SOURCE)

warnings: golden_companion_reconstruction.pdf
	@if grep -En 'Overfull|Underfull|LaTeX Warning|Package .* Warning|undefined references|Citation .* undefined' golden_companion_reconstruction.log; then \
		exit 1; \
	fi

clean:
	$(LATEXMK) -c $(SOURCE)

distclean:
	$(LATEXMK) -C $(SOURCE)
