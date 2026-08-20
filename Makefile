export SOURCE_DATE_EPOCH = 1767225600
export FORCE_SOURCE_DATE = 1

TEXSHELL ?= nix develop .\#manuscript --command
LATEXMK ?= $(TEXSHELL) latexmk
LATEXMK_FLAGS ?= -xelatex -interaction=nonstopmode -halt-on-error
PYTHON ?= nix shell nixpkgs\#python3 -c python3
SOURCE := chordal_conference_reconstruction.tex

.PHONY: all check evidence manuscript warnings clean distclean

all: manuscript

check: evidence manuscript warnings

evidence:
	$(PYTHON) verification/evidence/paper_ii_chordal_axis.py --check
	$(PYTHON) verification/evidence/conference_node_completeness.py --check

manuscript: $(SOURCE)
	$(LATEXMK) $(LATEXMK_FLAGS) $(SOURCE)

warnings: manuscript
	@if grep -En 'Overfull|Underfull|LaTeX Warning|Package .* Warning|undefined references|Citation .* undefined' chordal_conference_reconstruction.log; then \
		exit 1; \
	fi

clean:
	$(LATEXMK) -c $(SOURCE)

distclean:
	$(LATEXMK) -C $(SOURCE)
