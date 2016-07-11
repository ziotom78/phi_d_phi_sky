FPC := fpc
NOWEAVE := noweave
NOTANGLE := notangle
TEXI2PDF := lualatex
CPIF := cpif
PLANCKLIBPATH := ../planck
FPCOPTS = -g -gl -gv -gw -B -Co -Cr -vi -vh
FIGURES = quantile_processing.pdf quantile_mem_layout.pdf

.phony: all

all: phisky phid phiskyandd.pdf

phisky: phisky.pas
	$(FPC) $(FPCOPTS) -Fo$(PLANCKLIBPATH) -Fu$(PLANCKLIBPATH) \
		-k-lmpi -k-ldl -k-lhwloc $<

phid: phid.pas
	$(FPC) $(FPCOPTS) -Fo$(PLANCKLIBPATH) -Fu$(PLANCKLIBPATH) \
		-k-lmpi -k-ldl -k-lhwloc $<

phiskyandd.pdf: phiskyandd.tex $(FIGURES)
	$(TEXI2PDF) $<

phisky.pas: phiskyandd.nw
	$(NOTANGLE) -R$@ -L'{line %L "%F"}' $^ | $(CPIF) $@

phid.pas: phiskyandd.nw
	$(NOTANGLE) -R$@ -L'{line %L "%F"}' $^ | $(CPIF) $@

%.tex: %.nw
	$(NOWEAVE) -n -delay -index $< | $(CPIF) $@

%.pdf: %.svg
	inkscape --export-pdf=$@ $<
