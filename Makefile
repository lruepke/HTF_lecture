# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# make cases as archieve files
casesPackage = mkdir -p ${BUILDDIR}/$@/cases \
				&& zip -r ${BUILDDIR}/$@/cases/$2.zip source/lectures/$1/cases/$2 -x source/lectures/$1/cases/$2/*00\* source/lectures/$1/cases/$2/jupyter/data*\* \

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	@if [ $@ = "html" ]; then\
		$(call casesPackage,L04,Jupp_Schultz) ;\
	fi