# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = diffpy.structure
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Manual publishing to the gh-pages branch

GITREPOPATH = $(shell cd $(CURDIR) && git rev-parse --git-dir)
GITREMOTE = origin
GITREMOTEURL = $(shell git config --get remote.$(GITREMOTE).url)
GITLASTCOMMIT = $(shell git rev-parse --short HEAD)

publish:
	@test -d build/html || \
	    ( echo >&2 "Run 'make html' first!"; false )
	git show-ref --verify --quiet refs/heads/gh-pages || \
	    git branch --track gh-pages $(GITREMOTE)/gh-pages
	test -d build/gh-pages || \
	    git clone -s -b gh-pages $(GITREPOPATH) build/gh-pages
	cd build/gh-pages && \
	    git pull $(GITREMOTEURL) gh-pages
	rsync -acv --delete --exclude=.git --exclude=.rsync-exclude \
	    --exclude-from=build/gh-pages/.rsync-exclude \
	    --link-dest=$(CURDIR)/build/html build/html/ build/gh-pages/
	cd build/gh-pages && \
	    git add --all . && \
	    git diff --cached --quiet || \
	    git commit -m "Sync with the source at $(GITLASTCOMMIT)."
	cd build/gh-pages && \
	    git push origin gh-pages
