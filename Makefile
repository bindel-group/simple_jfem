LDOC= 	ldoc

DOCS=	docs/quadrules.qmd \
	docs/shapes.qmd \
	docs/mesh.qmd \
	docs/assemble.qmd \
	docs/fem.qmd \
	docs/element.qmd 

.PHONY: all doc test clean

all:

docs/%.qmd: src/%.jl
	$(LDOC) -l sh -highlight julia -p quarto -o $@ $^

docs/index.html: docs/index.qmd $(DOCS)
	( cd docs ; quarto render index.qmd --to html )

docs/index.pdf: docs/index.qmd $(DOCS)
	( cd docs ; quarto render index.qmd --to pdf )

doc: docs/index.pdf docs/index.html

publish: docs/index.pdf docs/index.html
	( cd docs ; quarto publish gh-pages index.qmd )

test:
	( for f in test/test*.jl ; do julia $$f ; done )

clean:
	rm -f src/*~ test/*~ Makefile.bak
	rm -f $(DOCS)
