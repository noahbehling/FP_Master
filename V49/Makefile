all: build/main.pdf

## hier Python-Skripte:
#build/plot.pdf: plot.py matplotlibrc header-matplotlib.tex | build
#	TEXINPUTS=$$(pwd): python plot.py
#
	
## hier weitere Abhängigkeiten für build/main.pdf deklarieren:
build/eh.pdf: diff.py matplotlibrc header-matplotlib.tex | build 
	TEXINPUTS=$$(pwd): python diff.py
build/uh.pdf: t1.py matplotlibrc header-matplotlib.tex | build 
	TEXINPUTS=$$(pwd): python t1.py	

build/main.pdf: build/eh.pdf build/uh.pdf

build/main.pdf: FORCE | build
	  TEXINPUTS=build: \
	  BIBINPUTS=build: \
	  max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=nonstopmode \
	  --halt-on-error \
	main.tex


build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean
