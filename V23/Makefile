all: build/main.pdf

## hier Python-Skripte:
#build/plot.pdf: plot.py matplotlibrc header-matplotlib.tex | build
#	TEXINPUTS=$$(pwd): python plot.py
#


build/Polar_04.pdf: PolarplotsWasserstoffatom.py matplotlibrc header-matplotlib.tex | build 
	TEXINPUTS=$$(pwd): python PolarplotsWasserstoffatom.py

build/ausgleichsgerade.pdf: WasserstoffatomRing.py matplotlibrc header-matplotlib.tex | build 
	TEXINPUTS=$$(pwd):	python WasserstoffatomRing.py

build/resonanz1.pdf: PolarRing9mm.py matplotlibrc header-matplotlib.tex | build 
	TEXINPUTS=$$(pwd):	python PolarRing9mm.py	

build/frequenz_blende.pdf: WasserstoffmolekülBlende.py matplotlibrc header-matplotlib.tex | build 
	TEXINPUTS=$$(pwd):	python WasserstoffmolekülBlende.py	

build/blende_polar3.pdf: PolarWasserstoffmolekülBlende.py matplotlibrc header-matplotlib.tex | build 
	TEXINPUTS=$$(pwd):	python PolarWasserstoffmolekülBlende.py		
## hier weitere Abhängigkeiten für build/main.pdf deklarieren:
build/main.pdf: build/Polar_04.pdf build/ausgleichsgerade.pdf build/resonanz1.pdf build/frequenz_blende.pdf build/blende_polar3.pdf

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
