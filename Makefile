
CXX=g++

LIB_FLAGS = -lGoptical -larmadillo $(EXTRA_LIB_FLAGS)

OPT = -O2

CXXFLAGS = $(DEBUG) $(FINAL) $(OPT) $(EXTRA_OPT)

EYEFILES = \
	src/eye_eye.cpp \
	src/eye_analysis.cpp \
	src/main.cpp

GENERATED_FILES = eye

all: $(GENERATED_FILES)

eye: $(EYEFILES)
	$(CXX) $(CXXFLAGS)  -o $@ $^ $(LIB_FLAGS)

install: $(GENERATED_FILES)
	sudo rm /usr/local/bin/$(GENERATED_FILES)
	sudo cp eye /usr/local/bin/$(GENERATED_FILES)

python: eye cython/setup.py
	./cython/setup.py build_ext --inplace

plots: src/EyePlots.py
	./src/EyePlots.py -m $(args)

.PHONY: clean

clean:
	rm -f eye

