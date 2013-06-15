
CXX=g++

LIB_FLAGS = -l Goptical $(EXTRA_LIB_FLAGS)

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

plots: src/EyePlots.py
	./src/EyePlots.py

.PHONY: clean

clean:
	rm -f eye

