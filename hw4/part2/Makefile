MPICXX := mpicxx
CXXFLAGS := -std=c++17 -Wall

EXE := matmul

all: CXXFLAGS += -O3
all: $(EXE)

debug: CXXFLAGS += -g -Wextra -Wno-deprecated-declarations -fsanitize=address -fdiagnostics-color=always
debug: $(EXE)

$(EXE): main.cc matrix.cc
	$(MPICXX) $(CXXFLAGS) $^ -o $@

.PHONY: clean
clean:
	rm -f $(EXE)