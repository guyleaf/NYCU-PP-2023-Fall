SANITIZE := 0

ifndef PANEL_WIDTH
	PANEL_WIDTH := 32
endif

MPICXX := mpicxx
CXXFLAGS := -std=c++17 -Wall -DPANEL_WIDTH=$(PANEL_WIDTH)

EXE := matmul

all: CXXFLAGS += -O3
all: $(EXE)

debug: CXXFLAGS += -g -DDEBUG -Wextra -Wno-deprecated-declarations
ifeq ($(SANITIZE),1)
        CXXFLAGS += -fsanitize=address -fdiagnostics-color=always
endif
debug: $(EXE)

$(EXE): main.cc matrix.cc
	$(MPICXX) $(CXXFLAGS) $^ -o $@

.PHONY: clean
clean:
	rm -f $(EXE)
