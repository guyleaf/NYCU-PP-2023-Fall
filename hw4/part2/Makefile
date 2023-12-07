SANITIZE := 0

MPICXX := mpicxx
CXXFLAGS := -std=c++17 -Wall

EXE := matmul

ifdef PANEL_WIDTH
	CXXFLAGS += -DPANEL_WIDTH=$(PANEL_WIDTH)
endif

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
