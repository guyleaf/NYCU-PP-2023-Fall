MPICXX := mpicxx
CXXFLAGS := -Wall

EXES := hello.out pi_block_linear.out pi_block_tree.out pi_nonblock_linear.out pi_gather.out pi_reduce.out pi_one_side.out

all: CXXFLAGS += -O3
all: $(EXES)

debug: CXXFLAGS += -g -Wextra -Wno-deprecated-declarations -fsanitize=address -fdiagnostics-color=always
debug: $(EXES)

%.out:src/%.cc
	$(MPICXX) $(CXXFLAGS) $^ -o $@

.PHONY: clean
clean:
	rm -f $(EXES)