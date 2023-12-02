MPICXX := mpicxx
CXXFLAGS := -Wall

PART1_EXES := hello.out pi_block_linear.out pi_block_tree.out pi_nonblock_linear.out pi_gather.out pi_reduce.out pi_one_side.out
PART1_EXES := $(addprefix part1/, $(PART1_EXES))

all: part1 part2

part1: $(PART1_EXES)
part2:
	$(MAKE) -C $@ $(MAKECMDGOALS)

%.out:%.cc
	$(MPICXX) $(CXXFLAGS) $^ -o $@

.PHONY: clean part1 part2
clean: part2
	rm -f $(PART1_EXES)
