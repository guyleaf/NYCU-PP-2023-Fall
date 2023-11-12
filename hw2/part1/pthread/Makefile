C := clang
CFLAGS := -m64 -std=gnu17 -Iinclude -Wall

EXE := pi.out
SRCS := $(wildcard src/*.c)
OBJS := $(SRCS:%.c=%.o)

ifeq (/usr/bin/clang-11, $(wildcard /usr/bin/clang-11*))
    C=clang-11
endif

all: CFLAGS += -O3 -mavx2 -march=native -ffast-math
all: $(EXE)

debug: CFLAGS += -g -Wextra -Wno-deprecated-declarations -fsanitize=address -fdiagnostics-color=always
debug: $(EXE)

$(EXE):$(OBJS)
	$(C) $(CFLAGS) $^ -o $@ -lpthread

%.o:%.c
	$(C) $(CFLAGS) -c $^ -o $@

.PHONY: clean
clean:
	rm -f $(EXE) $(OBJS) *.out