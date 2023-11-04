C := gcc
CFLAGS := -m64 -std=gnu17 -Iinclude -Wall

EXE := pi.out
SRCS := $(wildcard src/*.c)
OBJS := $(SRCS:%.c=%.o)

ifeq (/usr/bin/gcc-10, $(wildcard /usr/bin/gcc-10*))
    C=gcc-10
endif

all: CFLAGS += -O3 -mavx2
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