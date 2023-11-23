VERBOSE := 0
SANITIZE := 0
CPPFLAGS := -I../ -std=c++17 -fopenmp -O3

all: default grade

debug: CPPFLAGS += -g -Wextra -Wno-deprecated-declarations
ifeq ($(VERBOSE),1)
	CPPFLAGS += -DDEBUG -DVERBOSE
endif
ifeq ($(SANITIZE),1)
	CPPFLAGS += -fsanitize=address -fdiagnostics-color=always
endif
debug: default grade

default: main.cpp bfs.cpp
	g++ $(CPPFLAGS) -o bfs.out main.cpp bfs.cpp ../common/graph.cpp ref_bfs.a
grade: grade.cpp bfs.cpp
	g++ $(CPPFLAGS) -o bfs_grader.out grade.cpp bfs.cpp ../common/graph.cpp ref_bfs.a
clean:
	rm -rf bfs_grader.out bfs.out  *~ *.*~