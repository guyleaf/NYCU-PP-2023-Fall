CPPFLAGS := -I../ -std=c++17 -fopenmp -O3

all: default grade

debug: CPPFLAGS += -g -DDEBUG -DVERBOSE
debug: default grade

default: main.cpp bfs.cpp
	g++ $(CPPFLAGS) -o bfs.out main.cpp bfs.cpp ../common/graph.cpp ref_bfs.a
grade: grade.cpp bfs.cpp
	g++ $(CPPFLAGS) -o bfs_grader.out grade.cpp bfs.cpp ../common/graph.cpp ref_bfs.a
clean:
	rm -rf bfs_grader.out bfs.out  *~ *.*~
