CPPFLAGS := -I../ -std=c++17 -fopenmp -O3

all: default grade

debug: CPPFLAGS += -g
debug: default grade

default: page_rank.cpp main.cpp
	g++ $(CPPFLAGS) -o pr.out main.cpp page_rank.cpp ../common/graph.cpp ref_pr.a
grade: page_rank.cpp grade.cpp
	g++ $(CPPFLAGS) -o pr_grader.out grade.cpp page_rank.cpp ../common/graph.cpp ref_pr.a
clean:
	rm -rf pr.out pr_grader.out *~ *.*~
