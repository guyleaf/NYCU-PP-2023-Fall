CXX := g++
CXXFLAGS := -I./common -std=c++17 -Wall

ifeq (/usr/bin/g++-10,$(wildcard /usr/bin/g++-10*))
    CXX=g++-10
endif

all: CXXFLAGS += -O3
all: myexp

debug: CXXFLAGS += -g -Wextra -Wno-deprecated-declarations -fsanitize=address -fdiagnostics-color=always
debug: myexp

logger.o: logger.cpp logger.h PPintrin.h PPintrin.cpp def.h
	$(CXX) $(CXXFLAGS) -c logger.cpp

PPintrin.o: PPintrin.cpp PPintrin.h logger.cpp logger.h def.h
	$(CXX) $(CXXFLAGS) -c PPintrin.cpp

myexp: PPintrin.o logger.o main.cpp serialOP.cpp vectorOP.cpp 
	g++ $(CXXFLAGS) logger.o PPintrin.o main.cpp serialOP.cpp vectorOP.cpp -o myexp

clean:
	rm -f *.o *.s myexp *~
