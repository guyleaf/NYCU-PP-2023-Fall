CXX := g++
CXXFLAGS := -m64 -std=c++17 -Wall

APP_NAME=pi.out

ifeq (/usr/bin/g++-10,$(wildcard /usr/bin/g++-10*))
    CXX=g++-10
endif

all: CXXFLAGS += -O3
all: $(APP_NAME)

debug: CXXFLAGS += -g -Wextra -Wno-deprecated-declarations -fsanitize=address -fdiagnostics-color=always
debug: $(APP_NAME)

$(APP_NAME): utility.hpp pi.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) utility.hpp pi.cpp -o $(APP_NAME) -lpthread

clean:
	rm -f *.o *.out