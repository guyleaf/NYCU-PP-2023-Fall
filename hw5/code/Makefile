CXX = g++-8
NVCC = nvcc -ccbin=$(CXX)
CXXFLAGS = -I./common -Iobjs/ -O3 -std=c++17 -Wall -g -fPIC -lm

KERNELS = kernel1 kernel2 kernel3 kernel4
APP_NAME = mandelbrot
OBJDIR = objs
COMMONDIR = ./common

CUDA_LINK_FLAGS = -rdc=true -gencode=arch=compute_61,code=sm_61 -Xcompiler '-fPIC'
CUDA_COMPILE_FLAGS = --device-c -gencode=arch=compute_61,code=sm_61 -Xcompiler '-fPIC' -g -O3

PPM_CXX = $(COMMONDIR)/ppm.cpp
PPM_OBJ = $(addprefix $(OBJDIR)/, $(subst $(COMMONDIR)/,, $(PPM_CXX:.cpp=.o)))

EXES = $(addsuffix .out,$(addprefix $(APP_NAME)_,$(KERNELS)))

default: $(EXES)

.PHONY: dirs clean

dirs:
		/bin/mkdir -p $(OBJDIR)/

clean:
		/bin/rm -rf $(OBJDIR) *.ppm *~ $(EXES)

OBJS=$(OBJDIR)/main.o $(OBJDIR)/mandelbrotSerial.o $(OBJDIR)/mandelbrotThread.o $(PPM_OBJ)

$(APP_NAME)_%.out: dirs $(OBJS) $(OBJDIR)/%.o
		$(NVCC) ${CUDA_LINK_FLAGS} -o $@ $(OBJS) $(OBJDIR)/$*.o mandelbrotThreadRef.a

$(OBJDIR)/kernel%.o : kernel%.cu kernel.h
	${NVCC} ${CUDA_COMPILE_FLAGS} -c kernel$*.cu -o $@

$(OBJDIR)/%.o: %.cpp
		$(CXX) $< $(CXXFLAGS) -c -o $@

$(OBJDIR)/%.o: $(COMMONDIR)/%.cpp
	$(CXX) $< $(CXXFLAGS) -c -o $@

$(OBJDIR)/main.o: $(COMMONDIR)/CycleTimer.h kernel.h
